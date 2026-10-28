#!/usr/bin/env python3
"""Fusion V3: Token-wise Cross-Attention Ortho Concat vs Gated MoE Add."""
import json, math, time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy import stats
try:
    from sklearn.metrics import roc_auc_score
except ImportError:  # foundry env lacks sklearn; only per_bin_analysis needs it
    roc_auc_score = None
import torch, torch.nn as nn

DEVICE = torch.device("cuda")
DM, NH, BS = 256, 8, 64
EPOCHS, LR, WD = 300, 1e-4, 0.05
FAESM_DIM, FINAL_DIM = 1280, 128
SPLITS, INITS = 5, 3
SPLIT_SEEDS = [100, 101, 102, 103, 104]
DROPOUT = 0.1
CKPT_DIR = Path("J:/conda_envs/foundry/DMS_Project/logs/ckpts")
CKPT_DIR.mkdir(parents=True, exist_ok=True)

class PosEnc(nn.Module):
    def __init__(self,d,m=500):
        super().__init__(); pe=torch.zeros(m,d); pos=torch.arange(m).float().unsqueeze(1)
        div=torch.exp(torch.arange(0,d,2).float()*(-math.log(10000)/d))
        pe[:,0::2]=torch.sin(pos*div); pe[:,1::2]=torch.cos(pos*div)
        self.register_buffer('pe',pe.unsqueeze(0))
    def forward(self,x): return x+self.pe[:,:x.size(1)]

class Encoder(nn.Module):
    def __init__(self,d_in):
        super().__init__()
        self.in_proj=nn.Sequential(nn.Linear(d_in,DM),nn.GELU(),nn.Dropout(DROPOUT))
        self.pe=PosEnc(DM)
        el=nn.TransformerEncoderLayer(d_model=DM,nhead=NH,dim_feedforward=DM*4,dropout=DROPOUT,batch_first=True)
        self.enc=nn.TransformerEncoder(el,num_layers=3)
        self.out_proj=nn.Sequential(nn.Linear(DM,DM*2),nn.GELU(),nn.Dropout(DROPOUT),
                                     nn.Linear(DM*2,FINAL_DIM),nn.GELU())
    def forward(self,x):
        x=self.in_proj(x); x=self.pe(x); x=self.enc(x)
        return self.out_proj(x)

class AttnPool(nn.Module):
    def __init__(self,d):
        super().__init__(); self.w=nn.Linear(d,1)
    def forward(self,x):  # (B,L,D)
        a=torch.softmax(self.w(x),dim=1)  # (B,L,1)
        return (x*a).sum(dim=1)  # (B,D)

class CrossAttnOrthoConcatFusion(nn.Module):
    """
    1. 2-token Cross-Attention per residue
    2. Split and apply Orthogonal Penalty
    3. Concat in feature dimension
    """
    def __init__(self):
        super().__init__()
        self.z_enc = Encoder(128)
        self.f_enc = Encoder(FAESM_DIM)
        
        # 逐残基的 2-token 交互模块 (128维)
        interact_el = nn.TransformerEncoderLayer(d_model=FINAL_DIM, nhead=NH//2, dim_feedforward=FINAL_DIM*4, dropout=DROPOUT, batch_first=True)
        self.interaction = nn.TransformerEncoder(interact_el, num_layers=1)
        
        self.pe = PosEnc(DM) # 拼接后变成 256维 (DM)
        el = nn.TransformerEncoderLayer(d_model=DM, nhead=NH, dim_feedforward=DM*4, dropout=DROPOUT, batch_first=True)
        self.joint_enc = nn.TransformerEncoder(el, num_layers=2)
        self.pool = AttnPool(DM)
        self.head = nn.Sequential(nn.Linear(DM, DM//2), nn.GELU(), nn.Dropout(DROPOUT), nn.Linear(DM//2, 1))

    def forward(self, xz, xf):
        z = self.z_enc(xz)  # (B, L, 128)
        f = self.f_enc(xf)  # (B, L, 128)
        B, L, D = z.shape
        
        # 沿第2维度堆叠，并压平Batch和Length -> (B*L, 2, 128)
        # 让每个位点上的 z_i 和 f_i 自己做注意力交流
        tokens = torch.stack([z, f], dim=2).view(B*L, 2, D)
        interacted = self.interaction(tokens)
        interacted = interacted.view(B, L, 2, D)
        
        # 交流完毕，重新切分开
        z_out = interacted[:, :, 0, :]
        f_out = interacted[:, :, 1, :]
        
        # 将真正互补的特征拼接 -> (B, L, 256)
        j = torch.cat([z_out, f_out], dim=-1)
        j = self.pe(j); j = self.joint_enc(j); j = self.pool(j)
        pred = self.head(j).squeeze(-1)
        
        # 返回预测值以及切分后的Token，供训练循环施加正交惩罚
        return pred, z_out, f_out

class GatedAddFusion(nn.Module):
    """
    Token-wise Feature Mixture-of-Experts (Gating)
    """
    def __init__(self):
        super().__init__()
        self.z_enc = Encoder(128)
        self.f_enc = Encoder(FAESM_DIM)
        
        # 逐残基动态路由门控 (看一眼256维组合特征，输出128维的混合比例)
        self.gater = nn.Sequential(
            nn.Linear(FINAL_DIM * 2, FINAL_DIM),
            nn.GELU(),
            nn.Linear(FINAL_DIM, FINAL_DIM),
            nn.Sigmoid()
        )
        
        self.pe = PosEnc(FINAL_DIM) # 混合后依然是 128维 (FINAL_DIM)
        el = nn.TransformerEncoderLayer(d_model=FINAL_DIM, nhead=NH//2, dim_feedforward=FINAL_DIM*4, dropout=DROPOUT, batch_first=True)
        self.joint_enc = nn.TransformerEncoder(el, num_layers=2)
        self.pool = AttnPool(FINAL_DIM)
        self.head = nn.Sequential(nn.Linear(FINAL_DIM, FINAL_DIM//2), nn.GELU(), nn.Dropout(DROPOUT), nn.Linear(FINAL_DIM//2, 1))

    def forward(self, xz, xf):
        z = self.z_enc(xz)  # (B, L, 128)
        f = self.f_enc(xf)  # (B, L, 128)
        
        # 动态算出每个维度应该相信结构还是相信进化
        gate = self.gater(torch.cat([z, f], dim=-1))
        j = gate * z + (1.0 - gate) * f
        
        j = self.pe(j); j = self.joint_enc(j); j = self.pool(j)
        return self.head(j).squeeze(-1)

class SinglePredictorV2(nn.Module):
    def __init__(self,d_in):
        super().__init__()
        self.enc=Encoder(d_in)
        self.pool=AttnPool(FINAL_DIM)
        self.head=nn.Sequential(nn.Linear(FINAL_DIM,FINAL_DIM//2),nn.GELU(),nn.Dropout(DROPOUT),nn.Linear(FINAL_DIM//2,1))
    def forward(self,x):
        x=self.enc(x); x=self.pool(x)
        return self.head(x).squeeze(-1)

def _load_one_faesm(args):
    m,faesm_dir,L=args; fp=Path(faesm_dir)/f"{m['mutant_id']}.pt"
    if not fp.exists(): return None
    try:
        d=torch.load(str(fp),map_location='cpu',weights_only=True)
        emb=torch.cat([d[(32,h)].float() for h in range(20)],dim=-1)
        return emb[0,:L,:]
    except: return None

def load_data(protein,L):
    zii_dir=Path(f"K:/output_heads/{protein}/zii"); faesm_dir=Path(f"K:/output_heads/{protein}/faesm")
    with open(f"K:/output_heads/{protein}/metadata.jsonl") as f: meta=[json.loads(l) for l in f if l.strip()]
    zii_files=set(p.stem.replace('_zii','') for p in zii_dir.glob("*.pt"))
    zii_mutants=[m for m in meta if m['mutant_id'] in zii_files]
    Xz=[]
    for m in zii_mutants:
        z=torch.load(str(zii_dir/f"{m['mutant_id']}_zii.pt"),map_location='cpu',weights_only=True)
        if isinstance(z,dict): z=z.get('zii_latest',list(z.values())[0])
        Xz.append(z[0,:L].float().numpy())
    Xz=np.stack(Xz)
    lk='bind_avg' if 'rbd' in protein else 'score'
    y=np.array([m.get(lk,m.get('fitness',m.get('score',0))) for m in zii_mutants])
    sites=np.array([m.get('site_rbd',m.get('position',m.get('site',0))) for m in zii_mutants])
    args_list=[(m,faesm_dir,L) for m in zii_mutants]
    Xf_list=[]; valid_idx=[]
    with ProcessPoolExecutor(max_workers=8) as ex:
        for i,r in enumerate(ex.map(_load_one_faesm,args_list)):
            if r is not None: Xf_list.append(r.numpy()); valid_idx.append(i)
    Xf=np.stack(Xf_list); n=min(len(Xz),len(Xf)); vi=np.array(valid_idx[:n])
    return Xz[vi],Xf[:n],y[vi],sites[vi]

def _train_single(Xt,yt,Xv,yv,model_class,d_in,tag):
    best_val=-1; best_ep=0; best_preds_final=None
    for init_i in range(INITS):
        model=model_class(d_in).to(DEVICE)
        opt=torch.optim.AdamW(model.parameters(),lr=LR,weight_decay=WD)
        sched=torch.optim.lr_scheduler.CosineAnnealingLR(opt,EPOCHS)
        best=-1
        for ep in range(EPOCHS):
            model.train(); idx=torch.randperm(len(Xt))
            for bi in range(0,len(Xt),BS):
                b=idx[bi:bi+BS]; loss=nn.MSELoss()(model(Xt[b].to(DEVICE)),yt[b].to(DEVICE))
                opt.zero_grad(); loss.backward(); opt.step()
            sched.step()
            model.eval()
            with torch.no_grad():
                batches=[model(Xv[bi:min(bi+BS,len(Xv))].to(DEVICE)) for bi in range(0,len(Xv),BS)]
                pv=torch.cat(batches).cpu(); va,_=stats.spearmanr(pv.numpy(),yv.numpy())
                if np.isnan(va): va=-1
                if va>best: best=va; best_ep=ep; best_preds=pv.numpy()
            del batches; torch.cuda.empty_cache()
        print(f"      {tag}.init{init_i+1}: ρ={best:.4f} @ ep{best_ep+1}")
        if best>best_val: best_val=best; best_preds_final=best_preds
    return best_val, best_preds_final

def _train_fusion(Xtz,Xtf,yt,Xvz,Xvf,yv,model_class,tag,cos_penalty=0.05):
    best_val=-1; best_ep=0; best_preds_final=None
    for init_i in range(INITS):
        model=model_class().to(DEVICE)
        opt=torch.optim.AdamW(model.parameters(),lr=LR,weight_decay=WD)
        sched=torch.optim.lr_scheduler.CosineAnnealingLR(opt,EPOCHS)
        best=-1
        for ep in range(EPOCHS):
            model.train(); idx=torch.randperm(len(Xtz))
            for bi in range(0,len(Xtz),BS):
                b=idx[bi:bi+BS]
                out=model(Xtz[b].to(DEVICE),Xtf[b].to(DEVICE))
                if isinstance(out,tuple):
                    pred, z_out, f_out = out
                    # 逐 Token 计算正交惩罚：将 (B, L, D) 压平为 (B*L, D) 进行约束
                    z_n = z_out.reshape(-1, z_out.size(-1))
                    f_n = f_out.reshape(-1, f_out.size(-1))
                    cos = torch.abs((z_n * f_n).sum(-1) / (z_n.norm(dim=-1) * f_n.norm(dim=-1) + 1e-8)).mean()
                    loss = nn.MSELoss()(pred, yt[b].to(DEVICE)) + cos_penalty * cos
                else:
                    loss = nn.MSELoss()(out, yt[b].to(DEVICE))
                opt.zero_grad(); loss.backward(); opt.step()
            sched.step()
            model.eval()
            with torch.no_grad():
                batches=[(model(Xvz[bi:min(bi+BS,len(Xvz))].to(DEVICE),Xvf[bi:min(bi+BS,len(Xvf))].to(DEVICE))) for bi in range(0,len(Xvz),BS)]
                pv=torch.cat([b[0] if isinstance(b,tuple) else b for b in batches]).cpu()
                va,_=stats.spearmanr(pv.numpy(),yv.numpy())
                if np.isnan(va): va=-1
                if va>best: best=va; best_ep=ep; best_preds=pv.numpy()
            del batches; torch.cuda.empty_cache()
        print(f"      {tag}.init{init_i+1}: ρ={best:.4f} @ ep{best_ep+1}")
        if best>best_val: best_val=best; best_preds_final=best_preds
    return best_val, best_preds_final

def per_bin_analysis(preds,true):
    bounds=np.percentile(true,np.linspace(0,100,6))
    within,aucs=[],[]
    for i in range(5):
        m=(true>=bounds[i])&(true<bounds[i+1]) if i<4 else (true>=bounds[i])
        if m.sum()<5: within.append(np.nan)
        else: rho,_=stats.spearmanr(preds[m],true[m]); within.append(rho)
    for i in range(4):
        lo=(true>=bounds[i])&(true<bounds[i+1])
        hi=(true>=bounds[i+1])&(true<bounds[i+2]) if i<3 else (true>=bounds[i+1])
        if lo.sum()<5 or hi.sum()<5: aucs.append(np.nan)
        else:
            y=np.concatenate([np.zeros(lo.sum()),np.ones(hi.sum())])
            p=np.concatenate([preds[lo],preds[hi]])
            aucs.append(roc_auc_score(y,p))
    return within, aucs

def run(protein,L,lk):
    print(f"\n{'='*60}\n{protein.upper()} L={L}\n{'='*60}")
    Xz,Xf,y,sites=load_data(protein,L)
    ys=(y-y.mean())/max(y.std(),1e-8)
    all_pos=np.unique(sites); total=len(all_pos)
    all_rz,all_rf,all_cat,all_add=[],[],[],[]
    all_b_z,all_b_f,all_b_cat,all_b_add=[],[],[],[]
    
    for sd in range(SPLITS):
        seed=SPLIT_SEEDS[sd]
        # Per-position leak: 80/20 within each position
        np.random.seed(seed); torch.manual_seed(seed)
        tm_list, vm_list = [], []
        for pos in all_pos:
            pos_idx = np.where(sites == pos)[0]
            np.random.shuffle(pos_idx)
            n_val = max(1, int(len(pos_idx) * 0.2))
            vm_list.extend(pos_idx[:n_val])
            tm_list.extend(pos_idx[n_val:])
        tm = np.array(tm_list); vm = np.array(vm_list)
        print(f"\n  Split {sd+1}: train={len(tm)} val={len(vm)} (within-position leak)")
        print(f"    train_idx: {sorted(tm[:20].tolist())}{'...' if len(tm)>20 else ''}")
        print(f"    val_idx:   {sorted(vm[:20].tolist())}{'...' if len(vm)>20 else ''}")
        Xtz=torch.tensor(Xz[tm],dtype=torch.float32); Xvz=torch.tensor(Xz[vm],dtype=torch.float32)
        Xtf=torch.tensor(Xf[tm],dtype=torch.float32); Xvf=torch.tensor(Xf[vm],dtype=torch.float32)
        yt=torch.tensor(ys[tm],dtype=torch.float32); yv=torch.tensor(ys[vm],dtype=torch.float32)
        
        b_z,pz=_train_single(Xtz,yt,Xvz,yv,SinglePredictorV2,128,"Z_II")
        b_f,pf=_train_single(Xtf,yt,Xvf,yv,SinglePredictorV2,FAESM_DIM,"FAESM")
        # Concat is now the CrossAttnOrthoConcatFusion
        b_cat,pcat=_train_fusion(Xtz,Xtf,yt,Xvz,Xvf,yv,CrossAttnOrthoConcatFusion,"Concat",cos_penalty=0.05)
        # Add is now the GatedAddFusion
        b_add,padd=_train_fusion(Xtz,Xtf,yt,Xvz,Xvf,yv,GatedAddFusion,"Add")
        
        print(f"    -> Z_II={b_z:.4f}  FAESM={b_f:.4f}  Concat={b_cat:.4f}  Add={b_add:.4f}")
        # Save split results
        torch.save({'split':sd+1,'seed':seed,'val_idx':vm.tolist(),'train_idx':tm.tolist(),
                    'Z_II':b_z,'FAESM':b_f,'Concat':b_cat,'Add':b_add},
                   CKPT_DIR/f"split{sd+1}_seed{seed}_results.pt")
        all_rz.append(b_z); all_rf.append(b_f); all_cat.append(b_cat); all_add.append(b_add)
        yv_np=yv.numpy()
        for pred,bins_list in [(pz,all_b_z),(pf,all_b_f),(pcat,all_b_cat),(padd,all_b_add)]:
            w,a=per_bin_analysis(pred,yv_np); bins_list.append((w,a))
    
    rz=np.array(all_rz); rf=np.array(all_rf); rcat=np.array(all_cat); radd=np.array(all_add)
    print(f"\n{'='*60}\n{protein} SUMMARY")
    print(f"  {'Method':<10s} {'μ':>8s} {'σ':>8s} {'Best':>8s}")
    for name,arr in [("Z_II",rz),("FAESM",rf),("Concat",rcat),("Add",radd)]:
        print(f"  {name:<10s} {np.mean(arr):>8.4f} {np.std(arr):>8.3f} {np.max(arr):>8.4f}")
    
    print(f"\n  Per-Bin ρ (averaged):")
    print(f"  {'Bin':<12s} {'Z_II':>8s} {'FAESM':>8s} {'Concat':>8s} {'Add':>8s} {'MaxS':>8s}")
    for bi,bn in enumerate(["Lethal","Bad","Neutral","Good","Beneficial"]):
        wz=np.nanmean([b[0][bi] for b in all_b_z])
        wf=np.nanmean([b[0][bi] for b in all_b_f])
        wcat=np.nanmean([b[0][bi] for b in all_b_cat])
        wadd=np.nanmean([b[0][bi] for b in all_b_add])
        mxs=max(wz,wf)
        c="OK" if wcat>=mxs-0.02 else "  "; a="OK" if wadd>=mxs-0.02 else "  "
        print(f"  {bn:<12s} {wz:>+8.4f} {wf:>+8.4f} {wcat:>+8.4f}{c:<2s} {wadd:>+8.4f}{a:<2s} {mxs:>+8.4f}")
    
    print(f"\n  Per-Bin AUC (averaged):")
    print(f"  {'Pair':<16s} {'Z_II':>8s} {'FAESM':>8s} {'Concat':>8s} {'Add':>8s}")
    for pi,pn in enumerate(["LethalvsBad","BadvsNeutral","NeutralvsGood","GoodvsBeneficial"]):
        az=np.nanmean([b[1][pi] for b in all_b_z])
        af=np.nanmean([b[1][pi] for b in all_b_f])
        acat=np.nanmean([b[1][pi] for b in all_b_cat])
        aadd=np.nanmean([b[1][pi] for b in all_b_add])
        print(f"  {pn:<16s} {az:>8.4f} {af:>8.4f} {acat:>8.4f} {aadd:>8.4f}")

if __name__=="__main__":
    import sys
    LOG_DIR=Path("J:/conda_envs/foundry/DMS_Project/logs")
    LOG_DIR.mkdir(exist_ok=True)
    log_file=LOG_DIR/f"fusion_v3_rbd_{int(time.time())}.log"
    class Tee:
        def __init__(self,*files): self.files=files
        def write(self,obj):
            for f in self.files: f.write(obj); f.flush()
        def flush(self):
            for f in self.files: f.flush()
    f=open(log_file,'w')
    sys.stdout=Tee(sys.stdout,f)
    sys.stderr=Tee(sys.stderr,f)
    print(f"Logging to {log_file}")
    print(f"Fusion V3: CrossAttn Ortho Concat vs Gated Add"); print(f"Device: {torch.cuda.get_device_name(0)}")
    run("rbd",201,"bind_avg")
    f.close()

