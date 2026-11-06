# PLOS Computational Biology 初投目录

目标期刊：PLOS Computational Biology  
类型：Research Article  
Section：Biological Macromolecules

本目录包含完整初投文件及投稿辅助材料，但**尚未在 Editorial Manager 正式提交**。

## 上传文件

| 文件 | 投稿系统中的用途 |
| --- | --- |
| `Manuscript.pdf` | 主文；14 页，含正文、2 张主图、参考文献和 SI captions |
| `Cover_Letter.pdf` | 独立 Cover Letter；1 页 |
| `Fig1.tif`、`Fig2.tif` | 主图；600 dpi、RGB、LZW 压缩 TIFF |
| `S1_Text.pdf` | Supporting Text；9 页扩展方法与质量控制 |
| `S1_Fig.tif`、`S2_Fig.tif` | Supporting Figures；600 dpi、RGB、LZW 压缩 TIFF |
| `S1_Table.xlsx` | DMS 去偏、block-CV、轻量 probe；13 个 worksheets |
| `S2_Table.xlsx` | 13-entry distogram panel、5 个模板 steering targets、GA88/GB88 E49b |
| `S3_Table.xlsx` | 30 个坐标采样 cell 汇总 + 全部 480 条结构级 RMSD |
| `S4_Table.xlsx` | 融合、梯度干预、D1、16 次 E429 复测、CASP3 界面逐突变体数据 |

SI  captions 已置于主文参考文献之后。Excel 附件保留原始浮点数值，不以正文四舍五入后的表格作为数据源。

## 不上传的辅助材料

- `Submission_Metadata.md`：标题、摘要、声明和投稿字段，供复制到系统。
- `source/`：主文、S1 Text、Cover Letter 的 LaTeX 源文件及矢量图，留作返修。
- `SHA256SUMS.txt`：目录内文件完整性哈希。
- `validation.json`：自动核验结果。

初投不必上传 `source/`。PLOS 允许初投使用含正文和主图的单个 PDF；本目录同时提供独立主图和 SI，便于完整提交及后续返修。

## 提交前仍需人工完成

1. 作者审阅并批准正文、Cover Letter 和所有声明。
2. 在系统中选定合适的 Academic Editor，并填写至少 4 名无利益冲突的潜在审稿人。
3. 确认是否申请 Publication Fee Assistance；只能在首次投稿时申请。
4. 按 `Submission_Metadata.md` 填写投稿字段并确认没有相关稿件在其他期刊审理。
5. 上传后逐页检查系统生成的合并 PDF，再执行最终 Submit。

## 官方要求

- 投稿入口：https://www.editorialmanager.com/pcompbiol
- 稿件要求：https://journals.plos.org/ploscompbiol/s/submission-guidelines
- Supporting Information：https://journals.plos.org/ploscompbiol/s/supporting-information
- 主图要求：https://journals.plos.org/ploscompbiol/s/figures

## 可复现构建

在仓库根目录：

```bash
python scripts/make_figures.py
python scripts/make_supporting_figures.py
python scripts/make_supporting_tables.py
python scripts/validate_repository.py
python scripts/assemble_submission.py
python scripts/validate_submission.py
```

LaTeX 通过项目的 LaTeX 编译工具完成；不要用另一份 PDF 覆盖 `Manuscript.pdf` 而不重新执行校验。
