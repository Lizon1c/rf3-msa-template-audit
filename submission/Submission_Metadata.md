# PLOS Computational Biology — Submission Metadata

本文件用于复制到投稿系统，不是需要上传的稿件附件。英文声明在提交前须由作者再次确认。

## Article

- **Article type:** Research Article
- **Section:** Biological Macromolecules
- **Full title:** An input-channel audit of RoseTTAFold3: MSA writing, mutation-effect prediction, and conformational-state selection
- **Short title:** Input-channel audit of RoseTTAFold3
- **Keywords:** RoseTTAFold3; AlphaFold3; multiple sequence alignment; structural templates; deep mutational scanning; conformational selection; pair representation; protein structure prediction
- **Submission portal:** https://www.editorialmanager.com/pcompbiol

## Author and corresponding author

- **Name:** Cheng Ding
- **Affiliation:** School of Life Science and Technology, ShanghaiTech University, Shanghai, China
- **Email:** dingcheng2024@shanghaitech.edu.cn
- **ORCID:** 0009-0001-7221-9066
- **Author Contributions:** Conceptualization, Data curation, Formal analysis, Investigation, Methodology, Software, Validation, Visualization, Writing – original draft, Writing – review & editing

## Abstract

AlphaFold3-style structure predictors discard the explicit multiple sequence alignment (MSA) representation after a compact module within each recycle, leaving the pair representation to carry evolutionary information. Here we audit this input channel in RoseTTAFold3 across representation geometry, deep mutational scanning (DMS), and conformational-state selection. Across three proteins, MSA input reduced effective rank by 16–42% while preserving mutant-conditioned directions and improved cross-position DMS prediction. De-biased gains were +0.14–0.15 Spearman ρ on RBD and +0.14–0.18 on CASP3; contiguous-block cross-validation yielded gains of +0.256 and +0.264, respectively. A single MSA recycle installed approximately 75% of the geometric and 85% of the functional effect, with only about 25% of the geometric effect erased by nine subsequent MSA-free recycles. On RBD, conservation contributed approximately two-thirds and covariation one-third of the alignment-content gain, which approached saturation near eight rows. Distogram analysis of 13 two-state entries and explicit sampling of 480 structures showed that MSA effects on conformational selection were target-dependent and consistent with family consensus. Structural templates constituted a distinct prior: they impaired RBD mutation-effect prediction and strongly steered sampled conformations, but no universal template–MSA hierarchy was observed. Native MSAs preserved or restored the native fold in GA88 and produced a mixed ensemble in GB88 when opposed by an alternative-fold template. Several template-conditioned distogram preferences disagreed with sampled coordinates. These results establish a multi-level account of MSA and template priors in RoseTTAFold3 and show that conformational audits should include explicit coordinate generation.

## Author Summary

Proteins can adopt different shapes to perform their functions, and computational models are increasingly used to predict these shapes and the effects of mutations. These models can use information from many related protein sequences, but it has been unclear how this information is retained inside newer model architectures. We investigated RoseTTAFold3 and found that evolutionary information leaves a lasting imprint after a single processing step. This imprint improves the prediction of mutation effects and can influence which of several possible protein shapes the model selects. Supplying an example structure, known as a structural template, affects the model in a different way: it can strongly guide the predicted shape but reduce sensitivity to mutations. The outcome of combining evolutionary information and structural templates depends on the protein. We also found that an internal distance-based score can disagree with the three-dimensional structures that the model actually generates. Our findings provide a unified way to examine how different inputs guide structure prediction and highlight the importance of checking generated structures directly.

## Financial Disclosure Statement

The author(s) received no specific funding for this work.

该句采用 PLOS 对无经费研究的建议表述。若投稿前获得已承诺的经费支持，必须同步修改此处、Cover Letter 和 PFA 信息。

## Competing Interests

The author declares that no competing interests exist.

## Data Availability Statement

Analysis code and the numerical outputs underlying all reported summary statistics are available at https://github.com/Lizon1c/rf3-msa-template-audit. The 480 sampled coordinate files and analysis-ready feature archives are available at https://doi.org/10.5281/zenodo.23859152. Publicly available structural and DMS input datasets are identified in Materials and Methods.

## Generative AI disclosure

Generative AI tools were used to assist in unifying and editing the manuscript text from the author's experimental archive and fact pack; the author verified all numerical values, citations, and scientific statements and takes full responsibility for the content of the manuscript.

该披露已在 Acknowledgments 和 Cover Letter 中体现；若投稿系统另有生成式 AI 问题，按实际情况一致填写。

## Related manuscripts

No related or duplicate manuscript is currently under consideration or accepted for publication elsewhere.

提交前请作者确认该声明仍成立。若另有相关稿件进入投稿状态，需要向期刊申报并上传相关稿件。

## Editor and reviewer suggestions — author action required

- **Academic Editor:** 待作者选定；应熟悉蛋白质结构预测、构象系综或计算结构生物学。
- **Potential reviewers:** 至少 4 人，逐人记录姓名、单位、邮箱、专业方向及利益冲突核查结果。
- 不得选择近期合作者、同单位人员或存在其他利益冲突的人员。
- **Opposed reviewers:** 无默认名单；如确有需要，由作者另附理由。

## Publication Fee Assistance — author decision required at initial submission

PFA 只能在首次投稿时申请。当前 Financial Disclosure 为无专项经费，但是否申请 PFA 仍须根据真实支付能力决定，并如实说明已调查的单位、图书馆、政府、基金会、导师或其他科研经费渠道。不要把尚未承诺的科创基金或导师支持写成已有经费，也不要在 Cover Letter 中写 PFA 申请理由。

## Final author checks

- 核对作者名单、单位、邮箱、ORCID 和所有声明。
- 核对 Cover Letter 的原创性、未一稿多投、无利益冲突和无经费表述。
- 在系统中填写 Abstract、Financial Disclosure、Competing Interests 和完整 Data Availability Statement。
- 上传后检查系统生成的合并 PDF：主图、SI captions、链接、页序和特殊字符必须完整。
