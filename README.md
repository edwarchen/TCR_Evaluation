# TCR_pro 项目说明

本项目主要用于评估引物对 V/J 参考模板的覆盖率，并输出覆盖统计、覆盖图与引物命中位点信息。

## 目录说明

- `data/`
  - 存放参考序列（FASTA）。
  - 当前包含：`TRBV.fasta`、`TRBJ.fasta`、`IGHV.fasta`、`IGHJ.fasta`、`IGKV.fasta`、`IGKJ.fasta`、`IGLV.fasta`、`IGLJ.fasta`。
  - 所有参考序列已经去除序列中的`.`等gap

- `data/primer_set/`
  - 存放引物文件（Excel 和 FASTA）。将Excel中的引物序列转换为了FASTA格式；
  - 主要输入文件在这里，例如：`TRBV_primers.fasta`、`TRBJ_primers.fasta`、`IGHV_primers.fasta`、`IGHJ_primers.fasta`。

- `scripts/`
  - 存放评估与辅助脚本（R/Python）。
  
- `results/`
  - 存放运行输出结果。
  - `coverage_tables/`：覆盖统计表、位点表。
  - `coverage_plots/`：覆盖图（png）。
  - `primer_reference_heatmaps/`：引物与参考序列对齐热图。
  - `uncovered_templates/`：未覆盖模板表及其分组可视化图。

## 脚本作用

- `scripts/TRB_eval.R`
  - 评估 `TRBV/TRBJ` 引物覆盖率。
  - 输出 `TRBV_cov.csv`、`TRBJ_cov.csv` 和对应 `*_cov.png`。

- `scripts/*_eval.R`
  - 评估 `IGHV/IGHJ` 引物覆盖率。
  - 支持 J 引物裁剪参数（`trim`）。
  - 额外输出引物在模板上的匹配位点：`IGHV_binding_sites.csv`、`IGHJ_binding_sites.csv`。
  - 脚本需要自行指定允许错配数（`max_mismatch`）。
  - 能够输出未覆盖到的引物，用于后续的补充引物设计。

- `scripts/Primer_coverage_heatmap.R`
  - 根据 `*_binding_sites.csv` 和对应参考 FASTA 绘制“参考序列 vs 引物序列”对齐热图。
  - 图中上方是完整参考序列，下方是按命中 `start/end` 对齐后的引物序列，未覆盖位置用 `-` 补齐。
  - 脚本中需要自行指定基因家族名称（如 `IGHV`、`IGHJ`）。

- `scripts/Uncovered_plot.R`
  - 将未覆盖模板的 CSV 画成碱基热图。
  - 当前支持传入一个目录，批量处理该目录下所有 CSV。

- `scripts/split_primers.py`
  - 从引物表（Excel/CSV）拆分出 V/J 引物 FASTA。
  - 自动给 V 引物加 `_fw`，J 引物加 `_rev`（按脚本当前逻辑）。

- `scripts/step_by_step_eval.R`
  - 精简版示例脚本，用于快速验证覆盖率流程（便于调试）。

- `scripts/official_pipeline.R`
  - openPrimeR 官方示例流程（参考用），展示模板读取、约束检查和覆盖率统计的标准调用方式。

## 结果文件说明

- 覆盖统计：`results/coverage_tables/*_cov.csv`
- 覆盖图：`results/coverage_plots/*_cov.png`
  - 显示各模板/分组的覆盖情况。
- 位点结果（部分脚本）：`results/coverage_tables/*_binding_sites.csv`
  - 记录引物命中模板的具体区间（`start/end`）。
- 引物-参考对齐图：`results/primer_reference_heatmaps/**.png`
  - 每张图显示一个 `template_id + primer_id` 的两行对齐关系。
- 未覆盖模板：`results/uncovered_templates/*_uncovered.csv`
  - 历史分析产物，可继续按家族拆分或作图。

## 最小使用方式

在项目根目录运行：

```bash
Rscript scripts/TRB_eval.R
Rscript scripts/IGH_eval.R
Rscript scripts/IGK_eval.R
Rscript scripts/IGL_eval.R
```

若要画未覆盖模板热图，先在 `scripts/Uncovered_plot.R` 中设置目录，再运行：

```bash
Rscript scripts/Uncovered_plot.R
```

## 备注

- 各评估脚本里 `cfg` 段定义了输入参考与引物路径，是你更换数据时最常改的部分。
- 若运行报错，先检查：
  - 文件路径是否存在；
  - 引物 ID 后缀是否与方向一致（`_fw` / `_rev`）；
  - 参考序列是否为预期清洗版本（是否需要 `*_clean.fasta`）。
