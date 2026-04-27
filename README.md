# TCR/BCR Primer Evaluation Pipeline 🧬

![License](https://img.shields.io/badge/License-MIT-blue.svg)
![Python](https://img.shields.io/badge/Python-3.8%2B-blue)
![R](https://img.shields.io/badge/R-4.0%2B-blue)
![openPrimeR](https://img.shields.io/badge/openPrimeR-Supported-green)

A comprehensive, automated bioinformatics pipeline for evaluating and visualizing the coverage of T-cell receptor (TCR) and B-cell receptor (BCR) multiplex PCR primers. 

This tool seamlessly integrates primer formatting, in silico PCR evaluation (via `openPrimeR`), and high-quality ggplot2 visualizations to help researchers design and optimize immune repertoire sequencing (Rep-Seq) primers.

---

## 🌟 Key Features

- **Automated Workflow**: Run the entire evaluation process—from raw Excel/CSV primer lists to final publication-ready heatmaps—with a single command.
- **High-Resolution Visualization**:
  - **Primer Binding Heatmaps**: Visualize exactly where and how each primer binds to the reference sequences, including mismatch positions.
  - **Uncovered Sequence Analysis**: Automatically extracts and plots reference sequences that failed to amplify, highlighting sequence features to guide primer redesign.
- **Multi-Target Support**: Easily switch between different immune loci (e.g., `TRB`, `IGH`, `IGK`, `IGL`).
- **Mismatch Tolerance**: Accurately simulates PCR conditions by allowing up to 3 mismatches (configurable).

---

## 🏗️ Pipeline Architecture

The pipeline consists of four integrated modules, orchestrated by `main.py`:

1. **`split_primers.py`**: Parses raw Excel/CSV primer sets and splits them into V-forward and J-reverse FASTA files.
2. **`evaluation.R`**: The core analytical engine. Uses `openPrimeR` to evaluate primer coverage against reference templates, extracting binding sites and identifying uncovered sequences.
3. **`Primer_coverage_heatmap.R`**: Generates detailed nucleotide-level alignment heatmaps for covered regions.
4. **`Uncovered_plot.R`**: Generates nucleotide-level heatmaps for sequences missed by the primer set, grouped by gene family.

---

## 🚀 Getting Started

### Prerequisites

Ensure you have the following installed on your system:

- **Python 3.8+**
  - `pandas`
  - `openpyxl`
- **R 4.0+**
  - `openPrimeR` (Bioconductor)
  - `Biostrings` (Bioconductor)
  - `ggplot2`
  - `tidyr`

#### Installing R Dependencies
```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("openPrimeR", "Biostrings"))
install.packages(c("ggplot2", "tidyr"))
```

#### Installing Python Dependencies
```bash
pip install pandas openpyxl argparse
```

### Directory Structure

Before running, ensure your project directory is structured as follows:

```text
├── data/
│   ├── IGHV.fasta              # Reference V sequences for IGH
│   ├── IGHJ.fasta              # Reference J sequences for IGH
│   ├── TRBV.fasta              # Reference V sequences for TRB
│   ├── TRBJ.fasta              # Reference J sequences for TRB
│   └── primer_set/
│       └── tcr_primers.xlsx    # Your input primer list
├── main.py                     # The main orchestrator script
├── split_primers.py            # Python formatting script
├── evaluation.R                # R evaluation script
├── Primer_coverage_heatmap.R   # R visualization script
└── Uncovered_plot.R            # R visualization script
```

*Note: Reference FASTA headers must follow the IMGT format: `ACCESSION|GROUP|SPECIES|FUNCTION`.*

---

## 💻 Usage

Run the complete pipeline using the `main.py` wrapper. 

### Basic Command

```bash
python main.py -i data/primer_set/tcr_primers.xlsx -t TRB
```

### Arguments

| Argument | Short | Description | Default |
| :--- | :---: | :--- | :--- |
| `--input` | `-i` | Path to the input Excel or CSV file containing primers. | **Required** |
| `--target` | `-t` | Target locus name (e.g., `TRB`, `IGH`, `IGK`, `IGL`). The script will automatically look for `{TARGET}V` and `{TARGET}J` references. | **Required** |
| `--v_pattern`| | String pattern used to identify V region primers in the input file. | `V` |
| `--j_pattern`| | String pattern used to identify J region primers in the input file. | `J` |

### Input File Format

The input Excel (`.xlsx`) or `.csv` file must contain at least two columns:
1. **Primer ID**: (e.g., `TRBV1`, `IGHJ4`)
2. **Sequence**: (e.g., `GCTACTTCGGAGCCTCGG`)

---

## 📊 Output Directory Structure

Upon successful execution, the pipeline generates a `results/` folder containing comprehensive reports and plots:

```text
results/
├── coverage_tables/
│   ├── TRBV_binding_sites.csv       # Exact binding coordinates & mismatches
│   └── TRBV_cov.csv                 # Overall coverage statistics
├── coverage_plots/
│   └── TRBV_cov.png                 # Summary coverage plot
├── primer_reference_heatmaps/
│   └── TRBV/
│       ├── TRBV1_covered_region.png # High-res nucleotide alignment heatmap
│       └── ...
└── uncovered_templates/
    └── TRBV_unc/
        ├── TRBV_family1.csv         # Uncovered sequences by family
        └── TRBV_family1_heatmap.png # Nucleotide heatmap for redesign
```

---

## 💡 Use Cases for Primer Optimization

1. **Identify Blind Spots**: Review the `uncovered_templates` heatmaps to see exactly which conserved regions are missing coverage.
2. **Optimize Mismatches**: Use the `primer_reference_heatmaps` to visualize where mismatches occur. If a mismatch is consistently at the 3' end, the primer must be redesigned.
3. **Patent & IP Documentation**: The outputs generated by this pipeline (binding site CSVs, comprehensive heatmaps) serve as excellent technical disclosure materials for patent applications and software copyright registrations.

---

## 📄 License

This project is licensed under the MIT License - see the LICENSE file for details.

---
*Developed for advanced Immune Repertoire Sequencing (Rep-Seq) primer design and evaluation.*
