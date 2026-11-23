# High-fidelity synthetic data to boost under-represented patients and create virtual trial cohorts: the RATE-AF case study

This repository contains the complete codebase to **reproduce the experiments** described in:

> Draghi, B., Attal, D., Kotecha, D. et al. Synthetic data to boost under-represented patients and create virtual trial cohorts: the RATE-AF case study. Sci Rep 16, 36678 (2026). > [https://doi.org/10.1038/s41598-026-60252-z](https://doi.org/10.1038/s41598-026-60252-z)

The scripts here implement the **controlled**, **exploratory**, and **replication** boosting experiments, building on the original **RATE-AF randomized clinical trial**  
> Kotecha, D. et al. *Effect of Digoxin vs Bisoprolol for Heart Rate Control in Atrial Fibrillation on Patient-Reported Quality of Life.*  
> *JAMA* 324(24):2497–2508. [https://doi.org/10.1001/jama.2020.23138](https://doi.org/10.1001/jama.2020.23138)

---

##  Background

The **RATE-AF trial** compared *low-dose digoxin* vs *bisoprolol* for heart-rate control in older patients with permanent atrial fibrillation and heart-failure symptoms.  
This repository provides a reproducible synthetic-data framework that:

1. Reconstructs the baseline RATE-AF population structure.  
2. Generates **synthetic boosted datasets** to recover under-represented patterns.  
3. Evaluates **utility, fidelity, and bias** across all experimental conditions.  
4. Demonstrates how **synthetic data can enhance trial inclusivity** and support *virtual trial replication*.

---

## 📂 Repository Structure

```
synthetic-RATEAF/
│
├── data/
│ ├── RATAF-dummy.csv # schema-only example (no real data available)
│
├── main/
│ ├── controlled_boosting_experiment.R # main script for controlled boosting experiment
│ ├── exploratory_boosting_experiment.R # main script for exploratory boosting experiment
│ ├── replication_experiment.R # main script for full synthetic replication experiment
│
├── utils/
│ ├── functions.R # helper functions for preprocessing & variable mapping
│ ├── synthetic_data_generation.R # conditional data generation and boosting logic
│ ├── utility_evaluation.R # model-based utility and performance metrics
│ ├── fidelity_evaluation.R # log-probability and distributional fidelity measures
│ ├── bias_analysis.R # subgroup fairness and bias detection utilities
│ └── evaluation/ # internal scripts for aggregating and plotting results
│
├── .gitignore
├── .Rproj.user/
├── .Rhistory
└── README.md # you are here

```

---

## ⚙️ How to Run

### Requirements

```
- R version: ≥ 4.2  
- Tested environment: R 4.3.1 on Windows 10  
```

All scripts use base R and tidyverse functionality only.
No external datasets are provided. The included RATEAF-dummy.csv allows you to test that the code executes successfully.

### Step 1 – Download the Code
Download or clone this repository and open it in **RStudio**.

### Step 2 – Working Directory Handling
Each main script automatically detects its own location when executed in **RStudio** or via **Rscript**, and sets the correct working directory.  
No manual path edits are required under normal conditions.

However, if automatic detection fails (e.g., running from an unusual environment),  
you can manually set your working directory to the repository’s `/main` folder using your full local path:

```r
# Example (Windows)
setwd("H:/Projects/synthetic-RATEAF/main")

# Example (Mac/Linux)
setwd("~/Projects/synthetic-RATEAF/main")

```
💡 Tip:
Each script (controlled_boosting_experiment.R, exploratory_boosting_experiment.R, replication_experiment.R)
includes an internal routine that automatically sets the working directory and sources the required utilities from ../utils/.
You can disable or adjust this behavior if integrating into another workflow.

### Step 3 – Data Preparation

Place your dataset in the `data/` folder.  

- If you have **approved access to the real RATE-AF dataset** (e.g., `RATE-AF.dta`),  
  simply place the file in this folder — the scripts are already configured to recognize and process it automatically.  

- If you **do not have access** to the real data, you can use the provided [`RATEAF-dummy.csv`](data/RATEAF-dummy.csv) —  
  a **two-row illustrative file** that serves as a reference for variable names and expected structure to create your dummy data or to be able to adapt the code to your data.
  (*no real clinical information included*).  

Please see the **Data Policy** section for further clarification on data governance and access.


### Step 4 – Run Experiments

Each script in `/main` corresponds to one of the analyses presented in the NPJ manuscript.

```
# Controlled boosting (targeted subgroup augmentation)
source("main/controlled_boosting_experiment.R")

# Exploratory boosting (coverage-based augmentation)
source("main/exploratory_boosting_experiment.R")

# Replication (full virtual cohort reconstruction)
source("main/replication_experiment.R")
```
### Outputs

Each experiment generates a combination of **summary tables** and **diagnostic plots** under `/results/`.  
The main output categories include:

| Output Type | Description |
|--------------|-------------|
| **Fidelity metrics** | Quantifies similarity between real and synthetic data (e.g., marginal and conditional log-density, correlation structure). |
| **Utility metrics** | Evaluates predictive consistency and downstream model performance between real and synthetic datasets. |
| **Bias analysis** | Examines subgroup representation, coverage, and fairness before and after boosting. |
| **Plots** | Includes density overlays, conditional coverage maps, and improvement trajectories. |

Together, these outputs reproduce the figures and numerical results discussed in the NPJ manuscript.

---

## 🔒 Data Policy

This repository contains **only code and synthetic examples**.  
No identifiable patient data are included or distributed.  
All analyses involving real RATE-AF data must be conducted under approved governance and data-sharing agreements with **CPRD / MHRA**.  
Researchers wishing to access the original RATE-AF clinical trial data should contact the Chief Investigator **Prof Dipak Kotecha (University of Birmingham)**.

### Example Dataset – `RATEAF-dummy.csv`

The file [`data/RATEAF-dummy.csv`](data/RATEAF-dummy.csv) is a **two-row illustrative dataset** provided as a **guideline** for expected structure and variable formatting.  
It contains *no real clinical information* and is intended solely to demonstrate correct metadata and column definitions.

Researchers interested in testing or extending the framework can:

1. **Generate a dummy dataset** that follows the same metadata schema defined in this file, **or**  
2. **Adapt the framework** to work with alternative datasets, following the same logic and preprocessing steps applied here, using the CSV structure as reference.

---

## Citation

If you use this repository or build upon this codebase, please cite:
> Draghi, B. et al. *High-fidelity synthetic data to boost under-represented patients and create virtual trial cohorts: the RATE-AF case study.*  
> Submitted to *NPJ Digital Medicine – Digital Health Equity & Access*, 2025.


---

## Reproducibility Checklist

| Item | Details |
|------|----------|
| **R version** | 4.3.1 (Windows 10) |
| **Random seed** | `set.seed(123)` used consistently across scripts |
| **Execution time** | ≈ 5–10 minutes per experiment on a standard laptop (Intel i7, 16 GB RAM) |
| **Output directory** | `/results/` (auto-generated) |
| **Dependencies** | tidyverse, purrr, readr, ggplot2, dplyr |
| **Data** | `RATEAF-dummy.csv` included for schema testing only |


### Adaptability and Reuse

While some preprocessing and postprocessing functions in this repository are **specific to the RATE-AF dataset**  
(e.g., variable mappings, NYHA/EHRA conversions, and clinical endpoint formatting),  
the overall pipeline has been designed for **general applicability**.

The core modules within `/utils/` (particularly `synthetic_data_generation.R`, `utility_evaluation.R`, and `fidelity_evaluation.R`)  
can be easily adapted to handle **other clinical or tabular datasets** with minimal modification.

Researchers interested in:
- extending the methodology to a different trial or registry dataset,  
- adapting the boosting logic for new variable structures, or  
- integrating the evaluation metrics into their own pipeline  

are warmly encouraged to get in touch with the author.


---

## 📧 Contact

**Barbara Draghi**  
Clinical Practice Research Datalink (MHRA) / Brunel University London

📧 barbara.draghi@mhra.gov.uk  / barbara.draghi@brunel.ac.uk
🌐 [https://barbaradraghi.github.io](https://barbaradraghi.github.io)


---

*This repository is shared for transparency and educational reuse.  
It does not contain or distribute any real patient-level data.*
