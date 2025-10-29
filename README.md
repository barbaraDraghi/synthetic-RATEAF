# Synthetic RATE-AF – High-Fidelity Synthetic Data Experiments

This repository contains the complete codebase to **reproduce the experiments** described in:

> **Draghi, B. et al. (submitted)**  
> *High-fidelity synthetic data to boost under-represented patients and create virtual trial cohorts: the RATE-AF case study*  
> Submitted to *NPJ Digital Medicine – Digital Health Equity & Access* Collection.

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
│ ├── dummy_rateaf.csv # schema-only example (no real data available)
│ └── README.md # describes data fields and formats
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
No external datasets are provided. The included dummy_rateaf.csv allows you to test that the code executes successfully.


### Data Preparation

Place your dataset in the `data/` folder.
If no data are available, you can use the provided dummy_rateaf.csv (structure-only, no clinical meaning).


### Run Experiments

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
All analyses involving real RATE-AF data must be conducted under approved governance and data-sharing agreements with CPRD / MHRA.

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
| **Data** | `dummy_rateaf.csv` included for schema testing only |


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
