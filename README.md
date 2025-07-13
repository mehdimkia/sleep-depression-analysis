# Sleep Duration & Fragmentation Predict Incident Depressive Symptoms
**R survival analysis · Maastricht Study · MSc thesis project**

## TL;DR
- Objective: Does objectively measured sleep predict 4-year depression risk?
- Methods: Cox models, interaction terms, spline regression, sensitivity analyses.
- Key result: High sleep fragmentation (Q4 vs Q1) ↑ risk by 54 %.

## Reproduce
1. Clone repo
2. `Rscript scripts/01_descriptives.R` … etc.
3. Figures land in `Outputs/`.

## Project structure
```txt
sleep-depression-analysis/
├── LICENSE
├── Outputs/
│   ├── SplineModeling/
│   ├── CoxModels/
│   ├── InteractionTests/
│   └── Descriptives/
└── Scripts/
    ├── SplineModeling/
    ├── InteractionTests/
    ├── Descriptives/
    ├── CoxModels/
    └── SensitivityAnalysis/

