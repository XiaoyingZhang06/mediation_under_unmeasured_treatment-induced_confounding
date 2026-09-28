# Identification and Semiparametric Inference of Interventional Effects under Unmeasured Treatment-Induced Confounding

This repository contains the R code for the simulation studies and the real-data application in the paper.

## Structure

```
application/
  application.R        # Application: effect of discrimination on well-being among older adults (HRS 2020)
  wellbeing_data.rds   # Preprocessed analysis data (n = 4090)
simulation/
  biM_para.R           # Binary mediator, parametric models
  conM_para.R          # Continuous mediator, parametric models
  biM_ML.R             # Binary mediator, nonlinear setting with kernel methods and cross-fitting
  conM_ML.R            # Continuous mediator, nonlinear setting with kernel methods and cross-fitting
```

## Usage

Run each script from its own directory, for example:

```bash
cd simulation && Rscript biM_para.R
cd application && Rscript application.R
```

- Parametric simulations use n = 1000 with 1000 replications and 200 bootstrap samples.
- ML simulations use n = 250, 500, and 1000, with 1000 replications each.

## Data

`wellbeing_data.rds` is derived from the RAND HRS 1992–2020 file and the HRS 2020 Leave-Behind questionnaire. Variables were selected and missing values were imputed with random forests. The raw-data processing code is kept as comments at the top of `application.R`. The raw data must be obtained from HRS.
