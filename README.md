# Data and code for Lai et al. (2026) Oikos

Welcome! This is a repository for the data and models in our paper:

> Lai HR, Chong KY, Neo L, Wardle DA, Vesk PA (2026). 
> Asymmetric transferability of trait–environment relationships between native
> and exotic plant species in Singaporean forests differing in land-use history.
> *Oikos*

## Usage

The main scripts are in the `R` folder:

- `1_prep_inputs.R` wrangles the data for analyses (see also Notes below)
- `2a_greta.R` fits the statistical models
- `3_posteriors.R` examines the model results, plots some figures etc.
- `cross_val_greta.R` calculates the prediction accuracies from cross validation
- `KL_div.R` calculates the Kullback--Leibler divergences
- `horizon.R` computes and plots the prediction horizons
- Auxillary analyses:
  - `ess.R` calculates the effective sample sizes of the models
  - `geweke.R` computes the Geweke diagnostics to assess model convergence
  - `plooic.R` computes the parameter-to-data ratio

## Notes

- The filenames have a mixed naming conventions due to legacy, they don't mean
  much in terms of ordering
- Over the course of working this into a manuscript, we have packaged the data
  into a proper `R` package [`novelforestSG`](https://github.com/hrlai/novelforestSG).
  But our codes have not fully transitioned, some parts still reads data
  directly from the `Data` folder. Please consider using the `novelforestSG`
  package in your future work.

## Contact

Please feel free to post your issue in the
[Issue page](https://github.com/hrlai/Lai_et_al_2026_Oikos/issues).
