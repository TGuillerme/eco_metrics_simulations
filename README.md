# The what, how and why of trait-based analyses in ecology

Authors: [Thomas Guillerme](https://github/TGuillerme), [Pedro Cardoso](https://github.com/cardosopmb), [Maria Wagner Jørgensen](https://github.com/MariaWagnerJ), [Stefano Mammola](https://github.com/StefanoMammola), [Thomas J. Matthews](https://github.com/txm676).

This repository contains all the code and data used in [this paper](preprint).

[![DOI](zenodo)](zenodo)

<!-- ## Supplementary material -->

 <!-- * The supplementary vignette is available [here]() (or [here in html]()) -->

## Data

The empirical data is available [here](https://github.com/TGuillerme/eco_metrics_simulations/tree/master/Data/Raw).
The processed data is available [here](https://github.com/TGuillerme/eco_metrics_simulations/tree/master/Data/Processed).

## Reproducing the whole analysis

To reproduce the entire paper, you can follow these steps:

 1. **Data simulation visualisation (figure 1)**: the script `01_Simulating_spaces_and_stressors.Rmd` guides you through the steps to simulate the data for the paper. It also allows to reproduce figure 1.
 2. **Simulated data pipeline**: the script `02_Simulation_pipeline.Rmd` allows to generate the individual scripts to reproduce the simulations for each stressor and different numbers of dimensions. This generates a function factory that in turns is rand through the replicating function in order to get an optimal number of replicates. Note that reproducing a single run takes several minutes while reproducing all of them takes several hours. This script generates the files `sim_<STRESSOR_NAME>_<DIMENSIONS>d.rda`.
 3. **Empirical data pipeline**: the script `03_Empirical_pipeline.Rmd` calculates the metrics values the same way as the `02_Simulation_pipeline.Rmd` one. It generates the file `empirical_data.rda`.
 4. **Summarising the results (figure 2 and 3)**: the script `04_Summarise_results.Rmd` allows to read in the different results, scale them, analyse them and summarise them. It can be used to reproduce figure 2 and 3 as well as the supplementary figures and tables.


<!-- ## Citing this work -->

<!-- To cite the paper, please use: -->
 <!-- * Preprint -->

<!-- To cite this repository, please use:
 * Zenodo
 -->
