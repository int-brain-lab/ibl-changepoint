# Change-point detection and ideal-observer analyses

This repository contains a toolbox for model fitting and plotting of ideal observer models for the IBL task.
For the moment, the code is only available in Matlab.

### Installation guide

1. [Download](https://github.com/int-brain-lab/ibl-changepoint/archive/master.zip) or clone the repository on your local computer.
2. Add the base folder of the repository (which contains the `ibl_change_add2path.m` file) to your Matlab path. 
   **Warning:** It is not recommended to permanently add the entire repository tree to your Matlab path, since this might cause function name clashes with other projects.
3. If you do not already have them, install *Bayesian Adaptive Direct Search* (BADS; an optimization toolbox) from [here](https://github.com/lacerbi/bads) and *Variational Bayesian Monte Carlo* (VBMC; an approximate posterior inference toolbox) from [here](https://github.com/lacerbi/vbmc).
4. Download the IBL data for the exemplar mice to CSV files in the `data` folder using [this Python script](https://github.com/int-brain-lab/ibl-changepoint/blob/master/matlab/data/fetch_data.ipynb) (you will need IBL access credentials).

### Basic overview

As a simple example, we are going to fit an "omniscient" Bayesian observer to a mouse stable-sessions data, assuming a biased lapse (see below for an explanation). The exemplar mice are `CSHL_005`,`CSHL_007`,`IBL-T1`,`IBL-T4`,`ibl_witten_04`,`ibl_witten_05`; we are going here to fit the first one:

```
changepoint_add2path;           % Add the changepoint analysis toolbox to the Matlab path for this session
model_name = 'omniscient_biasedlapse';
mouse_name = 'CSHL_005';
data = read_data_from_csv(mouse_name);
Nopts = 1;                      % Perform only one optimization (but you should use multiple starting points)
params = fit_model(model_name,data,Nopts);
plot_fit(params,data);
```
