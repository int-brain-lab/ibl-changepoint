# Change-point detection and ideal-observer analyses

This repository contains a toolbox for model fitting and plotting of ideal observer models for the IBL task.
For the moment, the code is only available in Matlab.

### Installation guide

1. [Download](https://github.com/int-brain-lab/ibl-changepoint/archive/master.zip) or clone the repository on your local computer.
2. Add the base folder of the repository (which contains the `ibl_changepoint_add2path.m` file) to your Matlab path. 
   **Warning:** It is not recommended to permanently add the entire repository tree to your Matlab path, since this might cause function name clashes with other projects.
3. If you do not already have them, install *Bayesian Adaptive Direct Search* (BADS; an optimization toolbox) from [here](https://github.com/lacerbi/bads) and *Variational Bayesian Monte Carlo* (VBMC; an approximate posterior inference toolbox) from [here](https://github.com/lacerbi/vbmc).
4. Download the IBL data for the exemplar mice to CSV files in the `data` folder using [this Python script](https://github.com/int-brain-lab/ibl-changepoint/blob/master/data/fetch_data.ipynb) (you will need IBL access credentials).

### Basic overview

As a simple example, we are going to fit an "omniscient" Bayesian observer to a mouse stable-sessions data, assuming a biased lapse (see below for an explanation). The exemplar mice are `CSHL_005`,`CSHL_007`,`IBL-T1`,`IBL-T4`,`ibl_witten_04`,`ibl_witten_05`; we are going here to fit the first one:

```
ibl_changepoint_add2path;           % Add the changepoint analysis toolbox to the Matlab path for this session
model_name = 'omniscient_contrastnoise_biasedlapse';
mouse_name = 'CSHL_005';
data = read_data_from_csv(mouse_name);
Nopts = 1;                      % Perform only one optimization (but you should use multiple starting points)
params = fit_model(model_name,data,Nopts);
plot_fit(params,data);
```

The available models are:
- `psychofun`: simple psychometric function fit (separately for each probability level).
- `omniscient`: "omniscient" Bayesian observer that knows the true stimulus probability in each block.
- `changepoint`: change-point Bayesian observer that tracks the stimulus probability from trial to trial.

Specific features can be added to the base observer models, separated by subscripts:
- `_lapse`: adds a probability `lapse_rate` of a random response with equal probability.
- `_biasedlapse`: as above, but lapses have a `lapse_bias` probability of responding "Left" (`lapse_bias = 0.5` is unbiased lapse).
- `_softmax`: adds a softmax probabilistic step to the decision rule.
- `_runlength`: adds flexibility to run-length related parameters (`changepoint` observer only).
- `_freesym`: adds flexibility to the a-priori block probabilities for biased blocks, assuming symmetry (`changepoint` observer only).

Specific noise models are added as:
- `contrastnoise`: a simple noise model based on a noisy measurement of the contrast level.
- `nakarusthon`: a noise model inspired by the [Naka-Rushton model](https://www.jneurosci.org/content/17/21/8621) of contrast perception.

Multiple datasets and models can be fitted in batch by using the `batch_model_fit.m` function.

If a dataset name is suffixed with `_unbiased`, only unbiased (50/50) blocks are loaded.

### Predicting optimal bias shifts

The `fit_all_unbiased.m` script contains the entire pipeline to fit psychometric functions and several "omniscient" models to the unbiased data only, in order to make predictions about the mice behavior in the full sessions that include change-point blocks.
The `scripts` folder contains bash scripts to launch several jobs on a computer cluster (which uses the Slurm job manager).

Note that in addition to obtaining maximum-likelihood estimates of the parameters (via [BADS](https://github.com/lacerbi/bads)), we also compute (approximate) Bayesian posteriors using [VBMC](https://github.com/lacerbi/vbmc).

### Troubleshooting

The IBL change-point toolbox is work in progress, and I plan to extend the documentation where needed.
In the meanwhile, for any question please drop me an email at <luigi.acerbi@unige.ch>.

### License

The IBL change-point toolbox is released under the terms of the [MIT license](https://github.com/int-brain-lab/ibl-changepoint/blob/master/LICENSE).

