# Merge Chains executable

This executable takes multiple .root chains and merges them into one. There are options for thinning, burn-in, excluding branches or including branches into the output merged chain.

The default output is `merged_chain.root`, can be overwritten with `--output` argument. No burn-in or thinning applied by default.

For the full list of available arguments, run
```bash
./merge_chains --help
```


**It is best to run MCMC diagnostics before merging the chains**: the diagnostics will tell you how much burn-in to remove, and whether it is appropriate to apply any thinning or not (e.g. if the autocorrelations are long, thinning won't have much impact on the posterior).

## Default chain merging

To merge chains without any burn-in removed (e.g. appropriate for NUTS MCMC):

```bash
./merge_chains /location/of/your/chains
```

## Adding and removing branches

If running the above, cowards the end of the cout you'll see two printout lines, one with the all available branches/parameters listed, and one line with final list of branches to save in the output merged file.

In my example the these lines are:

```
All available branches (28): ['sin2th_12', 'sin2th_23', 'sin2th_13', 'delm2_12', 'delm2_23', 'delta_cp', 'Eps_ee', 'Eps_emu', 'Eps_etau', 'Eps_mumu', 'Eps_mutau', 'Eps_tautau', 'Delta_emu', 'Delta_etau', 'Delta_mutau', 'ElecCoup', 'UpCoup', 'DownCoup', 'ProdH', 'xsec_0', 'LogL', 'accProb', 'step', 'stepTime', 'LogL_sample_0', 'LogL_sample_1', 'LogL_systematic_osc_cov', 'LogL_systematic_xsec_cov']
Final list of branches to keep (7): ['sin2th_12', 'sin2th_23', 'sin2th_13', 'delm2_12', 'delm2_23', 'delta_cp', 'LogL']
```

Lets say I want to get rid of `'sin2th_12', 'delm2_12'` from the saved output, and add `'step', 'ProdH'` from the available branches. Simply run:

```bash
./merge_chains /location/of/your/chains --ignore-branches sin2th_12 delm2_12 --keep-branches step ProdH
```

## Common usage

Most commonly you will want to apply burn-in, perhaps some thinning, maybe keep a few branches for diagnostics and ignore others. Example:

```bash
./merge_chains /folder/to/your/chains --burn-in 100000 --thin 10 --output merged_mcmc_stan_ana2024prod5.1_realdatafit.root --keep-branches Calibration RelativeCalib --ignore-branches logprob step
```
