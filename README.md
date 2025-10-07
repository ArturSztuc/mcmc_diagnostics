# mcmc_diagnostics

A small package for running MCMC diagnostics on Aria, Stan or MaCh3 (experimental) chains. Can be used for:

1. Removing the burin-in
2. Getting the effective sample size
3. Determining the thinning factor to use (if any)
4. Checking if the chains have converged, and whether we have enough MCMC steps
5. Identifying various degeneracies related to MCMC sampling
6. Merging MCMC chains into one

## Installing dependencies

Create a virtual python environment and install the requirements via pip as usual:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Once the requirements are installed inside of the virtual environment, `source .venv/bin/activate` can be run at the next login time to re-activate it.

## Running

There are two executables, one for MCMC diagnostics, and another for merging chains into one, with burn-in and thinning options.

1. [README for diagnostics](docs/diagnose_mcmc.md)
2. [README for chain merging](docs/merge_chains.md)