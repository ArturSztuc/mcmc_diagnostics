import uproot

import numpy as np

from tqdm import tqdm

from colorama import Fore, Back

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def autocorr_fft_padded(x, lags):
    n=len(x)
    # pad 0s to 2n-1
    ext_size=2*n-1
    # nearest power of 2
    fsize=2**np.ceil(np.log2(ext_size)).astype('int')

    # subtract mean and calculate variance
    xp=x-np.mean(x)
    var=np.var(x)

    # do fft and ifft
    cf=np.fft.fft(xp,fsize)
    sf=cf.conjugate()*cf
    corr=np.fft.ifft(sf).real
    corr=corr/var/n

    return corr[:len(lags)]

def autocorr_naive(x, lags):
  mean = x.mean()
  var= np.sum((x - mean)**2)

  autocors = []
  for lag in lags:
    lhs = (x[:len(x) - lag] - mean)
    rhs = (x[lag:] - mean)
    autocors.append(np.sum(lhs * rhs) / var)
  
  return np.asarray(autocors)

def plot_autocorrelation(autocorrelation, key, pdf):

  plt.plot(autocorrelation)
  
  plt.title(f"Autocorrelations per chain for {key}")

  plt.axhline(0, color='black', linestyle='--')
  plt.xlabel("Lag")
  plt.ylabel("Autocorrelation")

  plt.tight_layout()
  pdf.savefig()
  plt.close()

def plot_all_autocorrelations_one_plot(autocorrelations, branch_keywords, pdf):
  first_systematic = True
  for [key, autocorrelation] in autocorrelations.items():
    interesting = False
    for oscpar in branch_keywords:
      if oscpar in key:
        interesting = True
        break
    if interesting:
      newkey = key
      if newkey.startswith("_"):
        newkey = newkey[1:]

      plt.plot(autocorrelation, label=f"{newkey}")
    else:
      if first_systematic:
        plt.plot(autocorrelation, label="systematic", color="black", linewidth=0.1)
        first_systematic = False
      else:
        plt.plot(autocorrelation, color="black", linewidth=0.1)

  plt.title("Autocorrelation per parameter")
  plt.axhline(0, color='black', linestyle='--')
  plt.xlabel("Lag")
  plt.ylabel("Autocorrelation")
  plt.legend()
  #plt.tight_layout()
  pdf.savefig()
  plt.close()

def get_autocorrelations(metadata, max_lag=100, burn_in=0):
  autocorrelations = {}

  if burn_in == 0:
    print(Back.RED + "Warning: burn-in is 0, this may lead to incorrect autocorrelation results" + Back.RESET)

  # Get the key names
  with uproot.open(metadata.files[0]) as f:
      example_chain = f[metadata.ttree_location]
      keys = [key for key in example_chain.keys() if key not in metadata.ignored_branches]

  # Initialise the autocorrelations
  for key in keys:
    autocorrelations[key] = np.zeros(max_lag)

  # Iterate over the files
  for filename in tqdm(metadata.files, desc="Getting autocorrelations"):
    chain = uproot.open(filename)[metadata.ttree_location]

    # Add the autocorrelations to the total for each key
    for key in keys:
      data = np.asarray(chain[key].array()[burn_in:])
      autocorrelations[key] += autocorr_fft_padded(data, range(max_lag))
  
  # Average over the chains
  for key in autocorrelations:
    autocorrelations[key] /= len(metadata.files)
  
  return autocorrelations

def make_autocorrelation_plots(metadata, max_lag=100, burn_in=0, output_file="autocorrelations.pdf"):
  autocorrelations = get_autocorrelations(metadata, max_lag, burn_in)

  with open(output_file, "wb") as f:
    pdf = PdfPages(f)
    # Plot all autocorrelations in one plot
    plot_all_autocorrelations_one_plot(autocorrelations, metadata.key_branches, pdf)

    # Plot each autocorrelation separately
    for key in autocorrelations:
      plot_autocorrelation(autocorrelations[key], key, pdf)
    pdf.close()
    print(f"Autocorrelations saved to {output_file}")
