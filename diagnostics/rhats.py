from colorama import Fore, Back
import uproot
from tqdm import tqdm
import numpy as np

from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt

def calculate_gelman_rubin(x):
  m, n = x.shape
  chain_means = np.mean(x, axis=1)
  W = np.mean(np.var(x, axis=1, ddof=1))
  B = n / (m - 1) * np.sum((chain_means - np.mean(chain_means))**2)
  var = (n - 1) / n * W + B / n
  return np.sqrt(var / W)
  
def get_rhat_naive(files, burn_in, ignore_branches_list):
  example_chain = uproot.open(files[0])["samples/samples"]
  keys = [key for key in example_chain.keys() if key not in ignore_branches_list]

  dict_rhat = {}
  dict_within_chain_rhat = {}
  for key in tqdm(keys, desc="Calculating Rhats"):
    chains = []
    within_chain_rhats = []
    for file in files:
      chain = uproot.open(file)["samples/samples"]

      # Get the within-chain rhats
      chainhalf_size = len(chain[key].array()[burn_in:]) // 2

      rhat_inside = calculate_gelman_rubin(np.asarray([chain[key].array()[burn_in:chainhalf_size + burn_in], chain[key].array()[chainhalf_size + burn_in:]]))

      within_chain_rhats.append(rhat_inside)
      chains.append(chain[key].array()[burn_in:])
    
    # Get the between-chain rhat
    chains = np.asarray(chains)
    rhat = calculate_gelman_rubin(chains)

    newkey = key
    if newkey.startswith("_"):
      newkey = newkey[1:]

    dict_rhat[newkey] = rhat
    dict_within_chain_rhat[newkey] = within_chain_rhats
  
  return dict_rhat, dict_within_chain_rhat

def get_rhat(metadata, burn_in):
  if burn_in == 0:
    print(Back.RED + "Warning: burn-in is 0, this may lead to incorrect rhat results" + Back.RESET)

  # Open a sample file to get parameter names
  with uproot.open(metadata.files[0]) as f:
      example_chain = f[metadata.ttree_location]
      keys = [key for key in example_chain.keys() if key not in metadata.ignored_branches]
  
  # Data structures to store the statistics
  M = len(metadata.files) 
  sum_means = {}
  sum_vars = {}
  chain_means = {}
  for key in keys:
      newkey = key
      if newkey.startswith("_"):
          newkey = newkey[1:]
      sum_means[newkey] = 0.0
      sum_vars[newkey] = 0.0
      chain_means[newkey] = []

  dict_rhat = {}
  dict_within_chain_rhat = {}

  for file in tqdm(metadata.files, desc="Calculating Rhats (Gellman-Rubin)"):
      with uproot.open(file) as f:
          chain = f[metadata.ttree_location]

          for key in keys:
              # Read only the required portion of the chain
              data = chain[key].array(library="np")[burn_in:]

              # Compute within-chain R̂ using split chains
              half = len(data) // 2
              chain1, chain2 = data[:half], data[half:]
              within_rhat = calculate_gelman_rubin(np.stack([chain1, chain2]))

              newkey = key
              if newkey.startswith("_"):
                newkey = newkey[1:]
              # Store within-chain rhat
              dict_within_chain_rhat.setdefault(newkey, []).append(within_rhat)

              # Update statistics for between-chain rhat
              chain_mean = np.mean(data)
              chain_var = np.var(data, ddof=1)

              sum_means[newkey] += chain_mean
              sum_vars[newkey] += chain_var
              chain_means[newkey].append(chain_mean)

  # Compute final between-chain rhat
  for key in keys:
      newkey = key
      if newkey.startswith("_"):
        newkey = newkey[1:]
      global_mean = sum_means[newkey] / M
      W = sum_vars[newkey] / M
      B = np.sum((np.array(chain_means[newkey]) - global_mean) ** 2) * (len(data) / (M - 1))
      var = (len(data) - 1) / len(data) * W + B / len(data)
      dict_rhat[newkey] = np.sqrt(var / W)

  return dict_rhat, dict_within_chain_rhat

def plot_rhat_matrix(dict_between_rhats, dict_within_rhats):
  import seaborn as sns
  import pandas as pd

  df_between = pd.DataFrame(dict_between_rhats, index=["Total"])
  df_within = pd.DataFrame(dict_within_rhats)

  df = pd.concat([df_between, df_within])

  df = df.T

  cmap = ListedColormap(["white", "lightcoral", "darkred"])

  bounds = [0, 1.01, 1.05, 1.1]

  norm = plt.cm.colors.BoundaryNorm(bounds, cmap.N)

  plt.figure(figsize=(12, 8))
  ax = sns.heatmap(df, cmap=cmap, norm=norm, cbar=True, linewidths=0.5, linecolor="black")
  ax.set_title("Rhats")
  ax.set_xticks(np.arange(df.shape[1]) + 0.5)
  ax.set_yticks(np.arange(df.shape[0]) + 0.5)

  ax.set_xticklabels(df.columns, rotation=45, ha="right", fontsize=8)
  ax.set_yticklabels(df.index, rotation=0, fontsize=8)

  plt.tight_layout()
  plt.savefig("rhats.pdf")
  plt.close()

def make_rhat_plots(metadata, burn_in):
  dict_between_rhats, dict_within_rhats = get_rhat(metadata, burn_in)
  plot_rhat_matrix(dict_between_rhats, dict_within_rhats)