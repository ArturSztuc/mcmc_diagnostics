import uproot
import numpy as np

from tqdm import tqdm

from colorama import Fore, Back

def get_step_acceptances(metadata, burn_in):
  if burn_in == 0:
    print(Back.RED + "Warning: burn-in is 0, this may lead to incorrect step acceptances printed" + Back.RESET)

  acceptances = []
  
  # Go over each file and calculate the step acceptance
  for file in tqdm(metadata.files, desc="Getting step acceptances"):
    diffs = []
    with uproot.open(file) as f:
      example_chain = f[metadata.ttree_location]
      keys = [key for key in example_chain.keys() if key not in metadata.ignored_branches]

      # Get the differences between consecutive samples for the first non-ignored branch
      data = np.asarray(example_chain[keys[0]].array()[burn_in:])
      data_diff = np.diff(data)
      diffs.append(data_diff)

      # Calculate the percentage of accepted steps for the chain
      acceptance = np.sum(data_diff != 0) / len(data_diff)
      acceptance *= 100.0;
      acceptances.append(acceptance)

  # Calculate the total acceptance across all files
  diffs = np.concatenate(diffs)
  total_acceptance = np.sum(diffs != 0) / len(diffs)
  total_acceptance *= 100.0
  acceptances.insert(0, total_acceptance)

  return acceptances

def print_step_acceptance(metadata, burn_in=0):
  acceptances = get_step_acceptances(metadata,  burn_in)

  print(f"Step acceptances for chains from {metadata.sampler_name} sampler:")

  foreground = Fore.GREEN

  # Calculate the distance from the perfect acceptance as a percentage
  distance_from_perfect = np.abs(acceptances[0] - metadata.perfect_acceptance) / metadata.perfect_acceptance * 100.0
  if distance_from_perfect > 25.0:
    foreground = Fore.RED

  # Print the acceptances for each file and the total step acceptance
  for i, file in enumerate(metadata.files):
    print(f"  - {file}: {foreground}{acceptances[i + 1]:.2f}%{Fore.RESET}")
  print(f"  - Total accepted steps: {foreground}{acceptances[0]:.2f}%{Fore.RESET} (perfect acceptance for {metadata.sampler_name}: {metadata.perfect_acceptance:.2f}%)")
  if distance_from_perfect > 25.0:
    print(Back.RED + Fore.WHITE + f"Warning: total step acceptance far away from the perfect acceptance of {metadata.sampler_name}!" + Back.RESET + Fore.RESET)
    if acceptances[0] < metadata.perfect_acceptance:
      print(Back.RED + Fore.WHITE + "Step-sizes need to be decreased." + Back.RESET + Fore.RESET)
    else:
      print(Back.RED + Fore.WHITE + "Step-sizes need to be increased." + Back.RESET + Fore.RESET)
  else:
    print(Back.GREEN + f"Total step acceptance is close to the perfect acceptance of {metadata.sampler_name}!" + Back.RESET)