import uproot
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.backends.backend_pdf import PdfPages
from itertools import combinations
from tqdm import tqdm

# Calculate the quantiles
def get_quantile_thresholds(arr, quantile_levels):
    assert all(0 <= f <= 1 for f in quantile_levels)

    arr = arr.flatten()
    total = np.sum(arr)
    breakpoint_sums = quantile_levels * total
    sorting_indices = np.argsort(arr)[::-1]
    cdf = np.cumsum(arr[sorting_indices])

    return (arr[sorting_indices[np.searchsorted(cdf, qlevel)]] for qlevel in breakpoint_sums)


def update_histograms(data, key, hist_dict, xedges_dict, split_type):
    hist, _ = np.histogram(data, bins=xedges_dict[key])
    hist_dict[key][split_type] += hist
 
def update_2d_histograms(data_x, data_y, pair, hist_dict, xedges_dict, split_type):
    hist, _, _ = np.histogram2d(data_x, data_y, bins=(xedges_dict[pair[0]], xedges_dict[pair[1]]))
    hist_dict[pair][split_type] += hist

def make_split_posteriors(metadata, n_bins, burn_in, output_file="split_posterior.pdf"):
    pdf = PdfPages(output_file)
    example_chain = uproot.open(metadata.files[0])[metadata.ttree_location]
    keys = [key for key in example_chain.keys() if key not in metadata.ignored_branches]
    keys_important = [key for key in keys if any(samplerkey in key for samplerkey in metadata.key_branches)]
    pairs_important = list(combinations(keys_important, 2))

    quantile_thresholds = np.array([0.6827, 0.9545, 0.9973])
    sigma_styles = {"1s": "solid", "2s": "dashdot", "3s": "dotted"}
    chain_names = ["Full posterior", "Left side of chains", "Right side of chains", "First half of chains", "Second half of chains"]
    chain_colors = {"full": "black", "left": "darkred", "right": "lightcoral", "first": "blue", "second": "cornflowerblue"}
    
    histograms = {key: {split: np.zeros(n_bins) for split in chain_colors.keys()} for key in keys}
    histograms_2d = {pair: {split: np.zeros((n_bins, n_bins)) for split in chain_colors.keys()} for pair in pairs_important}
    
    xedges_dict = {}
    nfiles = len(metadata.files)

    for file_idx, file in enumerate(tqdm(metadata.files, desc="Processing MCMC files")):
        chain = uproot.open(file)[metadata.ttree_location]
        
        for key in keys:
            data = chain[key].array(library="np")[burn_in:]
            half = len(data) // 2
            if "32" in key:
                data = np.abs(data)
            
            if key not in xedges_dict:
                xedges_dict[key] = np.linspace(np.min(data), np.max(data), n_bins + 1)

            update_histograms(data[:half], key, histograms, xedges_dict, "left")
            update_histograms(data[half:], key, histograms, xedges_dict, "right")
            histsum = histograms[key]["left"] + histograms[key]["right"]
            histograms[key]["full"] += histsum
            
            if file_idx / nfiles < 0.5:
                histograms[key]["first"] += histsum
            else:
                histograms[key]["second"] += histsum

        for pair in pairs_important:
            data_x = chain[pair[0]].array(library="np")[burn_in:]
            data_y = chain[pair[1]].array(library="np")[burn_in:]
            half = len(data_x) // 2
            
            update_2d_histograms(data_x[:half], data_y[:half], pair, histograms_2d, xedges_dict, "left")
            update_2d_histograms(data_x[half:], data_y[half:], pair, histograms_2d, xedges_dict, "right")
            histsum = histograms_2d[pair]["left"] + histograms_2d[pair]["right"]
            histograms_2d[pair]["full"] += histsum
            
            if file_idx / nfiles < 0.5:
                histograms_2d[pair]["first"] += histsum
            else:
                histograms_2d[pair]["second"] += histsum

    def plot_contours(pair):
        plt.figure(figsize=(10, 8))
        for chain, color in chain_colors.items():
            contour_breakpoints = np.array(list(get_quantile_thresholds(histograms_2d[pair][chain.lower()], quantile_thresholds)))[::-1]
            plt.contour(xedges_dict[pair[0]][:-1], xedges_dict[pair[1]][:-1],
                        histograms_2d[pair][chain.lower()].T, levels=contour_breakpoints,
                        colors=color, linestyles=[sigma_styles["3s"], sigma_styles["2s"], sigma_styles["1s"]], linewidths=1)

        all_handles = [mlines.Line2D([], [], color=color, visible=False, label=f"{chain_name}:")
                       for chain_name in chain_names]

        legend_handles = [mlines.Line2D([], [], color=color, linestyle=sigma_styles[sigma], label=sigma.replace('s', r'$\sigma$')) 
                         for sigma in ["1s", "2s", "3s"] for chain, color in chain_colors.items() ]
        
        legend_handles = all_handles + legend_handles

        plt.title(f"Split posterior for {pair[0]} vs {pair[1]}")
        plt.xlabel(pair[0])
        plt.ylabel(pair[1])
        plt.legend(handles=legend_handles, title="Credible intervals", loc="upper right", ncol=4)
        plt.tight_layout()
        pdf.savefig()
        plt.close()

    for pair in tqdm(pairs_important, desc="Generating 2D posterior plots"):
        plot_contours(pair)

    for key in tqdm(keys, desc="Generating split posteriors"):
        plt.figure(figsize=(10, 8))

        for (split, color), label in zip(chain_colors.items(), chain_names):
            if np.sum(histograms[key][split]) > 0:
                histograms[key][split] /= np.sum(histograms[key][split])
                bincenters = 0.5 * (xedges_dict[key][1:] + xedges_dict[key][:-1])
                plt.step(bincenters, histograms[key][split], where='mid', color=color, linestyle='--' if split != "full" else '-', label=label)

        plt.title(f"Split posterior for {key}")
        plt.xlabel(key)
        plt.ylabel("Posterior probability density")
        plt.legend()
        plt.tight_layout()
        pdf.savefig()
        plt.close()

    pdf.close()
