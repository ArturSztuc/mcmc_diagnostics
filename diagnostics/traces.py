import uproot
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from tqdm import tqdm

from .sampler_metadata import SamplerMetadata

def make_trace_plots(metadata, xbins=1000, ybins=100, output_file="trace_plots.pdf"):
    pdf = PdfPages(f"trace_plots.pdf")

    # Open one file to get parameter names
    with uproot.open(metadata.files[0]) as f:
        example_chain = f[metadata.ttree_location]
        keys = [key for key in example_chain.keys() if key not in metadata.ignored_branches]
        length = len(example_chain[keys[0]].array())

    # Initialise histograms for each parameter
    histograms = {key: np.zeros((xbins, ybins)) for key in keys}
    xedges = np.linspace(0, length, xbins+1)#len(files) * 100_000, bins + 1)  # Approximate iteration range
    yedges_dict = {}

    # Process one file at a time, adding to the histograms
    for file in tqdm(metadata.files, desc="Processing MCMC chains for trace heatmaps"):
        with uproot.open(file) as f:
            chain = f[metadata.ttree_location]

            for key in keys:
                data = chain[key].array(library="np")
                # Take abs of dm32, the bimodality is difficult to look at
                if "32" in key:
                    data = np.abs(data)
                iterations = np.arange(len(data))

                # Set y-bins dynamically (first file defines the range)
                if key not in yedges_dict:
                    yedges_dict[key] = np.linspace(np.min(data), np.max(data), 101)

                # Update histogram incrementally (no need to store full chains!)
                hist, _, _ = np.histogram2d(iterations, data, bins=(xedges, yedges_dict[key]))
                histograms[key] += hist  # Accumulate counts

    # Generate plots
    for key in tqdm(keys, desc="Generating trace heatmaps"):
        plt.figure(figsize=(10, 4))

        plt.imshow(histograms[key].T, aspect="auto", origin="lower",
                   extent=[xedges[0], xedges[-1], yedges_dict[key][0], yedges_dict[key][-1]],
                   cmap="inferno", interpolation="nearest")

        plt.colorbar(label="Density (number of samples)")
        plt.title(f"Heatmap trace plot for {key}")
        plt.xlabel("Iteration")
        plt.ylabel(key)
        plt.tight_layout()
        pdf.savefig()
        plt.close()
    pdf.close()