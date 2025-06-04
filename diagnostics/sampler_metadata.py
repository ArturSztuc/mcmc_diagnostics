import uproot

IGNORE_BRANCHES_STAN = ["accept_stat__", "stepsize__", "treedepth__", "n_leapfrog__", "divergent__", "energy__", "stepnum"]
BRANCH_KEYWORDS_STAN = ["logprob", "Th13", "Th23", "dCP", "DmSq32"]

IGNORE_BRANCHES_ARIA = ["MH", "stepnum"]
BRANCHES_KEYWORDS_ARIA = ["logprob", "th13", "th23", "delta(pi)", "dmsq32"]

class SamplerMetadata:
    """
    Class to hold metadata for a sampler.
    """

    def __init__(self, files: list[str]):
        """
        Initialize the SamplerMetadata object.

        Args:
            name (str): The name of the sampler.
            description (str): A description of the sampler.
        """

        self.files = [f for f in files if f.endswith(".root")]
        self.__fill_metadata()

    def __init(self, files: list[str], 
               sampler_name: str, 
               ignored_branches: list[str],
               key_branches: list[str],
               ttree_location: str):
              
        """
        Initialize the SamplerMetadata object.

        Args:
            name (str): The name of the sampler.
            description (str): A description of the sampler.
        """

        self.files = [f for f in files if f.endswith(".root")]
        self.sampler_name = sampler_name
        self.ignored_branches = ignored_branches
        self.key_branches = key_branches
        self.ttree_location = ttree_location

    def print_metadata(self):
        print(f"Sampler metadata initialized:")
        print(f"  - Sampler name: {self.sampler_name}")
        print(f"  - Ignored branches: {self.ignored_branches}")
        print(f"  - Key branches: {self.key_branches}")
        print(f"  - TTree location: {self.ttree_location}")
        print(f"  - Number of files: {len(self.files)}")

    def get_default_maxlag(self):
        if self.sampler_name == "stan":
            return 200
        elif self.sampler_name == "aria":
            return 20_000
        else:
            raise ValueError("Unknown sampler type. Please check the files.")

    def __fill_metadata(self):

        with uproot.open(self.files[0]) as f:
            if "samples/samples" in f:
                self.sampler_name = "stan"
                self.ignored_branches = IGNORE_BRANCHES_STAN
                self.key_branches = BRANCH_KEYWORDS_STAN
                self.ttree_location = "samples/samples"
            elif "run/samples" in f:
                self.sampler_name = "aria"
                self.ignored_branches = IGNORE_BRANCHES_ARIA
                self.key_branches = BRANCHES_KEYWORDS_ARIA
                self.ttree_location = "run/samples"
            else:
                raise ValueError("Unknown sampler type. Please check the files.")

    def __repr__(self):
        return f"SamplerMetadata(name={self.name}, description={self.description})"
