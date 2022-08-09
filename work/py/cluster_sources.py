from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage

# read and translate R functions
r_file_source = open("./rscripts/cluster_sources.R", "r").read()

cluster_sources = SignatureTranslatedAnonymousPackage(r_file_source, "cluster_sources")

results_dir: str = "/Users/joshl/docker/madrid/local_files/results"
context_names: list[str] = ["immNK", "naiveB"]
source_type: str = "zFPKM"
use_trna: bool = True
use_mrna: bool = True
binarize_data: bool = True
default_bin: int = -3
cluster_sources.main(
    results_dir=results_dir,
    context_names=context_names,
    source_type=source_type,
    use_trna=use_trna,
    use_mrna=use_mrna,
    binarize_data=binarize_data,
    default_bin=default_bin
)
