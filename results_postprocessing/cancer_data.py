import sys, os
from typing import Dict

from SamplesDB import SamplesDB, Sample
from enums import MSI_CLASSIFICATION




def load_cancer_type_metadata(fp="C:/Users/avrah/MaruvkaLab/Texas_samples_organization/tcga_metadata_upd.txt") -> Dict[str, str]:
    if sys.platform.startswith("linux"):
        # fp="/mnt/c/Users/avrah/MaruvkaLab/Texas_samples_organization/tcga_metadata_upd.txt"
        # fp="/mnt/c/Users/avrah/MaruvkaLab/Texas_samples_organization/tcga_metadata_upd.txt"
        fp = "/home/avraham/MaruvkaLab/msmutect_postprocessing/data/tcga_metadata_upd.txt"
    cancer_type_lookup_table = dict()
    with open(fp, 'r') as metadata_file:
        metadata_lines = metadata_file.readlines()
    for metadata in metadata_lines:
        name, cancer_type = metadata.split("\t")[:2]
        cancer_type_lookup_table[name] = cancer_type
    return cancer_type_lookup_table


def load_dir_wclassfication(directory: str, classification: str, extension: str = ".called.filt.mut.tsv.gz", forced_cancer_type: str = None, ):
    ret = []
    cancer_type_set = load_cancer_type_metadata()
    for file in os.listdir((directory)):
        if not file.endswith(extension):
            continue

        if classification == MSI_CLASSIFICATION.NEGATIVE_CONTROL:
            cancer_type = "GIB"
        elif forced_cancer_type is not None:
            cancer_type = forced_cancer_type
        else:
            cancer_type = cancer_type_set[Sample.patient_name_from_filename(file)]
        ret.append(Sample(os.path.join(directory, file), cancer_type, classification))
    return ret


def load_tcga_samples(msi_dir: str, mss_dir: str, mss_only: bool, forced_cancer_type: str = None) -> SamplesDB:
    if not mss_only:
        msi_samples = load_dir_wclassfication(msi_dir, MSI_CLASSIFICATION.MSI, forced_cancer_type=forced_cancer_type)
    else:
        msi_samples = []
    mss_samples = load_dir_wclassfication(mss_dir, MSI_CLASSIFICATION.MSS, forced_cancer_type=forced_cancer_type)
    ret = SamplesDB()
    ret.add_multiple_samples(msi_samples+mss_samples)
    return ret


def data_directory():
    if sys.platform.startswith("linux"):
        return "/home/avraham/MaruvkaLab/msmutect_postprocessing/data"
    else:
        return "C:/Users/avrah/MaruvkaLab/post_processing_code_for_australians/data"


def results_directory():
    if sys.platform.startswith("linux"):
        return "/home/avraham/MaruvkaLab/msmutect_postprocessing/results"
    else:
        return "C:/Users/avrah/MaruvkaLab/post_processing_code_for_australians/results"


def graphs_directory():
    return "C:/Users/avrah/MaruvkaLab/post_processing_code_for_australians/graphs"

def load_all_samples(mss_only: bool = False) -> SamplesDB:
    return load_tcga_samples(os.path.join(data_directory(), "msi"), os.path.join(data_directory(), "mss"), mss_only)


def load_test_samples(mss_only: bool = False) -> SamplesDB:
    return load_tcga_samples(os.path.join(data_directory(), "msi_example_files"),
                             os.path.join(data_directory(), "mss_example_files"), mss_only, forced_cancer_type="TEST")


def load_all_samples_gib() -> SamplesDB:
    negative_controls = load_dir_wclassfication(os.path.join(data_directory(), "gib_files_filtered"), MSI_CLASSIFICATION.NEGATIVE_CONTROL, extension=".filtered.full.mut.tsv")
    ret = SamplesDB()
    ret.add_multiple_samples(negative_controls)
    return ret
