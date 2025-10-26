import glob, os


def filtered_mutation_files_list(directory: str):
    all_files = glob.glob(os.path.join(directory, "*.tsv"))
    return all_files
    # return [os.path.join(directory, f) for f in all_files]

