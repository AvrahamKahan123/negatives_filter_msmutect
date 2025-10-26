import os, time
import numpy as np
from typing import List
import pandas as pd
import multiprocessing as mp


def get_filename_no_ext(fp: str):
    fp = os.path.basename(fp)
    period_loc = fp.find(".")
    return fp[:period_loc]


def analyze_single_file(fp: str, log_thresholds: List[float], column_name: str ) -> List[int]:
    """

    :param fp: file path of  file to analyze
    :param log_thresholds: log10 thresholds to test number of mutations for
    :return: len(thresholds)+1 list of number of mutations for [original, threshold_0, threshold_1, ..., threshold_n]
    """
    mutations_df = pd.read_csv(fp, delimiter="\t")
    thresholds = list(10**np.array(log_thresholds))
    mutations_per_threshold = [len(mutations_df)]+[(mutations_df[column_name] <= thresh).sum() for thresh in thresholds]
    return [get_filename_no_ext(fp)] + mutations_per_threshold


def analyze_list_of_files(files: List[str], log_thresholds: List[float], original_column_name: str = "LOG10 P_VAL=-1.5", column_name="FISHER_TEST_P_VALUE") -> pd.DataFrame:
    other_cols_prefix: str = "LOG10 P_VAL="
    results = [analyze_single_file(fp, log_thresholds, column_name) for fp in files]
    return pd.DataFrame(results, columns=["CASE"]+[original_column_name]+[f"{other_cols_prefix}{p}" for p in log_thresholds])


class BatchManager:
    def __init__(self, src_dir: str, batch_size: int):
        self.src_dir = src_dir
        self.batch_size = batch_size
        self.location = 0
        self.all_files = [os.path.join(src_dir, f) for f in os.listdir(src_dir)]

    def get_next_batch(self):
        ret = self.all_files[self.location: self.location+self.batch_size]
        self.location+=self.batch_size
        return ret

    def done(self):
        return self.location>len(self.all_files)

    def get_all_batches(self) -> List[List[str]]:
        ret = []
        while not self.done():
            ret.append(self.get_next_batch())
        return ret


class ProcessManager:
    def __init__(self, num_processes: int):
        self.num_processes = num_processes
        self.active_processes = []
        self.finished_processes = []
        self.queues = []

    def done(self):
        return len(self.active_processes)==0

    def has_room_for_new_process(self):
        return len(self.active_processes) < self.num_processes

    def add_process(self, func, args):
        new_process = mp.Process(target=func, args=args)
        new_process.start()
        self.active_processes.append(new_process)

    def reap(self):
        for process in self.active_processes:
            if not process.is_alive():
                self.finished_processes.append(process)

    def extract_results(self) -> List[pd.DataFrame]:
        pass


def main(src_dir: str = "data/MSMuTect_called_mut_filt_fixed", output_file: str="all_results.csv", batch_size: int = 100, metric="KS"):
    file_manager = BatchManager(src_dir, batch_size)
    all_batches = file_manager.get_all_batches()
    # thresholds = list(10**log_thresholds)
    if metric=="FISHER":
        log_thresholds = np.flip(np.arange(-5, -1.5, 0.25), axis=0)
        arguments = [(batch, log_thresholds) for batch in all_batches]
    elif metric== "KS":
        log_thresholds = np.flip(np.arange(-7, 0, 0.25), axis=0)
        arguments = [(batch, log_thresholds, "ORIGINAL", "KS_TEST_PVALUE") for batch in all_batches]
    else:
        raise RuntimeError("Unknown metric")

    with mp.Pool(processes=12) as pool:
        results = pool.starmap(analyze_list_of_files, arguments)
        # process_manager.reap()
        # if process_manager.done() and file_manager.done():
        #     results = process_manager.extract_results()
        #     break
        # elif process_manager.has_room_for_new_process():
        #     new_batch = file_manager.get_next_batch()
        #     if new_batch is not None:
        #         process_manager.add_process(new_batch)
        # else:
        #     time.sleep(2)
    combined = pd.concat(results, ignore_index=True)
    combined.to_csv(output_file, sep="\t", index=False)


if __name__ == '__main__':
    # main("data/hg003_filtered")
    main(src_dir="C:/Users/avrah/MaruvkaLab/post_processing_code_for_australians/data/gib_files", output_file="../results/gib_ks.tsv")
    main(src_dir="/data/msi", output_file="../results/msi_ks.tsv")
    main(src_dir="/data/mss", output_file="../results/mss_ks.tsv")

