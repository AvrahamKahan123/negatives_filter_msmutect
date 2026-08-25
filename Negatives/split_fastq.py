#!/Local/bfe_maruvka/anaconda3/bin/python3 -u

import time, random, os, sys
import numpy as np
from typing import List
from dataclasses import dataclass
from collections import OrderedDict


class NumPicker_300x:
    # efficient way to pick number with weighted probabilities
    # all my BAM files are 300x anyways, so might as well optimize for that
    def __init__(self, relative_depths: List[int]):
        if sum(relative_depths)>300:
            raise RuntimeError("Too many reads requested")
        self.map = self.initialize_map(relative_depths)

    def initialize_map(self, relative_depths: List[int]) -> np.ndarray:
        map = np.zeros((300), np.int32)  # build a map to map the random numbers to choices of which file
        current_increment = 0
        for i, increment in enumerate(relative_depths):
            new_increment = current_increment + increment
            map[current_increment:new_increment] = i
            current_increment = new_increment
        return map

    def pick_num(self):
        random_num = random.random()
        return self.map[int(random_num*300)]


class ReadBucket:
    def __init__(self, name: str, dest_dir: str, is_null: bool):
        if not is_null:
            self.fastq_1 = os.path.join(dest_dir, name+"_1.fastq")
            self.opened_fastq_1 = open(self.fastq_1, 'w')
            self.fastq_2 = os.path.join(dest_dir, name+"_2.fastq")
            self.opened_fastq_2 = open(self.fastq_2, 'w')
        self.is_null = is_null

    def add_reads(self, read_1: str, read_2: str):
        if self.is_null:
            return
        else:
            self.opened_fastq_1.write(read_1)
            self.opened_fastq_2.write(read_2)

    def close(self):
        if self.is_null:
            return
        self.opened_fastq_1.close()
        self.opened_fastq_2.close()


def next_4_lines(f) -> str:
    lines = [next(f, None) for _ in range(4)]
    if any(line is None for line in lines):
        return None
    return "".join(lines)


def pick_500million_nums(ra):
    res = [0 for _ in range(7)]
    for i in range(500_000_000):
        x=ra.pick_num()
        res[x]+=1
    print(res)
    exit()


def split_300x_into_3_pseudopatients(fastq_1_fp: str, fastq_2_fp: str, destination_dir: str):

    x75_depth_buckets = [ReadBucket(f"x75_{i}", destination_dir, is_null=False) for i in range(2)]
    x30_depth_buckets = [ReadBucket(f"x30_{i}", destination_dir, is_null=False) for i in range(5)]

    all_buckets_wdepths = OrderedDict([(bucket, 75) for bucket in x75_depth_buckets]+[(bucket, 30) for bucket in x30_depth_buckets])
    all_buckets = list(all_buckets_wdepths.keys())
    all_bucket_depths = list(all_buckets_wdepths.values())
    opened_fastq_1 = open(fastq_1_fp, 'r')
    opened_fastq_2 = open(fastq_2_fp, 'r')
    read_assigner = NumPicker_300x(all_bucket_depths)

    # pick_500million_nums(read_assigner)

    f1_lines = next_4_lines(opened_fastq_1)
    f2_lines = next_4_lines(opened_fastq_2)
    while True:
        if f1_lines is None or f2_lines is None:
            for bucket in all_buckets:
                bucket.close()
                exit()
        else:
            bucket_number = read_assigner.pick_num()
            all_buckets[bucket_number].add_reads(f1_lines, f2_lines)
            f1_lines = next_4_lines(opened_fastq_1)
            f2_lines = next_4_lines(opened_fastq_2)


if __name__ == '__main__':
    split_300x_into_3_pseudopatients(sys.argv[1], sys.argv[2], sys.argv[3])
