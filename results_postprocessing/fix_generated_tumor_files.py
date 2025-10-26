import os.path, re
from typing import List

import pandas as pd


def break_line(line: str) -> List[str]:
    return line.split("\t")


def rejoin_line(broken_line: List[str]) -> str:
    return "\t".join(broken_line)

def num_columns_in_line(line):
    return len(break_line(line))


def flatten_file(input_file: str, output_file: str):
    pattern_column = 3
    allowed_number_of_motifs = 6
    first_motif_column = 6
    output_lines = []
    with open(input_file, 'r') as f:
        first_line = f.readline()
        num_columns_in_df = len(first_line.split("\t"))
        output_lines.append(first_line)
        for line in f:
            if num_columns_in_line(line) > num_columns_in_df: # improper line
                broken_line = break_line(line)
                num_columns_in_current_line = len(broken_line)
                if num_columns_in_current_line > (num_columns_in_df + 3) and not re.search(r'\d', broken_line[num_columns_in_df+pattern_column-1]): # duplicate_line
                    output_lines.append(rejoin_line(broken_line[:num_columns_in_df]))
                    output_lines.append(rejoin_line(broken_line[num_columns_in_df:]))
                else: # just has too many motif repeats
                    output_lines.append(rejoin_line(broken_line[:first_motif_column]+
                                                    broken_line[first_motif_column:first_motif_column+allowed_number_of_motifs]+
                                                    broken_line[-allowed_number_of_motifs:]))

            else:
                output_lines.append(line)
    with open(output_file, 'w+') as output:
        output.write("\n".join(output_lines))



def replace_floats_with_ints(input_file: str, output_file: str):
    with open(input_file, "r", encoding="utf-8") as f:
        content = f.read()
    updated_content = content.replace(".0", "")
    with open(output_file, "w", encoding="utf-8") as f:
        f.write(updated_content)
    # df = pd.read_csv(input_file, delimiter="\t")
    # df = df.apply(lambda col: col.astype('Int64') if col.dtype == 'float' else col)
    # df.to_csv(output_file, index=False, sep="\t")


def write_na_string(input_file: str, output_file: str):
    max_tabs = 15
    with open(input_file, "r", encoding="utf-8") as f:
        content = f.read()
    content = content.replace("\t\n", "\tNA\n")
    def replacer(match):
        tabs = match.group(0)
        num_tabs = len(tabs)
        return ("\tNA" * (num_tabs - 1)) + "\t"
    pattern = r"\t{2," + str(max_tabs) + r"}"
    updated_content = re.sub(pattern, replacer, content)
    with open(output_file, "w", encoding="utf-8") as f:
        f.write(updated_content)


def fix_fake_tumor_file(input_file: str, output_file: str):
    flatten_file(input_file, output_file)
    replace_floats_with_ints(output_file, output_file)
    write_na_string(output_file, output_file)


def main():
    src_dir = "C:/Users/avrah/MaruvkaLab/post_processing_code_for_australians/data/positives"
    fix_fake_tumor_file(os.path.join(src_dir, "TCGA-AP-A05N_1x_0.35purity.hist.tsv"),
                        os.path.join(src_dir, "tumor_0.35purity.hist.tsv"))



if __name__ == '__main__':
    main()