import glob, os, re


def count_files(src_dir: str):
    all_files = glob.glob(os.path.join(src_dir, "*.npy"))
    sorted_all_files =  None

def extract_number(fp: str, num: int):
    numbers = re.findall(r"\d+", fp)
    return int(numbers[num])

def verify_completion(completed_npys_list_file: str):
    with open(completed_npys_list_file, 'r') as f:
        completed_files = [s.rstrip() for s in f.readlines()]
    sorted_completed_files = list(sorted(completed_files, key=lambda s: extract_number(s, 0)))
    with open("c2.txt", 'w+') as croc:
        croc.writelines([s+"\n" for s in sorted_completed_files])
    last = extract_number(sorted_completed_files[0], 1)
    for f in sorted_completed_files[1:]:
        if extract_number(f, 0) != last :#and last!=1110:
            print((f"Missing file: {last}"))
            #raise RuntimeError(f"Missing file: {last}")

        last = extract_number(f, 1)



def main():
    verify_completion("complete.txt")


if __name__ == '__main__':
    main()