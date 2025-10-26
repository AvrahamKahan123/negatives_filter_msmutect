from typing import Dict

def num_motif_repeats():
    return 6

def fill_with_default(lst: list, default_val: object, desired_length: int):
    if desired_length < len(lst):
        return lst[:desired_length]
    else:
        return lst + [default_val for _ in range(desired_length - len(lst))]


def dict_to_csv_representation(repeats_dict: Dict[int, int]) -> str:
    num_repeats = num_motif_repeats()
    sorted_repeat_lengths = sorted(repeats_dict, key=repeats_dict.get, reverse=True)
    sorted_repeat_lengths_str_form = [str(repeat) for repeat in sorted_repeat_lengths]
    repeat_lengths_final = fill_with_default(sorted_repeat_lengths_str_form, "NA", num_repeats)
    sorted_repeat_supports_str_form = [str(repeats_dict[repeat]) for repeat in sorted_repeat_lengths]
    repeat_lengths_supports = fill_with_default(sorted_repeat_supports_str_form, "NA", num_repeats)
    return ",".join(repeat_lengths_final) +"," + ",".join(repeat_lengths_supports)

def return_na_if_none_else_return_input(x):
    if x is None:
        return "NA"
    else:
        return x


if __name__ == '__main__':
    print(dict_to_csv_representation({5: 3, 7:2, 1:10, 32: 2}))
