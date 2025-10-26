import os


original_dir = "/storage/bfe_maruvka/gaiafr/Research_project/WGS_SNVs_indels_analysis_project/MS_Analysis/Google_Cloud_MSMuTect_final_output_new_May2025_run/MSMuTect_called_mut_filt/"
with open("relevant_tsv.txt", "r") as tsvs_list:
    lines = tsvs_list.readlines()
# print(lines)
for l in lines:
    tcga_id = l[23:35]
    suffix = l[36:]
    purity = suffix[:suffix.find("p")]
    input_tumor_file = l.rstrip()
    input_normal_file = tcga_id+".called.filt.hist.tsv"

    # input_normal_file = os.path.join(original_dir, tcga_id+".called.filt.mut.tsv")
    output_file = f"ouput_files/{tcga_id}_{purity}"
    print(f"{input_tumor_file} {input_normal_file} {output_file}")