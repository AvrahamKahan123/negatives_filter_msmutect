All code should be written in python files in this directory. It may import from other files as needed
All results and downloaded data should be saved in local subdirectories. Care should be taken to not completely fill up the disk

The code in this directory is meant to verify the current version of MSMuTect marks as mutations a subset of the mutations
that the old version did. Here are the steps to do so:
Use the following samples to verify:
TCGA-A6-5661
TCGA-AJ-A3BH
TCGA-AP-A05N
TCGA-FI-A2D4
TCGA-OR-A5LB

The current filtered samples to do so are at 
/storage/bfe_maruvka/gaiafr/Research_project/WGS_SNVs_indels_analysis_project/MS_Analysis/Google_Cloud_MSMuTect_final_output_new_May2025_run/MSMuTect_called_mut_filt_fixed
with the following file pattern:
TCGA-ZU-A8S4.called.filt.mut.tsv.gz

The current full run files are at 
/storage/bfe_maruvka/gaiafr/Research_project/WGS_SNVs_indels_analysis_project/MS_Analysis/Google_Cloud_MSMuTect_final_output_new_May2025_run
with the files having the following file pattern: TCGA-ZU-A8S4.full.mut.tsv.gz

step 1: write a script to scp the data to local
step 2: run all the filtered files through yossi_filter() in ../../results_postprocessing/analyze_full_tcga_mut_set.py
step 3: Run MSMuTect with the --from_file mode on the outputs from the full output file. An example of how to do this can be seen at
../../run_all_GIB_files_through_new_msmutect/run_msmutect.py




