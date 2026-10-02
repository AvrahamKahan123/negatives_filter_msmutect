Write all code here: do not change in other directories without asking me
I can provide a full file list of directories if neccesary
The idea is to run MSMuTect over the all our TCGA files using from_file, using the currently existing histograms
The run environment is the university cluster. Jobs are submitted with condor. I can only write to directory /storage/bfe_maruvka/avrahamk and under that
The ideal way to run most of these jobs is queue_from with an inputs file

1. All files are in ls /storage/bfe_maruvka/gaiafr/Research_project/WGS_SNVs_indels_analysis_project/MS_Analysis/Google_Cloud_MSMuTect_final_output_new_May2025_run/*.gz
with names like TCGA-W5-AA2K.full.mut.tsv.gz. Save only lines with a call of "M" (in the call field) to a new file with the same sample name in /storage/bfe_maruvka/avrahamk/run_over_whole_tcga_again/filtered. This should be one script+condor submit file
2. Add germline annotation with code in converting_from_old_format_to_new. Again, condor queue from. Save results to /storage/bfe_maruvka/avrahamk/run_over_whole_tcga_again/filtered_wgermline 
3. Run msmutect from file on them. Save results to /storage/bfe_maruvka/avrahamk/run_over_whole_tcga_again/msmutect_results