We want to create a "mixed histograms" file for each of the following 5 samples:
TCGA-A6-5661
TCGA-AJ-A3BH
TCGA-AP-A05N
TCGA-FI-A2D4
TCGA-OR-A5LB

For the same filtered files, write a bash+condor script to run the mixed vs the normal from_file
Normal have names like: /storage/bfe_maruvka/avrahamk/positives_simulation/mixed_histograms/TCGA-OR-A5LB_1x.normal.hist.tsv
pesudo-tumors have names like: /storage/bfe_maruvka/avrahamk/positives_simulation/mixed_histograms/TCGA-A6-5661_1x_0.15purity.hist.tsv
for all values from 0.1 to 1, in jumps of 0.05

Use a queue from file for all these values to run them all in parralel