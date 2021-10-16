###########################
#### Section7_TF_motif ####
###########################


#### We used HOMER/v4.11 (Heinz et al. 2010) to discover DNA sequence motifs 
#### that are specifically enriched in the core promoters of key genes (-300 
#### to +50 bp of annotated TSSs). Running the findMotifs.pl function with 
#### default parameters, each test set was defined as the core promoters of 
#### either upregulated or downregulated key genes at specific cell edges 
#### (excluding sets with fewer than 5 key genes), and the background as core 
#### promoters of key genes from all edges not in the test set. Motif quality 
#### was evaluated based on a q-value, which was calculated for each motif by 
#### 100 iterations of randomizing data labels and re-running HOMER. Motifs were 
#### aligned to known motif binding sequences based on the JASPAR and internal 
#### HOMER databases with default parameters as well (Khan et al. 2018). 
#### Mapping of specific motif positions around the TSS was assessed with the 
#### HOMER function annotatePeaks.pl using the following parameters: tss mm10 -hist 10 -ghist.

#### Please email either Chengxiang Qiu (cxqiu@uw.edu) or Tony Li (tli824@uw.edu)
#### if you have any problems to perform this analysis.