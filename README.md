# tome_code
The scripts used for the analyses in "Systematic reconstruction of cellular trajectories during mouse embryogenesis"

Mammalian embryogenesis is characterized by rapid cellular proliferation and diversification. Within a few weeks, a single cell zygote gives rise to millions of cells expressing a panoply of molecular programs, including much of the diversity that will subsequently be present in adult tissues. Although intensively studied, a comprehensive delineation of the major cellular trajectories that comprise mammalian development in vivo remains elusive. Here we set out to integrate several single cell RNA-seq datasets (scRNA-seq) that collectively span mouse gastrulation and organogenesis. To bridge technologies, these datasets were supplemented with new, intensive profiling of around 150,000 nuclei from  a series of E8.5 embryos in 1-somite increments with an improved combinatorial indexing protocol. Overall, we define cell states at each of 19 successive stages spanning E3.5 to E13.5, heuristically connect them to their pseudo-ancestors and pseudo-descendants, and for a subset of stages, deconvolve their approximate spatial distributions. Despite being constructed through automated procedures, the resulting trajectories of mammalian embryogenesis (TOME) are largely consistent with our contemporary understanding of mammalian development. We leverage TOME to nominate transcription factors (TF) and TF motifs as key regulators of each branch point at which a new cell type emerges. Finally, we apply the same procedures to single cell datasets of zebrafish and frog embryogenesis, and nominate “cell type homologs” based on shared regulators and transcriptional states.

Preprint: https://www.biorxiv.org/content/10.1101/2021.06.08.447626v1.abstract

Data website: http://tome.gs.washington.edu

Section 1: Intensive scRNA-seq of individual embryos during early somitogenesis
The scripts under this section were used to process the sci-RNA-seq data step by step, and make the plots showed in Fig. 1

Section 2: Systematic reconstruction of the cellular trajectories of mouse embryogenesis from E3.5 to E13.5
The scripts under this section were used to create the cellular trajectories step by step, and make the plots showed in Fig. 2

Section 3: RNA-velocity analysis
The scripts under this section were used to perform the RNA-velocity analysis, and make the plots showed in Fig. 3

Section 4: Inference of the approximate spatial locations of cell states during mouse gastrulation
The scripts under this section were used to infer the spatial locations of individual cell states, and make the plots showed in Fig. 4

Section 5: Inferring the molecular histories of individual cell types
The scripts under this section were used to perform the pseudobulk analysis, and make the plots showed in Supplementary Fig. 17-19

Section 6: Systematic nomination of key transcription factors for cell type specification
The scripts under this section were used to identify key transcription factors for cell type specification, and make the plots showed in Fig. 5

Section 7: Identification of core promoter cis-regulatory motifs involved in in vivo cell type specification
The scripts under this section were used to identify core promoter cis-regulatory motifs, and make the plots showed in Supplementary Fig. 21

Section 8: Comparison of the cellular trajectories of mouse, zebrafish and frog embryogenesis
The scripts under this section were used to perform the "cross-species" comparison analysis, and make the plots showed in Fig. 6-7

help_code: I posted all the relative codes or data profiles which are necessary in this folder, and I have indicated where to import them in individual script. For example, 
Cao_NNLS_regression.R is used to perform the NNLS regression analysis as showed in Extended Data Fig. 10
help_code.R includes some general function which was used in multiple scripts (simply used by source()).
run_scrublet.py is used to perform the doublets removing step of section 1.
Mus_musculus_TF.txt is the TF list of mouse (downloaded from http://bioinfo.life.hust.edu.cn/AnimalTFDB/) used in the section 6.
Danio_rerio_TF.txt is the TF list of mouse (downloaded from http://bioinfo.life.hust.edu.cn/AnimalTFDB/) used in the section 8.
Xenopus_tropicalis_TF.txt is the TF list of mouse (downloaded from http://bioinfo.life.hust.edu.cn/AnimalTFDB/) used in the section 8.
