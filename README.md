# Data Processing and Quality Control
This project contains reproduced results of Marisa.et.al. as the submission of my Master's final collaborative project(role:programmer) for class BF528 at Boston University.


This analysis is focussed only on reproducing the results from the comparison of the C3 and C4 tumor subtypes from the studies of **[Marisa et al](https://pubmed.ncbi.nlm.nih.gov/23700391/)**. The research was conducted in a two-phase design, where an initial set of “discovery” samples was used to identify patterns among the samples, and a separate set of “validation” samples was used to test if the results from the discovery set were robust. For this analysis, the discovery and validation set samples have been combined into a single dataset that will be used. There are 134 samples in total.

This project contains the R Script that shows data preprocesisng and quality control by : 
  - 1. Normalization of all of microarrays together using the Robust Multiarray Averaging (RMA) algorithm 
  - 2. Computation of standard quality control metrics on the normalized data (Batch effect correction using ComBat)
  - 3. Visualization of the distribution of samples using Principal Component Analysis (PCA).
