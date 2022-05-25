# Partitioning temperature sensitivity README

This file accompanies the paper entitled "Partitioning the apparent temperature sensitivity: revisiting the difference between autotrophic and heterotrophic protists" by Chen et al. submitted to The American Naturalist and introduces each file (including the source data file) used for statistical analysis and plotting. 

# Metadata for essential files

## Data files
1. Microzoo_23Mar2020.csv: Original data of growth rate ~ temperature of microzooplankton.

2. NEWphy24Mar2020.Rdata: merged  data of phytoplankton growth rate ~ temperature. This dataset was expanded from Chen & Laws (2017) and Kremer et al. (2017).

3. lno10523-sup-0008-suppinfo8.csv: original data from Kremer et al. (2017).

4. Phyto_19Mar2020.csv: phytoplankton data of growth rate ~ temperature expanded from Chen and Laws (2017). 

5. Data_source_references.docx: the word document containing the references lists of the data sources.

## R scripts

1. prep_data.R: R script used for reading the source data files ('Microzoo_23Mar2020.csv' & 'NEWphy24Mar2020.Rdata') and choose appropriate data that fit the selection criteria (see text for details). It also calculates the relevant information such as mean Boltzmann temperature for each taxon.

2. Phy_merge.R: R script to merge two phytoplankton datasets from Chen & Laws (2017) and Kremer et al. (2017) and remove the duplicates.

3. Fig1Eapp.R: R script to generate Fig. 1 in the paper. 

4. Fig2.R: R script to generate Fig. 2 in the paper.

5. Table2.R: R script to generate Table 2 in the main paper. It also contains the code to generate Table S1-S3, Fig. S1 and S2 in the supplements.

## References

1. Chen, B., and E. A. Laws. 2017. Is there a difference of temperature sensitivity between marine phytoplankton and heterotrophs? Limnology and Oceanography 62:806–817.

2. Kremer, C. T., M. K. Thomas, and E. Litchman. 2017. Temperature- and size-scaling of phytoplankton population growth rates: Reconciling the Eppley curve and the metabolic theory of ecology. Limnology and Oceanography 62:1658–1670.

