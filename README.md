# Partitioning temperature sensitivity README

This file accompanies the paper entitled "Partitioning the apparent temperature sensitivity: revisiting the difference between autotrophic and heterotrophic protists" by Chen et al. submitted to The American Naturalist and introduces each file (including the source data file) used for statistical analysis and plotting.

## Authors

Bingzhang Chen, David J. S. Montagnes, Qing Wang, Hongbin Liu, Susanne Menden-Deuer

## Contact email

[bingzhang.chen\@strath.ac.uk](mailto:bingzhang.chen@strath.ac.uk){.email}

## Brief summary of the study

This study revisits the difference in thermal sensitivity between autotrophs and heterotrophs by developing an analytic method to separate the within-taxa and across-taxa temperature responses.

## Contributor

Bingzhang Chen is responsible for collecting the data and writing the R code.

## LICENSE

All the codes are covered by the MIT license and the data are covered by the CC0 license. Please see the **LICENSE** file for details.

# Metadata

## Software and packages

All the codes are written in R 4.2.0. The following R packages are used: **foreach** (version 1.5.2), **nlme** (version 3.1-157), **plyr** (version 1.8.7), **dplyr** (version 1.0.9). The  detailed package versions and environment can be found in the file **Session_Info.png**.  

## How to run the code

1.  Download the source code and data. The most updated code and data are available at <https://github.com/BingzhangChen/ActivationEnergy.git> which can be obtained either by using git clone or directly downloaded from Github.

2.  To run the main analysis (**Table 2** in the main paper), run **Table2.R** in R. This script also generates Fig. S1 and S2.

3.  To reproduce **Fig. 1** in the main paper, run **Fig1Eapp.R**.

4.  To reproduce **Fig. 2** in the main paper, run **Fig2.R**.

## Data files

1.  HProtist.csv: Raw data of growth rate \~ temperature of heterotrophic protists. The **Reference** column indicates the original data source with the detailed reference shown in **Data_source_references.docx**. The **Group** column indicates the broad functional group (Ciliates, Flagellates, Dinoflagellates, or Amoebae) of the taxa. The **Genus** column indicates the genus of the taxa. The **Species** column indicates the species name of the taxa. The **ID** column is a tag to indicate an independent temperature-growth experiment; there can be multiple **ID**s for a same species. **Habitat** indicates whether it is a freshwater or marine species. **Volume** indicates the cell volume (unit: $\mu m^3$) of the taxa. Most of the times we used a same value of cell volume for a given taxon cultured at different temperatures unless different values of cell volume had been reported by the original study. The **Temperature** column indicates the experimental temperature (degree in Celsius). The **Growth** column indicates the measured *per capita* growth rate (unit: $d^{-1}$).

2.  Merged_PHY.Rdata: merged data of phytoplankton growth rate \~ temperature generated by **Phy_merge.R**. This dataset was expanded from Chen and Laws (2017) and Kremer et al. (2017) with some additions. If this .Rdata file is loaded in R, a dataframe called **newdat** can be found with the column headings the same as those in **HProtist.csv** except that we do not provide the original references which can be found in Chen and Laws (2017) and Kremer et al. (2017). The **ID** names starting with **Chen** indicate that the data were from Chen and Laws (2017) and those starting with **Kremer** indicate that the data were from Kremer et al. (2017).

3.  Kremer2017.csv: original data from Kremer et al. (2017) used for generating Merged_PHY.Rdata. It is not needed for the main analysis. The **isolate.code** column indicates the ID of an independent temperature-growth experiment. **genus** indicates the genus of the taxa. **species** and **name** indicate the species name of the taxa. The **group** indicates the functional group of the taxa. **environment** indicates whether it is a freshwater or marine species. **bv** indicates the log10 cell volume (unit: $\mu m^3$) of the taxa. The **temperature** column indicates the experimental temperature (degree in Celsius). The **r** column indicates the measured *per capita* growth rate (unit: $d^{-1}$).

4.  Chen2017.csv: phytoplankton data of growth rate \~ temperature expanded from Chen and Laws (2017) used for generating Merged_PHY.Rdata. It is not needed for the main analysis. The **Reference** column indicates the original reference. The **ID** column indicates the ID of an independent temperature-growth experiment. The columns of **Phylum**, **Class**, **Order**, **Family**, **Genus**, **Species** indicate the hierarchical taxonomic classification of the taxa. **Habitat** indicates whether it is a freshwater or marine species.**Temp** indicates the experimental temperature (degree in Celsius). **Volume** indicates the log10 cell volume (unit: $\mu m^3$) of the taxa. **Growth** indicates the measured *per capita* growth rate (unit: $d^{-1}$). **Lat** and **Lon** indicate the latitude and longitude of the location where the taxa were isolated.

5.  Insects.csv: insect data of per capita growth rate \~ temperature from Rezende and Bozinovic (2019). The **Species** column indicates the species name of the taxa. The **ID** column indicates the ID of an independent temperature-growth experiment. **Temperature** indicates the experimental temperature (degree in Celsius). **Performance** indicates the measured *per capita* growth rate (unit: $d^{-1}$).

6.  Smith2019Bac.csv: bacterial data of per capita growth rate \~ temperature from Smith et al. (2019). The **StandardisedTraitName** column indicates Standardised name of trait after conversion to SI units. The **StandardisedTraitValue** indicates Trait value after standardisation to SI units. **StandardisedTraitUnit** indicates Units of the standardised trait. **AmbientTemp** indicates Temperature at trait recording. **AmbientTempUnit** indicates Units the temperature was measured in. **Latitude** indicates Latitude of the location. **Longitude** indicates Longitude of the location. **ConKingdom** indicates the taxon's kingdom. **ConPhylum** indicates the taxon's phylum. **ConClass** indicates the taxon's class. **ConOrder** indicates the taxon's order. **ConFamily** indicates the taxon's family. **ConGenus** indicates the taxon's genus. **ConSpecies** indicates the taxon's species. **ConTrophic** indicates the taxon's trophic mode (autotrophic vs. heterotrophic). **ConThermy** indicates whether the taxon is ecotothermic or endothermic. **ConTemp** indicates the temperature of the taxon at trait measurement. **ConTempUnit** indicates units of the temperature recording. **Citation** indicates reference for the original trait recording and **DOI** indicates its DOI.

7.  Data_source_references.docx: the word document containing the references lists of the data sources of the autotrophic and heterotrophic protists.

## R scripts

1.  prep_data.R: R script used for reading the source data files (**HProtist.csv** and **Merged_PHY.Rdata**) and choose appropriate data that fit the selection criteria (see the main manuscript for details). It also calculates the relevant information such as mean Boltzmann temperature for each taxon.

2.  Phy_merge.R: R script to merge two phytoplankton datasets from Chen & Laws (2017) and Kremer et al. (2017) and remove the duplicates.

3.  Fig1Eapp.R: R script to generate Fig. 1 in the paper.

4.  Fig2.R: R script to generate Fig. 2 in the paper.

5.  Table2.R: R script to generate Table 2 in the main paper. It also contains the code to generate Table S1-S3, Fig. S1 and S2 in the supplements.

## Other supplemental materials

1.  AmNat60700Suppl.pdf: Supplemental files including Supplement 1 (Derivations of Eq. 1 and Eq. 2 in the main text and additional analysis results of autotrophic and heterotrophic prokaryotes as well as insects (Table S1 and S2)) and Supplement 2 (Estimations of $E_{app}$ by incorporating cell size (Table S3)).

2.  TableS4_Autotrophic_protists_summary.pdf: Summaries of autotrophic protists used in this study (Table S4).

3.  TableS5_Heterotrophic_protists_summary.pdf: Summaries of heterotrophic protists used in this study (Table S5).

4.  AProtist_bytaxon_linear.pdf: Individual fits of growth rate \~ temperature of each taxon of autotrophic protists using Ordinary Least-Square linear regression models.

5.  HProtist_bytaxon_linear.pdf: Individual fits of growth rate \~ temperature of each taxon of heterotrophic protists using Ordinary Least-Square linear regression models.

6. Session_Info.png: the output of the command **sessionInfo()** showing the detailed software package versions as well as dependencies.

## Funding

This study has been supported by Southern Marine Science and Engineering Guangdong Laboratory (Guangzhou) (SMSEGL20SC02), a FILAMO mobility grant provided by the University of Bergen, Norway, a Leverhulme Trust Research Project Grant (RPG-2020-389), and the National Science Foundation, Biological Oceanography OCE-1736635.

## References

1.  Chen, B., and E. A. Laws. 2017. Is there a difference of temperature sensitivity between marine phytoplankton and heterotrophs? Limnology and Oceanography 62:806--817.

2.  Kremer, C. T., M. K. Thomas, and E. Litchman. 2017. Temperature- and size-scaling of phytoplankton population growth rates: Reconciling the Eppley curve and the metabolic theory of ecology. Limnology and Oceanography 62:1658--1670.

3.  Rezende, E. L., and F. Bozinovic. 2019. Thermal performance across levels of biological organization. Philosophical Transactions of the Royal Society B: Biological Sciences 374:20180549.

4.  Smith, T. P., T. J. H. Thomas, B. García-Carreras, S. Sal, G. Yvon-Durocher, T. Bell, and S. Pawar. 2019. Community-level respiration of prokaryotic microbes may rise with global warming. Nature Communications 10:5124.