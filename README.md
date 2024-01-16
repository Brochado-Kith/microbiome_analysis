# Microbiome analysis
Thank you for reach these scripts, hope you find them usefull.

The following codes are focused on microbiome analysis with R.
Alpha, Beta and Relative abundance differences (RA).

These scripts are made with shinny code that allows interactivity.
You will be able to choose which groups will be on the analysis, percentaje of prevalence of Clusters, etc.

First of all, it is important to you to know that the analysis in the scripts starts from a phyloseq object.
If you do not have a phyloseq object, you could start from [this script](https://github.com/Brochado-Kith/microbiome_analysis/blob/main/Scripts/creating_phyloseq_object.rmd).

Phyloseq objects are a compound of tables with information interconected.

> otu_table: Table including the counts of our clusters (OTUs, ASV...) of our sequencing esperiment. This table will have samples as columns and clusters as rows. [example otu table]()

> sample_data: Table with the clinical data of the samples included in the study. This table will have the samples as rows and the clinical variables as columns. [example sample data]()

> tax_table: Table where taxonomic information is stored. This table will have clusters as rows and taxonomic ranks as columns. [example tax_table]()

optional

> phy_tree: This object will have a phylogenetic tree that will be usefull for subsequent analysis. example [phy tree]()

## Alpha diversity

Once you have created the phyloseq object you can start the analysis.
First we will calculate the alpha diversity in our data. Alpha diversity denominates the diviersity within each sample. 

There are different estimators of alpha diversity:

> Chao1: Estimates the diversity by using the total number of species(N) and the number of singletons(S) (species observed only once) and doubletons(D) (species observed two times). The formula will be **Chao1=N+S^2/(2xD)**.

> Shannon: Estimates the heterogeneicity of the sample taking into account the number of species and the abundance of them. It is explained as the probability of taking a certain species above others. If one specie is predominant, the probability will be higher and the diversity will be lower.

![f837a019109e054043577a3396c2cca5a2c3ae50](https://github.com/Brochado-Kith/microbiome_analysis/assets/135698696/efd0d792-ace9-457d-b06e-b1b94d67cc74)

> Simpson:

