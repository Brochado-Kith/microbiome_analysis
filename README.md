# Microbiome analysis
Thank you for reach to these scripts, hope you find them usefull.

The following codes are focused on microbiome analysis with R.
Alpha, Beta and Relative abundance differences (RA) will be explored along this pipeline.
You can choose to run each script separatedly or to use the [app](https://github.com/Brochado-Kith/microbiome_analysis/blob/main/App/App_metagenomics.r) I made joining all the statistical scripts (you will need to create the phyloseq object).

These scripts are made with shinny code that allows interactivity.
You will be able to choose which groups will be on the analysis, percentaje of prevalence of Clusters, etc.

First of all, it is important to you to know that the analysis in the scripts starts from a phyloseq object.
If you do not have a phyloseq object, you could start from [this script](https://github.com/Brochado-Kith/microbiome_analysis/blob/main/Scripts/creating_phyloseq_object.rmd).

Phyloseq objects are a compound of tables with information interconected.

> otu_table: Table including the counts of our clusters (OTUs, ASV...) of our sequencing esperiment. This table will have samples as columns and clusters as rows. [example otu table](https://github.com/Brochado-Kith/microbiome_analysis/blob/main/data/OTU_table.txt)

> sample_data: Table with the clinical data of the samples included in the study. This table will have the samples as rows and the clinical variables as columns. [example sample data](https://github.com/Brochado-Kith/microbiome_analysis/blob/main/data/Sample_data.txt)

> tax_table: Table where taxonomic information is stored. This table will have clusters as rows and taxonomic ranks as columns. [example tax_table](https://github.com/Brochado-Kith/microbiome_analysis/blob/main/data/taxa_data.txt)

optional

> phy_tree: This object will have a phylogenetic tree that will be usefull for subsequent analysis. example [phy tree]()



NOW YOU CAN TRY THE APP INSTEAD OF THE SEPARATED SCRIPTS!!!
Download it and give it a try.



## Alpha diversity

Once you have created the phyloseq object you can start the analysis.
First we will calculate the [alpha diversity](https://github.com/Brochado-Kith/microbiome_analysis/blob/main/Scripts/Alpha_diversity.rmd) in our data. Alpha diversity denominates the diversity within each sample. 

There are different estimators of alpha diversity:

> Chao1 index: Estimates the diversity by using the total number of species(N) and the number of singletons(S) (species observed only once) and doubletons(D) (species observed two times). The formula will be **Chao1=N+S^2/(2xD)**.

> Shannon index: Estimates the heterogeneicity of the sample taking into account the number of species and the abundance of them. It is explained as the probability of taking a certain species above others. If one specie is predominant, the probability will be higher and the diversity will be lower.

![f837a019109e054043577a3396c2cca5a2c3ae50](https://github.com/Brochado-Kith/microbiome_analysis/assets/135698696/efd0d792-ace9-457d-b06e-b1b94d67cc74)

> Simpson index: Estimates the probability that when randomly choosing two species these two species correspond to the same cluster (OTU, ASV...).

![d2880c237ac6db179f9c52a1ba06512a3965af2a](https://github.com/Brochado-Kith/microbiome_analysis/assets/135698696/cdc91246-19ec-4303-9a18-85e8b9fe8510)

## Beta diversity
[Beta diversity](https://github.com/Brochado-Kith/microbiome_analysis/blob/main/Scripts/Beta_diversity.rmd) estimates the community differences within groups.
First thing we will do is to create the distance matrix. This matrix is calculated by estimating the distance between each pair of invdividuals. The next step will be to calculate the beta dispersion. At the final step, we will carry out a permanova to see the differences in beta diversity.
As with the alpha diversity stimation, we will have differents ways to calculate the distance matrix.

> Bray curtis: Bray-Curtis dissimilarity takes into account those species that are shared between samples and sum the lesser value of eacch shared cluster (Cij) and the total value of the sum of each sample (Si and Sj).

![imagen](https://github.com/Brochado-Kith/microbiome_analysis/assets/135698696/e12c9801-ddcf-4ae4-9c00-1c18c9a430cc)


> Jaccard: Jaccard index is estimated as the number of clusters in one of the samples (a). The number of clusters in the other sample (b). The number of clusters shared by both samples.

![imagen](https://github.com/Brochado-Kith/microbiome_analysis/assets/135698696/3f5ab013-a3c8-414b-9f17-6e704992a0aa)


> Aitchison: One of the particularities of the Aitchison distance is that the data have to be centered log ratio (CLR) transformed. This distance estimation corresponds to Euclidean distances but with the CLR transformation.

> Unifrac: Unifrac distances will take into account the phylogenetic tree to calculate the distance matrix. There are two tipes of unifrac distance estimation, weigthed and unweighted. Weigthed unifrac will consider the abundance of unshared clusters when estimating the distance. Unweigthed unifrac will only consider presence absence of unshared clusters.

## Relative abundance differences

In this pipeline we will perform a generalize linear model (glm) under a negative binomial distribution by using the DESeq2 package and a glm under a zero inflated negative binomial ([ZINB](https://github.com/Brochado-Kith/microbiome_analysis/blob/main/Scripts/Relative_abundance/Relative_abundance_ZINB)) distribution.



