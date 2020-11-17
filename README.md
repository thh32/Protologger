![logo](/images/Protologger-logo.png)


# What is Protologger?

Protologger is an all-in-one genome description tool, aimed at simplifying the process of gathering the data required for writing protologues. This includes providing; taxonomic, functional and ecological insights, described in detail below;

Taxonomic placement
- 16S rRNA identity values
- 16S rRNA gene phylogenetic tree
- Taxonomic assignment based on GTDB (GTDB-Tk r89)
- Genomic tree (phylophlan3)
- ANI values (FastANI)
- POCP values

Functional analysis
- KEGG based pathway reconstruction
- CAZyme profiling
- Antibiotic resistance profiling (CARD)

Ecological analysis
- Prevalence and average relative abundance across 1,000 samples for 19 unique environments (IMNGS)
- Occurence based on MASH comparison to a database of >50,000 MAGs from thousands of samples

# What are protologues and why are they needed?
According to the latest version of the International Code on the Nomenclature of Prokaryotes (ICNP), publication of a novel taxa must include a description of that taxas' features. The format of this information is termed a 'protologue', of which some examples are provided on the example page. Whilst the ICNP are vague on the specifics of what should be included, generally protologues include; the functional features, isolation source and taxonomic placement of the species in relation to existing validly named taxa.
Therefore, Protologger provides all the necessary information for writing protologues, reducing the burden on cultivation experts for the validation of names for novel taxa.

# Dataset
Within the publication we provide the Protologger output for four distinct datasets; the [HBC](https://www.nature.com/articles/s41587-018-0009-7), the [BIO-ML collection](https://www.nature.com/articles/s41591-019-0559-3), the [Hungate1000](https://www.nature.com/articles/nbt.4110) and the [iMGMC](https://www.sciencedirect.com/science/article/pii/S2211124720301972?via%3Dihub). 

Additionally, Protologger has been applied to characterise isolates from the chicken gut (CHiBAC) which are available [here](LONK)

All datasets are downloadable [here](https://drive.google.com/file/d/1abNuXifhd2mH8txxkVhUO9MOZLcMETTb/view?usp=sharing)

# Citation
We ask that anyone who uses Protologger cites not only our publication but the list of publications below which provide tools and databases which are integral for Protologgers working;

## Tools
- Living Tree Project (https://www.sciencedirect.com/science/article/abs/pii/S072320200800060X?via%3Dihub)
- UCHIME (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3150044/pdf/btr381.pdf)
- LPSN (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3965054/pdf/gkt1111.pdf)
- GTDB-Tk (https://academic.oup.com/bioinformatics/article/36/6/1925/5626182)
- MUSCLE (https://academic.oup.com/nar/article/32/5/1792/2380623)
- PROKKA (https://academic.oup.com/bioinformatics/article/30/14/2068/2390517)
- PRODIGAL (https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-119)
- MASH (https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0997-x)
- blastN (https://www.sciencedirect.com/science/article/pii/S0022283605803602?via%3Dihub)
- DIAMOND (https://www.nature.com/articles/nmeth.3176)
- CheckM (https://genome.cshlp.org/content/25/7/1043)
- fastANI (https://www.nature.com/articles/s41467-018-07641-9)


## Databases 
- CARD (https://academic.oup.com/nar/article/48/D1/D517/5608993)
- CAZyme (https://academic.oup.com/nar/article/37/suppl_1/D233/1003505)
- KEGG (https://academic.oup.com/nar/article/27/1/29/1238108)
- IMNGS amplicon samples (https://www.nature.com/articles/srep33721)


