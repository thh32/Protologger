![logo](/images/Protologger-logo.png)


# What is Protologger?

Protologger is an all-in-one genome description tool, aimed at simplifying the process of gathering the data required for writing protologues. This includes providing; taxonomic, functional and ecological insights, described in detail below;

### Taxonomic placement
- 16S rRNA identity values
- 16S rRNA gene phylogenetic tree
- Taxonomic assignment based on GTDB (GTDB-Tk r89)
- Genomic tree (phylophlan3)
- ANI values (FastANI)
- POCP values

### Functional analysis
- KEGG based pathway reconstruction
- CAZyme profiling
- Antibiotic resistance profiling (CARD)

### Ecological analysis
- Prevalence and average relative abundance across 1,000 samples for 19 unique environments (IMNGS)
- Occurence based on MASH comparison to a database of >50,000 MAGs from thousands of samples


# Reference
Hitch, T.C.A., Riedel, T., Oren, A. et al. Automated analysis of genomic sequences facilitates high-throughput and comprehensive description of bacteria. ISME Communications. 1, 16 (2021). https://doi.org/10.1038/s43705-021-00017-z


# What are protologues and why are they needed?
According to the latest version of the International Code on the Nomenclature of Prokaryotes (ICNP), publication of a novel taxa must include a description of that taxas' features. The format of this information is termed a 'protologue', of which some examples are provided on the example page. Whilst the ICNP are vague on the specifics of what should be included, generally protologues include; the functional features, isolation source and taxonomic placement of the species in relation to existing validly named taxa.
Therefore, Protologger provides all the necessary information for writing protologues, reducing the burden on cultivation experts for the validation of names for novel taxa.

# Installation

## Web-server

If you have single isolates that you wish to study, why not use our Galaxy web-server which, freely available at; http://protologger.de/

Included on our website is an implementation of GAN (The Great Automatic Nomenclator) (https://github.com/telatin/gan) which accepts Protologger output to generate ecological and functionally informed names.


## Conda environment

Protologger has been re-written in python3 for easy conda installation.

There are four steps that are required to get Protologger working;

1. Create a python3 environment using the command; `conda create -n protologger python=3.7 prokka`
2. Install Protologger into this environment; `conda install -c thitch protologger`
3. Once installed, the databases must be downloaded using the following command; `setup-protologger.sh`
4. Make sure you have Usearch installed (version 5.2.32 is the tested version) and is in your $PATH

When finished (~5 hours depending on your internet speed), Protologger will be ready to run. Additionally, the command `protologger-update.sh` can be run to download the latest validation list which is updated monthly.

# Usage

Protologger.py requires three inputs to be run on the commandline, detailed below.

| Input flag | File type | Description                                                                                              |
| ------------- | ---------| -------------------------------------------------------------------------------------------------------- |
| -r     | Nucleotide FASTA file | Provide a file containing the 16S rRNA gene sequence for your species of interest           |
| -g    | Nucleotide FASTA file | Provide the genome file of your species of interest                                                    |
| -p       | String | Provide the name of your project which will be used to name your input and output folders                                               |
| -q       | NA | This option activates 'quick' mode which ignores both GTDB-Tk and PhyloPhlAn analysis, meaning Protologger can be run on a desktop PC                                               |



# Datasets
Within the publication we provide the Protologger output for four distinct datasets; the [HBC](https://www.nature.com/articles/s41587-018-0009-7), the [BIO-ML collection](https://www.nature.com/articles/s41591-019-0559-3), the [Hungate1000](https://www.nature.com/articles/nbt.4110) and the [iMGMC](https://www.sciencedirect.com/science/article/pii/S2211124720301972?via%3Dihub). 

The Protologger output from all four datasets are downloadable [here](https://drive.google.com/file/d/1abNuXifhd2mH8txxkVhUO9MOZLcMETTb/view?usp=sharing).


### ChiBAC - Chicken Bacterial Collection
Protologger has been applied to characterise isolates from the chicken gut (CHiBAC), for which the entire Protologger outputs are available [here](https://drive.google.com/file/d/19icEV9xH6PBiW1Mekd2duG6IJP3jf9HY/view?usp=sharing)


## MAG datasets

We always aim to expand the MAG database used within Protologger. If there is an additional dataset you wish included, please contact us at; admin@protologger.de

The list of currently included datasets (numbers represent Bacterial and Archaeal numbers) is as follows;

| Publication  | Number of MAGs | Description                                                                                              |
| ------------- | ---------| -------------------------------------------------------------------------------------------------------- |
| Parks et al (2020)     | 3,397 | Generic           |
| Woodcraft et al (2018)    | 568 | Permafrost                                                   |
| Anantharaman et al (2016)       | 303 | Groundwater                                               |
| Crits-Cristoph et al (2018)       | 225 | Soil                                               |
| Dombrowski et al (2018)       | 36 | Hydro-thermal sediment                                               |
| Tully et al (2018)       | 339 | Ocean                                               |
| Lesker et al (2020)       | 831 | Mouse gut                                            |
| Wylensek et al (2020)       | 589 | Pig gut                                               |
| Stewart et al (2018)       | 488 | Bovine rumen                                               |
| Almeida et al (2019)       | 39,891 | Human gut                                               |
| Pasolli et al (2019)       | 225 | Human stool/vagina/skin/oral cavity                                               |
| Manara et al (2019)       | 1,008 | soil                                               |


The list of datasets currently undergoing integration are;

| Publication  | Number of MAGs | Description                                                                                              |
| ------------- | ---------| -------------------------------------------------------------------------------------------------------- |
| Wilkinson et al (2020)     | 3,397 | African Boran rumen           |
| Chen et al (2021)     | 1,358 | Pig gut           |
| Albanese et al (2021)     | 263 | Cryptoendolithic community in Antarctica           |
| Levin et al (2021)     | 1,209 | Gut of 180  wild species           |
| Jegousse et al (2021)     | 219 | Icelandic marine water           |
| Robbins et al (2021)     | ~1200 | Coral sponge           |
| Wibowo et al (2021)     | 498 | Ancient human gut            |
| Collins et al (2021)     | 111 | Deep sea fish gut            |
| Nayfach et al (2021)     | 52,515 | Earth associated microbiota           |
| Becraft et al (2021)     | 126 | Deep subsurface            |
| Krüger et al (2019)     | 3101 | Algae blooms             |
| Lavrinienko et al (2020)     | 254 | Bank vole gut            |


# Citation
Hitch, T.C.A., Riedel, T., Oren, A. et al. Automated analysis of genomic sequences facilitates high-throughput and comprehensive description of bacteria. ISME COMMUN. 1, 16 (2021). https://doi.org/10.1038/s43705-021-00017-z

We ask that anyone who uses Protologger cites not only our publication but the list of publications below which provide tools and databases which are integral for Protologgers working;

### Tools
- [Living Tree Project](https://www.sciencedirect.com/science/article/abs/pii/S072320200800060X?via%3Dihub)
- [UCHIME](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3150044/pdf/btr381.pdf)
- [LPSN](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3965054/pdf/gkt1111.pdf)
- [GTDB-Tk](https://academic.oup.com/bioinformatics/article/36/6/1925/5626182)
- [MUSCLE](https://academic.oup.com/nar/article/32/5/1792/2380623)
- [PROKKA](https://academic.oup.com/bioinformatics/article/30/14/2068/2390517)
- [PRODIGAL](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-119)
- [MASH](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0997-x)
- [blastN](https://www.sciencedirect.com/science/article/pii/S0022283605803602?via%3Dihub)
- [DIAMOND](https://www.nature.com/articles/nmeth.3176)
- [CheckM](https://genome.cshlp.org/content/25/7/1043)
- [FastANI](https://www.nature.com/articles/s41467-018-07641-9)
- [PhyloPhlAn3](https://www.nature.com/articles/s41467-020-16366-7)
- [FastTree](https://academic.oup.com/mbe/article/26/7/1641/1128976)

### Databases 
- [CARD](https://academic.oup.com/nar/article/48/D1/D517/5608993)
- [CAZyme](https://academic.oup.com/nar/article/37/suppl_1/D233/1003505)
- [KEGG](https://academic.oup.com/nar/article/27/1/29/1238108)
- [IMNGS amplicon samples](https://www.nature.com/articles/srep33721)


