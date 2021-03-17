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

# What are protologues and why are they needed?
According to the latest version of the International Code on the Nomenclature of Prokaryotes (ICNP), publication of a novel taxa must include a description of that taxas' features. The format of this information is termed a 'protologue', of which some examples are provided on the example page. Whilst the ICNP are vague on the specifics of what should be included, generally protologues include; the functional features, isolation source and taxonomic placement of the species in relation to existing validly named taxa.
Therefore, Protologger provides all the necessary information for writing protologues, reducing the burden on cultivation experts for the validation of names for novel taxa.

# Installation

## Web-server

If you have single isolates that you wish to study, why not use our Galaxy web-server which, freely available at; http://protologger.de/

Included on our website is an implementation of GAN (The Great Automatic Nomenclator) (https://github.com/telatin/gan) which accepts Protologger output to generate ecological and functionally informed names.


## Conda environment

Due to the use of both Python 2.7 and Python 3 within Protologger a Conda package has not been made, instead we offer an alternative solution via a conda environment.

For this you must first have both `conda` and `git` installed.

Download the setup-protologger-env.sh file, [here](https://github.com/thh32/Protologger/blob/master/scripts/setup-protologger-env.sh).

Move this script to the folder in which you wish to host the Protologger databases (~50Gb).

Run the command `. setup-protologger-env.sh` and wait, the installation process takes ~2-3 hours, depending on your internet connection as multiple databases must be downloaded.

This script creates a conda environment called `protologger` that has all the required tools pre-installed, apart from Usearch (v5.2.32), which must be installed manually. Make sure that the usearch executable is in your $PATH so it can be called using `usearch`.

<b> IT IS INSTALLED!</b> You now have access to your own installation of Protologger, but be aware, Protologger is very RAM hungry and our own version is run on a machine with 500Gb RAM and 52 cores.

You can update the databases used by Protologger at any point using the command; `protologger-update.sh`.

Protologger currently utilises r89 of the genome taxonomy database (GTDB) and the Living Tree Project (LTP v132) as base databases and all code is designed around these. Further updates will be implemented for GTDB r95 soon.

# Usage

Protologger requires three inputs to be run on the commandline, detailed below.

| Option command| File type | Description                                                                                              |
| ------------- | -------------------------------------------------------------------------------------------------------- |
| -s     | Nucleotide FASTA file | Provide a file containing the 16S rRNA gene sequence for your species of interest           |
| -g    | Nucleotide FASTA file | Provide the genome file of your species of interest                                                    |
| -p       | String | Provide the name of your project which will be used to name your input and output folders                                               |



# Datasets
Within the publication we provide the Protologger output for four distinct datasets; the [HBC](https://www.nature.com/articles/s41587-018-0009-7), the [BIO-ML collection](https://www.nature.com/articles/s41591-019-0559-3), the [Hungate1000](https://www.nature.com/articles/nbt.4110) and the [iMGMC](https://www.sciencedirect.com/science/article/pii/S2211124720301972?via%3Dihub). 

The Protologger output from all four datasets are downloadable [here](https://drive.google.com/file/d/1abNuXifhd2mH8txxkVhUO9MOZLcMETTb/view?usp=sharing).

### ChiBAC - Chicken Bacterial Collection
Protologger has been applied to characterise isolates from the chicken gut (CHiBAC), for which the entire Protologger outputs are available [here](https://drive.google.com/file/d/19icEV9xH6PBiW1Mekd2duG6IJP3jf9HY/view?usp=sharing)


# Citation
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


