#! /bin/bash

echo 'Starting....'
echo $PWD
git clone https://github.com/thh32/Protologger
chmod +x $PWD/Protologger/scripts/*
export PATH="$PWD/Protologger/scripts:$PATH"

echo "1. Preparation of Protologger's internal datasets..."
wget -O "bin.tar.gz" "https://zenodo.org/record/4561557/files/bin.tar.gz"
mv bin.tar.gz $PWD/Protologger/bin.tar.gz
cd Protologger/
tar -xf bin.tar.gz 
rm bin.tar.gz
cd ../
export PROTOLOGGER_DATA_DIR=$PWD/Protologger/bin/
echo 'Protologger data is stored at; ' $PROTOLOGGER_DATA_DIR
echo "1. Internal datasets downloaded and ready."


# Create conda environment
echo "2. Creating Conda environment"
conda create -c bioconda -n protologger python=2.7 gtdbtk muscle hmmer prokka mash diamond checkm-genome clustalw htseq


# Download GTDB database
source ~/anaconda3/etc/profile.d/conda.sh
source ~/miniconda3/etc/profile.d/conda.sh
conda activate protologger

echo "2. Conda environment created, use the command; conda activate protologger"

echo "3. Downloading GTDB-Tk database"
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/release89/89.0/gtdbtk_r89_data.tar.gz -P ${GTDBTK_DATA_PATH}
tar xvzf ${GTDBTK_DATA_PATH}/gtdbtk_r89_data.tar.gz -C ${GTDBTK_DATA_PATH} --strip 1
rm ${GTDBTK_DATA_PATH}/gtdbtk_r89_data.tar.gz
echo "3. GTDB-Tk database downloaded and installed"

# Download and install phylophlan
echo "4. Downloading and installing PhyloPhlAn"
git clone https://github.com/biobakery/phylophlan
cd phylophlan
sudo python3.6 setup.py install
#python3.6 phylophlan/phylophlan.py -h
export phylophlan=$PWD/phylophlan/phylophlan.py 
cd ../
echo "4. PhyloPhlAn has been installed"



# Prep GTDB data
echo '5. Preparation of genome database'
#cp -r ${GTDBTK_DATA_PATH}/fastani/database ${PROTOLOGGER_DATA_DIR}/Genome-database
wget -O "bac120_metadata_r89.tsv" "https://data.ace.uq.edu.au/public/gtdb/data/releases/release89/89.0/bac120_metadata_r89.tsv"

mkdir -p ${PROTOLOGGER_DATA_DIR}/GTDB-TK
mv bac120_metadata_r89.tsv ${PROTOLOGGER_DATA_DIR}/GTDB-TK/
echo '5. The genome database has been correctly placed'

echo "6. Downloading CheckM database"
mkdir -p ${PROTOLOGGER_DATA_DIR}/CheckM
cd ${PROTOLOGGER_DATA_DIR}/CheckM
wget -O "checkm_data_2015_01_16.tar.gz" "https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz"
tar xvzf checkm_data_2015_01_16.tar.gz 
checkm data setRoot .
cd ../../

echo "7. Download Uniprot database for KEGG analysis"
mkdir -p ${PROTOLOGGER_DATA_DIR}/KEGG-analysis
cd ${PROTOLOGGER_DATA_DIR}/KEGG-analysis
wget -O "idmapping.dat.gz" "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz"
gzip -dc idmapping.dat.gz | awk '{if($2=="KO") print ​$1,$3}' OFS="\t" | gzip > idmapping_KO.tab.gz
rm -r idmapping.dat.gz
cd ../../
echo "7. Uniprot database integrated."


echo 'Protologger installation finalised.'
