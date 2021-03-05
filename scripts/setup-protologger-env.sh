#! /bin/bash

echo 'Starting....'
echo $PWD
git clone https://github.com/thh32/Protologger
export PATH="$PWD/protologger/scripts:$PATH"
mkdir -p $PWD/protologger/bin
export PROTOLOGGER_DATA_DIR=$PWD/protologger/bin/
echo 'Protologger data is stored at; ' $PROTOLOGGER_DATA_DIR


# Create conda environment
conda create -c bioconda -n protologger python=2.7 gtdbtk muscle hmmer prokka mash diamond checkm-genome clustalw htseq


# Download GTDB database
source ~/anaconda3/etc/profile.d/conda.sh
source ~/miniconda3/etc/profile.d/conda.sh
conda activate protologger
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/release89/89.0/gtdbtk_r89_data.tar.gz -P ${GTDBTK_DATA_PATH}
tar xvzf ${GTDBTK_DATA_PATH}/gtdbtk_r89_data.tar.gz -C ${GTDBTK_DATA_PATH} --strip 1
rm ${GTDBTK_DATA_PATH}/gtdbtk_r89_data.tar.gz

# Download and install phylophlan
git clone https://github.com/biobakery/phylophlan
cd phylophlan
sudo python3.6 setup.py install
python3.6 phylophlan/phylophlan.py 
export phylophlan=$PWD/phylophlan/phylophlan.py 

cd ../
