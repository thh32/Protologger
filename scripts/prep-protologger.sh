echo "1. Preparation of Protologger's internal datasets..."
wget -O "bin.tar.gz" "https://zenodo.org/record/4561557/files/bin.tar.gz"
mv bin.tar.gz $PROTOLOGGER_DATA_DIR/bin.tar.gz
tar -xf ${PROTOLOGGER_DATA_DIR}bin.tar.gz
rm ${PROTOLOGGER_DATA_DIR}bin.tar.gz
echo "1. Internal datasets moved."

# Prep GTDB data
echo '2. Preparation of GTDB-Tk dataset...'
cp $GTDBTK_DATA_PATH/fastani/database $PROTOLOGGER_DATA_DIR/Genome-database
wget -O "bac120_metadata_r89.tsv" "https://data.ace.uq.edu.au/public/gtdb/data/releases/release89/89.0/bac120_metadata_r89.tsv"

mkdir -p $PROTOLOGGER_DATA_DIRbin/GTDB-TK
mv bac120_metadata_r89.tsv $PROTOLOGGER_DATA_DIR/GTDB-Tk/
echo '2. GTDB-Tk dependancies moved.'
