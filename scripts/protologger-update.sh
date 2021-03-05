echo 'Downloading latest list of valid species names...'
wget --no-check-certificate --content-disposition https://raw.githubusercontent.com/thh32/Protologger/master/DSMZ-latest.tab -O ${PROTOLOGGER_DATA_DIR}DSMZ-valid-list/DSMZ-August-2019.tab


echo 'File updated and now contains the following valid species names; '
${PROTOLOGGER_DATA_DIR}DSMZ-valid-list/DSMZ-August-2019.tab | wc -l

echo 'Closing.' 
