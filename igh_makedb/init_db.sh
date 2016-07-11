#!/bin/bash

DB=$1
WDIR=$(dirname "${2}")

# Create SQL query
QUERY=$WDIR/temp.sql

echo "USE $1;" > $QUERY
echo >> $QUERY

TABLES=( patients libs clusterings sequences sequence_strings quality_scores lineages )

for TABLE in "${TABLES[@]}"
do
    echo "LOAD DATA LOCAL INFILE '$WDIR/$TABLE.txt.temp' REPLACE INTO TABLE $TABLE" >> $QUERY
    echo "FIELDS TERMINATED BY '\t' ESCAPED BY '' LINES TERMINATED BY '\n';" >> $QUERY
    echo >> $QUERY
done

