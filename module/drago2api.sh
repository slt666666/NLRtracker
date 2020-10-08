#!/bin/bash

fileType=$(file $1)

if [[ $fileType =~ "ASCII text" ]]; then
	reader="cat"
elif [[ $fileType =~ "gzip compressed data" ]]; then
	reader="zcat"
else
	echo "Unknown format" >&2
	exit 1
fi

#preparing input
temp_file=$(mktemp)
$reader $1 | awk '{if ($0~/^>/) {gsub(/[^A-Za-z0-9_\.:\-> ]/, "_") ; print $1} else {gsub(/[\* ]/,"") ; print}}' \
| sed $'s/\r//' > $temp_file

cat $temp_file \
| curl -X POST --data-binary @- http://prgdb.org/prgdb/drago2/pipe/ \
| sed 's/\\t/\t/g;s/\\n/\n/g;s/"//g'
RETURN_CODE=$(echo ${PIPESTATUS[@]} | awk '($0~/[1-9]+/)' | wc -l)
if [[ $RETURN_CODE != 0 ]] ; then echo "Curl request failed" >&2 ; exit 1 ; fi
