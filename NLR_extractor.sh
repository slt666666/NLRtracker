#!/bin/bash

# log
LOG_DATE=`date '+%Y-%m-%d-%H:%M:%S'`
exec 1> >(tee -a log/${LOG_DATE}_out.log)
exec 2> >(tee -a log/${LOG_DATE}_err.log)

# help
function usage {
  cat <<EOM
Usage: $(basename "$0") [OPTION]...

  -h               Display help

  (required)
  -a Filepath      File path to amino acid seqence file (.fasta)
  -o String        Directory name to save output 
  
  (optional)
  -i Filepath      Result of Interproscan (.tsv)
  -d Filepath      Result of DRAGO2 API (.txt)
EOM
  exit 2
}

# check options
while getopts ":a:i:d:o:h" optKey; do
  case "$optKey" in
    a)
      echo -e "\n-------------------------------input--------------------------------";
	  if [ -f ${OPTARG} ]; then
	  	echo "Fasta file             = ${OPTARG}"
	  	fasta=${OPTARG}
	  	FLG_A=1
	  fi
      ;;
    i)
      if [ -f ${OPTARG} ]; then
	  	echo "result of interproscan = ${OPTARG}"
	  	interpro_result=${OPTARG}
	  	FLG_I=1
	  else
	  	echo "${OPTARG} does not exits. Run interproscan in this pipeline."
	  fi
	  ;;
	d)
	  if [ -f ${OPTARG} ]; then
	  	echo "result of DRAGO2       = ${OPTARG}"
	  	DRAGO2_result=${OPTARG}
	  	FLG_D=1
	  else
	  	echo "${OPTARG} does not exits. Run DORAGO2 in this pipeline."
	  fi
	  ;;
    o)
      FLG_O=1
      echo "output directory       = ${OPTARG}"
      outdir=${OPTARG}
      echo -e "-------------------------------input--------------------------------\n";
      ;;
    '-h'|'--help'|* )
      usage
      ;;
  esac
done

# check fasta file
if [ -z $FLG_A ]; then
  echo -e "$(basename $0):「-a」option is required\n"
  usage
  exit 1
fi

# check header
if [ -z $FLG_O ]; then
  echo -e "$(basename $0):「-o」option is required\n"
  usage
  exit 1
fi

# Main pipeline
mkdir $outdir

# 1. Interproscan
if [ -z $FLG_I ]; then
  echo -e "\nRun Interproscan"
  interproscan.sh -i $fasta -f gff3 -o "${outdir}/interpro_result.gff"
  interpro_result="${outdir}/interpro_result.gff"
    
else
  echo "Pass Interproscan (Use $interpro_result as output of Interproscan)"
fi
  
  
# 2. NLR_extractor.R
if [ -f $interpro_result ]; then
  echo -e "\nRun NLR_extractor"
  Rscript ./module/NLR_extractor.R $interpro_result $outdir

else
  echo "Interproscan file error."
  exit 1
fi


# 3. DRAGO2
if [ -z $FLG_D ]; then
  echo -e "\nRun DRAGO2"
  ./module/drago2api.sh $fasta > ${outdir}/DRAGO2_out.txt
  DRAGO2_result="${outdir}/DRAGO2_out.txt"

else
  echo "Pass DRAGO2 (Use $DRAGO2_result as output of DRAGO2)"
fi


# extract DRAGO2 sepcific ids & domain check
if [ -f $DRAGO2_result ]; then
  echo -e "\nextract DRAGO2 sepcific ids"
  Rscript ./module/DRAGO2_extract.R $DRAGO2_result $outdir $interpro_result

else
  echo "DRAGO2 extract error."
  exit 1
fi

# check output
if [ -f "$outdir/DRAGO2_specific_Structure.tsv" ]; then
  echo "Finish NLR_extractor!"
fi