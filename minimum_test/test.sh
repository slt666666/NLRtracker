# minimum test

echo "normal run"
bash ../NLR_extractor.sh -s sample_data/sample.fasta -o test_out

echo "no -s option"
bash ../NLR_extractor.sh -o test_out2

echo "-s option file doesn't exist"
bash ../NLR_extractor.sh -s no_path.fasta -o test_out3

echo "no -o option"
bash ../NLR_extractor.sh -s sample_data/sample.fasta

echo "check -i option"
bash ../NLR_extractor.sh -s sample_data/sample.fasta -o test_out4 -i sample_data/test_interpro.gff3

echo "check -f option"
bash ../NLR_extractor.sh -s sample_data/sample.fasta -o test_out5 -f sample_data/test_fimo.gff

echo "check -i & -f option"
bash ../NLR_extractor.sh -s sample_data/sample.fasta -o test_out6 -i sample_data/test_interpro.gff3 -f sample_data/test_fimo.gff

echo "check -c set"
bash ../NLR_extractor.sh -s sample_data/sample_cds.fasta -o test_out7 -t n

echo "check -c set"
bash ../NLR_extractor.sh -s sample_data/sample.fasta -o test_out8 -c 10

echo "check -m set"
bash ../NLR_extractor.sh -s sample_data/sample.fasta -o test_out9 -m sample_data/test.xml

echo "check -d set"
bash ../NLR_extractor.sh -s sample_data/sample.fasta -o test_out10 -m sample_data/test.list







