# Bennett / Panama
esearch -db sra -query PRJNA523634 | \
    efetch -format native -mode xml | \
    xtract -pattern EXPERIMENT_PACKAGE \
    -SRA RUN@accession \
    -NAME SAMPLE@alias \
    -group SAMPLE_ATTRIBUTES \
    -block SAMPLE_ATTRIBUTES \
    -if VALUE -contains aegypti \
    -element "&SRA","&NAME",VALUE > bennett.tsv

# Hery / Guadeloupe / French Guiana
esearch -db sra -query PRJNA600474 | \
    efetch -format native -mode xml | \
    xtract -pattern EXPERIMENT_PACKAGE \
    -SRA RUN@accession \
    -NAME SAMPLE@alias \
    -group SAMPLE_ATTRIBUTES \
    -block SAMPLE_ATTRIBUTES \
    -element "&SRA","&NAME",VALUE > hery.tsv

# Hernandez / Mexico
esearch -db sra -query PRJNA836318 | \
    efetch -format native -mode xml | \
    xtract -pattern EXPERIMENT_PACKAGE \
    -SRA RUN@accession \
    -NAME SAMPLE@alias \
    -group SAMPLE_ATTRIBUTES \
    -block SAMPLE_ATTRIBUTES \
    -element "&SRA","&NAME",VALUE > hernandez.tsv

# Rodpai / Thailand
esearch -db sra -query PRJNA919511 | \
    efetch -format native -mode xml | \
    xtract -pattern EXPERIMENT_PACKAGE \
    -SRA RUN@accession \
    -NAME SAMPLE@alias \
    -group SAMPLE_ATTRIBUTES \
    -block SAMPLE_ATTRIBUTES \
    -if "&NAME" -contains Aey \
    -element "&SRA","&NAME",VALUE > rodpai.tsv

# Prefetch and download SRA data
threads=5
for study in *.tsv; do
    #mkdir -p ${study%%.tsv}
    cd ${study%%.tsv}
    #prefetch --max-size u --option-file <(cut -f 1 ../$study)
    cut -f 1 ../$study | xargs -n 1 fasterq-dump -O fastq_files -e $threads
    pigz fastq_files/*.fastq
    cd ..
done