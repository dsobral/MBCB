# MBCB - A brief (bioinformatics) history of the monkeypox outbreak in Portugal

Reference: 
https://www.ncbi.nlm.nih.gov/nuccore/NC_063383


[First case](https://doi.org/10.1038/s41591-022-01907-y):

docker pull ncbi/sra-tools:3.0.1

ONT data of first case; [ERR9769166](ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR976/006/ERR9769166/ERR9769166.fastq.gz)
fasterq-dump ERR9769166

Illumina data for case 5 (shotgun metagenomics): ERR9769167 (ONT); ERR9769171 (Illumina) 
fasterq-dump ERR9769171

fasterq-dump ERR9769167

Case 10 (shotgun metagenomics)
fasterq-dump ERR9769176

[A retrospective overview of the 2022 outbreak](https://doi.org/10.1038/s41591-023-02542-x):

Routine surveillance of monkeypox 

[amplicon sequencing]
PT400 (B.1): ERR10513212
fasterq-dump ERR10513212

PT428 (A2.3): ERR10513231
fasterq-dump ERR10513231


docker pull biocontainers/bwa:v0.7.17-3-deb_cv1

docker pull biocontainers/samtools:v1.9-4-deb_cv1

docker pull biocontainers/freebayes:v1.2.0-2-deb_cv1




