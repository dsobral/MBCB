# MBCB - Genetic Alterations and Functional Impact

For this session, we will analyse some data from the 2022 mpox (formerly known as Monkeypox) outbreak in Portugal.

We will start with data from the [first cases](https://doi.org/10.1038/s41591-022-01907-y) in early May 2022.

Briefly, total DNA was extracted from the clinical sample, and shotgun metagenomics sequencing was performed (more details available in the paper). Reads were human-depleted using BMTagger and subsequently mapped to the reference genome [MPXV-UK_P2 MT903344.1](https://www.ncbi.nlm.nih.gov/nuccore/MT903344.1).

We are going to analyse data from the sample [Monkeypox/PT0005/2022](https://www.ncbi.nlm.nih.gov/nuccore/ON585037.1). Due to ethical reasons, only the reads mapping to the monkeypox virus reference are submitted to public repositories, and that is what we will use in our analysis.

We first need to map the raw sequencing files to the reference genome. As you saw in the previous session, for this, we need the fastq files of our sample, as well as the fasta file of the reference genome. Then, we will use bwa to align the reads and samtools to transform the SAM file that is the output of bwa into a sorted, indexed, BAM file that will be used for further analysis.

-  **TASK**: Download the fasta file of the reference from https://www.ncbi.nlm.nih.gov/nuccore/MT903344.1 
<details><summary>Click Here to see a hint</summary><p>  
  Click on "Send to" > File > Format (FASTA) > Create File 
</p></details>

Now, we need to create an index of the reference genome to be able to use bwa to align the reads against the reference. As in the last session, we will make use of [docker](https://www.docker.com/) images to facilitate reproducible installation in (almost) any environment. The following tasks assume the environment you are working in has already docker installed and ready to use.

If you have not done so before, pull a docker image for bwa eg.:
  > docker pull biocontainers/bwa:v0.7.17-3-deb_cv1

-  **TASK**: Generate a bwa index of the fasta file you just downloaded
<details><summary>Click Here to see a hint</summary><p>  
> docker run -v /yourfolder:/data biocontainers/bwa:v0.7.17-3-deb_cv1 bwa index MT903344.1.fasta
</p></details>

  
https://www.ncbi.nlm.nih.gov/nuccore/NC_063383


docker pull ncbi/sra-tools:3.0.1

ONT data of first case; [ERR9769166](ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR976/006/ERR9769166/ERR9769166.fastq.gz)
fasterq-dump ERR9769166

Illumina data for case 5 (shotgun metagenomics): ERR9769167 (ONT); ERR9769171 (Illumina) 
fasterq-dump ERR9769171
It has a deletion at NC_063383:11,326-12,238
fasterq-dump ERR9769167

breseq -r NC_063383.1.gb -j 4 -n ERR9769171 -o ERR9769171_breseq ERR9769171_1.fastq.gz ERR9769171_2.fastq.gz


Case 10 (shotgun metagenomics)
fasterq-dump ERR9769176

[A retrospective overview of the 2022 outbreak](https://doi.org/10.1038/s41591-023-02542-x):

Routine surveillance of monkeypox 

[amplicon sequencing]
PT400 (B.1): ERR10513212
fasterq-dump ERR10513212

PT428 (A2.3): ERR10513231
fasterq-dump ERR10513231



docker pull staphb/minimap2:2.26

docker pull staphb/samtools:1.19

docker pull biocontainers/freebayes:v1.2.0-2-deb_cv1

docker pull jysgro/breseq:ub2304_py3114_R422_br0381_bt245


