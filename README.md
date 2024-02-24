# MBCB - Genetic Alterations and Functional Impact

For this session, we will analyse some data from the 2022 mpox (formerly known as Monkeypox) outbreak in Portugal.

We will start with data from the [first cases](https://doi.org/10.1038/s41591-022-01907-y) in early May 2022.

Briefly, total DNA was extracted from the clinical sample, and shotgun metagenomics sequencing was performed (more details available in the paper). Reads were human-depleted using BMTagger and subsequently mapped to the reference genome [MPXV-UK_P2 MT903344.1](https://www.ncbi.nlm.nih.gov/nuccore/MT903344.1).

We are going to analyse data from the sample [Monkeypox/PT0005/2022](https://www.ncbi.nlm.nih.gov/nuccore/ON585037.1). Due to ethical reasons, only the reads mapping to the monkeypox virus reference are submitted to public repositories, and that is what we will use in our analysis.

We first need to map the raw sequencing files to the reference genome. As you saw in the previous session, for this, we need the fastq files of our sample, as well as the fasta file of the reference genome. Then, we will use bwa to align the reads and samtools to transform the SAM file that is the output of bwa into a sorted, indexed, BAM file that will be used for further analysis.

**TASK**: Download the fasta file of the reference from https://www.ncbi.nlm.nih.gov/nuccore/MT903344.1 
<details><summary>Click Here to see a hint</summary><p>  
  Click on: Send to > File > Format (FASTA) > Create File 
</p></details>
<br/>

As in the last session, we will make use of [docker](https://www.docker.com/) images to facilitate reproducible installation of different bioinformatic tools in (almost) any environment. The following tasks assume the environment you are working in has already docker installed and ready to use.

 **TASK**: Download the paired fastq files using sra-tools fasterq-dump with the repository id ERR9769171  
<details><summary>Click Here to see a hint</summary><p>  

  If you have not done so before, pull a docker image for sra-tools eg.:
> docker pull ncbi/sra-tools:3.0.1

Next, run fasterq-dump using the sra-tools docker image:
> docker run --rm -v $PWD:/data ncbi/sra-tools:3.0.1 fasterq-dump --outdir /data ERR9769171

Since these reads are already the reads that map to the monkeypox genome, they have been already quality processed, and thus can be used directly as they are for subsequent analyses.

</p></details>
<br/>


**TASK**: Generate an indexed BAM file with the alignments of the monkeypox sample reads against the MT903344.1 reference
<details><summary>Click Here to see a hint</summary><p>

If you have not done so before, pull a docker image for bwa and samtools eg.:
> docker pull biocontainers/bwa:v0.7.17-3-deb_cv1

> docker pull biocontainers/samtools:v1.9-4-deb_cv1

We first need to create an index of the reference genome to be able to use bwa to align the reads against the reference. 
> docker run --rm -v $PWD:/data biocontainers/bwa:v0.7.17-3-deb_cv1 bwa index MT903344.1.fasta

Next we use bwa to generate alignments:
> docker run --rm -v $PWD:/data biocontainers/bwa:v0.7.17-3-deb_cv1 bwa mem MT903344.1.fasta ERR9769171_1.fastq ERR9769171_2.fastq > ERR9769171.sam 

Finally, we use samtools to convert the sam to bam, sort it by position, and index it.

> docker run --rm -v $PWD:/data biocontainers/samtools:v1.9-4-deb_cv1 samtools view -Sb ERR9769171.sam > ERR9769171.bam

> docker run --rm -v $PWD:/data biocontainers/samtools:v1.9-4-deb_cv1 samtools sort -o ERR9769171.sorted.bam ERR9769171.bam

> docker run --rm -v $PWD:/data biocontainers/samtools:v1.9-4-deb_cv1 samtools index ERR9769171.sorted.bam

Note: you should now also have the file ERR9769171.sorted.bam.bai

</p></details>
<br/>

Now that we obtained our alignments, let's try to obtain variants, ie., differences from the reference. For this we will use [freebayes](https://github.com/freebayes/freebayes), a tool that is relatively simple to use, and has been used in several different contexts. Note that there are several programs to perform variant calling, and with varying performance in different situations. 

Let's pull a docker image for freebayes (note that this is not the most recent version):
> docker pull biocontainers/freebayes:v1.2.0-2-deb_cv1

We will start by looking at the available options to run the software
> docker run --rm biocontainers/freebayes:v1.2.0-2-deb_cv1 freebayes -h

As you can see, freebayes has several options, although it can used even without explicitly providing any of them. One important important parameter is ploidy, as it determines how many possible haplotypes we should expect have at a given locus. 

**Question**: Knowing that this is a virus, do you think using the default ploidy of 2 is a good idea?
<details><summary>Click Here to see a hint</summary><p>

</p></details>
<br/>


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


