# Data Generation Script for VarFish Course

This repository contains scripts for generating Data for the VarFish course.

## Build `bamsurgeon` Singularity Image

First, create the Docker image.

```
# git clone https://github.com/adamewing/bamsurgeon.git
# cd bamsurgeon
# docker build -t bamsurgeon:local-build .
# docker save bamsurgeon:local-build -o bamsurgeon.tar
# singularity build bamsurgeon.sif docker-archive://bamsurgeon.tar
## now scp it next to this `README.md` file ;-)
```

## Create `conda` environment for Snakemake

```
# mamba create -y -n varfish-course-scripts samtools=1.9 htslib=1.9 wget bedtools samtools
```

## Used Data

- IGSR has 30x WGS data aligned to GRCh38 here:
    - https://www.internationalgenome.org/data-portal/data-collection/30x-grch38
- Raw Data:
    - https://ftp-trace.ncbi.nlm.nih.gov/1000genomes/ftp/1000G_2504_high_coverage/additional_698_related/

```
# wget \
    https://ftp-trace.ncbi.nlm.nih.gov/1000genomes/ftp/1000G_2504_high_coverage/additional_698_related/1000G_698_related_high_coverage.sequence.index \
    https://ftp-trace.ncbi.nlm.nih.gov/1000genomes/ftp/1000G_2504_high_coverage/additional_698_related/20130606_g1k_3202_samples_ped_population.txt \
    https://ftp-trace.ncbi.nlm.nih.gov/1000genomes/ftp/1000G_2504_high_coverage/1000G_2504_high_coverage.sequence.index
# cat >IBS029.ped <<"EOF"
IBS029 HG01679 0 0 2 IBS EUR
IBS029 HG01680 0 0 1 IBS EUR
IBS029 HG01681 HG01680 HG01679 2 IBS EUR
EOF
# mkdir -p IBS029 && cd IBS029
# wget $(for url in $(egrep -h 'HG01679|HG01680|HG01681' *.index | cut -f 1); do echo $url $url.crai; done)
```

```
# REFERENCE=/fast/projects/cubit/current/static_data/reference/GRCh38/hs38/hs38.fa
# samtools view -@ 4 -T /fast/projects/cubit/current/static_data/reference/GRCh38/hs38DH/hs38DH.fa -b -o HG01679.final.bam HG01679.final.cram
```
