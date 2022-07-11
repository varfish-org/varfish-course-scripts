configfile: "config.yaml"

PATH_REF = config["path_ref"]

def input_default(_wildcards):
    result = []
    for pedigree in config["pedigrees"]:
        ped_name = pedigree["name"]
        for sample in pedigree["samples"]:
            sample_name = sample["name"]
            result += [
                f"output/{sample_name}/{sample_name}_1.fastq.gz",
                f"output/{sample_name}/{sample_name}_2.fastq.gz",
            ]
    return result


rule default:
    input: input_default


rule download:
    output:
        cram="work/download/{sample_name}/{sample_name}.cram",
        crai="work/download/{sample_name}/{sample_name}.crai",
    resources:
        mem="2G",
        time="02:00:00",
        partition="critical",
    run:
        download_url = None
        for pedigree in config["pedigrees"]:
            for sample in pedigree["samples"]:
                if sample["name"] == wildcards.sample_name:
                    download_url = sample["download_url"]
        shell(r"""
        set -x
        set -euo pipefail

        URL=%s
        wget -O {output.cram} $URL
        wget -O {output.crai} $URL.crai
        """ % download_url)


def get_pedigree(wildcards):
    for pedigree in config["pedigrees"]:
        for sample in pedigree["samples"]:
            if sample["name"] == wildcards.sample_name:
                return pedigree["name"]
    return None


rule cram_to_bam:
    input:
        cram="work/download/{sample_name}/{sample_name}.cram",
    output:
        bam="work/cram_to_bam/{sample_name}/{sample_name}.bam",
        bai="work/cram_to_bam/{sample_name}/{sample_name}.bam.bai",
    resources:
        mem="2G",
        time="08:00:00",
        partition="critical",
    threads: 4
    shell:
        r"""
        samtools view -@ 4 -T {PATH_REF} -O BAM -o {output.bam} {input.cram}
        samtools index {output.bam}
        """

rule bamsurgeon_snv:
    input:
        bam="work/cram_to_bam/{sample_name}/{sample_name}.bam",
        bai="work/cram_to_bam/{sample_name}/{sample_name}.bam.bai",
    params:
        pedigree=get_pedigree,
    output:
        bam="work/bamsurgeon-snv/{sample_name}/{sample_name}.bam",
        bai="work/bamsurgeon-snv/{sample_name}/{sample_name}.bam.bai",
    resources:
        mem="16G",
        time="24:00:00",
        partition="critical",
    shell:
        r"""
        export TMPDIR=$PWD/tmp
        mkdir -p $TMPDIR
        export TMPDIR=$(mktemp -d)
        trap "rm -rf $TMPDIR" ERR EXIT

        if [[ -e {params.pedigree}/{wildcards.sample_name}.GRCh38.snv.txt ]]; then
            singularity run --env PATH=$PATH:/bin/velvet-1.2.10 --bind /data,/fast \
                bamsurgeon.sif \
                /usr/bin/python3.8 \
                    /bamsurgeon/bin/addsnv.py \
                    --reference {PATH_REF} \
                    -f {input.bam} \
                    -o $TMPDIR/unsorted.bam \
                    --tmpdir $TMPDIR \
                    -v {params.pedigree}/{wildcards.sample_name}.GRCh38.snv.txt
            samtools sort -O BAM -o {output.bam} $TMPDIR/unsorted.bam
            samtools index {output.bam}
        else
            ln -sr {input.bam} {output.bam}
            ln -sr {input.bai} {output.bai}
        fi
        """

rule bamsurgeon_indel:
    input:
        bam="work/bamsurgeon-snv/{sample_name}/{sample_name}.bam",
        bai="work/bamsurgeon-snv/{sample_name}/{sample_name}.bam.bai",
    output:
        bam="work/bamsurgeon-indel/{sample_name}/{sample_name}.bam",
        bai="work/bamsurgeon-indel/{sample_name}/{sample_name}.bam.bai",
    resources:
        mem="16G",
        time="24:00:00",
        partition="critical",
    params:
        pedigree=get_pedigree,
    shell:
        r"""
        export TMPDIR=$PWD/tmp
        mkdir -p $TMPDIR
        export TMPDIR=$(mktemp -d)
        trap "rm -rf $TMPDIR" ERR EXIT

        if [[ -e {params.pedigree}/{wildcards.sample_name}.GRCh38.indel.txt ]]; then
            singularity run --env PATH=$PATH:/bin/velvet-1.2.10 --bind /data,/fast \
                bamsurgeon.sif \
                /usr/bin/python3.8 \
                    /bamsurgeon/bin/addindel.py \
                    --reference {PATH_REF} \
                    -f {input.bam} \
                    -o $TMPDIR/unsorted.bam \
                    --tmpdir $TMPDIR \
                    -v {params.pedigree}/{wildcards.sample_name}.GRCh38.indel.txt
            samtools sort -O BAM -o {output.bam} $TMPDIR/unsorted.bam
            samtools index {output.bam}
        else
            ln -sr {input.bam} {output.bam}
            ln -sr {input.bai} {output.bai}
        fi
        """

rule bamsurgeon_sv:
    input:
        bam="work/bamsurgeon-indel/{sample_name}/{sample_name}.bam",
        bai="work/bamsurgeon-indel/{sample_name}/{sample_name}.bam.bai",
    output:
        bam="work/bamsurgeon-sv/{sample_name}/{sample_name}.bam",
        bai="work/bamsurgeon-sv/{sample_name}/{sample_name}.bam.bai",
    resources:
        mem="16G",
        time="24:00:00",
        partition="critical",
    params:
        pedigree=get_pedigree,
    shell:
        r"""
        export TMPDIR=$PWD/tmp
        mkdir -p $TMPDIR
        export TMPDIR=$(mktemp -d)
        trap "rm -rf $TMPDIR" ERR EXIT

        if [[ -e {params.pedigree}/{wildcards.sample_name}.GRCh38.sv.txt ]]; then
            singularity run --env PATH=$PATH:/bin/velvet-1.2.10 --bind /data,/fast \
                bamsurgeon.sif \
                /usr/bin/python3.8 \
                    /bamsurgeon/bin/addsv.py \
                    --reference {PATH_REF} \
                    -f {input.bam} \
                    -o $TMPDIR/unsorted.bam \
                    --tmpdir $TMPDIR \
                    -v {params.pedigree}/{wildcards.sample_name}.GRCh38.sv.txt
            samtools sort -O BAM -o {output.bam} $TMPDIR/unsorted.bam
            samtools index {output.bam}
        else
            ln -sr {input.bam} {output.bam}
            ln -sr {input.bai} {output.bai}
        fi
        """

rule bam_to_fastq:
    input:
        bam="work/bamsurgeon-sv/{sample_name}/{sample_name}.bam",
    output:
        fastq_1="output/{sample_name}/{sample_name}_1.fastq.gz",
        fastq_2="output/{sample_name}/{sample_name}_2.fastq.gz",
        md5_1="output/{sample_name}/{sample_name}_1.fastq.gz.md5",
        md5_2="output/{sample_name}/{sample_name}_2.fastq.gz.md5",
    resources:
        mem="2G",
        time="24:00:00",
        partition="critical",
    shell:
        r"""
        export TMPDIR=$(mktemp -d)
        trap "rm -rf $TMPDIR" EXIT ERR

        samtools sort -n -o $TMPDIR/qname-sort.bam {input.bam}
        bedtools bamtofastq -i $TMPDIR/qname-sort.bam \
            -fq >(gzip -c >{output.fastq_1}) \
            -fq2 >(gzip -c >{output.fastq_2})

        pushd $(dirname {output.fastq_1})
        md5sum $(basename {output.fastq_1}) > $(basename {output.fastq_1}).md5
        md5sum $(basename {output.fastq_2}) > $(basename {output.fastq_2}).md5
        """
