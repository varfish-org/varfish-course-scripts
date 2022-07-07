#!/usr/bin/bash

while read -r line; do
    IFS='|' read -r -a array <<< "$(echo $line | tr -d '[:space:]')"
    pedigree=$(echo ${array[1]} | sed -e 's|/.*||g')
    sample=${array[2]}
    snv1=${array[5]}
    snv2=${array[7]}
    sv=${array[9]}
    #echo $pedigree $sample $role $snv1 $snv2 $sv
    rm -f $pedigree/$sample.GRCh38.snv.txt $pedigree/$sample.GRCh38.sv.txt
    if [[ ! -z "$snv1" ]]; then
        echo "$snv1" \
        | tr ',' '\t' \
        >> $pedigree/$sample.GRCh38.snv.txt
        >&2 echo $pedigree/$sample.GRCh38.snv.txt
    fi
    if [[ ! -z "$snv2" ]]; then
        echo "$snv2" \
        | tr ',' '\t' \
        >> $pedigree/$sample.GRCh38.snv.txt
        >&2 echo $pedigree/$sample.GRCh38.snv.txt
    fi
    if [[ ! -z "$sv" ]]; then
        echo "$sv" \
        | tr ',' '\t' \
        >> $pedigree/$sample.GRCh38.sv.txt
        >&2 echo $pedigree/$sample.GRCh38.sv.txt
    fi
done < <(perl -p -e 's/\t/|/g' variants/variants.tsv | grep -v '^#')
