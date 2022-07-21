#!/usr/bin/bash

while read -r line; do
    IFS='|' read -r -a array <<< "$(echo $line | tr -d '[:space:]')"
    pedigree=$(echo ${array[1]} | sed -e 's|/.*||g')
    sample=${array[2]}
    snv1=${array[5]}
    snv2=${array[7]}
    sv=${array[9]}
    #echo $pedigree $sample $role $snv1 $snv2 $sv
    rm -f $pedigree/$sample.GRCh38.snv.txt $pedigree/$sample.GRCh38.sv.txt  # cruft
    rm -f $pedigree/$sample.GRCh38.yaml
    any=0
    if [[ ! -z "$snv1" ]]; then
        [[ "$any" -eq 0 ]] && echo "var_specs:" >$pedigree/$sample.GRCh38.yaml
        any=1
        cat >>$pedigree/$sample.GRCh38.yaml <<EOF
  - var_type: Snv
    chromosome: "$(echo $snv1 | cut -d , -f 1)"
    start: $(($(echo $snv1 | cut -d , -f 2) + 1))
    end: $(echo $snv1 | cut -d , -f 3)
    reference: "."
    alternative: "$(echo $snv1 | cut -d , -f 5)"
    aaf: $(echo $snv1 | cut -d , -f 4)
EOF
        >&2 echo $pedigree/$sample.GRCh38.yaml
    fi
    if [[ ! -z "$snv2" ]]; then
        [[ "$any" -eq 0 ]] && echo "var_specs:" >$pedigree/$sample.GRCh38.yaml
        any=1
        cat >>$pedigree/$sample.GRCh38.yaml <<EOF
  - var_type: Snv
    chromosome: "$(echo $snv2 | cut -d , -f 1)"
    start: $(($(echo $snv2 | cut -d , -f 2) + 1))
    end: $(echo $snv2 | cut -d , -f 3)
    reference: "."
    alternative: "$(echo $snv2 | cut -d , -f 5)"
    aaf: $(echo $snv2 | cut -d , -f 4)
EOF
        >&2 echo $pedigree/$sample.GRCh38.yaml
    fi
    if [[ ! -z "$sv" ]]; then
        [[ "$any" -eq 0 ]] && echo "var_specs:" >$pedigree/$sample.GRCh38.yaml
        any=1
        cat >>$pedigree/$sample.GRCh38.yaml <<EOF
  - var_type: "$(echo $sv | cut -d , -f 4)"
    chromosome: "$(echo $sv | cut -d , -f 1)"
    start: $(($(echo $sv | cut -d , -f 2) + 1))
    end: $(echo $sv | cut -d , -f 3)
    aaf: $(echo $sv | cut -d , -f 5)
EOF
        >&2 echo $pedigree/$sample.GRCh38.yaml
    fi
done < <(perl -p -e 's/\t/|/g' variants/variants.tsv | grep -v '^#')
