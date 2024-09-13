#!/bin/bash

max_jobs=50  # 设置最大并行任务数

cat sample.txt | while read id; do
    mkdir -p "$id"
    for vcf in ALL.chr*.vcf.gz; do
        while [ $(jobs -r | wc -l) -ge $max_jobs ]; do
            sleep 1  # 等待有空闲的任务槽
        done
        {
            bcftools view -s "$id" -g het  -v snps -o "$id/${id}_${vcf%.vcf.gz}_hetsnps.vcf" -O z "$vcf" --threads 30
        } &
    done
    wait  # 等待所有后台进程完成
done
