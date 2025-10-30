import os
import pandas as pd
from snakemake.utils import min_version

#using
# - snakemake v6.9.1
# - using conda v22.9.0

# #ex:
# snakemake \
# --printshellcmds \
# --cluster 'qsub -V -cwd -j y -S /bin/bash -pe smp {threads} -q all.q -o {log} -e {log}' \
# --jobs 200 \
# --latency-wait 200 \
# --keep-going \
# --rerun-incomplete \
# --snake code/filter_map_genotype.snake \
# --use-conda \
# --jobname snk.NGR.{name}.jid{jobid}


##### set minimum snakemake version #####
min_version("6.9.1")

#set main project dir and work from there
proj_dir    = "/master/nplatt/sch_hae_its-nigeria"
data_dir    = f"{proj_dir}/data"
results_dir = f"{proj_dir}/results"
envs_dir    = f"{proj_dir}/envs"
logs_dir    = f"{results_dir}/logs"
seq_dir     = f"{data_dir}/seq_data"
tmp_dir     = f"{results_dir}/scratch"
gatk_bin    = f"{proj_dir}/bin/gatk-4.2.0.0/gatk"

# ============================
# SET REFERENCE GENOMES
# ============================
genomes = {
    "SH_V3": f"{data_dir}/SH_V3.fa"
}

# ============================
# GET SAMPLES
# ============================
df = pd.read_csv("its-nigeria_samplesheet.csv", sep=",", header=0)
samples =list(df["wgs_id"])
# samples=["Sb_NG_au_1.1", "c_Sh_NG_bo_2_2"]

# ============================
# DEFINE INTERVALS FOR SCATTERING
# ============================
# Split each genome into intervals for parallel processing
num_intervals = 10  # Number of intervals to split the genome into
interval_ids = [f"{i:04d}" for i in range(num_intervals)]  # Generate interval IDs (e.g., '0001', '0002')

localrules: 
    all,

# Define Snakemake rules
rule all:
    input:
        expand("{dir}/filtered_reads/{id}_filtered_R{read}.fq.gz", dir=results_dir, id=samples, read=["1", "2", "X"]),
        expand("{dir}/mapped_reads/{genome}/{id}_flagstat.txt", dir=results_dir, genome=genomes.keys(), id=samples),
        expand("{dir}/mapped_reads/{genome}/{id}_processed.cram", dir=results_dir, genome=genomes.keys(), id=samples),
        expand("{dir}/mapped_reads/{genome}/{id}_processed.cram.crai", dir=results_dir, genome=genomes.keys(), id=samples),
        expand("{dir}/intervals/{genome}/{interval_id}-scattered.interval_list",   dir=results_dir, genome=genomes.keys(), interval_id=interval_ids),
        expand("{dir}/haplotype_caller/{genome}/{sample}_{interval_id}.vcf.gz",  dir=results_dir, genome=genomes.keys(), sample=samples, interval_id=interval_ids),
        expand("{dir}/haplotype_caller/{genome}/{sample}.merged.vcf.gz",           dir=results_dir, genome=genomes.keys(), sample=samples),
        expand("{dir}/haplotype_caller/{genome}/{sample}.merged.vcf.gz.tbi",       dir=results_dir, genome=genomes.keys(), sample=samples),
        expand("{dir}/genotype_intervals/{genome}/{interval_id}-scattered.interval_list", dir=results_dir, genome=genomes.keys(), interval_id=interval_ids),
        expand("{dir}/haplotype_caller/{genome}/{sample}.merged.vcf.gz", dir=results_dir, genome=genomes.keys(), sample=samples),
        expand("{dir}/genotype/{genome}/unfiltered.vcf.gz", dir=results_dir, genome=genomes.keys()),
        expand("{dir}/filter_genotypes/{genome}/coarse_filtered.vcf.gz", dir=results_dir, genome=genomes.keys()),
        expand("{dir}/filter_genotypes/{genome}/coarse_filtered.indv_freq.tbl", dir=results_dir, genome=genomes.keys()),
        expand("{dir}/filter_genotypes/{genome}/coarse_filtered.site_freq.tbl", dir=results_dir, genome=genomes.keys()),
        expand("{dir}/mosdepth/{genome}/{sample}.mosdepth.summary.txt", dir=results_dir, genome=genomes.keys(), sample=samples)


rule filter_reads:
    input:
        r1           = seq_dir + "/{id}_R1.fq.gz",
        r2           = seq_dir + "/{id}_R2.fq.gz",
        adapter_file = data_dir + "/adapters.fas"
    output:
        r1_pe = results_dir + "/filtered_reads/{id}_filtered_R1.fq.gz",
        r2_pe = results_dir + "/filtered_reads/{id}_filtered_R2.fq.gz",
        rx_se = results_dir + "/filtered_reads/{id}_filtered_RX.fq.gz",
        r1_se = results_dir + "/filtered_reads/{id}_filtered_SE_R1.fq.gz",
        r2_se = results_dir + "/filtered_reads/{id}_filtered_SE_R2.fq.gz"
    threads:
        12
    log:
        logs_dir + "/filter_reads.{id}.log"
    conda:
        envs_dir + "/trimmomatic.yml"
    shell:
        """
        trimmomatic \
            PE \
            -threads {threads} \
            -phred33 \
            {input.r1} \
            {input.r2} \
            {output.r1_pe} \
            {output.r1_se} \
            {output.r2_pe} \
            {output.r2_se} \
            LEADING:10 \
            TRAILING:10 \
            SLIDINGWINDOW:4:15 \
            MINLEN:36 \
            ILLUMINACLIP:{input.adapter_file}:2:30:10:1:true

        zcat {output.r1_se} {output.r2_se} | gzip >{output.rx_se} 
        """

rule bbmap:
    input:
        r1_read_fq = rules.filter_reads.output.r1_pe,
        r2_read_fq = rules.filter_reads.output.r2_pe,
        ref = genomes["SH_V3"],
    output:
        sam   = temp(f"{results_dir}/mapped_reads/{genome}/{id}.sam")
    threads:
        64
    log:
        logs_dir + "/pe_bbmap.{id}"
    conda:
        envs_dir + "/bbmap.yml"
    shell:
        """
        bbmap.sh \
            -ref={input.ref} \
            -nodisk \
            -in={input.r1_read_fq} \
            -in2={input.r2_read_fq} \
            -threads={threads} \
            ambig=toss \
            interleaved=false \
            -Xmx128g \
            -eoom \
            -out={output.sam} \
            minid=0.9
        """

rule sort_bam:
    input:
        sam = rules.bbmap.output.sam
    output:
        bam = f"{results_dir}/mapped_reads/{genome}/{id}_sorted.bam",
    threads:
        12
    log:
        logs_dir + "/sort_bam.{id}"
    conda:
        envs_dir + "/samtools.yml"
    shell:
        """
        samtools view -Sb {input.sam} | samtools sort -o {output.bam}
        """

rule add_readgroups:
    input:
        bam = rules.sort_bam.output.bam,
    output:
        bam = f"{results_dir}/mapped_reads/{genome}/{id}_rg.bam"
    threads:
        1
    log:
        logs_dir + "/add_readgroups.{id}"
    shell:
        """
        {gatk_bin} AddOrReplaceReadGroups \
            --INPUT {input.bam} \
            --OUTPUT {output.bam} \
            --RGPU unk \
            --RGLB library1 \
            --RGPL illumina \
            --RGSM {wildcards.id} \
            --RGID {wildcards.id}
        """

rule mark_duplicates_in_bam:
    input:
        bam = rules.add_readgroups.output.bam
    output:
        metrics = f"{results_dir}/mapped_reads/{genome}/{id}_dupmetrics.log",
        bam     = f"{results_dir}/mapped_reads/{genome}/{id}_processed.bam"
    threads:
        2
    log:
        logs_dir + "/mark_duplicates_in_bam.{id}"
    shell:
        """
        {gatk_bin} MarkDuplicates \
            --INPUT {input.bam} \
            --OUTPUT {output.bam} \
            --METRICS_FILE {output.metrics} \
            --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 900 \
            --ASSUME_SORT_ORDER coordinate
        """

rule index_bam:
    input:
        bam = rules.mark_duplicates_in_bam.output.bam
    output:
        bai = results_dir + "/mapped_reads/{genome}/{id}_processed.bam.bai"
    threads:
        1
    resources:
        cores=lambda wildcards, threads: threads
    log:
        logs_dir + "/index_bam.{genome}.{id}"
    conda:
        envs_dir + "/samtools.yml"
    shell:
        """
        samtools index {input.bam}
        """

rule flagstat:
    input:
        bam = rules.sort_bam.output.bam
    output:
        txt = results_dir + "/mapped_reads/{genome}/{id}_flagstat.txt"
    threads:
        1
    resources:
        cores=lambda wildcards, threads: threads
    log:
        logs_dir + "/flagstat.{genome}.{id}"
    conda:
        envs_dir + "/samtools.yml"
    shell:
        """
        samtools flagstat {input.bam} > {output.txt}
        """

rule mosdepth:
    input:
        bam = rules.mark_duplicates_in_bam.output.bam,
        bai = rules.index_bam.output.bai,
    output:
        results_dir + "/mosdepth/{genome}/{id}.mosdepth.global.dist.txt",
        results_dir + "/mosdepth/{genome}/{id}.mosdepth.summary.txt",
        results_dir + "/mosdepth/{genome}/{id}.per-base.bed.gz",
        results_dir + "/mosdepth/{genome}/{id}.per-base.bed.gz.csi"
    params:
        prefix = results_dir + "/mosdepth/{genome}/{id}"
    threads:
        4
    resources:
        cores=lambda wildcards, threads: threads
    log:
        logs_dir + "/mosdepth.{genome}.{id}"
    conda:
        envs_dir + "/mosdepth.yml"
    shell:
        """
        mosdepth --threads 4 {params.prefix} {input.bam}
        """

rule bam_to_cram:
    input:
        bam = rules.mark_duplicates_in_bam.output.bam,
    output:
        cram = results_dir + "/mapped_reads/{genome}/{id}_processed.cram"
    params:
        ref = lambda wildcards: genomes[wildcards.genome]  # Reference genome for CRAM compression
    threads:
        12
    resources:
        cores=lambda wildcards, threads: threads
    log:
        logs_dir + "/bam_to_cram.{genome}.{id}.log"
    conda:
        envs_dir + "/samtools.yml"
    shell:
        """
        samtools view -@ {threads} -C -T {params.ref} {input.bam} -o {output.cram}
        """

rule index_cram:
    input:
        cram = rules.bam_to_cram.output.cram
    output:
        crai = results_dir + "/mapped_reads/{genome}/{id}_processed.cram.crai"
    threads:
        1
    resources:
        cores=lambda wildcards, threads: threads
    log:
        logs_dir + "/index_cram.{genome}.{id}.log"
    conda:
        envs_dir + "/samtools.yml"
    shell:
        """
        samtools index {input.cram}
        """

# ============================
# RULE: SPLIT INTERVALS
# ============================
# This rule splits the genome into smaller intervals for parallel processing.
# By dividing the genome, the subsequent variant calling steps can be run faster.

rule split_intervals:
    input:
        ref = lambda wildcards: genomes[wildcards.genome]  # Reference genome
    output:
        interval_files = expand(
            results_dir + "/intervals/{genome}/{interval_id}-scattered.interval_list",
            genome="{genome}",
            interval_id=interval_ids
        )  # Explicitly list all interval files as outputs
    params:
        tmp_dir = tmp_dir,
        num_intervals = num_intervals,
        out_dir = results_dir + "/intervals/{genome}/"  # Directory for GATK output
    threads:
        12
    resources:
        cores=lambda wildcards, threads: threads
    log:
        logs_dir + "/split_intervals.{genome}.log"  # Log per genome
    shell:
        """
        mkdir -p {params.out_dir}
        {gatk_bin} SplitIntervals \
            --scatter-count {params.num_intervals} \
            -R {input.ref} \
            -O {params.out_dir}
        """

# ============================
# RULE: RUN HAPLOTYPECALLER
# ============================
# This rule runs GATK's HaplotypeCaller for each genome and interval, creating VCF files for each sample and interval.

rule haplotype_caller:
    input:
        crai = results_dir + "/mapped_reads/{genome}/{sample}_processed.cram.crai",
        cram = results_dir + "/mapped_reads/{genome}/{sample}_processed.cram",
        ref = lambda wildcards: genomes[wildcards.genome],
        interval = results_dir + "/intervals/{genome}/{interval_id}-scattered.interval_list"
    output:
        vcf = temp(results_dir + "/haplotype_caller/{genome}/{sample}_{interval_id}.vcf.gz")
    params:
        tmp_dir = tmp_dir
    threads:
        4
    resources:
        cores=lambda wildcards, threads: threads
    log:
        logs_dir + "/haplotype_caller.{genome}.{sample}.{interval_id}"
    shell:
        """

        echo "Running on compute node: $(hostname)"
        
        {gatk_bin} HaplotypeCaller --java-options "-Xmx8g" \
            -R {input.ref} \
            -I {input.cram} \
            -L {input.interval} \
            -O {output.vcf} \
            --native-pair-hmm-threads {threads} \
            --tmp-dir {params.tmp_dir} \
            -ERC GVCF \
            >{log}.stdout 2>{log}.stderr
        """

# ============================
# RULE: MERGE HAPLOTYPECALLER OUTPUTS
# ============================
# After running HaplotypeCaller for each interval, this rule merges the resulting VCFs into a single VCF for each sample.

rule merge_haplotype_caller:
    input:
        vcfs = expand(results_dir + "/haplotype_caller/{genome}/{sample}_{interval_id}.vcf.gz",  genome="{genome}", 
                                                                                                 sample="{sample}", 
                                                                                                 interval_id=interval_ids)
    output:
        merged_vcf = results_dir + "/haplotype_caller/{genome}/{sample}.merged.vcf.gz",  # Merged VCF output
        vcfs_list = temp(results_dir + "/haplotype_caller/{genome}/{sample}.hc_vcfs.list")
    params:
        tmp_dir = tmp_dir
    log:
        logs_dir + "/merge_haplotype_caller.{genome}.{sample}"  # Log file
    threads:
        2
    resources:
        cores=lambda wildcards, threads: threads
    shell:
        """
        ls {input.vcfs} >{output.vcfs_list}
        
        {gatk_bin} GatherVcfs \
            -O {output.merged_vcf} \
            -I {output.vcfs_list} \
            --TMP_DIR {params.tmp_dir} \
            >{log}.stdout 2>{log}.stderr
        """

# ============================
# RULE: INDEX MERGED HAPLOTYPECALLER VCFS
# ============================
# After merging HaplotypeCaller intervals from each individual, then index the VCF files for downstream processing

rule tabix_vcf:
    input:
        merged_vcf = results_dir + "/haplotype_caller/{genome}/{sample}.merged.vcf.gz",  # Merged VCF output
    output:
        tabix_index = results_dir + "/haplotype_caller/{genome}/{sample}.merged.vcf.gz.tbi"
    params:
        tmp_dir = tmp_dir
    log:
        logs_dir + "/tabix_vcf.{genome}.{sample}"  # Log file
    threads:
        1
    conda:
        envs_dir + "/vcftools.yml"
    resources:
        cores=lambda wildcards, threads: threads
    shell:
        """
        tabix -p vcf {input.merged_vcf}
        """

# ============================
# RULE: SPLIT INTERVALS
# ============================
# This rule splits a reference genome into multiple intervals
# to enable parallel processing in subsequent steps, speeding up
# variant calling and genotyping.
rule split_genotype_intervals:
    input:
        ref = lambda wildcards: genomes[wildcards.genome]  # Reference genome
    output:
        interval_files = expand(
            results_dir + "/genotype_intervals/{genome}/{interval_id}-scattered.interval_list",
            genome="{genome}",
            interval_id=interval_ids
        )
    params:
        tmp_dir = tmp_dir,
        num_intervals = num_intervals,
        out_dir = results_dir + "/genotype_intervals/{genome}/"  # Output directory for intervals
    threads:
        12
    resources:
        cores=lambda wildcards, threads: threads
    log:
        logs_dir + "/split_genotype_intervals.{genome}.log"
    shell:
        """
        mkdir -p {params.out_dir}
        {gatk_bin} SplitIntervals \
            --scatter-count {params.num_intervals} \
            -R {input.ref} \
            -O {params.out_dir}
        """

# ============================
# RULE: CREATE VCF LIST FOR GDBI
# ============================
# Generates a list of VCF files for input into GenomicsDBImport.
# The list is used to import multiple VCFs into GenomicsDB.
rule input_gdbi_vcf_list:
    input:
        vcfs = expand(results_dir + "/haplotype_caller/{genome}/{sample}.merged.vcf.gz", genome="{genome}", sample=samples)
    output:
        vcfs_list = results_dir + "/gdbimport/{genome}/vcfs.list"
    params:
        tmp_dir = tmp_dir
    log:
        logs_dir + "/input_gdbi_vcf_list.{genome}.log"
    threads:
        1
    resources:
        cores=lambda wildcards, threads: threads
    shell:
        """
        ls {input.vcfs} > {output.vcfs_list}
        """

# ============================
# RULE: GENOMICSDBIMPORT
# ============================
# Imports VCF files into a GenomicsDB for joint genotyping.
# Each genome's intervals are processed separately.
rule gdbimport:
    input:
        vcfs = expand(results_dir + "/haplotype_caller/{genome}/{sample}.merged.vcf.gz", genome="{genome}", sample=samples),
        tbis = expand(results_dir + "/haplotype_caller/{genome}/{sample}.merged.vcf.gz.tbi", genome="{genome}", sample=samples),
        interval_list = results_dir + "/genotype_intervals/{genome}/{interval_id}-scattered.interval_list",
        vcfs_list = results_dir + "/gdbimport/{genome}/vcfs.list"
    output:
        callset = results_dir + "/gdbimport/{genome}/{interval_id}/callset.json"
    params:
        tmp_dir = tmp_dir,
        interval_dir = results_dir + "/gdbimport/{genome}/{interval_id}"
    threads:
        12
    resources:
        cores=lambda wildcards, threads: threads
    log:
        logs_dir + "/gdbimport.{genome}.{interval_id}.log"
    shell:
        """
        [ -d {params.interval_dir} ] && rm -rf {params.interval_dir}
        {gatk_bin} GenomicsDBImport \
            -V {input.vcfs_list} \
            --genomicsdb-workspace-path {params.interval_dir} \
            -L {input.interval_list} \
            --reader-threads {threads} \
            --batch-size {threads} \
            --tmp-dir {params.tmp_dir} \
            --overwrite-existing-genomicsdb-workspace
        """

# ============================
# RULE: GENOTYPE GVCFS
# ============================
# Performs joint genotyping for each interval using GenotypeGVCFs.
# Outputs individual VCF files for each interval.
rule genotype:
    input:
        results_dir + "/gdbimport/{genome}/{interval_id}/callset.json",
        ref = lambda wildcards: genomes[wildcards.genome]
    output:
        vcf = temp(results_dir + "/genotype/{genome}/{interval_id}.vcf.gz"),
        tbi = temp(results_dir + "/genotype/{genome}/{interval_id}.vcf.gz.tbi")
    params:
        tmp_dir = tmp_dir,
        gdb_dir = results_dir + "/gdbimport/{genome}/{interval_id}"
    threads:
        8
    resources:
        cores=lambda wildcards, threads: threads
    log:
        logs_dir + "/genotype.{genome}.{interval_id}.log"
    shell:
        """
        {gatk_bin} GenotypeGVCFs \
            -R {input.ref} \
            -V gendb://{params.gdb_dir} \
            -O {output.vcf} \
            -new-qual \
            --tmp-dir {params.tmp_dir}
        """  

# ============================
# RULE: SORT GENOTYPED GVCFS
# ============================
# Sort each of the GVCFs before merging

rule sort_genotype_gvcfs:
    input:
        vcf = results_dir + "/genotype/{genome}/{interval_id}.vcf.gz",
        tbi = results_dir + "/genotype/{genome}/{interval_id}.vcf.gz.tbi"
    output:
        vcf = temp(results_dir + "/genotype/{genome}/{interval_id}.sorted.vcf.gz"),
        tbi = temp(results_dir + "/genotype/{genome}/{interval_id}.sorted.vcf.gz.tbi")
    params:
        tmp_dir = tmp_dir,
        picard_bin = "/master/nplatt/sch_hae_nigeria/bin/gatk-4.2.0.0/gatk-package-4.2.0.0-local.jar",
        memory=112
    threads:
        24
    resources:
        cores=lambda wildcards, threads: threads
    log:
        logs_dir + "/sort_genotype_gvcfs.{genome}.{interval_id}log"
    shell:
        """
        java -Xmx{params.memory}g -Djava.io.tmpdir={params.tmp_dir} \
            -jar {params.picard_bin} SortVcf \
            -I {input.vcf} \
            -O {output.vcf}
        """

# ============================
# RULE: MERGE INTERVAL VCFs
# ============================
# Merges all interval VCF files for a genome into a single VCF.
rule merge_gvcf:
    input:
        vcfs = expand(results_dir + "/genotype/{genome}/{interval_id}.sorted.vcf.gz", genome="{genome}", interval_id=interval_ids),
        tbis = expand(results_dir + "/genotype/{genome}/{interval_id}.sorted.vcf.gz.tbi", genome="{genome}", interval_id=interval_ids)
    output:
        vcf_list = temp(results_dir + "/genotype/{genome}/gdbi.list"),
        vcf = temp(results_dir + "/genotype/{genome}/raw.vcf.gz"),
        tbi = temp(results_dir + "/genotype/{genome}/raw.vcf.gz.tbi")
    params:
        tmp_dir = tmp_dir
    threads:
        12
    resources:
        cores=lambda wildcards, threads: threads
    log:
        logs_dir + "/merge_gvcf.{genome}.log"
    shell:
        """
        ls {input.vcfs} > {output.vcf_list}
        
        {gatk_bin} MergeVcfs \
            --MAX_RECORDS_IN_RAM 50000000 \
            -I {output.vcf_list} \
            -O {output.vcf} \
            --TMP_DIR {params.tmp_dir}
        """

# ============================
# RULE: SORT MERGED VCF
# ============================
# Sorts the merged VCF file to prepare for annotation or filtering.
rule sort_vcf:
    input:
        vcf = results_dir + "/genotype/{genome}/raw.vcf.gz",
        tbi = results_dir + "/genotype/{genome}/raw.vcf.gz.tbi"
    output:
        vcf = temp(results_dir + "/genotype/{genome}/sorted.vcf.gz"),
        tbi = temp(results_dir + "/genotype/{genome}/sorted.vcf.gz.tbi")
    params:
        tmp_dir = tmp_dir,
        picard_bin = "/master/nplatt/sch_hae_nigeria/bin/gatk-4.2.0.0/gatk-package-4.2.0.0-local.jar"
    threads:
        12
    resources:
        cores=lambda wildcards, threads: threads
    log:
        logs_dir + "/sort_vcf.{genome}.log"
    shell:
        """
        java -Xmx256g -Djava.io.tmpdir={params.tmp_dir} \
            -jar {params.picard_bin} SortVcf \
            -I {input.vcf} \
            -O {output.vcf}
        """

# ============================
# RULE: ANNOTATE VARIANTS
# ============================
# Annotates VCF with unique IDs for each variant based on "CHROM:POS".
# Compresses and indexes the output VCF.
rule annotate_vcf:
    input:
        vcf = results_dir + "/genotype/{genome}/sorted.vcf.gz",
        tbi = results_dir + "/genotype/{genome}/sorted.vcf.gz.tbi"
    output:
        vcf = results_dir + "/genotype/{genome}/unfiltered.vcf.gz",
        idx = results_dir + "/genotype/{genome}/unfiltered.vcf.gz.tbi"
    threads:
        12
    resources:
        cores=lambda wildcards, threads: threads
    conda:
        envs_dir + "/vcftools.yml"
    log:
        logs_dir + "/annotate_vcf.{genome}.log"
    shell:
        """
        bcftools annotate \
            --set-id +'%CHROM\:%POS' \
            {input.vcf} | \
        bgzip -c > {output.vcf}
        
        tabix -p vcf {output.vcf}
        """

# ============================
# RULE: SOFT FILTER VCF
# ============================
# This rule applies "soft" filters to the merged VCF.

rule soft_filter_vcf:
    input:
        vcf = results_dir + "/genotype/{genome}/unfiltered.vcf.gz",  # Input sorted VCF
        idx = results_dir + "/genotype/{genome}/unfiltered.vcf.gz.tbi",
        ref = lambda wildcards: genomes[wildcards.genome]  # Reference genome
    output:
        vcf = temp(results_dir + "/filter_genotypes/{genome}/soft_filtered.vcf.gz"),  # Output soft filtered VCF
        tbi = temp(results_dir + "/filter_genotypes/{genome}/soft_filtered.vcf.gz.tbi")  # Output soft filtered VCF
    threads:
        12
    resources:
        cores=lambda wildcards, threads: threads
    log:
        logs_dir + "/soft_filter_vcf.{genome}.log"  # Log file
    shell:
        """
        {gatk_bin} VariantFiltration \
                -R {input.ref} \
                -V {input.vcf} \
                -O {output.vcf} \
                --filter-name "QD_lt_2_or_missing"            \
                --filter-expression '!vc.hasAttribute("QD") || QD < 2.0' \
                --filter-name "MQ_lt_30_or_missing"           \
                --filter-expression '!vc.hasAttribute("MQ") || MQ < 30.0' \
                --filter-name "FS_gt_60_or_missing"           \
                --filter-expression '!vc.hasAttribute("FS") || FS > 60.0' \
                --filter-name "SOR_gt_3_or_missing"           \
                --filter-expression '!vc.hasAttribute("SOR") || SOR > 3.0' \
                --filter-name "MQRankSum_lt_-12.5_or_missing" \
                --filter-expression '!vc.hasAttribute("MQRankSum") || MQRankSum < -12.5' \
                --filter-name "ReadPosRankSum_lt_-8_or_missing" \
                --filter-expression '!vc.hasAttribute("ReadPosRankSum") || ReadPosRankSum < -8.0'
        """

# ============================
# RULE: SOFT FILTER VCF
# ============================
# This rule applies "soft" filters to the merged VCF.

rule coarse_gvcf_filter:
    input:
        vcf = results_dir + "/filter_genotypes/{genome}/soft_filtered.vcf.gz",
        tbi = results_dir + "/filter_genotypes/{genome}/soft_filtered.vcf.gz.tbi"  # Output soft filtered VCF
    output:
        vcf = temp(results_dir + "/filter_genotypes/{genome}/coarse_filtered.recode.vcf")
    params:
        min_gq      = 30,
        min_q       = 30,
        min_dp      = 8,
        min_alleles = 2,
        max_alleles = 2,
        out_prefix = results_dir + "/filter_genotypes/{genome}/coarse_filtered"
    threads:
        1
    conda:
        envs_dir + "/vcftools.yml"
    resources:
        cores=lambda wildcards, threads: threads
    log:
        logs_dir + "/initial_gvcf_filter.{genome}.log"  # Log per genome
    shell:
        """
        vcftools \
            --gzvcf {input.vcf} \
            --minGQ {params.min_gq} \
            --minDP {params.min_dp} \
            --minQ {params.min_q} \
            --min-alleles {params.min_alleles} \
            --max-alleles {params.max_alleles} \
            --remove-indels \
            --recode \
            --recode-INFO-all \
            --remove-filtered-all \
            --remove-filtered-geno-all \
            --out {params.out_prefix}
        """

# ============================
# RULE: COMPRESS AND INDEX GVCF
# ============================
# This rule applies "soft" filters to the merged VCF.

rule index_coarse_gvcf:
    input:
        vcf = results_dir + "/filter_genotypes/{genome}/coarse_filtered.recode.vcf"
    output:
        vcf = results_dir + "/filter_genotypes/{genome}/coarse_filtered.vcf.gz",
        tbi = results_dir + "/filter_genotypes/{genome}/coarse_filtered.vcf.gz.tbi"
    params:
    threads:
        1
    conda:
        envs_dir + "/vcftools.yml"
    resources:
        cores=lambda wildcards, threads: threads
    log:
        logs_dir + "/index_coarse_gvcf.{genome}.log"  # Log per genome
    shell:
        """
        bgzip -c {input.vcf} >{output.vcf}
        tabix -p vcf {output.vcf} 
        """


# ============================
# RULE: CALCULATE GENOTYPING RATE PER INDIVIDUAL
# ============================
# This rule applies calcualtes the genotyping rage per individual for downstream filtering

rule initial_genotyping_rate_per_individual:
    input:
        vcf = results_dir + "/filter_genotypes/{genome}/coarse_filtered.vcf.gz"
    output:
        tbl = results_dir + "/filter_genotypes/{genome}/coarse_filtered.indv_freq.tbl"
    threads:
        1
    conda:
        envs_dir + "/vcftools.yml"
    resources:
        cores=lambda wildcards, threads: threads
    log:
        logs_dir + "/initial_genotyping_rate_per_individual.{genome}.log"  # Log per genome
    shell:
        """
        vcftools \
            --gzvcf {input.vcf} \
            --missing-indv \
            --stdout \
            >{output.tbl}
        """

# ============================
# RULE: CALCULATE GENOTYPING RATE PER SITE
# ============================
# This rule applies calcualtes the genotyping rage per site for downstream filtering

rule initial_genotyping_rate_per_site:
    input:
        vcf = results_dir + "/filter_genotypes/{genome}/coarse_filtered.vcf.gz"
    output:
        tbl = results_dir + "/filter_genotypes/{genome}/coarse_filtered.site_freq.tbl"
    threads:
        1
    conda:
        envs_dir + "/vcftools.yml"
    resources:
        cores=lambda wildcards, threads: threads
    log:
        logs_dir + "/initial_genotyping_rate_per_site.{genome}.log"  # Log per genome
    shell:
        """
        vcftools \
            --gzvcf {input.vcf} \
            --missing-site \
            --stdout \
            >{output.tbl}
        """
