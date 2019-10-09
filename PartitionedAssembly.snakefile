import os
import tempfile
import subprocess
import os.path
# Snakemake and working directories
SD = os.path.dirname(workflow.snakefile)

# Config
configfile: "partitioned_assembly.json"

wd=config["working_dir"]
sample=config["sample"]
haps=["1", "2"]
ref=config["ref"]

fai = open(ref + ".fai")
chroms = [l.strip().split()[0] for l in fai ]

rule all:
    input:
 #       hapBams=expand(wd + "/{chrom}.{hap}.bam", chrom=chroms,hap=haps),
#        convert_part=expand(wd + "/{chrom}.{hap}.bam.fasta", chrom=chroms,hap=haps),
#        assemble_part=expand(wd + "/{chrom}.{hap}.assembly.fasta", chrom=chroms,hap=haps),
#        assemble_remap=expand(wd + "/{chrom}.{hap}.assembly.fasta.bam", chrom=chroms,hap=haps),
#        asmBam=expand(wd + "/{chrom}.{hap}.assembly.fasta.bam", chrom=chroms,hap=haps),
        cons=expand("{chrom}.{hap}.assembly.consensus.fasta", chrom=chroms,hap=haps)


rule SplitBam:
    input:
        bam=config["bam"],
        vcf=config["vcf"]
    output:
        hapBams=temp(expand(wd + "/{{chrom}}.{hap}.bam", hap=haps)),
    params:
        grid_opts=config["grid_medium"],
        sd=SD,
        sample=config["sample"],
        wd=wd,
        ref=config["ref"]
    shell:"""
samtools view -h {input.bam} {wildcards.chrom} | {params.sd}/pbgreedyphase/partitionByPhasedSNVs  --sample {params.sample} --vcf {input.vcf} --sam=/dev/stdin --h1={output.hapBams[0]} --h2={output.hapBams[1]} --phaseStats={params.wd}/{wildcards.chrom}.phase_stats --ref={params.ref}
"""

rule MakeFasta:
    input:
        hapBam=wd + "/{chrom}.{hap}.bam"
    output:
        hapFasta=temp(wd + "/{chrom}.{hap}.bam.fasta")
    params:
        grid_opts=config["grid_small"]
    shell:"""
samtools fasta {input.hapBam} -F 2304 > {output.hapFasta}
"""

rule AssembleFasta:
    input:
        hapFasta=wd + "/{chrom}.{hap}.bam.fasta"
    output:
        asm=temp(wd + "/{chrom}.{hap}.raw_assembly.fasta")
    params:
        grid_opts=config["grid_large"],
        ref=config["ref"]
    shell:"""
d=$(dirname {input.hapFasta})
mkdir -p $d
gs=`cat {params.ref}.fai | awk -vc={wildcards.chrom} '{{ if ($1 == c) print $2;}}'`
/home/cmb-16/mjc/shared/software_packages/wtdbg2/wtdbg2 -x sq -t 16 -g$gs  -i {input.hapFasta} -fo {output.asm} -L5000
/home/cmb-16/mjc/shared/software_packages/wtdbg2/wtpoa-cns -t 16 -i {output.asm}.ctg.lay.gz -o {output.asm}
samtools faidx {output.asm}
"""

rule RemapBam:
    input:
        asm=wd + "/{chrom}.{hap}.raw_assembly.fasta",
        hapBam=wd + "/{chrom}.{hap}.bam"
    output:
        asmBam=temp(wd + "/{chrom}.{hap}.assembly.fasta.bam"),
    params:
        grid_opts=config["grid_large"]
    shell:"""
. /home/cmb-16/mjc/mchaisso/projects/phasedsv_dev/phasedsv/dep/build/bin/activate pacbio
pbmm2 index {input.asm} {input.asm}.mmi
pbmm2 align {input.asm}.mmi {input.hapBam} -j 16 | samtools sort -T $TMPDIR/{wildcards.chrom}.{wildcards.hap} -m4G -@2 -o {output.asmBam}
pbindex {output.asmBam}
samtools index {output.asmBam}
"""

rule CallConsensus:
    input:
        asm=wd + "/{chrom}.{hap}.raw_assembly.fasta",
        asmBam=wd + "/{chrom}.{hap}.assembly.fasta.bam",
    output:
        cons="{chrom}.{hap}.assembly.consensus.fasta",
    params:
        grid_opts=config["grid_manycore"]
    shell:"""
. /home/cmb-16/mjc/mchaisso/projects/phasedsv_dev/phasedsv/dep/build/bin/activate pacbio
arrow -j 16 --referenceFilename {input.asm} --noEvidenceConsensusCall lowercasereference -o {output.cons} {input.asmBam}
rm -rf {input.asm}.*
"""
