import os
import tempfile

#
# Locate the tempdir for grid processing
#
configfile: "phase.json"

if "TMPDIR" in os.environ:
    TMPDIR = os.environ['TMPDIR']
elif "TMPDIR" in config:
    TMPDIR = config['TMPDIR']
else:
    TMPDIR = tempfile.gettempdir()

SD = os.path.dirname(workflow.snakefile)

faiFile = open(config["ref"] +".fai")
chroms = [l.split()[0] for l in faiFile ]

wd="/staging/mjc/mchaisso/phasing"

sample=config["sample"]

shell.prefix("set +eu ")

rule all:
    input:
        calledVCF=sample+".unphased.vcf",
        indexedVcf=sample+".wh.vcf.gz"


rule CallVariantsByChrom:
    input:
        srBam=config["srbam"]
    output:
        calledChrom=sample+".unphased.{chrom}.vcf"
    resources:
        mem_gb=48,
        threads=4
    params:
        grid_opts="sbatch -c 4 --mem=48G --time=24:00:00 -p cmb"
    shell:"""
which java

/home/cmb-16/mjc/shared/software_packages/gatk-4.0.7.0/gatk HaplotypeCaller --input {input.srBam} --output {output.calledChrom} -L {wildcards.chrom} --reference /home/cmb-16/mjc/shared/references/hg38_noalts/hg38.no_alts.fasta 
"""

rule CombineUnphasedByChrom:
    input:
        chroms=expand(sample+".unphased.{chrom}.vcf", chrom=chroms)
    output:
        phasedGenome=sample+".unphased.vcf"
    params:
        grid_opts="sbatch -c 1 --mem=2G --time=1:00:00 -p cmb"
    resources:
        mem_gb=4,
        threads=1
    shell:"""
grep "^#" {input.chroms[0]} | grep -v "CHROM" > {output.phasedGenome}
echo "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tunknown" >> {output.phasedGenome}
cat {input.chroms} | grep -v "^#" >> {output.phasedGenome}
"""
        
rule PhaseByChromosome:
    input:
        lrBam=config["bam"],
        srVCF=sample+".unphased.vcf"
    output:
        phasedChromosome=sample+".phased.{chrom}.vcf"
    resources:
        mem_gb=16,
        threads=1
    params:
        ref=config["ref"],
        grid_opts="sbatch -c 1 --mem=24G --time=12:00:00 -p cmb"
    shell:"""
which whatshap
whatshap phase  {input.srVCF} {input.lrBam} --ignore-read-groups   --chromosome {wildcards.chrom} --reference {params.ref} -o /dev/stdout | egrep -w "^#|^{wildcards.chrom}" > {output.phasedChromosome}
"""


rule CombinePhasedChromosomes:
    input:
        chroms=expand(sample+".phased.{chrom}.vcf", chrom=chroms)
    output:
        phasedGenome=sample+".wh.vcf"
    params:
        grid_opts="sbatch -c 1 --mem=2G --time=1:00:00 -p cmb"
    resources:
        mem_gb=4,
        threads=1
    shell:"""
grep "^#" {input.chroms[0]} > {output.phasedGenome}
echo "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tunknown" >> {output.phasedGenome}
cat {input.chroms} | grep -v "^#" >> {output.phasedGenome}
"""


rule IndexVCF:
    input:
        vcf="{input}.wh.vcf"
    output:
        gzvcf="{input}.wh.vcf.gz"
    params:
        grid_opts="sbatch -c 1 --mem=4G --time=1:00:00 -p cmb"
    resources:
        mem_gb=4,
        threads=1
    shell:"""
bgzip {input.vcf}
tabix {input.vcf}.gz
"""
