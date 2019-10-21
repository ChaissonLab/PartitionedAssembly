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
        indexedVcf=sample+".wh.vcf.gz"


rule PhaseByChromosome:
    input:
        lrBam=config["bam"],
        srVCF=config["srvcf"]
    output:
        phasedChromosome=sample+".{chrom}.vcf"
    params:
        ref=config["ref"]
    shell:"""
whatshap phase  {input.srVCF} {input.lrBam} --ignore-read-groups   --chromosome {wildcards.chrom} --reference {params.ref} -o /dev/stdout | egrep -w "^#|^{wildcards.chrom}" > {output.phasedChromosome}
"""


rule CombinePhasedChromosomes:
    input:
        chroms=expand(sample+".{chrom}.vcf", chrom=chroms)
    output:
        phasedGenome=sample+".wh.vcf"
    shell:"""
grep "^#" {input.chroms[0]} > {output.phasedGenome}
cat {input.chroms} | grep -v "^#" >> {output.phasedGenome}
"""


rule IndexVCF:
    input:
        vcf="{input}.wh.vcf"
    output:
        gzvcf="{input}.wh.vcf.gz"
    shell:"""
bgzip {input.vcf}
tabix {input.vcf}.gz
"""
