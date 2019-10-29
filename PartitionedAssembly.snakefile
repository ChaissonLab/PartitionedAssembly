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
allChroms = [l.strip().split()[0] for l in fai ]
chroms = []
allowedAssemblers = ["wtdbg2", "flye"]

if "assembler" not in config:
    config["assembler"] = "wtdbg2"
elif config["assembler"] not in allowedAssemblers:
    print("Configuration error. The assembler must be one of " + ",".join(allowedAssemblers))
    sys.exit(1)

for chrom in allChroms:
    if chrom != "chrM" and chrom != "chrY":
        chroms.append(chrom)


shell.prefix("set +eu ")

rule all:
    input:
 #       hapBams=expand(wd + "/{chrom}.{hap}.bam", chrom=chroms,hap=haps),
#        convert_part=expand(wd + "/{chrom}.{hap}.bam.fasta", chrom=chroms,hap=haps),
#        assemble_part=expand(wd + "/{chrom}.{hap}.assembly.fasta", chrom=chroms,hap=haps),
#        assemble_remap=expand(wd + "/{chrom}.{hap}.assembly.fasta.bam", chrom=chroms,hap=haps),
#        asmBam=expand(wd + "/{chrom}.{hap}.assembly.fasta.bam", chrom=chroms,hap=haps),
        cons=expand("{chrom}.{hap}.assembly.consensus.fasta", chrom=chroms,hap=haps),
        combined=expand("assembly.{hap}.consensus.fasta", hap=haps)


rule SplitBam:
    input:
        bam=config["bam"],
        vcf=config["vcf"]
    output:
        hapBams=temp(expand(wd + "/{{chrom}}.{hap}.bam", hap=haps)),
    resources:
        mem_gb=4,
        threads=14
    params:
        grid_opts=config["grid_medium"],
        sd=SD,
        sample=config["sample"],
        wd=wd,
        ref=config["ref"]
    shell:"""
samtools merge - {input.bam} -R {wildcards.chrom} | samtools view -h - | {params.sd}/pbgreedyphase/partitionByPhasedSNVs  --sample {params.sample} --vcf {input.vcf} --sam=/dev/stdin --h1={output.hapBams[0]} --h2={output.hapBams[1]} --phaseStats={params.wd}/{wildcards.chrom}.phase_stats --ref={params.ref}
"""

rule MakeFasta:
    input:
        hapBam=wd + "/{chrom}.{hap}.bam"
    output:
        hapFasta=temp(wd + "/{chrom}.{hap}.bam.fasta")
    params:
        grid_opts=config["grid_small"]
    resources:
        mem_gb=1,
        threads=1
    shell:"""
samtools fasta {input.hapBam} -F 2304 > {output.hapFasta}
"""

rule AssembleFasta:
    input:
        hapFasta=wd + "/{chrom}.{hap}.bam.fasta"
    output:
        asm=temp(wd + "/{chrom}.{hap}.raw_assembly.fasta")
    resources:
        mem_gb=lambda wildcards, attempt: (attempt+1)*16,
        threads=16
    params:
        grid_opts=config["grid_large"],
        ref=config["ref"],
        read_type=config["read-type"],
        assembler=config["assembler"],
        working_directory=wd
    shell:"""

d=$(dirname {input.hapFasta})
mkdir -p $d

if [ "{params.assembler}" == "wtdbg2" ]; then
  if [ "{params.read_type}" == "CCS" ]; then
     preset=ccs
  else
     preset=sq
  fi

  gs=`cat {params.ref}.fai | awk -vc={wildcards.chrom} '{{ if ($1 == c) print $2;}}'`
  /home/cmb-16/mjc/shared/software_packages/wtdbg2/wtdbg2 -x $preset -t 16 -g$gs  -i {input.hapFasta} -fo {output.asm} -L5000
  /home/cmb-16/mjc/shared/software_packages/wtdbg2/wtpoa-cns -t 16 -i {output.asm}.ctg.lay.gz -o {output.asm}

  rm -f {output.asm}.*.gz {output.asm}.*.frg.*
fi
if [ "{params.assembler}" == "flye" ]; then
  
  if [ "{params.read_type}" == "CCS" ]; then
     readType="--pacbio-corr"
  else
     readType="--pacbio-raw"
  fi

  gs=`cat {params.ref}.fai | awk -vc={wildcards.chrom} '{{ if ($1 == c) print $2;}}'`
  
. /home/cmb-16/mjc/mchaisso/projects/phasedsv_dev/phasedsv/dep/build/bin/activate python2
echo "read type " $readType
/home/cmb-16/mjc/mchaisso/software/Flye/bin/flye $readType {input.hapFasta} --genome-size $gs -o {params.working_directory}/{wildcards.chrom}.{wildcards.hap}  -t 16 -i 1
mv {params.working_directory}/{wildcards.chrom}.{wildcards.hap}/assembly.fasta {output.asm}
fi

samtools faidx {output.asm}
"""

rule RemapBam:
    input:
        asm=wd + "/{chrom}.{hap}.raw_assembly.fasta",
        hapBam=wd + "/{chrom}.{hap}.bam"
    output:
        asmBam=temp(wd + "/{chrom}.{hap}.assembly.fasta.bam"),
    resources:
        mem_gb=lambda wildcards, attempt: attempt*16,
        threads=16
    params:
        grid_opts=config["grid_large"]
    shell:"""
set +e
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
        fasta=wd + "/{chrom}.{hap}.bam.fasta",
    output:
        cons="{chrom}.{hap}.assembly.consensus.fasta",
    resources:
        mem_gb=lambda wildcards, attempt: attempt*16,
        threads=16
    params:
        grid_opts=config["grid_manycore"],
        read_type=config["read-type"],
        assembler=config["assembler"]
    shell:"""
if [ ! -e {input.asmBam}.bai ]; then
    samtools index {input.asmBam}
fi
if [ "{params.read_type}" == "CCS" ]; then
    samtools view -h {input.asmBam} | gzip -c > {input.asmBam}.sam.gz
    racon -t 16 {input.fasta} {input.asmBam}.sam.gz {input.asm} > {output.cons}
    true
else
. /home/cmb-16/mjc/mchaisso/projects/phasedsv_dev/phasedsv/dep/build/bin/activate pacbio

if [ ! -e {input.asmBam}.pbi ]; then 
    pbindex {input.asmBam}
fi

arrow -j 16 --referenceFilename {input.asm} --noEvidenceConsensusCall lowercasereference -o {output.cons} {input.asmBam}
fi
rm -rf {input.asm}.*

"""

rule CombineChromosomes:
    input:
        cons=expand("{chrom}.{{hap}}.assembly.consensus.fasta", chrom=chroms),
    output:
        combined="assembly.{hap}.consensus.fasta",
    resources:
        threads=1,
        mem_gb=4
    params:
        allChroms=chroms,
        grid_opts=config["grid_small"]
    shell:"""
for chrom in {params.allChroms}; do 
  cat $chrom.{wildcards.hap}.assembly.consensus.fasta | awk -vchrom=$chrom '{{ if (substr($1,0,1) == ">") {{ print ">"chrom"/"substr($1,2); }} else {{ print;}} }}'
done > {output.combined}
"""
