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
allowedAssemblers = ["wtdbg2", "flye", "falcon", "canu"]
allowedReadTypes= ["raw", "ccs", "ont"]


if "read-type" not in config:
    config["read-type"] = "raw"
elif config["read-type"] not in allowedReadTypes:
    print("Configuration error. read-type must be set in the json to one of " + ",".join(allowedReadTypes))
    sys.exit(1)
if "assembler" not in config:
    config["assembler"] = "wtdbg2"
elif config["assembler"] not in allowedAssemblers:
    print("Configuration error. The assembler must be one of " + ",".join(allowedAssemblers))
    sys.exit(1)

for chrom in allChroms:
    if chrom != "chrM" and chrom != "chrY":
        chroms.append(chrom)


#shell.prefix("set +eu ")

if config["assembler"] == "falcon":
    localrules: AssembleFasta


rule all:
    input:
        hapBams=expand(wd + "/{chrom}.{hap}.bam", chrom=chroms,hap=haps),
        convert_part=expand(wd + "/{chrom}.{hap}.bam.fasta", chrom=chroms,hap=haps),
        assemble_part=expand(wd + "/{chrom}.{hap}.raw_assembly.fasta", chrom=chroms,hap=haps),
        assemble_remap=expand(wd + "/{chrom}.{hap}.assembly.fasta.bam", chrom=chroms,hap=haps),
        asmBam=expand(wd + "/{chrom}.{hap}.assembly.fasta.bam", chrom=chroms,hap=haps),
        cons=expand("{chrom}.{hap}.assembly.consensus.fasta", chrom=chroms,hap=haps),
        combined=expand("assembly.{hap}.consensus.fasta", hap=haps)


rule SplitBam:
    input:
        bam=config["bam"],
        vcf=config["vcf"]
    output:
        hapBams=temp(expand(wd + "/{{chrom}}.{hap}.bam", hap=haps)),
#        hapBams=temp(expand(wd + "/{{chrom}}.{hap}.bam", hap=haps)),
    resources:
        mem_gb=4,
        threads=3
    params:
        grid_opts=config["grid_medium"],
        sd=SD,
        sample=config["sample"],
        wd=wd,
        ref=config["ref"],
        node_constraint=" "
    shell:"""
mkdir -p  {params.wd}
samtools merge - {input.bam} -R {wildcards.chrom} | samtools view -h - | {params.sd}/pbgreedyphase/partitionByPhasedSNVs  --sample {params.sample} --vcf {input.vcf} --sam=/dev/stdin --h1={output.hapBams[0]} --h2={output.hapBams[1]} --phaseStats={params.wd}/{wildcards.chrom}.phase_stats --ref={params.ref}
"""

rule MakeFasta:
    input:
        hapBam=wd + "/{chrom}.{hap}.bam"
    output:
#        hapFasta=temp(wd + "/{chrom}.{hap}.bam.fasta")
        hapFasta=wd + "/{chrom}.{hap}.bam.fasta"
    params:
        grid_opts=config["grid_medium"],
        node_constraint="",
        sd=SD,
    resources:
        mem_gb=lambda wildcards, attempt: 20,
        threads=1
    shell:"""

samtools fasta {input.hapBam} -F 2304 | {params.sd}/shuffleReads /dev/stdin {output.hapFasta} || true

"""

def GetSubmission(method):
    if method == "falcon":
        # falcon jobs are a pipeline that needs to run locally
        print("On Grid submit, using the shell " + method)
        gridSubmit=SD + "/GridWrapper.sh shell "
    else:
        print("On Grid submit, using the cluster " + method)
        # non-falcon run distributed
        gridSubmit=SD + "/GridWrapper.sh grid " + config["grid_large"]
    return gridSubmit


def GetThreads(assembler):
    if assembler == "falcon" or assembler == "canu":
        return 1
    else:
        return 16

def GetMemGB(assembler, attempt):
    if assembler=="falcon":
        return 32
    elif assembler == "flye":
        return min(96, 64+16*(attempt-1))
    elif assembler == "canu":
        return 8
    else:
        return min(48,16*(attempt+1))

rule AssembleFasta:
    input:
        hapFasta=wd + "/{chrom}.{hap}.bam.fasta"
    output:
        asm=wd + "/{chrom}.{hap}.raw_assembly.fasta"
#        asm=temp(wd + "/{chrom}.{hap}.raw_assembly.fasta")
    resources:
        mem_gb=lambda wildcards, attempt: GetMemGB(config["assembler"], attempt),
        threads=GetThreads(config["assembler"])
    params:
        grid_opts=config["grid_large"],
        ref=config["ref"],
        read_type=config["read-type"],
        assembler=config["assembler"],
        working_directory=wd,
        sd=SD,
        node_constraint="",
        partition=config["partition"]
    shell:"""
mkdir -p  {params.working_directory}
mkdir -p asm_{wildcards.chrom}_{wildcards.hap}
cwd=`pwd`
date > $cwd/asm_{wildcards.chrom}_{wildcards.hap}/timestamp.start
gs=`cat {params.ref}.fai | awk -vc={wildcards.chrom} '{{ if ($1 == c) print $2;}}'`
if [ "{params.assembler}" == "wtdbg2" ]; then
  if [ "{params.read_type}" == "ccs" ]; then
     preset=ccs
  else
     preset=sq
  fi

  cd asm_{wildcards.chrom}_{wildcards.hap} && \
  /home/cmb-16/mjc/shared/software_packages/wtdbg2/wtdbg2 -x $preset -t 16 -g$gs  -i {input.hapFasta} -fo {output.asm} -L5000 && \
  /home/cmb-16/mjc/shared/software_packages/wtdbg2/wtpoa-cns -t 16 -i {output.asm}.ctg.lay.gz -o {output.asm} && \
   rm {output.asm}.* && \
  cd .. && \
  rm -f {output.asm}.*.gz {output.asm}.*.frg.*
fi

if [ "{params.assembler}" == "flye" ]; then
  
  if [ "{params.read_type}" == "ccs" ]; then
     readType="--pacbio-corr"
  elif ["{params.read_type"} == "raw" ]; then
     readType="--pacbio-raw"
  else
     readType="--nano-raw"
  fi

  
. /home/cmb-16/mjc/mchaisso/projects/phasedsv_dev/phasedsv/dep/build/bin/activate python2
echo "read type " $readType
/home/cmb-16/mjc/mchaisso/software/Flye/bin/flye $readType {input.hapFasta} --genome-size $gs -o {params.working_directory}/{wildcards.chrom}.{wildcards.hap}  -t 16 -i 1
mv {params.working_directory}/{wildcards.chrom}.{wildcards.hap}/assembly.fasta {output.asm}
fi

if [ "{params.assembler}" == "falcon" ]; then

  if [ "{params.read_type}" == "ccs" ]; then
     cfgVer="ccs"
  else
     cfgVer="raw"
  fi
  fileName={input.hapFasta}
  if [ "{params.read_type}" == "ont" ]; then
    if [ ! -e {input.hapFasta}.ont.fasta ]; then
       cat '{input.hapFasta}' | awk 'BEGIN {{ i=1; }} {{ if (substr($1,0,1) != ">") {{ l=length($1); print ">ont_run/"i"/1_"l; print $1; i+=1; }} }}' > {input.hapFasta}.ont.fasta
    fi
    fileName={input.hapFasta}.ont.fasta
  fi
  wd=$PWD
  echo "Running falcon assembly $cfgVer"
  mkdir -p asm_{wildcards.chrom}_{wildcards.hap}
  cd asm_{wildcards.chrom}_{wildcards.hap} && \
  echo $fileName > input.fofn  && \
  cat {params.sd}/falcon_cfg/cfg.$cfgVer.1 | sed "s/PARTITION/{params.partition}/g" > falcon.cfg && \
  echo "genome_size = $gs" >> falcon.cfg && \
  cat {params.sd}/falcon_cfg/cfg.$cfgVer.2 >> falcon.cfg && \
  fc_run.py falcon.cfg && \
  cp 2-asm-falcon/p_ctg.fasta {output.asm} 
  if [ -e {output.asm} ]; then 
     rm -rf 0-rawreads 1-preads_ovl 2-asm-falcon
  fi
fi    

if [ "{params.assembler}" == "canu" ]; then
  
  rm -f canu_asm_{wildcards.chrom}_{wildcards.hap}/status.txt
  if [ ! -e canu_asm_{wildcards.chrom}_{wildcards.hap}/asm.contigs.fasta ]; then
  if [ "{params.read_type}" == "ccs" ]; then
     /home/cmb-16/mjc/shared/software_packages/canu/Linux-amd64/bin/canu -p asm -d canu_asm_{wildcards.chrom}_{wildcards.hap} genomeSize=$gs correctedErrorRate=0.015 OvlMerThreshold=200 batOptions="-eg 0.01 -eM 0.01 -dg 6 -db 6 -dr 1 -ca 50 -cp 5" gridOptions="-p {params.partition} --time=24:00:00" onSuccess={params.sd}/canu/OnSuccess.sh onFailure={params.sd}/canu/OnFailure.sh -pacbio-hifi {input.hapFasta} >& submit_canu_{wildcards.chrom}_{wildcards.hap}.txt 
  elif [ "{params.read_type}" == "raw" ]; then
    /home/cmb-16/mjc/shared/software_packages/canu/Linux-amd64/bin/canu -p asm -d canu_asm_{wildcards.chrom}_{wildcards.hap}  genomeSize=$gs gridOptions=" -p {params.partition} --time=24:00:00 " OvlMerThreshold=200 onSuccess={params.sd}/canu/OnSuccess.sh onFailure={params.sd}/canu/OnFailure.sh -pacbio-raw {input.hapFasta} >& submit_canu_{wildcards.chrom}_{wildcards.hap}.txt 
  else
     /home/cmb-16/mjc/shared/software_packages/canu/Linux-amd64/bin/canu -p asm -d canu_asm_{wildcards.chrom}_{wildcards.hap}  genomeSize=$gs gridOptions=" -p {params.partition} --time=24:00:00 " OvlMerThreshold=200 onSuccess={params.sd}/canu/OnSuccess.sh onFailure={params.sd}/canu/OnFailure.sh stopOnLowCoverage=7 -nanopore-raw {input.hapFasta} >& submit_canu_{wildcards.chrom}_{wildcards.hap}.txt 
  fi
  
  while [ ! -e canu_asm_{wildcards.chrom}_{wildcards.hap}/status.txt ]; do
    sleep 60
    jobid=`grep "Submitted batch job" submit_canu_{wildcards.chrom}_{wildcards.hap}.txt | awk '{{ print $4; }}'`
    sacct -j $jobid | grep -q "FAILED"
    status=$?
    if [ $status -eq 0 ]; then
       echo "failed" > canu_asm_{wildcards.chrom}_{wildcards.hap}/status.txt 
    fi
  done
  else
    echo "success" > canu_asm_{wildcards.chrom}_{wildcards.hap}/status.txt
  fi
  status=`cat canu_asm_{wildcards.chrom}_{wildcards.hap}/status.txt`;
  if [ "$status" != "success" ]; then
    echo "Status wasn't success "
    exit 1
  fi  
  echo "Made it to copy command."
  cp $cwd/canu_asm_{wildcards.chrom}_{wildcards.hap}/asm.contigs.fasta {output.asm}
  samtools faidx {output.asm}
  
fi

date > $cwd/asm_{wildcards.chrom}_{wildcards.hap}/timestamp.end
exit 0

"""


rule RemapBam:
    input:
        asm=wd + "/{chrom}.{hap}.raw_assembly.fasta",
        hapBam=wd + "/{chrom}.{hap}.bam"
    output:
        asmBam=temp(wd + "/{chrom}.{hap}.assembly.fasta.bam"),
#        asmBam=temp(wd + "/{chrom}.{hap}.assembly.fasta.bam"),
    resources:
        mem_gb=lambda wildcards, attempt: min(32,attempt*16 + 8),
        threads=16
    params:
        grid_opts=config["grid_large"],
        readtype=config["read-type"],
        node_constraint="--constraint=\"[E5-2640v3]\""
    shell:"""

if [ "{params.readtype}" = "ont" ]; then
  samtools fastq {input.hapBam} | minimap2 {input.asm} - -t 16 -a | samtools sort -T $TMPDIR/{wildcards.chrom}.{wildcards.hap} -m4G -@2 -o {output.asmBam}
else
   . /home/cmb-16/mjc/mchaisso/projects/phasedsv_dev/phasedsv/dep/build/bin/activate pacbio
   pbmm2 index {input.asm} {input.asm}.mmi
   pbmm2 align {input.asm}.mmi {input.hapBam} -j 16 | samtools sort -T $TMPDIR/{wildcards.chrom}.{wildcards.hap} -m4G -@2 -o {output.asmBam}
   pbindex {output.asmBam}
fi
fs=`stat --printf="%s" {output.asmBam}`
if [ $fs -lt 1000000 ]; then
   rm {output.asmBam}
   exit 1
fi

samtools index {output.asmBam}
"""

def GetConstraint(method):
    if method == "racon":
        return " --constraint=\"[E5-2640v3]\""
    else:
        return ""

rule CallConsensus:
    input:
        asm=wd + "/{chrom}.{hap}.raw_assembly.fasta",
        asmBam=wd + "/{chrom}.{hap}.assembly.fasta.bam",
        fasta=wd + "/{chrom}.{hap}.bam.fasta",
    output:
        cons="{chrom}.{hap}.assembly.consensus.fasta",
    resources:
        mem_gb=lambda wildcards, attempt: min(3,attempt)*16,
        threads=16
    params:
        grid_opts=config["grid_manycore"],
        consensus=config["consensus"],
        assembler=config["assembler"],
        node_constraint=GetConstraint(config["consensus"])
    shell:"""
set +e

if [ ! -e {input.asmBam}.bai ]; then
    samtools index {input.asmBam}
fi
if [ "{params.consensus}" == "racon" ]; then
    if [ ! -e {input.asmBam}.sam.gz ]; then
        samtools view -h {input.asmBam} | gzip -c > {input.asmBam}.sam.gz
    fi
    cat /proc/cpuinfo | grep sse4_1
    grep name /proc/cpuinfo
    /home/cmb-16/mjc/shared/software_packages/racon/build/bin/racon -t 16 {input.fasta} {input.asmBam}.sam.gz {input.asm} > {output.cons}
    true
else
  . /home/cmb-16/mjc/mchaisso/projects/phasedsv_dev/phasedsv/dep/build/bin/activate pacbio

  if [ ! -e {input.asmBam}.pbi ]; then 
      pbindex {input.asmBam}
  fi
  if [ ! -e {input.asm}.fai ]; then
      samtools faidx {input.asm}
  fi
  arrow -j 16 --referenceFilename {input.asm} --noEvidenceConsensusCall lowercasereference -o {output.cons} {input.asmBam}
fi

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
        grid_opts=config["grid_medium"],
        allChroms=chroms,
        node_constraint=""
    shell:"""
set +e
for chrom in {params.allChroms} ; do 
  cat $chrom.{wildcards.hap}.assembly.consensus.fasta | awk -vchrom=$chrom '{{ if (substr($1,0,1) == ">") {{ print ">"chrom"/"substr($1,2); }} else {{ print;}} }}'
done > {output.combined}
"""
