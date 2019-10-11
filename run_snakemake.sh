#!/usr/bin/env bash
snakemake -p -s /home/cmb-16/mjc/mchaisso/projects/PartitionedAssembly/PartitionedAssembly.snakefile -j 18 --restart-times 3 --cluster {params.grid_opts} -c {resources.threads} --mem={resources.mem_gb}G -k --rerun-incomplete
