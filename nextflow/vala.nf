#!/usr/bin/env nextflow

expr = file(params.expr)
genes = file(params.genes)
targets = params.targets
if (targets != '' ) {
  targets = file(params.targets)
  targets = '-t ' + targets
}

if (params.clr) {

  process clr {
    errorStrategy params.error_strategy
    cpus params.clr_settings.cores
    memory params.clr_settings.memory
    publishDir params.out + '/networks/clr'
    queue params.slurm_partition
    clusterOptions '-n ' + params.clr_settings.tasks + '-A ' + params.slurm_account
    time params.clr_settings.ptime

    input:
    file expr
    file genes
    val targets

    output:
    file 'clr_network.tsv' into clr_network_raw
    
    script:
    if (params.executor == 'slurm')
    {
      """
      srun \
      mi -m CLR -o clr_network.tsv -B ${params.clr_settings.batchsize} \
         -b ${params.clr_settings.bins} -s ${params.clr_settings.spline} \
         -i ${expr} -g ${genes} -O ${params.clr_settings.cores} ${targets}
      """
    }
    else
    {
      """
      mpirun -np ${params.clr_settings.tasks} \
      mi -m CLR -o clr_network.tsv -B ${params.clr_settings.batchsize} \
         -b ${params.clr_settings.bins} -s ${params.clr_settings.spline} \
         -i ${expr} -g ${genes} -O ${params.clr_settings.cores} ${targets}
      """
    }
  }

  process clr_import {
    errorStrategy params.error_strategy
    validExitStatus 0,3
    cpus params.clr_settings.importcores
    memory params.clr_settings.importmem
    publishDir params.out + '/networks/clr'
    queue params.slurm_partition
    clusterOptions '-n 1 -A ' + params.slurm_account
    time params.clr_settings.itime

    input:
    file genes
    file clr_net from clr_network_raw
    output:
    file 'clr_network.sf' into clr_sf

    script:
    if (targets == '')
    {
      """
      seidr import -F lm -u -r -z -n ${params.clr_settings.importname} \
                   -i ${clr_net} -g ${genes} -o clr_network.sf \
                   -O ${params.clr_settings.importcores}
      """
    }
    else
    {
      """
      seidr import -F el -u -r -z -n ${params.clr_settings.importname} \
                   -i ${clr_net} -g ${genes} -o clr_network.sf \
                   -O ${params.clr_settings.importcores}
      """
    }
  }
}

if (params.aracne) {

  process aracne {
    errorStrategy params.error_strategy
    cpus params.aracne_settings.cores
    memory params.aracne_settings.memory
    publishDir params.out + '/networks/aracne'
    queue params.slurm_partition
    clusterOptions '-n ' + params.aracne_settings.tasks + '-A ' + params.slurm_account
    time params.aracne_settings.ptime

    input:
    file expr
    file genes
    val targets

    output:
    file 'aracne_network.tsv' into aracne_network_raw

    script:
    if (params.executor == 'slurm')
    {
      """
      srun \
      mi -m ARACNE -o aracne_network.tsv -B ${params.aracne_settings.batchsize} \
         -b ${params.aracne_settings.bins} -s ${params.aracne_settings.spline} \
         -i ${expr} -g ${genes} -O ${params.aracne_settings.cores} ${targets}
      """
    }
    else
    {
      """
      mpirun -n ${params.aracne_settings.tasks} \
      mi -m ARACNE -o aracne_network.tsv -B ${params.aracne_settings.batchsize} \
         -b ${params.aracne_settings.bins} -s ${params.aracne_settings.spline} \
         -i ${expr} -g ${genes} -O ${params.aracne_settings.cores} ${targets}
      """
   }
  }

  process aracne_import {
    errorStrategy 'finish'
    validExitStatus 0,3
    cpus params.aracne_settings.importcores
    memory params.aracne_settings.importmem
    publishDir params.out + '/networks/aracne'
    queue params.slurm_partition
    clusterOptions '-n 1 -A ' + params.slurm_account
    time params.aracne_settings.itime

    input:
    file genes
    file aracne_net from aracne_network_raw

    output:
    file 'aracne_network.sf' into aracne_sf

    script:
    if (targets == '')
    {
      """
      seidr import -F lm -u -r -z -n ${params.aracne_settings.importname} \
                   -i ${aracne_net} -g ${genes} -o aracne_network.sf \
                   -O ${params.aracne_settings.importcores}
      """
    }
    else
    {
      """
      seidr import -F el -u -r -z -n ${params.aracne_settings.importname} \
                   -i ${aracne_net} -g ${genes} -o aracne_network.sf \
                   -O ${params.aracne_settings.importcores}
      """
    }
  }
}

if (params.anova) {

  process anova {
    errorStrategy params.error_strategy
    cpus params.anova_settings.cores
    memory params.anova_settings.memory
    publishDir params.out + '/networks/anova'
    queue params.slurm_partition
    clusterOptions '-n ' + params.anova_settings.tasks + '-A ' + params.slurm_account
    time params.anova_settings.ptime

    input:
    file expr
    file genes
    file params.anova_settings.meta_file
    val targets

    output:
    file 'anova_network.tsv' into anova_network_raw

    """
    export OMP_NUM_THREADS=1
    anoverence -i ${expr} -g ${genes} -e ${params.anova_settings.meta_file} \
               -w ${params.anova_settings.weight} ${targets}
    """
  }

  process anova_import {
    errorStrategy params.error_strategy
    validExitStatus 0,3
    cpus params.anova_settings.importcores
    memory params.anova_settings.importmem
    publishDir params.out + '/networks/anova'
    queue params.slurm_partition
    clusterOptions '-n 1 -A ' + params.slurm_account
    time params.anova_settings.itime

    input:
    file genes
    file anova_net from anova_network_raw

    output:
    file 'anova_network.sf' into anova_sf

    script:
    if (targets == '')
    {
      """
      seidr import -F lm -u -r -z -n ${params.anova_settings.importname} \
                   -i ${aracne_net} -g ${genes} -o anova_network.sf \
                   -O ${params.anova_settings.importcores}
      """
    }
    else
    {
      """
      seidr import -F el -u -r -z -n ${params.anova_settings.importname} \
                   -i ${aracne_net} -g ${genes} -o anova_network.sf \
                   -O ${params.anova_settings.importcores}
      """
    }
  }
}

if (params.pearson) {

  process pearson {
    errorStrategy params.error_strategy
    cpus params.pearson_settings.cores
    memory params.pearson_settings.memory
    publishDir params.out + '/networks/pearson'
    queue params.slurm_partition
    clusterOptions '-n ' + params.pearson_settings.tasks + '-A ' + params.slurm_account
    time params.pearson_settings.ptime

    input:
    file expr
    file genes
    val targets

    output:
    file 'pearson_network.tsv' into pearson_network_raw

    """
    export OMP_NUM_THREADS=1
    correlation -m pearson -i ${expr} -g ${genes} -o pearson_network.tsv \
                ${params.pearson_settings.scale} \
                ${params.pearson_settings.absolute} ${targets}
    """
  }

  process pearson_import {
    errorStrategy params.error_strategy
    validExitStatus 0,3
    cpus params.pearson_settings.importcores
    memory params.pearson_settings.importmem
    publishDir params.out + '/networks/pearson'
    queue params.slurm_partition
    clusterOptions '-n 1 -A ' + params.slurm_account
    time params.pearson_settings.itime

    input:
    file genes
    file pearson_net from pearson_network_raw

    output:
    file 'pearson_network.sf' into pearson_sf

    script:
    if (targets == '')
    {
      """
      seidr import -F lm -A -u -r -z -n ${params.pearson_settings.importname} \
                   -i ${pearson_net} -g ${genes} -o pearson_network.sf \
                   -O ${params.pearson_settings.importcores}
      """
    }
    else
    {
      """
      seidr import -F el -A -u -r -z -n ${params.pearson_settings.importname} \
                   -i ${pearson_net} -g ${genes} -o pearson_network.sf \
                   -O ${params.pearson_settings.importcores}
      """
    }
  }
}

if (params.spearman) {

  process spearman {
    errorStrategy params.error_strategy
    cpus params.spearman_settings.cores
    memory params.spearman_settings.memory
    publishDir params.out + '/networks/spearman'
    queue params.slurm_partition
    clusterOptions '-n ' + params.spearman_settings.tasks + '-A ' + params.slurm_account
    time params.spearman_settings.ptime

    input:
    file expr
    file genes
    val targets

    output:
    file 'spearman_network.tsv' into spearman_network_raw

    """
    export OMP_NUM_THREADS=1
    correlation -m spearman -i ${expr} -g ${genes} -o spearman_network.tsv \
                ${params.spearman_settings.scale} \
                ${params.spearman_settings.absolute} ${targets}
    """
  }

  process spearman_import {
    errorStrategy params.error_strategy
    validExitStatus 0,3
    cpus params.spearman_settings.importcores
    memory params.spearman_settings.importmem
    publishDir params.out + '/networks/spearman'
    queue params.slurm_partition
    clusterOptions '-n 1 -A ' + params.slurm_account
    time params.spearman_settings.itime

    input:
    file genes
    file spearman_net from spearman_network_raw

    output:
    file 'spearman_network.sf' into spearman_sf

    script:
    if (targets == '')
    {
      """
      seidr import -F lm -A -u -r -z -n ${params.spearman_settings.importname} \
                   -i ${spearman_net} -g ${genes} -o spearman_network.sf \
                   -O ${params.spearman_settings.importcores}
      """
    }
    else
    {
      """
      seidr import -F el -A -u -r -z -n ${params.spearman_settings.importname} \
                   -i ${spearman_net} -g ${genes} -o spearman_network.sf \
                   -O ${params.spearman_settings.importcores}
      """
    }
  }
}

if (params.elnet) {

  process elnet {
    errorStrategy params.error_strategy
    cpus params.elnet_settings.cores
    memory params.elnet_settings.memory
    publishDir params.out + '/networks/elnet'
    queue params.slurm_partition
    clusterOptions '-n ' + params.elnet_settings.tasks + '-A ' + params.slurm_account
    time params.elnet_settings.ptime

    input:
    file expr
    file genes
    val targets

    output:
    file 'elnet_network.tsv' into elnet_network_raw


    script:
    if (params.executor == 'slurm')
    {
      """
      srun  \
      el-ensemble -o elnet_network.tsv -B ${params.elnet_settings.batchsize} \
         -l ${params.elnet_settings.min_lambda} -a ${params.elnet_settings.alpha} \
         -n ${params.elnet_settings.nlambda} \
         -X ${params.elnet_settings.max_experiment_size} \
         -x ${params.elnet_settings.min_experiment_size} \
         -P ${params.elnet_settings.max_predictor_size} \
         -p ${params.elnet_settings.min_predictor_size} \
         -e ${params.elnet_settings.ensemble} \
         -i ${expr} -g ${genes} \
         -O ${params.elnet_settings.cores} ${targets}
      """
    }
    else
    {
      """
      mpirun -np ${params.elnet_settings.tasks}  \
      el-ensemble -o elnet_network.tsv -B ${params.elnet_settings.batchsize} \
         -l ${params.elnet_settings.min_lambda} -a ${params.elnet_settings.alpha} \
         -n ${params.elnet_settings.nlambda} \
         -X ${params.elnet_settings.max_experiment_size} \
         -x ${params.elnet_settings.min_experiment_size} \
         -P ${params.elnet_settings.max_predictor_size} \
         -p ${params.elnet_settings.min_predictor_size} \
         -e ${params.elnet_settings.ensemble} \
         -i ${expr} -g ${genes} \
         -O ${params.elnet_settings.cores} ${targets}
      """
   }
  }

  process elnet_import {
    errorStrategy params.error_strategy
    validExitStatus 0,3
    cpus params.elnet_settings.importcores
    memory params.elnet_settings.importmem
    publishDir params.out + '/networks/elnet'
    queue params.slurm_partition
    clusterOptions '-n 1 -A ' + params.slurm_account
    time params.spearman_settings.itime

    input:
    file genes
    file elnet_net from elnet_network_raw
    output:
    file 'elnet_network.sf' into elnet_sf

    script:
    if (targets == '')
    {
      """
      seidr import -F m -r -z -n ${params.elnet_settings.importname} \
                   -i ${elnet_net} -g ${genes} -o elnet_network.sf \
                   -O ${params.elnet_settings.importcores}
      """
    }
    else
    {
      """
      seidr import -F el -r -z -n ${params.elnet_settings.importname} \
                   -i ${elnet_net} -g ${genes} -o elnet_network.sf \
                   -O ${params.elnet_settings.importcores}
      """
    }
  }
}

if (params.svm) {

  process svm {
    errorStrategy params.error_strategy
    cpus params.svm_settings.cores
    memory params.svm_settings.memory
    publishDir params.out + '/networks/svm'
    queue params.slurm_partition
    clusterOptions '-n ' + params.svm_settings.tasks + '-A ' + params.slurm_account
    time params.svm_settings.ptime

    input:
    file expr
    file genes
    val targets

    output:
    file 'svm_network.tsv' into svm_network_raw

    script:
    if (params.executor == 'slurm')
    {
      """
      srun \
      svm-ensemble -o svm_network.tsv -B ${params.svm_settings.batchsize} \
         -X ${params.svm_settings.max_experiment_size} \
         -x ${params.svm_settings.min_experiment_size} \
         -P ${params.svm_settings.max_predictor_size} \
         -p ${params.svm_settings.min_predictor_size} \
         -e ${params.svm_settings.ensemble} \
         -i ${expr} -g ${genes} \
         -O ${params.svm_settings.cores} ${targets}
      """
    }
    else
    {
      """
      mpirun -np ${params.svm_settings.tasks} \
      svm-ensemble -o svm_network.tsv -B ${params.svm_settings.batchsize} \
         -X ${params.svm_settings.max_experiment_size} \
         -x ${params.svm_settings.min_experiment_size} \
         -P ${params.svm_settings.max_predictor_size} \
         -p ${params.svm_settings.min_predictor_size} \
         -e ${params.svm_settings.ensemble} \
         -i ${expr} -g ${genes} \
         -O ${params.svm_settings.cores} ${targets}
      """
    }
  }

  process svm_import {
    errorStrategy params.error_strategy
    validExitStatus 0,3
    cpus params.svm_settings.importcores
    memory params.svm_settings.importmem
    publishDir params.out + '/networks/svm'
    queue params.slurm_partition
    clusterOptions '-n 1 -A ' + params.slurm_account
    time params.spearman_settings.itime

    input:
    file genes
    file svm_net from svm_network_raw

    output:
    file 'svm_network.sf' into svm_sf

    script:
    if (targets == '')
    {
      """
      seidr import -F m -r -z -n ${params.svm_settings.importname} \
                   -i ${svm_net} -g ${genes} -o svm_network.sf \
                   -O ${params.svm_settings.importcores}
      """
    }
    else
    {
      """
      seidr import -F el -r -z -n ${params.svm_settings.importname} \
                   -i ${svm_net} -g ${genes} -o svm_network.sf \
                   -O ${params.svm_settings.importcores}
      """
    }
  }
}

if (params.llr) {

  process llr {
    errorStrategy params.error_strategy
    cpus params.llr_settings.cores
    memory params.llr_settings.memory
    publishDir params.out + '/networks/llr'
    queue params.slurm_partition
    clusterOptions '-n ' + params.llr_settings.tasks + '-A ' + params.slurm_account
    time params.llr_settings.ptime

    input:
    file expr
    file genes
    val targets

    output:
    file 'llr_network.tsv' into llr_network_raw

    script:
    if (params.executor == 'slurm')
    {
      """
      srun \
      llr-ensemble -o llr_network.tsv -B ${params.llr_settings.batchsize} \
         -X ${params.llr_settings.max_experiment_size} \
         -x ${params.llr_settings.min_experiment_size} \
         -P ${params.llr_settings.max_predictor_size} \
         -p ${params.llr_settings.min_predictor_size} \
         -e ${params.llr_settings.ensemble} \
         -i ${expr} -g ${genes} \
         -O ${params.llr_settings.cores} ${targets}
      """
    }
    else
    {
      """
      mpirun -np ${params.llr_settings.tasks} \
      llr-ensemble -o llr_network.tsv -B ${params.llr_settings.batchsize} \
         -X ${params.llr_settings.max_experiment_size} \
         -x ${params.llr_settings.min_experiment_size} \
         -P ${params.llr_settings.max_predictor_size} \
         -p ${params.llr_settings.min_predictor_size} \
         -e ${params.llr_settings.ensemble} \
         -i ${expr} -g ${genes} \
         -O ${params.llr_settings.cores} ${targets}
      """
    }
  }

  process llr_import {
    errorStrategy params.error_strategy
    validExitStatus 0,3
    cpus params.llr_settings.importcores
    memory params.llr_settings.importmem
    publishDir params.out + '/networks/llr'
    queue params.slurm_partition
    clusterOptions '-n 1 -A ' + params.slurm_account
    time params.spearman_settings.itime

    input:
    file genes
    file llr_net from llr_network_raw

    output:
    file 'llr_network.sf' into llr_sf

    script:
    if (targets == '')
    {
      """
      seidr import -F m -r -z -n ${params.llr_settings.importname} \
                   -i ${llr_net} -g ${genes} -o llr_network.sf \
                   -O ${params.llr_settings.importcores}
      """
    }
    else
    {
      """
      seidr import -F el -r -z -n ${params.llr_settings.importname} \
                   -i ${llr_net} -g ${genes} -o llr_network.sf \
                   -O ${params.llr_settings.importcores}
      """
    }
  }
}

if (params.pcor) {

  process pcor {
    errorStrategy params.error_strategy
    cpus params.pcor_settings.cores
    memory params.pcor_settings.memory
    publishDir params.out + '/networks/pcor'
    queue params.slurm_partition
    clusterOptions '-n ' + params.pcor_settings.tasks + '-A ' + params.slurm_account
    time params.pcor_settings.ptime

    input:
    file expr
    file genes
    val targets

    output:
    file 'pcor_network.tsv' into pcor_network_raw

    """
    pcor -i ${expr} -g ${genes} -o pcor_network.tsv \
         ${params.pcor_settings.absolute} ${targets}
    """
  }

  process pcor_import {
    errorStrategy params.error_strategy
    validExitStatus 0,3
    cpus params.pcor_settings.importcores
    memory params.pcor_settings.importmem
    publishDir params.out + '/networks/pcor'
    queue params.slurm_partition
    clusterOptions '-n 1 -A ' + params.slurm_account
    time params.spearman_settings.itime

    input:
    file genes
    file pcor_net from pcor_network_raw

    output:
    file 'pcor_network.sf' into pcor_sf

    script:
    if (targets == '')
    {
      """
      seidr import -F lm -A -u -r -z -n ${params.pcor_settings.importname} \
                   -i ${pcor_net} -g ${genes} -o pcor_network.sf \
                   -O ${params.pcor_settings.importcores}
      """
    }
    else
    {
      """
      seidr import -F el -A -u -r -z -n ${params.pcor_settings.importname} \
                   -i ${pcor_net} -g ${genes} -o pcor_network.sf \
                   -O ${params.pcor_settings.importcores}
      """
    }
  }
}

if (params.narromi) {

  process narromi {
    errorStrategy params.error_strategy
    cpus params.narromi_settings.cores
    memory params.narromi_settings.memory
    publishDir params.out + '/networks/narromi'
    queue params.slurm_partition
    clusterOptions '-n ' + params.narromi_settings.tasks + '-A ' + params.slurm_account
    time params.narromi_settings.ptime

    input:
    file expr
    file genes
    val targets

    output:
    file 'narromi_network.tsv' into narromi_network_raw

    script:
    if (params.executor == 'slurm')
    {
      """
      srun \
      narromi -o narromi_network.tsv -B ${params.narromi_settings.batchsize} \
         -a ${params.narromi_settings.alpha} \
         -m ${params.narromi_settings.method} \
         -i ${expr} -g ${genes} \
         -O ${params.narromi_settings.tasks} ${targets}
      """
    }
    else
    {
      """
      mpirun -np ${params.narromi_settings.tasks} \
      narromi -o narromi_network.tsv -B ${params.narromi_settings.batchsize} \
         -a ${params.narromi_settings.alpha} \
         -m ${params.narromi_settings.method} \
         -i ${expr} -g ${genes} \
         -O ${params.narromi_settings.tasks} ${targets}
      """
    }
  }

  process narromi_import {
    errorStrategy params.error_strategy
    validExitStatus 0,3
    cpus params.narromi_settings.importcores
    memory params.narromi_settings.importmem
    publishDir params.out + '/networks/narromi'
    queue params.slurm_partition
    clusterOptions '-n 1 -A ' + params.slurm_account
    time params.spearman_settings.itime

    input:
    file genes
    file narromi_net from narromi_network_raw

    output:
    file 'narromi_network.sf' into narromi_sf

    script:
    if (targets == '')
    {
      """
      seidr import -F m -r -z -n ${params.narromi_settings.importname} \
                   -i ${narromi_net} -g ${genes} -o narromi_network.sf \
                   -O ${params.narromi_settings.importcores}
      """
    }
    else
    {
      """
      seidr import -F el -r -z -n ${params.narromi_settings.importname} \
                   -i ${narromi_net} -g ${genes} -o narromi_network.sf \
                   -O ${params.narromi_settings.importcores}
      """
    }
  }
}

if (params.tigress) {

  process tigress {
    errorStrategy params.error_strategy
    cpus params.tigress_settings.cores
    memory params.tigress_settings.memory
    publishDir params.out + '/networks/tigress'
    queue params.slurm_partition
    clusterOptions '-n ' + params.tigress_settings.tasks + '-A ' + params.slurm_account
    time params.tigress_settings.ptime

    input:
    file expr
    file genes
    val targets

    output:
    file 'tigress_network.tsv' into tigress_network_raw

    script:
    if (params.executor == 'slurm')
    {
      """
      srun \
      tigress -o tigress_network.tsv -B ${params.tigress_settings.batchsize} \
         -n ${params.tigress_settings.nlambda} \
         -l ${params.tigress_settings.min_lambda} \
         ${params.tigress_settings.scale} \
         -i ${expr} -g ${genes} \
         -O ${params.tigress_settings.cores} ${targets}
      """
    }
    else
    {
      """
      mpirun -np ${params.tigress_settings.tasks} \
      tigress -o tigress_network.tsv -B ${params.tigress_settings.batchsize} \
         -n ${params.tigress_settings.nlambda} \
         -l ${params.tigress_settings.min_lambda} \
         ${params.tigress_settings.scale} \
         -i ${expr} -g ${genes} \
         -O ${params.tigress_settings.cores} ${targets}
      """ 
    }
  }

  process tigress_import {
    errorStrategy params.error_strategy
    validExitStatus 0,3
    cpus params.tigress_settings.importcores
    memory params.tigress_settings.importmem
    publishDir params.out + '/networks/tigress'
    queue params.slurm_partition
    clusterOptions '-n 1 -A ' + params.slurm_account
    time params.spearman_settings.itime

    input:
    file genes
    file tigress_net from tigress_network_raw

    output:
    file 'tigress_network.sf' into tigress_sf

    script:
    if (targets == '')
    {
      """
      seidr import -F m -r -z -n ${params.tigress_settings.importname} \
                   -i ${tigress_net} -g ${genes} -o tigress_network.sf \
                   -O ${params.tigress_settings.importcores}
      """
    }
    else
    {
      """
      seidr import -F el -r -z -n ${params.tigress_settings.importname} \
                   -i ${tigress_net} -g ${genes} -o tigress_network.sf \
                   -O ${params.tigress_settings.importcores}
      """
    }
  }
}

if (params.genie3) {

  process genie3 {
    errorStrategy params.error_strategy
    cpus params.genie3_settings.cores
    memory params.genie3_settings.memory
    publishDir params.out + '/networks/genie3'
    queue params.slurm_partition
    clusterOptions '-n ' + params.genie3_settings.tasks + '-A ' + params.slurm_account
    time params.genie3_settings.ptime

    input:
    file expr
    file genes
    val targets

    output:
    file 'genie3_network.tsv' into genie3_network_raw

    script:
    if (params.executor == 'slurm')
    {
      """
      srun \
      genie3 -o genie3_network.tsv -B ${params.genie3_settings.batchsize} \
         -p ${params.genie3_settings.min_prop} \
         -a ${params.genie3_settings.alpha} \
         -N ${params.genie3_settings.min_node_size} \
         -m ${params.genie3_settings.mtry} \
         -n ${params.genie3_settings.ntree} \
         ${params.genie3_settings.scale} \
         -i ${expr} -g ${genes} \
         -O ${params.genie3_settings.cores} ${targets}
      """
    }
    else
    {
      """
      mpirun -np ${params.genie3_settings.tasks} \
      genie3 -o genie3_network.tsv -B ${params.genie3_settings.batchsize} \
         -p ${params.genie3_settings.min_prop} \
         -a ${params.genie3_settings.alpha} \
         -N ${params.genie3_settings.min_node_size} \
         -m ${params.genie3_settings.mtry} \
         -n ${params.genie3_settings.ntree} \
         ${params.genie3_settings.scale} \
         -i ${expr} -g ${genes} \
         -O ${params.genie3_settings.cores} ${targets}
      """
    }
  }

  process genie3_import {
    errorStrategy params.error_strategy
    validExitStatus 0,3
    cpus params.genie3_settings.importcores
    memory params.genie3_settings.importmem
    publishDir params.out + '/networks/genie3'
    queue params.slurm_partition
    clusterOptions '-n 1 -A ' + params.slurm_account
    time params.spearman_settings.itime

    input:
    file genes
    file genie3_net from genie3_network_raw

    output:
    file 'genie3_network.sf' into genie3_sf

    script:
    if (targets == '')
    {
      """
      seidr import -F m -r -z -n ${params.genie3_settings.importname} \
                   -i ${genie3_net} -g ${genes} -o genie3_network.sf \
                   -O ${params.genie3_settings.importcores}
      """
    }
    else
    {
      """
      seidr import -F el -r -z -n ${params.genie3_settings.importname} \
                   -i ${genie3_net} -g ${genes} -o genie3_network.sf \
                   -O ${params.genie3_settings.importcores}
      """
    }
  }
}

if (params.plsnet) {

  process plsnet {
    errorStrategy params.error_strategy
    cpus params.plsnet_settings.cores
    memory params.plsnet_settings.memory
    publishDir params.out + '/networks/plsnet'
    queue params.slurm_partition
    clusterOptions '-n ' + params.plsnet_settings.tasks + '-A ' + params.slurm_account
    time params.plsnet_settings.ptime

    input:
    file expr
    file genes
    val targets

    output:
    file 'plsnet_network.tsv' into plsnet_network_raw

    script:
    if (params.executor == 'slurm')
    {
      """
      srun \
      plsnet -o plsnet_network.tsv -B ${params.plsnet_settings.batchsize} \
         -c ${params.plsnet_settings.components} \
         -p ${params.plsnet_settings.predictors} \
         -e ${params.plsnet_settings.ensemble} \
         ${params.plsnet_settings.scale} \
         -i ${expr} -g ${genes} \
         -O ${params.plsnet_settings.cores} ${targets}
      """
    }
    else
    {
      """
      mpirun -np ${params.plsnet_settings.tasks} \
      plsnet -o plsnet_network.tsv -B ${params.plsnet_settings.batchsize} \
         -c ${params.plsnet_settings.components} \
         -p ${params.plsnet_settings.predictors} \
         -e ${params.plsnet_settings.ensemble} \
         ${params.plsnet_settings.scale} \
         -i ${expr} -g ${genes} \
         -O ${params.plsnet_settings.cores} ${targets}
      """
    }
  }

  process plsnet_import {
    errorStrategy params.error_strategy
    validExitStatus 0,3
    cpus params.plsnet_settings.importcores
    memory params.plsnet_settings.importmem
    publishDir params.out + '/networks/plsnet'
    queue params.slurm_partition
    clusterOptions '-n 1 -A ' + params.slurm_account
    time params.spearman_settings.itime
    
    input:
    file genes
    file plsnet_net from plsnet_network_raw

    output:
    file 'plsnet_network.sf' into plsnet_sf

    script:
    if (targets == '')
    {
      """
      seidr import -F m -r -z -n ${params.plsnet_settings.importname} \
                   -i ${plsnet_net} -g ${genes} -o plsnet_network.sf \
                   -O ${params.plsnet_settings.importcores}
      """
    }
    else
    {
      """
      seidr import -F el -r -z -n ${params.plsnet_settings.importname} \
                   -i ${plsnet_net} -g ${genes} -o plsnet_network.sf \
                   -O ${params.plsnet_settings.importcores}
      """
    }
  }
}