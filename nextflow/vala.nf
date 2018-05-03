#!/usr/bin/env nextflow

expr = file(params.expr)
genes = file(params.genes)
targets = file(params.targets)
if (targets != '' ) {
  targets = '-t ' + targets
}

if (params.clr) {

  process clr {
    cpus params.clr_settings.cores
    memory params.clr_settings.memory
    publishDir params.out + '/networks/clr'
    input:
    file expr
    file genes
    val targets
    output:
    file 'clr_network.tsv' into clr_network_raw

    """
    export OMP_NUM_THREADS=1
    mpirun -x OMP_NUM_THREADS --oversubscribe -np ${params.clr_settings.cores} \
    mi -m CLR -o clr_network.tsv -B ${params.aracne_settings.batchsize} \
       -b ${params.clr_settings.bins} -s ${params.clr_settings.spline} \
       -i ${expr} -g ${genes} ${targets}
    """
  }

  process clr_import {
    validExitStatus 0,3
    cpus params.clr_settings.importcores
    memory params.clr_settings.importmem
    publishDir params.out + '/networks/clr'
    input:
    file genes
    file clr_net from clr_network_raw
    output:
    file 'clr_network.sf' into clr_sf

    script:
    if (targets == '')
    {
      """
      export OMP_NUM_THREADS=${params.clr_settings.importcores}
      seidr import -F lm -u -r -z -n ${params.clr_settings.importname} \
                   -i ${clr_net} -g ${genes} -o clr_network.sf
      """
    }
    else
    {
      """
      export OMP_NUM_THREADS=${params.clr_settings.importcores}
      seidr import -F el -u -r -z -n ${params.clr_settings.importname} \
                   -i ${clr_net} -g ${genes} -o clr_network.sf
      """
    }
  }
}

if (params.aracne) {

  process aracne {
    cpus params.aracne_settings.cores
    memory params.aracne_settings.memory
    publishDir params.out + '/networks/aracne'
    input:
    file expr
    file genes
    val targets
    output:
    file 'aracne_network.tsv' into aracne_network_raw

    """
    export OMP_NUM_THREADS=1
    mpirun -x OMP_NUM_THREADS --oversubscribe -np ${params.aracne_settings.cores} \
    mi -m ARACNE -o aracne_network.tsv -B ${params.aracne_settings.batchsize} \
       -b ${params.aracne_settings.bins} -s ${params.aracne_settings.spline} \
       -i ${expr} -g ${genes} ${targets}
    """
  }

  process aracne_import {
    validExitStatus 0,3
    cpus params.aracne_settings.importcores
    memory params.aracne_settings.importmem
    publishDir params.out + '/networks/aracne'
    input:
    file genes
    file aracne_net from aracne_network_raw
    output:
    file 'aracne_network.sf' into aracne_sf

    script:
    if (targets == '')
    {
      """
      export OMP_NUM_THREADS=${params.aracne_settings.importcores}
      seidr import -F lm -u -r -z -n ${params.aracne_settings.importname} \
                   -i ${aracne_net} -g ${genes} -o aracne_network.sf
      """
    }
    else
    {
      """
      export OMP_NUM_THREADS=${params.aracne_settings.importcores}
      seidr import -F el -u -r -z -n ${params.aracne_settings.importname} \
                   -i ${aracne_net} -g ${genes} -o aracne_network.sf
      """
    }
  }
}

if (params.anova) {

  process anova {
    cpus params.anova_settings.cores
    memory params.anova_settings.memory
    publishDir params.out + '/networks/anova'
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
    validExitStatus 0,3
    cpus params.anova_settings.importcores
    memory params.anova_settings.importmem
    publishDir params.out + '/networks/anova'
    input:
    file genes
    file anova_net from anova_network_raw
    output:
    file 'anova_network.sf' into anova_sf

    script:
    if (targets == '')
    {
      """
      export OMP_NUM_THREADS=${params.anova_settings.importcores}
      seidr import -F lm -u -r -z -n ${params.anova_settings.importname} \
                   -i ${aracne_net} -g ${genes} -o anova_network.sf
      """
    }
    else
    {
      """
      export OMP_NUM_THREADS=${params.anova_settings.importcores}
      seidr import -F el -u -r -z -n ${params.anova_settings.importname} \
                   -i ${aracne_net} -g ${genes} -o anova_network.sf
      """
    }
  }
}

if (params.pearson) {

  process pearson {
    cpus params.pearson_settings.cores
    memory params.pearson_settings.memory
    publishDir params.out + '/networks/pearson'
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
    validExitStatus 0,3
    cpus params.pearson_settings.importcores
    memory params.pearson_settings.importmem
    publishDir params.out + '/networks/pearson'
    input:
    file genes
    file pearson_net from pearson_network_raw
    output:
    file 'pearson_network.sf' into pearson_sf

    script:
    if (targets == '')
    {
      """
      export OMP_NUM_THREADS=${params.pearson_settings.importcores}
      seidr import -F lm -A -u -r -z -n ${params.pearson_settings.importname} \
                   -i ${pearson_net} -g ${genes} -o pearson_network.sf
      """
    }
    else
    {
      """
      export OMP_NUM_THREADS=${params.pearson_settings.importcores}
      seidr import -F el -A -u -r -z -n ${params.pearson_settings.importname} \
                   -i ${pearson_net} -g ${genes} -o pearson_network.sf
      """
    }
  }
}

if (params.spearman) {

  process spearman {
    cpus params.spearman_settings.cores
    memory params.spearman_settings.memory
    publishDir params.out + '/networks/spearman'
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
    validExitStatus 0,3
    cpus params.spearman_settings.importcores
    memory params.spearman_settings.importmem
    publishDir params.out + '/networks/spearman'
    input:
    file genes
    file spearman_net from spearman_network_raw
    output:
    file 'spearman_network.sf' into spearman_sf

    script:
    if (targets == '')
    {
      """
      export OMP_NUM_THREADS=${params.spearman_settings.importcores}
      seidr import -F lm -A -u -r -z -n ${params.spearman_settings.importname} \
                   -i ${spearman_net} -g ${genes} -o spearman_network.sf
      """
    }
    else
    {
      """
      export OMP_NUM_THREADS=${params.spearman_settings.importcores}
      seidr import -F el -A -u -r -z -n ${params.spearman_settings.importname} \
                   -i ${spearman_net} -g ${genes} -o spearman_network.sf
      """
    }
  }
}

if (params.elnet) {

  process elnet {
    cpus params.elnet_settings.cores
    memory params.elnet_settings.memory
    publishDir params.out + '/networks/elnet'
    input:
    file expr
    file genes
    val targets
    output:
    file 'elnet_network.tsv' into elnet_network_raw

    """
    export OMP_NUM_THREADS=1
    mpirun -x OMP_NUM_THREADS --oversubscribe -np ${params.elnet_settings.cores} \
    el-ensemble -o elnet_network.tsv -b ${params.elnet_settings.batchsize} \
       -l ${params.elnet_settings.min_lambda} -a ${params.elnet_settings.alpha} \
       -n ${params.elnet_settings.nlambda} \
       -X ${params.elnet_settings.max_experiment_size} \
       -x ${params.elnet_settings.min_experiment_size} \
       -P ${params.elnet_settings.max_predictor_size} \
       -p ${params.elnet_settings.min_predictor_size} \
       -e ${params.elnet_settings.ensemble} \
       -i ${expr} -g ${genes} ${targets}
    """
  }

  process elnet_import {
    validExitStatus 0,3
    cpus params.elnet_settings.importcores
    memory params.elnet_settings.importmem
    publishDir params.out + '/networks/elnet'
    input:
    file genes
    file elnet_net from elnet_network_raw
    output:
    file 'elnet_network.sf' into elnet_sf

    script:
    if (targets == '')
    {
      """
      export OMP_NUM_THREADS=${params.elnet_settings.importcores}
      seidr import -F m -r -z -n ${params.elnet_settings.importname} \
                   -i ${elnet_net} -g ${genes} -o elnet_network.sf
      """
    }
    else
    {
      """
      export OMP_NUM_THREADS=${params.elnet_settings.importcores}
      seidr import -F el -r -z -n ${params.elnet_settings.importname} \
                   -i ${elnet_net} -g ${genes} -o elnet_network.sf
      """
    }
  }
}

if (params.svm) {

  process svm {
    cpus params.svm_settings.cores
    memory params.svm_settings.memory
    publishDir params.out + '/networks/svm'
    input:
    file expr
    file genes
    val targets
    output:
    file 'svm_network.tsv' into svm_network_raw

    """
    export OMP_NUM_THREADS=1
    mpirun -x OMP_NUM_THREADS --oversubscribe -np ${params.svm_settings.cores} \
    svm-ensemble -o svm_network.tsv -b ${params.svm_settings.batchsize} \
       -X ${params.svm_settings.max_experiment_size} \
       -x ${params.svm_settings.min_experiment_size} \
       -P ${params.svm_settings.max_predictor_size} \
       -p ${params.svm_settings.min_predictor_size} \
       -e ${params.svm_settings.ensemble} \
       -i ${expr} -g ${genes} ${targets}
    """
  }

  process svm_import {
    validExitStatus 0,3
    cpus params.svm_settings.importcores
    memory params.svm_settings.importmem
    publishDir params.out + '/networks/svm'
    input:
    file genes
    file svm_net from svm_network_raw
    output:
    file 'svm_network.sf' into svm_sf

    script:
    if (targets == '')
    {
      """
      export OMP_NUM_THREADS=${params.svm_settings.importcores}
      seidr import -F m -r -z -n ${params.svm_settings.importname} \
                   -i ${svm_net} -g ${genes} -o svm_network.sf
      """
    }
    else
    {
      """
      export OMP_NUM_THREADS=${params.svm_settings.importcores}
      seidr import -F el -r -z -n ${params.svm_settings.importname} \
                   -i ${svm_net} -g ${genes} -o svm_network.sf
      """
    }
  }
}

if (params.llr) {

  process llr {
    cpus params.llr_settings.cores
    memory params.llr_settings.memory
    publishDir params.out + '/networks/llr'
    input:
    file expr
    file genes
    val targets
    output:
    file 'llr_network.tsv' into llr_network_raw

    """
    export OMP_NUM_THREADS=1
    mpirun -x OMP_NUM_THREADS --oversubscribe -np ${params.llr_settings.cores} \
    llr-ensemble -o llr_network.tsv -b ${params.llr_settings.batchsize} \
       -X ${params.llr_settings.max_experiment_size} \
       -x ${params.llr_settings.min_experiment_size} \
       -P ${params.llr_settings.max_predictor_size} \
       -p ${params.llr_settings.min_predictor_size} \
       -e ${params.llr_settings.ensemble} \
       -i ${expr} -g ${genes} ${targets}
    """
  }

  process llr_import {
    validExitStatus 0,3
    cpus params.llr_settings.importcores
    memory params.llr_settings.importmem
    publishDir params.out + '/networks/llr'
    input:
    file genes
    file llr_net from llr_network_raw
    output:
    file 'llr_network.sf' into llr_sf

    script:
    if (targets == '')
    {
      """
      export OMP_NUM_THREADS=${params.llr_settings.importcores}
      seidr import -F m -r -z -n ${params.llr_settings.importname} \
                   -i ${llr_net} -g ${genes} -o llr_network.sf
      """
    }
    else
    {
      """
      export OMP_NUM_THREADS=${params.llr_settings.importcores}
      seidr import -F el -r -z -n ${params.llr_settings.importname} \
                   -i ${llr_net} -g ${genes} -o llr_network.sf
      """
    }
  }
}

if (params.pcor) {

  process pcor {
    cpus params.pcor_settings.cores
    memory params.pcor_settings.memory
    publishDir params.out + '/networks/pcor'
    input:
    file expr
    file genes
    val targets
    output:
    file 'pcor_network.tsv' into pcor_network_raw

    """
    export OMP_NUM_THREADS=1
    pcor -i ${expr} -g ${genes} -o pcor_network.tsv \
         ${params.pcor_settings.absolute} ${targets}
    """
  }

  process pcor_import {
    validExitStatus 0,3
    cpus params.pcor_settings.importcores
    memory params.pcor_settings.importmem
    publishDir params.out + '/networks/pcor'
    input:
    file genes
    file pcor_net from pcor_network_raw
    output:
    file 'pcor_network.sf' into pcor_sf

    script:
    if (targets == '')
    {
      """
      export OMP_NUM_THREADS=${params.pcor_settings.importcores}
      seidr import -F lm -A -u -r -z -n ${params.pcor_settings.importname} \
                   -i ${pcor_net} -g ${genes} -o pcor_network.sf
      """
    }
    else
    {
      """
      export OMP_NUM_THREADS=${params.pcor_settings.importcores}
      seidr import -F el -A -u -r -z -n ${params.pcor_settings.importname} \
                   -i ${pcor_net} -g ${genes} -o pcor_network.sf
      """
    }
  }
}

if (params.narromi) {

  process narromi {
    cpus params.narromi_settings.cores
    memory params.narromi_settings.memory
    publishDir params.out + '/networks/narromi'
    input:
    file expr
    file genes
    val targets
    output:
    file 'narromi_network.tsv' into narromi_network_raw

    """
    export OMP_NUM_THREADS=1
    mpirun -x OMP_NUM_THREADS --oversubscribe -np ${params.narromi_settings.cores} \
    narromi -o narromi_network.tsv -b ${params.narromi_settings.batchsize} \
       -a ${params.narromi_settings.alpha} \
       -m ${params.narromi_settings.method} \
       -i ${expr} -g ${genes} ${targets}
    """
  }

  process narromi_import {
    validExitStatus 0,3
    cpus params.narromi_settings.importcores
    memory params.narromi_settings.importmem
    publishDir params.out + '/networks/narromi'
    input:
    file genes
    file narromi_net from narromi_network_raw
    output:
    file 'narromi_network.sf' into narromi_sf

    script:
    if (targets == '')
    {
      """
      export OMP_NUM_THREADS=${params.narromi_settings.importcores}
      seidr import -F m -r -z -n ${params.narromi_settings.importname} \
                   -i ${narromi_net} -g ${genes} -o narromi_network.sf
      """
    }
    else
    {
      """
      export OMP_NUM_THREADS=${params.narromi_settings.importcores}
      seidr import -F el -r -z -n ${params.narromi_settings.importname} \
                   -i ${narromi_net} -g ${genes} -o narromi_network.sf
      """
    }
  }
}

if (params.tigress) {

  process tigress {
    cpus params.tigress_settings.cores
    memory params.tigress_settings.memory
    publishDir params.out + '/networks/tigress'
    input:
    file expr
    file genes
    val targets
    output:
    file 'tigress_network.tsv' into tigress_network_raw

    """
    export OMP_NUM_THREADS=1
    mpirun -x OMP_NUM_THREADS --oversubscribe -np ${params.tigress_settings.cores} \
    tigress -o tigress_network.tsv -b ${params.tigress_settings.batchsize} \
       -n ${params.tigress_settings.nlambda} \
       -l ${params.tigress_settings.min_lambda} \
       ${params.tigress_settings.scale} \
       -i ${expr} -g ${genes} ${targets}
    """
  }

  process tigress_import {
    validExitStatus 0,3
    cpus params.tigress_settings.importcores
    memory params.tigress_settings.importmem
    publishDir params.out + '/networks/tigress'
    input:
    file genes
    file tigress_net from tigress_network_raw
    output:
    file 'tigress_network.sf' into tigress_sf

    script:
    if (targets == '')
    {
      """
      export OMP_NUM_THREADS=${params.tigress_settings.importcores}
      seidr import -F m -r -z -n ${params.tigress_settings.importname} \
                   -i ${tigress_net} -g ${genes} -o tigress_network.sf
      """
    }
    else
    {
      """
      export OMP_NUM_THREADS=${params.tigress_settings.importcores}
      seidr import -F el -r -z -n ${params.tigress_settings.importname} \
                   -i ${tigress_net} -g ${genes} -o tigress_network.sf
      """
    }
  }
}

if (params.genie3) {

  process genie3 {
    cpus params.genie3_settings.cores
    memory params.genie3_settings.memory
    publishDir params.out + '/networks/genie3'
    input:
    file expr
    file genes
    val targets
    output:
    file 'genie3_network.tsv' into genie3_network_raw

    """
    export OMP_NUM_THREADS=1
    mpirun -x OMP_NUM_THREADS --oversubscribe -np ${params.genie3_settings.cores} \
    genie3 -o genie3_network.tsv -b ${params.genie3_settings.batchsize} \
       -p ${params.genie3_settings.min_prop} \
       -a ${params.genie3_settings.alpha} \
       -N ${params.genie3_settings.min_node_size} \
       -m ${params.genie3_settings.mtry} \
       -n ${params.genie3_settings.ntree} \
       ${params.genie3_settings.scale} \
       -i ${expr} -g ${genes} ${targets}
    """
  }

  process genie3_import {
    validExitStatus 0,3
    cpus params.genie3_settings.importcores
    memory params.genie3_settings.importmem
    publishDir params.out + '/networks/genie3'
    input:
    file genes
    file genie3_net from genie3_network_raw
    output:
    file 'genie3_network.sf' into genie3_sf

    script:
    if (targets == '')
    {
      """
      export OMP_NUM_THREADS=${params.genie3_settings.importcores}
      seidr import -F m -r -z -n ${params.genie3_settings.importname} \
                   -i ${genie3_net} -g ${genes} -o genie3_network.sf
      """
    }
    else
    {
      """
      export OMP_NUM_THREADS=${params.genie3_settings.importcores}
      seidr import -F el -r -z -n ${params.genie3_settings.importname} \
                   -i ${genie3_net} -g ${genes} -o genie3_network.sf
      """
    }
  }
}

if (params.plsnet) {

  process plsnet {
    cpus params.plsnet_settings.cores
    memory params.plsnet_settings.memory
    publishDir params.out + '/networks/plsnet'
    input:
    file expr
    file genes
    val targets
    output:
    file 'plsnet_network.tsv' into plsnet_network_raw

    """
    export OMP_NUM_THREADS=1
    mpirun --oversubscribe -np ${params.plsnet_settings.cores} \
    plsnet -o plsnet_network.tsv -b ${params.plsnet_settings.batchsize} \
       -c ${params.plsnet_settings.components} \
       -p ${params.plsnet_settings.predictors} \
       -e ${params.plsnet_settings.ensemble} \
       ${params.plsnet_settings.scale} \
       -i ${expr} -g ${genes} ${targets}
    """
  }

  process plsnet_import {
    validExitStatus 0,3
    cpus params.plsnet_settings.importcores
    memory params.plsnet_settings.importmem
    publishDir params.out + '/networks/plsnet'
    input:
    file genes
    file plsnet_net from plsnet_network_raw
    output:
    file 'plsnet_network.sf' into plsnet_sf

    script:
    if (targets == '')
    {
      """
      export OMP_NUM_THREADS=${params.plsnet_settings.importcores}
      seidr import -F m -r -z -n ${params.plsnet_settings.importname} \
                   -i ${plsnet_net} -g ${genes} -o plsnet_network.sf
      """
    }
    else
    {
      """
      export OMP_NUM_THREADS=${params.plsnet_settings.importcores}
      seidr import -F el -r -z -n ${params.plsnet_settings.importname} \
                   -i ${plsnet_net} -g ${genes} -o plsnet_network.sf
      """
    }
  }
}