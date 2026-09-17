#!/usr/bin/env nextflow
/*
 * Nextflow port of run.py.
 *
 * Each step below is a thin wrapper around the corresponding scripts/N.*.py
 * script. The scripts themselves read paths out of the defaults/experiment
 * config files and read/write data on the shared filesystem (not through
 * Nextflow channels), so steps are chained purely to enforce execution
 * order 0 -> 7, matching the sequential behavior of run.py.
 */

nextflow.enable.dsl = 2

import groovy.json.JsonSlurper

params.defaults_config   = "config/defaults.yaml"
params.experiment_config = "config/experiment.json"
params.recipe_dir        = "${projectDir}"

def defaultsConfig   = file(params.defaults_config).toAbsolutePath()
def experimentConfig = file(params.experiment_config).toAbsolutePath()

if (!defaultsConfig.exists()) {
    exit 1, "Defaults config not found: ${defaultsConfig}"
}
if (!experimentConfig.exists()) {
    exit 1, "Experiment config not found: ${experimentConfig}"
}

def experiment = new JsonSlurper().parse(experimentConfig)

// Directories that must be visible inside a container at their original host
// path, since the scripts read/write absolute paths taken from the configs
// rather than through Nextflow's own file staging. Only used when a
// container engine (e.g. -profile docker) is enabled; harmless otherwise.
def mountDirs = [
    params.recipe_dir,
    defaultsConfig.parent.toString(),
    experimentConfig.parent.toString(),
] as Set
if (experiment.file_location) {
    mountDirs << file(experiment.file_location).toAbsolutePath().toString()
}
if (experiment.match_to_library && experiment.library_location) {
    mountDirs << file(experiment.library_location).toAbsolutePath().parent.toString()
}
def dockerMounts = mountDirs.collect { "-v ${it}:${it}" }.join(' ')

process IMAGE_QC {
    tag "0.image-qc"
    containerOptions dockerMounts

    input:
    val ready
    val defaults_config
    val experiment_config

    output:
    val true

    script:
    """
    cd ${params.recipe_dir}
    python ${params.recipe_dir}/scripts/0.image-qc.py ${defaults_config} ${experiment_config}
    """
}

process PROCESS_SBS {
    tag "1.process-SBS"
    containerOptions dockerMounts

    input:
    val ready
    val defaults_config
    val experiment_config

    output:
    val true

    script:
    """
    cd ${params.recipe_dir}
    python ${params.recipe_dir}/scripts/1.process-SBS.py ${defaults_config} ${experiment_config}
    """
}

process MERGE_SINGLE_CELLS {
    tag "2.merge-single-cells"
    containerOptions dockerMounts

    input:
    val ready
    val defaults_config
    val experiment_config

    output:
    val true

    script:
    """
    cd ${params.recipe_dir}
    python ${params.recipe_dir}/scripts/2.merge-single-cells.py ${defaults_config} ${experiment_config}
    """
}

process SUMMARIZE_SBS {
    tag "3.summarize-SBS"
    containerOptions dockerMounts

    input:
    val ready
    val defaults_config
    val experiment_config

    output:
    val true

    script:
    """
    cd ${params.recipe_dir}
    python ${params.recipe_dir}/scripts/3.summarize-SBS.py ${defaults_config} ${experiment_config}
    """
}

process AGGREGATE {
    tag "4.aggregate"
    containerOptions dockerMounts

    input:
    val ready
    val defaults_config
    val experiment_config

    output:
    val true

    script:
    """
    cd ${params.recipe_dir}
    python ${params.recipe_dir}/scripts/4.aggregate.py ${defaults_config} ${experiment_config}
    """
}

process NORMALIZE {
    tag "5.normalize"
    containerOptions dockerMounts

    input:
    val ready
    val defaults_config
    val experiment_config

    output:
    val true

    script:
    """
    cd ${params.recipe_dir}
    python ${params.recipe_dir}/scripts/5.normalize.py ${defaults_config} ${experiment_config}
    """
}

process FEATURE_SELECT {
    tag "6.feature-select"
    containerOptions dockerMounts

    input:
    val ready
    val defaults_config
    val experiment_config

    output:
    val true

    script:
    """
    cd ${params.recipe_dir}
    python ${params.recipe_dir}/scripts/6.feature-select.py ${defaults_config} ${experiment_config}
    """
}

process EXPLORE {
    tag "7.explore"
    containerOptions dockerMounts

    input:
    val ready
    val defaults_config
    val experiment_config

    output:
    val true

    script:
    """
    cd ${params.recipe_dir}
    python ${params.recipe_dir}/scripts/7.explore.py ${defaults_config} ${experiment_config}
    """
}

workflow {
    def defaults_config   = defaultsConfig.toString()
    def experiment_config = experimentConfig.toString()

    def ready = Channel.value(true)

    if (experiment.perform_process_qc) {
        ready = IMAGE_QC(ready, defaults_config, experiment_config)
    }
    if (experiment.perform_process_SBS) {
        ready = PROCESS_SBS(ready, defaults_config, experiment_config)
    }
    if (experiment.perform_process_single_cells) {
        ready = MERGE_SINGLE_CELLS(ready, defaults_config, experiment_config)
    }
    if (experiment.perform_summarize_SBS) {
        ready = SUMMARIZE_SBS(ready, defaults_config, experiment_config)
    }
    if (experiment.perform_aggregate) {
        ready = AGGREGATE(ready, defaults_config, experiment_config)
    }
    if (experiment.perform_normalize) {
        ready = NORMALIZE(ready, defaults_config, experiment_config)
    }
    if (experiment.perform_feature_select) {
        ready = FEATURE_SELECT(ready, defaults_config, experiment_config)
    }
    if (experiment.perform_explore) {
        ready = EXPLORE(ready, defaults_config, experiment_config)
    }
}
