nextflow.enable.dsl=2

/*
 * AF3 + OTR metrics + clean OTR scoring pipeline.
 *
 * This keeps the original main.nf AF3 data/inference stages, then runs:
 *   1. af3_otr_pipeline.py
 *   2. compute_otr_clean.py
 */

params.inputs    = null
params.out_base  = "/home/biolab/Projects/DesCARTes_wd/results/OTRtest/af3_outputs_99"
params.af3_root  = null
params.db_dir    = null
params.model_dir = null
params.c_db_dir    = "/root/public_databases"
params.c_model_dir = "/root/models"

params.run_af3 = true

params.report_dir = "${projectDir}"
params.metrics_csv = "af3_otr_metrics_99_freesasa_cdr.csv"
params.clean_output_csv = "otr_clean_scored_99.csv"
params.summary_csv = "otr_clean_summary_99.csv"
params.prodigy_cache = "prodigy_delta_g_cache_99.csv"

params.chain_roles = "${projectDir}/../data/chain_roles_99.csv"
params.cdr_ranges = "${projectDir}/../data/cdr_ranges_kabat_default_99.csv"
params.affinity_csv = null

params.use_freesasa = true
params.recursive = true
params.run_prodigy_if_missing = true
params.prodigy_bin = "prodigy"
params.prodigy_timeout = 120
params.structure_col = "job_dir"
params.group_column = null
params.use_bsa_penalty = false

params.contact_cutoff = 4.5
params.proxy_variable_n = 120
params.pae_good = 5.0
params.pae_bad = 12.0
params.affinity_good_dg = -15.0
params.affinity_bad_dg = 0.0

params.help = false

def boolParam(value) {
    if (value instanceof Boolean) {
        return value
    }
    return value?.toString()?.toLowerCase() in ['1', 'true', 'yes', 'y']
}

def shellQuote(value) {
    return "'" + value.toString().replace("'", "'\"'\"'") + "'"
}

def basename(value) {
    return file(value.toString()).getName()
}

def resolveProjectPath(value) {
    if (!value) {
        return null
    }
    def text = value.toString()
    return text.startsWith('/') ? text : "${projectDir}/${text}"
}

def metricCsvName = basename(params.metrics_csv)
def metricReadmeName = metricCsvName.replaceFirst(/\.csv$/, '.README.txt')
def cleanCsvName = basename(params.clean_output_csv)
def summaryCsvName = params.summary_csv ? basename(params.summary_csv) : null
def prodigyCacheName = params.prodigy_cache ? basename(params.prodigy_cache) : null
def af3Root = params.af3_root ?: params.out_base

workflow {
    if (boolParam(params.help)) {
        log.info """
        AF3 OTR clean pipeline

        Full run:
          nextflow run af3_otr_clean.nf \\
            --inputs '/home/biolab/Projects/DesCARTes_wd/data/phase2_OTRtestPDBbind/*_json/*.json' \\
            --db_dir /mnt/ssd0/alphafold3_db \\
            --model_dir /home/biolab/alphafold3/models \\
            -resume

        Score existing AF3 outputs only:
          nextflow run af3_otr_clean.nf \\
            --run_af3 false \\
            --db_dir /mnt/ssd0/alphafold3_db \\
            --model_dir /home/biolab/alphafold3/models \\
            -resume

        Key outputs are published to:
          ${params.report_dir}
        """.stripIndent()
        return
    }

    Channel
        .value(file("${projectDir}/af3_otr_pipeline.py"))
        .set { af3_otr_script_ch }

    Channel
        .value(file("${projectDir}/compute_otr_clean.py"))
        .set { compute_clean_script_ch }

    if (boolParam(params.run_af3)) {
        if (!params.inputs) {
            error "Missing --inputs. Provide input JSON glob or run scoring only with --run_af3 false."
        }

        Channel
            .fromPath(params.inputs, checkIfExists: true)
            .map { json_file ->
                def id = json_file.baseName
                def parent_dir = json_file.parent.name
                def target_subdir = parent_dir.replace('input', 'output')
                tuple(id, target_subdir, json_file)
            }
            .set { input_ch }

        data_ch = AF3_DATA_PIPELINE(input_ch)
        AF3_INFERENCE(data_ch)
        SYNC_AF3_OUTPUTS(AF3_INFERENCE.out.job_outputs.collect())
        metrics_trigger_ch = SYNC_AF3_OUTPUTS.out.done
    } else {
        Channel
            .value("using_existing_af3_outputs")
            .set { metrics_trigger_ch }
    }

    AF3_OTR_METRICS(metrics_trigger_ch, af3_otr_script_ch)
    OTR_CLEAN_SCORE(AF3_OTR_METRICS.out.metrics_csv, compute_clean_script_ch)
}

process AF3_DATA_PIPELINE {
    tag { id }

    input:
        tuple val(id), val(sub_dir), path(json_file)

    output:
        tuple val(id), val(sub_dir), path("${id}/**/*_data.json")

    script:
    """
    set -euo pipefail
    mkdir -p ${id}

    python /app/alphafold/run_alphafold.py \\
        --norun_inference \\
        --json_path ${json_file} \\
        --output_dir ${id} \\
        --db_dir ${params.c_db_dir} \\
        --model_dir ${params.c_model_dir}
    """
}

process AF3_INFERENCE {
    tag { id }

    input:
        tuple val(id), val(sub_dir), path(data_json)

    output:
        tuple val(id), val(sub_dir), path("*"), emit: job_outputs

    script:
    """
    set -euo pipefail
    export XLA_PYTHON_CLIENT_PREALLOCATE=false

    python /app/alphafold/run_alphafold.py \\
        --norun_data_pipeline \\
        --json_path ${data_json} \\
        --output_dir . \\
        --db_dir ${params.c_db_dir} \\
        --model_dir ${params.c_model_dir} \\
        --num_diffusion_samples 1
    """
}

process SYNC_AF3_OUTPUTS {
    tag "sync_af3_outputs"
    container null

    input:
        val completed_jobs

    output:
        path "af3_outputs.synced", emit: done

    script:
    def copyCommands = completed_jobs.collect { row ->
        def id = row[0].toString()
        def sub_dir = row[1].toString()
        def outputs = row[2] instanceof Collection ? row[2] : [row[2]]
        def dest = "${params.out_base}/${sub_dir}/${id}"
        def lines = ["mkdir -p ${shellQuote(dest)}"]
        lines.addAll(outputs.collect { out_path ->
            "cp -a ${shellQuote(out_path)} ${shellQuote(dest + '/')}"
        })
        lines.join('\n')
    }.join('\n')

    """
    set -euo pipefail
    mkdir -p ${shellQuote(params.out_base)}
    ${copyCommands}
    touch af3_outputs.synced
    """
}

process AF3_OTR_METRICS {
    tag "af3_otr_metrics"
    container null
    publishDir "${params.report_dir}", mode: 'copy', overwrite: true

    input:
        val trigger
        path af3_otr_script

    output:
        path "${metricCsvName}", emit: metrics_csv
        path "${metricReadmeName}", emit: metrics_readme

    script:
    def chainRolesPath = resolveProjectPath(params.chain_roles)
    def cdrRangesPath = resolveProjectPath(params.cdr_ranges)
    def affinityPath = resolveProjectPath(params.affinity_csv)
    def chainRolesOpt = chainRolesPath ? "--chain_roles ${shellQuote(chainRolesPath)}" : ""
    def cdrRangesOpt = cdrRangesPath ? "--cdr_ranges ${shellQuote(cdrRangesPath)}" : ""
    def affinityOpt = affinityPath ? "--affinity_csv ${shellQuote(affinityPath)}" : ""
    def freesasaOpt = boolParam(params.use_freesasa) ? "--use_freesasa" : ""
    def recursiveOpt = boolParam(params.recursive) ? "--recursive" : ""

    """
    set -euo pipefail

    python ${af3_otr_script} \\
        --af3_root ${shellQuote(af3Root)} \\
        --output_csv ${shellQuote(metricCsvName)} \\
        ${chainRolesOpt} \\
        ${cdrRangesOpt} \\
        ${affinityOpt} \\
        --contact_cutoff ${params.contact_cutoff} \\
        --proxy_variable_n ${params.proxy_variable_n} \\
        --pae_good ${params.pae_good} \\
        --pae_bad ${params.pae_bad} \\
        --affinity_good_dg ${params.affinity_good_dg} \\
        --affinity_bad_dg ${params.affinity_bad_dg} \\
        ${freesasaOpt} \\
        ${recursiveOpt}
    """
}

process OTR_CLEAN_SCORE {
    tag "otr_clean_score"
    container null
    publishDir "${params.report_dir}", mode: 'copy', overwrite: true

    input:
        path metrics_csv
        path compute_clean_script

    output:
        path "${cleanCsvName}", emit: scored_csv
        path "${summaryCsvName}", optional: true, emit: summary_csv
        path "${prodigyCacheName}", optional: true, emit: prodigy_cache

    script:
    def affinityPath = resolveProjectPath(params.affinity_csv)
    def cachePath = resolveProjectPath(params.prodigy_cache)
    def affinityOpt = affinityPath ? "--affinity_csv ${shellQuote(affinityPath)}" : ""
    def summaryOpt = summaryCsvName ? "--summary_csv ${shellQuote(summaryCsvName)}" : ""
    def groupOpt = params.group_column ? "--group_column ${shellQuote(params.group_column)}" : ""
    def bsaPenaltyOpt = boolParam(params.use_bsa_penalty) ? "--use_bsa_penalty" : ""
    def runProdigyOpt = boolParam(params.run_prodigy_if_missing) ? "--run_prodigy_if_missing" : ""
    def cacheOpt = prodigyCacheName ? "--prodigy_cache ${shellQuote(prodigyCacheName)}" : ""
    def cachePrep = cachePath ? """
    if [ -f ${shellQuote(cachePath)} ]; then
        cp ${shellQuote(cachePath)} ${shellQuote(prodigyCacheName)}
    fi
    """.stripIndent() : ""

    """
    set -euo pipefail

    ${cachePrep}

    python ${compute_clean_script} \\
        --metrics_csv ${metrics_csv} \\
        --output_csv ${shellQuote(cleanCsvName)} \\
        ${summaryOpt} \\
        ${affinityOpt} \\
        ${groupOpt} \\
        ${bsaPenaltyOpt} \\
        ${runProdigyOpt} \\
        --prodigy_bin ${shellQuote(params.prodigy_bin)} \\
        ${cacheOpt} \\
        --structure_col ${shellQuote(params.structure_col)} \\
        --prodigy_timeout ${params.prodigy_timeout}
    """
}
