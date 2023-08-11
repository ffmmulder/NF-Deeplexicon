#!/usr/bin/env nextflow
nextflow.enable.dsl=2

def helpMessage() {
    // Log colors ANSI codes
    c_reset = "\033[0m";
    c_yellow = "\033[0;33m";
    c_blue = "\033[0;34m";

    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run NF-Deeplexicon/main.nf --fast5_path <fast5_dir> --out_dir <output_dir> -c <path/to/analysis.config --email <email address>
${c_blue}    Mandatory arguments: ${c_reset}
${c_yellow}        --fast5_path [str] ${c_reset}             Path to a directory containing fast5 files.
                                                    Files should be named in the following format: <samplename>_pass_<id>_<chunk>.fast5
${c_yellow}        --out_dir [str] ${c_reset}                 The output directory where the results will be saved.
${c_yellow}        --email [str] ${c_reset}                   The email address to send workflow summary and MultiQC report to.
${c_yellow}        --account [str] ${c_reset}                   The cluster account to be used for the computation and jobs.

${c_blue}    Standard options: ${c_reset}
        --profile [str]                 Configuration profile to use, leave empty to run locally.
                                              Available: slurm, SGE, singularity.
        --fastq_path [str] ${c_reset}             Path to a directory containing passed fastq files. If not specified the parent path of '--fast5_path' is assumed with 'fast5_pass' being substituted for 'fastq_pass'
    """.stripIndent()
}

// Show help message and exit.
if(params.help){
  helpMessage()
  exit 0
}
// Minimal required parameters.
if (!params.out_dir) {
   exit 1, "Output directory not found, please provide the correct path! (--out_dir)"
}
if (!params.fast5_path) {
   exit 1, "Fast5 directory not found, please provide the correct path! (--fast5_path)"
}

if (!params.fastq_path) {
   params.fastq_path=params.fast5_path.replace("fast5_pass", "fastq_pass")
   print "No fastq path provided, assuming it is ${params.fastq_path}..."
}

// Utils modules
include { ExportParams as Workflow_ExportParams } from params.nextflowmodules_path+'/Utils/workflow.nf'

include { extractAllFast5FromDir } from params.nextflowmodules_path+'/Utils/fast5.nf' params(params)  
include { extractAllNanoporeFastqFromDir } from params.nextflowmodules_path+'/Utils/fastq.nf' params(params)  

// Deeplexicon modules
include { dmux as DeeplexiconDmux } from params.nextflowmodules_path+'/Deeplexicon/1.2.0/dmux.nf' //params(gatk_path: "$params.gatk_path", genome:"$params.genome", optional: "--intervals $params.dxtracks_path/$params.fingerprint_target --output_mode EMIT_ALL_SITES")
include { split as DeeplexiconSplit } from params.nextflowmodules_path+'/Deeplexicon/1.2.0/split.nf' // params(genome:"$params.genome", optional: '-c 100 -M')
include { deeplexicon_concatTsv as DeeplexiconConcatTsv } from params.nextflowmodules_path+'/Deeplexicon/1.2.0/tsv.nf'
include { deeplexicon_concatFastq } from params.nextflowmodules_path+'/Deeplexicon/1.2.0/fastq.nf'

workflow {

    def sample = "dummy_sample"
    def fast5_by_id
    def fastq_files
    def run_name = ""
    if ( params.run_name) {
        run_name = params.run_name
    }

    //dmux fast5  files
    if ( ! params.tsv_path) {
        //collect fast5 files, map based on first element (=id)
        def fast5_files = extractAllFast5FromDir(params.fast5_path).map { [it[0],it[1],it[2]]}
        
        //first dmux, creates one tsv per fast5 (fastest due to parallel processing)
        DeeplexiconDmux(fast5_files)

        //concat resuulting tsv files per id
        def tsv_files = DeeplexiconDmux.out.groupTuple(by:[0]).map{id,chunk,tsv -> [id, chunk, tsv]}
        DeeplexiconConcatTsv(tsv_files)
    } else {
        println "Skipping fast5 demux, change to use TSV path: ${params.tsv_path}"

        include { extractAllTsvFromDir } from params.nextflowmodules_path+'/Utils/tsv.nf'
        include { deeplexicon_copyTsv } from params.nextflowmodules_path+'/Deeplexicon/1.2.0/tsv.nf'

        def tsv_files = extractAllTsvFromDir(params.tsv_path+"/*/").map { [it[0],it[1],it[2]]}
        deeplexicon_copyTsv(tsv_files)
        def tsv_by_id = deeplexicon_copyTsv.out.groupTuple(by:[0]).map{id,chunk,tsv -> [id, chunk, tsv]}

        DeeplexiconConcatTsv(tsv_by_id)
    }
    //collect fastq files, map based on first token (=id)
    fastq_files = extractAllNanoporeFastqFromDir(params.fastq_path).map { [it[0],it[1],it[2]]}

    //Also concat fastq files into one for each (sample)id
    def fastq_by_id = fastq_files.groupTuple(by:[0]).map{id,chunk,fastq -> [id, chunk, fastq]}
    deeplexicon_concatFastq(fastq_by_id)

    tsv_fastq_by_id = DeeplexiconConcatTsv.out.join(deeplexicon_concatFastq.out, by:[0]).view()
    //split concatenated fastq into the seperate fastqs based on found barcodes
    DeeplexiconSplit(tsv_fastq_by_id)

    // Create log files: Repository versions and Workflow params
    VersionLog()
    Workflow_ExportParams()
}

// Workflow completion notification
workflow.onComplete {
    // HTML Template
    def template = new File("$baseDir/assets/workflow_complete.html")
    def binding = [
        runName: run_name,
        workflow: workflow
    ]
    def engine = new groovy.text.GStringTemplateEngine()
    def email_html = engine.createTemplate(template).make(binding).toString()

    // Send email
    if (workflow.success) {
        def subject = "Deeplexicon Workflow Successful: ${run_name}"
        sendMail(to: params.email.trim(), subject: subject, body: email_html)
//        sendMail(to: params.email.trim(), subject: subject, body: email_html, attach: "${params.outdir}/log/${analysis_id}_multiqc_report.html")
    } else {
        def subject = "Deeplexicon Workflow Failed: ${run_name}"
        sendMail(to: params.email.trim(), subject: subject, body: email_html)
    }
}

//from other pipeline, modify for this
process VersionLog {
    // Custom process to log repository versions
    tag {"VersionLog ${analysis_id}"}
    label 'VersionLog'
    shell = ['/bin/bash', '-eo', 'pipefail']
    cache = false  //Disable cache to force a new version log when restarting the workflow.

    output:
        path('repository_version.log', emit: log_file)

    script:
        """
        echo 'Deeplexicon' > repository_version.log
        git --git-dir=${workflow.projectDir}/.git log --pretty=oneline --decorate -n 2 >> repository_version.log

        echo 'NextflowModules' >> repository_version.log
        git --git-dir=${params.nextflowmodules_path}/.git log --pretty=oneline --decorate -n 2 >> repository_version.log
        """
}