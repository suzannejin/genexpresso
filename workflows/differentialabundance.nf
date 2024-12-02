/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def checkPathParamList = [ params.input ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

def exp_meta = [ "id": params.study_name  ]
if (params.input) { ch_input = Channel.of([ exp_meta, file(params.input, checkIfExists: true) ]) } else { exit 1, 'Input not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { VALIDATE_RNASEQ                                   } from '../subworkflows/local/validate_rnaseq/main'
include { VALIDATE_AFFY_ARRAY                               } from '../subworkflows/local/validate_affy_array/main'
include { VALIDATE_GEO_SOFT_FILE                            } from '../subworkflows/local/validate_geo_soft_file/main'
include { VALIDATE_MAXQUANT                                 } from '../subworkflows/local/validate_maxquant/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ABUNDANCE_DIFFERENTIAL_FILTER as DIFFERENTIAL     } from '../subworkflows/nf-core/abundance_differential_filter/main'

include { SHINYNGS_VALIDATEFOMCOMPONENTS as VALIDATOR       } from '../modules/nf-core/shinyngs/validatefomcomponents/main'
include { CUSTOM_MATRIXFILTER                               } from '../modules/nf-core/custom/matrixfilter/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DIFFERENTIALABUNDANCE {

    take:
    ch_tools

    main:

    ch_versions = Channel.empty()

    // ----------------------------------------
    // PARSE INPUTS
    // ----------------------------------------

    // Parse the matrix and features

    if (params.study_type == "rnaseq") {
        VALIDATE_RNASEQ(ch_input)
        ch_matrix = VALIDATE_RNASEQ.out.matrix
        ch_features = VALIDATE_RNASEQ.out.features
        ch_versions = ch_versions.mix(VALIDATE_RNASEQ.out.versions)

    } else if (params.study_type == 'affy_array') {
        VALIDATE_AFFY_ARRAY(ch_input)
        ch_matrix = VALIDATE_AFFY_ARRAY.out.matrix
        ch_features = VALIDATE_AFFY_ARRAY.out.features
        ch_versions = ch_versions.mix(VALIDATE_AFFY_ARRAY.out.versions)

    } else if (params.study_type == 'geo_soft_file') {
        VALIDATE_GEO_SOFT_FILE(ch_input)
        ch_matrix = VALIDATE_GEO_SOFT_FILE.out.matrix
        ch_features = VALIDATE_GEO_SOFT_FILE.out.features
        ch_versions = ch_versions.mix(VALIDATE_GEO_SOFT_FILE.out.versions)

    } else if (params.study_type == 'maxquant') {
        VALIDATE_MAXQUANT(ch_input)
        ch_matrix = VALIDATE_MAXQUANT.out.matrix
        ch_features = VALIDATE_MAXQUANT.out.features
        ch_versions = ch_versions.mix(VALIDATE_MAXQUANT.out.versions)

    } else {
        error("Study type not recognized")
    }

    // Parse contrast file

    ch_contrasts = Channel.from([[exp_meta, file(params.contrasts)]])

    // Parse transcript lenght and control features - optional - for normalization

    if (params.transcript_length_matrix) { 
        ch_transcript_lengths = Channel.of([ exp_meta, file(params.transcript_length_matrix, checkIfExists: true)]).first() 
    } else { 
        ch_transcript_lengths = [[],[]] 
    }

    if (params.control_features) {
        ch_control_features = Channel.of([ exp_meta, file(params.control_features, checkIfExists: true)]).first() 
    } else { 
        ch_control_features = [[],[]] 
    }

    // Parse gene sets files - optional - for enrichment analysis

    if (params.gene_sets_files) {
        gene_sets_files = params.gene_sets_files.split(",")
        ch_gene_sets = Channel.of(gene_sets_files).map { file(it, checkIfExists: true) }
    }

    // Parse other files

    report_file = file(params.report_file, checkIfExists: true)
    logo_file = file(params.logo_file, checkIfExists: true)
    css_file = file(params.css_file, checkIfExists: true)
    citations_file = file(params.citations_file, checkIfExists: true)

    // ----------------------------------------
    // VALIDATE INPUT FORMATS
    // ----------------------------------------

    /* It validates the consistency of features and samples annotations with matrices and contrasts
        and returns the files with the proper formats.

        Example:
        If a matrix has a column for gene_id, another column for gene_name, and the
        rest of columns for sample expressions, it will remove the gene_name column and return the
        matrix with only the gene_id column and the sample expressions columns
    */

    VALIDATOR(
        ch_input.join(ch_matrix),
        ch_features,
        ch_contrasts
    )

    ch_matrix = VALIDATOR.out.assays
    ch_input = VALIDATOR.out.sample_meta

    // Split the contrasts up so we can run differential analyses and
    // downstream plots separately.
    // Replace NA strings that might have snuck into the blocking column

    ch_contrasts = VALIDATOR.out.contrasts
        .map{it[1]}
        .splitCsv ( header:true, sep:'\t' )
        .map{
            it.blocking = it.blocking.replace('NA', '')
            if (!it.id){
                it.id = it.values().join('_')
            }
            tuple(it, it.variable, it.reference, it.target)
        }

    // ----------------------------------------
    // FILTER INPUT MATRIX
    // ----------------------------------------

    /* Filter out the lowly expressed features and samples from the matrix */

    CUSTOM_MATRIXFILTER(
        ch_matrix,
        ch_input
    )
    ch_matrix = CUSTOM_MATRIXFILTER.out.filtered

    // ----------------------------------------
    // RUN DIFFERENTIAL ANALYSIS
    // ----------------------------------------

    // combine the abundance matrix with the tools channel, to dictate which methods to run
    ch_input_to_differential = ch_matrix
        .combine(ch_tools)
        .map { meta_matrix, matrix, pathway_name, differential_map, correlation_map, enrichment_map ->
            def meta = meta_matrix.clone() + ['pathway_name' : pathway_name['pathway_name']] + differential_map.clone()
            def fc_threshold = differential_map['differential_min_fold_change'] ? differential_map['differential_min_fold_change'] : params.differential_min_fold_change
            def pvalue_threshold = differential_map['differential_max_qval'] ? differential_map['differential_max_qval'] : params.differential_max_qval
            [ meta, matrix, differential_map['diff_method'], fc_threshold, pvalue_threshold ]
        }

    ch_input_to_differential.view()

    // // run the differential subworkflow
    // DIFFERENTIAL(
    //     ch_input_to_differential,
    //     ch_input,
    //     ch_transcript_lengths,
    //     ch_control_features,
    //     ch_contrasts
    // )

}   

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/