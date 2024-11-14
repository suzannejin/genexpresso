/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def checkPathParamList = [ params.input ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
def exp_meta = [ "id": params.study_name  ]
if (params.input) { ch_samplesheet = Channel.of([ exp_meta, file(params.input, checkIfExists: true) ]) } else { exit 1, 'Input samplesheet not specified!' }

// abundance matrix channel
abundance_file = file(params.matrix, checkIfExists: true)
ch_abundance = Channel.of([ exp_meta, abundance_file])

// features channel
matrix_as_anno_filename = "${workflow.workDir}/${abundance_file.getBaseName()}_as_anno.${abundance_file.getExtension()}"
ch_features = ch_abundance
    .map{ meta, matrix ->
        matrix_copy = file(matrix_as_anno_filename)
        matrix_copy.exists() && matrix.getText().md5().equals(matrix_copy.getText().md5()) ?: matrix.copyTo(matrix_as_anno_filename)
        [ meta, file(matrix_as_anno_filename) ]
    }

// contrasts file channel
ch_contrasts = Channel.from([[exp_meta, file(params.contrasts)]])

// optional files channels
if (params.transcript_length_matrix) { ch_transcript_lengths = Channel.of([ exp_meta, file(params.transcript_length_matrix, checkIfExists: true)]).first() } else { ch_transcript_lengths = [[],[]] }
if (params.control_features) { ch_control_features = Channel.of([ exp_meta, file(params.control_features, checkIfExists: true)]).first() } else { ch_control_features = [[],[]] }

// gene sets
gene_sets_files = params.gene_sets_files.split(",")
ch_gene_sets = Channel.of(gene_sets_files).map { file(it, checkIfExists: true) }

//other files
report_file = file(params.report_file, checkIfExists: true)
logo_file = file(params.logo_file, checkIfExists: true)
css_file = file(params.css_file, checkIfExists: true)
citations_file = file(params.citations_file, checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { samplesheetToList } from 'plugin/nf-schema'

def preprocess_channel(ch_input, method_type, method_name) {
    /* This function filters the input channel by method,
        and adds the pathway_name and method arguments to the meta.
        Then the method arguments can be used in the analysis steps through ext.args

        Example:
        ch_input = [ [pathway_name: 'deseq2_with_gsea', diff_method: 'deseq2', args_diff: '--deseq2_param1: x', enr_method: 'gsea', args_enr: '']
                    [id: 'data1'],
                    data
                    ]
        method_type = 'diff'
        method_name = 'deseq2'
        return [ [id: 'data1', pathway_name: 'deseq2_with_gsea', args: '--deseq2_param1: x'], data ]
    */

    method_field_name = method_type + "_method"
    args_field_name = "args_" + method_type

    return ch_input
            .filter{ it[0][method_field_name] == method_name }
            .map { it ->
                def tools = it[0]
                def meta = it[1] + [ 'pathway_name': tools[pathway_name], 'args': tools[args_field_name] ]
                def data = it[2..-1]
                return [meta, data].flatten()
            }
}

def postprocess_channel(ch_input, ch_tools, method_type, method_name) {
    /* This function add the tools info back to the input channel.
        Then the returned channel can be used in the next analysis step.

        Example:
        ch_input = [ [id: 'data1', pathway_name: 'deseq2_with_gsea', args: '--deseq2_param1: x'], data ]
        method_type = "diff"
        method_name = "deseq2"
        ch_tools = [ [pathway_name: 'deseq2_with_gsea', diff_method: 'deseq2', args_diff: '--deseq2_param1: x', enr_method: 'gsea', args_enr: ''] ]
        return [ [pathway_name: 'deseq2_with_gsea', diff_method: 'deseq2', args_diff: '--deseq2_param1: x', enr_method: 'gsea', args_enr: ''],
                [id: 'data1'],
                deseq2_results
                ]
    */

    method_field_name = method_type + "_method"
    args_field_name = "args_" + method_type

    return ch_tools
            .filter { it[0][method_field_name] == method_name }
            .combine(ch_input)
            .filter { it[0]['pathway_name'] == it[1]['pathway_name'] }
            .map { it ->
                def meta = it[1] - ['pathway_name': it[1]['pathway_name'], 'args': it[1]['args']]
                return [it[0], meta, it[2..-1]].flatten()
            }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { SHINYNGS_VALIDATEFOMCOMPONENTS as VALIDATOR       } from '../modules/nf-core/shinyngs/validatefomcomponents/main'
include { CUSTOM_MATRIXFILTER                               } from '../modules/nf-core/custom/matrixfilter/main'
include { GPROFILER2_GOST                                   } from '../modules/nf-core/gprofiler2/gost/main'

include { ABUNDANCE_DIFFERENTIAL_FILTER as DIFF_DESEQ2      } from '../subworkflows/nf-core/abundance_differential_filter/main'
include { ABUNDANCE_DIFFERENTIAL_FILTER as DIFF_LIMMA       } from '../subworkflows/nf-core/abundance_differential_filter/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DIFFERENTIALABUNDANCE {

    main:

    ch_results_genewise = Channel.empty()
    ch_results_genewise_filtered = Channel.empty()
    ch_versions = Channel.empty()

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
        ch_samplesheet.join(ch_abundance),
        ch_features,
        ch_contrasts
    )

    ch_abundance = VALIDATOR.out.assays
    ch_samplesheet = VALIDATOR.out.sample_meta

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
        ch_abundance,
        ch_samplesheet
    )
    ch_abundance = CUSTOM_MATRIXFILTER.out.filtered

    // ----------------------------------------
    // PREPARE TOOLS CHANNEL
    // ----------------------------------------

    /* Convert the toolsheet.csv in a channel with the proper format

        This tools channel is used to define which methods will be performed in each analysis step.
        And how to combine the different outputs from different methods at different steps.
        The schema functionalities are used for the parsing.

        The tools channel is a dictionary with the following keys (corresponding to the toolsheet.csv columns):
            - pathway_name
            - diff_method
            - args_diff
            - cor_method
            - args_cor
            - enr_method
            - args_enr
    */

    // parse the toolsheet.csv file
    ch_tools = Channel.fromList(samplesheetToList(params.tools, './assets/schema_tools.json'))

    // Filter the tools to the pathway(s) of interest (as indicated by --pathway flag),
    // or run everything if requested (if --pathway all)
    if (params.pathway == "all") {
        ch_tools
            .set{ ch_tools }
    } else {
        ch_tools
            .filter{
                it["pathway_name"] in params.pathway.tokenize(',')
            }
            .set{ ch_tools }
    }

    //  Prepare the inputs for the analysis section
    // The inputs are combined with the tools information to be used in the analysis steps.

    ch_samplesheet_with_method = ch_tools.combine(ch_samplesheet)
    ch_abundance_with_method = ch_tools.combine(ch_abundance)
    ch_features_with_method = ch_tools.combine(ch_features)
    ch_contrasts_with_method = ch_tools.combine(ch_contrasts)

    ch_transcript_lengths_with_method = ch_tools.combine(ch_transcript_lengths)
    ch_control_features_with_method = ch_tools.combine(ch_control_features)
    ch_gene_sets_with_method = ch_tools.combine(ch_gene_sets)

    // ----------------------------------------
    // RUN DIFFERENTIAL ABUNDANCE ANALYSIS
    // ----------------------------------------

    // run deseq2

    DIFF_DESEQ2(
        preprocess_channel(ch_abundance_with_method, "diff", "deseq2"),
        preprocess_channel(ch_transcript_lengths_with_method, "diff", "deseq2"),
        preprocess_channel(ch_control_features_with_method, "diff", "deseq2"),
        preprocess_channel(ch_samplesheet_with_method, "diff", "deseq2"),
        preprocess_channel(ch_contrasts_with_method, "diff", "deseq2"),
        "deseq2",
        params.differential_min_fold_change,
        params.differential_max_qval
    )

    ch_results_genewise = postprocess_channel(DIFF_DESEQ2.out.results_genewise, ch_tools, "diff", "deseq2")
        .mix(ch_results_genewise)
    ch_results_genewise_filtered = postprocess_channel(DIFF_DESEQ2.out.results_genewise_filtered, ch_tools, "diff", "deseq2")
        .mix(ch_results_genewise_filtered)

    // run limma

    DIFF_LIMMA(
        preprocess_channel(ch_abundance_with_method, "diff", "limma"),
        preprocess_channel(ch_transcript_lengths_with_method, "diff", "limma"),
        preprocess_channel(ch_control_features_with_method, "diff", "limma"),
        preprocess_channel(ch_samplesheet_with_method, "diff", "limma"),
        preprocess_channel(ch_contrasts_with_method, "diff", "limma"),
        "limma",
        params.differential_min_fold_change,
        params.differential_max_qval
    )

    ch_results_genewise = postprocess_channel(DIFF_LIMMA.out.results_genewise, ch_tools, "diff", "limma").mix(ch_results_genewise)
    ch_results_genewise_filtered = postprocess_channel(DIFF_LIMMA.out.results_genewise_filtered, ch_tools, "diff", "limma").mix(ch_results_genewise_filtered)

    // ----------------------------------------
    // RUN ENRICHMENT ANALYSIS
    // ----------------------------------------

    // run gsea

    // run gprofiler2

    GPROFILER2_GOST(
        preprocess_channel(ch_results_genewise_filtered, "enr", "gprofiler2"),
        preprocess_channel(ch_gene_sets_with_method, "enr", "gprofiler2"),
        preprocess_channel(ch_abundance_with_method, "enr", "gprofiler2").map { meta, abundance -> abundance }
    )

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
