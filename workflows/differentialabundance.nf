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

// gene sets
ch_gene_sets = Channel.of([exp_meta, file(params.gene_sets_files, checkIfExists: true)])

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

include { ABUNDANCE_DIFFERENTIAL_FILTER as DIFFERENTIAL     } from '../subworkflows/nf-core/abundance_differential_filter/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DIFFERENTIALABUNDANCE {

    main:

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
        The nf-schema functionalities are used for the parsing.
    */

    // parse the toolsheet.csv file using nf-schema
    // and organize the channel into the proper format
    ch_tools = Channel.fromList(samplesheetToList(params.tools, './assets/schema_tools.json'))
                    .map {
                        it ->
                            def pathway_name     = it[0].subMap(["pathway_name"])
                            def differential_map = it[0].subMap(["diff_method","args_diff"])
                            def correlation_map  = it[0].subMap(["cor_method","args_cor"])
                            def enrichment_map   = it[0].subMap(["enr_method","args_enr"])
                            [ pathway_name, differential_map, correlation_map, enrichment_map ]
                    }.unique()

    // Filter the tools to the pathway(s) of interest (as indicated by --pathway flag),
    // or run everything if requested (if --pathway all)
    if (params.pathway != "all") {
        ch_tools
            .filter{
                it[0]["pathway_name"] in params.pathway.tokenize(',')
            }
            .set{ ch_tools }
    }

    // remove the pathway_name information if there is only one pathway required
    // This is done to enable resuming, for example, deseq2 analysis in the following scenario:
    //    - run 1: deseq2 + gprofiler2
    //    - run 2: deseq2 + gsea
    // Providing the same meta (with pathway_name ==  null) will enable resume work properly
    ch_tools.count()
        .combine(ch_tools)
            .map { n, pathway_name, differential_map, correlation_map, enrichment_map ->
                def new_pathway_name = n == 1 ? ['pathway_name': null] : ['pathway_name': pathway_name['pathway_name']]
                [ new_pathway_name, differential_map, correlation_map, enrichment_map ]
            }
            .set{ ch_tools }

    // ----------------------------------------
    // RUN DIFFERENTIAL ABUNDANCE ANALYSIS
    // ----------------------------------------

    // ch_input_to_differential = ch_abundance
    //     .combine()

}   

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/