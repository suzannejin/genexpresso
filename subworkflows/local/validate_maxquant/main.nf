include { GUNZIP as GUNZIP_GTF                 } from '../../../modules/nf-core/gunzip/main'
include { PROTEUS_READPROTEINGROUPS as PROTEUS } from '../../../modules/nf-core/proteus/readproteingroups/main'

workflow VALIDATE_MAXQUANT {
    take:
    ch_input      // [meta, input_file]

    main:

    exp_meta = ch_input.first()[0]

    // ========================================
    // Get matrix
    // ========================================

    // Make channel for proteus

    matrix_file = file(params.matrix, checkIfExists: true)
    ch_matrix = Channel.of([ exp_meta, matrix_file])

    ch_proteus_in = ch_matrix
        .join(ch_input)
        .map { meta, matrix, samplesheet ->
            [samplesheet, matrix]}

    // We'll be running Proteus once per unique contrast variable to generate plots
    // TODO: there should probably be a separate plotting module in proteus to simplify this

    ch_contrasts_file = Channel.from([[exp_meta, file(params.contrasts)]])

    ch_contrast_variables = ch_contrasts_file
        .splitCsv ( header:true, sep:(params.contrasts.endsWith('tsv') ? '\t' : ','))
        .map{ it.tail().first() }
        .map{
            tuple('id': it.variable)
        }
        .unique()   // uniquify to keep each contrast variable only once (in case it exists in multiple lines for blocking etc.)

    // Run proteus to get protein abundance matrix

    PROTEUS(
        ch_contrast_variables.combine(ch_proteus_in)
    )
    ch_matrix = PROTEUS.out.raw_tab // or norm_tab?
        .first()
        .map{ meta, matrix -> tuple(exp_meta, matrix) }
    ch_versions = ch_versions
        .mix(PROTEUS.out.versions)

    // ========================================
    // Get features
    // ========================================

    // If user has provided a feature annotation table, use that
    if (params.features){
        ch_features = Channel.of([ exp_meta, file(params.features, checkIfExists: true)])
    } else if (params.gtf){

        // Get feature annotations from a GTF file, gunzip if necessary

        file_gtf_in = file(params.gtf)
        file_gtf = [ [ "id": file_gtf_in.simpleName ], file_gtf_in ]

        if ( params.gtf.endsWith('.gz') ){
            GUNZIP_GTF(file_gtf)
            file_gtf = GUNZIP_GTF.out.gunzip
            ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
        }

        // Get a features table from the GTF and combine with the matrix and sample
        // annotation (fom = features/ observations/ abundance)

        GTF_TO_TABLE( file_gtf, [[ "id":""], []])
        ch_features = GTF_TO_TABLE.out.feature_annotation
            .map{
                tuple( exp_meta, it[1])
            }
        ch_versions = ch_versions
            .mix(GTF_TO_TABLE.out.versions)
    }
    else{

        // Otherwise we can just use the matrix input; save it to the workdir so that it does not
        // just appear wherever the user runs the pipeline
        matrix_file = ch_matrix.first()[1]
        matrix_as_anno_filename = "${workflow.workDir}/${matrix_file.getBaseName()}_as_anno.${matrix_file.getExtension()}"
        ch_features = ch_matrix
            .map{ meta, matrix ->
                matrix_copy = file(matrix_as_anno_filename)
                matrix_copy.exists() && matrix.getText().md5().equals(matrix_copy.getText().md5()) ?: matrix.copyTo(matrix_as_anno_filename)
                [ meta, file(matrix_as_anno_filename) ]
            }
    }

    emit:
    matrix   = ch_matrix
    features = ch_features
    versions = ch_versions
}