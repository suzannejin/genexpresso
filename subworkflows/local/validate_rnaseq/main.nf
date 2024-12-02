include { GUNZIP as GUNZIP_GTF } from '../../../modules/nf-core/gunzip/main'


workflow VALIDATE_RNASEQ {

    take:
    ch_input      // [meta, input_file]

    main:

    exp_meta = ch_input.first()[0]

    // ========================================
    // Get matrix
    // ========================================

    matrix_file = file(params.matrix, checkIfExists: true)
    ch_matrix = Channel.of([ exp_meta, matrix_file])

    // ========================================
    // Get features
    // ========================================

    // If user has provided a feature annotation table, use that
    if (params.features){
        ch_features = Channel.of([ exp_meta], file(params.features, checkIfExists: true)])
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
    }else{

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