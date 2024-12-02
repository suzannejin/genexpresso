include { UNTAR             } from '../../../modules/nf-core/untar/main.nf'
include { AFFY_JUSTRMA      } from '../../../modules/nf-core/affy/justrma/main'

workflow VALIDATE_AFFY_ARRAY {
    take:
    ch_input      // [meta, input_file]

    main:

    exp_meta = ch_input.first()[0]

    // ========================================
    // Get matrix
    // ========================================

    // channel for CEL files - mandatory

    if (params.affy_cel_files_archive) {
        ch_celfiles = Channel.of([ exp_meta, file(params.affy_cel_files_archive, checkIfExists: true) ])
    } else {
        error("CEL files archive not specified!")
    }

    // Uncompress the CEL files archive

    UNTAR ( ch_celfiles )

    ch_affy_input = ch_input
        .join(UNTAR.out.untar)

    // Run affy to derive the abundance matrix and features

    AFFY_JUSTRMA (
        ch_affy_input,
        [[],[]]
    )

    ch_matrix = AFFY_JUSTRMA.out.expression

    ch_versions = ch_versions
            .mix(AFFY_JUSTRMA.out.versions)

    // ========================================
    // Get features
    // ========================================

    if (params.features) {
        ch_features = Channel.of([ exp_meta, file(params.features, checkIfExists: true)])
    } else {
        ch_features = AFFY_JUSTRMA.out.annotation
    }

    emit:
    matrix   = ch_matrix
    features = ch_features
    versions = ch_versions
}