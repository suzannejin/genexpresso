workflow VALIDATE_GEO_SOFT_FILE {

    take:
    ch_input      // [meta, input_file]

    main:

    exp_meta = ch_input.first()[0]

    // ========================================
    // Get matrix
    // ========================================

    // Parse the channel to query GSE

    if (params.querygse && params.features_metadata_cols) {
        ch_querygse = Channel.of([exp_meta, params.querygse])
    } else {
        error("Query GSE not specified or features metadata columns not specified")
    }

    // Query the matrix

    GEOQUERY_GETGEO(ch_querygse)

    ch_matrix = GEOQUERY_GETGEO.out.expression

    ch_versions = ch_versions
        .mix(GEOQUERY_GETGEO.out.versions)

    // ========================================
    // Get features
    // ========================================
    
    if (params.features) {
        ch_features = Channel.of([ exp_meta, file(params.features, checkIfExists: true)])
    } else {
        ch_soft_features = GEOQUERY_GETGEO.out.annotation
    }

    emit:
    matrix   = ch_matrix
    features = ch_features
    versions = ch_versions
}