//
// Perform enrichment analysis
//
include { MYGENE } from "../../../modules/nf-core/mygene/main.nf"
include { PROPR_GREA as GREA } from "../../../modules/local/propr/grea/main.nf"
include { GPROFILER2_GOST } from "../../../modules/nf-core/gprofiler2/gost/main.nf"

workflow ENRICHMENT {
    take:
    ch_counts                       // [ meta, counts] with meta keys: method, args_cor
    ch_results_genewise             // [ meta, results] with meta keys: method, args_cor
    ch_results_genewise_filtered    // [ meta, results] with meta keys: method, args_cor
    ch_adjacency                    // [ meta, adj_matrix] with meta keys: method, args_cor
    // TODO: add ch_gm when provided by user, etc.

    main:

    // initialize empty results channels
    ch_enriched = Channel.empty()
    ch_gmt      = Channel.empty()
    ch_versions = Channel.empty()

    ch_adjacency
        .branch {
            grea: it[0]["method"] == "grea"
        }
        .set { ch_adjacency }

    // ----------------------------------------------------
    // Construct gene set database
    // ----------------------------------------------------

    // TODO this should be optional, only run when there is no gene set data provided by user

    MYGENE(ch_counts.take(1))  // only one data is provided to this pipeline
    ch_gmt      = MYGENE.out.gmt
    ch_versions = ch_versions.mix(MYGENE.out.versions)

    // ----------------------------------------------------
    // Perform enrichment analysis with GREA
    // ----------------------------------------------------

    GREA(ch_adjacency.grea.unique(), ch_gmt.collect())
    ch_enriched = ch_enriched.mix(GREA.out.results)
    ch_versions = ch_versions.mix(GREA.out.versions)

    // ----------------------------------------------------
    // Perform enrichment analysis with GSEA
    // ----------------------------------------------------

    // todo: add gsea here

    // ----------------------------------------------------
    // Perform enrichment analysis with gprofiler2
    // ----------------------------------------------------

    // todo: add gprofiler2 here

    // Define background file
    if (!params.gprofiler2_background_file) {
        // If deactivated, use empty list as "background"
        ch_background = []
    } else if (params.gprofiler2_background_file == "auto") {
        // If auto, use input matrix as background
        ch_background = ch_counts.map { meta, counts -> counts }
    } else {
        ch_background = Channel.from(file(params.gprofiler2_background_file, checkIfExists: true))
    }

    // rearrange channel for GPROFILER2_GOST process
    ch_gmt = ch_gmt.map { meta, gmt -> gmt }

    ch_results_genewise_filtered
        .branch {
            grea: it[0]["method"] == "gprofiler2"
        }
        .set { ch_results_genewise_filtered }

    GPROFILER2_GOST(ch_results_genewise_filtered, ch_gmt, ch_background)
    ch_versions = ch_versions.mix(GPROFILER2_GOST.out.versions)

    emit:
    enriched = ch_enriched  // channel: [ tsv ]
    versions = ch_versions  // channel: [ versions.yml ]
}
