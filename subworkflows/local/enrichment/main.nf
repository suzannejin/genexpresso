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

    ch_adjacency
        .branch {
            grea: it[0]["method"] == "grea"
            gsea: it[0]["method"] == "gsea"
        }
        .set { ch_adjacency }

    // ----------------------------------------------------
    // Construct gene set database
    // ----------------------------------------------------

    // TODO this should be optional, only run when there is no gene set data provided by user

    // empty counts channel of ch_adjacency is empty to skip unnecessary MYGENE computations
    ch_counts
        .combine(ch_adjacency.grea)
        .map{ meta_counts, counts, meta_adjacency, adjacency -> [meta_counts, counts]}
        .unique()
        .set{ch_counts}

    MYGENE(ch_counts.take(1))  // only one data is provided to this pipeline
    ch_gmt = MYGENE.out.gmt

    // ----------------------------------------------------
    // Perform enrichment analysis with GREA
    // ----------------------------------------------------

    GREA(ch_adjacency.grea.unique(), ch_gmt.collect())
    ch_enriched = ch_enriched.mix(GREA.out.results)

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

    // rearrage channel for GPROFILER2_GOST process
    ch_gmt = ch_gmt.map { meta, gmt -> gmt }

    ch_results_genewise_filtered
        .branch {
            grea: it[0]["method"] == "gprofiler2"
        }
        .set { ch_results_genewise_filtered }

    GPROFILER2_GOST(ch_results_genewise_filtered, ch_gmt, ch_background)

    emit:
    enriched = ch_enriched
}
