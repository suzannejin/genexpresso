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

    // ----------------------------------------------------
    // Construct gene set database
    // ----------------------------------------------------

    // TODO this should be optional, only run when there is no gene set data provided by user

    MYGENE(ch_counts.take(1))  // only one data is provided to this pipeline
    ch_gmt = MYGENE.out.gmt

    // ----------------------------------------------------
    // Perform enrichment analysis with GREA
    // ----------------------------------------------------

    ch_adjacency
        .filter { it[0]["method"] == "grea" }
        .unique()
        .set { ch_adjacency_to_grea }

    GREA(ch_adjacency_to_grea, ch_gmt.collect())
    ch_enriched = ch_enriched.mix(GREA.out.results)

    // ----------------------------------------------------
    // Perform enrichment analysis with GSEA
    // ----------------------------------------------------

    // todo: add gsea here

    // ----------------------------------------------------
    // Perform enrichment analysis with gprofiler2
    // ----------------------------------------------------

    // parse input channels
    // TODO we need to find a way to combine these information with also args coming from toolsheet and modules.config

    if (!params.gprofiler2_background_file) {  // If deactivated, use empty list as "background"
        ch_background = []
    } else if (params.gprofiler2_background_file == "auto") {  // If auto, use input matrix as background
        ch_background = ch_counts.map { meta, counts -> counts }
    } else {
        ch_background = Channel.from(file(params.gprofiler2_background_file, checkIfExists: true))
    }

    ch_gmt = ch_gmt.map { meta, gmt -> gmt }

    ch_results_genewise_filtered
        .filter { it[0]["method"] == "gprofiler2" }
        .unique()
        .set { ch_for_gprofiler2 }

    // run gprofiler2

    GPROFILER2_GOST(
        ch_for_gprofiler2,
        ch_gmt,
        ch_background
    )

    // collect results

    ch_enriched = ch_enriched.mix(GPROFILER2_GOST.out.all_enrich)

    // TODO also collect the files for report purposes
    emit:
    enriched = ch_enriched
}
