//
// Perform enrichment analysis
//
include { MYGENE } from "../../../modules/nf-core/mygene/main.nf"
include { PROPR_GREA as GREA } from "../../../modules/local/propr/grea/main.nf"

workflow ENRICHMENT {
    take:
    ch_counts
    ch_results_genewise
    ch_results_genewise_filtered
    ch_adjacency
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

    emit:
    enriched = ch_enriched
}
