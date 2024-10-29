//
// Perform enrichment analysis
//
include { MYGENE } from "../../../modules/nf-core/mygene/main.nf"
include { PROPR_GREA as GREA } from "../../../modules/local/propr/grea/main.nf"
include { GPROFILER2_GOST } from "../../../modules/nf-core/gprofiler2/gost/main.nf"

workflow ENRICHMENT {
    take:
    ch_tools        // [ pathway_name, enrichment_map ]
    ch_counts
    ch_results_genewise
    ch_results_genewise_filtered
    ch_adjacency
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
        .map { meta, matrix -> [meta.subMap(["pathway_name"]), meta, matrix] }
        .join(ch_tools, by: [0])
        .map {
            pathway_name, meta, matrix, meta_tools ->
                def new_meta = meta.clone() + meta_tools.clone()
                [ new_meta, matrix ]
            }
        .branch {
            grea:  it[0]["enr_method"] == "grea"
            gsea: it[0]["enr_method"] == "gsea"
        }
        .set { ch_adjacency }

    GREA(ch_adjacency.grea, ch_gmt.collect()) //currently, ch_gmt.collect() returns an empty channel, so this does not run
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
        .map { meta, matrix -> [meta.subMap(["pathway_name"]), meta, matrix] }
        .join(ch_tools, by: [0])
        .map {
            pathway_name, meta, matrix, meta_tools ->
                def new_meta = meta.clone() + meta_tools.clone()
                [ new_meta, matrix ]
        }
        .filter {
            it[0].enr_method == "gprofiler2"
        }
        .set { ch_results_genewise_filtered }

    GPROFILER2_GOST(ch_results_genewise_filtered, ch_gmt, ch_background)

    emit:
    enriched = ch_enriched
}