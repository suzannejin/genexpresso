//
// Perform correlation analysis
//
include {PROPR_PROPR } from "../../../modules/local/propr/propr/main.nf"

workflow CORRELATION {
    take:
    ch_counts       // [ meta, counts] with meta keys: method, args_cor

    main:

    // initialize empty results channels
    ch_correlation = Channel.empty()
    ch_adjacency   = Channel.empty()

    // branch tools to select the correct correlation analysis method
    ch_counts
        .branch {
            propr:  it[0]["method"] == "propr"
        }
        .set { ch_counts }

    // ----------------------------------------------------
    // Perform correlation analysis with PROPR
    // ----------------------------------------------------

    PROPR_PROPR(ch_counts.propr.unique())
    ch_correlation = PROPR_PROPR.out.matrix.mix(ch_correlation)
    ch_adjacency = PROPR_PROPR.out.adjacency.mix(ch_adjacency)

    // TODO: divide propr module into cor, propr, pcor, pcorbshrink, etc.

    emit:
    correlation = ch_correlation
    adjacency   = ch_adjacency
}
