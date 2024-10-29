//
// Run experimental analysis
//
include { DIFFERENTIAL } from '../differential/main.nf'
include { CORRELATION } from '../correlation/main.nf'
include { ENRICHMENT } from '../enrichment/main.nf'

include { CUSTOM_MATRIXFILTER } from '../../../modules/nf-core/custom/matrixfilter/main'
include { SHINYNGS_STATICEXPLORATORY as PLOT_EXPLORATORY } from '../../../modules/nf-core/shinyngs/staticexploratory/main.nf'

workflow EXPERIMENTAL {
    take:
    ch_contrasts    // [ meta, contrast_variable, reference, target ]
    ch_samplesheet  // [ meta, samplesheet ]
    ch_features     // [ meta, features ]
    ch_counts       // [ meta, counts]
    ch_tools        // [ pathway_name, differential_map, correlation_map, enrichment_map ]

    main:

    // split toolsheet into channels
    ch_tools
        .multiMap{
            pathway_name, differential_map, correlation_map, enrichment_map ->
                diff: [ pathway_name, differential_map ]
                corr: [ pathway_name, correlation_map ]
                enr:  [ pathway_name, enrichment_map ]
        }
        .set{ ch_tools }


    // initialize empty results channels
    ch_results_pairwise = Channel.empty()               // differential results for pairwise analysis - it should be a table
    ch_results_pairwise_filtered = Channel.empty()      // differential results for pairwise analysis - filtered - it should be a table
    ch_results_genewise = Channel.empty()               // differential results for genewise analysis - it should be a table
    ch_results_genewise_filtered = Channel.empty()      // differential results for genewise analysis - filtered - it should be a table
    ch_adjacency = Channel.empty()                      // adjacency matrix showing the connections between the genes, with values 1|0
    ch_matrix = Channel.empty()                         // correlation matrix
    ch_enriched = Channel.empty()                       // output table from enrichment analysis

    // ----------------------------------------------------
    // DATA PREPROCESSING
    // ----------------------------------------------------

    // filter out low-abundance features
    CUSTOM_MATRIXFILTER(
        ch_counts,
        ch_samplesheet
    )
    ch_counts_filtered = CUSTOM_MATRIXFILTER.out.filtered

    // ----------------------------------------------------
    // DIFFERENTIAL ANALYSIS BLOCK
    // ----------------------------------------------------

    DIFFERENTIAL(
        ch_tools.diff,
        ch_counts_filtered,
        ch_samplesheet,
        ch_contrasts
    )
    ch_results_pairwise = ch_results_pairwise.mix(DIFFERENTIAL.out.results_pairwise)
    ch_results_pairwise_filtered = ch_results_pairwise_filtered.mix(DIFFERENTIAL.out.results_pairwise_filtered)
    ch_results_genewise = ch_results_genewise.mix(DIFFERENTIAL.out.results_genewise)
    ch_results_genewise_filtered = ch_results_genewise_filtered.mix(DIFFERENTIAL.out.results_genewise_filtered)
    ch_adjacency = ch_adjacency.mix(DIFFERENTIAL.out.adjacency)

    // ----------------------------------------------------
    // CORRELATION ANALYSIS BLOCK
    // ----------------------------------------------------

    CORRELATION(
        ch_tools.corr,
        ch_counts_filtered
    )
    ch_matrix = ch_matrix.mix(CORRELATION.out.matrix)
    ch_adjacency = ch_adjacency.mix(CORRELATION.out.adjacency)

    // ----------------------------------------------------
    // FUNCTIONAL ENRICHMENT BLOCK
    // ----------------------------------------------------

    ENRICHMENT(
        ch_tools.enr,
        ch_counts_filtered,
        ch_results_genewise,
        ch_results_genewise_filtered,
        ch_adjacency
    )
    ch_enriched = ch_enriched.mix(ENRICHMENT.out.enriched)

    // ----------------------------------------------------
    // VISUALIZATION BLOCK
    // ----------------------------------------------------

    // compute exploratory plots for data
    // TODO currently it is visualizing the raw counts, but it should also visualize the differently normalized data

    ch_contrasts
        .map {
            [ "id": it[1] ]
        }
        .unique()
        .set{ ch_contrast_variables }

    ch_counts_filtered
        .join(ch_samplesheet)
        .join(ch_features)
        .combine(ch_contrast_variables)
        .map { meta_counts, counts, samplesheet, features, meta_contrast ->
            def meta = meta_counts.clone() + meta_contrast.clone()
            [ meta, samplesheet, features, counts ]
        }
        .set { ch_to_plot_exploratory}
    ch_to_plot_exploratory.view()
    PLOT_EXPLORATORY(ch_to_plot_exploratory)

    // TODO: add other visualization stuff here
}
