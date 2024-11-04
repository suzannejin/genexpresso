//
// Perform differential analysis
//
include { PROPR_PROPD } from "../../../modules/local/propr/propd/main.nf"
include { DESEQ2_DIFFERENTIAL as DESEQ2_NORM } from "../../../modules/nf-core/deseq2/differential/main"
include { DESEQ2_DIFFERENTIAL } from '../../../modules/nf-core/deseq2/differential/main'
include { LIMMA_DIFFERENTIAL } from '../../../modules/nf-core/limma/differential/main'
include { FILTER_DIFFTABLE as FILTER_DIFFTABLE_DESEQ2 } from '../../../modules/local/filter_difftable'
include { FILTER_DIFFTABLE as FILTER_DIFFTABLE_LIMMA } from '../../../modules/local/filter_difftable'

workflow DIFFERENTIAL {
    take:
    ch_counts             // [ meta_exp, counts ] with meta keys: method, args_diff
    ch_samplesheet        // [ meta_exp, samplesheet ]
    ch_contrasts          // [ meta_contrast, contrast_variable, reference, target ]

    main:

    // initialize empty results channels
    // NOTE that ch_results pairwise and adjacency are a special case of results, which stores pairwise DE results (see propd)
    // whereas ch_results stores gene-wise DE results (traditional methods like deseq2 and limma only provide gene-wise results)

    ch_results_genewise          = Channel.empty()
    ch_results_genewise_filtered = Channel.empty()
    ch_results_pairwise          = Channel.empty()
    ch_results_pairwise_filtered = Channel.empty()
    ch_adjacency                 = Channel.empty()
    ch_norm                      = Channel.empty()  // channel to store the normalized data
    ch_norm_for_plotting         = Channel.empty()  // channel to store the normalized data for plotting

    // branch the data channel to the correct differential analysis method, based on the method key specified in the meta data
    ch_counts
        .branch {
            propd:  it[0]["method"] == "propd"
            deseq2: it[0]["method"] == "deseq2"
            limma:  it[0]["method"] == "limma"
        }
        .set { ch_counts }

    // ----------------------------------------------------
    // Perform differential analysis with PROPD
    // ----------------------------------------------------

    // TODO propd currently don't support blocking, so we should not run propd with same contrast_variable, reference and target,
    // but different blocking variable, since it will simply run the same analysis again.

    ch_counts.propd
        .combine(ch_samplesheet)
        .filter{ meta_counts, counts, meta_samplesheet, samplesheet -> meta_counts.subMap(meta_samplesheet.keySet()) == meta_samplesheet }
        .combine(ch_contrasts)
        .map {
            meta_data, counts, meta_samplesheet, samplesheet, meta_contrast, contrast_variable, reference, target ->
                def meta = meta_data.clone() + ['contrast': meta_contrast.id]
                return [ meta, counts, samplesheet, contrast_variable, reference, target ]
        }
        .unique()
        .set { ch_propd }

    PROPR_PROPD(ch_propd)

    ch_results_genewise          = PROPR_PROPD.out.connectivity.mix(ch_results_genewise)
    ch_results_genewise_filtered = PROPR_PROPD.out.hub_genes.mix(ch_results_genewise_filtered)
    ch_results_pairwise          = PROPR_PROPD.out.results.mix(ch_results_pairwise)
    ch_results_pairwise_filtered = PROPR_PROPD.out.results_filtered.mix(ch_results_pairwise_filtered)
    ch_adjacency                 = PROPR_PROPD.out.adjacency.mix(ch_adjacency)

    // ----------------------------------------------------
    // Perform differential analysis with DESEQ2
    // ----------------------------------------------------

    // parse input channels for deseq2 modules

    if (params.transcript_length_matrix) { ch_transcript_lengths = Channel.of([ exp_meta, file(params.transcript_length_matrix, checkIfExists: true)]).first() } else { ch_transcript_lengths = Channel.of([[],[]]) }
    if (params.control_features) { ch_control_features = Channel.of([ exp_meta, file(params.control_features, checkIfExists: true)]).first() } else { ch_control_features = Channel.of([[],[]]) }
    ch_counts.deseq2
        .combine(ch_samplesheet)
        .combine(ch_contrasts)
        .combine(ch_transcript_lengths)
        .combine(ch_control_features)
        .unique()
        .multiMap {
            meta_counts, counts, meta_samplesheet, samplesheet, meta_contrast, contrast_variable, reference, target, meta_lengths, lengths, meta_control, control ->
                def meta = meta_counts.clone() + ['contrast': meta_contrast.id]
                contrast: [ meta, contrast_variable, reference, target ]
                samples_and_matrix: [ meta, samplesheet, counts ]
                control_features:   [ meta, control ]
                transcript_lengths: [ meta, lengths ]
        }
        .set { ch_deseq2 }

    // normalize the data using deseq2
    // NOTE that for the moment the output of this module is only needed for plot_exploratory,
    // and as input for GSEA (for the moment GSEA cannot take the preranked DE results from DESEQ2_DIFFERENTIAL)

    DESEQ2_NORM (
            ch_deseq2.contrast.first(),
            ch_deseq2.samples_and_matrix,
            ch_deseq2.control_features,
            ch_deseq2.transcript_lengths
        )

    ch_norm = DESEQ2_NORM.out.normalised_counts.mix(ch_norm)
    ch_norm_for_plotting = DESEQ2_NORM.out.normalised_counts
                            .join(DESEQ2_NORM.out.rlog_counts)
                            .join(DESEQ2_NORM.out.vst_counts) // CHECK if this is correct, otherwise add only if these outputs are present
                            .map{ it.tail() }
                            .mix(ch_norm_for_plotting)

    // run DE analysis using deseq2

    DESEQ2_DIFFERENTIAL (
            ch_deseq2.contrast,
            ch_deseq2.samples_and_matrix,
            ch_deseq2.control_features,
            ch_deseq2.transcript_lengths
        )

    ch_results_genewise = DESEQ2_DIFFERENTIAL.out.results.mix(ch_results_genewise)

    // filter results

    // TODO modify the module to accept these parameters as meta/ext.args in the same way how propd does
    ch_logfc_deseq2 = Channel.value([ "log2FoldChange", params.differential_min_fold_change ])
    ch_padj_deseq2 = Channel.value([ "padj", params.differential_max_qval ])

    FILTER_DIFFTABLE_DESEQ2(
        DESEQ2_DIFFERENTIAL.out.results,
        ch_logfc_deseq2,
        ch_padj_deseq2
    )

    ch_results_genewise_filtered = FILTER_DIFFTABLE_DESEQ2.out.filtered.mix(ch_results_genewise_filtered)

    // ----------------------------------------------------
    // Perform differential analysis with limma
    // ----------------------------------------------------

    // parse input channels for limma
    // TODO provide the normalized data to limma

    ch_counts.limma
        .combine(ch_samplesheet)
        .filter{ meta_counts, counts, meta_samplesheet, samplesheet -> meta_counts.subMap(meta_samplesheet.keySet()) == meta_samplesheet }
        .combine(ch_contrasts)
        .unique()
        .multiMap {
            meta_counts, counts, meta_samplesheet, samplesheet, meta_contrast, contrast_variable, reference, target ->
                def meta = meta_counts.clone() + meta_contrast.clone()
                input1:  [ meta, contrast_variable, reference, target ]
                input2:  [ meta, samplesheet, counts ]
        }
        .set { ch_limma }

    // run limma

    LIMMA_DIFFERENTIAL(ch_limma.input1, ch_limma.input2)

    ch_results_genewise = LIMMA_DIFFERENTIAL.out.results.mix(ch_results_genewise)

    // filter results

    // note that these are column names specific for limma output table
    // TODO modify the module to accept these parameters as meta/ext.args in the same way how propd does
    ch_logfc_limma = Channel.value([ "logFC", params.differential_min_fold_change ])
    ch_padj_limma = Channel.value([ "adj.P.Val", params.differential_max_qval ])

    FILTER_DIFFTABLE_LIMMA(
        LIMMA_DIFFERENTIAL.out.results,
        ch_logfc_limma,
        ch_padj_limma
    )

    ch_results_genewise_filtered = FILTER_DIFFTABLE_LIMMA.out.filtered.mix(ch_results_genewise_filtered)

    emit:
    results_genewise          = ch_results_genewise
    results_genewise_filtered = ch_results_genewise_filtered
    results_pairwise          = ch_results_pairwise
    results_pairwise_filtered = ch_results_pairwise_filtered
    adjacency                 = ch_adjacency
    ch_norm                   = ch_norm
    ch_norm_for_plotting      = ch_norm_for_plotting
}
