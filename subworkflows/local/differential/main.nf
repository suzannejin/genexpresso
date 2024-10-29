//
// Perform differential analysis
//
include { PROPR_PROPD as PROPD } from "../../../modules/local/propr/propd/main.nf"
include { DESEQ2_DIFFERENTIAL } from '../../../modules/nf-core/deseq2/differential/main'
include { LIMMA_DIFFERENTIAL } from '../../../modules/nf-core/limma/differential/main'
include { FILTER_DIFFTABLE as FILTER_DIFFTABLE_LIMMA } from '../../../modules/local/filter_difftable'

def correct_meta_data = { meta, data, pathway ->
    def meta_clone = meta.clone() + pathway
    meta_clone.remove('diff_method')
    meta_clone.remove('args_diff')
    return [meta_clone, data]
}

workflow DIFFERENTIAL {
    take:
    ch_tools              // [ pathway_name, differential_map ]
    ch_counts             // [ meta_exp, counts ]
    ch_samplesheet        // [ meta_exp, samplesheet ]
    ch_contrasts          // [ meta_contrast, contrast_variable, reference, target ]

    main:

    // initialize empty results channels
    ch_results_pairwise          = Channel.empty()
    ch_results_pairwise_filtered = Channel.empty()
    ch_results_genewise          = Channel.empty()
    ch_results_genewise_filtered = Channel.empty()
    ch_adjacency                 = Channel.empty()

    // branch tools to select the correct differential analysis method
    ch_tools
        .branch {
            propd:  it[1]["diff_method"] == "propd"
            deseq2: it[1]["diff_method"] == "deseq2"
            limma: it[1]["diff_method"] == "limma"
        }
        .set { ch_tools_single }

    // ----------------------------------------------------
    // Perform differential analysis with propd
    // ----------------------------------------------------

    // TODO propd currently don't support blocking, so we should not run propd with same contrast_variable, reference and target,
    // but different blocking variable, since it will simply run the same analysis again.

    ch_counts
        .join(ch_samplesheet)
        .combine(ch_contrasts)
        .combine(ch_tools_single.propd)
        .multiMap {
            meta_data, counts, samplesheet, meta_contrast, contrast_variable, reference, target, pathway, meta_tools ->
                def meta = meta_data.clone() + ['contrast': meta_contrast.id] + meta_tools.clone()
                input:   [ meta, counts, samplesheet, contrast_variable, reference, target ]
                pathway: [ meta, pathway ]
        }
        .set { ch_propd }

    PROPD(ch_propd.input.unique())
    ch_results_pairwise          = PROPD.out.results
                                        .join(ch_propd.pathway).map(correct_meta_data).mix(ch_results_pairwise)
    ch_results_pairwise_filtered = PROPD.out.results_filtered
                                        .join(ch_propd.pathway).map(correct_meta_data).mix(ch_results_pairwise_filtered)
    ch_results_genewise          = PROPD.out.connectivity
                                        .join(ch_propd.pathway).map(correct_meta_data).mix(ch_results_genewise)
    ch_results_genewise_filtered = PROPD.out.hub_genes
                                        .join(ch_propd.pathway).map(correct_meta_data).mix(ch_results_genewise_filtered)
    ch_adjacency                 = PROPD.out.adjacency
                                        .join(ch_propd.pathway).map(correct_meta_data).mix(ch_adjacency)

    // ----------------------------------------------------
    // Perform differential analysis with DESeq2
    // ----------------------------------------------------

    // ToDo: In order to use deseq2 the downstream processes need to be updated to process the output correctly
    // if (params.transcript_length_matrix) { ch_transcript_lengths = Channel.of([ exp_meta, file(params.transcript_length_matrix, checkIfExists: true)]).first() } else { ch_transcript_lengths = [[],[]] }
    // if (params.control_features) { ch_control_features = Channel.of([ exp_meta, file(params.control_features, checkIfExists: true)]).first() } else { ch_control_features = [[],[]] }

    // ch_samplesheet
    //     .join(ch_counts)
    //     .first()
    //     .combine(ch_tools_single.deseq2)
    //     .set { ch_counts_deseq2 }

    // DESEQ2_DIFFERENTIAL (
    //     ch_contrasts,
    //     ch_counts_deseq2,
    //     ch_control_features,
    //     ch_transcript_lengths
    // )
    // ch_results = ch_results
    //     .mix(DESEQ2_DIFFERENTIAL.out.results)

    // ----------------------------------------------------
    // Perform differential analysis with limma
    // ----------------------------------------------------

    // combine the input channels with the tools information
    // in this way, limma will only be run if the user have specified it, as informed by ch_tools
    ch_counts
        .join(ch_samplesheet)
        .combine(ch_contrasts)
        .combine(ch_tools_single.limma)
        .unique()
        .multiMap {
            meta_data, counts, samplesheet, meta_contrast, contrast_variable, reference, target, pathway, meta_tools ->
                def meta = meta_data.clone() + meta_contrast.clone() + meta_tools.clone()
                input1:  [ meta, contrast_variable, reference, target ]
                input2:  [ meta, samplesheet, counts ]
                pathway: [ meta, pathway ]
        }
        .set { ch_limma }

    // run limma
    LIMMA_DIFFERENTIAL(ch_limma.input1, ch_limma.input2)

    // filter results
    ch_logfc = Channel.value([ params.differential_fc_column, params.differential_min_fold_change ])
    ch_padj = Channel.value([ params.differential_qval_column, params.differential_max_qval ])
    FILTER_DIFFTABLE_LIMMA(
        LIMMA_DIFFERENTIAL.out.results,
        ch_logfc,
        ch_padj
    )

    // collect results
    ch_results_genewise = LIMMA_DIFFERENTIAL.out.results
                            .join(ch_limma.pathway).map(correct_meta_data).mix(ch_results_genewise)
    ch_results_genewise_filtered = FILTER_DIFFTABLE_LIMMA.out.results
                            .join(ch_limma.pathway).map(correct_meta_data).mix(ch_results_genewise_filtered)

    emit:
    results_pairwise          = ch_results_pairwise
    results_pairwise_filtered = ch_results_pairwise_filtered
    results_genewise          = ch_results_genewise
    results_genewise_filtered = ch_results_genewise_filtered
    adjacency                 = ch_adjacency
}
