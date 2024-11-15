//
// Perform differential analysis
//
include { PROPR_PROPD as PROPD } from "../../../modules/local/propr/propd/main.nf"
include { DESEQ2_DIFFERENTIAL } from '../../../modules/nf-core/deseq2/differential/main'
include { DESEQ2_DIFFERENTIAL as DESEQ2_NORM } from "../../../modules/nf-core/deseq2/differential/main"
include { LIMMA_DIFFERENTIAL } from '../../../modules/nf-core/limma/differential/main'
include { FILTER_DIFFTABLE as FILTER_DIFFTABLE_LIMMA } from '../../../modules/local/filter_difftable'
include { FILTER_DIFFTABLE as FILTER_DIFFTABLE_DESEQ2 } from '../../../modules/local/filter_difftable'

workflow DIFFERENTIAL {
    take:
    ch_counts             // [ meta_exp, counts ] with meta keys: method, args_diff
    ch_samplesheet        // [ meta_exp, samplesheet ]
    ch_contrasts          // [ meta_contrast, contrast_variable, reference, target ]
    ch_transcript_lengths
    ch_control_features

    main:

    // initialize empty results channels
    ch_results_pairwise          = Channel.empty()
    ch_results_pairwise_filtered = Channel.empty()
    ch_results_genewise          = Channel.empty()
    ch_results_genewise_filtered = Channel.empty()
    ch_adjacency                 = Channel.empty()
    ch_versions                  = Channel.empty()

    // branch tools to select the correct differential analysis method
    ch_counts
        .branch {
            propd:  it[0]["method"] == "propd"
            deseq2: it[0]["method"] == "deseq2"
            limma:  it[0]["method"] == "limma"
        }
        .set { ch_counts }

    // ----------------------------------------------------
    // Perform differential analysis with propd
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
        .set { ch_propd }

    PROPD(ch_propd.unique())
    ch_results_pairwise          = PROPD.out.results.mix(ch_results_pairwise)
    ch_results_pairwise_filtered = PROPD.out.results_filtered.mix(ch_results_pairwise_filtered)
    ch_results_genewise          = PROPD.out.connectivity.mix(ch_results_genewise)
    ch_results_genewise_filtered = PROPD.out.hub_genes.mix(ch_results_genewise_filtered)
    ch_adjacency                 = PROPD.out.adjacency.mix(ch_adjacency)
    ch_versions                  = PROPD.out.versions.mix(ch_versions)

    // ----------------------------------------------------
    // Perform differential analysis with DESeq2
    // ----------------------------------------------------

    ch_counts.deseq2
        .combine(ch_samplesheet)
        .filter{ meta_counts, counts, meta_samplesheet, samplesheet -> meta_counts.subMap(meta_samplesheet.keySet()) == meta_samplesheet }
        .combine(ch_contrasts)
        .combine(ch_transcript_lengths)
        .combine(ch_control_features)
        .multiMap {
            meta_data, counts, meta_samplesheet, samplesheet, meta_contrast, contrast_variable, reference, target, meta_lengths, lengths, meta_control, control ->
                def meta = meta_data.clone() + meta_contrast.clone() + meta_lengths.clone() + meta_control.clone()
                contrast: [ meta, contrast_variable, reference, target ]
                samples_and_matrix: [ meta, samplesheet, counts ]
                control_features:   [ meta, control ]
                transcript_lengths: [ meta, lengths ]
        }
        .set { ch_deseq2 }

    // do we need this process DESEQ2_NORM?
    DESEQ2_NORM (
            ch_deseq2.contrast.first(),
            ch_deseq2.samples_and_matrix,
            ch_deseq2.control_features,
            ch_deseq2.transcript_lengths
        )

    DESEQ2_DIFFERENTIAL (
            ch_deseq2.contrast,
            ch_deseq2.samples_and_matrix,
            ch_deseq2.control_features,
            ch_deseq2.transcript_lengths
        )

    ch_norm_deseq2         = DESEQ2_NORM.out.normalised_counts
    ch_differential_deseq2 = DESEQ2_DIFFERENTIAL.out.results
    ch_model_deseq2        = DESEQ2_DIFFERENTIAL.out.model
    ch_versions            = DESEQ2_DIFFERENTIAL.out.versions.mix(ch_versions)

    ch_processed_matrices = ch_norm_deseq2
    if ('rlog' in params.deseq2_vs_method){
        ch_processed_matrices = ch_processed_matrices.join(DESEQ2_NORM.out.rlog_counts)
    }
    if ('vst' in params.deseq2_vs_method){
        ch_processed_matrices = ch_processed_matrices.join(DESEQ2_NORM.out.vst_counts)
    }
    ch_processed_matrices = ch_processed_matrices
        .map{ it.tail() }

    // TODO modify the module to accept these parameters as meta/ext.args in the same way how propd does
    ch_logfc_deseq2 = Channel.value([ "log2FoldChange", params.differential_min_fold_change ])
    ch_padj_deseq2 = Channel.value([ "padj", params.differential_max_qval ])

    FILTER_DIFFTABLE_DESEQ2(
        ch_differential_deseq2,
        ch_logfc_deseq2,
        ch_padj_deseq2
    )

    ch_results_genewise          = DESEQ2_DIFFERENTIAL.out.results.mix(ch_results_genewise)
    ch_results_genewise_filtered = FILTER_DIFFTABLE_DESEQ2.out.filtered.mix(ch_results_genewise_filtered)
    ch_versions                  = FILTER_DIFFTABLE_DESEQ2.out.versions.mix(ch_versions)

    // ----------------------------------------------------
    // Perform differential analysis with limma
    // ----------------------------------------------------

    // combine the input channels with the tools information
    // in this way, limma will only be run if the user have specified it, as informed by ch_tools
    ch_counts.limma
        .combine(ch_samplesheet)
        .filter{ meta_counts, counts, meta_samplesheet, samplesheet -> meta_counts.subMap(meta_samplesheet.keySet()) == meta_samplesheet }
        .combine(ch_contrasts)
        .unique()
        .multiMap {
            meta_data, counts, meta_samplesheet, samplesheet, meta_contrast, contrast_variable, reference, target ->
                def meta = meta_data.clone() + meta_contrast.clone()
                input1:  [ meta, contrast_variable, reference, target ]
                input2:  [ meta, samplesheet, counts ]
        }
        .set { ch_limma }


    // run limma
    LIMMA_DIFFERENTIAL(ch_limma.input1, ch_limma.input2)
    ch_versions = LIMMA_DIFFERENTIAL.out.versions.mix(ch_versions)

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

    // collect results
    ch_results_genewise          = LIMMA_DIFFERENTIAL.out.results.mix(ch_results_genewise)
    ch_results_genewise_filtered = FILTER_DIFFTABLE_LIMMA.out.filtered.mix(ch_results_genewise_filtered)
    ch_versions                  = FILTER_DIFFTABLE_LIMMA.out.versions.mix(ch_versions)

    emit:
    results_pairwise          = ch_results_pairwise           // channel: [ tsv ]
    results_pairwise_filtered = ch_results_pairwise_filtered  // channel: [ tsv ]
    results_genewise          = ch_results_genewise           // channel: [ tsv ]
    results_genewise_filtered = ch_results_genewise_filtered  // channel: [ tsv ]
    adjacency                 = ch_adjacency                  // channel: [ tsv ]
    versions                  = ch_versions                   // channel: [ versions.yml ]
}
