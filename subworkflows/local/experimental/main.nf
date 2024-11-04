/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { DIFFERENTIAL }        from '../differential/main.nf'
include { CORRELATION }         from '../correlation/main.nf'
include { ENRICHMENT }          from '../enrichment/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DEFINE AUXILIARY FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// function used to preprocess the input of a subworkflow
// Basically, given ch_input and ch_tools_args,
// it updates the meta data of ch_input with the required method value and args

def preprocess_subworkflow_input( ch_input, ch_tools_args, method_field_name) {
    return ch_input
        .combine(ch_tools_args)
        // filter the tools_args to match the pathway_name
        // if no pathway_name is provided, then it will match all
        // NOTE that no pathway_name is provided when ch_tools only have one element. By doing this, it avoids the recomputation of processes that use the same tool/args but belong to different pathways
        .filter{  meta, input, pathway, arg_maps -> meta["pathway_name"] ? meta["pathway_name"] == pathway["pathway_name"] : true }
        // update meta with method value and args
        .map{ meta, input, pathway, arg_map ->
            def meta_clone = meta.clone() + pathway + arg_map.clone()
            meta_clone["method"] = meta_clone.remove(method_field_name)
            return [meta_clone, input]
        }
}

// function used to postprocess the output of a subworkflow
// Basically, given ch_results and field_name,
// it removes the field_name from the meta data

def postprocess_subworkflow_output( ch_results, field_name ) {
    return ch_results
        .map{ meta, data ->
            def meta_clone = meta.clone()
            meta_clone.removeAll{it.key in field_name}
            return [meta_clone, data]
        }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow EXPERIMENTAL {
    take:
    ch_contrasts    // [ meta, contrast_variable, reference, target ]
    ch_samplesheet  // [ meta, samplesheet ]
    ch_counts       // [ meta, counts]
    ch_tools        // [ pathway_name, differential_map, correlation_map, enrichment_map ]

    main:

    // split toolsheet into channels: diff, corr, enr
    // NOTE that when only one tools pathway is provided, the pathway_name is not needed, and hence removed
    // by doing so, it avoids the recomputation of processes that use the same tool/args but belong to different pathways

    ch_tools.count()
        .combine(ch_tools)
        .multiMap{
            n, pathway, differential_map, correlation_map, enrichment_map ->
                def pathway_name = n == 1 ? ["pathway_name":""] : pathway
                diff: [ pathway_name, differential_map ]
                corr: [ pathway_name, correlation_map ]
                enr:  [ pathway_name, enrichment_map ]
        }
        .set{ ch_tools }

    // initialize empty results channels

    ch_results_genewise = Channel.empty()               // differential results for genewise analysis - it should be a table
    ch_results_genewise_filtered = Channel.empty()      // differential results for genewise analysis - filtered - it should be a table
    ch_results_pairwise = Channel.empty()               // differential results for pairwise analysis - it should be a table
    ch_results_pairwise_filtered = Channel.empty()      // differential results for pairwise analysis - filtered - it should be a table
    ch_adjacency = Channel.empty()                      // adjacency matrix showing the connections between the genes, with values 1|0
    ch_correlation = Channel.empty()                    // correlation matrix
    ch_enriched = Channel.empty()                       // output table from enrichment analysis

    // ----------------------------------------------------
    // DIFFERENTIAL ANALYSIS BLOCK
    // ----------------------------------------------------

    // preprocess the input of the subworkflow, by adding 'method' and 'args_diff'
    // NOTE that here we only preprocess ch_counts, to use it as the carrier of the method information
    // since the rest of data channels would be combined with the ch_counts correspondingly, we don't need to preprocess them all
    // This MUST be changed, if this is not the case (eg. if one wants to use a specific ch_contrasts element for a specific method, etc)

    preprocess_subworkflow_input(ch_counts, ch_tools.diff, "diff_method")
        .set{ ch_counts_diff }

    // run differential subworkflow

    DIFFERENTIAL(
        ch_counts_diff,
        ch_samplesheet,
        ch_contrasts
    )

    // collect and postprocess the output of the subworkflow, by removing 'method' and 'args_diff'

    ch_results_genewise          = postprocess_subworkflow_output(DIFFERENTIAL.out.results_genewise,["method", "args_diff"]).mix(ch_results_genewise)
    ch_results_genewise_filtered = postprocess_subworkflow_output(DIFFERENTIAL.out.results_genewise_filtered,["method", "args_diff"]).mix(ch_results_genewise_filtered)
    ch_results_pairwise          = postprocess_subworkflow_output(DIFFERENTIAL.out.results_pairwise,["method", "args_diff"]).mix(ch_results_pairwise)
    ch_results_pairwise_filtered = postprocess_subworkflow_output(DIFFERENTIAL.out.results_pairwise_filtered,["method", "args_diff"]).mix(ch_results_pairwise_filtered)
    ch_adjacency                 = postprocess_subworkflow_output(DIFFERENTIAL.out.adjacency,["method", "args_diff"]).mix(ch_adjacency)

    // ----------------------------------------------------
    // CORRELATION ANALYSIS BLOCK
    // ----------------------------------------------------

    // preprocess the input of the subworkflow, by adding 'method' and 'args_cor'

    preprocess_subworkflow_input(ch_counts, ch_tools.corr, "cor_method")
        .set{ ch_counts_corr }

    // run correlation subworkflow

    CORRELATION(
        ch_counts_corr
    )

    // collect and postprocess the output of the subworkflow, by removing 'method' and 'args_cor'

    ch_correlation = postprocess_subworkflow_output(CORRELATION.out.matrix,["method", "args_cor"]).mix(ch_correlation)
    ch_adjacency   = postprocess_subworkflow_output(CORRELATION.out.adjacency,["method", "args_cor"]).mix(ch_adjacency)

    // ----------------------------------------------------
    // FUNCTIONAL ENRICHMENT BLOCK
    // ----------------------------------------------------

    // preprocess the input of the subworkflow, by adding 'method' and 'args_enr'
    // this is done by matching the 'pathway_name' in ch_tools

    preprocess_subworkflow_input(ch_counts, ch_tools.enr, "enr_method")
        .set{ ch_counts_enr }
    preprocess_subworkflow_input(ch_results_genewise, ch_tools.enr, "enr_method")
        .set{ ch_results_genewise_enr }
    preprocess_subworkflow_input(ch_results_genewise_filtered, ch_tools.enr, "enr_method")
        .set{ ch_results_genewise_filtered_enr }
    preprocess_subworkflow_input(ch_adjacency, ch_tools.enr, "enr_method")
        .set{ ch_adjacency_enr }

    // run enrichment subworkflow
    // TODO don't run if no enrichment analysis is needed

    ENRICHMENT(
        ch_counts_enr,
        ch_results_genewise_enr,
        ch_results_genewise_filtered_enr,
        ch_adjacency_enr
    )

    // collect the output of the subworkflow

    ch_enriched = ch_enriched.mix(ENRICHMENT.out.enriched)

    // ----------------------------------------------------
    // VISUALIZATION BLOCK
    // ----------------------------------------------------

    // TODO: call visualization stuff here
}
