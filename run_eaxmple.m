clear 

addpath(genpath('code'));

path.SC = fullfile('input', 'datasets');
dataset_name = 'Pollen';
ds_path = fullfile(path.SC, dataset_name);

path.output = fullfile('output', 'example');
if( ~exist(path.output, 'dir') )
    system(sprintf('mkdir %s -p', path.output));
end

%% Read Pollen dataset
    tic; [expression, sample_names, gene_names] = my_tblread(fullfile(ds_path, 'expression.txt')); toc

    mask = (sum(expression, 2) == 0);
    expression(mask, :) = [];
    gene_names(mask) = [];

    sample_annotations = my_dlmread(fullfile(ds_path, 'sample_annotations.txt'));
    Labels = sample_annotations(:, end);

    n = size(expression, 2);

%% Identify archetypes -- functional identities  
    [ adjusted_expressions, C, archetypes, S, cell_type_no, cell_types, log ] = processCells( expression, 'max_k', 50, 'experiment_name', dataset_name, 'annotations', sample_annotations(:, end), 'ER_pval_threshold', 0.05, 'debug', true );
    
    save(fullfile(path.output, dataset_name, '.mat'), 'adjusted_expressions', 'C', 'archetypes', 'S', 'G', 'cell_types', 'node_labels', 'log');

%% Assess enrichment of cell type labels around each archetype.
    [ arch_annotations ] = archLabelEnrichment( S, Labels, 'pval_threshold', 1e-3 );
    [ arch_annotations_full ] = archEnrichment( S, sample_annotations, 'pval_threshold', 1e-5 );
    close(child_handles(mask));
    
    
%% Project functional identities onto 2D plane
    UL = unique(Labels);
    [~, l] = ismember(Labels, UL);
    
    perplexity = 50;
    [ fig, ydata, iid_points ] = plotContinuum(S, 'dim', 2, 'perplexity', perplexity, 'Labels', Labels, 'arch_labels', {});


    export_fig(fullfile(path.output, sprintf('%s_Fiedler_perp=%d.eps', dataset_name, perplexity)), '-eps', '-transparent', '-painters'); 

%% Construct TRN
    RegNet = my_dlmread(fullfile('input', 'TRN', 'trrust_rawdata.txt'), '\t');
    negative_edges = strncmp(RegNet(:, 3), 'Repr', 4);
    RegNet(negative_edges, :)= [];
    TRN_edgeList = RegNet(:, [1, 2]);


    selfloop_mask = arrayfun(@(row) strcmp(TRN_edgeList(row, 1), TRN_edgeList(row, 2)), 1:size(TRN_edgeList, 1));
    TRN_edgeList(selfloop_mask, :) = [];

    [residual_expression, sorted_genes, TFs, mHG_pvals, mHG_pvals_corrected, selected_genes, TF_table, subTRN_EdgeList, subTRN_NodeAnnotations] = dissectTRN( archetypes, gene_names, TRN_edgeList, 1:size(archetypes, 2), 'pval_threshold', 0.05);

    save(fullfile(path.output, dataset_names, '_TRN.mat'), 'residual_expression', 'sorted_genes', 'TFs', 'mHG_pvals', 'mHG_pvals_corrected', 'selected_genes', 'TF_table', 'subTRN_EdgeList', 'subTRN_NodeAnnotations');
