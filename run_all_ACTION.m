clear 

addpath(genpath('code'));


path.SC = fullfile('input', 'datasets');
path.results = fullfile('results', 'ACTION');

ds_list = dir(sprintf('%s/*', path.SC)); 
ds_list(~[ds_list.isdir]) = [];
ds_list(1:2) = [];
ds_no = numel(ds_list);   

dataset_names = {ds_list(:).name};


KK = [18, 8, 29, 12];

batch_attr = { [], [], [], [], [], [] };

if( ~exist(path.results, 'dir') )
    system(sprintf('mkdir %s -p', path.results));
end

%%
ARI = zeros(numel(dataset_names), 1);
NMI = zeros(numel(dataset_names), 1);
actual_k = zeros(numel(dataset_names), 1);
Time = zeros(numel(dataset_names), 1);

adjusted_expressions = cell(numel(dataset_names), 1);
C = cell(numel(dataset_names), 1);
archetypes = cell(numel(dataset_names), 1);
S = cell(numel(dataset_names), 1);
estimated_k = zeros(numel(dataset_names), 1);
all_labels = cell(numel(dataset_names), 1);
log = cell(numel(dataset_names), 1);

arch_annotations = cell(numel(dataset_names), 1);
arch_annotations_full = cell(numel(dataset_names), 1);

for ds_id = 4%1:numel(dataset_names)
    ds_path = fullfile(path.SC, dataset_names{ds_id});
    fprintf('Reading SC dataset %s\n', dataset_names{ds_id});    

    tic; [expression, sample_names, gene_names] = my_tblread(fullfile(ds_path, 'expression.txt')); toc
    unexpressed_mask = sum(logical(expression), 2) <= 0; % a minimum of 10 samples
    expression(unexpressed_mask, :) = [];
    gene_names(unexpressed_mask) = [];

    sample_annotations = my_dlmread(fullfile(ds_path, 'sample_annotations.txt'));

    n = size(expression, 2);

    Labels = sample_annotations(:, end);
    actual_k(ds_id) = numel(unique(Labels));
    UL = unique(Labels);
    [~, l] = ismember(Labels, UL);


    fprintf('Number of true cell types = %d\n', actual_k(ds_id));
    display(UL);  

    tic
    [ adjusted_expressions{ds_id}, C{ds_id}, archetypes{ds_id}, S{ds_id}, estimated_k(ds_id), all_labels{ds_id}, log{ds_id} ] = processCells( expression, 'max_k', 50, 'experiment_name', dataset_names{ds_id}, 'batch_annotations', sample_annotations(:, batch_attr{ds_id}), 'annotations', sample_annotations(:, end), 'debug', true);
    C{ds_id} = sparse(C{ds_id});
    S{ds_id} = sparse(S{ds_id});

%     [ adjusted_expressions{ds_id}, C{ds_id}, archetypes{ds_id}, S{ds_id}, estimated_k(ds_id), all_labels{ds_id}, log{ds_id} ] = processCells( expression, 'k', KK(ds_id), 'experiment_name', dataset_names{ds_id}, 'batch_annotations', sample_annotations(:, batch_attr{ds_id}), 'annotations', sample_annotations(:, end), 'visualize', false);
    Time(ds_id) = toc;  

    
    labels = all_labels{ds_id};
    save(fullfile(path.results, sprintf('%s.mat', dataset_names{ds_id})), 'labels');    


    [ arch_annotations{ds_id} ] = archLabelEnrichment( S{ds_id}, Labels, 'pval_threshold', 1e-3 );
    [ arch_annotations_full{ds_id} ] = archEnrichment( S{ds_id}, sample_annotations, 'pval_threshold', 1e-5 );


    NMI(ds_id) = nmi(l, labels);
    ARI(ds_id) = adjustedrand(l, labels'); 
end


save(fullfile(path.results, 'ACTION_full_results.mat'), 'C', 'archetypes', 'S', 'estimated_k', 'all_labels', 'log', 'arch_annotations', 'arch_annotations_full', 'NMI', 'ARI', 'Time');
