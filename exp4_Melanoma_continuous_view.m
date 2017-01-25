clear 

addpath(genpath('code'));


path.results = 'results';
path.SC = fullfile('input', 'datasets');
path.output = fullfile('output', 'continuous_view');

ds_list = dir(fullfile(path.SC, '*')); 
ds_list(~[ds_list.isdir]) = [];
ds_list(1:2) = [];
ds_no = numel(ds_list);   

dataset_names = {ds_list(:).name};


if( ~exist(path.output, 'dir') )
    system(sprintf('mkdir %s -p', path.output));
end


Melanoma_ds_id = 2;
ds_path = fullfile(path.SC, dataset_names{Melanoma_ds_id});


load(fullfile(path.results, 'ACTION', 'ACTION_full_results.mat'), 'S', 'archetypes', 'arch_annotations', 'arch_annotations_full');


%% use Fiedler/t-SNE to project results on to 2D space
    fprintf('Processing %s\n', dataset_names{Melanoma_ds_id});

    sample_annotations = my_dlmread(fullfile(ds_path, 'sample_annotations.txt'));
    Labels = sample_annotations(:, end);
    UL = unique(Labels);
    [~, l] = ismember(Labels, UL);
    
    perplexity = 50;
    [ fig, ydata, iid_points ] = plotContinuum(S{Melanoma_ds_id}, 'dim', 2, 'perplexity', perplexity, 'Labels', Labels, 'arch_labels', {});
%     [ fig, ydata, iid_points ] = plotContinuum(S{Melanoma_ds_id}, 'dim', 3, 'perplexity', perplexity, 'Labels', Labels, 'arch_labels', {}); % If 3D needed


    export_fig(fullfile(path.output, sprintf('%s_Fiedler_perp=%d.eps', dataset_names{Melanoma_ds_id}, perplexity)), '-eps', '-transparent', '-painters');        
%     close all;      
