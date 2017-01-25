clear 

addpath(genpath('code'));


path.results = 'results';
path.SC = fullfile('input', 'datasets');
path.output = fullfile('output', 'kernel');

ds_list = dir(sprintf('%s/*', path.SC)); 
ds_list(~[ds_list.isdir]) = [];
ds_list(1:2) = [];
ds_no = numel(ds_list);   

dataset_names = {ds_list(:).name};
methods = {'ACTION', 'Isomap', 'MDS', 'SIMLR'};


if( ~exist(path.output, 'dir') )
    system(sprintf('mkdir %s -p', path.output));
    mkdir(fullfile(path.output, 'Mat'));    
    mkdir(fullfile(path.output, 'Heatmap'));    
end

path.results = fullfile('results', 'kernel');

sample_no = 100; % for Kernel-Kmeans
seed = 0;

%%

    NMI = zeros(numel(dataset_names), numel(methods));
    ARI = zeros(numel(dataset_names), numel(methods));
    Timing = zeros(numel(dataset_names), numel(methods));
    dataset_size = zeros(numel(dataset_names), 1);
    all_labels = cell(sample_no, 1);
    energy = ones(sample_no, 1);
    
    ARI_samples_cell = cell(numel(methods), 1);
    NMI_samples_cell = cell(numel(methods), 1);
    for ds_id = 1:numel(dataset_names)
        
        ds_path = fullfile(path.SC, dataset_names{ds_id});
        fprintf('Reading SC dataset %s\n', dataset_names{ds_id});


        tic; [expression, sample_names, gene_names] = my_tblread(fullfile(ds_path, 'expression.txt')); toc
        unexpressed_mask = sum(expression, 2) == 0;
        expression(unexpressed_mask, :) = [];
        gene_names(unexpressed_mask) = [];

        sample_annotations = my_dlmread(fullfile(ds_path, 'sample_annotations.txt'));
        Labels = sample_annotations(:, end);
        k = numel(unique(Labels));
        UL = unique(Labels);
        [~, true_labels] = ismember(Labels, UL);

        dataset_size(ds_id) = size(sample_annotations, 1);

        for method_index = 1:numel(methods)        
            method_name = methods{method_index};
            fprintf('Running %s ...\n', method_name);
            exp_name = strcat(dataset_names{ds_id}, '_', method_name);

            % Compute kernel
            rng(seed)
            if(strcmp(method_name, 'ACTION'))
                S = getSimilarity(expression, 'ACTION');
                save(fullfile(path.results, sprintf('%s_%s.mat', dataset_names{ds_id}, method_name)), 'S');
%                 load(fullfile(path.results, sprintf('%s_%s.mat', dataset_names{ds_id}, method_name)));
                
                NMI_samples = zeros(sample_no, 1);
                ARI_samples = zeros(sample_no, 1);
                for j = 1:sample_no
                    fprintf('\tsample %d/%d\n', j, sample_no);
                    labels = k2means(S, k);
                    NMI_samples(j) = nmi(labels, true_labels);
                    ARI_samples(j) = adjustedrand(labels, true_labels);
                end
                NMI_samples_cell{method_index}{ds_id, 1} = NMI_samples;
                ARI_samples_cell{method_index}{ds_id, 1} = ARI_samples;

                NMI(ds_id, method_index) = mean(NMI_samples);
                ARI(ds_id, method_index) = mean(ARI_samples);   
                
                NMI_idx(ds_id, method_index) = 1;
                ARI_idx(ds_id, method_index) = 1;                
            else                
                for K = 5:5:50
                    fprintf('K = %d\n', K);
                    S = getSimilarity(expression, method_name, 'dim', K, 'k', numel(UL));
                    save(sprintf('%s_%s_k=%d.mat', dataset_names{ds_id}, method_name, K), 'S');
%                     load(fullfile(path.results, sprintf('%s_%s_k=%d.mat', dataset_names{ds_id}, method_name, K)));


                    NMI_samples = zeros(sample_no, 1);
                    ARI_samples = zeros(sample_no, 1);
                    for j = 1:sample_no
                        fprintf('\tsample %d/%d\n', j, sample_no);
                        labels = k2means(S, k);
                        NMI_samples(j) = nmi(labels, true_labels);
                        ARI_samples(j) = adjustedrand(labels, true_labels);
                    end
                    NMI_samples_cell{method_index}{ds_id, K} = NMI_samples;
                    ARI_samples_cell{method_index}{ds_id, K} = ARI_samples;
                                        
                    NMI_K{method_index}(ds_id, K) = mean(NMI_samples);
                    ARI_K{method_index}(ds_id, K) = mean(ARI_samples);           
                    
                    fprintf('k = %d, NMI = %.2f, ARI = %.2f\n', K, NMI_K{method_index}(ds_id, K), ARI_K{method_index}(ds_id, K));
                end                                
                [NMI(ds_id, method_index), NMI_idx(ds_id, method_index)] = max(NMI_K{method_index}(ds_id, :));
                [ARI(ds_id, method_index), ARI_idx(ds_id, method_index)] = max(ARI_K{method_index}(ds_id, :));
            end            
        end
    end

%%
    for ds_id = 1:numel(dataset_names)
        [~, perm] = sort(NMI(ds_id, :), 'descend');
        
        for j = 1:numel(perm)-1
            pop1 = NMI_samples_cell{perm(j)}{ds_id, NMI_idx(ds_id, perm(j))};
            pop2 = NMI_samples_cell{perm(j+1)}{ds_id, NMI_idx(ds_id, perm(j+1))};
            [~, NMI_pvals(ds_id, j)] = ttest2(pop1, pop2, 'tail', 'right');
        end
        
        [~, perm] = sort(ARI(ds_id, :), 'descend');
        for j = 1:numel(perm)-1
            pop1 = ARI_samples_cell{perm(j)}{ds_id, ARI_idx(ds_id, perm(j))};
            pop2 = ARI_samples_cell{perm(j+1)}{ds_id, ARI_idx(ds_id, perm(j+1))};
            [~, ARI_pvals(ds_id, j)] = ttest2(pop1, pop2, 'tail', 'right');            
        end
    end
%%



    for ds_id = 1:numel(dataset_names)
        
        for method_index = 1:numel(methods)        
            method_name = methods{method_index};
            if(strcmp(method_name, 'ACTION'))
                NMI(ds_id, method_index) = mean(NMI_samples_cell{method_index}{ds_id, 1});
                ARI(ds_id, method_index) = mean(ARI_samples_cell{method_index}{ds_id, 1});               
            else                
                for K = 5:5:50                                        
                    NMI_K{method_index}(ds_id, K) = mean(NMI_samples_cell{method_index}{ds_id, K});
                    ARI_K{method_index}(ds_id, K) = mean(ARI_samples_cell{method_index}{ds_id, K});           
                end                                
                [NMI(ds_id, method_index), NMI_idx(ds_id, method_index)] = max(NMI_K{method_index}(ds_id, :));
                [ARI(ds_id, method_index), ARI_idx(ds_id, method_index)] = max(ARI_K{method_index}(ds_id, :));
            end            
        end
    end


%%
sample_range = 5:5:50;
idx = -1;
x_NMI = [];
x_ARI = [];
for ds_id = 1:4
    for method_index = 2:4
        x_NMI(method_index-1, 1:10) = NMI_K{method_index}(ds_id, sample_range);
        x_ARI(method_index-1, 1:10) = ARI_K{method_index}(ds_id, sample_range);

        mean_rank = tiedrank(x_NMI(method_index-1, :)) + tiedrank(x_ARI(method_index-1, :));
        [~, idx] = max(mean_rank);        
        best_k(ds_id, method_index) = sample_range(idx);
    end
    
    plot(sample_range, x_NMI', 'Linewidth', 2);
    set(gca,'FontSize', 12, 'FontWeight','bold');   
    xlabel('Dimension (k)');
    ylabel('NMI');
    legend({'Isomap', 'MDS', 'SIMLR'}, 'Location', 'SouthWest');
    export_fig(fullfile(path.output, strcat(dataset_names{ds_id}, '_dim_NMI.eps')), '-eps', '-transparent', '-painters');    
    close
    
    
    plot(sample_range, x_ARI', 'Linewidth', 2);
    set(gca,'FontSize', 12, 'FontWeight','bold');   
    xlabel('Dimension (k)');
    ylabel('ARI');
    legend({'Isomap', 'MDS', 'SIMLR'}, 'Location', 'SouthWest');
    export_fig(fullfile(path.output, strcat(dataset_names{ds_id}, '_dim_ARI.eps')), '-eps', '-transparent', '-painters');    
    close
end
%%
    exp_name = 'kernelComp_mean';
    
    CC = cbrewer('seq', 'BuGn', 20);

    figure
    data_raw = (round(NMI*100)/100);
    data = (round(zscore(NMI, 0, 2)*1000)/1000);
    heatmap(data_raw, methods, dataset_names, data_raw, 'colormap', CC, 'TextColor', [0.0, 0.0, 0.0], 'FontSize', 12);
    set(gca,'FontSize', 12, 'FontWeight','bold');   
    export_fig(fullfile(path.output, strcat(exp_name, '_NMI.eps')), '-eps', '-transparent', '-painters');

    figure
    data_raw = (round(ARI*100)/100);
    data = (round(zscore(ARI, 0, 2)*100)/100);
    heatmap(data_raw, methods, dataset_names, data_raw, 'colormap', CC, 'TextColor', [0.0, 0.0, 0.0], 'FontSize', 12);
    set(gca,'FontSize', 12, 'FontWeight','bold');   
    export_fig(fullfile(path.output, strcat(exp_name, '_ARI.eps')), '-eps', '-transparent', '-painters');
      
    
    
%%
    save(fullfile(path.output, strcat('kernel_experiment', '.mat')));
