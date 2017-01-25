function [ arch_annotations ] = archEnrichment( M, annotations, varargin )
    params = inputParser;
    params.addParamValue('pval_threshold', 0.05, @(x) isscalar(x) & x >= 0 & x <= 1); % pvalue threshold
    params.addParamValue('overlapping', false, @(x) islogical(x)); 

    params.parse(varargin{:});
    par = params.Results;

    
    if(par.overlapping)
        [X, X_perm] = sort(M, 'descend');
        cut_size = round(sum(X).^2./sum(X.^2))';

        selected_node_clusters = arrayfun(@(c) X_perm(1:cut_size(c), c), 1:size(X, 2), 'UniformOutput', false)';         
        
        MM = zeros(size(M));
        for i = 1:size(MM, 2)
            MM(selected_node_clusters{i}, i) = 1;
        end
        clusters = Membership2Cluster(MM');        
    else        
        [~, l] = max(M);
        clusters = Labels2Clusters(l);    
    end
    cluster_no = numel(clusters);
    
    sizes = arrayfun(@(col) numel(unique(annotations(:, col))), 1:size(annotations, 2));
    selected_attributes = find( (sizes < size(annotations, 1)) & (1 < sizes) );
    attr_no = numel(selected_attributes);

    Attributes = annotations(:, selected_attributes);
    
    uniqAttr = cell(attr_no, 1);
    uniqAttr_counts = cell(attr_no, 1);
    for j = 1:numel(selected_attributes)
        uniqAttr{j} = unique(Attributes(:, j));
        counts = zeros(numel(uniqAttr{j}), 1);
        for k = 1:numel(counts)
            counts(k) = nnz(strcmp(Attributes(:, j), uniqAttr{j}{k}));
        end
        uniqAttr_counts{j} = counts;
    end
    
    pop_size = size(annotations, 1);
    arch_annotations = cell(cluster_no, 1);
    for i = 1:cluster_no
        samples = clusters{i};        
        
        labels = cell(1, numel(selected_attributes));
        for j = 1:numel(selected_attributes)
            sig_vals = {};
            for k = 1:numel(uniqAttr{j})
                success = nnz(strcmp(Attributes(samples, j), uniqAttr{j}{k}));
                
                pval = hygecdf(success-1, pop_size, numel(samples), uniqAttr_counts{j}(k), 'upper');
                if(pval*numel(uniqAttr{j})*numel(uniqAttr{j}) < par.pval_threshold)
                    sig_vals = union(sig_vals, uniqAttr{j}{k});
                end                
            end
            labels{j} = sprintf('[%s]', strjoin(sig_vals, ','));
        end
        arch_annotations{i} = strjoin(labels, '-');
    end
    
    
    cluster_labels = cell(cluster_no+1, 1);
    cluster_labels{1} = '-';
    for i = 2:cluster_no+1
        cluster_labels{i} = sprintf('%d) %s', i, arch_annotations{i-1});
    end          

end

