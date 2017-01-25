function [ arch_annotations, coverage ] = archLabelEnrichment( M, Labels, varargin )
    params = inputParser;
    params.addParamValue('pval_threshold', 0.05, @(x) isscalar(x) & x >= 0 & x <= 1); % pvalue threshold
    
    params.parse(varargin{:});
    par = params.Results;
    
    [~, l] = max(M);
    clusters = Labels2Clusters(l);    
    cluster_no = numel(clusters);
        
    UL = unique(Labels);
    Label_counts = arrayfun(@(i) nnz(strcmp(Labels, UL{i})), 1:numel(UL));
    
    pop_size = numel(Labels);
    arch_annotations = cell(cluster_no, 1);
    annotation_mask = zeros(cluster_no, numel(UL));
    for i = 1:cluster_no
        sig_labels = {};
        samples = clusters{i};        
        for k = 1:numel(UL)
            success = nnz(strcmp(Labels(samples), UL{k}));
            pval = hygecdf(success-1, pop_size, numel(samples), Label_counts(k), 'upper');
            if(pval < par.pval_threshold)
                sig_labels = union(sig_labels, UL{k});
                annotation_mask(i, k) = 1;
            end                            
        end        
        arch_annotations{i} = strjoin(sig_labels, ' || ');        
    end
    coverage = [UL, num2cell([Label_counts; sum(annotation_mask, 1)]')];
end

