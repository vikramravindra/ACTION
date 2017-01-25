function [ weights ] = genesSpecificity( X, varargin )
    params = inputParser;    
    params.addParamValue('order_limit', 2, @(x) isscalar(x) & x >= 0); 
    params.parse(varargin{:});
    par = params.Results;

    [~, cell_count] = size(X);
    
    X(X < 0) = 0;
    row_min = min(X, [], 2);
    X_shifted = bsxfun(@minus, X, row_min);

    row_sum = sum(X_shifted, 2);
    P = sparse(bsxfun(@rdivide, X_shifted, row_sum));
    
    I = spfun(@(x) -log2(x), P); % Information content of the relative expression of each gene   
    H = sum(P .* I, 2); % Shannon entropy of the relative expression -- overall tissue-specificity of genes. It has unit of "bits" and is between 0 -> completely selective) and log2(k) -> Uniform/HK gene    
    uniformity  = H/log2(cell_count); % Scale between 0-1

    % Identify the cutting point to distinguish specific genes
    v = uniformity - min(uniformity);
    v = sort(v);    
    cut_idx= cut(v, 'selection_method', 'participation_ratio_stringent', 'selection_half', 'bottom');
    transition_val = v(cut_idx);
    transition_val = transition_val + min(uniformity);
    weights = full(transition_val ./ uniformity);

    % Limit weights to ensure stability
    weights(weights < 10^(-par.order_limit)) = 10^(-par.order_limit);
    weights(weights > 10^(par.order_limit)) = 10^(par.order_limit);    
end

