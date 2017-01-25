function S = getSimilarity(X, method_name, varargin)

params = inputParser;
params.addParamValue('dim', 20, @(x) isscalar(x) & x > 0);
params.addParamValue('cluster_no', -1, @(x) isscalar(x) & x > 0);
params.addParamValue('kappa', 1, @(x) isscalar(x) & x > 0);
params.parse(varargin{:});
par = params.Results;

switch(method_name)
    
    case 'SIMLR'
        [~, S] = SIMLR(X, par.cluster_no, par.dim, 0, 1);        
        
    case 'Isomap'
        [mappedX] = compute_mapping(X', 'Isomap', par.dim);
        S = corr(mappedX');

    case 'MDS'
        [mappedX] = compute_mapping(X', 'MDS', par.dim);
        S = corr(mappedX');

        
    case 'ACTION'
        X = normalize(X, 'pnorm', 1);
        Y = OPA(X); 
        
        w = genesSpecificity(Y);
        Z = bsxfun(@times, Y, w);     
        S = Z'*Z;     
        
    otherwise
        error('Unknown method %s', method_name);            
end
end