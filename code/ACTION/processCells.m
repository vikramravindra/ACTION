function [ Z, final_C, archetypes, final_S, k, cell_types, trace ] = processCells( X, varargin )
    rng('default') % Only useful if we don't use SPA for initialization
    [m, n] = size(X);

    params = inputParser;    
    params.addParamValue('large_scale', false, @(x) islogical(x)); % should we sacrifice accuracy for speed?
    
    
    params.addParamValue('batch_annotations', {}, @(x) iscell(x)); % Annotations for each cell to correct for confounding factors
    params.addParamValue('annotations', repmat({'1'}, n, 1), @(x) iscell(x)); % trusted annotations, for labeling archetypes
    
    params.addParamValue('debug', false, @(x) islogical(x));    

    params.addParamValue('k', 10, @(x) isscalar(x) & x > 0 & x <= min(m, n)); % Set k (#archetypes) to a given value
    params.addParamValue('delta', 0, @(x) isscalar(x) & x >= 0); % margin of error for PCHA
    params.addParamValue('max_k', -1, @(x) isscalar(x) & x > 0 & x <= min(m, n)); % Maximum k to consider for automatic k identification (equivalent to the maximum number of estimated tasks/functions)

    params.addParamValue('precond', 'minvol', @(x) ischar(x)); % methods to use for preconditioning: none, prewhitening, minvol

    params.addParamValue('ER_pval_threshold', 1e-3, @(x) isscalar(x) & x >= 0 & x <= 1); 
    

    
    % Only for post-processing label assignments
    params.addParamValue('unassign_noisy_cells', false, @(x) islogical(x)); % Should we filter final cells according to their error in their measurements?
    params.addParamValue('max_acceptable_cell_error', 1, @(x) isscalar(x)); % how much relative residual error is accepted in estimating X(:, i) - W(i, :)*H(i, :)' 
    params.addParamValue('min_membership_score', 0, @(x) isscalar(x) & x <= 1 & x >= 0); % if value of membership matrix is less than this, then we ignore it
    
    params.addParamValue('output_path', 'output', @(x) ischar(x)); % Folder to store intermediate results
    params.addParamValue('experiment_name', 'default', @(x) ischar(x)); % Folder to store intermediate results
    
    params.parse(varargin{:});
    par = params.Results;

    if( ~exist(par.output_path, 'dir') )
        system( sprintf('mkdir %s -p', par.output_path) ); % make sure it works properly in Windows
    end

    
    fprintf('Processing %d cells each having %d genes ... \n', n, m);

    % Normalize expression profiles such that each column have norm_1=1    
    X = normalize(X, 'pnorm', 1, 'dim', 1);
    
    tic
    fprintf('\tRunning OPA (%d confounding factors) ... ', size(par.batch_annotations, 2));
    Y = OPA(X, 'factors', par.batch_annotations); % Y'Y ./ m; This should be equivalent to partialcorr(X, par.hk), if there is no confounding factor;    
    clear X;
    fprintf('Done (dt = %.2f s)\n', toc);
    
    
    tic
    fprintf('\tComputing expression specificity of genes ... ');
    w = genesSpecificity(Y);
    fprintf('Done (dt = %.2f s)\n', toc);

    tic
    fprintf('\tComputing final Kernel ... ');
    Z = full(bsxfun(@times, Y, w));    
    Ker = Z'*Z;
    fprintf('Done (dt = %.2f s)\n', toc);

    
%     if(par.k > 0 && par.k <= n)
%         fprintf('\tSetting archetype count (k) = %d\n', par.k);
%         k = par.k;
%         trace = [];
%         G = [];
%         node_labels = [];
%     else
%         fprintf('\tUnknown archetype count (k)\n'); 
%         if(~par.large_scale)
%             [ k, trace, G, node_labels ] = estimateK(Z, 'max_k', par.max_k, 'annotations', par.annotations, 'visualize', par.visualize, 'output_path', par.output_path, 'experiment_name', par.experiment_name);        
%         else
%             % use eigenpairs to estimate k. More efficient but much less effective
%             k = round(log(n));
%             trace = [];
%             G = [];
%             node_labels = {};
%         end
% 
%         fprintf('\tEstimated k = %d\n', k);
%     end
    

    tic
    fprintf('\tSolving near-separable NMF ... ');    
    if(par.max_k == -1)
%         [trace.K] = adaptiveSPA_naive(Z, 'r', par.k, 'locality_pval', 0.05, 'debug', true, 'window_size', 0, 'true_labels', par.annotations, 'precond', 'minvol', 'scale', true);        
%         [trace.K, trace.best_k, trace.smoothed_logPvals, trace.log] = adaptiveSPA_mutual_trace(Z, 'r', par.k, 'locality_pval', 0.05, 'debug', true, 'window_size', 0, 'true_labels', par.annotations, 'precond', 'prewhitening', 'scale', true);        
%         [~, trace.K] = VCA(Z, 'Endmembers', par.k, 'verbose', 'on'); 
%           trace.K = FastConicalHull(Z,par.k); trace.K = trace.K'; 
        [trace.K] = adaptiveSPA_naive(Z, 'r', par.k, 'debug', true, 'true_labels', par.annotations, 'precond', par.precond, 'scale', true);
        k = par.k;
    else
        [trace.K, trace.log] = adaptiveSPA_localityER(Z, 'r', par.max_k, 'scale', true, 'precond', par.precond, 'true_labels', par.annotations, 'debug', par.debug, 'pval_threshold', par.ER_pval_threshold); 
        k = trace.log.best_k;
    end
    K = trace.K(1:k);
    fprintf('\tSetting k = %d\n',  k);



%     k = par.k;
%     [K] = PrecSPA(Z,k,1,1,1); 
%     trace = [];
    


    
%     if(par.prewhitening && ~par.large_scale)
%         [~,~,vh] = svds(Z, k); % pre-whitening    
%         K = FastSepNMF(vh',k); % Initial points on the convext hull of data points using Successive Projection Algorithm (SPA)
% %         K = FastSepNMF(Z,k); % Initial points on the convext hull of data points using Successive Projection Algorithm (SPA)
%     else
%         % SVD is expensive! Ignore whitenning
%         K = FastSepNMF(Z,k); % Initial points on the convext hull of data points using Successive Projection Algorithm (SPA)
%     end
    fprintf('Done (dt = %.2f s)\n', toc);

    
    
    
%     final_C = sparse(K, 1:k, 1, n, k);    
%     archetypes = Z*final_C;    
    
%     cvx_begin quiet
%     variable H(k, n)
%     minimize( norm(archetypes*H - Z, 'fro') )
%     subject to
%         0 <= H <= 1
%         ones(1, k)*H == ones(1, n)
%     cvx_end


%     H = zeros(k, n);
%     lsqlin_options = optimset('LargeScale','off', 'Display', 'off');            
%     for sample_id=1:n
%         fprintf('%d/%d\n', sample_id, n);
%         H(:, sample_id) = lsqlin(archetypes, Z(:, sample_id), [], [], ones(1, k), 1, zeros(k, 1), ones(k, 1), zeros(k, 1), lsqlin_options);
%     end 

%      H = nnlsm_activeset(Z*final_C, Z); final_S = normalize(final_S, 'pnorm', 1, 'dim', 1);
            
%     final_S = H;
%     [~, cell_types] = max(final_S);
                
    
    
    fprintf('\tAdjusting archetypes ... \n');
    opts.C = sparse(K, 1:k, 1, n, k);    
    opts.verbose = 1;
    [final_S, final_C] = PCHAkernel(Ker, k, 1:n, 1:n, par.delta, opts);    

    archetypes = Z*final_C;    

    unassigned_cells = find(sum(final_S) == 0);
    [~, idx] = max(cell2mat(arrayfun(@(curr_cell) Z(:, curr_cell)'*Z(:, setdiff(1:n, curr_cell)), unassigned_cells, 'uniformoutput', false)'), [], 2);
    final_S(:, unassigned_cells) = final_S(:, idx);
    
    [~, cell_types] = max(final_S);
        
    
    % Fitering unreliable assignments
%     final_S(:, max(final_S) < par.min_membership_score) = 0;   % Remove cell type associations for cells that are not really close to any archetype (if requested) 
%     Err = abs(X - X*final_C*final_S);
%     total_err = sum(sum(Err.^2)); % Squared norm frob.
%     if(par.unassign_noisy_cells)
%         cell_err = sum(Err, 1) ./ sum(abs(X));
%         final_S(:, cell_err <= par.max_acceptable_cell_error) = 0;
%     end    
%     final_S = normalize(final_S, 'pnorm', 1, 'dim', 1);

        

%     archetypes_adjusted = 1 ./ (1 + exp(-archetypes));
end

