function [K, trace] = adaptiveSPA_localityER(X, varargin) 
    [~, n] = size(X);
    
    params = inputParser;    
    params.addParamValue('pval_threshold', 0.001, @(x) isscalar(x) & x >= 0 & x<= 1); 
    params.addParamValue('z_threshold', 1.96, @(x) isscalar(x) & x >= 0 & x<= 1); 
    params.addParamValue('r', min(10, n), @(x) isscalar(x) & x > 0 & x <= n);
    params.addParamValue('true_labels', ones(n, 1), @(x) isvector(x));
    params.addParamValue('debug', false, @(x) islogical(x));    
    params.addParamValue('precond', 'none', @(x) ischar(x));    
    params.addParamValue('scale', true, @(x) islogical(x));        
    params.parse(varargin{:});
    par = params.Results;

    Ker = X'*X;
    
    if(iscell(par.true_labels))
        UL = unique(par.true_labels);
        [~, par.true_labels] = ismember(par.true_labels, UL);
    end

    switch(par.precond)
        case 'none'
            M = X;

        case 'TnS'
            tic; [~, S, V] = svd(X, 'econ'); toc
            M = S*V'; % U'*X
            
        case 'prewhitening'
            [~, ~, V_r] = svds(X, par.r); 
            M = V_r';
            
        case 'minvol'
            [~, S_r, V_r] = svds(X, par.r); 
            Xr = S_r*V_r';
            A_star = minvolell(Xr, 1e-6, par.r, par.debug); 
            P = chol(A_star);
            M = P*Xr;
            
        otherwise
            error('Unknown preconditioning %s', par.precond);
    end            
    
    if(par.scale)
        M = normalize(M, 'pnorm', 1, 'dim', 1);
    end
    
    normM = sum(M.^2); 
    nM = max(normM); 
    normM1 = normM;
    
    i = 1; 
    trace.meta_pval = ones(par.r, 1);
    K = zeros(par.r, 1);
    U = zeros(size(M, 1), par.r);    
    trace.er_pval = ones(par.r, 1);
    
    while i <= par.r && max(normM)/nM > 0
        % Select the column of M with largest l2-norm
        [a,b] = max(normM); 

        % Check ties up to 1e-6 precision
        b = find((a-normM)/a <= 1e-2); 

        if length(b) > 1
            if( i <= 1 )
                [~,d] = max(normM1(b)); 
            else                
                CC = corr(X(:, b), X(:, K(1:(i-1))));
                bestMatch = max(CC, [], 2);
                [~, d] = min(bestMatch);
            end
            b = b(d); 
        end

        % Update the index set, and extracted column
        K(i) = b; U(:,i) = M(:,b);         

        if(i > 1)                        
            subKer = Ker(K(1:i), K(1:i));
            ind = find(tril(ones(size(subKer)), -1));
            vals = subKer(ind);
            z = zscore(vals);
            
            mask = (z > par.z_threshold);
            G = zeros(size(subKer));
            G(ind(mask)) = 1;
            G = max(G, G');
            density = nnz(G) / (i*(i-1));

            
            [cc] = components(G);
            comp_no = max(cc);
            nontrivial_comps = find(arrayfun(@(c_id) nnz(cc == c_id), 1:comp_no) > 1);
            if(~isempty(nontrivial_comps))
                comp_pvals = ones(numel(nontrivial_comps), 1);
                for j = 1:numel(nontrivial_comps)
                    mc_vertices = find(cc == nontrivial_comps(j));
                    mc_density = nnz(G(mc_vertices, mc_vertices)) / (numel(mc_vertices)*(numel(mc_vertices)-1));

                    comp_pvals(j) = ER_DensityPval(i, density, mc_density, numel(mc_vertices));

                end
                trace.er_pval(i) = min(comp_pvals*numel(nontrivial_comps));
            end            
            
            

        else
            trace.meta_pval(i) = 1;
        end
        

        
        
        % Compute (I-u_{i-1}u_{i-1}^T)...(I-u_1u_1^T) U(:,i), that is,
        % R^(i)(:,K(i)), where R^(i) is the ith residual (with R^(1) = M).
        for j = 1 : i-1
            U(:,i) = U(:,i) - U(:,j)*(U(:,j)'*U(:,i));
        end

        % Normalize U(:,i)
        U(:,i) = U(:,i)/norm(U(:,i));     

        % Compute v = u_i^T(I-u_{i-1}u_{i-1}^T)...(I-u_1u_1^T)
        v = U(:,i); 
        for j = i-1 : -1 : 1
            v = v - (v'*U(:,j))*U(:,j); 
        end

        % Update the norm of the columns of M after orhogonal projection using
        % the formula ||par.r^(i)_k||^2 = ||par.r^(i-1)_k||^2 - ( v^T m_k )^2 for all k. 
        normM = normM - (v'*M).^2;   
        normM(b) = -inf;
        
        i = i+1;
    end

    if(i < par.r)
        K(i:end) = [];
    end
    
   trace.er_pval_adj = min(1, (i-1)*trace.er_pval);
   trace.best_k = find(trace.er_pval_adj < par.pval_threshold, 1, 'first')-1;
   
   if(par.debug)
       plot(-log10(trace.er_pval), 'LineWidth', 2);
       xlabel('k');
       ylabel('-log(pval)');
   end

end
    