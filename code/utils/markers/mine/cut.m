function [ num_selected_entries ] = cut( v, varargin )
    params = inputParser;
    params.addParamValue('selection_method'           , 'participation_ratio',@(x) ischar(x) ); % selection method. Options:
    params.addParamValue('selection_half'           , 'top',@(x) ischar(x) ); % Should we select the top half or the bottom half of scores?
    params.addParamValue('cut_threshold'           , 100      ,@(x) isscalar(x) & x <= 100 & x > 0); % Percentage of top-ranked scores to choose    
    params.addParamValue('debug'           , false      ,@(x) islogical(x)); 
        
    params.parse(varargin{:});
    par = params.Results;

    badV = isnan(v) | isinf(v);
    v(badV) = [];
    n = numel(v);
    
    switch(par.selection_method)
        case 'Elbow'
            v = sort(v);                
            [ ~, num_selected_entries ] = findElbow_maxDist( v );            
            if(  (num_selected_entries < n/2 && strcmp(par.selection_half, 'top') ) || ... % Right-elbow 
                (num_selected_entries > n/2 && strcmp(par.selection_half, 'bottom') )  ) % Left-elbow
                    num_selected_entries = n - num_selected_entries;
            end
        case 'participation_ratio'
            num_selected_entries = round(sum(abs(v))^2/sum(v.^2)); % participation ratio relaxed (allows more nonzeros)
            if( strcmp(par.selection_half, 'bottom') )
                num_selected_entries = n - num_selected_entries; % Participation ratio returns #nnz=#elements in Top half. If we want bottom half, we have to reverse it.
            end
        case 'participation_ratio_stringent'
            num_selected_entries = round(sum(v.^2)^2/sum(v.^4)); % participation ratio 
            if( strcmp(par.selection_half, 'bottom') )
                num_selected_entries = n - num_selected_entries; % Participation ratio returns #nnz=#elements in Top half. If we want bottom half, we have to reverse it.
            end
        case 'MeanDiff'
            q = nan(length(v)-1, 1);
            for i = 1:length(v)-1
                mu1 = mean(v(1:i)); mu2 = mean(v(i+1:end));
                sigma1 = std(v(1:i)); sigma2 = std(v(i+1:end));
                q(i) = (mu2 - mu1) / (sigma1 + sigma2);
            end
            [~, cut_idx] = min(q);       
            
            if( strcmp(par.selection_half, 'bottom') )
                num_selected_entries = cut_idx; % Participation ratio returns #nnz=#elements in Top half. If we want bottom half, we have to reverse it.
            else
                num_selected_entries = n - cut_idx;
            end   

        case 'ttest'
            tstat = nan(n-4, 1);
            pval = nan(n-4, 1);
            for k = 3:n-2
                [~,p,~,stats] = ttest2(v(1:k), v(k+1:n),'Vartype','unequal');
                pval(k-2) = -log10(p);
                tstat(k-2) = stats.tstat;
                
%                 
%                 mu1(k-2) = mean(v(1:k)); mu2(k-2) = mean(v(k+1:n));
%                 sigma1_sqr = var(v(1:k)); sigma2_sqr = var(v(k+1:n));                
%                 sigma(k) = sqrt( (sigma1_sqr / k) + (sigma2_sqr / (n-k)) );
%                 
%                 delta_mean(k) = mu2(k-2)-mu1(k-2);
%                 
%                 tstat(k-2) = (delta_mean(k)) / sigma(k);
            end
            if(par.debug)
                plot(n-2:-1:3, pval');
            end
            [~, cut_idx] = max(pval);       
            
            if( strcmp(par.selection_half, 'bottom') )
                num_selected_entries = cut_idx+1; 
            else
                num_selected_entries = n - cut_idx-1;
            end  
            
        
        case 'signrank'
            pval = nan(n-4, 1);
            for k = 3:n-2
                [~,p] = signrank(v(1:k), v(k+1:n),'Tail','Smalle');
                pval(k-2) = p;                
            end
            if(par.debug)
                plot(n-2:-1:3, -log10(pval)');
            end
            [~, cut_idx] = min(pval);       
            
            if( strcmp(par.selection_half, 'bottom') )
                num_selected_entries = cut_idx+1; 
            else
                num_selected_entries = n - cut_idx-1;
            end              

        case 'ranksum'
            pval = nan(n-4, 1);
            for k = 3:n-2
                p = ranksum(v(1:k), v(k+1:n));
                pval(k-2) = p;                
            end
            if(par.debug)
                plot(n-2:-1:3, -log10(pval)');
            end

            [~, cut_idx] = min(pval);       
            
            if( strcmp(par.selection_half, 'bottom') )
                num_selected_entries = cut_idx+1; 
            else
                num_selected_entries = n - cut_idx-1;
            end              
            
        case 'Gaussian'
%             pisqrt = sqrt(2*pi);
            L = nan(n-2, 1);
            for k = 2:n-1
                mu1 = mean(v(1:k)); mu2 = mean(v(k+1:n));
                sigma1 = std(v(1:k)); sigma2 = std(v(k+1:n));
                
                L1 = sum(log(normpdf(v(1:k), mu1, sigma1)));
                L2 = sum(log(normpdf(v(k+1:n), mu2, sigma2)));
                L(k-1) = L1 + L2;
%                 z1 =(v(1:k) - mu1).^2 / (2*sigma1^2); z2 = (v(k+1:n) - mu2).^2 / (2*sigma2^2);
%                 
%                 L1 = k* log(1 / (sigma1*pisqrt))*sum( z1 );
%                 L2 = (k-n)*log(1 / (sigma2*pisqrt))*sum( z2 );
%                 L(k-1) = -(L1 + L2);
            end
            if(par.debug)
                plot(2:n-1, L');
            end
            [~, cut_idx] = max(L);       
            
            if( strcmp(par.selection_half, 'bottom') )
                num_selected_entries = cut_idx+1; 
            else
                num_selected_entries = n - cut_idx-1;
            end   

            
            
        case 'AIC' % Unstable: needs work
            v = sort(v);                
            num_selected_entries = AIC_outlierDet(v, 0.9);
        case 'fix_threshold' 
            num_selected_entries = nnz(v >= par.cut_threshold);
        case 'fix_perc'
            num_selected_entries = round(par.cut_threshold*numel(v)/100);
        otherwise 
            error('Cell type selection method unknown: %s\n', par.selection_method);
    end  
end

