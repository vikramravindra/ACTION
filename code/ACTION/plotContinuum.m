function [ fig, ydata, iid_points ] = plotContinuum( X, varargin )
    [n, m] = size(X);

    params = inputParser;    
    params.addParamValue('Labels', repmat({'none'}, m, 1), @(x) iscell(x)); % Cell array of text labels, one for each cell (used for coloring the plot)
    params.addParamValue('arch_labels', {}, @(x) iscell(x)); % Cell array of text labels, one for each archetype (used to create figure legend)
    params.addParamValue('fiedler', true, @(x) islogical(x));  % Should we use Fiedler embedding as initialization? Highly suggested to be true
    params.addParamValue('dim', 2, @(x) x <= 3 & x > 1);  % Dimension of projection space (2D or 3D)
    params.addParamValue('perplexity', 30, @(x) isscalar(x)); % parameter to tSNE
    params.addParamValue('tsne_iters', 1000, @(x) isscalar(x) & x > 0);  % How many iterations of tSNE should we use?
    params.addParamValue('coloring', 'convex', @(x) ischar(x));  % either 'convex' (combination of archetype colors) or 'labels' (known cell types from par.Labels)
    params.parse(varargin{:});
    par = params.Results;
            
    UL = unique(par.Labels);
    [~, par.labels] = ismember(par.Labels, UL);
    
    par.labels = [par.labels; (max(par.labels)+1)*ones(n,1)];
    rng('default'); % Useless when we use Fiedler initialization
    
    X = normalize(X, 'pnorm', 1, 'dim', 1); % Just to make sure .. !
    A = [X, eye(n)];

    Dist = squareform(pdist(A', 'euclidean'));

    iid_points = logical(sum(Dist(1:size(X, 2), size(X, 2)+1:end) == 0, 2));
    Dist(iid_points, :)= [];
    Dist(:, iid_points)= [];
    labels = [par.labels(~iid_points); ones(n, 1)*max(par.labels)+1];
    m = size(Dist, 1)-n;
    X(:, iid_points) = [];

    if(par.fiedler)
        W = affinityMatrix(Dist, round(m/10), 1);
    %     W = d2p(Dist .^ 2, par.perplexity, 1e-5); 

        L = diag(sum(W)) - W;    
        [V, eigenval] = eigs(L, par.dim+2, 'sm'); % first one is zero, last one is only used for computing eigen-gap
        %eigen_gap = eigenval(1,1)-eigenval(2,2);    
        U = V(:, 2:par.dim+1);
        F = bsxfun(@rdivide, U, 1 ./ sqrt(diag(eigenval(2:par.dim+1, 2:par.dim+1)))');
        scatter(F(:, 1), F(:, 2), 10, labels);

        ydata = tsne_d(Dist, labels, F, par.perplexity, par.tsne_iters);
    else
        ydata = tsne_d(Dist, labels, par.dim, par.perplexity, par.tsne_iters);
    end
    

    fig = figure;
    set(fig, 'Position', [0 0 1000 1000])
    hold all
    xlabel('X','FontSize', 16, 'FontWeight','bold');
    ylabel('Y','FontSize', 16, 'FontWeight','bold');
    if(par.dim == 3)
        zlabel('Z','FontSize', 16, 'FontWeight','bold');
    end        

    if(strcmp(par.coloring,'convex'))
        C = cbrewer('qual', 'Paired', n);
    %     C = cbrewer('qual', 'Dark2', n);

        for i = 1:n
            if(par.dim == 2)
                scatter(ydata(m+i, 1), ydata(m+i, 2), 75, C(i, :), 'filled');
            else
                scatter3(ydata(m+i, 1), ydata(m+i, 2), ydata(m+i, 3), 75, C(i, :), 'filled');
            end                
        end


        if(~isempty(par.arch_labels))
            legend(arrayfun(@(i) sprintf('Archetype %d: %s', i, par.arch_labels{i}), 1:n, 'Uniformoutput', false), 'Location', 'SouthOutside');
        end

        for i = 1:m
            CC = 0;
            for j = 1:n
                CC = CC + X(j, i) * C(j, :);
            end
            CC = max(min(1, CC), 0);
            if(par.dim == 2)
                scatter(ydata(i, 1), ydata(i, 2), 15, CC, 'filled');
            else
                scatter3(ydata(i, 1), ydata(i, 2), ydata(i, 3), 15, CC, 'filled');
            end
        end    
    elseif(strcmp(par.coloring, 'labeled'))        
        C = cbrewer('qual', 'Set1', numel(UL));          
        
        for ll = 1:numel(UL)            
            if(par.dim == 2)
                scatter(ydata(labels == ll, 1), ydata(labels == ll, 2), 15, C(ll, :), 'filled');
            else
                scatter3(ydata(labels == ll, 1), ydata(labels == ll, 2), ydata(labels == ll, 3), 15, C(ll, :), 'filled');
            end                
        end         
        legend(UL, 'EdgeColor', 'k');             
    end
    
    if(par.dim == 3)
        view(3);
    end

    for i = 1:n
        if(par.dim == 2)
            text(ydata(m+i, 1), ydata(m+i, 2), sprintf('A%d', i),'FontSize', 14, 'FontWeight','bold')            
        else
            text(ydata(m+i, 1), ydata(m+i, 2), ydata(m+i, 3), sprintf('A%d', i),'FontSize', 14, 'FontWeight','bold')
        end                
    end        
    
    set(gca,'FontSize', 12, 'FontWeight','bold'); 
    
    
    Min = min(ydata);
    Max = max(ydata);
    R = (Max - Min) / 7;
    set(gca, 'XLim', [Min(1)-R(1), Max(1)+R(1)])
    set(gca, 'YLim', [Min(2)-R(2), Max(2)+R(2)])
    if(par.dim == 3)
        set(gca, 'ZLim', [Min(3)-R(3), Max(3)+R(3)])
    end
  
end

