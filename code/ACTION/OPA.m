function [ Y, X_z ] = OPA( X, varargin )
%% Input) X: Input expression matrix.
% <'housekeeping', h_vec>: h_vec is the signature vector of housekeeping genes
% <'fators', factor_table>: factor_table is a cell array with each column
% representing a categorical confounding factor 
% Output) Y: Adjusted transcriptional signatures

    X = full(X);

    params = inputParser;    
    params.addParamValue('housekeeping', mean(X, 2), @(x) isvector(x)); % vector of gene expression for housekeeping genes
    params.addParamValue('factors', {}, @(x) iscell(x)); % table of confounding factors   
    params.parse(varargin{:});
    par = params.Results;

    [~, cell_count] = size(X);
%    
%     [u, ~, ~] = svds(X, 3);
%     par.housekeeping = u;
%     
%     if(sign(sum(par.housekeeping(:, 1))) ~= sign(sum(sum(X))))
%         par.housekeeping = -par.housekeeping;
%     end

    
    if(~isempty(par.factors))
        % Construct basis for the confounding subspace
        sizes = arrayfun(@(col) numel(unique(par.factors(:, col))), 1:size(par.factors, 2));
        selected_attributes = find( (sizes < round(cell_count/10)) & (1 < sizes) ); % Ignore overly specific factors or factors with no specificity

        attr_no = numel(selected_attributes);
        [~, perm] = sort(sizes(selected_attributes));
        selected_attributes = selected_attributes(perm); % From the most common attribute to more specific ones

        correction_basis = cell(1, numel(selected_attributes));
        for i = 1:attr_no
            curr_attr = selected_attributes(i);
            unique_attr = unique(par.factors(:, curr_attr));


            [~, Labels] = ismember(par.factors(:, curr_attr), unique_attr);
            correction_basis{i} = cell2mat(arrayfun(@(label) trimmean(X(:, Labels == label), 25, 2), 1:numel(unique_attr), 'UniformOutput', false));        
        end

        % Combine confounding and housekeeping subspaces
        Basis_z = my_zscore([par.housekeeping, cell2mat(correction_basis)]);
%         Basis_z = zscore([par.housekeeping, cell2mat(correction_basis)]);
    else
        Basis_z = my_zscore(par.housekeeping);
%         Basis_z = zscore(par.housekeeping);
    end
    
    % Project to the orthogonal subspace
    X_z = my_zscore(X); % zscore(X); % Raw signatures
    clear X;
    
    X_r = orthoProject( X_z, Basis_z );        
    clear X_z Basis_z
    
    % Standardize
    Y = my_zscore(X_r);%    zscore(X_r); % Adjusted signatures 
    clear X_r;
end

