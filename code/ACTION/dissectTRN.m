function [residual_expression, sorted_genes, TFs, mHG_pvals, mHG_pvals_corrected, selected_genes, TF_table, subTRN_EdgeList, subTRN_NodeAnnotations] = dissectTRN( archetypes, gene_names, TRN_edgeList, Arch_range, varargin )
    params = inputParser;
    params.addParamValue('min_targets', 5, @(x) isscalar(x));
    params.addParamValue('pval_threshold', 0.05, @(x) isscalar(x));
    params.addParamValue('z_threshold', 1.65, @(x) isscalar(x));
    params.addParamValue('orthogonalize', true, @(x) islogical(x));           
    params.parse(varargin{:});
    par = params.Results;

    arch_no = size(archetypes, 2);

%     % Get rid of ribosomal proteins
%     Ribo_mask = (strncmpi(gene_names, 'RPL', 3) | strncmpi(gene_names, 'RPS', 3));
%     archetypes(Ribo_mask, :) = [];
%     gene_names(Ribo_mask, :) = [];

%     archetypes_adjusted = archetypes;

%     archetypes_adjusted(archetypes_adjusted < 0) = 0;
%     archetypes_adjusted = size(archetypes, 1)*normalize(archetypes_adjusted, 'pnorm', 1, 'dim', 1);
%     archetypes_adjusted(archetypes_adjusted < 0) = 0;    
%     if(max(archetypes) > 32)        
%         archetypes_adjusted = log2(archetypes_adjusted);
%         archetypes_adjusted(isinf(archetypes_adjusted)) = 0;
%     end

    archetypes(archetypes < 0) = 0;
    
    if(par.orthogonalize)
        residual_expression = zeros(size(archetypes));
        for i = 1:arch_no
            x = archetypes(:, i);
            A = archetypes(:, setdiff(1:arch_no, i));
            residual_expression(:, i) = x - A/(A'*A)*(A'*x);
    %         v = mean(A, 2);
    %         residual_expression(:, i) = x - v/(v'*v)*(v'*x);
        end

%     residual_expression = zscore(residual_expression);
%     residual_expression = OPA(archetypes_adjusted);
%     w = genesSpecificity(residual_expression);
%     residual_expression = bsxfun(@times, residual_expression, w);
%     

%     residual_expression = archetypes_adjusted;
%     residual_expression = tiedrank(archetypes) + tiedrank(residual_expression);    
    else
        residual_expression = archetypes;
    end
    
    
    [~, perm] = sort(residual_expression, 'descend');
    sorted_genes_cell = arrayfun(@(col) gene_names(perm(:, col)), 1:size(perm, 2), 'Uniformoutput', false);   
    sorted_genes = [sorted_genes_cell{:}];
   

    
    TFs = unique(TRN_edgeList(:, 1));    
    


    
    mHG_pvals = ones(arch_no, numel(TFs));
    selected_genes = cell(arch_no, numel(TFs));
    for arch_idx = 1:numel(Arch_range)
        i = Arch_range(arch_idx);
        
        fprintf('Archetype %d (%d/%d)\n', i, arch_idx, numel(Arch_range));
%         [~, perm] = sort(residual_expression(:, i), 'descend');
%         sorted_gene_names = gene_names(perm);
        perc = 0;
        for k = 1:numel(TFs)   
            if(round(100*k/numel(TFs)) > perc)
                perc = round(100*k/numel(TFs));
                fprintf('%d %%\n', perc);
            end
            rows = strcmp(TRN_edgeList(:, 1), TFs(k));
            targets = unique(TRN_edgeList(rows, 2));
%             if(numel(targets) < par.min_targets)
%                 continue;
%             end
            labels = ismember( sorted_genes(:, i), targets);
%             [v1, p1] = mHG(labels, 'exact');
%             [v2, p2] = mHG(flipud(labels), 'exact');
%             if (p1 <= p2)
%                 HGT_vec = v1;
%                 mHG_pvals(i, k) = p1;
%             else
%                 HGT_vec = v2;
%                 mHG_pvals(i, k) = p2;
%             end
            pos_count = size(residual_expression, 1); 
%             pos_count = nnz(residual_expression(:, i) > 0);
%             pos_count = round(0.5*size(residual_expression, 1));
%             pos_count = 1000;
            [HGT_vec, mHG_pvals(i, k)] = mHG(labels, 'exact', pos_count);
            [~, idx] = min(HGT_vec);
            selected_genes{i, k} = intersect(sorted_genes(1:idx, i), unique(TRN_edgeList(rows, 2)));            

%             plot(-log10(HGT_vec), 'LineWidth', 3)
%             set(gca,'FontSize', 14, 'FontWeight','bold');            
%             xlabel('Rank','FontSize', 18, 'FontWeight','bold');
%             ylabel('Enrichment','FontSize', 18, 'FontWeight','bold');
%             
%             ylim([0, ceil(-log10(min(HGT_vec)))]);
%             set(gca, 'xscale', 'log');
% 
%             line([idx, idx], [0, -log10(HGT_vec(idx))], 'color', [0.5, 0.5, 0.5], 'Linewidth', 3, 'LineStyle', '--')            
%             line([1,idx], [-log10(HGT_vec(idx)), -log10(HGT_vec(idx))], 'color', [0.5, 0.5, 0.5], 'Linewidth', 3, 'LineStyle', '--')            
        end        
    end           
    mHG_pvals(mHG_pvals < 0) = 0;
    

    
%     mHG_pvals_corrected = mHG_pvals;

%     [~, ~, padj] = fdr(mHG_pvals(:));
%     mHG_pvals_corrected = reshape(padj, size(mHG_pvals, 1), size(mHG_pvals, 2));
%     
    mHG_pvals_corrected = ones(size(mHG_pvals));
    for arch_idx = 1:numel(Arch_range)
        i = Arch_range(arch_idx);
        
        pvals = mHG_pvals(i, :);
        [~, ~, mHG_pvals_corrected(i, :)] = fdr(pvals);
    end


    [sorted_mHG_pvals, sorted_mHG_pvals_idx] = sort(mHG_pvals_corrected, 2);
    cut_size = sum(sorted_mHG_pvals < par.pval_threshold, 2);
    TF_table = arrayfun(@(row) [TFs(sorted_mHG_pvals_idx(row, 1:cut_size(row))), num2cell(sorted_mHG_pvals(row, 1:cut_size(row)))'], 1:numel(cut_size), 'UniformOutput', false);


    TF_nodes = arrayfun(@(row) TFs(sorted_mHG_pvals_idx(row, 1:cut_size(row))), 1:numel(cut_size), 'UniformOutput', false);
    TG_nodes = arrayfun(@(row) selected_genes(row, (sorted_mHG_pvals_idx(row, 1:cut_size(row))))', 1:numel(cut_size), 'UniformOutput', false);
    TG_nodes = cellfun(@(x) unique(vertcat(x{:})), TG_nodes, 'UniformOutput', false);

    subTRN_EdgeList = cell(size(mHG_pvals, 1), 1);
    subTRN_NodeAnnotations = cell(size(mHG_pvals, 1), 1);
    for arch_idx = 1:numel(Arch_range)
        i = Arch_range(arch_idx);
        
        if(isempty(TF_nodes{i}) || isempty(TG_nodes{i}))
            continue;
        end
        row_mask = ismember(TRN_edgeList(:, 1), TF_nodes{i}) & ismember(TRN_edgeList(:, 2), TG_nodes{i});
        subTRN_EdgeList{i} = TRN_edgeList(row_mask, :);

        all_nodes = union(subTRN_EdgeList{i}(:, 1), subTRN_EdgeList{i}(:, 2));
        
        
        [~, ii] = ismember(subTRN_EdgeList{i}(:, 1), all_nodes);
        [~, jj] = ismember(subTRN_EdgeList{i}(:, 2), all_nodes);
        Z = unique([ii, jj], 'rows');        
        subTRN_EdgeList{i} = [all_nodes(Z(:, 1)), all_nodes(Z(:, 2))];
        
        
        
        subTRN_NodeAnnotations{i} = [all_nodes, repmat({'TG'}, numel(all_nodes), 1), num2cell(zeros(numel(all_nodes), 1)), num2cell(zeros(numel(all_nodes), 1))];

        subTRN_NodeAnnotations{i}(ismember(all_nodes, subTRN_EdgeList{i}(:, 1)), 2) = {'TF'}; 




        [~, perm] = ismember(subTRN_NodeAnnotations{i}(ismember(all_nodes, subTRN_EdgeList{i}(:, 1)), 1), TF_table{i}(:, 1));
        subTRN_NodeAnnotations{i}(ismember(all_nodes, subTRN_EdgeList{i}(:, 1)), 3) = num2cell(-log10(cell2mat(TF_table{i}(perm, 2))));

        [~, perm] = ismember(subTRN_NodeAnnotations{i}(ismember(all_nodes, subTRN_EdgeList{i}(:, 2)), 1), gene_names);
        subTRN_NodeAnnotations{i}(ismember(all_nodes, subTRN_EdgeList{i}(:, 2)), 4) = num2cell(residual_expression(perm, i));
        v = cell2mat(subTRN_NodeAnnotations{i}(:, 3));
        v = v ./ max(v);
        u = cell2mat(subTRN_NodeAnnotations{i}(:, 4));
        u = u ./ max(u);
        subTRN_NodeAnnotations{i}(:, 5) = num2cell(u+v); % node size
    end    
end

