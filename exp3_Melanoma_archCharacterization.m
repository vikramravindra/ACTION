clear 
addpath(genpath('code'));


load 'results/ACTION/ACTION_full_results.mat'

path.output = fullfile('output', 'archCharacterization');
if(~exist(path.output, 'dir'))
    system(sprintf('mkdir %s -p', path.output));
end

archs = archetypes{2};
arch_mat = zscore(archs, 0, 2);

S = S{2};
[~, labels] = max(S);    
UL = unique(labels);

arch_annotations = arch_annotations_full{2};
dlmcell(fullfile(path.output, 'arch_annotations.txt'), arch_annotations);



%% Load expression datasets
ds_path = fullfile('input', 'datasets', 'Melanoma');
tic; [expression, sample_names, gene_names] = my_tblread(fullfile(ds_path, 'expression.txt')); toc
unexpressed_mask = sum(expression, 2) == 0;
expression(unexpressed_mask, :) = [];
gene_names(unexpressed_mask) = [];
Z = zscore(expression, 0, 2);


%%
sample_annotations = my_dlmread(fullfile(ds_path, 'sample_annotations.txt'));
UP = unique(sample_annotations(:, 2));
tumor_counts = cellfun(@(patient) nnz(strcmp(sample_annotations(:, 2), patient) & strcmp(sample_annotations(:, 4), 'malignant')), UP)

[~, perm] = sort(tumor_counts, 'descend');

tumor_count_report = [UP(perm), num2cell(tumor_counts(perm))]


%% Preprocess markers from the original paper
CellTypes = {'Melanoma', 'T-cells', 'B-cells', 'Macrophages', 'Endothelial cells', 'CAFs'};
Melanoma_Markers = {{'MIA', 'TYR', 'SLC45A2', 'CDH19', 'PMEL', 'SLC24A5', 'MAGEA6', 'GJB1', 'PLP1', 'PRAME', 'CAPN3', 'ERBB3', 'GPM6B', 'S100B', 'FXYD3', 'PAX3', 'S100A1', 'MLANA', 'SLC26A2', 'GPR143', 'CSPG4', 'SOX10', 'MLPH', 'LOXL4', 'PLEKHB1', 'RAB38', 'QPCT', 'BIRC7', 'MFI2', 'LINC00473', 'SEMA3B', 'SERPINA3', 'PIR', 'MITF', 'ST6GALNAC2', 'ROPN1B', 'CDH1', 'ABCB5', 'QDPR', 'SERPINE2', 'ATP1A1', 'ST3GAL4', 'CDK2', 'ACSL3', 'NT5DC3', 'IGSF8', 'MBP'},
{'CD2', 'CD3D', 'CD3E', 'CD3G', 'CD8A', 'SIRPG', 'TIGIT', 'GZMK', 'ITK', 'SH2D1A', 'CD247', 'PRF1', 'NKG7', 'IL2RB', 'SH2D2A', 'KLRK1', 'ZAP70', 'CD7', 'CST7', 'LAT', 'PYHIN1', 'SLA2', 'STAT4', 'CD6', 'CCL5', 'CD96', 'TC2N', 'FYN', 'LCK', 'TCF7', 'TOX', 'IL32', 'SPOCK2', 'SKAP1', 'CD28', 'CBLB', 'APOBEC3G', 'PRDM1'},
{'CD19', 'CD79A', 'CD79B', 'BLK', 'MS4A1', 'BANK1', 'IGLL3P', 'FCRL1', 'PAX5', 'CLEC17A', 'CD22', 'BCL11A', 'VPREB3', 'HLA-DOB', 'STAP1', 'FAM129C', 'TLR10', 'RALGPS2', 'AFF3', 'POU2AF1', 'CXCR5', 'PLCG2', 'HVCN1', 'CCR6', 'P2RX5', 'BLNK', 'KIAA0226L', 'POU2F2', 'IRF8', 'FCRLA', 'CD37'}, 
{'CD163', 'CD14', 'CSF1R', 'C1QC', 'VSIG4', 'C1QA', 'FCER1G', 'F13A1', 'TYROBP', 'MSR1', 'C1QB', 'MS4A4A', 'FPR1', 'S100A9', 'IGSF6', 'LILRB4', 'FPR3', 'SIGLEC1', 'LILRA1', 'LYZ', 'HK3', 'SLC11A1', 'CSF3R', 'CD300E', 'PILRA', 'FCGR3A', 'AIF1', 'SIGLEC9', 'FCGR1C', 'OLR1', 'TLR2', 'LILRB2', 'C5AR1', 'FCGR1A', 'MS4A6A', 'C3AR1', 'HCK', 'IL4I1', 'LST1', 'LILRA5', 'CSTA', 'IFI30', 'CD68', 'TBXAS1', 'FCGR1B', 'LILRA6', 'CXCL16', 'NCF2', 'RAB20', 'MS4A7', 'NLRP3', 'LRRC25', 'ADAP2', 'SPP1', 'CCR1', 'TNFSF13', 'RASSF4', 'SERPINA1', 'MAFB', 'IL18', 'FGL2', 'SIRPB1', 'CLEC4A', 'MNDA', 'FCGR2A', 'CLEC7A', 'SLAMF8', 'SLC7A7', 'ITGAX', 'BCL2A1', 'PLAUR', 'SLCO2B1', 'PLBD1', 'APOC1', 'RNF144B', 'SLC31A2', 'PTAFR', 'NINJ1', 'ITGAM', 'CPVL', 'PLIN2', 'C1orf162', 'FTL', 'LIPA', 'CD86', 'GLUL', 'FGR', 'GK', 'TYMP', 'GPX1', 'NPL', 'ACSL1'},
{'PECAM1', 'VWF', 'CDH5', 'CLDN5', 'PLVAP', 'ECSCR', 'SLCO2A1', 'CCL14', 'MMRN1', 'MYCT1', 'KDR', 'TM4SF18', 'TIE1', 'ERG', 'FABP4', 'SDPR', 'HYAL2', 'FLT4', 'EGFL7', 'ESAM', 'CXorf36', 'TEK', 'TSPAN18', 'EMCN', 'MMRN2', 'ELTD1', 'PDE2A', 'NOS3', 'ROBO4', 'APOLD1', 'PTPRB', 'RHOJ', 'RAMP2', 'GPR116', 'F2RL3', 'JUP', 'CCBP2', 'GPR146', 'RGS16', 'TSPAN7', 'RAMP3', 'PLA2G4C', 'TGM2', 'LDB2', 'PRCP', 'ID1', 'SMAD1', 'AFAP1L1', 'ELK3', 'ANGPT2', 'LYVE1', 'ARHGAP29', 'IL3RA', 'ADCY4', 'TFPI', 'TNFAIP1', 'SYT15', 'DYSF', 'PODXL', 'SEMA3A', 'DOCK9', 'F8', 'NPDC1', 'TSPAN15', 'CD34', 'THBD', 'ITGB4', 'RASA4', 'COL4A1', 'ECE1', 'GFOD2', 'EFNA1', 'PVRL2', 'GNG11', 'HERC2P2', 'MALL', 'HERC2P9', 'PPM1F', 'PKP4', 'LIMS3', 'CD9', 'RAI14', 'ZNF521', 'RGL2', 'HSPG2', 'TGFBR2', 'RBP1', 'FXYD6', 'MATN2', 'S1PR1', 'PIEZO1', 'PDGFA', 'ADAM15', 'HAPLN3', 'APP'},
{'FAP', 'THY1', 'DCN', 'COL1A1', 'COL1A2', 'COL6A1', 'COL6A2', 'COL6A3', 'CXCL14', 'LUM', 'COL3A1', 'DPT', 'ISLR', 'PODN', 'CD248', 'FGF7', 'MXRA8', 'PDGFRL', 'COL14A1', 'MFAP5', 'MEG3', 'SULF1', 'AOX1', 'SVEP1', 'LPAR1', 'PDGFRB', 'TAGLN', 'IGFBP6', 'FBLN1', 'CA12', 'SPOCK1', 'TPM2', 'THBS2', 'FBLN5', 'TMEM119', 'ADAM33', 'PRRX1', 'PCOLCE', 'IGF2', 'GFPT2', 'PDGFRA', 'CRISPLD2', 'CPE', 'F3', 'MFAP4', 'C1S', 'PTGIS', 'LOX', 'CYP1B1', 'CLDN11', 'SERPINF1', 'OLFML3', 'COL5A2', 'ACTA2', 'MSC', 'VASN', 'ABI3BP', 'C1R', 'ANTXR1', 'MGST1', 'C3', 'PALLD', 'FBN1', 'CPXM1', 'CYBRD1', 'IGFBP5', 'PRELP', 'PAPSS2', 'MMP2', 'CKAP4', 'CCDC80', 'ADAMTS2', 'TPM1', 'PCSK5', 'ELN', 'CXCL12', 'OLFML2B', 'PLAC9', 'RCN3', 'LTBP2', 'NID2', 'SCARA3', 'AMOTL2', 'TPST1', 'MIR100HG', 'CTGF', 'RARRES2', 'FHL2'}};

markers_idx = cell(numel(Melanoma_Markers), 1);
for j = 1:numel(Melanoma_Markers)
    [mask, markers_idx{j}] = ismember(Melanoma_Markers{j}, gene_names);
    markers_idx{j}(~mask) = [];
end




col_labels = arrayfun(@(x) sprintf('A%d', x), 1:numel(arch_annotations), 'UniformOutput', false);
CC = cbrewer('seq', 'BuGn', 128);

%% Node annotations
    node_names = arrayfun(@(c) sprintf('A%d', c), 1:numel(UL), 'UniformOutput', false)';
    Marker_expressions = cell(numel(UL)+1, numel(Melanoma_Markers)+3);
    Marker_expressions(2:end, 1) = node_names;
    Marker_expressions(1, :) = [{'ID', 'Size', 'annotations'}, CellTypes];
    Marker_expressions(2:end, 3) = arch_annotations;
    
    
    for i = 1:numel(UL)
        samples = find(labels == UL(i));
        Marker_expressions{i+1, 2} = numel(samples);
        for j = 1:numel(Melanoma_Markers)
            subExpr = arch_mat(markers_idx{j}, i);            
            Marker_expressions{i+1, j+3} = trimmean((nonzeros(subExpr)), 10);
        end
    end

    
    dlmcell(fullfile(path.output, 'Melanoma_MarkerExpressions.txt'), Marker_expressions);

%%
MarkerEnrichments = cell2mat(Marker_expressions(2:end, 4:end));

heatmap(MarkerEnrichments',col_labels, CellTypes, [], 'colormap', CC);
set(gcf, 'Position', [0 0 100*numel(col_labels) 100*numel(CellTypes)])
pause(1);

set(gca, 'XTick', 1:numel(col_labels));
set(gca, 'XTickLabel', col_labels);

set(gca, 'YTick', 1:numel(CellTypes));
set(gca, 'YTickLabel', CellTypes);
 
colorbar();

set(gca,'FontSize', 14, 'FontWeight','bold'); 
export_fig(fullfile(path.output, 'marker_enrichment.eps'), '-eps', '-transparent', '-painters');
close
%% Immune Subtyping    
%     LM22_names = {'B cells naive', 'B cells memory', 'Plasma cells', 'T cells CD8', 'T cells CD4 naive', 'T cells CD4 memory activated', 'NK cells activated', 'Monocytes', 'Dendritic cells activated', 'Eosinophils', 'Neutrophils'};


    
    [LM22_marker_expressions, LM22_immune_cellTypes, LM22_marker_genes] = my_tblread(fullfile('input', 'markers', 'LM22.txt'));
    [mask, LM22_marker_idx] = ismember(LM22_marker_genes, gene_names);
    LM22_marker_idx(~mask) = [];
    LM22_marker_expressions(~mask, :) = [];    

    PP = normalize(LM22_marker_expressions, 'pnorm', 1, 'dim', 2);
    marker_mask = (PP > (1/size(PP, 2)));     
    
    LM22_enrichment = zeros(numel(LM22_immune_cellTypes), numel(UL));
    for i = 1:numel(UL)
        samples = find(labels == UL(i));
        for j = 1:numel(LM22_immune_cellTypes)
            rows = LM22_marker_idx(marker_mask(:, j));

            subExpr = Z(rows, samples);  
%             subExpr = arch_mat(rows, i);
            LM22_enrichment(j, i) = trimmean((nonzeros(subExpr)), 1);
        end
    end    
    

row_labels = LM22_immune_cellTypes;

LM22_results = cell(size(LM22_enrichment)+1);
LM22_results(2:end, 2:end) = num2cell(LM22_enrichment);
LM22_results(2:end, 1) = row_labels;
LM22_results(1, 2:end) = col_labels;
dlmcell(fullfile(path.output, 'LM22_table.txt'), LM22_results);



beneign_archs = ~cellfun(@(arch) nnz(strcmp('[malignant]', strsplit(arch, '-')) | strncmp('[CAF', strsplit(arch, '-'), 4) | strncmp('[Endo', strsplit(arch, '-'), 5)), arch_annotations);

heatmap(zscore(LM22_enrichment(:, beneign_archs)), col_labels(beneign_archs), row_labels, [], 'colormap', CC);
set(gcf, 'Position', [0 0 50*numel(col_labels)+70 50*numel(row_labels)])
pause(1);

set(gca, 'XTick', 1:nnz(beneign_archs));
set(gca, 'XTickLabel', col_labels(beneign_archs));

set(gca, 'YTick', 1:numel(row_labels));
set(gca, 'YTickLabel', row_labels);
 
colorbar();

set(gca,'FontSize', 14, 'FontWeight','bold'); 
export_fig(fullfile(path.output, 'immuno_subtyping.eps'), '-eps', '-transparent', '-painters');
close
    
%%
[~, LM22_perm] = sort(LM22_enrichment, 'descend');
LM22_Dominant_Celltypes= arrayfun(@(col) LM22_immune_cellTypes(LM22_perm(:, col))', 1:size(LM22_enrichment, 2), 'UniformOutput', false);


CT_cell = [LM22_Dominant_Celltypes{:}];
CT_cell = [arrayfun(@(id) sprintf('A%d', id), 1:size(archs, 2), 'Uniformoutput', false); CT_cell];
CT_cell = [[{''}; arrayfun(@(id) sprintf('Rank %d', id), 1:numel(LM22_immune_cellTypes), 'Uniformoutput', false)'], CT_cell];

dlmcell(fullfile(path.output, 'Dominant_CellTypes.txt'), CT_cell);

