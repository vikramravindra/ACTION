clear 

addpath(genpath('code/'));


path.results = 'results';
path.SC = fullfile('input', 'datasets');
path.output = fullfile('output', 'TRN');

ds_list = dir(sprintf('%s/*', path.SC)); 
ds_list(~[ds_list.isdir]) = [];
ds_list(1:2) = [];
ds_no = numel(ds_list);   

dataset_names = {ds_list(:).name};


if( ~exist(path.output, 'dir') )
    system(sprintf('mkdir %s -p', path.output));
end


load(fullfile(path.results, 'ACTION', 'ACTION_full_results.mat'), 'adjusted_expressions', 'C', 'S', 'archetypes', 'arch_annotations', 'arch_annotations_full');


%% Read expression matrix -- wasteful since we only need gene_names here!
    Melanoma_ds_id = 2;
    ds_path = fullfile(path.SC, dataset_names{Melanoma_ds_id});

    tic; [expression, sample_names, gene_names] = my_tblread(fullfile(ds_path, 'expression.txt')); toc
    mask = (sum(expression, 2) == 0);
    expression(mask, :) = [];
    gene_names(mask) = [];

    sample_annotations = my_dlmread(fullfile(ds_path, 'sample_annotations.txt'));

    n = size(expression, 2);

    Labels = sample_annotations(:, end);
    UL = unique(Labels);
    [~, l] = ismember(Labels, UL);    
%% Construct background TRN

% RegNet = my_dlmread('input/TRN/RegNetwork.tsv', '\t');
% RegNet(strncmp('MI', RegNet(:, 2), 2), :) = []; % remove miRNA→TF and miRNA→gene interactions
% RegNet(strncmp('MI', RegNet(:, 4), 2), :) = []; % remove TF→miRNA interactions
% TRN_edgeList = RegNet(:, [1, 3]);


RegNet = my_dlmread('input/TRN/trrust_rawdata.txt', '\t');
negative_edges = strncmp(RegNet(:, 3), 'Repr', 4);
RegNet(negative_edges, :)= [];
TRN_edgeList = RegNet(:, [1, 2]);


selfloop_mask = arrayfun(@(row) strcmp(TRN_edgeList(row, 1), TRN_edgeList(row, 2)), 1:size(TRN_edgeList, 1));
TRN_edgeList(selfloop_mask, :) = [];

%% Construct archetype-specific TRNs
[residual_expression, sorted_genes, TFs, mHG_pvals, mHG_pvals_corrected, selected_genes, TF_table, subTRN_EdgeList, subTRN_NodeAnnotations] = dissectTRN( archetypes{Melanoma_ds_id}, gene_names, TRN_edgeList, 1:size(archetypes{Melanoma_ds_id}, 2), 'pval_threshold', 0.05);

save(fullfile(path.output, 'Melanoma_TRN.mat'), 'residual_expression', 'sorted_genes', 'TFs', 'mHG_pvals', 'mHG_pvals_corrected', 'selected_genes', 'TF_table', 'subTRN_EdgeList', 'subTRN_NodeAnnotations');

%% Export sorted genes, as well as their corresponding regulatory archetype/cell type-specific TFs to file
    
    dlmcell(fullfile(path.output, 'sorted_genes.tsv'), [arrayfun(@(id) sprintf('Celltype%d', id), 1:size(archetypes{Melanoma_ds_id}, 2), 'UniformOutput', false); sorted_genes]);
    
    for i = 1:size(archetypes{Melanoma_ds_id}, 2)
        dlmcell(fullfile(path.output, sprintf('Celltype%d_TFs.tsv', i)), [{'TF', 'pvalue'};TF_table{i}]);
    end


% %% p-value of patients with NRAS Q61L subclass that are enriched in archtype 4 (numbers may change if you run with different parameters!!)
% 
% NRAS_pval = hygecdf(4-1, 19, 5,4, 'upper') % Stats from Aviv's paper (19 patients, 4 of them ar NRAS Q61L mutants

%% Evaluate enrichment of markers for "proliferative" subclass. "Supposedly" this should be a good subclass, with these genes being good, but Cox says otherwise!
% 17 common genes out of 54 genes in the network (pval = 1.15e-12) [770 total]
MITF_edges = subTRN_EdgeList{4};
MITF_nodes = subTRN_NodeAnnotations{4};
MTF_TFs = TF_table{4};
total_genes = union(MITF_edges(:, 1), MITF_edges(:, 2));


proliferative_genes = {'MAD2L1BP', 'TXN2', 'ADAM10', 'NAGA', 'ACO2', 'MT-CO1', 'TRUB1', 'HMGCR', 'CCNC', 'TANC1', 'TRAK2', 'RPLP1', 'BCL2L13', 'SGSM3', 'KCP', 'POPDC3', 'TTL', 'PTCD3', 'ATE1', 'LMBRD1', 'RSAD1', 'GGCT', 'WDR11', 'RHBDD3', 'ATP6V1E1', 'ABCC5', 'TMEM33', 'NDUFAF4', 'TCTN3', 'DDI2', 'TAB3', 'CD320', 'MCOLN1', 'PTEN', 'CABIN1', 'GNPTG', 'FBXO7', 'CLN3', 'USP6NL', 'RP11-159F24.1', 'DUSP22', 'ESRRA', 'CLUH', 'PPA1', 'LOH12CR1', 'FASN', 'MINPP1', 'NPM3', 'GGA2', 'ZNF74', 'MOAP1', 'TLE1', 'ZNF706', 'FOXO4', 'WDFY1', 'DHTKD1', 'CNNM3', 'MMS19', 'PSEN2', 'GBAS', 'PCMT1', 'FARP1', 'FAM89A', 'MSI2', 'RPS12', 'TSPAN3', 'ANKRD46', 'DMXL2', 'HS1BP3', 'STX3', 'SREBF2', 'TMEM64', 'ABRACL', 'SSH2', 'SNHG16', 'TOM1', 'PPP2R5A', 'DTNB', 'MYLIP', 'OSBPL1A', 'MED15', 'FAM217B', 'TRAF4', 'DGCR6L', 'SEC11C', 'ZNF138', 'PEBP1', 'THEM6', 'AGPAT6', 'NECAB3', 'TYRO3', 'GYPC', 'SRP14-AS1', 'NHLRC3', 'ZFYVE27', 'GSTM4', 'BIN3', 'MYO10', 'CTNNB1', 'ATM', 'SLC2A11', 'TP53', 'SNX30', 'DHX33', 'FNBP1L', 'PPP1R3D', 'PAX3', 'MPC1', 'FAM63B', 'MOB1B', 'HPS1', 'PRMT7', 'LINC00263', 'EIF4EBP2', 'ACSL3', 'GM2A', 'SFXN4', 'CELSR2', 'ATP1A1', 'TYSND1', 'R3HCC1L', 'RASSF3', 'AMDHD2', 'IDI1', 'RHOQ', 'PACSIN2', 'TFAP2A', 'PIP4K2A', 'CHD7', 'FZD3', 'ALDH4A1', 'TMEM251', 'EPM2A', 'ATP6V0E2', 'FAM101B', 'TPT1-AS1', 'FAM161A', 'PLD1', 'MYO5A', 'ARHGAP19', 'SDCBP', 'INPP4B', 'MEGF9', 'RP11-1055B8.7', 'SLC19A2', 'WIPI1', 'KLHL22', 'C1orf85', 'FUOM', 'SLC6A8', 'DAAM1', 'ERMP1', 'LPCAT1', 'SMPD2', 'SLC39A11', 'POLR3G', 'OSER1-AS1', 'OAT', 'TULP4', 'ATP6V0A1', 'CYTH3', 'ZBTB10', 'TPD52', 'FAM210B', 'ARVCF', 'UBL3', 'LDLRAD3', 'SLC23A2', 'BACE2', 'AC083843.1', 'LONRF1', 'TTLL1', 'IVNS1ABP', 'ARRB1', 'STAT5A', 'RAB32', 'NANOS1', 'LIPA', 'BCL2', 'USP54', 'STXBP1', 'VAT1', 'SASH1', 'TMEM98', 'TRAPPC6A', 'ATP6V1C1', 'GRAMD4', 'FARP2', 'SERINC5', 'COMTD1', 'CHCHD6', 'IRS2', 'MYC', 'SLC25A23', 'ST3GAL5', 'MAP6D1', 'MFGE8', 'PPARGC1B', 'GLUD2', 'PI4KAP1', 'PARD6G', 'HSF4', 'FBXL7', 'GPR153', 'PFKFB2', 'PSAP', 'RP11-334C17.5', 'DIP2C', 'RAB3A', 'THNSL1', 'AVPI1', 'GNPTAB', 'LGALS3', 'TANGO2', 'ATP7B', 'GSTO1', 'SOX13', 'CLCN7', 'KREMEN1', 'SETDB2', 'PPFIBP2', 'ASB13', 'CERS4', 'CHCHD10', 'SCARB1', 'ARNT2', 'TBC1D14', 'EN2', 'QDPR', 'AK9', 'IGSF8', 'RGS20', 'KAZN', 'DSTYK', 'PODXL2', 'DFNB31', 'ZBED3-AS1', 'KIAA0930', 'FNIP2', 'GK', 'TPCN2', 'ANKRD44', 'KREMEN2', 'PLEKHG3', 'OSTM1', 'ATP11A', 'FAXDC2', 'NHSL1', 'MFAP3L', 'TCFL5', 'ITPK1', 'CHKA', 'TRIB1', 'TRIB2', 'PIK3C2B', 'MYH10', 'SLC25A16', 'CTD-3065J16.9', 'MAP3K1', 'GPNMB', 'IGSF3', 'HSPA12A', 'ADAM11', 'MARCKSL1', 'TATDN2P2', 'SOX5', 'C1orf198', 'ARSG', 'FMN1', 'RP11-356I2.4', 'FAHD2CP', 'ACSL1', 'GS1-358P8.4', 'GPRC5B', 'LYST', 'SCD', 'SLC18B1', 'LINC00473', 'CCDC78', 'SLC5A4', 'TPPP', 'LINC00920', 'SMPDL3A', 'FRMD3', 'SIRPA', 'TMEM170B', 'NXN', 'HPS4', 'TDRD3', 'GPR56', 'CTA-217C2.1', 'CRYL1', 'MICAL1', 'KIF1A', 'CLCN5', 'RHPN1', 'LEF1', 'RAB27A', 'RNF144B', 'LA16c-390E6.4', 'ADARB1', 'SEMA6D', 'CHST11', 'DNAJC12', 'REEP6', 'CDK18', 'GPR19', 'PDGFD', 'FN3K', 'RTN4R', 'SCAMP5', 'C10orf11', 'ZFYVE16', 'HMG20B', 'FAM149A', 'SHC4', 'FAM222A', 'SLC1A4', 'STX7', 'CUBN', 'HAGHL', 'TSPAN33', 'NEDD9', 'PLA2G6', 'CGN', 'FAM69B', 'ACACB', 'TRIM2', 'RP11-1055B8.4', 'PRKD3', 'RGS10', 'C21orf91', 'HSD17B14', 'TBC1D16', 'DGCR6', 'FCRLA', 'WDR91', 'GDPD5', 'SLC22A23', 'RP11-3L8.3', 'PDZRN3', 'NR4A1', 'PROS1', 'SS18L1', 'AP1S2', 'TSPAN10', 'BHLHE41', 'CCDC171', 'IL17D', 'PDK4', 'C2orf88', 'PRR5', 'CDK2', 'CTSH', 'MREG', 'FAM124A', 'RP11-558F24.4', 'ASAH1', 'HEY2', 'SDC3', 'RAP1GAP', 'CAPG', 'RAB6B', 'SORT1', 'MXI1', 'PCSK9', 'ST3GAL4', 'APOC1', 'ADAMTS17', 'EPHA5', 'HEY1', 'SLC12A7', 'SH3TC1', 'PAG1', 'RAB38', 'RP11-1055B8.2', 'SCIN', 'LINC00340', 'GAB2', 'CES3', 'DUSP15', 'CYP27A1', 'ATRNL1', 'FAXC', 'HES6', 'KLF15', 'DGCR5', 'AC005786.5', 'INPP5F', 'SGK1', 'NPAS1', 'CYGB', 'PIR', 'CERS1', 'AATK', 'SESN3', 'CDK5R1', 'SCUBE2', 'SLC17A9', 'SLC16A10', 'TMTC2', 'KIAA1598', 'RASSF2', 'MAST1', 'RP11-80F22.9', 'TBC1D7', 'PRKCH', 'PDE3B', 'C19orf71', 'SLAIN1', 'GJA3', 'FAM53B', 'RAB11FIP4', 'BAI1', 'C11orf96', 'SLC27A3', 'FGF13', 'BRSK2', 'EGLN3', 'GNAL', 'CABLES1', 'GPR137B', 'CXADR', 'SHC2', 'ST6GALNAC1', 'FBXL16', 'Z83851.1', 'ASRGL1', 'TNFRSF19', 'CEACAM1', 'SORL1', 'ANKRD6', 'ISG20', 'RIMS4', 'MYOM2', 'LAD1', 'ADRBK2', 'LZTS1', 'RNF125', 'TMEM255A', 'TRPM8', 'FAM20A', 'LONRF3', 'PPM1H', 'SPTBN2', 'TEX41', 'ST6GAL1', 'POU3F3', 'ADAM23', 'ANO4', 'MFSD12', 'RP11-390P2.4', 'ST6GALNAC2', 'RP11-137H2.6', 'OLFM2', 'TMCC2', 'GREB1', 'TTC39A', 'FAM213A', 'KCNS1', 'TNFRSF14', 'STXBP6', 'ITGA7', 'ALDH1A1', 'ZNF704', 'BAMBI', 'PGBD5', 'PRKCZ', 'IL6R', 'PLCL1', 'EGR3', 'ITPKB', 'NAT16', 'LRRC4', 'STOX2', 'OGDHL', 'PIK3AP1', 'PNLIPRP3', 'CNTN3', 'BAAT', 'COL25A1', 'CELF2', 'RASIP1', 'TMEM229B', 'PLEKHG1', 'PKNOX2', 'KRTAP19-1', 'SLC7A4', 'SLC24A4', 'ASB4', 'ST3GAL6-AS1', 'EFR3B', 'NKAIN1', 'CASKIN1', 'LDLRAD4', 'WNK2', 'TKTL1', 'RAB17', 'KNDC1', 'TESK2', 'CHN2', 'SLC7A8', 'FGD4', 'CRTAC1', 'PPARGC1A', 'RP11-527H14.2', 'FAM134B', 'RASEF', 'ACAN', 'CHST6', 'TFCP2L1', 'HMCN1', 'RP3-395M20.8', 'PNMAL1', 'RP3-527G5.1', 'PIP5K1B', 'IL12RB2', 'TENM1', 'RPS6KA2', 'CECR2', 'VGF', 'BCAN', 'ADCY1', 'RAB3C', 'CLDN14', 'RP11-481A20.11', 'MERTK', 'LINC00937', 'LINC00504', 'CCL18', 'RXRG', 'PHACTR1', 'CARD14', 'QPCT', 'GRASP', 'DLL3', 'PKLR', 'LRGUK', 'C1orf51', 'TEX15', 'B4GALNT3', 'KIAA1211', 'PELI2', 'POU3F2', 'PRKCB', 'HCG20', 'RP11-98L5.2', 'DISC1FP1', 'RP11-317M11.1', 'KBTBD11', 'LAMC3', 'RP11-557H15.4', 'MPZ', 'AC011294.3', 'RENBP', 'PRUNE2', 'LAMA1', 'B3GAT1', 'TINCR', 'NDN', 'TTYH2', 'CTB-151G24.1', 'TC2N', 'SLC16A6', 'RP11-143A12.3', 'LINC00518', 'MFI2', 'PRODH', 'GOLGA7B', 'POU3F1', 'CDH3', 'GPM6A', 'NR4A3', 'MAPT', 'LARGE', 'LRP2', 'LPL', 'TMPRSS5', 'RLBP1', 'GNG7', 'ACP5', 'RP11-93B14.5', 'FAM167B', 'EDNRB', 'TMC6', 'LINC00426', 'CA8', 'MYO16', 'ST3GAL6', 'CITED1', 'RASGEF1A', 'AC009784.3', 'DAPK1', 'KIT', 'GYG2', 'PLEKHB1', 'AP000479.1', 'TMPRSS13', 'NUP210', 'ALDH1A2', 'ZNF536', 'FAM19A5', 'RGS1', 'NKX2-5', 'OPLAH', 'TSPAN7', 'KDR', 'PLEKHH1', 'SBK1', 'FAM155B', 'ITGA9', 'BEST1', 'KCNAB2', 'RP11-2E17.1', 'MAGEB2', 'RP13-735L24.1', 'TUBB8P7', 'CPVL', 'DENND1C', 'MOB3B', 'ST8SIA6', 'WIPF3', 'PRSS33', 'CNTN1', 'CCDC64', 'APOD', 'NRG3', 'RP3-395M20.7', 'AC009499.1', 'PPP1R14C', 'CAPN3', 'SOX6', 'MBP', 'ISM1', 'MITF', 'CACNA1H', 'RP4-718J7.4', 'MGAT4A', 'CPB2-AS1', 'LRRC4B', 'EXTL1', 'SOX8', 'CRYAB', 'ROPN1B', 'SYT3', 'SGCA', 'NAT8L', 'RP11-509E16.1', 'PMP2', 'NMRK2', 'AC004988.1', 'GSTT1', 'LINC00589', 'RP11-189B4.6', 'TYRP1', 'FAM189A2', 'ALDH3B2', 'RP11-161M6.2', 'RAB33A', 'ZDHHC11B', 'ROBO2', 'SLC35F1', 'RRAGD', 'LINGO1', 'RP11-290F20.3', 'HSPB8', 'RP11-669N7.2', 'RP3-332B22.1', 'CACNA1D', 'SAMD5', 'GFPT2', 'ITGAX', 'CPN1', 'RP11-726G1.1', 'MSI1', 'TNRC6C-AS1', 'LGI3', 'MLIP', 'TUBB4A', 'COBL', 'LHFPL3-AS1', 'LINC00488', 'RP11-557H15.2', 'GALNT3', 'S100B', 'AC002511.1', 'SORBS1', 'IGF1', 'ADCY2', 'DLGAP1', 'SFTPC', 'HILS1', 'PLXNC1', 'MYH14', 'AC145110.1', 'PLA1A', 'GLB1L2', 'CHL1', 'GAPDHS', 'CTNNA2', 'RP11-599J14.2', 'SHE', 'CTD-2207A17.1', 'RP11-1055B8.3', 'RP11-104E19.1', 'AC096559.1', 'FXYD3', 'HPGD', 'MMP8', 'IRX6', 'RP11-347E10.1', 'GJB1', 'GPR143', 'PLP1', 'CDH19', 'IL16', 'MPPED2', 'RP11-252C15.1', 'CDH1', 'C10orf90', 'APOE', 'ITIH5', 'OCA2', 'MAF', 'FAM69C', 'RP4-529N6.1', 'SOX10', 'FRMD4B', 'ERBB3', 'SEMA6A', 'MCF2L', 'MAPK4', 'FAM174B', 'SLC38A8', 'MYO1D', 'LINC00520', 'NSG1', 'IGSF11', 'C20orf26', 'HRK', 'PAEP', 'SLC45A2', 'TRIM51', 'ATP10A', 'ROPN1', 'COL9A3', 'SCML4', 'RP11-429E11.2', 'GAS7', 'CA14', 'LINC00698', 'PTPRZ1', 'DCT', 'MID2', 'PMEL', 'SORCS1', 'ABCB5', 'PRDM7', 'GPM6B', 'LCP2', 'ENTHD1', 'MLANA', 'TYR', 'NKAIN4', 'IRF4', 'SGCD', 'SLC24A5', 'TRIM63', 'BIRC7', 'TRPM1'};

common_prolif_genes = MITF_nodes(ismember(MITF_nodes(:, 1), proliferative_genes), 1);
prolif_pval = hygecdf(numel(common_prolif_genes)-1, size(expression, 1), numel(proliferative_genes), numel(total_genes), 'upper');
fprintf('%d common genes out of %d genes in the network (pval = %.2e)\n', numel(common_prolif_genes), numel(total_genes), prolif_pval)

%% Compute shift in positive/negativeness of genes among identified genes (larger p-val:: 2.26-11)
% "A positive Cox coefficient indicates high expression of the gene
% increases the risk of death, while a negative Cox coefficient indicates the opposite."
Onco = readtable(fullfile('input', 'OncoLnc', 'SKCM_mrna.csv'));
background_dist = Onco.CoxCoefficient;

subOnco = Onco(ismember(Onco.UpdatedName, total_genes), :);
TRN_dist = subOnco.CoxCoefficient;


[~, Cox_pval_larger] = kstest2(background_dist, TRN_dist, 'Tail','larger');
[~, Cox_pval_smaller] = kstest2(background_dist, TRN_dist, 'Tail','smaller');

fprintf('Median of background is: %.2f, Median of TRN genes = %.2f (pval of being smaller = %.2e, pval of being larger = %.2e)\n', median(background_dist), median(TRN_dist), Cox_pval_smaller, Cox_pval_larger);

[f1, x1] = ksdensity(background_dist, linspace(min(background_dist), max(background_dist), 200));
[f2, x2] = ksdensity(TRN_dist, linspace(min(TRN_dist), max(TRN_dist), 200));

figure
% set(gcf,'units','normalized','outerposition',[0 0 1 1], 'PaperPositionMode', 'auto');
hold all

plot(x1, f1, 'Color', [0.7, 0.7, 0.7], 'LineWidth', 2);
plot(x2, f2, 'Color', [0 0 0.7], 'LineWidth', 2);  
xlabel('Cox coefficient','FontSize', 12, 'FontWeight','bold');
ylabel('pdf','FontSize', 12, 'FontWeight','bold');

title(sprintf('p-value = %.2e', Cox_pval_larger),'FontSize', 14, 'FontWeight','bold');
set(gca,'FontSize', 11, 'FontWeight','bold');   

legend({'Background', 'TRN genes'},'FontSize', 14, 'FontWeight','bold', 'Location', 'NorthWest');

%% Significant TFs for arch 4&5
total_TFs = union(TF_table{4}(:, 1), TF_table{5}(:, 1));
[~, idx] = ismember(total_TFs, TFs);
pvals_sub = mHG_pvals_corrected([4, 5], idx);
folds = -log10(pvals_sub(2, :) ./ pvals_sub(1, :));
[~, perm] = sort(folds);
geo_mean = sqrt(pvals_sub(2, :) .* pvals_sub(1, :));
Report =[total_TFs(perm), num2cell(pvals_sub(1, perm))', num2cell(pvals_sub(2, perm))', num2cell(folds(perm))', num2cell(-log10(geo_mean(perm))')];

dlmcell(fullfile(path.output, 'MITF_associated_TF_list.txt'), Report);

%% MITF targets for arch 4&5
targets4 = subTRN_EdgeList{4}(strcmp(subTRN_EdgeList{4}(:, 1), 'MITF'), 2);
targets5 = subTRN_EdgeList{5}(strcmp(subTRN_EdgeList{5}(:, 1), 'MITF'), 2);

Common_targets = intersect(targets4, targets5);
A4_targets = setdiff(targets4, targets5);
A5_targets = setdiff(targets5, targets4);

sorted_targets = [A4_targets; Common_targets; A5_targets];
MITF_target_report = [sorted_targets, num2cell(ismember(sorted_targets, targets4)), num2cell(ismember(sorted_targets, targets5))];

dlmcell(fullfile(path.output, 'MITF_targets.txt'), MITF_target_report);

%% Export TRN edge list/node annotation files to be imported in Cytoscape
[residual_expression, sorted_genes, TFs, mHG_pvals, mHG_pvals_corrected, selected_genes, TF_table, subTRN_EdgeList, subTRN_NodeAnnotations] = dissectTRN( archetypes{Melanoma_ds_id}, gene_names, TRN_edgeList, 1:size(archetypes{Melanoma_ds_id}, 2), 'pval_threshold', 1e-3);

save(fullfile(path.output, 'Melanoma_TRN_highSig.mat'), 'residual_expression', 'sorted_genes', 'TFs', 'mHG_pvals', 'mHG_pvals_corrected', 'selected_genes', 'TF_table', 'subTRN_EdgeList', 'subTRN_NodeAnnotations');

for i = 1:size(archetypes{Melanoma_ds_id}, 2)
    if(isempty(subTRN_EdgeList{i}))
        continue;
    end
    
    cnode_table = subTRN_NodeAnnotations{i};
    [mask, idx] = ismember(cnode_table(:, 1), Onco.UpdatedName);
    
    cox_coeff = Onco.CoxCoefficient(idx(mask));
    cox_coeff(Onco.BH_adjustedP_value(idx(mask)) > 0.05) = 0;
    cnode_table(mask, end+1) = num2cell(cox_coeff);
    cnode_table = [{'ID', 'Type', '-log(pval)', 'residual_expression', 'size', 'Cox'}; cnode_table];
    
    dlmcell(fullfile(path.output, sprintf('Arch%d_TRN_EdgeList.txt', i)), subTRN_EdgeList{i});
    dlmcell(fullfile(path.output, sprintf('Arch%d_TRN_NodeAnnotation.txt', i)), cnode_table);    
end

