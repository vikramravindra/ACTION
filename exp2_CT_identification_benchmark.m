clear 

addpath(genpath('code'));


path.results = 'results';
path.SC = fullfile('input', 'datasets');
path.output = fullfile('output', 'celltype_identification_benchmark');

ds_list = dir(sprintf('%s/*', path.SC)); 
ds_list(~[ds_list.isdir]) = [];
ds_list(1:2) = [];
ds_no = numel(ds_list);   

dataset_names = {ds_list(:).name};


if( ~exist(path.output, 'dir') )
    system(sprintf('mkdir %s -p', path.output));
end

methods = {'ACTION', 'ParTI', 'SCUBA', 'SNNCliq', 'TSCAN'};
method_count = numel(methods);


%%

ARI = zeros(numel(dataset_names), method_count);
NMI = zeros(numel(dataset_names), method_count);

for ds_id = 1:numel(dataset_names)
    ds_path = fullfile(path.SC, dataset_names{ds_id});
    fprintf('Reading SC dataset %s\n', dataset_names{ds_id});

    sample_annotations = my_dlmread(fullfile(ds_path, 'sample_annotations.txt'));
    Labels = sample_annotations(:, end);
       
    UL = unique(Labels);
    [~, l] = ismember(Labels, UL);
    CT_no(ds_id, 1) = numel(unique(l));
    
    display(UL);  

    for i = 1:method_count
        
        curr_method = methods{i};

        if(~exist(fullfile(path.results, curr_method, sprintf('%s.mat', dataset_names{ds_id})), 'file'))
            continue;
        end        
        clear labels
        
        load(fullfile(path.results, curr_method, sprintf('%s.mat', dataset_names{ds_id})), 'labels');

%         labels(labels == 0) = 1;
        NMI(ds_id, i) = nmi(l, labels);
        ARI(ds_id, i) = adjustedrand(l, labels);
        CT_count(ds_id, i) = numel(unique(labels));
    end
end


%%
C = cbrewer('qual', 'Set1', method_count);    
Edge_color = C;
Face_color= 0.7*C + 0.3*ones(size(C));


Data = NMI;
fig = figure;
set(gcf,'units','normalized','outerposition',[0 0 1 1], 'PaperPositionMode', 'auto');
handle = bar(Data,'grouped');
for i = 1:numel(handle)
    set(handle(i),'FaceColor',Face_color(i,:));
    set(handle(i),'EdgeColor',Edge_color(i,:));
end
set(gca, 'XTickLabel', dataset_names);
set(gca,'FontSize', 34, 'FontWeight','bold');     

ylabel('Normalized Mutual Information');
legend(methods, 'Interpreter', 'None', 'FontSize', 24,'fontWeight','bold', 'Location', 'NorthWest');
hold on


print(fig, fullfile(path.output, 'NMI.eps'), '-depsc2','-r300');
% close;  


Data = ARI;
fig = figure;
set(gcf,'units','normalized','outerposition',[0 0 1 1], 'PaperPositionMode', 'auto');
handle = bar(Data,'grouped');
for i = 1:numel(handle)
    set(handle(i),'FaceColor',Face_color(i,:));
    set(handle(i),'EdgeColor',Edge_color(i,:));
end
set(gca, 'XTickLabel', dataset_names);
set(gca,'FontSize', 34, 'FontWeight','bold');     


ylabel('Adjusted Rand Index');
legend(methods, 'Interpreter', 'None', 'FontSize', 24,'fontWeight','bold', 'Location', 'NorthWest');

hold on

print(fig, fullfile(path.output, 'ARI.eps'), '-depsc2','-r300');
% close;  
