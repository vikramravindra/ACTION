%% Write graph to a file in sif or dot format
% <adjmat> contains the adjacency matrix encoding the graph to be written
% <filename> is the file to be written to
% [filetype] is an optional input to determine the type of file to be
% written (options are 'sif'(default) or 'dot')
% sif: Cytoscape simple interaction format
% dot: GraphViz dot format
function [] = write_graph_from_adj(adjmat,filename,varargin)
    filetype='sif';
    if numel(varargin)>0
        switch varargin(1)
            case 'sif'
                filetype='sif';
            case 'dot'
                filetype='dot';
            otherwise
                sprintf('Invalid file format specifier: %s\n Assuming sif file.',varargin(1));
        end
    end
    file = fopen(filename,'w');
    if strcmp(filetype,'dot')
        fprintf(file,'digraph G {\n');
        for i=1:size(adjmat,1)
            for j=1:size(adjmat,2)
                if adjmat(i,j)
                    fprintf(file,'%i\t->\t%i\n',i,j);
                end
            end
        end
        fprintf(file,'}\n');
    end
    if strcmp(filetype,'sif')
        for i=1:size(adjmat,1)
            for j=1:size(adjmat,2)
                if adjmat(i,j)
                    fprintf(file,'%i\t->\t%i\n',i,j);
                end
            end
        end
    end
    fclose(file);
end