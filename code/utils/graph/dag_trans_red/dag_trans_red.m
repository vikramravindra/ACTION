%% Compute a transitive reduction of a DAG (directed ACYCLIC graph) 
% The DAG has to be specified in adjacency matrix form. 
% The algorithm: 
% 1. computes a topological ordering of the nodes
% 2. iteratively treats the nodes in their topological order starting with a node that
%    has no predecessors
% 3. finds indirect ancestors of node i and removes edges pointing from them
%    to i
% 4. adds the direct predecessors of node i to its ancestor list
% 
% Since the algorithm starts with the nodes highest in the topological
% ordering, step 2/4 are guaranteed to have a full list of indirect/direct ancestors
% and thus all transitive edges will be removed
function red_adjmat = dag_trans_red(adjmat)
    % % Check whether adjacency matrix is square
    V=size(adjmat,1);
    assert(V==size(adjmat,2),'Adjacency matrix has to be square');
    
    % % get topological ordering of the nodes
    node_indices = 1:V;
    node_order=[];
    start_nodes = node_indices(sum(adjmat,1)==0);%nodes having no incoming edges
    sort_adj=adjmat; % copy of the intitial matrix
    while ~isempty(start_nodes)
        n = start_nodes(1);
        start_nodes(1)=[];
        node_order=[node_order,n];
        successors = node_indices(logical(adjmat(n,:)));
        for mi=1:numel(successors)
            m=successors(mi);
            sort_adj(n,m)=0;
            if sum(sort_adj(:,m))==0
                start_nodes=[start_nodes,m];
            end
        end
    end
    % if there's still an edge left, the graph must have contained a cycle
    assert(sum(sum(sort_adj))==0,'The supplied graph has at least one cycle! This method works for DAGs only')
    % free memory of adjacency matrix copy used for finding the
    % toplogical order
    clear sort_adj;
    
    % % find transitive reduction
    % reorder adjacency matrix according to topological ordering
    red_adjmat = adjmat(node_order,node_order);
    ancs = zeros(V); % ancestors of each node
    for i=1:V
        % get indirect ancestors of node i
        ancs(i,:)=sum(ancs(logical(red_adjmat(:,i)),:),1)>0;
        % remove edges pointing directly to i from i's indirect ancestors
        red_adjmat(:,i)=red_adjmat(:,i) & ~ancs(i,:)';
        % add the direct predecessors of i to its ancestors
        ancs(i,:)=ancs(i,:)|red_adjmat(:,i)';
    end
    % restore initial matrix order
    re_order=1:V;
    for i=1:V
        re_order(node_order(i))=i;
    end
    red_adjmat = red_adjmat(re_order,re_order);
end