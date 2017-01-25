%% Generate a testcase for the dag_trans_red function:
% 1. compile a random DAG with up to <nodes> nodes and up to <edges> edges not containing any transitive edges
% 2. add up to <trans_edges> transitive edges that do not violate the DAG property to the random DAG created in step 1
% Return the created non-transitive DAG <adjm_nontrans> as well as the transitive DAG <adjm_trans> which is reducible to the non-transitive one as a test case for transitive reduction algorithms
function [adjm_nontrans,adjm_trans] = generate_testcase(nodes, edges, trans_edges)
    % % construct DAG without transitive edges
    node_indices=1:nodes;
    % initialize as empty graph with <nodes> nodes
    adjm_nontrans = zeros(nodes);
    ancs = zeros(nodes); % matrix of ancestors
    childs = zeros(nodes); % matrix of children
    % try <edges> times to add an edge
    for e=1:edges
        i=randi(nodes-1,1); % source node for new edge
        j=randi([i+1,nodes],1); % target node for new edge
        % only add edge if there's not yet a direct edge:
        % - from any ancestor of j to a child of i AND
        % - from any ancestor of i to a child of j AND
        % - from any ancestor of i to j AND
        % - from i to any child of j
        if sum(sum(adjm_nontrans(logical(ancs(j,:)),logical(childs(i,:)))))==0 && ...
           sum(sum(adjm_nontrans(logical(ancs(i,:)),logical(childs(j,:)))))==0 && ...
           sum(adjm_nontrans(logical(ancs(i,:)),j))==0 && ...
           sum(adjm_nontrans(i,logical(childs(j,:))))==0
            % add edge
            adjm_nontrans(i,j)=1;
            % update ancestor and children relationship in the graph
            ancs(j,i)=1;
            ancs(j,:)=ancs(j,:)|ancs(i,:);
            childs(i,j)=1;
            childs(i,:)=childs(i,:)|childs(j,:);
            childrenj=node_indices(logical(childs(j,:)));
            for ci=1:numel(childrenj)
                c=childrenj(ci);
                ancs(c,:)=ancs(c,:)|ancs(j,:);
            end
            ancestorsi=node_indices(logical(ancs(i,:)));
            for ai=1:numel(ancestorsi)
                a=ancestorsi(ai);
                childs(a,:)=childs(a,:)|childs(i,:);
            end
        end
    end
    
    % % create new graph by adding transitive edges
    adjm_trans=adjm_nontrans;
    % try to add <trans_edges> edges in total
    for e=1:trans_edges
        i=randi(nodes,1); % source node
        tmp=i;
        % make a random walk from source node
        while true
            succs = node_indices(logical(adjm_trans(tmp,:)));
            if isempty(succs)
                break
            end
            tmp=succs(randi(numel(succs),1));
            % stop random walk with probability 20%
            if rand(1)>.8
                break
            end
        end
        % if end node from random walk is different than source
        if tmp~=i
            % add an edge from source to end node
            adjm_trans(i,tmp)=1;
        end
    end
end