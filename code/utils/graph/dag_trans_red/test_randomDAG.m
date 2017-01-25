% generate a testcase:
% - <adj_nontrans> is a randomly generated non-transitive DAG
% - <adj_trans> is a DAG created by adding a few transitive edges to
%   <adj_nontrans>
% For large graphs the testcase generation may run for quite a while, 
% since the code is far from optimized for runtime
[adj_nontrans,adj_trans]=generate_testcase(50,100,50);

% compute transitive reduction of <adj_nontrans>
adj_reduced=dag_trans_red(adj_nontrans);

% test whether the computed transitive reduction is equal to <adj_nontrans>
sprintf('DAG transitive reduction matches true solution? %i',~any(any(adj_reduced~=adj_nontrans)))

% write created DAGs to sif files for visualization in Cytoscape 
% (hierarchical layout is best used to compare the graphs visually)
write_graph_from_adj(adj_nontrans,'true_solution.sif');
write_graph_from_adj(adj_trans,'full_graph.sif');
write_graph_from_adj(adj_reduced,'solution.sif');