# CoDO-mv

## Application
Given a background graph G, a set of sample subgraphs of various sizes, which has an overlap subgraph, 
this method uses both the size of the overlap set as well as its density to assess the statistical
significance of the overlap subgraph.

## Usage
src/main.c provides a simple example for the case with 2 sets. This case reduces to CoDO formulation. 
To use C implementation, use the following syntax:

**void CoDO(int x, int nL, int *L, int n, int m_overlap, int m_union, double *p);**

* x:         number of elements overlap between all subsets
* nL:        number of subsets
* L:         subset sizes
* n:         background size
* m_overlap: Number of edges in intersection
* m_union:   Number of edges in union
* p:         output probability

There is also a Matlab interface that is accessible using mex_CoDO function:

**mex_CoDO( subgraph_size_list, size_of_overlap, total_no_vertices, num_edges_in_overlap_subgraph, num_edges_in_union_of_sets);**


## Reference:
Combining Density and Overlap (CoDO): A New Method for Assessing the Significance of Overlap Among Subgraphs
https://arxiv.org/abs/1605.06167
