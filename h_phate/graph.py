import graphtools
import scprep
import graph_coarsening
import numpy as np
from scipy import sparse

class CoarseGraph:
    def __init__(self, K=10, r=0.9, method='variation_edges', algorithm='greedy', hvg_percentile=90):
        self.K = K
        self.r = r
        self.method = method
        self.algorithm = algorithm
        self.hvg_percentile = hvg_percentile
        
    def _build_map_refine(self, map_coarsen):
        """Builds the map from a coarse level to a fine level.
        
        Returns
        -------
        map_refine : dict
            map_refine[(i+1, i)][j] gives the indices in level i
            that map to index j in level i + 1
        """
        map_refine = dict()
        for i in range(self.n_levels-1):
            coarsener = self.map_coarsen[(i, i+1)]
            try:
                max_coarsening = coarsener.shape[1]
            except AttributeError:
                max_coarsening = np.max([len(c) for c in coarsener])
                coarsener = np.array([np.concatenate([c, np.repeat(-1, max_coarsening-len(c))]) for c in coarsener])

            first_nodes = coarsener[:,0]
            remaining_nodes = coarsener[:,1:]
            coarsener = coarsener[np.argsort(first_nodes)]
            nodes_deleted = np.unique(np.concatenate(remaining_nodes))
            nodes_retained = np.setdiff1d(np.arange(self.level_sizes[i]), nodes_deleted)
            refiner = np.full((len(nodes_retained), max_coarsening), -1)
            refiner[np.arange(len(nodes_retained)),0] = nodes_retained
            refiner[np.isin(nodes_retained, first_nodes)] = coarsener
            map_refine[(i+1, i)] = refiner
        return map_refine
    
    def _build_node_sizes(self):
        node_sizes = dict()
        node_sizes[0] = np.ones(self.level_sizes[0])
        for i in range(1, self.n_levels):
            map_refine = self.map_refine[(i, i-1)]
            node_sizes[i] = node_sizes[i-1][map_refine[:,0]]
            expand_idx = map_refine[:,1] != -1
            node_sizes[i][expand_idx] += node_sizes[i-1][map_refine[expand_idx,1]]
        return node_sizes
    
    def _build_gene_expression(self):
        gene_expression = dict()
        data_hvg = scprep.select.highly_variable_genes(self.data, percentile=self.hvg_percentile)
        self.gene_list = data_hvg.columns.to_numpy()
        gene_expression[0] = sparse.lil_matrix(scprep.utils.to_array_or_spmatrix(data_hvg))
        for i in range(1, self.n_levels):
            map_refine = self.map_refine[(i, i-1)]
            gene_expression[i] = gene_expression[i-1][map_refine[:,0]]
            expand_idx = map_refine[:,1] != -1
            gene_expression[i][expand_idx] += gene_expression[i-1][map_refine[expand_idx,1]]
            gene_expression[i][expand_idx] /= 2

        for i in gene_expression:
            gene_expression[i] = gene_expression[i].tocsr().T
            gene_expression[i] = scprep.normalize.library_size_normalize(gene_expression[i])
        
        return gene_expression
    
    def _build_adjacencies(self):
        A = dict()
        for level in range(self.n_levels):
            G = self.graphs[level]
            A[(level, level)] = G.W
            if level > 0:
                coarsener = self.coarseners[level-1]
                for from_level in range(0, level):
                    A[(level, from_level)] = coarsener @ A[(level-1, from_level)]
                    if level == from_level:
                        A[(level, from_level)] = graphtools.matrix.set_diagonal(A[(level, from_level)], 1)
        return A
        
    def fit(self, data):
        self.data = data
        G = graphtools.Graph(data, n_pca=100, use_pygsp=True)
        _, _, self.coarseners, self.graphs, hierarchy = graph_coarsening.coarsen(
            G, K=self.K, r=self.r, method=self.method, algorithm=self.algorithm
        )
        self.n_levels = len(self.graphs)
        self.level_sizes = [G.N for G in self.graphs] # largest to smallest
        self.map_coarsen = {(i, i+1) : h for i, h in enumerate(hierarchy)}
        self.map_refine = self._build_map_refine(self.map_coarsen)
        self.node_sizes = self._build_node_sizes()
        self.gene_expression = self._build_gene_expression()
        self.adjacency = self._build_adjacencies()
        return self
    
    def create_nodes_dict(self, initialized=False):
        level_nodes = [np.array([], dtype=int)
                            for i in range(self.n_levels)]
        if initialized:
            level_nodes[-1] = np.arange(self.level_sizes[-1])
        return level_nodes
    
    def refine_level(self, nodes, level):
        refined_nodes = self.map_refine[level, level-1][nodes]
        refined_nodes = np.setdiff1d(refined_nodes, [-1])
        return refined_nodes
    
    def refine(self, level_nodes):
        refined_level_nodes = self.create_nodes_dict()
        for i in range(1, self.n_levels):
            refined_level_nodes[i-1] = self.refine_level(level_nodes[i], i)
        refined_level_nodes[0] = np.union1d(refined_level_nodes[0], level_nodes[0])
        return refined_level_nodes