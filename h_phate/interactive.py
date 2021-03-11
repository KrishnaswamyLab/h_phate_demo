from .utils import setdiff1d

import numpy as np
from scipy import sparse
import graphtools
import phate
import scprep

class InteractiveCoarseGraph:
    def __init__(self, cg_op):
        self.cg_op = cg_op
        self.initialize()
        
    def initialize(self):
        self.level_nodes = self.cg_op.create_nodes_dict(initialized=True)
        
    def n_nodes(self):
        return np.sum([len(nodes) for nodes in self.level_nodes])

    def level_sizes(self):
        return [len(self.level_nodes[i]) for i in range(self.cg_op.n_levels)]
    
    def level_boundaries(self):
        return np.concatenate([[0], np.cumsum(self.level_sizes())])
    
    def select_nodes(self, select_idx):
        select_nodes = self.cg_op.create_nodes_dict()
        level_boundaries = self.level_boundaries()
        for i in range(self.cg_op.n_levels):
            level_select_idx = select_idx[(select_idx >= level_boundaries[i]) & 
                                          (select_idx < level_boundaries[i+1])]
            level_select_idx -= level_boundaries[i]
            select_nodes[i] = self.level_nodes[i][level_select_idx]
        return select_nodes
    
    def set_operate_nodes(self, nodes1, nodes2, operator):
        new_nodes = self.cg_op.create_nodes_dict()
        for i in range(self.cg_op.n_levels):
            new_nodes[i] = operator(nodes1[i], nodes2[i])
        return new_nodes
    
    def refine(self, select_idx, remove_deselected=False):
        select_nodes = self.select_nodes(select_idx)
        refined_nodes = self.cg_op.refine(select_nodes)
        if not remove_deselected:
            deselect_nodes = self.set_operate_nodes(self.level_nodes, select_nodes, setdiff1d)
            self.level_nodes = self.set_operate_nodes(refined_nodes, deselect_nodes, np.union1d)
        else:
            self.level_nodes = refined_nodes
    
    def remove(self, select_idx):
        select_nodes = self.select_nodes(select_idx)
        self.level_nodes = self.set_operate_nodes(self.level_nodes, select_nodes, setdiff1d)

    def node_sizes(self):
        return np.concatenate([self.cg_op.node_sizes[level][nodes] 
                               for level, nodes in enumerate(self.level_nodes)])

    def gene_expression(self):
        return sparse.hstack([self.cg_op.gene_expression[level][:,nodes]
                          for level, nodes in enumerate(self.level_nodes)])
    
    def adjacency(self):
        level_boundaries = self.level_boundaries()
        N = level_boundaries[-1]
        A = self.cg_op.adjacency
        A_combined = sparse.lil_matrix((N, N))
        for to_level in range(self.cg_op.n_levels):
            to_level_idx = np.arange(level_boundaries[to_level], level_boundaries[to_level+1])
            if len(to_level_idx) == 0:
                continue
            for from_level in range(0, to_level+1):
                from_level_idx = np.arange(level_boundaries[from_level], level_boundaries[from_level + 1])
                if len(from_level_idx) == 0:
                    continue
                A_submatrix = A[(to_level, from_level)][self.level_nodes[to_level]][:,self.level_nodes[from_level]]
                graphtools.matrix.set_submatrix(A_combined, to_level_idx, from_level_idx, A_submatrix)
                if from_level != to_level:
                    graphtools.matrix.set_submatrix(A_combined, from_level_idx, to_level_idx, A_submatrix.T)

        A_combined = graphtools.matrix.set_diagonal(A_combined, 1)
        return A_combined
    
    def embed(self):
        A = self.adjacency()
        node_sizes = self.node_sizes()
        Y = phate.PHATE(knn_dist='precomputed_affinity', verbose=0).fit_transform(A)
        return Y, node_sizes
    
    def embed_genes(self):
        X = self.gene_expression()
        gene_list = self.cg_op.gene_list
        gene_idx = np.arange(len(gene_list))
        X, gene_list, gene_idx = scprep.filter.filter_empty_cells(X, gene_list, gene_idx)
        Y = phate.PHATE(knn_max=15, verbose=0).fit_transform(X)
        return Y, gene_list, gene_idx