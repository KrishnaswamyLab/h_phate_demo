import plotly.graph_objs as go
import ipywidgets as widgets
import matplotlib as mpl
import numpy as np

from .utils import Lock

class PlotlyDashboard():
    
    def __init__(self, icg_op):
        self.icg_op = icg_op
        self.selected_idx = None
        self.lock = False
        self.initialize_figure()
        self.initialize_gene_figure()
        self.initialize_buttons()
        self.lock = Lock().lock
    
    def initialize_figure(self):
        Y, node_sizes = self.icg_op.embed()
        self.fig = go.FigureWidget(layout=dict(
            dragmode="lasso", hovermode="closest", 
            title=dict(text="Cells Graph"),
            xaxis=dict(showgrid=False, visible=False), 
            yaxis=dict(showgrid=False, visible=False)
        ))
        self.fig.add_scatter(x=Y[:,0], y=Y[:,1], mode='markers', 
                        marker=dict(
                            size=node_sizes, 
                            line=dict(width=0), 
                            opacity=0.4,
                            sizemin=1,
                            sizeref=2. * np.max(node_sizes) / (6 ** 2) ,
                        ))
        self.fig.data[0].on_selection(self.select)
        self.fig.data[0].on_click(self.click_cells_point)
    
    def initialize_gene_figure(self):
        Y, gene_list, self.gene_idx = self.icg_op.embed_genes()
        node_sizes = 10
        self.gene_fig = go.FigureWidget(layout=dict(
            dragmode="lasso", hovermode="closest", 
            title=dict(text="Genes Graph"),
            xaxis=dict(showgrid=False, visible=False), 
            yaxis=dict(showgrid=False, visible=False)
        ))
        self.gene_fig.add_scatter(x=Y[:,0], y=Y[:,1], text=gene_list, mode='markers', 
                        marker=dict(
                            size=node_sizes, 
                            line=dict(width=0), 
                            opacity=0.4,
                            sizemin=1,
                            sizeref=2. * np.max(node_sizes) / (6 ** 2) ,
                        ))
        self.gene_fig.data[0].on_click(self.click_genes_point)
        self.gene_fig.data[0].on_selection(self.click_genes_point)
    
    def initialize_buttons(self):
        self.expand_button = widgets.Button(description="Expand")
        self.zoom_button = widgets.Button(description="Zoom")
        self.filter_button = widgets.Button(description="Filter")
        self.reset_button = widgets.Button(description="Reset")
        self.genes_button = widgets.Button(description="Rebuild gene graph")
        self.expand_button.on_click(self.click_expand)
        self.zoom_button.on_click(self.click_zoom)
        self.filter_button.on_click(self.click_filter)
        self.reset_button.on_click(self.click_reset)
        self.genes_button.on_click(self.click_genes_button)
    
    def select(self, trace, points, selector):
        self.click_cells_point(trace, points, selector)
        self.selected_idx = np.array(points.point_inds, dtype=int)
    
    def rebuild(self):
        Y, node_sizes = self.icg_op.embed()
        with self.fig.batch_update():
            self.fig.data[0].x = Y[:,0]
            self.fig.data[0].y = Y[:,1]
            self.fig.data[0].marker['size'] = node_sizes
            self.fig.data[0].marker['sizeref'] = 2. * np.max(node_sizes) / (6 ** 2)
            self.fig.data[0].selectedpoints = None
            self.selected_idx = None
    
    def rebuild_genes(self):
        Y, gene_list, self.gene_idx = self.icg_op.embed_genes()
        with self.gene_fig.batch_update():
            self.gene_fig.data[0].x = Y[:,0]
            self.gene_fig.data[0].y = Y[:,1]
            self.gene_fig.data[0].text = gene_list
            self.gene_fig.data[0].selectedpoints = None
    
    def color_by_gene(self, gene_idx):
        expression = scprep.utils.toarray(self.icg_op.gene_expression().tocsr()[self.gene_idx[gene_idx]].sum(axis=0)).flatten()
        expression -= expression.min()
        expr_max = expression.max()
        if expr_max > 0:
            expression /= expr_max
        colors = [mpl.colors.to_hex(c) for c in mpl.cm.inferno(expression)]
        self.fig.data[0].marker.color = colors
    
    def color_by_cluster(self, cell_idx):
        expression = scprep.utils.toarray(self.icg_op.gene_expression().tocsc()[:,cell_idx].sum(axis=1)).flatten()
        expression -= expression.min()
        expr_max = expression.max()
        if expr_max > 0:
            expression /= expr_max
        colors = [mpl.colors.to_hex(c) for c in mpl.cm.inferno(expression)]
        self.gene_fig.data[0].marker.color = colors

    def expand(self, zoom=False, filter=False):
        if self.selected_idx is None:
            print("Select points with the Lasso tool to expand.")
            return
        with self.lock() as lock_fn:
            lock_fn()
            if filter:
                self.icg_op.remove(self.selected_idx)
            else:
                self.icg_op.refine(self.selected_idx, remove_deselected=zoom)

            self.rebuild()
        
    def click_expand(self, b):
        self.expand()

    def click_zoom(self, b):
        self.expand(zoom=True)

    def click_filter(self, b):
        self.expand(filter=True)

    def click_reset(self, b):
        with self.lock() as lock_fn:
            lock_fn()
            self.icg_op.initialize()
            self.rebuild()
            self.rebuild_genes()
    
    def click_genes_button(self, b):
        with self.lock() as lock_fn:
            lock_fn()
            self.rebuild_genes()
    
    def click_genes_point(self, trace, points, selector):
        with self.lock() as lock_fn:
            lock_fn()
            self.color_by_gene(points.point_inds)
    
    def click_cells_point(self, trace, points, selector):
        with self.lock() as lock_fn:
            lock_fn()
            self.color_by_cluster(points.point_inds)
        
    @property
    def dashboard(self):
        return widgets.VBox([
            widgets.HBox([
                self.fig, 
                self.gene_fig,
            ]),
            widgets.HBox([
                self.expand_button, 
                self.zoom_button, 
                self.filter_button, 
                self.reset_button,
                self.genes_button,
            ])
        ])