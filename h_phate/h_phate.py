from .graph import CoarseGraph
from .interactive import InteractiveCoarseGraph
from .dashboard import PlotlyDashboard

class H_PHATE():
    
    def __init__(self, K=10, r=0.9, method='variation_edges'):
        self.K = K
        self.r = r
        self.method = method
    
    def fit(self, data):
        self.cg_op = CoarseGraph(K=self.K, r=self.r, method=self.method).fit(data)
        self.icg_op = InteractiveCoarseGraph(self.cg_op)
        self.dashboard_op = PlotlyDashboard(self.icg_op)
        return self
    
    @property
    def dashboard(self):
        return self.dashboard_op.dashboard