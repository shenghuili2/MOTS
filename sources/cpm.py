import random
import networkx as nx
import matplotlib.pyplot as plt
class CPM(nx.DiGraph):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._dirty = True
        self._makespan = -1
        self._criticalPath = None

    def _forward(self):
        for n in nx.topological_sort(self):
            S = max([self.node[j]['C'] + self[j][n]['weight']
                     for j in self.predecessors(n)], default=0)
            self.node[n]['S']=S
            self.node[n]['C']=S + self.node[n]['weight']

    def _backward(self):
        for n in list(nx.topological_sort(self))[::-1]:
            Cp = min([self.node[j]['Sp'] - self[n][j]['weight']
                      for j in self.successors(n)], default=self._makespan)
            self.node[n]['Sp'] = Cp - self.node[n]['weight']
            self.node[n]['Cp'] = Cp

    def _computeCriticalPath(self):
        G = set()
        for n in self:
            if self.node[n]['C'] == self.node[n]['Cp']:
                G.add(n)
        self._criticalPath = self.subgraph(G)

    @property
    def makespan(self):
        if self._dirty:
            self._update()
        return self._makespan

    @property
    def criticalPath(self):
        if self._dirty:
            self._update()
        return self._criticalPath

    def _update(self):
        self._forward()
        self._makespan = max(nx.get_node_attributes(self, 'C').values())
        self._backward()
        self._computeCriticalPath()
        self._dirty = False

    def draw(self):
        for n in self:
            if n in self._criticalPath:
                self.node[n]['color'] = 'green'
            else:
                self.node[n]['color'] = 'red'
        pos = nx.random_layout(self)
        edge_lb = nx.get_edge_attributes(self, 'weight')
        node_lb = nx.get_node_attributes(self, 'weight')
        node_color = [self.node[n]['color'] for n in self.nodes()]
        nx.draw_networkx_nodes(self, pos=pos, label=node_lb, node_color=node_color)
        nx.draw_networkx_labels(self, pos=pos, labels=node_lb)
        nx.draw_networkx_edges(self, pos=pos)
        nx.draw_networkx_edge_labels(self.nodes(), pos=pos, edge_labels=edge_lb)
        plt.show()


if __name__ == "__main__":
    wf = nx.read_gexf("./dag/3.gexf")
    for (u, v) in wf.edges():
        wf[u][v]['weight'] = random.randint(5, 10)
    for n in wf.nodes():
        wf.node[n]['weight'] = random.randint(2, 8)

    cpmExample = CPM(wf)
    print(cpmExample.makespan)
    print(cpmExample.criticalPath.nodes())
    cpmExample.draw()
