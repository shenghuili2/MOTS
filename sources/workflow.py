import networkx as nx


class Workflow(nx.DiGraph):
    def __init__(self, graph):
        super().__init__(graph)
        self._add_sink_node()

    def _add_sink_node(self):
        src = 'src'
        des = 'des'
        self.add_node(src)
        self.add_node(des)

        nd = list(self.nodes())[0]
        cores = list(self.node[nd].keys())[:-1]

        for nd in self.nodes():
            if nd == src or nd == des:
                for c in cores:
                    self.node[nd][c] = 0
                continue

            if not list(self.predecessors(nd)):
                self.add_edge(src, nd, data=0)

            if not list(nx.descendants(self, nd)):
                self.add_edge(nd, des, data=0)




def load_dag(path):
    dag = nx.read_gexf(path)
    return Workflow(dag)