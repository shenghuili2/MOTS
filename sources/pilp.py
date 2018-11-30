from gurobipy import *
import networkx as nx
import matplotlib.pyplot as plt
from workflow import *

M = 100000
class ILPSolver(object):
    def __init__(self, wf, pfm):
        self.m = Model()
        self.DL = {}
        self.pfm = pfm
        self.nd_list = list(nx.topological_sort(wf))
        self.des_node = self.nd_list[-1]
        self.t_s = self.m.addVars(wf.nodes(), vtype=GRB.CONTINUOUS, name='t_s')
        self.t_e = self.m.addVars(wf.nodes(), vtype=GRB.CONTINUOUS, name='t_e')
        self.x_jk = self.m.addVars(wf.nodes(), pfm.nodes(), vtype=GRB.BINARY, name='x_jk')
        self.y_ij = self.m.addVars(wf.nodes(), wf.nodes(), vtype=GRB.BINARY, name='y_jk')
        self.rbt_assignment = {}
        # self.m.setParam("LogToConsole", False)

        self.ved_index = list(pfm.edges()) + [(u, v) for (v,u) in pfm.edges()
                                              ]+ [(n, n) for n in pfm.nodes()]
        self.f_ijkl = self.m.addVars(wf.edges(), self.ved_index,
                                     vtype=GRB.BINARY, name='f') 

        
        self.m.addConstrs((self.t_s[j] >= self.t_e[i] +
                           sum(self.f_ijkl[i, j, k, l] * int(wf[i][j]['data'] / pfm[k][l]['bandwidth'])
                               if k != l else 0
                               for (k, l) in self.ved_index)
                           for j in wf.nodes()
                           for i in wf.predecessors(j)), 'formula-5')

        
        self.m.addConstrs((self.t_e[j] == self.t_s[j] + sum(self.x_jk[j, k] * wf.node[j][k]
                                                            for k in pfm.nodes())
                           for j in wf.nodes()), 'formula-6')

        
        self.m.addConstrs(((self.f_ijkl[i, j, k, l] <= self.x_jk[i, k])
                           for (i, j) in wf.edges()
                           for (k, l) in self.ved_index),
                          'formula-9')
        
        self.m.addConstrs(((self.f_ijkl[i, j, k, l] <= self.x_jk[j, l])
                           for (i, j) in wf.edges()
                           for (k, l) in self.ved_index),
                          'formula-10')
        
        
        for (i, j) in wf.edges():
            for (k, l) in self.ved_index:
                self.m.addConstr(self.x_jk[i, k] + self.x_jk[j, l] - 1 <= self.f_ijkl[i, j, k, l], '')
       

        self.m.addConstrs(sum(self.x_jk[j, k] for k in pfm.nodes()) == 1
                          for j in wf.nodes())

        for i in self.nd_list:
            for j in self.nd_list:
                if i == j:
                    continue
                else:
                    self.m.addConstr(self.y_ij[i, j] + self.y_ij[j, i] == 1)

                    M = 100000
                    for k in pfm.nodes():
                        self.m.addConstr(self.t_s[j] +
                                         M * (3 - self.x_jk[i, k] - self.x_jk[j, k] -
                                              self.y_ij[i, j]) >= self.t_e[i])

        self.m.addConstr(self.t_e[self.des_node] >= 103)
        # self.m.addConstr(self.t_e[self.des_node] <= 1807)
        self.m.setObjective(self.t_e[self.des_node], GRB.MINIMIZE)

    def optimize(self):
        self.m.optimize()
        placement = [(i, k) for i in self.nd_list
                     for k in self.pfm if self.x_jk[i, k].X == 1]
        ts = [(i, self.t_s[i].X) for i in self.nd_list]
        te = [(i, self.t_e[i].X) for i in self.nd_list]
        return self.m.objVal, self.m.objBound

    def set_schedule(self, S):
        for id, i in enumerate(S):
            if id == len(S) - 1:
                break
            for j in S[id+1:]:
                self.m.addConstr(self.y_ij[i, j] == 1)
                self.m.addConstr(self.y_ij[j, i] == 0)

    def set_time(self, tl=20):
        self.m.setParam("TimeLimit", tl)

    def set_logoff(self):
        self.m.setParam("LogToConsole", False)


if __name__ == '__main__':
    num_nodes, num_cores, dag = 50, 10, 11
    workflow_path = "./SyntheticSettings/DAGs/%d nodes/%d Cores/dag%d_%dn_%dc.gexf"%\
                    (num_nodes, num_cores, dag, num_nodes, num_cores)
    system_path = "./SyntheticSettings/Systems/%dCoreSys.gexf" % num_cores

    dag = load_dag(workflow_path)
    platform = nx.read_gexf(system_path)

    for u, v in dag.edges():
        dag[u][v]['data'] = 0

    ilp_model = ILPSolver(dag, platform)
    result = ilp_model.optimize()
