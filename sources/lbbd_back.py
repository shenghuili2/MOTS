from gurobipy import *
import networkx as nx
import matplotlib.pyplot as plt
from workflow import *
import time
from itertools import combinations, chain
import math


M = 100000
TTL = 6000

class Master(object):
    def __init__(self, wf, pfm):
        self.m = Model()
        self.DL = {}
        self.nd_list = list(nx.topological_sort(wf))
        self.des_node = self.nd_list[-1]
        self.t_s = self.m.addVars(wf.nodes(), vtype=GRB.CONTINUOUS, name='t_s')
        self.t_e = self.m.addVars(wf.nodes(), vtype=GRB.CONTINUOUS, name='t_e')
        self.x_jk = self.m.addVars(wf.nodes(), pfm.nodes(), vtype=GRB.BINARY, name='x_jk')
        self.y_ij = self.m.addVars(wf.nodes(), wf.nodes(), vtype=GRB.BINARY, name='y_jk')
        self.rbt_assignment = {}
        # self.m.setParam("LogToConsole", False)

        self.ved_index = list(pfm.edges()) + [(n, n) for n in pfm.nodes()]
        self.f_ijkl = self.m.addVars(wf.edges(), self.ved_index,
                                     vtype=GRB.BINARY, name='f')

        self.m.addConstrs((self.t_s[j] >= self.t_e[i] +
                           sum(self.f_ijkl[i, j, k, l] * wf[i][j]['data'] / pfm[k][l]['bandwidth']
                               if k != l else 0
                               for (k, l) in self.ved_index)
                           for j in wf.nodes()
                           for i in wf.predecessors(j)), '')

        self.m.addConstrs((self.t_e[j] == self.t_s[j] + sum(self.x_jk[j, k] * wf.node[j][k]
                                                            for k in pfm.nodes())
                           for j in wf.nodes()), '')

        self.m.addConstrs(((self.f_ijkl[i, j, k, l] <= self.x_jk[i, k])
                           for (i, j) in wf.edges()
                           for (k, l) in self.ved_index), '')

        self.m.addConstrs(((self.f_ijkl[i, j, k, l] <= self.x_jk[j, l])
                           for (i, j) in wf.edges()
                           for (k, l) in self.ved_index), '')

        for (i, j) in wf.edges():
            for (k, l) in self.ved_index:
                self.m.addConstr(self.x_jk[i, k] + self.x_jk[j, l] - 1 <= self.f_ijkl[i, j, k, l], '')

        self.m.addConstrs(sum(self.x_jk[j, k] for k in pfm.nodes()) == 1
                          for j in wf.nodes())

        self.m.setObjective(self.t_e[self.des_node], GRB.MINIMIZE)

    def optimize(self):
        self.m.optimize()
        return self.m.objVal

    def update_dl(self, zs):
        print(zs)

    def write(self):
        self.m.write('master.lp')

    def solve_dependency(self, opt_x, opt_y):
        pass
        return 0, 0

class Slave(object):
    def __init__(self, wf, pfm):
        self.m = Model()
        self.DL = {}
        self.wf = wf
        self.pfm = pfm
        self.nd_list = list(nx.topological_sort(wf))
        self.des_node = self.nd_list[-1]
        self.t_s = self.m.addVars(wf.nodes(), vtype=GRB.CONTINUOUS, name='t_s')
        self.t_e = self.m.addVars(wf.nodes(), vtype=GRB.CONTINUOUS, name='t_e')
        self.x_jk = self.m.addVars(wf.nodes(), pfm.nodes(), vtype=GRB.BINARY, name='x_jk')
        self.y_ij = self.m.addVars(wf.nodes(), wf.nodes(), vtype=GRB.BINARY, name='y_jk')
        self.assignment = {}
        self.placement = {}

    def log_off(self):
        self.m.setParam("LogToConsole", False)
        self.rbt_assignment = {}

    def log_on(self):
        self.m.setParam("LogToConsole", True)

    def optimize(self):
        for r in self.pfm.nodes():
            self.assignment[r] = [t for t in self.nd_list if self.placement[t] == r]
            for (i, j) in combinations(self.assignment[r], 2):
                self.m.addConstr(self.y_ij[i, j] + self.y_ij[j, i] == 1, '')
                self.m.addGenConstrIndicator(self.y_ij[i, j], True, self.t_s[j] >= self.t_e[i])
                self.m.addGenConstrIndicator(self.y_ij[i, j], False, self.t_s[i] >= self.t_e[j])
        self.m.addConstrs(self.t_e[j] >= self.t_s[j] + self.wf.node[j][self.placement[j]] \
                          for j in self.nd_list)

        for j in self.nd_list:
            for i in self.wf.predecessors(j):
                if self.placement[i] == self.placement[j]:
                    self.m.addConstr(self.t_s[j] >= self.t_e[i])
                else:
                    k, l = self.placement[i], self.placement[j]
                    self.m.addConstr(self.t_s[j] >= self.t_e[i] +
                                     round(self.wf[i][j]['data'] / self.pfm[k][l]['bandwidth']))

        self.m.setObjective(self.t_e[self.des_node], GRB.MINIMIZE)

        #         self.log_off()
        #         self.m.setParam("TimeLimit", 600)
        self.m.optimize()
        self.rbt_occupy = {}
        for r in self.pfm.nodes():
            t_list = [k for (k, v) in self.placement.items()
                      if v == r]
            t_list.sort(key=lambda x: self.t_s[x].X)
            self.rbt_occupy[r] = t_list

        return self.m.objVal

    def set_task_placement(self, task, core):
        self.placement[task] = core

    def set_start(self):
        pass

    def write(self, name='sub_model.lp'):
        self.m.write(name)


def lbbd(dag, platform):
    ts = time.time()
    m_master = Master(dag, platform)
    lb = m_master.optimize()

    ub = 100000
    zs = (lb + ub) / 2
    obj = ub
    opt_x = {}
    opt_y = {}
    print('init master obj_val: ', lb)
    iter = 0
    while abs(ub - lb) > 2:
        print('ub, lb, zs: ', ub, lb, zs)
        m_master.m.reset()
        m_master.update_dl(zs)
        m_master.write()

        try:
            m_obj = m_master.optimize()
            print("master obj_val: ", m_obj)
        except:
            print("master infeasible")
            if abs(lb - zs) <= 0.5:
                lb += 1
            else:
                lb = zs
            zs = (lb + ub) / 2
            continue

        if abs(lb - zs) <= 0.5:
            ub = min(obj, ub)
        m_sub = Slave(dag, platform)
        for (k, v) in m_master.x_jk.items():
            if v.X >= 0.8:
                m_sub.set_task_placement(k[0], k[1])
        m_sub.set_start()

        z_sp = m_sub.optimize()
        print('sub problem objective value: ', z_sp)
        if z_sp < ub:
            m_sub.write()
            for (k, v) in m_master.x_jk.items():
                opt_x[k] = v.X
            for (k, v) in m_sub.y_ij.items():
                opt_y[k] = v.X

            ub = z_sp
            zs = (lb + ub) / 2
        else:
            zs = math.floor((lb + zs) / 2)

        obj = min(obj, z_sp)
        print('cur optimal obj: ', obj)
        ass_keys = [k for (k, v) in m_master.x_jk.items() if v.X >= 0.99]
        m_master.write()

        iter += 1
        print('current execution time: ', time.time() - ts)
        if time.time() - ts > TTL:
            status, obj = m_master.solve_dependency(opt_x, opt_y)
            print('final obj : ', obj)
            if status == GRB.OPTIMAL:
                return obj
            print('hybrid ttl reached! ttl =', TTL)
            return obj
    print('optimal solution found: %d, iteration: %d ' % (obj, iter))
    return round(obj)


if __name__ == '__main__':
    num_nodes, num_cores, dag = 30, 4, 1
    workflow_path = "./SyntheticSettings/DAGs/%d nodes/%d Cores/dag%d_%dn_%dc.gexf" % \
                    (num_nodes, num_cores, dag, num_nodes, num_cores)
    system_path = "./SyntheticSettings/Systems/%dCoreSys.gexf" % num_cores

    dag = load_dag(workflow_path)
    platform = nx.read_gexf(system_path)

    ilp_model = lbbd(dag, platform)
