from __future__ import print_function
from gurobipy import *
import networkx as nx
import copy
import math
from cpm import CPM
from itertools import combinations, chain
from greedy import *
from mows import *
import time
from ortools.sat.python import cp_model
import collections
import heft
import traceback
from gurobipy import *
import networkx as nx
import matplotlib.pyplot as plt
from workflow import *
from pilp import *
from greedy import *
# import lbbd
# import lbbd_
import hybrid
import heft

class MProblem(object):
    def __init__(self, crs, wf,
                 vType=GRB.BINARY):
        self.num_knp = 0
        self.crs = crs
        self.wf = wf
        self.m = Model()
        self.vType = vType
        self.keycut = 0
        self.nd_list = list(nx.topological_sort(wf))
        self.des_node = self.nd_list[-1]
        self.ub = 100000
        self.criticalPath = []
        self.DL = {}
        self.x_jk = self.m.addVars(wf.nodes(), crs.nodes(), vtype=self.vType, name='x_jk')
        self.t_s = self.m.addVars(wf.nodes(), vtype=GRB.CONTINUOUS, name='t_s')
        self.t_e = self.m.addVars(wf.nodes(), vtype=GRB.CONTINUOUS, name='t_e')
        self.ved_index = list(crs.edges()) + [(n, n) for n in crs.nodes()]

        self.f_ijkl = self.m.addVars(wf.edges(), self.ved_index,
                                     vtype=self.vType, name='f')

        self.n_cover_cut = 0
        self.n_ctc_cut = 0
        self.flag = True
        self.cp_cut_flag = True

        # formula-5
        self.m.addConstrs((self.t_s[j] >= self.t_e[i] +
                           sum(self.f_ijkl[i, j, k, l] *
                               round(wf[i][j]['data'] / crs[k][l]['bandwidth'])
                               if k != l else 0
                               for (k, l) in self.ved_index)
                           for j in wf.nodes()
                           for i in wf.predecessors(j)), 'formula-5')

        # formula-8
        self.m.addConstrs((self.t_e[j] == self.t_s[j] + sum(self.x_jk[j, k] * wf.node[j][k]
                                                            for k in crs.nodes())
                           for j in wf.nodes()), 'formula-6')

        # formula-3
        self.m.addConstrs(((self.f_ijkl[i, j, k, l] <= self.x_jk[i, k])
                           for (i, j) in wf.edges()
                           for (k, l) in self.ved_index),
                          'formula-9')
        # formula-3
        self.m.addConstrs(((self.f_ijkl[i, j, k, l] <= self.x_jk[j, l])
                           for (i, j) in wf.edges()
                           for (k, l) in self.ved_index),
                          'formula-10')

        # formula-4
        for (i, j) in wf.edges():
            for (k, l) in self.ved_index:
                self.m.addConstr(self.x_jk[i, k] + self.x_jk[j, l] - 1 <= self.f_ijkl[i, j, k, l], '')

        # formula-1
        self.m.addConstrs(sum(self.x_jk[j, k] for k in crs.nodes()) == 1
                          for j in wf.nodes())

        self.m.addConstr(self.x_jk['src', 'Core1'] == 1)
        self.m.addConstr(self.x_jk['des', 'Core1'] == 1)
        self.m.setObjective(self.t_e[self.des_node], GRB.MINIMIZE)

        # self.init_RL()
        # self.init_RL_downward()
        self.init_RL_down()
        self.init_di()
        self.log_off()

    def allow_cp_cuts(self, flag=True):
        self.cp_cut_flag = flag

    def log_off(self):
        self.m.setParam("LogToConsole", False)

    def log_on(self):
        self.m.setParam("LogToConsole", True)

    def optim(self):
        self.log_off()
        # self.m.addConstrs(sum(self.x_jk[i, r] * self.wf.node[i][r] for i in self.wf.nodes()) <= self.ub
        #                   for r in self.crs.nodes())
        self.m.optimize()
        self.rbt_occupy = {}
        for r in self.crs.nodes():
            t_list = [k[0] for (k, v) in self.x_jk.items()
                      if k[1] == r and v.X == 1]
            t_list.sort(key=lambda x: self.t_s[x].X)
            self.rbt_occupy[r] = t_list

        # max_single_core = max([sum(self.wf.node[i][r] for i in self.rbt_occupy[r])
        #                         for r in self.crs.nodes()])
        # print('max_single_core', max_single_core)
        self.tsk_placement = {}
        for t in self.wf.nodes():
            for r in self.crs:
                if self.x_jk[t, r].X == 1:
                    self.tsk_placement[t] = r
        return round(self.m.objVal)

    def set_init_xjk(self, xjk):
        for (key, value) in xjk.items():
            self.x_jk[key, value].start = 1

    def set_init_obj(self, obj):
        self.t_e[self.des_node].start = obj

    def _subsets(self, li):
        S = [li]
        li.sort()
        for i in range(1, len(li)):
            if i != i - 1:
                S.append(li[i:])
        return S

    def init_RL(self):
        self.RL = {self.nd_list[0]: 0}
        for u in self.nd_list[1:]:
            win = Model()
            win.setParam("LogToConsole", False)
            x_jk = win.addVars(self.wf.nodes(), self.crs.nodes(), vtype=GRB.BINARY, name='x_jk')
            z_u = win.addVar(vtype=GRB.CONTINUOUS, name='z_u')
            pred = list(self.wf.predecessors(u))
            win.addConstrs(sum(x_jk[j, r] for r in self.crs.nodes()) == 1
                           for j in pred)
            for S in self._subsets(pred):
                for r in self.crs.nodes():
                    R_s = min([self.RL[i] for i in S])
                    win.addConstr(R_s + sum(x_jk[j, r] * self.wf.node[j][r] for j in S) <= z_u)

            win.setObjective(z_u, GRB.MINIMIZE)
            win.write('init_RL.lp')
            win.optimize()
            self.RL[u] = z_u.X

            self.m.addConstr(self.t_s[u] >= self.RL[u])
        print("RL values: ", self.RL)

    def init_RL_downward(self):
        self.RL = {self.nd_list[0]: 0}
        for j in self.nd_list[1:]:
            pred = list(self.wf.predecessors(j))
            self.RL[j] = max([self.RL[i] +
                              min([self.wf.node[i][r] for r in self.crs.nodes()])
                              for i in pred])
        print("RL values: ", self.RL)

    def init_RL_down(self):
        self.RL = {self.nd_list[0]: 0}
        chose = []
        RL_iq = {self.nd_list[0]: {r: 0 for r in self.crs.nodes()}}
        for i in self.nd_list[1:]:
            pred = list(self.wf.predecessors(i))
            RL_i = {}
            for q in self.crs.nodes():
                max_pred = []

                for j in pred:
                    min_ij = []
                    for p in self.crs.nodes():
                        if p != q:
                            min_ij.append(
                                RL_iq[j][p] + self.wf.node[j][p] + self.wf[j][i]['data'] / self.crs[p][q]['bandwidth'])
                        else:
                            min_ij.append(RL_iq[j][p] + self.wf.node[j][p])
                    max_pred.append(min(min_ij))

                RL_i[q] = max(max_pred)
            self.RL[i] = min(RL_i.values())
            RL_iq[i] = RL_i
            chose.append(min(RL_i, key = RL_i.get))
        print(self.RL)

    def update_DL(self, zs):
        for t in self.nd_list:
            self.DL[t] = zs + self.di[t]

    def init_di(self):
        self.di = {self.nd_list[-1]: -1}
        for u in self.nd_list[-2::-1]:
            win = Model()
            win.setParam("LogToConsole", False)
            x_jk = win.addVars(self.wf.nodes(), self.crs.nodes(), vtype=GRB.BINARY, name='x_jk')
            z_u = win.addVar(vtype=GRB.CONTINUOUS, name='z_u')
            succ = list(self.wf.successors(u))
            win.addConstrs(sum(x_jk[j, r] for r in self.crs.nodes()) == 1
                           for j in succ)
            for S in self._subsets(succ):
                for r in self.crs.nodes():
                    R_s = min([self.di[i] for i in S])
                    win.addConstr(R_s + sum(x_jk[j, r] * self.wf.node[j][r] for j in S) <= z_u)

            win.setObjective(z_u, GRB.MINIMIZE)
            win.write('init_di.lp')
            win.optimize()
            self.di[u] = z_u.X

        for t in self.nd_list:
            self.di[t] = -self.di[t]

    def update_DL_upward(self, gz):
        self.DL = {self.nd_list[-1]: 1}
        RL_iq = {self.nd_list[-1]: {r: 0 for r in self.crs.nodes()}}
        for i in self.nd_list[-2::-1]:
            succ = list(self.wf.successors(i))
            RL_i = {}
            for q in self.crs.nodes():
                max_pred = []

                for j in succ:
                    min_ij = []
                    for p in self.crs.nodes():
                        if p != q:
                            min_ij.append(
                                RL_iq[j][p] + self.wf.node[j][p] + self.wf[i][j]['data'] / self.crs[p][q]['bandwidth'])
                        else:
                            min_ij.append(RL_iq[j][p] + self.wf.node[j][p])
                    max_pred.append(min(min_ij))

                RL_i[q] = max(max_pred)
            self.DL[i] = min(RL_i.values())
            RL_iq[i] = RL_i

        for t in self.nd_list:
            self.DL[t] = gz - self.DL[t]


    def updateDL(self, gz):
        nodes = list(nx.topological_sort(self.wf))#.reverse()
        nodes.reverse()
        self.DL = {nodes[0]: gz}
        nd = list(self.crs.nodes())
        a = nd[0]
        DL_iq = {nodes[0]: {n: gz for n in self.crs.nodes()}}
        for i in nodes[1:]:
            succ = list(self.wf.successors(i))
            DL_i = {}
            for q in self.crs.nodes():
                max_succ = []

                for j in succ:
                    min_ij = []
                    for p in self.crs.nodes():
                        if p != q:
                            min_ij.append(
                                DL_iq[j][p] - self.wf.node[j][p] - self.wf[i][j]['data'] / self.crs[p][q]['bandwidth'])
                        else:
                            min_ij.append(DL_iq[j][p] - self.wf.node[j][p])
                    max_succ.append(min(min_ij))

                DL_i[q] = max(max_succ)
            self.DL[i] = max(DL_i.values())
            DL_iq[i] = DL_i


    def rm_knap_cut(self):

        for i in range(self.num_knp):
            cst = self.m.getConstrByName('knp%d' % i)
            self.m.remove(cst)
        self.num_knp = 0

    def knapsack_cut(self):
        self.rm_knap_cut()
        for k in self.crs.nodes():
            # stage 1
            cand = [(i, j) for i in self.nd_list[1:-1]
                    for j in self.nd_list[1:-1]
                    if (self.RL[i] <= self.RL[j] and
                        self.DL[i] <= self.DL[j]) and i != j]
            for (i, j) in cand:
                S = [s for s in self.nd_list[1:-1] if (self.RL[s] >= self.RL[i]
                                                       and self.DL[s] <= self.DL[j])]
                if sum(self.wf.node[s][k] for s in S) >= self.DL[j] - self.RL[i]:
                    self.m.addConstr(sum(self.x_jk[s, k] * self.wf.node[s][k]
                                         for s in S) <= self.DL[j] - self.RL[i], 'knp%d' % self.num_knp)
                    self.num_knp += 1
            # stage 2
            for i in self.nd_list[1:-1]:
                S = [s for s in self.nd_list[1:-1] if (self.RL[i] <= self.RL[s] and
                                                       self.DL[s] <= self.DL[i]) and i != s]
                if len(S) <= 1:
                    continue

                A = min([self.RL[s] for s in S if s != i]) - self.RL[i]
                B = self.DL[i] - max([self.DL[s] for s in S if s != i])
                C = min(A, B)
                if sum(self.wf.node[s][k] for s in S) >= self.DL[i] - self.RL[i] - C:
                    self.m.addConstr(sum(self.x_jk[s, k] * self.wf.node[s][k]
                                         for s in S) <= self.DL[i] - self.RL[i] - C,
                                     'knp%d' % self.num_knp)
                    self.num_knp += 1
            # stage 3
            for i in self.nd_list[1:-1]:
                for j in self.nd_list[1:-1]:
                    if i == j:
                        continue
                    if self.wf.node[i][k] + self.wf.node[j][k] > \
                            max(self.DL[j] - self.RL[i], self.DL[i] - self.RL[j]):
                        self.m.addConstr(self.x_jk[i, k] + self.x_jk[j, k] <= 1, 'knp%d' % self.num_knp)
                        self.num_knp += 1

    def write(self):
        self.m.write('master.lp')

    def cover_cut(self, ass_keys):
        cnst = self.m.addConstr(sum([self.x_jk[k] for k in ass_keys])
                                <= len(self.nd_list) - 1, 'cover_cut%d' % self.n_cover_cut)
        self.n_cover_cut += 1

    def critical_path_cut(self, m_sub):
        # Find critical path
        sdag = copy.deepcopy(self.wf)
        for u in sdag.nodes():
            sdag.node[u]['weight'] = int(m_sub.t_e[u] - m_sub.t_s[u])

            if u == m_sub.des_node:
                break
            for v in list(sdag.successors(u)):
                ru = m_sub.tsk_placement[u]
                rv = m_sub.tsk_placement[v]
                if ru == rv:
                    sdag[u][v]['weight'] = 0
                else:
                    sdag[u][v]['weight'] = int(sdag[u][v]['data'] / self.crs[ru][rv]['bandwidth'])
        for r in self.crs.nodes():
            for (u, v) in zip(m_sub.rbt_occupy[r][:-1], m_sub.rbt_occupy[r][1:]):
                if u == 'src' or v == 'des' or u == 'des' or v == 'src':
                    continue
                sdag.add_edge(u, v, weight=0)

        sdag = CPM(sdag)

        # Multi critical paths
        path = sdag.criticalPath.nodes()
        self.critical(m_sub, sdag, path)
        for i in self.criticalPath:
            self.m.addConstr(sum(self.x_jk[t, m_sub.tsk_placement[t]]
                                for t in i
                                if t != 'src' and t != 'des'
                                )
                            <= len(i) - 3, 'path_cut%d' % self.n_ctc_cut)
            self.n_ctc_cut += 1



        # total = sum(self.wf.node[t][m_sub.tsk_placement[t]]
        #                          for t in path)
        # subsetlsit = [total - self.wf.node[t][m_sub.tsk_placement[t]]
        #               for t in path if t != 'src' and t != 'des']
        # print("total: ", total, [i for i in subsetlsit if i > self.ub])
        # sdag.draw()
        self.m.addConstr(sum(self.x_jk[t, m_sub.tsk_placement[t]]
                             for t in path
                             if t != 'src' and t != 'des'
                             )
                         <= len(path) - 3, 'path_cut%d' % self.n_ctc_cut)
        self.n_ctc_cut += 1

    def dfs_critical(self, sdag, path, p, i):
        if i == 'des':
            self.criticalPath.append(list(p))
        else:
            for j in sdag.successors(i):
                if j in path:
                    p.append(j)
                    self.dfs_critical(sdag, path, p, j)
                    p.remove(j)

    def critical(self, m_sub, sdag, path):
        self.criticalPath.clear()
        p, cost = ['src'], []
        self.dfs_critical(sdag, path, p, 'src')
        critic = self.criticalPath
        placement = m_sub.tsk_placement
        for i in range(len(critic)):
            for j in range(len(critic[i])):
                if critic[i][j] == 'src':
                    cost.append(0)
                    continue
                if placement[critic[i][j]] == placement[critic[i][j-1]]:
                    cost[i] = cost[i] + self.wf.node[critic[i][j]][placement[critic[i][j]]]
                else:
                    pred = critic[i][j-1]
                    succ = critic[i][j]

                    cost[i] = cost[i] + self.wf.node[succ][placement[succ]] + \
                              self.wf[pred][succ]['data'] / self.crs[placement[succ]][placement[pred]]['bandwidth']
        temp_path = []
        for i in range(len(cost)):
            if cost[i] == max(cost):
                temp_path.append(critic[i])
        self.criticalPath = temp_path

    def benders_cut(self, m_sub):
        if self.flag:
            self.flag = False
            return

        self.knapsack_cut()
        if self.cp_cut_flag:
            self.critical_path_cut(m_sub)
            ass_keys = [k for (k, v) in self.x_jk.items() if v.X >= 0.99]
            self.cover_cut(ass_keys)
        self.write()


class SCProblem(object):
    def __init__(self, crs, wf, ub=30000):
        self.tsk_placement = {}
        self.crs = crs
        self.wf = wf
        self.ub = int(ub)
        self.m = Model()
        self.DL = {}
        self.nd_list = list(nx.topological_sort(wf))
        self.des_node = self.nd_list[-1]
        self.t_s = {nd: 0 for nd in wf.nodes()}
        self.t_e = {nd: 0 for nd in wf.nodes()}
        self.x_jk = self.m.addVars(wf.nodes(), crs.nodes(), vtype=GRB.BINARY, name='x_jk')
        self.y_ij = self.m.addVars(wf.nodes(), wf.nodes(), vtype=GRB.BINARY, name='y_jk')
        self.m.setParam("LogToConsole", False)
        self.rbt_assignment = {}

    def log_off(self):
        self.m.setParam("LogToConsole", False)

    def log_on(self):
        self.m.setParam("LogToConsole", True)

    def cp_satoptim(self):
        # Instantiate a cp solver.
        model = cp_model.CpModel()
        t_ce = {}
        vars = []

        task_type = collections.namedtuple('task_type', 'start end interval')

        # Job Shop Scheduling
        # Variables
        horizon = self.ub
        all_tasks = {}
        for id in self.t_s.keys():
            core = self.tsk_placement[id]
            start_var = model.NewIntVar(0, horizon, 'start_%s_%s' % (id, core))
            duration = self.wf.node[id][core]
            end_var = model.NewIntVar(0, horizon, 'end_%s_%s' % (id, core))
            interval_var = model.NewIntervalVar(start_var, duration, end_var,
                                                'interval_%s_%s' % (id, core))
            all_tasks[id] = task_type(
                start=start_var, end=end_var, interval=interval_var)

        # Creates sequence variables and add disjunctive constraints.
        for core in self.crs.nodes():
            self.rbt_assignment[core] = [t for t in self.nd_list if self.tsk_placement[t] == core]
            intervals = []
            for tsk in self.rbt_assignment[core]:
                intervals.append(all_tasks[tsk].interval)
            model.AddNoOverlap(intervals)

        # Add precedence contraints.
        for j in self.nd_list:
            for i in self.wf.predecessors(j):
                if self.tsk_placement[i] == self.tsk_placement[j]:
                    model.Add(
                        all_tasks[j].start >= all_tasks[i].end)

                else:
                    k, l = self.tsk_placement[i], self.tsk_placement[j]
                    model.Add(
                        all_tasks[j].start >= all_tasks[i].end +
                        round(self.wf[i][j]['data'] / self.crs[k][l]['bandwidth']))

        model.Minimize(all_tasks[self.des_node].end)

        # Solve model.
        solver = cp_model.CpSolver()
        status = solver.Solve(model)

        if status == cp_model.OPTIMAL:
            # Print out makespan.
            self.optimal_obj = solver.ObjectiveValue()
            # print('Optimal Schedule Length: %i' % solver.ObjectiveValue())
            # print()

            for id in self.t_s.keys():
                # if collector.Value(best_solution, x[i][j]) == 1:
                # self.t_s[id] = collector.Value(best_solution,t_cs[id])
                self.t_s[id] = solver.Value(all_tasks[id].start)
                self.t_e[id] = solver.Value(all_tasks[id].end)
                # print("ts%s=%d, te%s=%d" % (id, self.t_s[id], id, self.t_e[id]))

            self.rbt_occupy = {}
            for r in self.crs.nodes():
                t_list = [k for (k, v) in self.tsk_placement.items()
                          if v == r]
                t_list.sort(key=lambda x: self.t_s[x])
                self.rbt_occupy[r] = t_list

        return self.optimal_obj

    def set_task_placement(self, k):
        self.tsk_placement[k[0]] = k[1]

    def set_start(self):
        pass

    def write(self, name='sub_model.lp'):
        self.m.write(name)


def lbbd(wf, crs):
    ts = time.time()

    #ub = 3121  # heft.heft(wf, crs)
    ub = heft.heft(wf, crs)

    m_master = MProblem(crs, wf)
    m_master.ub = ub
    lb = m_master.optim()

    # lb = 2011

    # ub = 5730#heft.heft(wf, crs)

    vgz = z_best = ub

    iter = 0
    print('ub, lb: ', ub, lb)

    m_sub = SCProblem(crs, wf)
    for (k, v) in m_master.x_jk.items():
        if v.X >= 0.8:
            m_sub.set_task_placement(k)

    # m_sub.set_start()
    # m_sub.log_on()
    z = m_sub.cp_satoptim()
    z_best = min(z, z_best)

    # z_best, lb = 1756, 1723
    while ub >= lb:


        ub = z_best
        MProblem.ub = ub
        # gz = int((lb + max(ub, vgz)) / 2)
        gz = math.ceil((lb + min(ub, vgz)) / 2)

        # m_master.updateDL(gz - 1)
        m_master.update_DL_upward(gz - 1)
        # m_master.update_DL(gz - 1)
        m_master.benders_cut(m_sub)
        iter += 1
        try:
            m_master.m.reset()
            m_obj = m_master.optim()
            vgz = gz

            m_sub = SCProblem(crs, wf)
            for (k, v) in m_master.x_jk.items():
                if v.X >= 0.9:
                    m_sub.set_task_placement(k)

            # m_sub.set_start()
            z = m_sub.cp_satoptim()
            z_best = min(z, z_best)
            m_master.allow_cp_cuts()

        except:
            # info = sys.exc_info()
            # print(info[0], ":", info[1])
            lb = gz + 1

            m_master.allow_cp_cuts(False)

            print('infeasible: ub: %d, lb: %d, z: %d, z_best: %d, gz: %d, vgz: %d, iter: %d'
                  % (ub, lb, z, z_best, gz, vgz, iter))
            # iter += 1
            if vgz < lb:
                vgz = lb + 1
            continue
        print('ub: %d, lb: %d, z: %d, z_best: %d, gz: %d, vgz: %d, iter: %d'
              % (ub, lb, z, z_best, gz, vgz, iter))
        # iter += 1
    print(iter)
    return ub, iter

if __name__ == "__main__":
    num_nodes, num_cores, dag = 20, 5, 11

    result_list = []
    time_list = []
    iter_list = []
    for i in range(dag, dag + 10):
        workflow_path = "./SyntheticSettings/DAGs/%d nodes/%d Cores/dag%d_%dn_%dc.gexf" % \
                        (num_nodes, num_cores, i, num_nodes, num_cores)
        system_path = "./SyntheticSettings/Systems/%dCoreSys.gexf" % num_cores

        dag = load_dag(workflow_path)
        # nx.draw(dag, with_labels=True)
        # plt.show()
        platform = nx.read_gexf(system_path)

        for u, v in dag.edges():
            pass
            # dag[u][v]['data'] = 0
            # dag[u][v]['data'] *= 160
        ts = time.time()
        ilp_model = ILPSolver(dag, platform)
        # result = ilp_model.optimize()
        # result = lbbd.lbbd(dag, platform)
        result = lbbd(dag, platform)
        duration = time.time() - ts
        # result = heft.heft(dag, platform)
        result_list.append(result)
        time_list.append(duration)
        # iter_list.append(iter)
        print("makespan:", result_list)
        print("time cost: ", time_list)
        print("iteration: ", iter_list)