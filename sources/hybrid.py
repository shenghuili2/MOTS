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
import threading
import multiprocessing
import greedy_ant
import matplotlib.pyplot as plt

cross_interval = 3600
def hybrid_callback(model, where):
    if where == GRB.Callback.MIP:
        model._shared_dict['objval'] = min(round(model.cbGet(GRB.Callback.MIP_OBJBST)),
                                          model._shared_dict['objval'])
        model._shared_dict['lb'] = max(round(model.cbGet(GRB.Callback.MIP_OBJBND)),
                                     model._shared_dict['lb'])
        if model._shared_dict['finish_flag']:
            model.terminate()
        if model.cbGet(GRB.Callback.MIP_SOLCNT) and not model._timer_flag:
            model._timer_flag = True
            model._ts = time.time()
        duration = time.time() - model._ts
        if duration > cross_interval and \
                model.cbGet(GRB.Callback.MIP_SOLCNT):
            model.terminate()


class MProblem(object):
    def __init__(self, crs, wf,
                 vType=GRB.BINARY):
        self.num_knp = 0
        self.crs = crs
        self.wf = wf
        self.lb = 0
        self.m = Model()
        self.vType = vType
        self.DL = {}
        self.ilp_on = True
        self.keycut = 0
        self.criticalPath = []
        self.nd_list = list(nx.topological_sort(wf))
        self.des_node = self.nd_list[-1]
        self.ub = 100000
        self.di = {self.nd_list[-1]: 0}
        self.RL = {(self.nd_list[0], k): 0 for k in self.crs.nodes()}
        self.x_jk = self.m.addVars(wf.nodes(), crs.nodes(), vtype=self.vType, name='x_jk')
        self.t_s = self.m.addVars(wf.nodes(), vtype=GRB.CONTINUOUS, name='t_s')
        self.t_e = self.m.addVars(wf.nodes(), vtype=GRB.CONTINUOUS, name='t_e')
        self.ved_index = list(crs.edges()) + [(u, v) for (v, u) in crs.edges()
                                              ] + [(n, n) for n in crs.nodes()]

        # self.f_ijkl = self.m.addVars(wf.edges(), self.ved_index,
        #                              vtype=self.vType, name='f')
        self.f_ijkl = self.m.addVars(wf.edges(), [(u, v) for u in self.crs.nodes()
                                                  for v in self.crs.nodes()],
                                     vtype=self.vType, name='f')

        self.y_ij = self.m.addVars(wf.nodes(), wf.nodes(), vtype=GRB.BINARY, name='y_jk')

        self.cover_cut_set = []
        self.n_cover_cut = 0
        self.n_ctc_cut = 0
        self.flag = True
        self.cp_cut_flag = True

        #
        self.m.addConstrs((self.t_s[j] >= self.t_e[i] +
                           sum(self.f_ijkl[i, j, k, l] *
                               int(wf[i][j]['data'] / crs[k][l]['bandwidth'])
                               if k != l else 0
                               for (k, l) in self.ved_index)
                           for j in wf.nodes()
                           for i in wf.predecessors(j)), 'formula-5')

        #
        self.m.addConstrs((self.t_e[j] == self.t_s[j] + sum(self.x_jk[j, k] * wf.node[j][k]
                                                            for k in crs.nodes())
                           for j in wf.nodes()), 'formula-6')

        for (i, j) in wf.edges():
            self.m.addConstr(sum(self.f_ijkl[i, j, k, l]
                                 for (k, l) in self.ved_index) == 1)

        self.m.addConstrs(((self.f_ijkl[i, j, k, l] <= self.x_jk[i, k])
                           for (i, j, k, l) in self.f_ijkl.keys()),
                          'formula-9')

        self.m.addConstrs(((self.f_ijkl[i, j, k, l] <= self.x_jk[j, l])
                           for (i, j, k, l) in self.f_ijkl.keys()),
                          'formula-10')

        for k in self.crs.nodes():
            for l in self.crs.nodes():
                if k != l and not self.crs.has_edge(k, l):
                    self.m.addConstrs(self.f_ijkl[i, j, k, l] == 0
                                      for i, j in self.wf.edges())

        self.m.addConstrs(self.x_jk[i, k] + self.x_jk[j, l] - 1 <= self.f_ijkl[i, j, k, l]
                          for (i, j, k, l) in self.f_ijkl.keys())

        self.m.setObjective(self.t_e[self.des_node], GRB.MINIMIZE)

        self.log_off()

    def lbbd_init(self):
        # self.init_RL()
        # self.init_RL_down()
        # self.init_di()
        # self.init_di()
        self.init_rl_ilp()
        self.init_di_ilp()
        # self.init_RL_down()
        # self.init_di()
        # self.m.setParam("TimeLimit", 10)

    def allow_cp_cuts(self, flag=True):
        self.cp_cut_flag = flag

    def log_off(self):
        self.m.setParam("LogToConsole", False)

    def log_on(self):
        self.m.setParam("LogToConsole", True)

    def ilp_init(self, pilp=False):
        # if not pilp:
        #     self.lbbd_init()
        for i in self.nd_list:
            for j in self.nd_list:
                if i == j:
                    continue
                else:
                    self.m.addConstr(self.y_ij[i, j] + self.y_ij[j, i] == 1)

                    M = 100000
                    for k in self.crs.nodes():
                        self.m.addConstr(self.t_s[j] +
                                         M * (3 - self.x_jk[i, k] - self.x_jk[j, k] -
                                              self.y_ij[i, j]) >= self.t_e[i])

    def ilp_call_back(self, model, where):
        if where == GRB.Callback.MIP:
            self.shared_dict['objval'] = min(round(model.cbGet(GRB.Callback.MIP_OBJBST)),
                                             self.shared_dict)
            self.shared_dict['lb'] = max(round(model.cbGet(GRB.Callback.MIP_OBJBND)),
                                         self.shared_dict)
            # global timer, bounds
            # duration = time.time() - timer
            # if duration > 10:
            #     bounds.append((int(model.cbGet(GRB.Callback.RUNTIME)),
            #                    model.cbGet(GRB.Callback.MIP_OBJBST),
            #                    model.cbGet(GRB.Callback.MIP_OBJBND)))
            #     timer = time.time()

    def ilp_optim(self, timer, shared_dict=None,
                  covercut_queue=None, with_knapsack_cut=False):
        self.log_on()
        self.shared_dict = shared_dict
        while True:
            for t, r in shared_dict['init_asmt']:
                self.x_jk[t, r].start = 1
            if with_knapsack_cut:
                self.m.addConstr(self.t_e[self.des_node] <= round(shared_dict['objval']) - 1)
                self.update_DL(round(shared_dict['objval']))
                self.knapsack_cut()

            if covercut_queue.qsize() >= 1:
                # self.rm_knap_cut()
                print('cover cut length: ', covercut_queue.qsize())
                while not covercut_queue:
                    i = covercut_queue.get()
                    self.m.addConstr(sum(self.x_jk[t[0], t[1]]
                                         for t in i
                                         )
                                     <= len(i) - 1, 'cover_cut%d' % self.n_cover_cut)
                    # self.n_cover_cut += 1

            self.m._shared_dict = shared_dict
            self.m._timer_flag = False
            self.m._ts = float('inf')
            self.m.optimize(hybrid_callback)

            if shared_dict['finish_flag']:
                exit()
            if self.m.status == GRB.Status.OPTIMAL or \
                self.m.status == GRB.Status.INFEASIBLE:
                shared_dict['lb'] = shared_dict['objval']
                shared_dict['ilp_done'] = True
                print('optimal sulotion found by ILP:', round(self.m.objVal))
                exit()
                

    def optim(self):
        ts = time.time()
        self.m.setParam("TimeLimit", 200)
        # self.log_on()
        self.m.optimize()
        print('assignment time: ', time.time() - ts)
        self.rbt_occupy = {}
        for r in self.crs.nodes():
            t_list = [k[0] for (k, v) in self.x_jk.items()
                      if k[1] == r and v.X == 1]
            t_list.sort(key=lambda x: self.t_s[x].X)
            self.rbt_occupy[r] = t_list

        self.tsk_placement = {}
        for t in self.wf.nodes():
            for r in self.crs:
                if self.x_jk[t, r].X == 1:
                    self.tsk_placement[t] = r
        return round(self.m.objVal)

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

            z_jk = win.addVars(self.wf.nodes(), self.crs.nodes(), vtype=GRB.BINARY, name='x_jk')
            z_u = win.addVar(vtype=GRB.CONTINUOUS, name='z_u')
            pred = list(self.wf.predecessors(u))
            win.addConstrs(sum(x_jk[j, r] for r in self.crs.nodes()) == 1
                           for j in pred)
            for S in self._subsets(pred):
                for r in self.crs.nodes():
                    R_s = min([self.RL[i] for i in S])
                    win.addConstr(R_s + sum(x_jk[j, r] * self.wf.node[j][r] for j in S) <= z_u)

            win.setObjective(z_u, GRB.MINIMIZE)
            win.setParam("TimeLimit", 10)
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

    def update_DL(self, zs):
        for t in self.nd_list:
            self.DL[t] = zs + self.di[t]

    def init_di(self):

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
        # print(self.di)

    def init_di_ilp(self):
        self.di = {(self.nd_list[-1], k): 0 for k in self.crs.nodes()}
        self.di['des'] = 0
        for u in self.nd_list[-2::-1]:
            succ = list(self.wf.successors(u))
            if succ == ['des']:
                self.di[u] = -1
                for k in self.crs.nodes():
                    self.di[u, k] = -1
                continue
            for k in self.crs.nodes():
                win = Model()
                # win.setParam("LogToConsole", False)

                x_jk = win.addVars(succ, self.crs.nodes(), vtype=GRB.BINARY, name='x_jk')
                tr_il = win.addVars(succ, self.crs.nodes(), vtype=GRB.CONTINUOUS, name='tr_il')
                te_il = win.addVars(succ, self.crs.nodes(), vtype=GRB.CONTINUOUS, name='te_il')
                ts_il = win.addVars(succ, self.crs.nodes(), vtype=GRB.CONTINUOUS, name='te_il')
                y_ij = win.addVars(succ, succ, vtype=GRB.BINARY, name='y_jk')
                z_k = win.addVar(vtype=GRB.INTEGER, name='z_k')

                win.addConstrs(sum(x_jk[j, r] for r in self.crs.nodes()) == 1
                               for j in succ)

                for i in succ:
                    for j in succ:
                        if i == j:
                            continue
                        else:
                            win.addConstr(y_ij[i, j] + y_ij[j, i] == 1)

                            win.addConstrs(tr_il[i, l] +
                                           10000 * (3 - x_jk[i, l] - x_jk[j, l] -
                                                    y_ij[i, j]) >= te_il[j, l]
                                           for l in self.crs.nodes())

                    win.addConstrs(ts_il[i, l] == max_(tr_il[i, l], self.di[i, l])
                                   for l in self.crs.nodes())

                    win.addConstrs(ts_il[i, l] + self.wf.node[i][l] <= te_il[i, l]
                                   + 10000 * (1 - x_jk[i, l])
                                   for l in self.crs.nodes())

                    for l in self.crs.nodes():
                        if l == k:
                            win.addConstr(te_il[i, l] <= z_k + 10000 * (1 - x_jk[i, l]))
                        else:
                            win.addConstr(te_il[i, l]
                                          + int(self.wf[u][i]['data']/self.crs[k][l]['bandwidth'])
                                          <= z_k + 10000 * (1 - x_jk[i, l]))
                win.setObjective(z_k, GRB.MINIMIZE)

                win.setParam("LogToConsole", False)
                win.setParam("TimeLimit", 200)
                win.optimize()
                self.di[u, k] = round(z_k.X)
            self.di[u] = min(self.di[u, k] for k in self.crs.nodes())

        for t in self.nd_list:
            self.di[t] = -self.di[t]
        print("di values: ", self.di)

    def init_rl_ilp(self):

        self.RL['src'] = 0
        for u in self.nd_list[1:]:
            pred = list(self.wf.predecessors(u))
            if pred == ['src']:
                self.RL[u] = 1
                for k in self.crs.nodes():
                    self.RL[u, k] = 1
                continue
            for k in self.crs.nodes():
                win = Model()
                # win.setParam("LogToConsole", False)

                x_jk = win.addVars(pred, self.crs.nodes(), vtype=GRB.BINARY, name='x_jk')
                tr_il = win.addVars(pred, self.crs.nodes(), vtype=GRB.INTEGER, name='tr_il')
                te_il = win.addVars(pred, self.crs.nodes(), vtype=GRB.INTEGER, name='te_il')
                ts_il = win.addVars(pred, self.crs.nodes(), vtype=GRB.INTEGER, name='te_il')
                y_ij = win.addVars(pred, pred, vtype=GRB.BINARY, name='y_jk')
                z_k = win.addVar(vtype=GRB.INTEGER, name='z_k')

                win.addConstrs(sum(x_jk[j, r] for r in self.crs.nodes()) == 1
                               for j in pred)

                for i in pred:
                    for j in pred:
                        if i == j:
                            continue
                        else:
                            win.addConstr(y_ij[i, j] + y_ij[j, i] == 1)

                        win.addConstrs(tr_il[i, l] +
                                       10000 * (3 - x_jk[i, l] - x_jk[j, l] -
                                                y_ij[i, j]) >= te_il[j, l]
                                       for l in self.crs.nodes())

                    win.addConstrs(ts_il[i, l] == max_(tr_il[i, l], self.RL[i, l])
                                   for l in self.crs.nodes())

                    win.addConstrs(ts_il[i, l] + self.wf.node[i][l] <= te_il[i, l]
                                   + 10000 * (1 - x_jk[i, l])
                                   for l in self.crs.nodes())

                    for l in self.crs.nodes():
                        if l == k:
                            win.addConstr(te_il[i, l] <= z_k + 10000 * (1 - x_jk[i, l]))
                        else:
                            win.addConstr(te_il[i, l]
                                          + int(self.wf[i][u]['data']/self.crs[l][k]['bandwidth'])
                                          <= z_k + 10000 * (1 - x_jk[i, l]))
                win.setObjective(z_k, GRB.MINIMIZE)

                win.setParam("LogToConsole", False)
                win.setParam("TimeLimit", 200)
                win.optimize()
                self.RL[u, k] = round(z_k.X)
            self.RL[u] = min(self.RL[u, k] for k in self.crs.nodes())
            self.m.addConstr(self.t_s[u] >= self.RL[u])
        print("RL values: ", self.RL)

    # Yufei version
    def init_RL_down(self):
        self.RL = {self.nd_list[0]: 0}
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
        # print(self.RL)

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

        # nx.draw(sdag, with_labels=True)
        # plt.show()
        # Multi critical paths
        path = sdag.criticalPath.nodes()
        self.critical(m_sub, sdag, path)
        # self.cover_cut_generation(m_sub, gz)

        for i in self.criticalPath:
            self.m.addConstr(sum(self.x_jk[t, m_sub.tsk_placement[t]]
                                 for t in i
                                 if t != 'src' and t != 'des'
                                 )
                             <= len(i) - 3, 'path_cut%d' % self.n_ctc_cut)
            self.n_ctc_cut += 1

        self.m.addConstr(sum(self.x_jk[t, m_sub.tsk_placement[t]]
                             for t in path
                             if t != 'src' and t != 'des'
                             )
                         <= len(path) - 3, 'path_cut%d' % self.n_ctc_cut)
        self.n_ctc_cut += 1

    def dfs_critical(self, sdag, path, p, i):
        if i == 'des' and len(path) - len(p) <= 2:
            self.criticalPath.append(list(p))
        else:
            for j in sdag.successors(i):
                if j in path:
                    p.append(j)
                    self.dfs_critical(sdag, path, p, j)
                    p.remove(j)

    def critical(self, m_sub, sdag, path):
        self.criticalPath.clear()

        dag_ = copy.deepcopy(sdag)
        for i in list(dag_.nodes):
            if i not in path:
                dag_.remove_node(i)

        for i in dag_.nodes:
            for j in dag_.successors(i):
                dag_[i][j]['weight'] = 1

        self.criticalPath = [nx.dag_longest_path(dag_)]

        # p, cost = ['src'], []
        # self.dfs_critical(sdag, path, p, 'src')
        # critic = self.criticalPath
        # placement = m_sub.tsk_placement
        # for i in range(len(critic)):
        #     for j in range(len(critic[i])):
        #         if critic[i][j] == 'src':
        #             cost.append(0)
        #             continue
        #         if placement[critic[i][j]] == placement[critic[i][j - 1]]:
        #             cost[i] = cost[i] + self.wf.node[critic[i][j]][placement[critic[i][j]]]
        #         else:
        #             pred = critic[i][j - 1]
        #             succ = critic[i][j]
        #
        #             cost[i] = cost[i] + self.wf.node[succ][placement[succ]] + \
        #                       self.wf[pred][succ]['data'] / self.crs[placement[succ]][placement[pred]]['bandwidth']
        # temp_path = []
        # for i in range(len(cost)):
        #     if cost[i] == max(cost):
        #         temp_path.append(critic[i])
        # self.criticalPath = temp_path

    def cover_cut_generation(self, m_sub, gz, covercut_queue):

        cover_cut_set = []
        for p in self.criticalPath:
            for (i, e) in enumerate(p):
                if m_sub.t_e[e] > self.ub:
                    cut = [(t, m_sub.tsk_placement[t])
                           for t in p[1:i + 1]]
                    cover_cut_set.append(cut)
                    break

        for cp in self.criticalPath:
            p = cp[::-1]
            for (i, e) in enumerate(p):
                if m_sub.t_e[p[0]] - m_sub.t_e[e] > self.ub:
                    cut = [(t, m_sub.tsk_placement[t])
                           for t in p[1:i + 1]]
                    cover_cut_set.append(cut)
                    break

        # self.cover_cut_set.extend(cover_cut_set)

        for i in cover_cut_set:
            covercut_queue.put(tuple(i))
            self.m.addConstr(sum(self.x_jk[t[0], t[1]]
                                 for t in i)
                             <= len(i) - 1, 'cover_cut%d' % self.n_cover_cut)
            self.n_cover_cut += 1

    # edited by cyf
    def cover_cut_computation(self, m_sub, gz):
        # schedule dag
        sdag = copy.deepcopy(self.wf)
        for u in sdag.nodes():
            for v in list(sdag.successors(u)):
                sdag[u][v]['data'] = 0
        for r in self.crs.nodes():
            for (u, v) in zip(m_sub.rbt_occupy[r][:-1], m_sub.rbt_occupy[r][1:]):
                if u == 'src' or v == 'des' or u == 'des' or v == 'src':
                    continue
                sdag.add_edge(u, v, data=0)

        cover_cut_set = []
        for i in self.criticalPath:
            zerodag = copy.deepcopy(sdag)

            ub = 1
            Xc = [('src', m_sub.tsk_placement['src']), (i[ub], m_sub.tsk_placement[i[ub]])]
            sub = SCProblem(self.crs, self.wf)
            while sub.cover_cut_optimize(Xc, self.DL, zerodag, m_sub) and ub < len(i):
                ub += 1
                if i[ub] == 'des':
                    break
                for j in Xc:
                    if j in list(sdag.predecessors(i[ub])):
                        zerodag[j][i[ub]]['data'] = sdag[j][i[ub]]['data']
                Xc.append((i[ub], m_sub.tsk_placement[i[ub]]))
            Xc.append(('des', m_sub.tsk_placement['des']))
            cover_cut_set.append(list(Xc)[1:-1])

        # reverse dag cover cut computation
        reverse_sdag = sdag.reverse()
        reverse_dag = self.wf.reverse()
        end = {}
        for i in self.nd_list:
            end[i] = gz - self.RL[i]

        for i in self.criticalPath:

            zerodag = copy.deepcopy(reverse_sdag)
            ub = len(i) - 2
            Xc = [('des', m_sub.tsk_placement['des']), (i[ub], m_sub.tsk_placement[i[ub]])]
            sub = SCProblem(self.crs, self.wf)

            while sub.cover_cut_optimize(Xc, end, zerodag, m_sub) and ub >= 0:
                ub -= 1
                if i[ub] == 'src':
                    break
                for j in Xc:
                    if j in list(sdag.predecessors(i[ub])):
                        if j not in list(reverse_dag.predecessors(i[ub])):
                            zerodag[j][i[ub]]['data'] = 0
                        else:
                            zerodag[j][i[ub]]['data'] = reverse_dag[j][i[ub]]['data']
                Xc.append((i[ub], m_sub.tsk_placement[i[ub]]))
            Xc.append(('src', m_sub.tsk_placement['src']))
            cover_cut_set.append(list(Xc))

        self.cover_cut_set.extend(cover_cut_set)
        for i in cover_cut_set:
            self.m.addConstr(sum(self.x_jk[t[0], t[1]]
                                 for t in i)
                             <= len(i) - 1, 'cover_cut%d' % self.n_cover_cut)
            self.n_cover_cut += 1

    def benders_cut(self, m_sub, gz, covercut_queue):
        if self.flag:
            self.flag = False
            return
        ts = time.time()
        self.knapsack_cut()
        print('knapsack time', time.time() - ts)
        if self.cp_cut_flag:

            self.critical_path_cut(m_sub)

            # self.cover_cut_computation(m_sub, gz)
            self.cover_cut_generation(m_sub, gz, covercut_queue)
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

    def cp_satoptim(self, ub=None):
        # Instantiate a cp solver.
        model = cp_model.CpModel()
        task_type = collections.namedtuple('task_type', 'start end interval')

        # Job Shop Scheduling
        # Variables
        horizon = 300000000
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
                    model.Add(all_tasks[j].start >= all_tasks[i].end)

                else:
                    k, l = self.tsk_placement[i], self.tsk_placement[j]
                    data_time = round(self.wf[i][j]['data'] / self.crs[k][l]['bandwidth'])
                    model.Add(all_tasks[j].start >= all_tasks[i].end + data_time)

        model.Minimize(all_tasks[self.des_node].end)

        # Solve model.
        solver = cp_model.CpSolver()
        status = solver.Solve(model)

        if status == cp_model.OPTIMAL:
            # Print out makespan.
            self.optimal_obj = solver.ObjectiveValue()
            if self.optimal_obj == 1:
                print('error')
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
        else:
            print(nx.is_directed_acyclic_graph(self.wf))
            return ub
        return self.optimal_obj

    def cp_satoptimize(self):
        # Instantiate a cp solver.
        t_begin = time.time()
        print('constraint programming begines')
        model = cp_model.CpModel()
        task_type = collections.namedtuple('task_type', 'start end interval')

        # Job Shop Scheduling
        # Variables
        horizon = 50000
        all_tasks = {}
        ts = {}
        te = {}
        durations = {}
        y = {}
        for id in self.t_s.keys():
            core = self.tsk_placement[id]
            start_var = model.NewIntVar(0, horizon, 'start_%s_%s' % (id, core))
            cost = self.wf.node[id][core]
            end_var = model.NewIntVar(0, horizon, 'end_%s_%s' % (id, core))
            ts[id] = start_var
            te[id] = end_var
            durations[id] = cost
            for id2 in self.t_s.keys():
                schedule_var = model.NewBoolVar('y_%s_%s' % (id, id2))
                y[id, id2] = schedule_var

        for core in self.crs.nodes():
            self.rbt_assignment[core] = [t for t in self.nd_list if self.tsk_placement[t] == core]
            intervals = []
            for i in self.rbt_assignment[core]:
                for j in self.rbt_assignment[core]:
                    if i != j:
                        if self.wf.has_edge(i, j):
                            model.Add(y[i, j] == 1)
                        model.Add(y[i, j] + y[j, i] == 1)
                        model.Add(ts[j] >= te[i]).OnlyEnforceIf(y[i, j])

        # Add precedence contraints.
        for j in self.nd_list:
            model.Add(te[j] == ts[j] + durations[j])
            for i in self.wf.predecessors(j):
                if self.tsk_placement[i] == self.tsk_placement[j]:
                    model.Add(ts[j] >= te[i])
                else:
                    k, l = self.tsk_placement[i], self.tsk_placement[j]
                    data_time = round(self.wf[i][j]['data'] / self.crs[k][l]['bandwidth'])
                    model.Add(ts[j] >= te[i] + data_time)

        model.Minimize(te[self.des_node])

        # Solve model.
        solver = cp_model.CpSolver()
        status = solver.Solve(model)

        # print('cp time: ', time.time() - t_begin)
        if status == cp_model.OPTIMAL:
            # Print out makespan.
            self.optimal_obj = solver.ObjectiveValue()
            # print('Optimal Schedule Length: %i' % solver.ObjectiveValue())
            # print()

            for id in self.t_s.keys():
                # if collector.Value(best_solution, x[i][j]) == 1:
                # self.t_s[id] = collector.Value(best_solution,t_cs[id])
                self.t_s[id] = solver.Value(ts[id])
                self.t_e[id] = solver.Value(te[id])
                # print("ts%s=%d, te%s=%d" % (id, self.t_s[id], id, self.t_e[id]))

            self.rbt_occupy = {}
            for r in self.crs.nodes():
                t_list = [k for (k, v) in self.tsk_placement.items()
                          if v == r]
                t_list.sort(key=lambda x: self.t_s[x])
                self.rbt_occupy[r] = t_list
        else:
            print('fuck')
        return self.optimal_obj

    def set_task_placement(self, k):
        self.tsk_placement[k[0]] = k[1]

    def set_start(self):
        pass

    def write(self, name='sub_model.lp'):
        self.m.write(name)

    # edited by cyf, cover cut computation's feasibility check
    def cover_cut_optimize(self, Xc, DL, sdag, m_sub):
        # Instantiate a cp solver.
        model = cp_model.CpModel()
        task_type = collections.namedtuple('task_type', 'start end interval')

        # Job Shop Scheduling
        # Variables
        horizon = 50000
        ts = {}
        te = {}
        durations = {}
        y = {}
        for id in m_sub.t_s.keys():
            core = m_sub.tsk_placement[id]
            start_var = model.NewIntVar(0, horizon, 'start_%s_%s' % (id, core))
            if id in Xc:
                cost = m_sub.wf.node[id][core]
            else:
                cost = 0
            end_var = model.NewIntVar(0, horizon, 'end_%s_%s' % (id, core))
            ts[id] = start_var
            te[id] = end_var
            durations[id] = cost
            for id2 in m_sub.t_s.keys():
                schedule_var = model.NewBoolVar('y_%s_%s' % (id, id2))
                y[id, id2] = schedule_var

        for core in self.crs.nodes():
            self.rbt_assignment[core] = [t for t in self.nd_list if m_sub.tsk_placement[t] == core]
            for i in self.rbt_assignment[core]:
                for j in self.rbt_assignment[core]:
                    if i != j:
                        model.Add(y[i, j] + y[j, i] == 1)
                        model.Add(ts[j] >= te[i]).OnlyEnforceIf(y[i, j])

        # Add precedence contraints.
        for j in self.nd_list:
            model.Add(te[j] == ts[j] + durations[j])
            model.Add(te[j] <= int(DL[j]))
            for i in sdag.predecessors(j):
                if m_sub.tsk_placement[i] == m_sub.tsk_placement[j]:
                    model.Add(ts[j] >= te[i])

                else:
                    k, l = m_sub.tsk_placement[i], m_sub.tsk_placement[j]
                    data_time = round(sdag[i][j]['data'] / self.crs[k][l]['bandwidth'])
                    model.Add(ts[j] >= te[i] + data_time)

        model.Minimize(te[self.des_node])
        # Solve model.
        solver = cp_model.CpSolver()
        status = solver.Solve(model)

        if status == cp_model.INFEASIBLE:
            return False
        else:
            return True


bounds = []
timer = None
def ilp_call_back(model, where):
    if where == GRB.Callback.MIP:
        global timer, bounds
        duration = time.time() - timer
        if duration > 10:
            bounds.append((int(model.cbGet(GRB.Callback.RUNTIME)),
                           model.cbGet(GRB.Callback.MIP_OBJBST),
                           model.cbGet(GRB.Callback.MIP_OBJBND)))
            timer = time.time()

def ilp(wf, crs, tl=3600):
    ilp_model = MProblem(crs, wf)
    ilp_model.log_on()
    ilp_model.ilp_init(pilp=True)

    global timer, bounds
    timer = time.time()
    ilp_model.m.setParam("TimeLimit", tl)
    ilp_model.log_on()
    ilp_model.m.optimize(ilp_call_back)
    print(bounds)
    return ilp_model.m.objVal, ilp_model.m.objBound





def hybrid(wf, crs, tl=3600, with_ilp=True):
    ts = time.time()
    global_timer = time.time()
    # tl = 400
    timer = cross_interval
    # ub, _ = greedy_ant.greedy_ant(wf, crs, 40)
    ub, assignment = heft.heft(wf, crs)
    # ub = heft.heft(wf, crs)

    manager = multiprocessing.Manager()
    shared_dict = manager.dict()
    shared_dict['objval'] = ub
    shared_dict['init_asmt'] = assignment
    shared_dict['optimal'] = ub
    shared_dict['ilp_done'] = False
    shared_dict['lb'] = 0
    shared_dict['ilp_flag'] = True
    shared_dict['cover_cut'] = []
    shared_dict['finish_flag'] = False
    covercut_queue = multiprocessing.Queue()


    ilp_model = MProblem(crs, wf)
    ilp_model.ilp_init()

    ilp_model.init_rl_ilp()
    ilp_model.init_di_ilp()
    # ilp_model.init_RL()
    # ilp_model.init_di()
    shared_dict['rl'] = ilp_model.RL
    shared_dict['di'] = ilp_model.di

    ilp_process = multiprocessing.Process(target=ilp_model.ilp_optim,
                                          args=(timer, shared_dict,
                                                covercut_queue, True))  #
    lbbd_process = multiprocessing.Process(target=lbbd,
                                           args=(wf, crs, shared_dict, tl, ub,
                                                 covercut_queue))  #

    lbbd_process.start()
    if with_ilp:
        ilp_process.start()
    while time.time()-ts < tl:
        if shared_dict['lb'] == shared_dict['objval'] or \
            shared_dict['finish_flag']:
            if lbbd_process.is_alive():
                lbbd_process.terminate()
            break
        time.sleep(2)
    if with_ilp:
        ilp_process.join()
    if lbbd_process.is_alive():
        lbbd_process.terminate()
    lbbd_process.join()

    print('final object value:', shared_dict['objval'])
    return shared_dict['objval'], shared_dict['lb']


def lbbd(wf, crs, shared_dict, tl=3600, ub=10000, covercut_queue=None):
    ts = time.time()
    global_timer = time.time()

    m_master = MProblem(crs, wf)
    m_master.di = shared_dict['di']
    m_master.RL = shared_dict['rl']

    lb = m_master.optim()
    vgz = z_best = ub
    iter = 0
    print('ub, lb: ', ub, lb)

    m_sub = SCProblem(crs, wf)
    for (k, v) in m_master.x_jk.items():
        if v.X >= 0.99:
            m_sub.set_task_placement(k)

    # z = m_sub.cp_satoptimize()
    z = m_sub.cp_satoptim()
    z_best = min(z, z_best)

    # edited by cyf
    gz_length = 0
    infeasible_count = 0

    bouds = []
    bound_timer = time.time()
    while ub > lb and time.time() - global_timer < tl and not shared_dict['ilp_done']:

        m_master.ub = ub
        shared_dict['objval'] = ub = min(shared_dict['objval'], ub, z_best)
        shared_dict['lb'] = lb = max(lb, shared_dict['lb'])
        MProblem.ub = ub
        gz = math.ceil((lb + min(ub, vgz)) / 2)

        m_master.update_DL(gz - 1)
        m_master.benders_cut(m_sub, gz - 1, covercut_queue)
        iter += 1
        try:

            m_master.m.reset()
            # m_master.log_on()
            m_obj = m_master.optim()

        except:
            lb = gz + 1
            m_master.allow_cp_cuts(False)

            # edited by cyf
            infeasible_count += 1
            if infeasible_count >= 3:
                gz_length = (ub - lb) // 5
                vgz += gz_length
            duration = time.time() - bound_timer
            if duration > 5:
                bouds.append((time.time() - global_timer, ub, lb))
                bound_timer = time.time()
            print(int(time.time() - global_timer),
                  'infeasible: ub: %d, lb: %d, z: %d, inf_num: %d, gz: %d, vgz: %d, iter: %d, gap: %.4f'
                  % (ub, lb, z, infeasible_count, gz, vgz, iter, (ub - lb) / ub))
            # print(bouds)
            if vgz < lb:
                vgz = lb
            continue

        # edited by cyf
        infeasible_count = 0
        gz_length //= 2
        vgz = gz + gz_length

        m_sub = SCProblem(crs, wf)
        for (k, v) in m_master.x_jk.items():
            if v.X >= 0.9:
                m_sub.set_task_placement(k)
        try:
            z = m_sub.cp_satoptim(ub)
        except:
            z = ub
        # z = m_sub.cp_satoptimize()
        if int(z) == int(ub):
            m_master.allow_cp_cuts(False)
        else:
            m_master.allow_cp_cuts(True)
        if int(z) < int(ub):
            asgmt = tuple([k
                 for (k, v) in m_master.x_jk.items()
                 if v.X >= 0.9])
            shared_dict['init_asmt'] = asgmt
        z_best = min(z, z_best)
        # m_master.allow_cp_cuts()
        bouds.append((time.time() - global_timer, ub, lb))
        bound_timer = time.time()
        print(int(time.time() - global_timer), 'ub: %d, lb: %d, z: %d, z_best: %d, gz: %d, vgz: %d, iter: %d, gap: %.4f'
              % (ub, lb, z, z_best, gz, vgz, iter, (ub - lb) / ub))

        if shared_dict['ilp_done']:
            lb = ub = shared_dict['objval']
    shared_dict['finish_flag'] = True
    print('objective value: ', ub)
    ub = min(ub, z_best, shared_dict['objval'])
    shared_dict['optimal'] = ub
    print(bouds)
    exit()