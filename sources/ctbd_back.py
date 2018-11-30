from gurobipy import *
import networkx as nx
import copy
import math
from cpm import CPM
from itertools import combinations, chain
from greedy import *
from mows import *
import time

class MProblem(object):
    def __init__(self, crs, wf,
                 vType=GRB.BINARY):
        self.num_knp = 0
        self.crs = crs
        self.wf = wf
        self.m = Model()
        self.vType = vType
        self.DL = {}
        self.keycut = 0
        self.nd_list = list(nx.topological_sort(wf))
        self.des_node = self.nd_list[-1]

        self.x_jk = self.m.addVars(wf.nodes(), crs.nodes(), vtype=self.vType, name='x_jk')
        self.t_s = self.m.addVars(wf.nodes(), vtype=GRB.CONTINUOUS, name='t_s')
        self.t_e = self.m.addVars(wf.nodes(), vtype=GRB.CONTINUOUS, name='t_e')
        self.ved_index = list(crs.edges()) + [(n, n) for n in crs.nodes()]

        self.f_ijkl = self.m.addVars(wf.edges(), self.ved_index,
                           vtype=self.vType, name='f')

        self.n_cover_cut = 0
        self.n_ctc_cut = 0
        # formula-3
        self.m.addConstrs((self.t_s[j] >= self.t_e[i] +
                           sum(self.f_ijkl[i, j, k, l] *
                                  round(wf[i][j]['data'] / crs[k][l]['bandwidth'])
                                             if k != l else 0
                                             for (k, l) in self.ved_index)
                      for j in wf.nodes()
                      for i in wf.predecessors(j)), 'formula-5')

        # formula-5
        self.m.addConstrs((self.t_e[j] == self.t_s[j] + sum(self.x_jk[j, k] * wf.node[j][k]
                                                            for k in crs.nodes())
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

        self.m.addConstrs(sum(self.x_jk[j, k] for k in crs.nodes()) == 1
                          for j in wf.nodes())

        self.m.setObjective(self.t_e[self.des_node], GRB.MINIMIZE)

        self.init_RL()
        self.init_di()
        self.log_off()

    def log_off(self):
        self.m.setParam("LogToConsole", False)

    def log_on(self):
        self.m.setParam("LogToConsole", True)

    def optim(self):
        self.m.optimize()
        self.rbt_occupy = {}
        for r in self.crs.nodes():
            t_list = [k[0] for (k, v) in self.x_jk.items()
                      if k[1]==r and v.X == 1]
            t_list.sort(key=lambda x: self.t_s[x].X)
            self.rbt_occupy[r] = t_list
            # print(r, t_list)

        self.tsk_placement = {}
        for t in self.wf.nodes():
            for r in self.crs:
                if self.x_jk[t, r].X == 1:
                    self.tsk_placement[t] = r
        return round(self.m.objVal)

    def set_init_xjk(self,xjk):
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
        self.RL = {self.nd_list[0] : 0}
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
        # print("RL values: ", self.RL)

    def update_DL(self, zs):
        for t in self.nd_list:
            self.DL[t] = zs + self.di[t]
        # print('DL values: ', self.DL)

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

    def rm_knap_cut(self):

        for i in range(self.num_knp):
            cst = self.m.getConstrByName('knp%d' % i)
            self.m.remove(cst)
        self.num_knp = 0

    def knapsack_cut(self):
        self.rm_knap_cut()
        for k in self.crs.nodes():
            # stage 1
            cand = [(i, j) for i in self.wf.nodes()
                            for j in self.wf.nodes()
                                if (self.RL[i] <= self.RL[j] and
                                    self.DL[i] <= self.DL[j])]
            for (i, j) in cand:
                S = [s for s in self.wf.nodes() if (self.RL[s] >= self.RL[i]
                                    and self.DL[s] <= self.DL[j] )]
                if sum(self.wf.node[s][k] for s in S) >= self.DL[j] - self.RL[i]:
                    self.m.addConstr(sum(self.x_jk[s, k]*self.wf.node[s][k]
                                         for s in S) <= self.DL[j] - self.RL[i], 'knp%d'%self.num_knp)
                    self.num_knp += 1
            # #stage 2
            for i in self.nd_list:
                S = [s for s in self.wf.nodes() if (self.RL[i] <= self.RL[s]  and
                                    self.DL[s] <= self.DL[i])]
                if len(S) <= 1:
                    continue

                A = min([self.RL[s] for s in S if s != i]) - self.RL[i]
                B = self.DL[i] - max([self.DL[s] for s in S if s != i])
                C = max(A, B)
                if sum(self.wf.node[s][k] for s in S) >= self.DL[i] - self.RL[i] - C:
                    self.m.addConstr(sum(self.x_jk[s, k] * self.wf.node[s][k]
                                         for s in S) <= self.DL[i] - self.RL[i] - C, 'knp%d'%self.num_knp)
                    self.num_knp += 1
            #stage 3
            for i in self.nd_list:
                for j in self.nd_list:
                    if i == j:
                        continue
                    if self.wf.node[i][k] + self.wf.node[j][k] > \
                                        max(self.DL[j]-self.RL[i], self.DL[i]-self.RL[j]):
                        self.m.addConstr(self.x_jk[i, k] + self.x_jk[j, k] <= 1, 'knp%d'%self.num_knp)
                        self.num_knp += 1

    def write(self):
        self.m.write('master.lp')

    def cover_cut(self, ass_keys):

        cnst = self.m.addConstr(sum([self.x_jk[k] for k in ass_keys])
                         <= len(self.nd_list) - 1, 'cover_cut%d'%self.n_cover_cut)
        self.n_cover_cut += 1

    def critical_path_cut(self, path, sub_model):
        self.m.addConstr(sum(self.x_jk[t, sub_model.tsk_placement[t]]
                                 for t in path) <= len(path) - 1, 'path_cut%d'%self.n_ctc_cut)
        self.n_ctc_cut += 1
        # self.m.getConstrByName()

    def solve_dependency(self, init_x, init_y):
        self.m.update()
        milp = self.m.copy()
        y_ij = milp.addVars(self.wf.nodes(), self.wf.nodes(), vtype=GRB.BINARY, name='y_ij')

        for i in range(self.num_knp):
            cst = milp.getConstrByName('knp%d' % i)
            milp.remove(cst)

        # for i in range(self.n_cover_cut):
        #     cst = milp.getConstrByName('cover_cut%d' % (i))
        #     milp.remove(cst)

        cst = milp.getConstrByName('cover_cut%d' % (self.keycut))
        milp.remove(cst)

        # for i in range(self.n_ctc_cut):
        cst = milp.getConstrByName('path_cut%d' % (self.keycut))
        milp.remove(cst)

        milp.update()
        for i in self.nd_list:
            for j in self.nd_list:
                if i == j:
                    continue
                else:
                    milp.addConstr(y_ij[i, j] + y_ij[j, i] == 1)
                    y_ij[i, j].start = init_y[i, j]
                    # milp.addConstr(y_ij[i, j] == init_y[i, j])
                    M = 1000

                    for k in self.crs.nodes():
                        ts = milp.getVarByName('t_s['+j+']')
                        te = milp.getVarByName('t_e['+i+']')
                        xik = milp.getVarByName('x_jk[' + i + ',' + k + ']')
                        xjk = milp.getVarByName('x_jk['+j + ',' + k+']')
                        milp.addConstr(ts + M * (3 - xik - xjk - y_ij[i, j]) >= te)
                        xik.start = init_x[i, k]
                        # milp.addConstr(xik == init_x[i, k])
        milp.setParam("LogToConsole", True)
        # milp.setParam("LogToConsole", False)
        milp.write('ilp.lp')
        # milp.setParam("TimeLimit", 20)
        milp.optimize()

        sta = milp.getAttr(GRB.Attr.Status)
        obj = milp.objVal
        print(obj)
        return sta, obj


class SCProblem(object):
    def __init__(self, crs, wf):
        self.tsk_placement = {}
        self.crs = crs
        self.wf = wf
        self.m = Model()
        self.DL = {}
        self.nd_list = list(nx.topological_sort(wf))
        self.des_node = self.nd_list[-1]
        self.t_s = self.m.addVars(wf.nodes(), vtype=GRB.CONTINUOUS, name='t_s')
        self.t_e = self.m.addVars(wf.nodes(), vtype=GRB.CONTINUOUS, name='t_e')
        self.x_jk = self.m.addVars(wf.nodes(), crs.nodes(), vtype=GRB.BINARY, name='x_jk')
        self.y_ij = self.m.addVars(wf.nodes(), wf.nodes(), vtype=GRB.BINARY, name='y_jk')
        self.m.setParam("LogToConsole", False)
        self.rbt_assignment = {}

    def log_off(self):
        self.m.setParam("LogToConsole", False)



    def log_on(self):
        self.m.setParam("LogToConsole", True)

    def optim(self):
        for r in self.crs.nodes():
            self.rbt_assignment[r] = [t for t in self.nd_list if self.tsk_placement[t] == r]
            for (i, j) in combinations(self.rbt_assignment[r], 2):
                self.m.addConstr(self.y_ij[i, j] + self.y_ij[j, i] == 1)
                self.m.addGenConstrIndicator(self.y_ij[i, j], True, self.t_s[j] >= self.t_e[i])
                self.m.addGenConstrIndicator(self.y_ij[i, j], False, self.t_s[i] >= self.t_e[j])
        self.m.addConstrs(self.t_e[j] >= self.t_s[j] + self.wf.node[j][self.tsk_placement[j]] \
                            for j in self.nd_list)

        for j in  self.nd_list:
            for i in self.wf.predecessors(j):
                if self.tsk_placement[i] == self.tsk_placement[j]:
                    self.m.addConstr(self.t_s[j] >= self.t_e[i])
                else:
                    k, l = self.tsk_placement[i], self.tsk_placement[j]
                    self.m.addConstr(self.t_s[j] >= self.t_e[i] +
                                     round(self.wf[i][j]['data'] / self.crs[k][l]['bandwidth']))

        self.m.setObjective(self.t_e[self.des_node], GRB.MINIMIZE)

        self.log_off()
        self.m.setParam("TimeLimit", 600)
        self.m.optimize()
        self.rbt_occupy = {}
        for r in self.crs.nodes():
            t_list = [k for (k, v) in self.tsk_placement.items()
                      if v == r]
            t_list.sort(key=lambda x: self.t_s[x].X)
            self.rbt_occupy[r] = t_list

        return round(self.m.objVal)

    def set_task_placement(self, k):
        self.tsk_placement[k[0]] = k[1]

    def set_start(self):
        pass

    def write(self, name='sub_model.lp'):
        self.m.write(name)


def ctbd(wf, crs):
    ts = time.time()
    m_master = MProblem(crs, wf)
    lb = m_master.optim()

    ub = greedy(wf, crs)
    zs = (lb + ub) / 2
    obj = ub
    opt_x = {}
    opt_y = {}
    # print('init master obj_val: ', lb)
    iter = 0
    while abs(ub - lb) > 2:
        print('ub, lb, zs: ' , ub,  lb, zs)
        m_master.m.reset()
        m_master.update_DL(zs)
        m_master.knapsack_cut()
        m_master.write()

        try:
            m_obj = m_master.optim()
        except:
            # print("master infeasible")
            if abs(lb - zs) <= 0.5:
                lb += 1
            else:
                lb = zs
            zs = (lb + ub) / 2
            continue

        if abs(lb - zs) <= 0.5:
            ub = min(obj, ub)
        m_sub = SCProblem(crs, wf)
        for (k, v) in m_master.x_jk.items():
            if v.X >= 0.8:
                m_sub.set_task_placement(k)

        m_sub.set_start()
        z_sp = m_sub.optim()
        # print('sub problem objective value: ', z_sp)
        if z_sp < ub:
            m_sub.write()
            m_master.keycut = m_master.n_cover_cut
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
        m_master.cover_cut(ass_keys)
        m_master.write()

        # Find critical path
        sdag = copy.deepcopy(wf)
        for u in sdag.nodes():
            sdag.node[u]['weight'] = int(m_sub.t_e[u].X - m_sub.t_s[u].X)

            if u  == m_sub.des_node:
                break
            for v in list(sdag.successors(u)):
                ru = m_sub.tsk_placement[u]
                rv = m_sub.tsk_placement[v]
                if ru == rv:
                    sdag[u][v]['weight'] = 0
                else:
                    sdag[u][v]['weight'] = int(sdag[u][v]['data'] / crs[ru][rv]['bandwidth'])
        for r in crs.nodes():
            for (u, v) in zip(m_sub.rbt_occupy[r][:-1], m_sub.rbt_occupy[r][1:]):
                sdag.add_edge(u, v, weight=0)

        sdag = CPM(sdag)
        path = sdag.criticalPath.nodes()
        m_master.critical_path_cut(path, m_sub)
        m_master.write()
        iter += 1
        # print('current execution time: ', time.time() - ts)
        if time.time() - ts > MOWS.TTL:

            status, obj = m_master.solve_dependency(opt_x, opt_y)
            print('final obj : ', obj)
            if status == GRB.OPTIMAL:
                return obj

            print('hybrid ttl reached! ttl =', MOWS.TTL)
            return obj
    print('optimal solution found: %d, iteration: %d '%(obj, iter))
    return round(obj)