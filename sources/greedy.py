import array, random
# from deap import creator, base, tools, algorithms
import networkx as nx
from functools import partial
import math
import numpy as np
import itertools
import copy
import csv
import pandas as pd
from gurobipy import *
import time


def heuristic(wf, pfm):
    tsks = list(nx.topological_sort(wf))
    des_node = tsks[-1]
    te = {}
    rbt = {}
    rbt_a = {k : 0 for k in pfm.nodes()}
    for u in tsks:
        tsk_finish = float('inf')
        choice_rbt = None
        for k in pfm.nodes():
            tsk_ready = max([0] +
                            [(te[v] +
                             round(wf[v][u]['data'] / pfm[rbt[v]][k]['bandwidth']))
                             if rbt[v] != k
                             else te[v]
                             for v in wf.predecessors(u)])
            tsk_start = max(tsk_ready, rbt_a[k])
            cur_finish = tsk_start + wf.node[u][k]
            if cur_finish < tsk_finish:
                tsk_finish = cur_finish
                rbt[u] = k
                te[u] = tsk_finish
                choice_rbt = k
        # print(choice_rbt, u)
        rbt_a[choice_rbt] = te[u]
    return te[des_node], rbt

class GreedySolver(object):
    def __init__(self, crs, wf):
        self.crs = crs
        self.wf = wf
    def run(self):
        res, _ = heuristic(self.crs, self.wf)
        return res

def greedy(dag, pfm):
    res, rbt  = heuristic(dag, pfm)
    # print('greedy objVal:', res, rbt)
    return res