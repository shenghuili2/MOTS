import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import random
import math
import pandas as pd
import copy
import re
import os
import crs
import networkx
class TaskGraph(nx.DiGraph):
    def __init__(self, graph_data, crs,
                 data_lb=10, data_ub=30, comp_lb = 10,
                 comp_ub = 20, speedup = 5, r_rate = 0.5):
        super().__init__(graph_data)
        self.crs = crs
        self.process_lb = comp_lb
        self.process_ub = comp_ub
        self.speedup = speedup
        self.r_rate = r_rate
        self.data_lb = data_lb
        self.data_ub = data_ub


    def __init_data_size(self):
        for ed in self.edges():
            self[ed[0]][ed[1]]['weight'] = \
                int(random.uniform(self.data_lb, self.data_ub))

    def set_local(self):
        for tsk in self.nodes():
            for rbt in self.crs.nodes():
                if rbt == 'cloud':
                        self.node[tsk][rbt] = 10000

    def init_cloud_process(self, flag='cloud'):
        for tsk in self.nodes():
            for rbt in self.crs.nodes():
                if rbt == 'cloud':
                    if flag == 'cloud':
                        self.node[tsk][rbt] = ((self.process_lb+self.process_ub) / 2) / self.speedup
                    else:
                        self.node[tsk][rbt] = 10000
    def init_local_process(self, flag='cloud'):
        for tsk in self.nodes():
            for rbt in self.crs.nodes():
                if rbt != 'cloud':
                    if flag == 'cloud':
                        self.node[tsk][rbt] = int(random.uniform(self.process_lb, self.process_ub))
                    elif flag == 'all':
                        if self.node[tsk]['r_value'] != crs.R_LOCAl:
                            self.node[tsk][rbt] = 1000

    def set_r_rate(self, rate):

        self.init_cloud_process()
        nr = round(len(self.nodes()) * rate)
        self.r_tasks = set(random.sample(self.nodes(), nr))
        for tsk in self.nodes():
            if tsk in self.r_tasks:
                self.node[tsk]['r_value'] = crs.R_LOCAl
                self.node[tsk]['cloud'] = 500
            else:
                self.node[tsk]['r_value'] = crs.R_CLOUD
        return self

    def reset_r_rate(self, rate):
        nr = round(len(self.nodes()) * rate)
        if nr > len(self.r_tasks):
            nb = nr - len(self.r_tasks)
            rest_ctsk = [n for n in self.nodes() if self.node[n]['r_value'] == crs.R_CLOUD]
            print(rest_ctsk)
            cand_r_tasks = random.sample(rest_ctsk, nb)
            for t in cand_r_tasks:
                self.node[t]['r_value'] = crs.R_LOCAl
                self.node[t]['cloud'] = 500
                self.r_tasks.add(t)
    def generate(self):
        self.__init_data_size()
        self.init_cloud_process()
        self.init_local_process()

        self.nd_list = list(self.nodes())
        self.set_r_rate(self.r_rate)
        return self

    def reset_data_size(self, data_lb, data_ub):
        self.data_lb = data_lb
        self.data_ub = data_ub
        for ed in self.edges():
            self[ed[0]][ed[1]]['weight'] = \
                int(random.uniform(self.data_lb, self.data_ub))

class TaskGraphGenerator(object):
    def __init__(self, crs, path='./dag.dot', di_graph=None,
                 data_lb=10, data_ub=30, comp_lb = 10, comp_ub = 20,
                    speedup = 5, r_rate = 0.5):

        self.crs = crs
        self.process_lb = comp_lb
        self.process_ub = comp_ub
        self.speedup = speedup
        self.r_rate = r_rate

        self.data_lb = data_lb
        self.data_ub = data_ub
        if di_graph == None:
            self.di_graph = self.__from_daggen(path)
        else:
            self.di_graph = di_graph

    def __from_daggen(self, path):
        with open(path, 'r') as f:

            raw_graph = f.read()
        raw_edges = re.findall('\d* -> \d*', raw_graph)
        edge_list = []
        for s in raw_edges:
            src = re.findall(r'\d+ ', s)[0]
            des = re.findall(r' \d+', s)[0]
            edge_list.append((str(src[:-1]), str(des[1:])))
        src = 'src'
        des = 'des'
        dag = nx.DiGraph(edge_list)
        dag.add_node(src)
        dag.add_node(des)
        for nd in dag.nodes():
            if nd == src or nd == des:
                continue
            if not list(dag.predecessors(nd)):
                dag.add_edge(src, nd)
            if not list(nx.descendants(dag, nd)):
                dag.add_edge(nd, des)
        # dag = TaskGraph(dag)
        return dag

    def generate(self):
        wf = TaskGraph(self.di_graph, crs=self.crs,
                            data_lb=self.data_lb, data_ub=self.data_ub,
                            comp_lb=self.process_lb, comp_ub=self.process_ub,
                            speedup=self.speedup, r_rate=self.r_rate)
        return  wf.generate()


class CRSFactory(object):
    def __init__(self, n_rbt):
        '''

        :param n_rbt: Number of robots
        '''
        self.n_rbt = n_rbt
        self.rbt_list = ['rbt'+str(i) for i in range(n_rbt)]
        self.crs = nx.DiGraph()
        self.node = 'cloud'


    def set_local_bandwidth(self, bd):
        '''

        :param bd: local bandwidth
        :return:
        '''
        self.localbd = bd


    def set_cloud_bandwidth(self, mean=3, var=1):
        '''
        R2C bandwidth, follows Gauss distribution
        :param mean: mean value of Gauss distribution
        :param var: variance value of Gauss distribution
        :return:
        '''
        self.cloudbdmean = mean
        self.cloudbdvar = var


    def create(self):
        self.crs.add_nodes_from(self.rbt_list)
        r2r_edges = []
        for i in self.rbt_list:
            for j in self.rbt_list:
                if i != j:
                    r2r_edges.append((i, j, self.localbd))
        self.crs.add_weighted_edges_from(r2r_edges)

        r2c_edges = []
        for i in self.rbt_list:
            bd = abs(random.gauss(self.cloudbdmean, self.cloudbdvar))
            r2c_edges.append(('cloud', i, bd))
            bd = abs(random.gauss(self.cloudbdmean, self.cloudbdvar))
            r2c_edges.append((i, 'cloud', bd))
        self.crs.add_weighted_edges_from(r2c_edges)

        return self.crs


    def write(self, path='./crs'):
        nx.write_gexf(self.crs, path)

def random_dag(nodes, den=0.7):
    n, den = nodes, den
    cmd =  './daggen/daggen ' + '-n ' + str(n) + ' --density ' + str(den) + ' fat 0.8 --dot --jump 2 > ./dag.dot' 
    os.system(cmd)
    with open('./dag.dot', 'r') as f:
        raw_graph = f.read()
    raw_edges = re.findall('\d* -> \d*', raw_graph)
    edge_list = []
    for s in raw_edges:
        src = re.findall (r'\d+ ', s)[0]
        des = re.findall(r' \d+', s)[0]
        edge_list.append((int(src[:-1]), int(des[1:])))
    src = 'src'
    des = 'des'
    dag = nx.DiGraph(edge_list)
    dag.add_node(src)
    dag.add_node(des)
    for nd in dag.nodes():
        if nd == src or nd == des:
            continue
        if not list(dag.predecessors(nd)):
            dag.add_edge(src, nd)
        if not list(nx.descendants(dag, nd)):
            dag.add_edge(nd, des)
    return dag



if __name__ == '__main__':

    crs = CRSFactory(4)
    crs.set_local_bandwidth(5)
    crs.set_cloud_bandwidth()
    crs = crs.create()
    nx.write_gexf(crs, './crs.gexf')

    tg = TaskGraph()
    tg.set_comp_cost(crs)
    wf = tg.create()
    nx.draw(wf, with_labels=True)
    plt.show()
    nx.write_gexf(wf, './wf.gexf')
    # generate_workflow(50)
    # cloudRoboticSystem()