from gurobipy import *
import networkx as nx
import matplotlib.pyplot as plt
from workflow import *
from pilp import *
import greedy
import greedy_ant
import lbbd
import hybrid
import heft
import time
import random


def main():
    num_nodes, dag =20, 1


    result_list = []
    time_list = []
    iter_list = []
    gap_list = []
    for i in range(dag, dag+10):
        workflow_path = './SyntheticSettings/default_dags/%d_nodes/dag_%d.gexf'%(num_nodes, i)
        system_path = './SyntheticSettings/default_dags/platform.gexf'

        dag = load_dag(workflow_path)
        platform = nx.read_gexf(system_path)
        # edge_list = list(platform.edges())
        # edge_list = random.Random(500).sample(edge_list, 2)
        # for (u, v) in edge_list:
        #     platform[u][v]['bandwidth'] = 100000
        # platform.remove_edges_from(edge_list[::2])
        #nx.draw(platform, with_labels=True)
        #plt.show()

        for u, v in platform.edges():
            pass
            # platform[u][v]['bandwidth'] /= 4
        for u, v in dag.edges():
            pass
            # dag[u][v]['data'] /= 2
            # dag[u][v]['data'] *= 4

        for nd in dag.nodes():
            pass
            # dag.node[nd]['Core6'] = dag.node[nd]['Core1'] // 5
        
        nooff_rate, localoff_rate = 0, 0
        noff_nodes = random.Random(401).sample(dag.nodes(), round(dag.number_of_nodes()*nooff_rate))
        localoff_nodes = random.Random(500).sample([i for i in dag.nodes() if not i in noff_nodes],
                                                   round(dag.number_of_nodes()*localoff_rate))
        for t in noff_nodes:
            for r in platform.nodes():
                if r != 'Core1':
                    dag.node[t][r] = 10000
                    # dag[u][v]['data'] *= 4

        for t in localoff_nodes:
                dag.node[t]['Core6'] = 10000
        
        #print('no offloading nodes:', noff_nodes)
        #print('local offloading nodes:', localoff_nodes)
        tl = 3600



        ts = time.time()
        # result,_ = lb, _ = greedy_ant.greedy_ant(dag, platform, 40)
        result, lb = hybrid.ilp(dag, platform, tl=tl)
        # result, lb = hybrid.hybrid(dag, platform, tl=tl, with_ilp=True)
        # result, lb = hybrid.hybrid(dag, platform, tl=tl, with_ilp=False)
        # result = lbbd.lbbd(dag,  platform)
        # result = lb = greedy.greedy(dag, platform)
        # result, _ = lb, _ = heft.heft(dag, platform)

        # result = hybrid.ilp(dag, platform, tl=tl)
        duration = time.time() - ts
        result_list.append(int(result))
        time_list.append(int(duration))
        gap = (result-lb)/result
        gap_list.append(round(gap, 5))
        # iter_list.append(iter)
        print("makespan = ", result_list)
        print("time_cost = ", time_list)
        print("gap = ", gap_list)
    print(result_list)
  

if __name__ == "__main__":
    main()