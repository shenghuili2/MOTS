from gurobipy import *
import networkx as nx
import matplotlib.pyplot as plt
from workflow import *
from pilp import *
from greedy import *
import lbbd
import lbbd_
import hybrid
import heft

# num_nodes, num_cores, dag =20, 5, 12
# num_nodes, num_cores, dag = 20, 5, 14
num_nodes, num_cores, dag = 50, 8, 1


result_list = []
time_list = []
iter_list = []
for i in range(dag, dag+10):
    workflow_path = "./SyntheticSettings/DAGs/%d nodes/%d Cores/dag%d_%dn_%dc.gexf"%\
                    (num_nodes, num_cores, i, num_nodes, num_cores)
    system_path = "./SyntheticSettings/Systems/%dCoreSys.gexf" % num_cores

    dag = load_dag(workflow_path)
    # nx.draw(dag, with_labels=True)
    # plt.show()
    platform = nx.read_gexf(system_path)

    for u, v in dag.edges():
        # if dag[u][v]['data'] == 96:
        # dag[u][v]['data'] /= 4
        # pass
        dag[u][v]['data'] = 0
        #dag[u][v]['data'] *= 160
    tl=3600
    ts = time.time()
    # result, lb = hybrid.ilp(dag, platform, tl=tl)
    result, lb = hybrid.hybrid(dag, platform, tl=tl, with_ilp=True)
    duration = time.time() - ts
    # result = heft.heft(dag, platform)
    result_list.append(result)
    time_list.append(duration)
    # iter_list.append(iter)
    print("makespan:", result_list)
    print("time cost: ", time_list)
    # print("iteration: ", iter_list)
