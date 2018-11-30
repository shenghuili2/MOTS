from workflow import *
import networkx as nx
import time


def init_avail(cores):
    avail = {}
    for i in range(1, cores + 1):
        avail["Core{}".format(i)] = 0
    return avail


def computation_cost(dag, platform):
    computation = {}
    nodes = list(dag.nodes())
    cores = list(platform.nodes())
    for nd in nodes:
        temp = 0
        for cr in cores:
            temp += dag.node[nd][cr]
        temp /= len(cores)
        computation[nd] = temp
    return computation


def communication_cost(dag, platform):
    communication = {}
    nodes = list(dag.nodes())
    cores = list(platform.nodes())
    aver = 0
    for cr1 in cores:
        for cr2 in cores:
            if cr1 != cr2:
                aver += platform[cr1][cr2]['bandwidth']
    aver = aver / (2 * len(cores))
    for nd in nodes:
        nd_i = {}
        for su in list(dag.successors(nd)):
            nd_i[su] = dag[nd][su]['data'] / aver
        communication[nd] = nd_i
    return communication


def all_su_visit(dag, visit, nd):
    flag = 1
    for i in list(dag.successors(nd)):
        if i not in visit:
            flag = 0
            break
    return flag


# formula-8,9
def compute_rank(dag, computation, communication):
    rank, li, visit = {}, [], ['des']
    rank['des'] = computation['des']
    for i in dag.predecessors('des'):
        li.append(i)
    while li:
        temp = []
        for i in li:
            res = [0]
            for su in dag.successors(i):
                res.append(communication[i][su] + rank[su])
            rank[i] = computation[i] + max(res)
            visit.append(i)
        for i in li:
            for pre in dag.predecessors(i):
                if all_su_visit(dag, visit, pre):
                    temp.append(pre)
        li = temp
    return rank


def find_min_EFT(EFT):
    res = ''
    minimum = float("inf")
    for i in EFT:
        if EFT[i] < minimum:
            minimum = EFT[i]
            res = i
    return res, minimum


def heft(dag, platform):
    computation = computation_cost(dag, platform)
    communication = communication_cost(dag, platform)
    avail = init_avail(nx.number_of_nodes(platform))
    rank, AFT, m = {}, {}, {}
    nodes = list(dag.nodes())
    rank = compute_rank(dag, computation, communication)
    del rank['src']
    m['src'] = 'Core1'
    AFT['src'] = 0
    schedule = sorted(rank.items(), key=lambda rank: rank[1], reverse=True)
    for i in range(0, len(schedule)):
        task_i = schedule[i][0]

        EST, EFT, ready = {}, {}, {}
        for j in list(platform.nodes()):
            temp = [0]
            for pre in list(dag.predecessors(task_i)):
                if m[pre] != j:
                    temp.append(AFT[pre] + dag[pre][task_i]['data'] / platform[m[pre]][j]['bandwidth'])
                else:
                    temp.append(AFT[pre])
            EST[j] = max(avail[j], max(temp))
            EFT[j] = EST[j] + dag.node[task_i][j]
        m[task_i], AFT[task_i] = find_min_EFT(EFT)
        avail[m[task_i]] = EST[m[task_i]] + dag.node[task_i][m[task_i]]
    print ("AFT".format(AFT['des']))
    return AFT['des'], tuple([(k, m[k]) for k in m.keys()])