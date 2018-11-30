import os
import re
import networkx as nx
def random_dag(nodes, den=0.4):
    n, den = nodes, den
    cmd =  './daggen/daggen ' + '-n ' + str(n) + ' --density ' + str(den) + ' fat 0.3 --dot --jump 3 > ./dag.dot'
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


if __name__ == "__main__":
    for i in range(3, 100):
        dag = random_dag(i)
        nx.write_gexf(dag, "./dag/%d.gexf" % i)
