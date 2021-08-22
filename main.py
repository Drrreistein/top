#!/usr/bin/python3.8

import numpy as np
from numpy.core.fromnumeric import _sum_dispatcher

from geneticalgorithm import geneticalgorithm as ga

import matplotlib.pyplot as plt
from fem import FEM_Solver
from plotter import window

fem = FEM_Solver()
original_edges = fem.conn

def sphere(x):
  return np.sum(x ** 2)

def nodes_connected(edges, nodes):

    res = np.zeros(len(nodes))
    ind = 0
    for n in nodes:
        for e in edges:
            if n==e[0] or n==e[1]:
                res[ind]=1
                break
        if res[ind]==0:
            return False
        ind += 1
    return True

def valid_edges(new_edges, support_nodes):

    nodes_edges = dict()
    for edge in new_edges:
        if nodes_edges.get(edge[0]) is None:
            nodes_edges[edge[0]] = set()
        nodes_edges[edge[0]].add(tuple(edge)) 

        if nodes_edges.get(edge[1]) is None:
            nodes_edges[edge[1]] = set()
        nodes_edges[edge[1]].add(tuple(edge)) 

    modify = False
    for node in nodes_edges.keys():
        if len(nodes_edges[node])==1 and node not in support_nodes:
            # print(f"single degree node")
            modify=True
            while True:
                edge = nodes_edges[node].pop()
                if node == edge[0]:
                    next_node = edge[1]
                else:
                    next_node = edge[0]
                nodes_edges[next_node].remove(edge)
                if len(nodes_edges[next_node])==1 and next_node not in support_nodes:
                    node = next_node
                else:
                    break
    if modify:
        res = set()
        for value in nodes_edges.values():
            res = res.union(value)
        new_edges = np.array(tuple(res))
    return new_edges

####################################### cost function #######################

def fem_cost(x):
    num_edges = int(np.sum(x))
    
    new_edges = np.zeros((num_edges, 2), dtype=int)

    ind=0
    for i in range(len(x)):
        if x[i]:
            new_edges[ind] = original_edges[i]
            ind += 1
    new_new_edges = valid_edges(new_edges, [1,2,3])
    
    if not nodes_connected(new_new_edges, [ 10]):
        return np.Inf

    try:
        all_dist = fem.get_displacement(new_new_edges)
        max_disp = np.max(np.abs(all_dist))
        # sum_dist = np.sum(np.abs(all_dist))
        all_force = fem.get_force(new_new_edges)
        max_force = np.max(np.abs(all_force))
        sum_force = np.sum(np.abs(all_force))
        if max_force==0:
            return np.Inf
    except:
        # print(f"singularity matrix")
        return np.Inf

    cost = fem.compute_total_edges_length(new_edges) 
    if max_disp>0.1:
        cost += 200 * (max_disp-0.1)

    if max_force/fem.A>300:
         cost += 2 * (max_force/fem.A-300)

    return cost * sum_force

n_agents = 100
n_variables = len(fem.conn)

####################################### geneticalgorithm #######################

algorithm_param = {'max_num_iteration': 1000,
                   'population_size':30,
                   'mutation_probability':0.1,
                   'elit_ratio': 0.1,
                   'crossover_probability': 0.5,
                   'parents_portion': 0.5,
                   'crossover_type':'uniform',
                   'max_iteration_without_improv':None}
model=ga(function=fem_cost,dimension=n_variables,variable_type='bool',\
    algorithm_parameters=algorithm_param)
model.run()
best_agent = model.best_variable

####################################### nlopt #######################

num_edges = int(np.sum(best_agent))
new_edges = np.zeros((num_edges, 2), dtype=int)
ind=0
for i in range(len(best_agent)):
    if best_agent[i]==1:
        new_edges[ind] = original_edges[i]
        ind += 1

win = window(fem)
win.plot_nodes(fem.node)
win.plot_edges(new_edges, fem.node)
win.plot_support(fem.isol, fem.node)
win.plot_forces(fem.F, fem.node)

new_edges = valid_edges(new_edges, [1,2,3])

print(f"displacement: \n {np.round(fem.get_displacement(new_edges),2)}")
stress = fem.get_force(new_edges)
print(f"forces: \n{np.round(stress,3)}")
# print(f"populations: \n{model.all_pops}")

# plot state after force loading

# win.plot_nodes(fem.node_after, color='blue',alpha=0.5)
colors = []
max_stress = np.max(np.abs(stress))
for f in stress:
    cl = np.interp(np.abs(f), [0, max_stress], [0.8,0])
    if f>0:
        colors.append(np.array([cl,cl,1]))
    elif f<0:
        colors.append(np.array([1,cl,cl]))
    else:
        colors.append(np.array([1,0.5,1]))
win.plot_edges(new_edges, fem.node_after, color=colors, alpha=1, linestyle='solid', linewidth=2)
win.plot_support(fem.isol, fem.node)
win.plot_forces(fem.F, fem.node_after)
plt.show()
