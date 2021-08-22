#!/usr/bin/python3.8
from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage, AnnotationBbox
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np

from fem import FEM_Solver

fem = FEM_Solver()

class window(object):
    def __init__(self, fem_solver:FEM_Solver):
        self.fem = fem_solver
        self.fig, self.ax = plt.subplots()
        self.init_plot()
    
    def init_plot(self):
        plt.axis('equal')
        # plt.axis('off')
        plt.xlabel(f"x")
        plt.ylabel(f"y")

    def plot_nodes(self, nodes, color='gray',alpha=0.2):
        plt.scatter(nodes[:,0], nodes[:,1], c=color, alpha=alpha)
        # plot index of nodes
        for i in range(len(nodes)):
            plt.text(nodes[i,0], nodes[i,1], f"{i+1}")

    def plot_edges(self, edges, nodes, color='gray', alpha=0.5, linestyle='dashed',linewidth=1):
        if isinstance(color, str):
            colors = np.repeat(color, len(edges))
        else:
            colors = color
        for i, e in enumerate(edges):
            node1 = e[0]-1
            node2 = e[1]-1
            plt.plot([nodes[node1][0], nodes[node2][0]], [nodes[node1][1], nodes[node2][1]], color=colors[i],\
                alpha=alpha, linestyle=linestyle, linewidth=linewidth)

    def plot_forces(self, forces, nodes):
        for n in range(len(nodes)*2):
            if forces[n]:
                force_node = nodes[int(n/2)]
                if np.mod(n,2)==0:
                    # horizontal force
                    plt.arrow(force_node[0], force_node[1], 0.2*np.sign(forces[n]), 0, \
                        width=0.001, head_width=0.05, head_length=0.05)
                    plt.text(force_node[0]+0.2*np.sign(forces[n]), force_node[1], f"{np.abs(forces[n])} N")
                else:
                    # veritcal force
                    plt.arrow(force_node[0], force_node[1], 0, 0.2*np.sign(forces[n]), \
                        width=0.001, head_width=0.05, head_length=0.05)
                    plt.text(force_node[0]+0, force_node[1]+0.2*np.sign(forces[n]), f"{np.abs(forces[n])} N")

    def plot_support(self, isol, nodes):
        arr_lena_h = mpimg.imread('./docs/triangle_h.png')
        arr_lena_v = mpimg.imread('./docs/triangle_v.png')

        for n in range(len(nodes)*2):
            if n not in isol:
                support_node = nodes[int(n/2)]
                if np.mod(n,2):
                    imagebox = OffsetImage(arr_lena_h, zoom=0.04 )
                    ab = AnnotationBbox(imagebox, (support_node[0], support_node[1]), \
                        pad=0, frameon=False, box_alignment=[0.5,1])
                else:
                    imagebox = OffsetImage(arr_lena_v, zoom=0.04 )
                    ab = AnnotationBbox(imagebox, (support_node[0], support_node[1]), \
                        pad=0, frameon=False, box_alignment=[1,0.5])
                self.ax.add_artist(ab)

    def plot_triangle(self):
        pass

if __name__=='__main__':
    print(f"hh")

    win = window(fem)
    # plot initial state
    win.plot_nodes(fem.node)
    win.plot_edges(fem.conn, fem.node)
    win.plot_support(fem.isol, fem.node)
    # win.plot_forces(fem.F, fem.node)

    # finite element analysis
    fem.get_displacement(fem.conn)
    fem.get_force(fem.conn)

    # plot state after force loading
    win.plot_nodes(fem.node_after, color='gray',alpha=1)
    win.plot_edges(fem.conn, fem.node_after, color='gray', alpha=1, linestyle='solid', linewidth=2)
    win.plot_forces(fem.F, fem.node_after)

    plt.show()  

