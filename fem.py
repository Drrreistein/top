#!/usr/bin/python3.8


import numpy as np
from copy import deepcopy

from numpy.core.arrayprint import dtype_is_implied
from numpy.lib.type_check import iscomplex

class FEM_Solver(object):
    """
    based on,
    Shubham Dhanale (2021). MATLAB program for 2D truss analysis (FEM) 
    (https://www.mathworks.com/matlabcentral/fileexchange/75983-matlab-program-for-2d-truss-analysis-fem), 
    MATLAB Central File Exchange. Retrieved August 19, 2021.
    """
    def __init__(self):
                
        # coordinate of all self.nodes
        # self.node = np.array([[0,0],[2,2],[3,0],[0,2]])
        self.node = np.array([[0,0],[0,1],[0,2],[1,0],[1,1],[1,2],[2,0],[2,1],[2,2],[3,0],[3,1],[3,2]])
        self.constraint = np.array([0,1,2,3,4,5])

        # self.node = np.array([[0,0],[0,1],[1,0],[1,1],[2,0],[2,1]])
        #         # fixed degree of freedom
        # self.constraint = np.array([0,1,2,3])

        # number of self.nodes
        self.nn=self.node.shape[0]
        # number of dof
        self.ndof=2*self.nn

        # self.connections of self.nodes - edges
        self.conn = self.full_combination(self.nn)
        # self.conn = np.array([[1, 4],[1, 5],[4, 5],[2, 4]])

        # number of edge
        self.ne=self.conn.shape[0]

        # area of cross section, /mm^2
        self.A=200
        # modulus of elasticity of the element, 
        self.E=29000

        self.isol = np.zeros(self.ndof-len(self.constraint), dtype=int)
        ind=0
        for i in range(self.ndof):
            if i not in self.constraint:
                self.isol[ind] = i
                ind += 1

        # force, /Newton
        self.F=np.zeros(2*self.nn)
        self.F[19]=-5e+4
        # self.F[9]=-1e+3

        # displacement of nodes
        self.d= None
        # force of edges
        self.force=None

    def full_combination(self, node_len):
        edge = np.zeros((int(node_len*(node_len-1)/2), 2), dtype=int)
        ind = 0
        for i in range(1,node_len+1):
            for j in range(i+1, node_len+1):
                edge[ind] = np.array([i, j])
                ind += 1
        return edge

    def get_displacement(self, edges):
        # Defines size of the global stiffness matrix 
        self.K = np.zeros((self.ndof,self.ndof))
        # Defines size of the displacement matrix    
        self.d = np.zeros(self.ndof)
        # calculate displacements of self.nodes
        used_nodes = set()
        for i in range(len(edges)):
            used_nodes.add(edges[i][0])
            used_nodes.add(edges[i][1])
            n0 = edges[i, 0] -1
            n1 = edges[i, 1] -1

            x0 = self.node[n0, 0]
            y0 = self.node[n0, 1]
            x1 = self.node[n1, 0]
            y1 = self.node[n1, 1]

            L = np.sqrt((x0-x1)**2+(y0-y1)**2)
            l = (x1-x0)/L
            m = (y1-y0)/L

            self.ke = self.A*self.E/L * np.array([[l*l,m*l,-l*l,-m*l],[m*l,m*m,-m*l,-m*m],[-l*l,-m*l,l*l,m*l],[-m*l,-m*m,m*l,m*m]])

            sctr=np.array([2*n0, 2*n0+1, 2*n1, 2*n1+1])
            self.K[np.ix_(sctr, sctr)] += self.ke
        for i in range(self.ndof):
            if self.F[i] != 0:
                used_nodes.add(np.floor(i/2)+1)
        isol = []
        for n in self.isol:
            if np.floor(n/2)+1 in used_nodes:
                isol.append(n)
        isol = np.array(isol)
        self.d[isol] = np.linalg.solve(self.K[np.ix_(isol, isol)], self.F[isol])
        # self.d[self.isol] = np.linalg.lstsq(self.K[np.ix_(self.isol, self.isol)], self.F[self.isol])[0]
        
        # position of nodes after force loaded
        self.node_after = np.zeros((self.nn, 2))
        for i in range(self.ndof):
            self.node_after[int(i/2),np.mod(i,2)] = self.node[int(i/2),np.mod(i,2)] + self.d[i]

        return self.d

    def get_force(self, edges):
        if self.d is None:
            print(f"calculate displacement first")
            return None
        self.ne = len(edges)
        self.strain = np.zeros(self.ne)
        self.stress = np.zeros(self.ne)
        self.force = np.zeros(self.ne)
        
        # calculate force and strain of edges
        for i in range(self.ne):
            n0 = edges[i, 0]-1
            n1 = edges[i, 1]-1

            x0 = self.node[n0, 0]
            y0 = self.node[n0, 1]
            x1 = self.node[n1, 0]
            y1 = self.node[n1, 1]

            L = np.sqrt((x0-x1)**2+(y0-y1)**2)
            xl = (x1-x0)/L
            yl = (y1-y0)/L

            B = 1/L * np.array([-xl, -yl, xl, yl])

            # globalized stiffness matrix
            sctr=[2*n0, 2*n0+1, 2*n1, 2*n1+1]

            self.strain[i]=B.dot(self.d[sctr])
            self.stress[i]=self.E*self.strain[i]
            self.force[i]=self.A*self.stress[i]

        return self.force

    def compute_total_edges_length(self, edges=None):
        ans = 0
        if edges is None:
            edges=self.conn

        for e in edges:
            x1 = self.node[e[0]-1][0]
            y1 = self.node[e[0]-1][1]

            x2 = self.node[e[1]-1][0]
            y2 = self.node[e[1]-1][1]

            ans += np.sqrt((x1-x2)**2+(y1-y2)**2)

        return ans

if __name__=='__main__':
    fem_solver = FEM_Solver()
    print(fem_solver.get_displacement(fem_solver.conn))
    print(fem_solver.get_force(fem_solver.conn))

