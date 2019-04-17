#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 15:16:29 2019

@author: paul
"""


import numpy as np
import matplotlib.pyplot as plt
import asciitable as asciitable
from scipy.spatial import Delaunay
import os
from scipy.sparse.csgraph import shortest_path
from mpl_toolkits.mplot3d import Axes3D

embryo = 4

if embryo<10:
    emb_name = "WT-EMB0"+str(embryo)
else:
    emb_name = "WT-EMB"+str(embryo)
file_name = "/Users/paul/Google Drive/New Lab/Ongoing Projects/cell packing/code/data/Galaxy35-Cell-lineage-and-position-WT-compressed/"+emb_name+".txt"
cell_lineage = np.loadtxt(fname = file_name,skiprows = 1,usecols = (0,2,3,4))

names_ = asciitable.read(file_name)

names = []
founder_cells = ['AB','C','D','E','EMS','MS','P1','P2','P3','P4','Z2','Z3']
seeds = ['AB','C','D','E','EMS','MS','P1','P2','P3','P4','Z2','Z3','ABal','ABar','ABpl','ABpr']
colors = np.zeros((len(names_),1))
positions = ['a','p','r','l','v','d']

for i in range(len(names_)):
    names.append(names_[i][1])
    n = names_[i][1]
    try:
        colors[i] = seeds.index(n[0:4])
    except ValueError:
        try:        
            colors[i] = seeds.index(n[0:2])
        except ValueError:
            colors[i] = seeds.index(n[0:1])
    
names_unique = list(set(names))
names_unique.sort()

all_cells = np.zeros((len(names_unique),4))

c = 0
for n in names_unique:
    ind = [i for i, x in enumerate(names) if x ==  n]
    all_cells[c,0] = c
    all_cells[c,1] = cell_lineage[ind[0],0]
    all_cells[c,2] = cell_lineage[ind[-1],0]
    if not n in founder_cells :
        all_cells[c,3] = names_unique.index(n[:-1])
    
    c = c+1    

    
fig0 = plt.figure()    
a, bin_edges = np.histogram(all_cells[:,2] - all_cells[:,1])
plt.hist(a, bins='auto')  
plt.show

size = []    

for t in range(int(max(cell_lineage[:,0]))):
    temp = np.where(cell_lineage[:,0] == t)
    temp1 = temp[0].shape
    size.append(temp1[0])
    
fig1 = plt.figure()
plt.plot(range(int(max(cell_lineage[:,0]))),size)
plt.ylabel('Number of Cells')
plt.xlabel('Time Steps')
plt.show


def plot_tri_2(ax, points, tri, color):
    edges = collect_edges(tri)
    x = np.array([])
    y = np.array([])
    z = np.array([])
    sizes = 100*np.ones((1,points.shape[0]))
    for (i,j) in edges:
        x = np.append(x, [points[i, 0], points[j, 0], np.nan])      
        y = np.append(y, [points[i, 1], points[j, 1], np.nan])      
        z = np.append(z, [points[i, 2], points[j, 2], np.nan])
    ax.plot3D(x, y, z, color='g', lw='0.1')

    ax.scatter(points[:,0], points[:,1], points[:,2],s = sizes, c=10*np.transpose(color)[0])

def collect_edges(tri):
    edges = set()

    def sorted_tuple(a,b):
        return (a,b) if a < b else (b,a)
    # Add edges of tetrahedron (sorted so we don't add an edge twice, even if it comes in reverse order).
    for (i0, i1, i2, i3) in tri.simplices:
        edges.add(sorted_tuple(i0,i1))
        edges.add(sorted_tuple(i0,i2))
        edges.add(sorted_tuple(i0,i3))
        edges.add(sorted_tuple(i1,i2))
        edges.add(sorted_tuple(i1,i3))
        edges.add(sorted_tuple(i2,i3))
    return edges


#os.mkdir('figures/'+emb_name+'/')
    
z_max = max(cell_lineage[:,3])
x_min = min(cell_lineage[:,1])
x_max = max(cell_lineage[:,1])
for t in np.arange(next(x[0] for x in enumerate(size) if x[1]>5),int(max(cell_lineage[:,0]))):
    indices = np.where(cell_lineage[:,0] == t)
    points = cell_lineage[indices[0],:]
    points = points[:,[1,2,3]]
    tri = Delaunay(points)
    fig = plt.figure(figsize = (16,12), dpi = 80, facecolor = 'w', edgecolor = 'k')
    ax = plt.axes(projection='3d')
    plot_tri_2(ax, points, tri, colors[indices[0]])        
#    plt.savefig('figures/'+emb_name+'/3D'+str(t)+'.png')
#    plt.close(fig)
    
    nb_sh = 10
    
    points_2D = points[:,[1,2]]
    
    for i in range(nb_sh):
        x_bound_min = x_min+i*(x_max-x_min)/nb_sh
        x_bound_max = x_min+(i+1)*(x_max-x_min)/nb_sh
        for j in range(len(points[:,0])):
            if (points[j,0] >= x_bound_min) & (points[j,0] < x_bound_max):
                if i%2:
                    points_2D[j,1] = (i+1)*z_max - points[j,2]
                else:
                    points_2D[j,1] = (i)*z_max + points[j,2]
     
    edges = collect_edges(tri)     
    
    fig = plt.figure(figsize = (16,12), dpi = 80, facecolor = 'w', edgecolor = 'k')
    
    for (i,j) in edges:
        plt.plot([points_2D[i,0],points_2D[j,0]],[points_2D[i,1],points_2D[j,1]],color = 'g',lw='0.1')
    sizes = 100*np.ones((1,points.shape[0]))
    plt.scatter(points_2D[:,0],points_2D[:,1],s = sizes,c = 10*np.transpose(colors[indices[0]])[0])     
#    plt.savefig('figures/'+emb_name+'/flat'+str(t)+'.png')
#    plt.close(fig)
    
    
    