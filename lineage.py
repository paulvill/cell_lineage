# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt
import asciitable as asciitable
from scipy.spatial import Delaunay
import os
from scipy.sparse.csgraph import shortest_path

for embryo in [1,4,12,13,14,15,16,19]:
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
    
    # we want to produce a table with the list of all the cells, their initial and final time points
    # ind = [i for i, x in enumerate(names) if x == "P1"]
    
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
    
#    fig = plt.figure()
#    ax = plt.axes(projection='3d')
#    
#    #for t in range(int(max(cell_lineage[:,0]))): 
#    #int(max(cell_lineage[:,0]))-100
#    indices = np.where(cell_lineage[:,0] == 30)
#    xdata = cell_lineage[indices,1]
#    ydata = cell_lineage[indices,2]
#    zdata = cell_lineage[indices,3]
#        
#    ax.scatter3D(xdata, ydata, zdata, c = 'g', cmap='Greens')
    
        
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
        
            
    #delaunay triangulation
    #for t in np.arange(next(x[0] for x in enumerate(size) if x[1]>5),int(max(cell_lineage[:,0]))):
    #    indices = np.where(cell_lineage[:,0] == t)
    #    points = cell_lineage[indices[0],:]
    #    points = points[:,[1,2,3]]
    #    tri = Delaunay(points)
    #    fig = plt.figure(figsize = (16,12), dpi = 80, facecolor = 'w', edgecolor = 'k')
    #    ax = plt.axes(projection='3d')
    #    plot_tri_2(ax, points, tri, colors[indices[0]])
    #    plt.savefig('figures/'+emb_name+'/'+str(t)+'.png')
    #    plt.close(fig)
    
    #t = 50
    #indices = np.where(cell_lineage[:,0] == t)
    #points = cell_lineage[indices[0],:]
    #points = points[:,[1,2,3]]
    #tri = Delaunay(points)
    #fig = plt.figure()
    #ax = plt.axes(projection='3d')
    #plot_tri_2(ax, points, tri)
    #
    #edges = collect_edges(tri)
    #
    #adjacency = np.zeros((points.shape[0],points.shape[0]))
    #
    #for i in edges:
    #    adjacency[i[0],i[1]] = 1
    #    adjacency[i[1],i[0]] = 1
    #    
    #imshow(adjacency)
    
    #adjacency = np.zeros((len(names_unique),len(names_unique)))
    #
    #c = 0
    #for t in np.arange(next(x[0] for x in enumerate(size) if x[1]>5),int(max(cell_lineage[:,0]))):
    #    indices = np.where(cell_lineage[:,0] == t)
    #    points = cell_lineage[indices[0],:]
    #    points = points[:,[1,2,3]]
    #    tri = Delaunay(points)
    #    edges = collect_edges(tri)
    #    for i in edges:
    #        i1 = names_unique.index(names[indices[0][i[0]]])
    #        i2 = names_unique.index(names[indices[0][i[1]]])
    #        adjacency[i1,i2] = adjacency[i1,i2]+1
    #        adjacency[i2,i1] = adjacency[i2,i1]+1
    #    c = c+1
    #          
    #fig = plt.figure()
    #imshow(adjacency)
    #plt.savefig('figures/'+emb_name+'/averaged_adjacency.png')
    
    #c = 0
    ##os.mkdir('figures/'+emb_name+'/adjacencies/')
    #
    #for t in np.arange(next(x[0] for x in enumerate(size) if x[1]>5),int(max(cell_lineage[:,0]))):
    #    adjacency = np.zeros((len(names_unique),len(names_unique)))
    #    indices = np.where(cell_lineage[:,0] == t)
    #    points = cell_lineage[indices[0],:]
    #    points = points[:,[1,2,3]]
    #    tri = Delaunay(points)
    #    edges = collect_edges(tri)
    #    for i in edges:
    #        i1 = names_unique.index(names[indices[0][i[0]]])
    #        i2 = names_unique.index(names[indices[0][i[1]]])
    #        adjacency[i1,i2] = adjacency[i1,i2]+1
    #        adjacency[i2,i1] = adjacency[i2,i1]+1
    #    c = c+1
    #
    #    fig = plt.figure()
    #    imshow(adjacency)
    #    plt.savefig('figures/'+emb_name+'/adjacencies/'+str(t)+'.png')
    #    plt.close(fig)
    
    # proportion of cells in each clone
    prop = []
    for i in range(len(seeds)):
        prop.append(np.where(colors == i)[0].shape[0]/colors.shape[0])
        
    # towards a genealogical distance versus geometrical distance plot
    genealogical_distance_static = np.zeros((len(names_unique),len(names_unique)))
    for i in range(len(names_unique)):
        for j in range(i+1,len(names_unique)):
            n = names_unique[i]
            m = names_unique[j]
            if len(m)>len(n):
                m1 = m
                n1 = n
            else:
                m1 = n
                n1 = m    
            for o in range(len(n1)):
                if n1[:o+1] == m1[:o+1]:
                    genealogical_distance_static[i,j] = genealogical_distance_static[i,j]+0
                else:
                    genealogical_distance_static[i,j] = genealogical_distance_static[i,j]+2
            genealogical_distance_static[i,j] = genealogical_distance_static[i,j] + (len(m1) - len(n1))
    genealogical_distance_static = genealogical_distance_static+np.transpose(genealogical_distance_static)
               
    
    #from matplotlib.image import NonUniformImage
    #import matplotlib as mpl
    
    os.mkdir('figures/'+emb_name+'/genealogy_distance/')
    for t in np.arange(next(x[0] for x in enumerate(size) if x[1]>5),int(max(cell_lineage[:,0]))):    
        adjacency = np.zeros((len(names_unique),len(names_unique)))
        indices = np.where(cell_lineage[:,0] == t)
        points = cell_lineage[indices[0],:]
        points = points[:,[1,2,3]]
        tri = Delaunay(points)
        edges = collect_edges(tri)
        
        G = np.zeros((points.shape[0],points.shape[0]))
        
        for e in edges:
            G[e[0],e[1]] = 1
            G[e[1],e[0]] = 1
        
        dist_matrix = shortest_path(G, return_predecessors=False, directed=False, unweighted=True)
        X = []
        Y = []
        dist_matrix_gen = np.zeros((points.shape[0],points.shape[0]))
        for i in range(points.shape[0]):
            for j in range(i+1,points.shape[0]):
                dist_matrix_gen[i,j] = genealogical_distance_static[names_unique.index(names[indices[0][i]]),names_unique.index(names[indices[0][j]])]
                X.append(dist_matrix[i,j])
                Y.append(dist_matrix_gen[i,j])
        dist_matrix_gen = dist_matrix_gen + np.transpose(dist_matrix_gen)
        
        plt.figure(figsize=(7, 7))
        xedges = np.array(range(0,11))+0.5
        yedges = np.array(range(0,21))+0.5
        plt.hist2d(X,Y, bins=(xedges, yedges))
        axes = plt.gca()
        axes.set_xlim([0.5,10.5])
        axes.set_ylim([0.5,20.5])
        plt.savefig('figures/'+emb_name+'/genealogy_distance/'+str(t)+'.png')
        plt.close(fig)

    os.mkdir('figures/'+emb_name+'/genealogy_distance/all')
    X = []
    Y = []
    for t in np.arange(next(x[0] for x in enumerate(size) if x[1]>5),int(max(cell_lineage[:,0]))):    
        adjacency = np.zeros((len(names_unique),len(names_unique)))
        indices = np.where(cell_lineage[:,0] == t)
        points = cell_lineage[indices[0],:]
        points = points[:,[1,2,3]]
        tri = Delaunay(points)
        edges = collect_edges(tri)
        
        G = np.zeros((points.shape[0],points.shape[0]))
        
        for e in edges:
            G[e[0],e[1]] = 1
            G[e[1],e[0]] = 1
        
        dist_matrix = shortest_path(G, return_predecessors=False, directed=False, unweighted=True)

        dist_matrix_gen = np.zeros((points.shape[0],points.shape[0]))
        for i in range(points.shape[0]):
            for j in range(i+1,points.shape[0]):
                dist_matrix_gen[i,j] = genealogical_distance_static[names_unique.index(names[indices[0][i]]),names_unique.index(names[indices[0][j]])]
                X.append(dist_matrix[i,j])
                Y.append(dist_matrix_gen[i,j])
        dist_matrix_gen = dist_matrix_gen + np.transpose(dist_matrix_gen)
        
    plt.figure(figsize=(7, 7))
    xedges = np.array(range(0,11))+0.5
    yedges = np.array(range(0,21))+0.5
    plt.hist2d(X,Y, bins=(xedges, yedges))
    axes = plt.gca()
    axes.set_xlim([0.5,10.5])
    axes.set_ylim([0.5,20.5])
    plt.savefig('figures/'+emb_name+'/genealogy_distance/all/all.png')
    plt.close(fig)
        