#  !/usr/bin/env python
#  -*- coding:utf-8 -*-

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

neunum_coef = 0.1
PN_num = int(90*neunum_coef)
LN_num = int(30*neunum_coef)
PN_stim_num = int(30*neunum_coef)
LN_stim_num = int(10*neunum_coef)

fig_name_list = \
    ['LN2PN_slow', 'LN2PN_GABA', 'LN2LN_GABA', 'PN2PN_nACH', 'PN2LN_nACH']
#          0             1             2             3             4
fig_name = fig_name_list[4]


def creat_graph(gName):
    # input: data file name, the main part
    # output: the graph build on data file
    # function:
    #     read coupling from txt data files, convert it into graphes,
    #     plot the graph, save a hard copy of the graph, and
    #     return the graph (as a networkx object)
    # code: ------
    # read from txt file
    gMat = np.loadtxt('mat_'+gName+'.txt')
    start, end = np.shape(gMat)  # start neuron & end neuron of each edge.
    grf = nx.MultiDiGraph()
    # add edges
    for i in range(start):
        for j in range(end):
            if gMat[i][j] != 0:
                grf.add_edge(gName[0:2]+'%2d' % i, gName[3:5]+'%2d' % j)
    # draw the graph and save a hard copy
    gPos = nx.circular_layout(grf)  # positions for all nodes
    nx.draw_networkx_nodes(grf, gPos, node_size=1000, node_color='g')
    nx.draw_networkx_edges(grf, gPos, width=1, edge_color='b')
    nx.draw_networkx_labels(grf, gPos, font_size=10, font_family='sans-serif')
    plt.axis('off')
    plt.savefig(gName+'.jpg')
    # return the graph
    return grf

if __name__ == '__main__':
    print fig_name
    creat_graph(fig_name)
