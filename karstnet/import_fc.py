#    Copyright (C) 2018 by
#    Philippe Renard <philippe.renard@unine.ch>
#    Pauline Collon <pauline.collon@univ-lorraine.fr>
#    All rights reserved.
#    MIT license.
#
"""
Karstnet
========

Karstnet is a Python package for the analysis of karstic networks.

License
-------
Released under the MIT license:
   Copyright (C) 2018 Karstnet Developers
   Philippe Renard <philippe.renard@unine.ch>
   Pauline Collon <pauline.collon@univ-lorraine.fr>
"""

#----External librairies importations
import numpy as np
import networkx as nx
import scipy.stats as st
import matplotlib.pyplot as plt
import sqlite3
import mplstereonet

#----Internal module dependancies
from karstnet.base import *



# *************************************************************
# -------------------Public Functions:-------------------------
# -------------------GRAPH GENERATORS--------------------------
# *************************************************************

def from_nxGraph(nxGraph, coordinates, properties={}):
    """
    Creates a Karst graph from a Networkx graph.

    Takes a graph in the Networkx format and the coordinates of the
    nodes to generate a karstic network.

    Parameters
    ----------
    nxGraph : networkx graph
        the input graph

    coordinates : dictionnary
        the coordinates of the node, keys are node names

    properties : dictionnary
        optional argument containing properties associated with the nodes

    Returns
    -------
    KGraph
        A KGraph object

    Examples
    --------
       >>> myKGraph = kn.from_nxGraph(G, coord)
       >>> myKGraph = kn.from_nxGraph(G, coord, prop)
    """
    # Initialization of the complete graph
    edges = nx.to_edgelist(nxGraph)
    Kg = KGraph(edges, coordinates, properties)

    return Kg


def from_nodlink_dat(basename):
    """
    Creates the Kgraph from two ascii files (nodes, and links).

    Parameters
    ----------
    basename : string
        The base name used for the input files.

        The input files are named using the following convention:
         - basename_nodes.dat: the matrix of 3D node coordinates
           one line per node, one column per coordinates

         - basename_links.dat: the list of edges

    Returns
    -------
    KGraph
        A KGraph object

    Examples
    --------
       >>> myKGraph = kn.from_nodlink_dat("MyKarst")
    """

    link_name = basename + '_links.dat'
    node_name = basename + '_nodes.dat'

    # Read data files if exist - otherwise return empty graph
    try:
        links = np.loadtxt(link_name).astype(int) - 1
    except OSError:
        print("IMPORT ERROR: Could not import {}".format(link_name))
        return

    try:
        nodes = np.loadtxt(node_name)
    except OSError:
        print("IMPORT ERROR: Could not import {}".format(node_name))
        return
    # Create the dictionnary of coordinates
    coord = dict(enumerate(nodes[:, :3].tolist()))

    if len(nodes[0] > 3):
        properties = dict(enumerate(nodes[:, 3:].tolist()))
    else:
        properties = {}

    Kg = KGraph(links, coord, properties)
    print("Graph successfully created from file !\n")
    return Kg


# modif PVernant 2019/11/25
# add a function to read form an SQL export of Therion

def from_therion_sql(basename):
    """
    Creates the Kgraph from on SQL file exported from a Therion survey file.

    Parameters
    ----------
    basename : string
        The base name used for the input files.

        The input file is named using the following convention:
         - basename.sql: the containing all the needed informations

    Returns
    -------
    KGraph
        A KGraph object

    Examples
    --------
       >>> myKGraph = kn.from_therion_sql("MyKarst")
    """

    sql_name = basename + '.sql'

    # Read data files if exist - otherwise return empty graph
    try:
        conn = sqlite3.connect(':memory:')
        conn.executescript(open(sql_name).read())
#    	conn.executescript(open('../data/g_huttes.sql').read())
    except OSError:
        print("IMPORT ERROR: Could not import {}".format(sql_name))
        return

    # Read the SQL file and extract nodes and links data
    c = conn.cursor()
    c.execute('select st.ID, st.NAME, FULL_NAME, X, Y, Z from STATION\
    st left join SURVEY su on st.SURVEY_ID = su.ID;')
    nodes_th = []
    stations_th = []
    stations_id = []
    for s in c.fetchall():
        nodes_th.append([s[3], s[4], s[5]])
        stations_th.append(s[1])
        stations_id.append(s[0])

    c.execute('select FROM_ID, TO_ID from SHOT;')
    links_th = []
    for s in c.fetchall():
        links_th.append([s[0], s[1]])

    # Remove the splay links
    T = [s != '.' for s in stations_th]
    links_th = np.asarray(links_th).astype(int) - 1
    stations_id = np.asarray(stations_id).astype(int)[T] - 1
    links_ok = np.isin(links_th, stations_id)
    links = links_th[np.logical_and(links_ok[:, 0], links_ok[:, 1])]

    # Create the dictionnary of coordinates
    nodes = np.asarray(nodes_th)
    coord = dict(enumerate(nodes[:, :3].tolist()))

    if len(nodes[0] > 3):
        properties = dict(enumerate(nodes[:, 3:].tolist()))
    else:
        properties = {}

    Kg = KGraph(links, coord, properties)
    print("Graph successfully created from file !\n")
    return Kg

# end of modif PVernant 2019/11/25
# --------------------------------

def from_pline(filename):
    """
    Creates a KGraph from a Pline (Gocad ascii object)

    The function reads the Pline ASCII file and manages the colocated
    vertices indicated by the mention "ATOM" in the  file.
    This version loads also the Properties stored on vertices.

    Parameters
    ----------
    filename : string
        The name of the GOCAD Pline ASCII file.

    Returns
    -------
    KGraph
        A KGraph object

    Examples
    --------
       >>> myKGraph = kn.from_pline("MyKarst.pl")
    """

    # Read data files if exist - otherwise return empty graph
    try:
        # Open the ascii file in reading mode
        f_pline = open(filename, 'r')
    except OSError:
        print("IMPORT ERROR: Could not import {}".format(filename))
        return

    #  To store 3D location
    coord = {}
    # To store properties
    prop = {}
    # To store list of edges (each edge is a tuple)
    edges = []

    # Counter of nodes: in pl format, nodes are duplicated when changing
    # iline (eq. for branch). This is symbolized by the word ATOM instead of
    # VRTX and a segment uses the atom index instead those of the vrtx.
    # To track correspondance between VRTX and ATOM and avoids duplicates,
    # we use a counter of nodes, a dictionnary of nodes and one of atoms
    cpt_nodes = 0
    dico_nodes = {}  # make the correspondance betwen vrtx index and node index
    dico_atom = {}  # to memorize the atom index and use it to write segments

    for line in f_pline:
        if 'VRTX' in line:
            cpt_nodes += 1
            # cle,num,x,y,z=ligne.split()

            # because we do not pressupose the number of properties
            data = line.rstrip().split(" ")

            dico_nodes[int(data[1])] = cpt_nodes  # vrtx index vs. node index

            # store 3D location (relating to node index)
            coord[cpt_nodes] = (float(data[2]), float(data[3]), float(data[4]))
            # store properties if exist (relating to node index)
            prop[cpt_nodes] = dict(enumerate(list(np.float_(data[5:]))))
        if 'ATOM ' in line:
            cle, num, ref = line.split()
            # Atom must link to node index, not the index of the VTRX
            dico_atom[int(num)] = dico_nodes[int(ref)]
        if 'SEG' in line:
            cle, refi, refj = line.split()
            i = int(refi)
            j = int(refj)
            # Treatment of i:
            if i in dico_atom:
                # Replace atom number by the corresponding node index
                i = dico_atom[i]
            else:
                # Replace vertex number by the corresponding node index
                i = dico_nodes[i]
            # Treatment of j:
            if j in dico_atom:
                # Replace atom number by the corresponding node index
                j = dico_atom[j]
            else:
                # Replace vertex number by the corresponding node index
                j = dico_nodes[j]
            # Add the edge with the correct node indices
            edges.append((i, j))

    Kg = KGraph(edges, coord, prop)

    print("Graph successfully created from file !\n")
    f_pline.close()

    return Kg
