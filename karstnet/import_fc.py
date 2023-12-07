#    Copyright (C) 2018-2023 by
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
   Copyright (C) 2018-2023 Karstnet Developers
   Philippe Renard <philippe.renard@unine.ch>
   Pauline Collon <pauline.collon@univ-lorraine.fr>
"""

# ----External libraries importations
import numpy as np
import networkx as nx
# import scipy.stats as st
# import matplotlib.pyplot as plt
import sqlite3

# ----Internal module dependent
from karstnet.base import KGraph
import karstnet as kn


# *************************************************************
# -------------------Public Functions:-------------------------
# -------------------GRAPH GENERATORS--------------------------
# *************************************************************

def from_nxGraph(nxGraph, coordinates, properties=None, verbose=True):
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
       >>> G = KGraph([],{})
       >>> myKGraph = kn.from_nxGraph(G, {})
       >>> myKGraph_with_prop = kn.from_nxGraph(G, {}, {})
    """
    # Initialization of the complete graph
    if properties is None:
        properties = dict()
    else:
        properties = properties
    edges = nx.to_edgelist(nxGraph)
    Kg = KGraph(edges, coordinates, properties, verbose=verbose)

    return Kg


def from_nodlink_dat(basename, verbose=True):
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

    # Determinate the number of columns of the intial sheet
    data_node = open(node_name, 'r')
    first_line = data_node.readline()
    param_number = len(first_line.split())
    convert = dict()
    for i in range(3):
        convert[i] = _check_type_coord
    for i in range(3, param_number):
        convert[i] = _check_type_prop

    # Read data files if exist - otherwise return empty graph
    try:
        links = np.loadtxt(link_name).astype(int) - 1
    except OSError:
        print("IMPORT ERROR: Could not import {}".format(link_name))
        return

    try:
        nodes = np.loadtxt(node_name, converters=convert)
    except OSError:
        print("IMPORT ERROR: Could not import {}".format(node_name))
        return
    # Create the dictionary of coordinates
    coord = dict(enumerate(nodes[:, :3].tolist()))

    if len(nodes[0] > 3):
        properties = dict(enumerate(nodes[:, 3:].tolist()))
    else:
        properties = {}

    Kg = KGraph(links, coord, properties, verbose=verbose)

    if verbose:
        print("Graph successfully created from file !\n")
    return Kg


# modif PVernant 2019/11/25
# add a function to read form an SQL export of Therion

def from_therion_sql(basename, verbose=True, properties=None):
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
    #    return

    # Read the SQL file
    c = conn.cursor()


    # load flags
    ############
    # flags will be used to ignore duplicate and surface points
    # duplicates are integreted in therion coordinates calculations,
    # but they should not be loaded for the graph.
    # this is the list of graphs 
    # 'dpl' = duplicate, 'srf' = surface shots
    # 'ent' = entrance, 'con' = continuation, 'fix' = fixed, 
    # 'spr' = spring, 'sin' = sink, 'dol' = doline, 'dig' = dig, 
    # 'air' =air-draught, 'ove' = overhang, 'arc' = arch attributes

    shot_flags_id = []
    c.execute('select SHOT_ID, FLAG from SHOT_FLAG where FLAG="srf" or FLAG="dpl"')
    for fl in c.fetchall():
        shot_flags_id.append(fl[0])

    station_flags_id = []
    c.execute('select STATION_ID, FLAG from STATION_FLAG where FLAG="srf" or FLAG="dpl"')
    for fl in c.fetchall():
        station_flags_id.append(fl[0])


    #import nodes without the splays
    #prevents extraction of anonymous survey point symbol (- or .) 
    #prevents extraction of surface and duplicate points (shot_flags_id)
    string_station_flags_id = ",".join(map(str,station_flags_id))
    c.execute('select st.ID, st.NAME, FULL_NAME, X, Y, Z from STATION \
                st left join SURVEY su on st.SURVEY_ID = su.ID \
                where st.NAME not in (".") and st.NAME not in ("-") and st.ID not in (%s)' % string_station_flags_id)
    #            where st.NAME not in (".","-") and st.ID not in (%s)' % string_station_flags_id)

    nodes_coord = []
    stations_id = []
    for s in c.fetchall():
        #extract x,y,z nodes coordinates
        nodes_coord.append([s[3], s[4], s[5]])
        #extract unique node id from Therion
        stations_id.append(s[0])
    #create dictionnary of the nodes coordinates
    coord = dict(zip(stations_id,nodes_coord))

    string_shot_flags_id = ",".join(map(str,shot_flags_id))
    #import links only for the nodes we exported
    string_id = ",".join(map(str,stations_id))
    c.execute('select FROM_ID, TO_ID from SHOT \
            where FROM_ID not in (%s) and FROM_ID in (%s) and TO_ID in (%s)' % (string_shot_flags_id,string_id,string_id))

    links = []
    for l in c.fetchall():
        links.append([l[0], l[1]])
    links = np.asarray(links).astype(int)

    Kg = kn.KGraph(links, coord, verbose=verbose, properties=properties)

    if verbose:
        print("Graph successfully created from file !\n")
#    return Kg

    return Kg

# end of modif PVernant 2019/11/25
# --------------------------------


def from_pline(filename, verbose=True):
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

    Kg = KGraph(edges, coord, prop, verbose=verbose)

    if verbose:
        print("Graph successfully created from file !\n")
    f_pline.close()

    return Kg

# **************************************************************
#
# -----------LH2022--NON public functions used by import_fc-----
#
# **************************************************************


def _check_type_prop(s):
    """
    Check if the value input is a positive float.

    Parameters
    ----------
        s : every type
            An object to check

    Returns:
    --------
        The object as float if it is a positive float initially
        As no definite object NAN if the input object is not a float,
        or a negative float
    """

    try:
        s = float(s)
        if s > 0.:
            return s
        return np.nan
    except ValueError:
        return np.nan


def _check_type_coord(s):
    """
    Check if the value input is a positive float.

    Parameters
    ----------
        s : every type
            An object to check

    Returns:
    --------
        The object as float if it is a postiv float intially
        As no definite object NAN if the input object is not a float,
        or a negativ float
    """

    try:
        return float(s)
    except ValueError:
        return np.nan

# **************************************************************
#
# -----------End Louise HOUBRE 2022----
#
# **************************************************************
