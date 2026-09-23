"""
Module `_io_func.pline`
-----------------------

The module `_io_func.pline` is accessible as `io_func`. It contains functions 
to import / export networkx graphs from / to Pline (Gocad ascii object).
"""

import numpy as np
import networkx as nx


# ===== Import functions based on Pline (Gocad ascii) ========================

# ----------------------------------------------------------------------------
def networkx_from_pline(
        filename, 
        pos_attr='pos',
        verbose=True):
    """
    Creates a networkx graph from a Pline (Gocad ascii object).

    The function reads the Pline ASCII file and manages the colocated
    vertices indicated by the mention "ATOM" in the  file.
    This version loads also the Properties stored on vertices.

    Parameters
    ----------
    filename : str
        the name of the GOCAD Pline ASCII file

    pos_attr : str, default: 'pos'
        name of the node attribute for position

    verbose : bool, default: True
        indicates if info is printed

    Returns
    -------
    G : networkx.Graph
        networkx graph

    Examples
    --------
        >>> G = kn.io.networkx_from_pline("MyKarst.pl")
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
    # we use a counter of nodes, a dictionary of nodes and one of atoms
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
            prop[cpt_nodes] = dict(enumerate(list(np.float64(data[5:]))))
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

    if len(prop):
        node_property = prop
        node_property_name = 'dict_of_prop'
    else:
        node_property = None
        node_property_name = None


    # Set output networkx graph
    G = nx.Graph()
    
    G.add_nodes_from(list(coord.keys()))
    
    nx.set_node_attributes(G, coord, pos_attr)
    
    if node_property is not None:
        nx.set_node_attributes(G, node_property, node_property_name)
    
    G.add_edges_from(edges)

    if verbose:
        print("Networkx graph successfully created from (pline) file !\n")

    f_pline.close()

    return G
# ----------------------------------------------------------------------------
