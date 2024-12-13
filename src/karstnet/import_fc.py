#    Copyright (C) 2018-2024 by
#    Philippe Renard <philippe.renard@unine.ch>
#    Pauline Collon <pauline.collon@univ-lorraine.fr>
#    All rights reserved.
#    MIT license.
#
"""
Karstnet Import FC
==================

Karstnet is a Python package for the analysis of karstic networks.

The Import FC module contains functions for import and exports
in various formats.

"""

# ----External librairies importations
import numpy as np
import networkx as nx
import scipy.stats as st
import matplotlib.pyplot as plt
import sqlite3
import mplstereonet

# ----Internal module dependancies
from karstnet.base import *


# *************************************************************
# -------------------Public Functions:-------------------------
# -------------------GRAPH GENERATORS--------------------------
# *************************************************************

def from_nxGraph(nxGraph, coordinates, properties={}, verbose=True):
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

    Kg = KGraph(links, coord, properties, verbose=verbose)

    if verbose:
        print("Graph successfully created from file !\n")
    return Kg


# modif PVernant 2019/11/25
# add a function to read form an SQL export of Therion

def from_therion_sql(basename, verbose=True):
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
    T = [((s != '.') & (s != '-')) for s in stations_th]
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

    Kg = KGraph(links, coord, properties, verbose=verbose)

    if verbose:
        print("Graph successfully created from file !\n")
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





########################################

def from_therion_sql_enhanced(  inputfile, 
                                cavename=None, 
                                crs=None, 
                                rights=None,
                                citation=None ):
    """ This function 
    1. loads all the data from Therion sql, 
    2. add flags on shots and stations,
    3. regroupe nodes with same exact geographic coordinates,
    4. rename nodes,

    5. optionnaly remove nodes with srf, dpl, rmv, art, or spl flags 
    6. optionnaly return a networkx graph or karstnet Kgraph

    #potential Therion node flags 
    # 'ent' = entrance, 'con' = continuation, 'fix' = fixed, 
    # 'spr' = spring, 'sin' = sink, 'dol' = doline, 'dig' = dig, 
    # 'air' =air-draught, 'ove' = overhang, 'arc' = arch attributes  
    #potential Therion shot flags:
    # 'dpl' = duplicate, 'srf' = surface shots 

    Parameters
    ----------
    inputfile : string
        path to the SQL file exported with Therion

    cavename : string, optional
        Name of the cave. Will be attached to the graph as metadata. By default None

    crs : string, optional
        Coordinate reference. Will be attached to the graph as metadata. By default None

    rights : string, optional
        Information about the rights related to the dataset. Example: CC-BY-NC-SA, ... . Will be attached to the graph as metadata. By default None

    citation : string, optional
        Description on how to cite the dataset. Will be attached to the graph as metadata. By default None

    Returns
    -------
    G : networkx graph 
        with optional properties on nodes and edges:
        Dictionnaries always present on node: 'fulladdress', 'idsql'
        Optional dictionnaries on node: 'flags', 'pos', 'splays'
        Optional dictionnaries on edge: 'flags'

    Example:
    --------

    >>> from_therion_sql_enhanced('inputfilepath.sql')

    Metadata can be accessed with:
    >>> G.graph #not to confuse with the Kg.graph from Karstnet
    List of attribute names attached to the nodes:
    >>> set([k for n in G.nodes for k in G.nodes[n].keys()])
    List of attribute names attached to the nodes:
    >>> set([k for n in G.edges for k in G.edges[n].keys()])
    Dictionnaries can be accessed with:
    >>> nx.get_node_attribute(G,'attribute_name')


    """   

    import time 
    import networkx as nx
    import sqlite3
    from sqlite3 import OperationalError
    import sys

    def list2dict(key_list, value_list):
        """Transform list to dictionnary by regouping values in list for identical keys. 
        Using dictionnary comprehension.

        Parameters
        ----------
        key_list : list
            Dictionnary keys. Usually a list of int
        value_list : list
            Dictionnary values. Can be a list of int, flot, array, or list.

        Returns
        -------
        dictionnary

        """
        
        return {key : [value_list[idx] 
                for idx in range(len(value_list)) if key_list[idx]== key]
                for key in set(key_list)}

    def read_sql_file(basename):
        """
        Parameters
        ----------
        basename : str
            name of the sql database. without the extension.

        Returns
        -------
        c : TYPE
            DESCRIPTION.
        """
        # sql_name = basename #+ '.sql'

        try:
            conn = sqlite3.connect(':memory:')
            conn.executescript(open(basename,  encoding='utf-8-sig').read())
        #    	conn.executescript(open('../data/g_huttes.sql').read())
        except OSError:
            print("IMPORT ERROR: Could not import {}".format(basename))
        #    return

        # Read the SQL file
        c = conn.cursor()
        return c
    
    def extract_flags(c,type, return_type='dictionnary'):
        """  Extract flags attached to the stations or the shots to the form of a dictionnary.
        Shot Flags currently present in Therion:
        'dpl' = duplicate, 'srf' = surface shots

        Parameters
        ----------
        c : sqlite3.Cursor
            SQL database request cursor. this is the output of the function read_sql_file(basename).
        type : string
            To select wether its a shot flag or a station flag

        Returns
        -------
        dict
            Dictionnary with flags as keys, and each key containing a list of id or list of tuples.
        """ 

        flag_list=[]

        if type == 'shot':
            #extract shot flags with from-to info
            try:
                c.execute('select SHOT_FLAG.FLAG, SHOT.FROM_ID, SHOT.TO_ID from SHOT, SHOT_FLAG  \
                                            where SHOT.ID = SHOT_FLAG.SHOT_ID')
            except OperationalError:
                print(f'Cannot find sql. Verify that .sql exists or that the path is correct ')

            keys = []
            values = []
            id_from = []
            id_to = []
            for s in c.fetchall():
                #create list of tuple from all the links (from to)
                keys.append((s[1],s[2]))
                id_from.append(s[1])
                id_to.append(s[2])
                #create list with all the flags as value
                values.append(s[0])
            
            if return_type == 'lists':
                return id_from,id_to,values

            if return_type == 'dictionnary':
                return list2dict(keys, values) #dict(zip(keys, values))

        elif type == 'station':
            try:
                c.execute('select STATION_ID, FLAG from STATION_FLAG')
            except OperationalError:
                print(f'Cannot find sql. Verify that .sql exists or that the path is correct ')
            keys = []
            values = []
            for s in c.fetchall():
                #create list of station id
                keys.append(s[0])
                #create list with all the flags as value
                values.append(s[1])

            if return_type == 'lists':
                return keys,values

            if return_type == 'dictionnary':
                return list2dict(keys, values) #dict(zip(keys, values))


    #########################################
    ########################################


    #to check running time
    start_time = time.time()

    #read the sql database
    c = read_sql_file(inputfile)
    
    # import all LINKS 
    ###############
    print(f'Therion Import -- Importing all links (including splays) -- {(time.time() - start_time)}s')
    try:
        c.execute('select FROM_ID, TO_ID from SHOT')
    except OperationalError as e:
        print(f'1. Cannot find sql here: {inputfile}\n verify that .sql exists or that the path is correct ')
        raise e
        
    links_all = []
    for l in c.fetchall():
        links_all.append(l)
   
    
    # import NODES
    ###############################################################
    #import all nodes
    #prevents extraction of anonymous survey point symbol (- or .)  
    print(f'Therion Import -- Importing all nodes data (including splays) -- {(time.time() - start_time)}s')

    try:   
        c.execute('select st.ID, st.NAME, st.SURVEY_ID, FULL_NAME, X, Y, Z from STATION st \
                left join SURVEY su on st.SURVEY_ID = su.ID') 
    except OperationalError:
        print(f'2. Cannot find sql here: {inputfile}\n verify that .sql exists or that the path is correct ')
    
                           
    nodes_coord = [] # this is all the coordinates, including the splays
    nodes_id = [] # this all the ids, including the splays. (rename??)
    nodes_tree_structure = []
    for s in c.fetchall():
        #extract x,y,z nodes coordinates. this is all the coordinates, including splays
        nodes_coord.append([s[4], s[5], s[6]])
        #extract unique node id from Therion. this all the ids, including splays
        nodes_id.append(s[0])
        #extract full stree structure from Therion. this is all the tree-structure, including the splays
        if s[3]=='':
            nodes_tree_structure.append(f'{s[1]}')
        else:
            address = '.'.join(s[3].split('.')[::-1])            
            nodes_tree_structure.append(f'{address}.{s[1]}')
        # nodes_tree_structure.append('%s@%s'%(s[1],s[3]))
    #create dictionnary of the nodes coordinates
    coord = dict(zip(nodes_id,nodes_coord))
    #save tree structure in the form of two list, to prevent data loss when combining nodes
    #only take the stations
    list_tree_oldi = []
    list_tree_values = []
    for i, tree in enumerate(nodes_tree_structure):
        if tree.startswith(('.','-'))==False:#tree[0].isdigit():  
            list_tree_values.append(nodes_tree_structure[i])
            list_tree_oldi.append(nodes_id[i])
        #if tree.startwith(-) or tree.startwith(.):


    
    #create graph with all the links
    #################################
    print(f'Therion Import -- Create initial graph with all the data points (including splays) -- {(time.time() - start_time)}s')
    G = nx.Graph(cavename=cavename, crs=crs, original_data_rights=rights, citation=citation)
    G.add_edges_from(links_all)
    # if the nodes attributes are the same for two combined nodes, it seems that it does not affect the combining
    nx.set_node_attributes(G, coord, 'pos')
    # nx.set_node_attributes(G, tree_structure, 'tree_structure')

    
    # Import splay leg 
    ###################################
    ##################################################################
    #remove nodes that are anonymous survey point symbol (- or .)
    
    splay_id = []  #this is the sql id of the splay itself
    splay_coord = []
    try:
        c.execute('select st.ID, st.NAME, X, Y, Z from STATION st \
                where st.NAME in (".","-") or st.NAME like "%splay%"' )
    except OperationalError:
        print(f'3. Cannot find sql here: {inputfile}\n verify that .sql exists or that the path is correct ')
    
    for s in c.fetchall():
        #extract x,y,z nodes coordinates
        splay_coord.append([s[2], s[3], s[4]])
        #extract unique node id from Therion
        splay_id.append(s[0])

    #!!! remove splays from the nodes. they will be imported later
    if splay_id:
        G.remove_nodes_from(splay_id)
    else: 
        print('no splays legs to remove')
      
    # Import splay leg shot info on nodes in the form of a list of coordinates of the end of the shot. 
    #create dictionnary of the nodes coordinates for each splays. 
    #the dictionnary key corresponds to the id for each splay in the sql database
    coord = dict(zip(splay_id,splay_coord))

    #import links only for the nodes we exported
    string_id = ",".join(map(str,splay_id))
    try:
        c.execute('select FROM_ID, TO_ID from SHOT \
                where TO_ID in (%s)' % (string_id))
    except OperationalError:
        print(f'4. Cannot find sql here: {inputfile}\n verify that .sql exists or that the path is correct ')
    links = []
    for l in c.fetchall():
        links.append([l[0], l[1]])

    #replace splay node id with the station id to which the splay is shot from
    #for example, if 2,3,4 are splay id, and attached to station 1, then all the id will be 1
    #splays_dict = defaultdict(list)
    #make two lists of splays 1. station of departure, 2. coordinates for arrival
    #(make a drawing to explain this)
    list_splays_oldi=[]
    list_splays_pos=[]
    for link in links:   
        if link:
            list_splays_oldi.append(link[0])
            list_splays_pos.append(coord[link[1]])
            #splays_dict[link[0]].extend([coord[link[1]]]) 
    
    #nx.set_node_attributes(G, splays_dict, 'splays')   
    
   
    #COMBINE IDENTIDAL STATIONS
    #Rename nodes and get ride of duplicate nodes with identical position
    ##############################################################################
    ##############################################################################
    # this rename nodes with identical position with the same id, 
    # which automatically regroup the nodes with identical name into one.
     
    #pos2d = {key: value[0:2] for key, value in nx.get_node_attributes(G,'coord').items()}
    # plt.figure()
    # nx.draw(G,pos=pos2d)
    
    #find nodes with duplicate positions:
    #create a list of lists of index where the coordinates are the same
    print(f'Therion Import -- Combine Stations with identical x,y,z -- {(time.time() - start_time)}s')
    unique_pos = [list(x) for x in set(tuple(x) for x in list(nx.get_node_attributes(G,'pos').values()))]
    # print(len(unique_pos))
    duplicates = []
    for i,position in enumerate(unique_pos):
        if i%1000 == 0:
            print(f'{i}/{len(unique_pos)} unique positions')
        #this could be sped up by inversing the dictionnary key and values??
        duplicates.append([key for key,coord in G.nodes('pos') if coord==position])
        #duplicates_fulladdress.append([])


        
    #rename nodes 
    ########################################################################
    #duplicate nodes are renamed with the same name
    #create new ids dictionnary to replace the initial indexes  
    #create new ids with repeating values for idential node position 
    #################################################################
    newis = []  
    print(f'Therion Import -- Rename nodes -- {(time.time() - start_time)}s') 
    # print(f'len(duplicates) = {len(duplicates)}')
    for i, index in enumerate(range(len(duplicates))):
        if i%1000 == 0:
            print(f'{i}/{len(duplicates)} nodes to rename')
        newis = newis + [index]*len(duplicates[i])
        
    #flatten the list of list of old ids
    #####################################
    print(f'Therion Import -- concatenate old ic in a dictionnary -- {(time.time() - start_time)}s')
    concat_oldi = [j for i in duplicates for j in i]   
    #the dictionnary has to be in the form of dict keys are the old keys, and the value is the new key
    index_dict = dict(zip(concat_oldi, newis ))
    # index_fulladdress = dict(zip(concat_fulladdress,newis))

    # #extract full tree info
    # #########################
    # concat_fulladdress = []
    # for index in concat_oldi:
    #     concat_fulladdress.append(G.nodes('tree_structure')[index])
    
    print(f'Therion Import -- Relabel nodes -- {(time.time() - start_time)}s')
    #rename nodes (nodes with same geographic posiion will be "merged" under the same name)
    G = nx.relabel_nodes(G,index_dict)
    #drop edges that link the node to themselves. happen because of the combining the nodes.
    print(f'Therion Import -- remove self links -- {(time.time() - start_time)}s') 
    G.remove_edges_from(list(nx.selfloop_edges(G)))


    #Add attributes to the graph with the new ids,
    ################################################################
    print(f'Therion Import --add dictionnaries to graph -- {(time.time() - start_time)}s')  
    #combines the information for nodes that are regrouped
    #this steps has to be mande after the nodes have been regrouped, otherwise, 
    #the networkx function just gets rid of attribute values is they exist in two or more combined nodes

    #SPLAYS
    #######
    print(f'Therion Import -- add splays -- {(time.time() - start_time)}s') 
    list_splays_newi = [index_dict.get(item, item)  for item in list_splays_oldi]
    dict_splays = list2dict(list_splays_newi, list_splays_pos)
    nx.set_node_attributes(G, dict_splays, 'splays') 

    #TREE
    #####
    print(f'Therion Import -- add fulladdress -- {(time.time() - start_time)}s')
    list_tree_newi = [index_dict.get(item, item)  for item in list_tree_oldi]
    dict_tree = list2dict(list_tree_newi, list_tree_values)
    nx.set_node_attributes(G, dict_tree, 'fulladdress') 


    #add potential node flags
    ###########################
     # 'ent' = entrance, 'con' = continuation, 'fix' = fixed, 
     # 'spr' = spring, 'sin' = sink, 'dol' = doline, 'dig' = dig, 
     # 'air' =air-draught, 'ove' = overhang, 'arc' = arch attributes
    #load the flags with the sql index
    print('Therion Import -- add flags') 
    list_node_flag_oldi, list_node_flag_values = extract_flags(c,'station', return_type='lists')
    list_node_flag_newi = [index_dict.get(item, item)  for item in list_node_flag_oldi]
    dict_node_flag = list2dict(list_node_flag_newi, list_node_flag_values)
    nx.set_node_attributes(G, dict_node_flag, 'flag') 
    
    #add potential edge flags
    ############################
    # Shot Flags
    # 'dpl' = duplicate, 'srf' = surface shots
    from_edge_flag_oldi, to_edge_flag_oldi, list_edge_flag_values = extract_flags(c,'shot', return_type='lists')   
    list_from_edge_flag_newi = [index_dict.get(item, item)  for item in from_edge_flag_oldi]
    list_to_edge_flag_newi = [index_dict.get(item, item)  for item in to_edge_flag_oldi]
    dict_edge_flag = list2dict(list(zip(list_from_edge_flag_newi,list_to_edge_flag_newi)), list_edge_flag_values)
    nx.set_edge_attributes(G, dict_edge_flag, 'flags') 
    

    #SQL IDs (oldi)
    #add old therion id name as a property
    ################  
    #has to be reversed from the oldi-newi dictionnary, 
    #but preserving the 
    print(f'Therion Import -- add sql ids -- {(time.time() - start_time)}s')
    sql_ids = {}
    for k, v in zip(newis, concat_oldi):
        sql_ids.setdefault(k, []).append(v)
    nx.set_node_attributes(G, sql_ids, 'idsql')
    
    #remove nodes that were isolated when removing the edges
    #not sure that this is still necessary
    print(f'Therion Import -- remove isolated nodes -- {(time.time() - start_time)}s')  
    G.remove_nodes_from(list(nx.isolates(G)))      

    # #remove unnecessary attributes
    # for (n,d) in G.nodes(data=True):
    #     del d["tree_structure"]




    return G

    # #can return either the Kgraph object or just the graph in the networkx format
    # if export_Kgraph == True:
    #     import karstnet as kn
    #     return kn.from_nxGraph(G, dict(G.nodes('pos')), properties=None, verbose=True)
    # else:
    # return H




        ##################
    # ADD MANUAL DATA
    ##################

def clean_therion_sql_graph(G,
                            flagged_edges_to_remove=['srf','dpl','rmv','art','spl'],
                            additional_edges = None, 
                            additional_flags_edges=None,
                            additional_flags_nodes=None):
    """_summary_

    Parameters
    ----------
    G : networkx graph 
        Graph produced with the function kn.from_therion_sql_enhanced


    flagged_edges_to_remove : list of strings, optional
        list of the flags for which edges should be removed, by default ['srf','dpl','rmv','art','spl']
        - 'dpl' : duplicate
        - 'srf' : surface
        - 'art' : artificial
        - 'rmv' : remove
        - 'spl' : splay (for example when a shot is made in a large room, star shots, ...)

    additional_edges : list of tuple or list of lists, optional
        list of edges to add to the graph. The edges will be created between already existing nodes, and a flag add will be added. 
        The naming convention for the node should be based on the full address of the point, by default None
        additional_edges = [['full_address.0','full_address.1],['full_address.3','full_address.10]]

    additional_flags_edges : dictionnary of list of tuple or list of lists, optional
        lists of edges to be flagged with corresponding flag to add. Any flag can be added. by default None
        However, for the duplicate edges and surface edges to be removed, it is necessary to use the correct flags.
        It is also possible to use the add edges with this dictionnary instead of using the 'add_edges' option.
        List of flagged edge that will be removed by default:  'dpl', 'srf', 'art', 'rmv', 'spl'

        Example of dictionnary: 
        ----------------------
        add_flags_edges = {'dpl':[['full_address.0','full_address.1],['full_address.3','full_address.10]], 'srf':[[full_address.3,full_address.2]]}   

    additional_flags_nodes : dictionnary of list, optional
        list of nodes to be flagged, by default None

        Example of dictionnary: 
        -----------------------
        additional_flags_nodes = {'ent':['full_address.0','full_address.1', 'full_address.4' ]


    Returns
    -------
    networkx graph
        Clean networkx graph
    """    ''''''

    import networkx as nx
    
    dict_address = dict(G.nodes('fulladdress'))
    inverse_dict_address = { v: k for k, l in dict_address.items() for v in l }

    # only adds missing edges. this is just a simpler option than using the dictionnaries. Not sure it is usefull..
    if additional_edges is not None:
        for i,edge in enumerate(additional_edges):
            edges[i][0] = [key for key, value in dict_address.items() if edge[0] in value ][0]
            edges[i][1] = [key for key, value in dict_address.items() if edge[1] in value ][0]            
            #check if there is flags already attached to the edges
            #if not, create a new edge with the flag value
            G.add_edge(edges[i][0],edges[i][1],flags=['add'])
    

    # add flags on nodes
    if additional_flags_nodes is not None: 
        print(f'Therion Import - adding manual node flags: {additional_flags_nodes.keys()}')
         
        for flag in additional_flags_nodes.keys():
            #loop throught the flags
            for fulladdress in additional_flags_nodes[flag]:
                id_node = inverse_dict_address[fulladdress]                      
                #check if node already has a flag
                #if not, create a new list of flag(s) attached to the node
                if G.nodes('flags')[id_node] is None:
                    # pass
                    #create flag on node
                    nx.set_node_attributes(G, {id_node:[flag]}, name='flags')
                #if yes, append the new flag to the list
                elif G.nodes('flags')[id_node] is not None:
                    #print(flag,id_node, fulladdress)
                    G.nodes[id_node]['flags'].append(flag)
    else:
        print('No manual node flags to add')

    if additional_flags_edges is not None: 
        print(f'Therion Import - adding manual edges flags: {additional_flags_edges.keys()}')
        #loop through all the flags
        for flag in additional_flags_edges.keys():
            edges = additional_flags_edges[flag]
            #loop through all the edges for each flag
            for i,addresses in enumerate(edges):
                edge_1 = inverse_dict_address[addresses[0]]
                edge_2 = inverse_dict_address[addresses[1]]
                #check if edge exists already (for example to add a duplicate flag on an exisiting edge)
                if G.has_edge(edge_1,edge_2):
                #check if there is flags already attached to the edges
                #if not, create a new edge with the dictionnary 'flags' and the flag value
                    if 'flags' not in G[edge_1][edge_2]:
                        nx.set_node_attributes(G,{'flags':[flag]})
                        
                    #if yes, append the flag to the list
                    elif 'flags' in G[edge_1][edge_2]:
                        G[edge_1][edge_2]['flags'].append(flag)
                #if edge does not exist yet, then just create a new edge with the appropriate flag name
                else:
                    G.add_edge(edge_1,edge_2,flags=[flag])

    else:
        print('Therion Import - No manual edge flags to add')

    if flagged_edges_to_remove != []:
        print(f'Therion Import - removing flagged_edges: {flagged_edges_to_remove}')



        H= G.copy()
        #edges_to_remove = list(dict(nx.get_edge_attributes(G,'flags')).keys()) 
        #extract flag unique values into a list
        flags = {x for l in list(nx.get_edge_attributes(G,'flags').values()) for x in l}
        #loop through the unique flags and 
        for flag in flags:
            #only remove edges with flag surface, duplicate, or remove
            if flag in flagged_edges_to_remove:
                
                list_edges = [edge for edge, action in nx.get_edge_attributes(G,'flags').items() if flag in action]
                #print('remove ', flag, list_edges)
                #remove edges    
                H.remove_edges_from(list_edges)
                #remove nodes that were isolated when removing the edges
                H.remove_nodes_from(list(nx.isolates(H)))
            else:
                pass
                #print('not removed', flag)


        print(f'Initial Graph size: {len(G)}, Graph size after removing flagged edges: {len(H)}')
        return H
    else:
        print("clean_therion_sql_graph - to remove flagged edges: flagged_edges_to_remove=['srf','dpl','rmv','art','spl']")