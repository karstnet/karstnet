"""
Module `_io_func.therion`
-------------------------

The module `_io_func.therion` is accessible as `io_func`. It contains functions
to import / export networkx graphs from / to therion format.
"""

import numpy as np
import networkx as nx
import sqlite3 # for import therion files
from sqlite3 import OperationalError
import karstnet as kn
import time 



def networkx_from_therion_sql(basename,
                              pos_attr='pos',
                              remove_flagged_edges=False,
                              verbose=True):
    """ This function loads all the data from a Therion SQL file,
    This function replace the initial function developped by PVernant, but fixes small issues with 
    node ids if they are not consecutives, and enables to remove links flagges as duplicate or surface.

    The function:
    - add flags on shots and stations,
    - regroupe nodes with same exact geographic coordinates,
    - rename nodes,
    - optionnaly remove nodes with srf, dpl, rmv, art, or spl flags ?? this is not true anymore. should I add this option??

    #potential Therion node flags 
    # 'ent' = entrance, 'con' = continuation, 'fix' = fixed, 
    # 'spr' = spring, 'sin' = sink, 'dol' = doline, 'dig' = dig, 
    # 'air' =air-draught, 'ove' = overhang, 'arc' = arch attributes  
    #potential Therion shot flags:
    # 'dpl' = duplicate, 'srf' = surface shots 

    Parameters
    ----------
    basename : string
        path to the SQL file exported with Therion

    pos_attr : str, default: 'pos'
        name of the node attribute for position

    remove_flagged_edges : bool
        If set to True, all the edges flagged with 'dpl' or 'srf' are removed from the dataset
        by default: False

    verbose : bool, default: True
        indicates if info is printed

    Returns
    -------
    G : networkx graph 
        with optional properties on nodes and edges:
        Dictionnaries always present on node: 'fulladdress', 'idsql'
        Optional dictionnaries on node: 'flags', 'pos', 'splays'
        Optional dictionnaries on edge: 'flags'

    Example:
    --------

    >>> networkx_from_therion_sql('cavename.sql')

    Metadata can be accessed with:
    >>> G.graph #not to confuse with the Kg.graph from Karstnet
    List of attribute names attached to the nodes:
    >>> set([k for n in G.nodes for k in G.nodes[n].keys()])
    List of attribute names attached to the nodes:
    >>> set([k for n in G.edges for k in G.edges[n].keys()])
    Dictionnaries can be accessed with:
    >>> nx.get_node_attribute(G,'attribute_name')


    """   




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
        if basename.endswith('.sql'):
            pass
        else:
            basename = basename+ '.sql'

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
                return kn.list2dict(keys, values) #dict(zip(keys, values))

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
                return kn.list2dict(keys, values) #dict(zip(keys, values))


    #########################################
    ########################################


    #to check running time
    start_time = time.time()

    #read the sql database
    c = read_sql_file(basename)
    
    # import all LINKS 
    ###############
    print(f'Therion Import -- Importing all links (including splays) -- {(time.time() - start_time)}s')
    try:
        c.execute('select FROM_ID, TO_ID from SHOT')
    except OperationalError as e:
        print(f'1. Cannot find sql here: {basename}\n verify that .sql exists or that the path is correct ')
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
        print(f'2. Cannot find sql here: {basename}\n verify that .sql exists or that the path is correct ')
    
                           
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
        # if tree.endswith(('.','-'))==True:
        #     # print(f'{tree} ends with . or -')
        if tree.endswith(('.','-'))==False:#tree[0].isdigit():  
            # print(f'{tree} does not ends with . or -')
            list_tree_values.append(nodes_tree_structure[i])
            list_tree_oldi.append(nodes_id[i])
        #if tree.startwith(-) or tree.startwith(.):


    
    #create graph with all the links
    #################################
    print(f'Therion Import -- Create initial graph with all the data points (including splays) -- {(time.time() - start_time)}s')
    G = nx.Graph()
    G.add_edges_from(links_all)
    # if the nodes attributes are the same for two combined nodes, it seems that it does not affect the combining
    nx.set_node_attributes(G, coord,  pos_attr)
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
        print(f'3. Cannot find sql here: {basename}\n verify that .sql exists or that the path is correct ')
    
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
        print(f'4. Cannot find sql here: {basename}\n verify that .sql exists or that the path is correct ')
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
        duplicates.append([key for key,coord in G.nodes( pos_attr) if coord==position])
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
    dict_splays = kn.utils.list2dict(list_splays_newi, list_splays_pos)
    nx.set_node_attributes(G, dict_splays, 'splays') 

    #TREE
    #####
    print(f'Therion Import -- add fulladdress -- {(time.time() - start_time)}s')
    list_tree_newi = [index_dict.get(item, item)  for item in list_tree_oldi]
    dict_tree = kn.utils.list2dict(list_tree_newi, list_tree_values)
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
    dict_node_flag = kn.utils.list2dict(list_node_flag_newi, list_node_flag_values)
    nx.set_node_attributes(G, dict_node_flag, 'flags') 
    
    #add potential edge flags
    ############################
    # Shot Flags
    # 'dpl' = duplicate, 'srf' = surface shots
    from_edge_flag_oldi, to_edge_flag_oldi, list_edge_flag_values = extract_flags(c,'shot', return_type='lists')   
    list_from_edge_flag_newi = [index_dict.get(item, item)  for item in from_edge_flag_oldi]
    list_to_edge_flag_newi = [index_dict.get(item, item)  for item in to_edge_flag_oldi]
    dict_edge_flag = kn.utils.list2dict(list(zip(list_from_edge_flag_newi,list_to_edge_flag_newi)), list_edge_flag_values)
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

    #optionally remove flagged edges
    if remove_flagged_edges==True:
        print(f'Therion Import -- remove_flagged_edges set to True. removing flagged edges -- {(time.time() - start_time)}s')  
        kn.clean.remove_flagged_edges(G)  


    if verbose:
        print("Networkx graph successfully created from sql file !\n")

    return G
