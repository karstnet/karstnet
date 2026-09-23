import networkx as nx
import numpy as np
import yaml
import re


def add_edges(G,additional_edges, dict_address=None, additional_flag = None, comment_string=None):
    """Add edges to a networkx cave graph, based on a list of edges
    Parameters
    ----------
    G : networkx graph 
        Graph produced with the function kn.io_func.networkx_from_therion_sql

    additional_edges : list of tuple or list of lists
        list of edges to add to the graph. The edges will be created between already existing nodes, and a flag add will be added. 

        example:
        -----------
        in the case that use_fulladdress == True: additional_edges = [['full_address.0','full_address.1],['full_address.3','full_address.10]]
        if the current networkx key is used: additional_edges = [[67,2],[110,232]]

    dict_address : dict, optional
        by default this function takes the current networkx keys of the graph G. 
        When dict_address is not None, then it uses other node identifiers (for example, it can be the original node name)
        The dictionnary 
        stored in the graph attributes G.nodes('fulladdress').
        This attribute is created with the therion import function.

        example:
        --------
        Each current node id can be associated to a series of older ids, in the case that stations where regrouped because they were at the same positions.
        dict_address = {id_0:['old_id_20','old_id_1','old_id_200'], id_1:['old_id_3']}

        When the cave is processed created with the import function therion, then the original node name is stored in the 
        graph attributes G.nodes('fulladdress'), in the form of a full path from the main folder. This is a therion standard:
        dict_address = {id_0:['full_address.0','full_address.1','full_address.4'], id_1:['full_address.0']}

    """
    
    #create dictionnary to find current id based on the fulladdress
    if dict_address is not None:
        inverse_dict_address = { v: k for k, l in dict_address.items() for v in l }

    for edge in additional_edges:
        if dict_address is not None:
            #find node id based on the full address and attach the edge
            edge_from = inverse_dict_address[edge[0]]
            edge_to = inverse_dict_address[edge[1]]
            # edges[i][0] = [key for key, value in dict_address.items() if edge[0] in value ][0]
            # edges[i][1] = [key for key, value in dict_address.items() if edge[1] in value ][0]  
            # create a new edge with the flag value
        else:
            edge_from = edge[0]
            edge_to = edge[1]
            

        #check if there is flags already attached to the edges
        if G.has_edge(edge_from,edge_to):
            print(f'edge {(edge_from,edge_to)} already exists. ignore add_edge')
        else:
            G.add_edge(edge_from,edge_to,flags=['add'])

            if additional_flag is not None:
                G.edges[(edge_from,edge_to)]['flags'].append(additional_flag)

            if comment_string is not None:
                nx.set_edge_attributes(G, {(edge_from,edge_to): {'comments': comment_string}})




                
            


def flag_nodes(G,flagged_nodes, dict_address=None):
    
    """Add a string in the node attribute 'flag' of the networkx cave graph.

    Parameters
    ----------
    G : networkx graph 
        Graph produced with the function kn.io_func.networkx_from_therion_sql
    
    flagged_nodes : dictionnary of list, optional
        list of nodes to be flagged, with associated flag.
        The dictionnary key is the flag, and the values associated with the key is the list of node ids to flag.

        Example of dictionnary: 
        -----------------------
        flagged_nodes = {'ent':['full_address.0','full_address.1','full_address.4']}
        flagged_nodes = {'ent':[old_id_1,old_id_4,old_id_400]}
        or 
        flagged_nodes = {'ent':[8,45,201]}

    dict_address : dict, optional
        by default this function takes the current networkx keys of the graph G. 
        When dict_address is not None, then it uses other node identifiers (for example, it can be the original node name)
        The dictionnary 
        stored in the graph attributes G.nodes('fulladdress').
        This attribute is created with the therion import function.
        
        example:
        --------
        Each current node id can be associated to a series of older ids, in the case that stations where regrouped because they were at the same positions.
        dict_address = {id_0:['old_id_20','old_id_1','old_id_200'], id_1:['old_id_3']}

        When the cave is processed created with the import function therion, then the original node name is stored in the 
        graph attributes G.nodes('fulladdress'), in the form of a full path from the main folder. This is a therion standard:
        dict_address = {id_0:['full_address.0','full_address.1','full_address.4'], id_1:['full_address.0']}


    """

    print(f'flag_nodes - adding manual node flags: {flagged_nodes.keys()}')

    #create dictionnary to find current id based on the fulladdress
    if dict_address is not None:
        inverse_dict_address = { v: k for k, l in dict_address.items() for v in l }

        
    for flag in flagged_nodes.keys():
        #loop throught the flags
        for node in flagged_nodes[flag]:
            
            #if the node is just a string or int, or... and not a list, then it means there is no comments attached
            if type(node)!=list:
                if dict_address is not None:
                    id_node = inverse_dict_address[node]   
                else:
                    id_node = node    

                #FIRST CHECK IF NODE EXISTS! IT MAY HAVE BEEN REMOVED DURING CLEANING   
                if G.has_node(id_node):              
                    #check if node already has a flag
                    #if not, create a new list of flag(s) attached to the node
                    if G.nodes('flags')[id_node] is None:
                        #create flag on node
                        nx.set_node_attributes(G, {id_node:[flag]}, name='flags')

                    #if yes, append the new flag to the list
                    elif G.nodes('flags')[id_node] is not None:
                        #to avoid repetition:
                        if flag not in G.nodes[id_node]['flags']:
                            G.nodes[id_node]['flags'].append(flag)
                else:
                    print(f'Node {id_node} does not exist anymore.')


            # Check if its a list, then it means there is a comment attached
            elif len(node)>1:
                # print(node, 'bigger than 1')
                if dict_address is not None:
                    id_node = inverse_dict_address[node[0]]   
                else:
                    id_node = node[0] 
                comment = node[1] 

                #FIRST CHECK IF NODE EXISTS! IT MAY HAVE BEEN REMOVED DURING CLEANING   
                if G.has_node(id_node):
                    #check if node already has a flag
                    #if not, create a new list of flag(s) attached to the node
                    if G.nodes('flags')[id_node] is None:
                        # pass
                        #create flag on node
                        nx.set_node_attributes(G, {id_node:[flag]}, name='flags')

                    #if yes, append the new flag to the list
                    elif G.nodes('flags')[id_node] is not None:
                        #print(flag,id_node, fulladdress)
                        if flag not in G.nodes[id_node]['flags']:
                            G.nodes[id_node]['flags'].append(flag)


                    #check if there is no comments already attached to the node
                    if G.nodes('comments')[id_node] is None:
                        nx.set_node_attributes(G, {id_node:[comment]}, name='comments')
                    # if there is already some comments attached to the node, add it to the list
                    #!!! this could be an issue if comments is not a list
                    elif G.nodes('comments')[id_node] is not None:
                        #to avoid repetition
                        if comment not in G.nodes[id_node]['comments']:
                            G.nodes[id_node]['comments'].append(comment)

                else:
                    print(f'Node {id_node} does not exist anymore.')

                    
                



def flag_edges(G, flagged_edges, dict_address = None):
    """Add a string in the edge attribute 'flag' of the networkx cave graph.

    Parameters
    ----------
    G : networkx graph 
        Graph produced with the function kn.io_func.networkx_from_therion_sql

    flagged_edges : dictionnary of list of tuple or list of lists, optional
        lists of edges to be flagged with corresponding flag to add. Any flag can be added. by default None
        However, for the duplicate edges and surface edges to be removed, it is necessary to use the correct flags.
        It is also possible to use the add edges with this dictionnary instead of using the 'add_edges' option.
        List of flagged edge that will be removed by default:  'dpl', 'srf', 'art', 'rmv', 'spl'

    use_fulladdress : bolean, default False
        by default this function takes the current networkx keys of the graph G. 
        When True, then it uses the "fulladdress" string stored in the graph attributes G.nodes('fulladdress').
        This attribute is created with the therion import function.

    Example of dictionnary: 
    ----------------------
    if use_fulladdress == True:
    flagged_edges = {'dpl':[['full_address.0','full_address.1],['full_address.3','full_address.10]], 'srf':[[full_address.3,full_address.2]]}   

    if we use the networkx key:
    flagged_edges = {'dpl':[[1,3],[11,5]], 'srf':[[5,11]]}   
    """
    print(f'flag_edges - adding manual edges flags: {flagged_edges.keys()}')

    #create dictionnary to find current id based on the fulladdress
    if dict_address is not None:
        inverse_dict_address = { v: k for k, l in dict_address.items() for v in l }

    
    #loop through all the flags
    for flag in flagged_edges.keys():
        print(flag)
        #loop through all the edges for each flag
        for edge in flagged_edges[flag]:
            if dict_address is not None:
                edge_from = inverse_dict_address[edge[0]]
                edge_to = inverse_dict_address[edge[1]]
            else:
                edge_from = edge[0]
                edge_to = edge[1]
            #check if edge exists already (for example to add a duplicate flag on an exisiting edge)
            if G.has_edge(edge_from,edge_to):
            #check if there is flags already attached to the edges
            #if not, create a new edge with the dictionnary 'flags' and the flag value
                if 'flags' not in G[edge_from][edge_to]:
                    # print(f'adding {flag} to edge {edge_from}-{edge_to}')
                    nx.set_edge_attributes(G,{(edge_from,edge_to):{'flags':['add',flag]}})
                #if yes, append the flag to the list
                elif 'flags' in G[edge_from][edge_to]:
                    # Supprime 'dpl' si on ajoute 'add'
                    if flag == 'add' and 'dpl' in G[edge_from][edge_to]['flags']:
                        G[edge_from][edge_to]['flags'].remove('dpl')
                    # print(f'appending {flag} to edge {edge_from}-{edge_to}')
                    #to avoid repetition
                    if 'add' not in G[edge_from][edge_to]['flags']:
                        # print('true')
                        G[edge_from][edge_to]['flags'].append('add')
                    if flag not in G[edge_from][edge_to]['flags']:
                        G[edge_from][edge_to]['flags'].append(flag)
                        

                #if there is a comment attached to the edge, create
                if len(edge)==3:
                    # if  the edge already exists, but there is no comments yet, create the comment attribute for this edge
                    if 'comments' not in G[edge_from][edge_to]:
                        nx.set_edge_attributes(G,{(edge_from,edge_to):{'comments':[edge[2]]}})
                        
                    #in the case that there is already comments on the edge, apped the current comment
                    elif 'comments' in G[edge_from][edge_to]:
                        if edge[2] not in G[edge_from][edge_to]['comments']:
                            G[edge_from][edge_to]['comments'].append(edge[2])
                        



            #if edge does not exist yet, then just create a new edge with the appropriate flag name            
            else:
                # add the flag add by default?? if yes, implement this feature (20 fev 2025)
                # print(f'creating edge and adding {flag} to edge {edge_from}-{edge_to}')
                G.add_edge(edge_from,edge_to,flags=['add', flag])
                # print(f'graph length: {len(G)}')
                if len(edge)==3: 
                    #add comments if there is any               
                    nx.set_edge_attributes(G,{(edge_from,edge_to):{'comments':[edge[2]]}})





# def remove_edges()

def remove_flagged_edges(G, flags_to_remove=['srf','dpl','rmv','art','spl'], attribute_name = 'flags'):
    """Remove edges flagged with certain strings.

    Parameters
    ----------
    G : networkx graph 
        Graph produced with the function kn.io_func.networkx_from_therion_sql

    flags_to_remove : list of strings, optional
        list of the flags for which edges should be removed, by default ['srf','dpl','rmv','art','spl']
        - 'dpl' : duplicate
        - 'srf' : surface
        - 'art' : artificial
        - 'rmv' : remove
        - 'spl' : splay (for example when a shot is made in a large room, star shots, ...)
    attribute_name : string
        name of the attribute attached to the graph 

    if the flag 'de' is present ont


    """
    #edges_to_remove = list(dict(nx.get_edge_attributes(G,'flags')).keys()) 
    #extract flag unique values into a list
    flags = {x for l in list(nx.get_edge_attributes(G,attribute_name).values()) for x in l}
    #loop through the unique flags and 


    for flag in flags:
        #only remove edges with flag surface, duplicate, or remove
        # the 'de' flags superseed the removing flags and prevents the points from being removed. its applied to edges that were fasly flagged
        if flag in flags_to_remove:
            list_edges = [edge for edge, action in nx.get_edge_attributes(G,attribute_name).items() if ((flag in action) & ('de' not in action))]
            #print('remove ', flag, list_edges)
            #remove edges    
            G.remove_edges_from(list_edges)

        else:
            pass
            #print('not removed', flag)
    #remove nodes that were isolated when removing the edges
    G.remove_nodes_from(list(nx.isolates(G)))



    # print(f'Initial Graph size: {len(G)}, Graph size after removing flagged edges: {len(H)}')
    # return H




# ----------------------------------------------------------------------------
def get_potential_connection(
        G,
        dist_horiz_max,
        dist_vert_max,
        node_deg=1,
        exclude_neighbors_up_to_edge=3,
        pos_attr='pos',
        return_dist=False,
        return_angle=False):
    """
    Retrieves list of potential new edges for a graph.

    For every node `u` of given degree `node_deg`:

    1. the nodes within a cylindrical box of horizontal radius `dist_horiz_max`
    and height `dist_vert_max` centered at `u` are retrieved in a list

    2. The nodes of this list at a distance in number of edges less than or
    equal to `exclude_neighbors_up_to_edge` are excluded.

    3. Then, the nearst node `v` to `u` (checking horizontal Euclidean distance)
    in this list is identified, and the tuple `(u, v)` is considered as a potential
    edge to be added to the graph (new connection)

    Parameters
    ----------
    G : networkx.Graph
        graph

    dist_horiz_max : float (positive)
        maximal horizontal distance (radius of the cylinder around checked nodes)

    dist_vert_max : float (positive)
        maximal vertical distance (height of the cylinder around checked nodes)

    exclude_neighbors_up_to_edge : int, default: 3
        distance in number of edges, nodes to a distance from a checked node at
        distance smaller than or equal to `exclude_neighbors_up_to_edge` are excluded
        from potential edge with the checked node

    node_deg : int, default: 1
        degree of the nodes to be checked as an extremity for
        potential new edge

    return_dist : bool, default: False
        if `True`, the length of the potential new edges (distance between the
        two extremities) are returned

    return_angle : bool, default: False
        if `True`, xxxxthe length of the potential new edges (distance between the
        two extremities) are returned

    Returns
    -------
    edges_list : list of 2-tuples
        list of potential new edges, each element is a 2-tuple (u, v), where
        u and v are the node ids of the two extremities of a potential new edge

    dist_list : list of floats, optional
        returned if `return_dist=True`, list of same length as `edges_list`, of the
        lengths of the potential new edges (distance between the two extremities)

    angle_list : list of lists of float(s), optional
        returned if `return_angle=True`, list of same length as `edges_list`, where
        `angle_list[i]` is the list of angles between the potential new edge `edge_list[i]`
        and the existing edge(s) whose one extremity is the node ``edge_list[i][0]`;
        each angle is in degree in the interval [0, 180]
    """
    # Set dictionary to convert node label (id) to node index, and vice versa
    node_label2index = {u:i for i, u in enumerate(G.nodes())}
    node_index2label = {i:u for i, u in enumerate(G.nodes())}

    pos = nx.get_node_attributes(G, pos_attr)
    pos = {k:np.asarray(v) for k, v in pos.items()} # convert tuple to array

    # Matrix of all positions
    pos_arr = np.asarray(list(pos.values()))

    rh2 = dist_horiz_max**2
    rv = dist_vert_max

    edges_list = []
    for u in G.nodes():
        if G.degree(u) != node_deg:
            continue

        # Get array (sel_ind_arr) of index of nodes within the cylindrical box centered at u:
        ind = node_label2index[u]
        lag = pos_arr - pos_arr[ind]
        disth2 = np.sum(lag[:,:2]**2, axis=1)
        distv = np.abs(lag[:,2])
        sel = np.all((disth2 <= rh2, distv <= rv), axis=0) # True at least for node u (at index ind)
        sel_ind_arr = np.where(sel)[0]
        if sel_ind_arr.size <= 1:
            continue

        # Get array (neigh_to_exclude_ind_arr) of neighbors index to exclude
        neigh_to_exclude_ind_arr = np.asarray([node_label2index[v] for v in list(nx.single_source_shortest_path_length(G, u, cutoff=exclude_neighbors_up_to_edge).keys())])

        # Update sel_ind_array
        sel_ind_arr = np.setdiff1d(sel_ind_arr, neigh_to_exclude_ind_arr)

        if sel_ind_arr.size == 0:
            continue

        # Get v the nearest node to u : potential edge (u, v)
        min_ind = sel_ind_arr[np.argmin(disth2[sel_ind_arr])]
        v = node_index2label[min_ind]
        edges_list.append((u, v))

    out = [edges_list]

    if return_dist:
        if len(edges_list):
            dist_list = [pos_arr[node_label2index[u]] - pos_arr[node_label2index[v]] for u, v in edges_list]
            dist_list = list(np.sqrt(np.sum(np.asarray(dist_list)**2, axis=1)))
        else:
            dist_list = []
        out.append(dist_list)

    if return_angle:
        angle_list = []
        for u, v in edges_list:
            pos_u = pos_arr[node_label2index[u]]
            pos_v = pos_arr[node_label2index[v]]
            uv = pos_v - pos_u
            uv_norm = np.sqrt(np.sum(uv**2))
            a = []
            for _, vi in G.edges(u):
                pos_vi = pos_arr[node_label2index[vi]]
                uvi = pos_vi - pos_u
                uvi_norm = np.sqrt(np.sum(uvi**2))
                a.append(np.rad2deg(np.arccos(np.sum(uv*uvi)/(uv_norm*uvi_norm))))
            angle_list.append(a)
        out.append(angle_list)

    if len(out) == 1:
        out = out[0]
    else:
        out = tuple(out)

    return out
# ----------------------------------------------------------------------------

def find_disconnected_node(G, H):
    """Identify node of degree one that were initially of degree 2 or more at an earlier stage of the cleaning process

    Parameters
    ----------
    G : networkx graph
        Graph exported from therion, containing all the data
    H : networkx graph
        Graph without the surface and duplicate shots

    Returns
    -------
    list
        list of disconnected node ids
    """    
    # G_raw = load_raw_therion_data(basename)    
    # G = load_therion_without_flagged_edges(basename)
    print( 'There is ', nx.number_connected_components(G), 'connected components in the original graph')
    print( 'There is ', nx.number_connected_components(H), 'connected components in the graph without flagged edges')
    
    closeby_all = []
    keys_disconnected_all =[]
    
    #cc_number = 0
    if nx.is_connected(H) == False:
        #iterate through the connectec components to find nodes where disconnection occured
        #search for the nodes that used to be degree >1 and are now degree 1.
        for i, subgraph_index in enumerate(nx.connected_components(H)):
            #print(i,subgraph_index )
            subgraph = nx.subgraph(H, subgraph_index)

            keys_disconnected_subgraph=[]   
            #find all the nodes where disconnection happened
            #look for all the nodes degree smaller in the cleaned file than in the original file 
            #keys_disconnected_subgraph = [k for k, v in dict(subgraph.degree()).items() if v == 1 and G_raw.degree()[k] >1]   
            for k in subgraph.nodes(): #dict(subgraph.degree()).items():
                # print(k)
                if subgraph.degree()[k]==1 and G.degree()[k] >1:
                    keys_disconnected_subgraph.append(k)

                
            keys_disconnected_all = keys_disconnected_all + keys_disconnected_subgraph
        return keys_disconnected_all
    else:
        print('There is no disconnected components, no need to merge')

# ----------------------------------------------------------------------------
def find_splays_via_letters(G, filepath=None):
    """Identifies splays connected to each station in a karst network graph and optionally exports the result to a YAML file.

    This function processes a networkx graph, where each node includes a `fulladdress` attribute.
    It identifies for each station the connected splays, whose last address element ends with a letter (e.g., 'A', 'b', etc.).

    The result is a dictionary where each key is the full address of a station and each value is a list of the full addresses
    of its connected splays. If a file path is provided, this dictionary is saved as a YAML file.


    Args:
        G (networkx graph): Graph exported from therion, containing all the data
        filepath (str, optional): Path to the YAML file. Defaults to None.

    Returns:
        dict: Dictionary with each station's full address and the list of its connected splays' full addresses.
         - Each key is the fulladdress of a node (station),
         - Each value is a list of fulladdresses of connected splays (nodes whose fulladdress ends with a letter).
    """    
    splay_dict = {} #Dictionnaire qui va contenir le résultat

    for node_id, data in G.nodes(data=True): 
        station_addr = data.get('fulladdress', [])
        if not station_addr:
            continue
        station_full = station_addr[-1]  # Dernière adresse de la station

        splay_list = [] # Liste temporaire pour stocker les splays connectés à cette station

        for neighbor in G.neighbors(node_id): # On parcourt tous les voisins (nœuds connectés) de cette station
            neighbor_data = G.nodes[neighbor] # On récupère les attributs du voisin
            neighbor_addr = neighbor_data.get('fulladdress', []) # Sa fulladdress
            if neighbor_addr:
                neighbor_full = neighbor_addr[-1]
                # Vérifie si le voisin est un splay (adresse terminée par une lettre)
                if re.match(r'.*[A-Za-z]$', neighbor_full):
                    splay_list.append(neighbor_full)

        if splay_list:
            splay_dict[station_full] = splay_list

    # Ajout automatique de l'extension si manquante
    if not filepath.endswith('.yaml'):
        filepath += '.yaml'

    if filepath is not None:
        # Sauvegarde du dictionnaire dans le fichier YAML
        with open(filepath, 'w') as file:
            yaml.dump(splay_dict, file, default_flow_style=False, sort_keys=False)

        print(f"Splay structure saved to {filepath}")
    
    return splay_dict
# ----------------------------------------------------------------------------
def find_splays_via_structure(G, filepath=None, degree_threshold=5, neighbors_degree_threshold=1, neighbors_count_threshold=4):
    """
    Identifies splays based on graph structure: high-degree stations connected to multiple terminal nodes.
    Optionally exports the result to a YAML file.

    This function processes a networkx graph where each node includes a `fulladdress` attribute.
    It identifies stations (nodes with high degree) that are connected to multiple neighbors of degree 1,
    and returns a dictionary mapping each station’s fulladdress to the list of connected splays (neighbor fulladdresses).

    If a file path is provided, the dictionary is saved as a YAML file.

    Args:
        G (networkx graph): Graph exported from therion, containing all the data
        filepath (str, optional): Path to the YAML file to export. Defaults to None.
        degree_threshold (int): Minimum degree for a node to be considered a potential central station.
        neighbors_degree_threshold (int): Maximum degree for a neighbor to be considered a splay.
        neighbors_count_threshold (int): Minimum number of splay-like neighbors required to retain the central station.

    Returns:
        dict: Dictionary with each central station’s fulladdress as key and the list of connected splays’ fulladdresses as value.
    """
    splay_dict = {}

    for node_id in G.nodes:
        node_data = G.nodes[node_id]
        station_addr = node_data.get('fulladdress', [])
        if not station_addr or G.degree(node_id) <= degree_threshold:
            continue

        station_full = station_addr[-1]
        splay_list = []

        for neighbor in G.neighbors(node_id):
            if G.degree(neighbor) <= neighbors_degree_threshold:
                neighbor_addr = G.nodes[neighbor].get('fulladdress', [])
                if neighbor_addr:
                    splay_list.append(neighbor_addr[-1])

        if len(splay_list) >= neighbors_count_threshold:
            splay_dict[station_full] = splay_list

    # Ajout automatique de l’extension si nécessaire
    if filepath:
        if not filepath.endswith('.yaml'):
            filepath += '.yaml'
        with open(filepath, 'w') as file:
            yaml.dump(splay_dict, file, default_flow_style=False, sort_keys=False)
        print(f"Splay structure saved to {filepath}")

    return splay_dict

# ----------------------------------------------------------------------------
def add_splay(G, splay_dict, dict_address=None):
    """
    Adds information about splay nodes to a NetworkX graph.

    For each station node listed in the input dictionary (splay_dict) :
    - Locates the corresponding node in the graph G (using fulladdress or node ID),
    - Locates each associated splay node,
    - Retrieves the 3D position of each splay (from the 'pos' attribute),
    - Adds a 'splays' attribute to the station node, containing the list of 3D positions,
    - Flags each splay node with 'spl' in its 'flags' attribute.

    This function assumes that splay_dict is of the form: {station_fulladdress: [splay_fulladdress1, splay_fulladdress2, ...]}

    Args:
        G (networkx.Graph): The input graph to be updated with splay data.
        splay_dict (dict): A dictionary that lists, for each station, the names (fulladdresses) of its splay nodes.
        dict_address (dict, optional): A dictionary that helps convert a fulladdress into a node ID in the graph.
            If provided, it's used to convert fulladdresses into node IDs.

    Returns:
        None. Updates the graph G.
    """
    #create dictionnary to find current id based on the fulladdress
    if dict_address is not None:
        inverse_dict_address = { v: k for k, l in dict_address.items() for v in l }

    # Pour chaque station dans le dictionnaire (YAML)
    for station_fulladdr, splay_addresses in splay_dict.items(): # Parcourt chaque station et sa liste de splays
        ## Trouver l'ID du noeud correspondant à cette station
        if dict_address is not None:
            station_node=inverse_dict_address.get(station_fulladdr)
        else: 
            station_node=station_fulladdr    

        if station_node is not None and station_node in G: # Si on a bien trouvé la station dans le graphe
            splay_positions = [] #Liste vide pour stocker les positions des splays

            for splay_addr in splay_addresses:  # Pour chaque splay associé à cette station
                # Trouver l'ID du splay correspondant à cette fulladdress
                if dict_address is not None:
                    splay_node=inverse_dict_address.get(splay_addr)
                else:
                    splay_node=splay_addr

                if splay_node is not None and splay_node in G:
                    pos = G.nodes[splay_node].get('pos')
                    if pos:
                        splay_positions.append(pos)
                    # Ajout du flag 'spl' au nœud splay
                    flags = G.nodes[splay_node].get('flags', [])
                    if 'spl' not in flags:
                        flags.append('spl')
                        G.nodes[splay_node]['flags'] = flags

            # Ajouter la liste des positions des splays à la station
            G.nodes[station_node]['splays'] = splay_positions
        else:
            print(f'{station_fulladdr}, does not exist in the graph - double check fulladdress or node ID')

    print("Splays successfully assigned to the corresponding stations.")
    print(f"{len(splay_dict)} stations updated with 'splays'.")


# ----------------------------------------------------------------------------
def transform_station_splay_dict(original_dict):
    """
    Transforms a dictionary of station-to-splays into a format suitable for edge flagging.

    This function takes a dictionary where each key is a station's fulladdress and each value is a list of splay fulladdresses connected to that station.
    It reformats this information into a dictionary with a single key 'spl', whose value is a list of [station, splay] pairs. This format is required by
    the flag_edges() function in the clean.py.

    Args:
        original_dict (dict): A dictionary that lists, for each station fulladdress, the splays associated with it.
        
    Returns:
        dict: A new dictionary with a single key 'spl' and a list of [station, splay] pairs as its value.
    """
    new_dict = {'spl': []}
    for station, splays in original_dict.items():
        for splay in splays:
            new_dict['spl'].append([station, splay])
    return new_dict

# ----------------------------------------------------------------------------
def find_shortest_connection_between_components(G, distance_threshold=15):
    """
    Finds the shortest Euclidean-distance connections between pairs of connected components
    in a graph, as long as they are within a given threshold.

    Parameters
    ----------
    G : networkx.Graph
        Graph exported from therion, containing all the data
    distance_threshold : float
        Maximum distance under which two nodes from different components are considered connectable.

    Returns
    -------
    List of tuples:
        Each tuple represents a connection and contains:
        (component_index_1, component_index_2, node_id_1, node_id_2, distance)
    """
    # Get the list of connected components in the graph
    components = list(nx.connected_components(G))

    # Extract node positions
    pos = nx.get_node_attributes(G, 'pos')

    shortest_links = []

    # Compare each unique pair of connected components
    for i in range(len(components)):
        for j in range(i + 1, len(components)):
            comp_i = components[i]
            comp_j = components[j]

            min_dist = float('inf')  # Initialize minimum distance
            best_pair = None         # Store the best pair of nodes

            # Search for the pair of nodes (one in each component) with the shortest distance
            for u in comp_i:
                for v in comp_j:
                    if u in pos and v in pos:
                        dist = np.linalg.norm(np.array(pos[u]) - np.array(pos[v]))
                        if dist < distance_threshold and dist < min_dist:
                            min_dist = dist
                            best_pair = (i, j, u, v, round(dist, 2))

            # If a valid pair was found under the threshold, store it
            if best_pair:
                shortest_links.append(best_pair)

    return shortest_links

#######################################
########################################

# def clean_graph()

