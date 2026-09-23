"""
The module `utils` contains basics functions.
"""


import numpy as np
import networkx as nx
import pandas as pd
import scipy
import os

# -------------------------------------------------------------------------
def rescale_array(arr, out_min, out_max, exponent=1, nan_val=None, in_min=None, in_max=None):
    """Rescales array values linearly.
    
    Parameters
    ----------
    arr : numpy array of floats
        input array

    out_min : float
        min value in the output array

    out_max : float
        max value in the output array

    exponent : float or int, default: 1
        positive value, polynomial exponent of the rescaling;
        by default (1): linear rescaling 
    
    nan_val : float (or int), optional 
        if specified (not `None`), the nan values in the input array
        are replaced by `nan_val` in the output array
    
    in_min : float, optional
        value used as the "input minimum";
        by default (`None`), minimal value of non nan entries in `arr` is used

    in_max : float, optional
        value used as the "input maximum";
        by default (`None`), maximal value of non nan entries in `arr` is used
    
    Returns
    -------
    out_arr : numpy array of floats
        array of same shape as `arr`, values of input array linearly rescaled
        between `out_min` and `out_max` (with nan replaced by `nan_val` if specified)
    """
    if in_min is None:
        in_min = np.nanmin(arr)
    
    if in_max is None:
        in_max = np.nanmax(arr)
    
    arr_out = np.full(arr.shape, np.nan)
    arr_out[~np.isnan(arr)] = out_min + np.pow((arr[~np.isnan(arr)] - in_min) / (in_max - in_min), exponent) * (out_max - out_min)

    if nan_val is not None:
        np.putmask(arr_out, np.isnan(arr_out), nan_val)
    
    return arr_out
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
def get_pos2d(G, pos_attr='pos'):
    """
    Gets node position in 2D (ignoring z-coordinate if it exists) of a networkx graph.

    Parameters
    ----------
    G : networkx.Graph
        graph with 2D or 3D position as node attribute

    node_attr : str, default: 'pos'
        name of the node attribute for position

    Returns
    -------
    pos2d : dict
        keys are node names, values are the 2D position of the nodes (array)
    """
    pos2d = {key: np.asarray(value[0:2]) for key, value in nx.get_node_attributes(G, pos_attr).items()}
    return pos2d
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
def get_pos3d(G, pos_attr='pos'):
    """
    Gets node position in 3D of a networkx graph.

    Parameters
    ----------
    G : networkx.Graph
        graph with 2D or 3D position as node attribute

    node_attr : str, default: 'pos'
        name of the node attribute for position

    Returns
    -------
    pos3d : dict
        keys are node names, values are the 3D position of the nodes (array)
    """
    pos = nx.get_node_attributes(G, pos_attr)
    try:
        if len(list(pos.values())[0]) == 2: # pos is 2D
            pos3d = {k:np.append(np.asarray(v[:2]), 0.0) for k, v in pos.items()}
        else: # pos is 3D
            pos3d = {k:np.asarray(v) for k, v in pos.items()}
            # pos3d = pos
    except:
        # e.g. when pos={}
        pos3d = {}

    return pos3d
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
def get_posz(G, pos_attr='pos'):
    """
    Gets z coordinates of node position of a networkx graph.

    Parameters
    ----------
    G : networkx.Graph
        graph with 2D or 3D position as node attribute

    node_attr : str, default: 'pos'
        name of the node attribute for position

    Returns
    -------
    posz : dict
        keys are node names, values are the z position of the nodes
    """
    pos = nx.get_node_attributes(G, pos_attr)
    try:
        if len(list(pos.values())[0]) == 2: # pos is 2D
            posz = {k:0.0 for k in pos.keys()}
        else: # pos is 3D
            posz = {k:v[2] for k, v in pos.items()}
            # pos3d = pos
    except:
        # e.g. when pos={}
        posz = {}

    return posz
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
def get_node_attribute_names(G):
    """
    Gets the set of node attribute names.

    Parameters
    ----------
    G : networkx.Graph
        graph

    Returns
    -------
    node_attr_list : list
        list of node attribute names: all attribute names (keys) attached 
        to at least one graph node
    """
    return list(set([k for n in G.nodes for k in G.nodes[n].keys()]))
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
def get_edge_attribute_names(G):
    """
    Gets the set of edge attribute names.

    Parameters
    ----------
    G : networkx.Graph
        graph

    Returns
    -------
    edge_attr_list: list
        list of edge attribute names: all attribute names (keys) attached 
        to at least one edge
    """
    # get list of node keys
    return list(set([k for n in G.edges for k in G.edges[n].keys()]))
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
def remove_node_attribute(G, attr):
    """
    Removes (if present) a node attribute (from all nodes).

    This function operates inplace on the graph `G`.

    Parameters
    ----------
    G : networkx.Graph
        graph

    attr : str
        name of the node attribute to be removed
    """
    for u in G.nodes():
        if attr in G.nodes[u].keys():
            del G.nodes[u][attr]
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
def remove_edge_attribute(G, attr):
    """
    Removes (if present) an edge attribute (from all edges).

    This function operates inplace on the graph `G`.

    Parameters
    ----------
    G : networkx.Graph
        graph

    attr : str
        name of the edge attribute to be removed
    """
    for e in G.edges():
        if attr in G.edges[e].keys():
            del G.edges[e][attr]
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
def copy_node_attribute(G, attr_src, attr_target):
    """
    Copies a node attribute.

    This function operates inplace on the graph `G`.

    Parameters
    ----------
    G : networkx.Graph
        graph

    attr_src : str
        name of the source node attribute

    attr_target : str
        name of the target node attribute
    """
    nx.set_node_attributes(G, nx.get_node_attributes(G, attr_src), attr_target)
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
def copy_edge_attribute(G, attr_src, attr_target):
    """
    Copies an edge attribute.

    This function operates inplace on the graph `G`.

    Parameters
    ----------
    G : networkx.Graph
        graph

    attr_src : str
        name of the source edge attribute

    attr_target : str
        name of the target edge attribute
    """
    nx.set_edge_attributes(G, nx.get_edge_attributes(G, attr_src), attr_target)
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
def get_node_label2index_dict(G):
    """
    Gets the correspondance from node label to node index.

    Parameters
    ----------
    G : networkx.Graph
        graph

    Returns
    -------
    node_label2index : dict
        dictionary where the keys are the node label and the value are
        the node index in the list `list(G.nodes())`
    """
    # Set dictionary to convert node label to node index
    node_label2index = {u:i for i, u in enumerate(G.nodes())}
    return node_label2index
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
def get_node_index2label_dict(G):
    """
    Gets the correspondance from node index to node label.

    Parameters
    ----------
    G : networkx.Graph
        graph

    Returns
    -------
    node_label2index : dict
        dictionary where the keys are the node index in the list 
        `list(G.nodes())` and the value are the node label
        
    """
    # Set dictionary to convert node index to node label
    node_index2label = {i:u for i, u in enumerate(G.nodes())}
    return node_index2label
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
def get_node_data_frame(
        G, 
        pos_attr='pos', 
        include_pos=True, 
        attr_list=None, 
        attr_ignore_list=None,
        round=None,
        verbose=True):
    """
    Gets the data frame of node attributes.

    Note: if an attribute (other than position) is an array (non-scalar), it is splitted 
    into several columns with names '<attr_name>_0', '<attr_name>_1', ..., 
    where '<attr_name>' is the name of the attribute; 
    position attribute is splitted into 'x', 'y'[, 'z'].

    Parameters
    ----------
    G : networkx.Graph
        graph

    pos_attr : str, default: 'pos'
        name of the node attribute for position

    include_pos : bool, default: True
        indicates if node position is included in the data frame

    attr_list : list of str, optional
        list of node attribute names (in addition to position) to be 
        included in the data frame;
        if `None` (default), all node attributes are included (except 
        those in `attr_ignore_list`)

    attr_ignore_list : list of str, optional
        list of node attribute names to be ignored in the data frame;
        if `None` (default), no node attribute is ignored

    round : int or dict-like, optional
        if `int`, number of decimals for rounding the node attribute values;
        if `dict-like`, keys are the node attribute names, and values are
        the number of decimals for rounding the corresponding node attribute values;
        if `None` (default), no rounding is performed

    verbose : bool, default: True
        indicates if some attributes can not be imported in the data frame 
        (e.g. inhomogeneous shape) 

    Returns
    -------
    df : pandas.DataFrame
        data frame of node attributes, with columns:
        - 'id' : node label
        - 'x', 'y'[, 'z'] : node position coordinates
        - other columns corresponding to other node attributes
    """
    # Get node list (ids)
    nodes_list = list(G.nodes())
    # nodes_list = [u for u in G.nodes()] # equivalent

    # Set data frame with node ids
    df = pd.DataFrame(nodes_list, columns=['id'])

    # Add node position
    if include_pos:
        pos_arr = np.asarray(list(nx.get_node_attributes(G, pos_attr).values()))
        if pos_arr.shape[1] == 2:
            pos_names = ['x', 'y']
        elif pos_arr.shape[1] == 3:
            pos_names = ['x', 'y', 'z']
        else:
            raise ValueError('Node positions must be 2D or 3D.')

        df[pos_names] = pos_arr
    
    # Add other node attributes
    attr_names = get_node_attribute_names(G)
    for attr in attr_names:
        if attr == pos_attr:
            continue

        if attr_list is not None and attr not in attr_list:
            continue

        if attr_ignore_list is not None and attr in attr_ignore_list:
            continue

        try:
            arr = np.asarray(list(nx.get_node_attributes(G, attr, default=np.nan).values()))
        except Exception as exc:
            if verbose:
                print(f'Warning: node attribute "{attr}" can not be imported in the data frame: {exc}')
            continue
    
        if arr.ndim == 1:
            # scalar attribute
            df[attr] = arr
        else: 
            # array attribute
            df[[attr + f'_{i}' for i in range(arr.shape[1])]] = arr

    # Rounding
    if round is not None:
        df = df.round(round)

    return df
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
def get_edge_data_frame(
        G, 
        edge_name=['from', 'to'], 
        attr_list=None, 
        attr_ignore_list=None,
        round=None,
        verbose=True):
    """
    Gets the data frame of edge attributes.

    Note: if an attribute is an array (non-scalar), it is splitted 
    into several columns with names '<attr_name>_0', '<attr_name>_1', ..., 
    where '<attr_name>' is the name of the attribute.

    Parameters
    ----------
    G : networkx.Graph
        graph

    edge_name : list of 2 str, default: ['from', 'to']
        names of the columns corresponding to the extremities of 
        the edges
    
    attr_list : list of str, optional
        list of edge attribute names (in addition to `edge_name`) to be 
        included in the data frame;
        if `None` (default), all edge attributes are included (except 
        those in `attr_ignore_list`)

    attr_ignore_list : list of str, optional
        list of edge attribute names to be ignored in the data frame;
        if `None` (default), no edge attribute is ignored

    round : int or dict-like, optional
        if `int`, number of decimals for rounding the edge attribute values;
        if `dict-like`, keys are the edge attribute names, and values are
        the number of decimals for rounding the corresponding edge attribute values;
        if `None` (default), no rounding is performed

    verbose : bool, default: True
        indicates if some attributes can not be imported in the data frame 
        (e.g. inhomogeneous shape) 
        
    Returns
    -------
    df : pandas.DataFrame
        data frame of edge attributes, with columns:
        - edge_name : node labels of the edge extremities
        - other columns corresponding to the edge attributes
    """
    # Get edge array
    edges_arr = np.asarray(list(G.edges()))
    # edges_arr = np.asarray([(u, v) for u, v in G.edges()]) # equivalent

    # Set data frame with edge extremities
    df = pd.DataFrame(edges_arr, columns=edge_name)

    # Add edge attributes
    attr_names = get_edge_attribute_names(G)
    for attr in attr_names:
        if attr_list is not None and attr not in attr_list:
            continue

        if attr_ignore_list is not None and attr in attr_ignore_list:
            continue

        try:
            arr = np.asarray(list(nx.get_edge_attributes(G, attr, default=np.nan).values()))
        except Exception as exc:
            if verbose:
                print(f'Warning: edge attribute "{attr}" can not be imported in the data frame: {exc}')
            continue

        if arr.ndim == 1:
            # scalar attribute
            df[attr] = arr
        else: 
            # array attribute
            df[[attr + f'_{i}' for i in range(arr.shape[1])]] = arr

    # Rounding
    if round is not None:
        df = df.round(round)

    return df
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
def get_edge_length_dict(G, pos_attr='pos', exponent=1):
    """
    Computes the length of every edge, raised to the given exponent.

    Parameters
    ----------
    G : networkx.Graph
        graph

    pos_attr : str, default: 'pos'
        name of the node attribute for position

    exponent : int or float, default: 1
        the edge length are raised to the power `exponent`

    Returns
    -------
    edge_length_dict : dict
        dictionary where the keys are the edges, and the values are
        the edge lengths raised to the power `exponent`
    """
    edge_length_dict = {e: float(np.power(np.sum((np.asarray(G.nodes[e[0]][pos_attr]) - np.asarray(G.nodes[e[1]][pos_attr]))**2), 0.5*exponent)) for e in G.edges()}

    return edge_length_dict
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
def mean_edge_length(G, pos_attr='pos'):
    """
    Computes the mean edge length.

    Parameters
    ----------
    G : networkx.Graph
        graph

    pos_attr : str, default: 'pos'
        name of the node attribute for position

    Returns
    -------
    mean_edge_length : float
        mean edge length        
    """
    if G.number_of_edges() == 0:
        return 0.0
    
    mean_edge_length = \
        float(np.mean(np.sqrt(np.sum(np.asarray(
            [np.asarray(G.nodes[u][pos_attr]) - np.asarray(G.nodes[v][pos_attr]) for u, v in G.edges()]
        )**2, axis=1))))

    return mean_edge_length
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
def euclidean_diameter(G, pos_attr='pos'):
    """
    Computes the Euclidean diameter.

    The Euclidean diameter is the maximum of the Euclidean distance
    between all pairs of nodes in the graph.

    Parameters
    ----------
    G : networkx.Graph
        graph

    pos_attr : str, default: 'pos'
        name of the node attribute for position

    Returns
    -------
    euclidean_diameter : float
        Eulidean diameter        
    """
    pos = np.asarray(list(nx.get_node_attributes(G, pos_attr).values()))
    
    try:
        euclidean_diameter = np.max(scipy.spatial.distance.pdist(pos))
    except:
        z = 0.0
        for i in range(len(pos)-1):
            z = max(z, np.sum((pos[i+1:, :] - pos[i, :])**2, axis=1).max())
        euclidean_diameter = np.sqrt(z)

    return euclidean_diameter
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
def outbox_diagonal(G, pos_attr='pos'):
    """
    Computes the length of the diagonal of the outbox including the graph.

    Parameters
    ----------
    G : networkx.Graph
        graph

    pos_attr : str, default: 'pos'
        name of the node attribute for position

    Returns
    -------
    outbox_diagonal : float
        outbox diagonal        
    """
    pos = np.asarray(list(nx.get_node_attributes(G, pos_attr).values()))
    outbox_diagonal = np.sqrt(np.sum(np.ptp(pos, axis=0)**2))

    return outbox_diagonal
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
def pca_of_node_position(G, pos_attr='pos'):
    """
    Performs the principal component analysis of the node position in 3D.

    Parameters
    ----------
    G : networkx.Graph
        graph with 2D or 3D position as node attribute

    pos_attr : str, default: 'pos'
        name of the node attribute for position

    Returns
    -------
    pca_axes : 2d numpy array of shape (3, 3)
        principal axes, in columns, i.e. :
        
        - pca_axes[:, i] : i-th principal axis

        the three principal axes is an orthonormal system, and 
        the axes are ordered such that the variance of the coordinates
        of the node position projected onto these axes are decreasing
    
    pca_var : 1d numpy array of shape (3, )
        variance of the coordinates of the node position projected
        onto each principal axe
    """
    pos = get_pos3d(G, pos_attr=pos_attr)
    pos = np.asarray(list(pos.values()))

    # Covariance matrix of point coordinates
    cov_mat = np.cov(pos.T)

    #  Diagonalization
    pca_var, pca_axes = np.linalg.eig(cov_mat)
    # -> pca_var: (1d-array) eigen values, variances along prinicpal axes
    # -> pca_axes: matrix whose columns are the eigen vectors (of norm 1), i.e. principal axes

    # Sort principal axes according to descending order for d
    ind = np.argsort(pca_var)[::-1]
    pca_var = pca_var[ind]
    pca_axes = pca_axes[:, ind]

    # Set the third principal axis such that det(pca_axes) = 1 (re-orientation if needed)
    pca_axes[:, 2] = np.cross(pca_axes[:, 0], pca_axes[:, 1])
    
    # We have
    #   pca_axes[:, i] : i-th principal axis
    #   pca_var[i]  : variance along the i-th principal axis
    #         
    return pca_axes, pca_var
# -------------------------------------------------------------------------

# ----------------------------------------------------------------------------
def extract_neighborhood(G, u, n=1):
    """
    Extracts a neighborhood of a graph node.

    Parameters
    ----------
    G : networkx.Graph
        input graph
    
    u : label
        label of a node in `G`
    
    n : int, default: 1
        all nodes in `G ` within a distance of (at most) `n`
        (edges) from `u` are considered in the neighborhood of `u`

    Returns
    -------
    Gsub : networkx.Graph
        output graph, subgraph of `G` consisting of the neighborhood of `u`
        only (see parameters above)
    """
    nodes_list = list(nx.single_source_shortest_path_length(G, u, cutoff=n).keys())
    Gsub = nx.subgraph(G, nodes_list)
    return Gsub
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
def set_connected_component_attribute(G, attr_name='connected_component_label'):
    """
    Sets / adds connected component label as node attributes to a graph.

    Labels are integers starting from 0.

    Note that inplace operations are done, i.e. the graph `G` is modified.

    Parameters
    ----------
    G : networkx.Graph
        graph

    attr_name : str, default: 'connected_component_label'
        name of the node attribute for connected component label

    Returns
    -------
    G : networkx.Graph
        updated networkx graph
    """
    node_list = []
    cc_id_list = []

    for i, cc in enumerate(nx.connected_components(G)):
        node_list += list(cc)
        cc_id_list += [i] * len(cc)
        
    cc_id_dict = dict(zip(node_list, cc_id_list))
    nx.set_node_attributes(G, cc_id_dict, attr_name)

    return G 
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
def get_list_of_connected_components(G, order='descending', renumbering_nodes=False):
    """
    Gets the list of connected component(s) of a graph.

    Parameters
    ----------
    G : networkx.Graph
        input graph

    order : str {'ascending', 'descending' (default), 'none'}    
        order of the connected component according to the number of nodes;
        note: any string other than 'ascending' and 'descending' is considered
        as 'none' (no sort)

    renumbering_nodes : bool, default: False
        if `True`, the nodes of every output graph are renumbered (label as 
        integers) starting from 0, otherwise (`False`) original node labels are
        kept

    Returns
    -------
    G_list : list of networkx.Graph
        list of graph(s) corresponding to the connected components of the
        input graph (ordered according to parameter `order`)
    """
    # Get the list of the set of nodes of each connected components
    G_cc_set_of_nodes = list(nx.connected_components(G))
    
    if order == 'ascending':
        isort = np.argsort([len(g_cc_s) for g_cc_s in G_cc_set_of_nodes])
    elif order == 'descending':
        isort = np.argsort([len(g_cc_s) for g_cc_s in G_cc_set_of_nodes])[::-1]
    else:
        isort = np.arange(len(G_cc_set_of_nodes))
    
    G_list = [G.subgraph(G_cc_set_of_nodes[i]) for i in isort]
        
    if renumbering_nodes:
        for i in range(len(G_list)):
            G_list[i] = nx.convert_node_labels_to_integers(G_list[i])

    return G_list
# ----------------------------------------------------------------------------

# def find_neighbors(G,key):   
# list(G.neighbors)
#     return [n for n in G.neighbors(key)]

# -------------------------------------------------------------------------
def replace_edge_key(edges,dict_replacement=None, inverse_dict=None):
    """Replace elements of a list of list based on 
    the corresponding key if element exists in the dictionnary item. 


    Parameters
    ----------
    edges : list of list
        list of edges 
    dict_replacement : dict, optional
        Dictionnary where the keys are the current node ids, and the corresponding item 
        is a list of the previous node ids. 
    inverse_dict : dict, optional
        Dictionnary where the key is the old node names, and the corresponding item is the 
        current node id. The item should be a single value.
    """
    if dict_replacement:
        inverse_dict = { v: k for k, l in dict_replacement.items() for v in l }
    edges_new = []
    for i,edge in enumerate(edges):
        #replace each value by the current ids of the graph
        edges_new.append((inverse_dict[edge[0]],inverse_dict[edge[1]]))

    return edges_new
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
def find_flagged_edges(G,attribute,value):
    """Finds all the edges where a certain value exist in a certain attribute attached to the graph.

    Parameters
    ----------
    G : networkx graph
        networkx graph with attributes on nodes
    attribute : string
        name of the attribute
    value : any
        value to find if present in the edge attribute

    Returns
    -------
    list
        list of tuple of edge
    """
    return [key for key, val in nx.get_edge_attributes(G,attribute).items() if value in val]
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
def find_flagged_nodes(G,attribute,value):
    """Finds all the node ids where a certain value exist in a certain attribute attached to the graph.

    Parameters
    ----------
    G : networkx graph
        networkx graph with attributes on nodes
    attribute : string
        name of the attribute
    value : any
        value to find if present in the node attribute

    Returns
    -------
    list
        list of node ids
    """
    return [key for key, val in nx.get_node_attributes(G,attribute).items() if value in val]
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
def find_value_in_node_attribute(G,attribute, value):
    """Finds all the node ids where a certain value exist in a certain attribute attached to the graph.

    Parameters
    ----------
    G : networkx graph
        networkx graph with attributes on nodes
    attribute : string
        name of the attribute
    value : any
        value to find if present in the node attribute

    Returns
    -------
    list
        list of node ids containing the value
    """
    return [i for i in dict(G.nodes(attribute)) if any(value in s for s in G.nodes(attribute)[i])] 
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
def find_elements(lst1, lst2):
    return list(map(lst1.__getitem__, lst2))
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
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
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
def make_filepath(outputpath,foldername):
    sep = '' if outputpath.endswith('/') else '/'
    filepath = f'{outputpath}{sep}{foldername}' if foldername.endswith('/') else f'{outputpath}{sep}{foldername}/'
    isExist = os.path.exists(filepath)
    if not isExist:
       os.makedirs(filepath)
    return filepath
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
def find_key_from_dict(dictionnary, value):
    #issue when the list is returned empty
    inverse = { v: k for k, l in dictionnary.items() for v in l } 
    return inverse[value]
    # [key for key,values in dictionnary.items() if value in values][0]
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
def find_key_from_fulladdress(G, value):
    fulladdress_dict = nx.get_node_attributes(G,'fulladdress')
    # node_id = [key for key,values in fulladdress_dict.items() if value in values]
    log_key = None
    for key,values in fulladdress_dict.items():

        if value in values:
            #print(value)
            log_key = key
            if log_key is None:
                print(f'{value} is missing in G - fulladdress')
            else:
                return(log_key)
# -------------------------------------------------------------------------
            
# -------------------------------------------------------------------------
# ATTRIBUTE TO DICTIONNARY
def attribute_dict_to_df(G, attribute_name, attribute_type):
    list_data = []
    
    if attribute_type == 'node':
        # attribute_name = str(attribute_name) #this prevent issues with idealized networks, where this value is a set
        attribute_dict = G.nodes(str(attribute_name))
        for key, values in attribute_dict:
            # print(key,values)
            if values is not None:
                #print(key)
                #for pos and width_thickness (list of 2 or 3 floats)
                #########################################################
                if attribute_name=='pos' or attribute_name=='csdim':
                    if type(values)==tuple:
                        list_data.append([key]+list(values))
                    else:
                        list_data.append([key]+values)
                #for splays (list of lists of arrays)
                ###########################################################                
                elif attribute_name=='splays':

                    for array in values: #!!! the toporobot export makes a list of list. simplify this in the toporobot to uniformise
                        list_data.append([key] + array)   
                                    
                    # for list_of_arrays in values:
                    #     print(list_of_arrays)
                    #     for array in list_of_arrays:
                    #         list_data.append([key] + array.tolist())
                        
                #for full_address_id, therion_sql_id, flags, ... (list of n strings or int)
                #########################################################  
                else:
                    #if attribute_name=='flags': print(key, values)
                    if type(values)==str or type(values)==int or type(values)==float: 
                        # print(key,'is string or int or float')
                        list_data.append([key,values])
                    elif type(values)==list:
                        for value in values:                           
                            list_data.append([key,value])
                    
    if attribute_type == 'edge':
        attribute_dict = dict(nx.get_edge_attributes(G,attribute_name))
        
        for key, values in attribute_dict.items():
            if values is not None:
                #if attribute_name=='flags': print(key, values)
                if type(values)==str or type(values)==int or type(values)==float: 
                    # print(key,'is string or int or float')
                    list_data.append(list(key) +[values])
                elif type(values)==list:
                    for value in values:                           
                        list_data.append([key[0],key[1],value])
                # if len(values.split(','))==1:
                #     list_data.append(list(key) +[values])
                # else:
                #     for value in values:
                #         list_data.append([key,value])

    return pd.DataFrame(list_data)
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
def graph_to_branches(G):
    """Breaks down networkx graph into individual branches, by creating new points with unique value index at interesections.
    This was conceived by Nina. for the purpose of the gocad export, which requires to loop through each branch.
    This creates as many connected components as there is branches.
    What it does is that it looks at each intersection (node degree >2) and disconnect (randomly) at the intersection.
    To disconnect, it take randomly a segment, remove it and recreate a segment at the same spot but with a different name
    For example, at an intersection of degree 3, it will only keep 2 segments attached and detach one segement. 

    Args:
        G (networkx Graph): 
            a graph containing values of position, connected_component_number, and intersection

    Returns:
        networkx Graph: separated in branches. 
    """    
    # def add_cc_attribute(G):
    #     """_summary_

    #     Parameters
    #     ----------
    #     G : networkx graph
    #         adds a number as an attribute for all the connected components, from 0 to i (i is the number of connected components)
    #         This is useful for gocad plotting for example
    #     """
    #     id_cc = []
    #     value_cc = []
        
    #     for i,cc in enumerate(nx.connected_components(G)):
    #         id_cc += list(cc)
    #         value_cc += [i] * len(cc)
        
    #     cc_number = dict(zip(id_cc, value_cc))
    #     nx.set_node_attributes(G, cc_number, 'connected_component_number') 

    #!!!here we should find a way to add any attrtibutes to the new create node number!!!
    #and create atoms to link points at branches intersection
    #BREAK THE GRAPH IN SINGLE BRANCHES
    max_value = max(G.nodes) 
    H = G.copy()

    # add connected component number to the graph
    H = set_connected_component_attribute(G, attr_name='connected_component_number')
    # add_cc_attribute(H)

    #loop through all nodes with degree larger than 2
    for intersection in [node for node,degree in dict(H.degree()).items() if degree >2]:
        #loop through all the neighbours minus the first two
        for node in [n for n in H.neighbors(intersection)][2:]:
            #set the new id for the new node
            max_value += 1  
            #create a new node at the intersection
            H.add_node(max_value,   pos=H.nodes('pos')[intersection], 
                                    connected_component_number=H.nodes('connected_component_number')[intersection],
                                    intersection=intersection)  
            #create an edge that connects the new node
            H.add_edge(max_value,node) 
            #remove the old ege
            H.remove_edge(intersection,node)      
    return H
# -------------------------------------------------------------------------


def find_closeby_node(G, dist, node_list_only=False, disconnected_only=False):
    """Find all the pair of nodes within a certain distance of each other.
    Creates a list of pair of nodes and the distance between the two points

    Parameters
    ----------
    G : networkx graph
        cave network graph with 2D or 3D coordinates
    dist : float
        distance in graph units between nodes (smaller or equal)
    node_list_only : bool
        if True, it returns a single list of nodes
        if False, it returns a list of list of pair of nodes and the distance between the two
    disconnected_only : bool
        if True, it only selects the pair of nodes without any connections

    returns
    -------
    list of values [node id 1, node id 2, distance between the two ]

    """
    import itertools
    pos = nx.get_node_attributes(G,'pos')
    list_closeby = []

    for a, b in itertools.combinations(G.nodes(), 2):
        # print(a,b)

        if disconnected_only == True:
            #if the edge does not exist, ignore
            if G.has_edge(a,b)==False:
                #calculate distance between the two edges
                distance = np.linalg.norm(np.array(pos[a])-np.array(pos[b]))
                #select point if the distance is smaller than
                if distance <= dist:
                    #to export a single list of nodes
                    
                    if node_list_only==True:
                        list_closeby.insert(-1,a)
                        list_closeby.insert(-1,b)  
                    #to export a list of list
                    else:
                        list_closeby.append([a,b,float(np.round(distance,2))])

        else:
            #calculate distance between the two edges
            distance = np.linalg.norm(np.array(pos[a])-np.array(pos[b]))
            #select point if the distance is smaller than
            if distance <= dist:
                #to export a single list of nodes
                if node_list_only==True:
                    list_closeby.insert(-1,a)
                    list_closeby.insert(-1,b)  
                #to export a list of list
                else:
                    list_closeby.append([a,b,float(np.round(distance,2))])

    if node_list_only == True:
        list_closeby = np.unique(list_closeby)

    return list_closeby


def set_graph_lengths(G):
    """from base.py
    Calculate edge length and attache it to the graph.
    It updates graph.
    """

    pos = nx.get_node_attributes(G, 'pos')

    # Creation of a dictionary to store the length of each edge
    length = {e:np.sqrt(np.sum((np.asarray(pos[e[1]])-np.asarray(pos[e[0]]))**2))
                for e in G.edges()}
    # length = {}
    # for e in self.graph.edges():
    #     dx = self.pos3d[e[0]][0] - self.pos3d[e[1]][0]
    #     dy = self.pos3d[e[0]][1] - self.pos3d[e[1]][1]
    #     dz = self.pos3d[e[0]][2] - self.pos3d[e[1]][2]
    #     length[e] = np.sqrt(dx ** 2 + dy ** 2 + dz ** 2)
    # Storing the length as an edge attribute
    nx.set_edge_attributes(G, length, 'length')