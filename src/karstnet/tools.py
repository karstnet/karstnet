"""
The module `tools` contains miscellaneous tools / algorithms.
"""


import numpy as np
import networkx as nx
import scipy
import matplotlib.pyplot as plt
import time
import multiprocessing

import geone

import karstnet as kn


# =============================================================================
# Function for simplifying graph
# =============================================================================

# ----------------------------------------------------------------------------
def remove_simple_cycles(
        G, 
        max_cycle_len=None,
        n_iter_max = None,
        node_attr_mode='mean',
        node_weight=None,
        node_attr_list=None,
        return_node_merging_count_dict=False,
        return_node_merging_list_dict=False,
        verbose=0):
    """
    Simplifies a graph by removing simple cycles.

    The function `networkx.simple_cyles` is used to identify the cycles.
    The nodes of any cycle whose length (in number of nodes) not exceeding 
    `max_cylcle_len` (if specified) are removed and replaced by one node, 
    and the edges issued from the removed nodes are removed and replaced 
    by new edges issued from the new node.

    The procedure is iterative: at each iteration, only cycles that
    don't share common nodes are treated.

    At each iteration, a new node has the label of the first node in the 
    cycle it comes from. In particular, the set of output nodes (labels) 
    is a subset of the set of input nodes (labels).

    The node attributes in the output graph are computed according to the 
    attributes of the nodes merged (aggregated) from the input graph, using 
    the mode(s)  specified in `node_attr_mode`.

    The edge attributes are not considered for new edges.

    The function returns an output graph (simplified as described above),
    and, optionally, the dictionary of merged nodes from the input graph,
    and / or the dictionary of the number of merged nodes from the input 
    graph.

    Note: this function should not be used if too many cycles are present
    (function `nx.simple_cycles` can be very time consuming in this case).

    Parameters
    ----------
    G : networkx.Graph
        input graph

    max_cycle_len : int, optional
        maximum cycle length (in number of nodes): cycles with more nodes
        are kept;
        by default (`None`): cycle length not limited
    
    n_iter_max : int, optional
        maximal number of iteration(s)
        by default (`None`): as many iterations as needed to remove all 
        cycles (not exceeding a length of `max_cycle_len`)

    node_attr_mode : str {'first', 'mean'} (or list of strs), default: 'mean'
        string or list of strings, the strings indicate how the node attributes
        (for nodes that are part of a cycle) are computed at each iteration
         of the  in the output graph (see above)

        - if string: the same operation (mode) is used for all considered \
        node attributes
        - if list of strings: the length of the list must be equal to the \
        the list `node_attr_list`, and `node_attr_mode[i]` indicates the \
        operation (mode) used for the node attribute `node_attr_list[i]`

    node_weight : str, optional
        name of the node attribute used as weight to compute the weighted
        mean if `node_attr_mode` is set to 'mean'

    node_attr_list : list of strs, optional
        list of node attributes of the input graph to be included in 
        the output graph;
        by default (`None`): the list of all nodes attributes is considered

    return_node_merging_count_dict : bool, default: False
        if `True`, the dictionary `node_merging_count_dict` is returned
        (see below)

    return_node_merging_list_dict : bool, default: False
        if `True`, the dictionary `node_merging_list_dict` is returned
        (see below)

    verbose : int, default: 0
        verbose mode, larger value implies more info printed

    Returns
    -------
    G_out : networkx.Graph
        output graph (simplifed as described above)

    node_merging_count_dict : dict
        return if `return_node_merging_count_dict=True`; dictionary with

        - key: node label in the input graph
        - value:
            - 0 if the node has been removed in the output graph
            - k > 0 if the node is in the output graph, k being the number of \
            nodes in the input graph that have been merged to form this node \
            in the output graph

    node_merging_list_dict : dict, optional
        return if `return_node_merging_list_dict=True`; dictionary with

        - key: node label in the input graph
        - value: the list of nodes in the input graph that have been merged to \
        form this node in the output graph
    """

    # Check node attributes of the input graph to be considered in the output graph
    # and corresponding operation (mode)

    # Get keys (attribute names) in the input graph
    # # - all keys
    # keys = [list(G.nodes[i].keys()) for i in G.nodes()] # list of lists
    # node_attr_all = np.unique([xij for xi in keys for xij in xi])
    # # - from one node
    # node_attr_all = G.nodes[list(G.nodes())[0]].keys()
    node_attr_all = kn.utils.get_node_attribute_names(G)

    # Check given node attributes name
    if node_attr_list is not None:
        if not np.all([k in node_attr_all for k in node_attr_list]):
            raise ValueError('Attribute name does not exist, check `attr_node_list` parameters')
    else:
        node_attr_list = node_attr_all

    # Check given mode
    if isinstance(node_attr_mode, list):
        if len(node_attr_mode) != len(node_attr_list):
            raise ValueError('Length of the list `node_attr_mode` is not valid')

        if np.any([s not in ('first', 'mean') for s in node_attr_mode]):
            raise ValueError('Entry of the list `node_attr_mode` is not valid')

    else: # `node_attr_mode` assumed to be a string
        if node_attr_mode not in ('first', 'mean'):
            raise ValueError('`node_attr_mode` is not valid')
        node_attr_mode = len(node_attr_list)*[node_attr_mode]

    # Initialize `node_merging_list_dict` (may be used even if not returned)
    node_merging_list_dict = {u:[u] for u in G.nodes()}
    
    # Copy graph to initialize output graph (G_out)
    G_out = G.copy()
    
    if n_iter_max is None:
        n_iter_max = G_out.number_of_nodes()

    for n_iter in range(n_iter_max):
        # Get simple cycles to be removed
        if max_cycle_len is not None:
            cycle_nodes_list = [c for c in nx.simple_cycles(G_out) if len(c) <= max_cycle_len]
        else:
            cycle_nodes_list = [c for c in nx.simple_cycles(G_out)]

        # Remove cycles so that two cycles do not have common nodes
        i0 = 0
        while i0 < len(cycle_nodes_list) - 1:
            ind = np.ones(len(cycle_nodes_list), dtype='bool')
            for i in range(i0+1, len(cycle_nodes_list)):
                for x in cycle_nodes_list[i]:
                    if x in cycle_nodes_list[i0]:
                        ind[i] = False
                        break
            cycle_nodes_list = [cycle_nodes_list[i] for i in np.where(ind)[0]]
            i0 = i0+1
            
        # Number of new nodes
        n_new_nodes = len(cycle_nodes_list)

        if verbose > 0:
            print(f'Iteration {n_iter+1} : number of cycles to simplify = {n_new_nodes}')

        if n_new_nodes == 0:
            break # exit for loop

        # Update `node_merging_list_dict`
        for c in cycle_nodes_list:
            for u in c[1:]:
                node_merging_list_dict[c[0]].extend(node_merging_list_dict[u])
                node_merging_list_dict[u] = []

        # Set list of nodes to remove (all except first node in each cycle)
        nodes_list_to_remove = [u for c in cycle_nodes_list for u in c[1:]]

        # Set new edges according to the edges involving nodes that will be removed
        # - Set dictionary with
        #   - key: current node label
        #   - value: list of indices of the cycles containing the node
        node2cycle = {u:[] for u in G_out.nodes()}
        for i, c in enumerate(cycle_nodes_list):
            for u in c:
                node2cycle[u].append(i)
        # -> node u s.t. node2cycle[u] is a non empty list will be removed

        new_edges = []
        for i, c in enumerate(cycle_nodes_list):
            for u in c[1:]:
                edges = []
                for v in G_out.neighbors(u):
                    if len(node2cycle[v]) == 0:
                        edges.append((c[0], v))
                    else:
                        for j in node2cycle[v]:
                            if j != i:
                                edges.append((c[0], cycle_nodes_list[j][0]))
                for e in edges:
                    if [e[0], e[1]] not in new_edges and [e[1], e[0]] not in new_edges:
                        new_edges.append(e)

        # Remove nodes in cycle
        G_out.remove_nodes_from(nodes_list_to_remove)

        # Add new edges
        G_out.add_edges_from(new_edges)

        # --- End for loop iter ---
                
    # Set node attributes (according to `node_attr_list` and `node_attr_mode`)
    for attr, mode in zip(node_attr_list, node_attr_mode):
        d = nx.get_node_attributes(G, attr, default=np.nan) # dictionary node_label:value_of_attribute in input graph
        if mode == 'first':
            for u, u_list in node_merging_list_dict.items():
                if len(u_list) <= 1:
                    continue
                G_out.nodes[u][attr] = d[u_list[0]][attr]

        elif mode == 'mean':
            v0 = list(d.values())[0] # value of one node in G
            if hasattr(v0, '__len__'):
                if node_weight is None:
                    for u, u_list in node_merging_list_dict.items():
                        if len(u_list) <= 1:
                            continue
                        G_out.nodes[u][attr] = np.nanmean(np.asarray([np.atleast_1d(d[v]) for v in u_list]), axis=0)
                else:
                    for u, u_list in node_merging_list_dict.items():
                        if len(u_list) <= 1:
                            continue
                        values = np.asarray([np.atleast_1d(d[v]) for v in u_list])
                        w = np.asarray([G.nodes[v][node_weight] for v in u_list]).reshape(-1, 1)
                        ind = np.all(~np.isnan(values), axis=1)
                        G_out.nodes[u][attr] = np.sum(values[ind]*w[ind], axis=0)/np.sum(w[ind])

            else:
                if node_weight is None:
                    for u, u_list in node_merging_list_dict.items():
                        if len(u_list) <= 1:
                            continue
                        G_out.nodes[u][attr] = np.nanmean(np.asarray([np.atleast_1d(d[v]) for v in u_list]), axis=0)[0]
                else:
                    for u, u_list in node_merging_list_dict.items():
                        if len(u_list) <= 1:
                            continue
                        values = np.asarray([np.atleast_1d(d[v]) for v in u_list])
                        w = np.asarray([G.nodes[v][node_weight] for v in u_list]).reshape(-1, 1)
                        ind = np.all(~np.isnan(values), axis=1)
                        G_out.nodes[u][attr] = (np.sum(values[ind]*w[ind], axis=0)/np.sum(w[ind]))[0]
 
    out = [G_out]
    if return_node_merging_count_dict:
        node_merging_count_dict = {u:len(u_list) for u, u_list in node_merging_list_dict.items()}
        out.append(node_merging_count_dict)

    if return_node_merging_list_dict:
        out.append(node_merging_list_dict)
    
    if len(out) == 1:
        out = out[0]
    else:
        out = tuple(out)
    
    return out
# ----------------------------------------------------------------------------


# =============================================================================
# Variogram on graph
# =============================================================================

# # -------------------------------------------------------------------------
# def variogram_cloud(
#         G, 
#         node_attr,
#         nsample_pair,
#         edge_length_attr=None,
#         exclude_nan=True,
#         hmax=None,
#         seed=None,
#         make_plot=True,
#         verbose=0,
#         **kwargs):
#     """
#     Computes the variogram cloud for the given node attribute `node_attr` 
#     which is assumed to be a scalar (float) or a sequence of floats.

#     For each sample pair of nodes i, j, let
    
#     - h(i, j) the shortest path length between nodes i and j
#     - gamma(i, j) = 0.5*(v[i]-v[j])**2 (gamma values for the attribute) \
#     where v[k] is the scalar (or vector) attribute attached to node k

#     This method retrieves all points (h(i, j), gamma(i, j)) for which h(i, j) is
#     less than or equal to `hmax` (if given), among the `nsample_pair` pairs
#     of nodes randomly chosen in the entire graph `G`.

#     Parameters
#     ----------
#     G : networkx.Graph
#         input graph

#     node_attr : str
#         name of the node attribute for which the variogram cloud is computed;
#         this node attribute should be a scalar (float) or a sequence of floats

#     nsample_pair : int
#         number of sample pairs of nodes (randomly chosen) in the original graph

#     edge_length_attr : str, optional
#         name of the edge attribute for length (used to compute the distance between
#         nodes);
#         by default (`None`): the edges have a length of one

#     hmax : float, optional
#         maximal distance between two nodes to be integrated in the variogram cloud

#     exclude_nan : bool, default: `True`
#         if `True`, any pair of nodes i, j within a branch, with at least 
#         one `nan` value for the considered attribute at these nodes, is 
#         excluded from the variogram cloud

#     seed : int, optional
#         seed to initialize the random number generator

#     make_plot : bool, default: `True`
#         indicates if the variogram cloud is plotted (in the current
#         figure axis)

#     verbose : int, default: 0
#         verbose mode, larger value implies more info printed

#     kwargs : dict
#         keyword arguments passed to the function `matplotlib.pyplot.plot`
#         (used if `make_plot=True`)

#     Returns
#     -------
#     h : 1d numpy array of floats
#         distance between pair of nodes (1st coordinate) of the variogram cloud points

#     gamma : 1d or 2d numpy array of floats
#         gamma values for the specified attribute of the variogram
#         cloud points:
        
#         - if the specified attribute (`node_attr`) is a scalar, \
#         `gamma` is a 1d array and the points (h, gamma) constitutes the variogram cloud \
#         for that attribute (may contains `nan` values if `exclude_nan=False`)
#         - if the specified attribute (`node_attr`) is a sequence, \
#         `gamma` is a 2d array and the points (h, gamma[:, i]) constitutes the variogram cloud \
#         for the component i of that attribute (may contains `nan` values if `exclude_nan=False`)
    
#     npair : int
#         number of points (pairs of data points considered) in the variogram cloud 
#         (length of `h` or `gamma`)

#     Examples
#     --------
#     >>> h, gamma, npair = myKGraph.variogram_cloud('node_prop', nsample_pair=1000, seed=235)
#     """

#     # Check node attribute
#     if node_attr not in kn.utils.get_node_attribute_names(G):
#         raise ValueError(f'{node_attr} is not a node attribute')

#     # Check edge length attribute
#     if edge_length_attr is not None and edge_length_attr not in kn.utils.get_edge_attribute_names(G):
#         raise ValueError(f'{edge_length_attr} is not an edge attribute')

#     # Set dictionary to convert node label (id) to node index, and vice versa
#     #node_label2index = {u:i for i, u in enumerate(G.nodes())}
#     node_index2label = {i:u for i, u in enumerate(G.nodes())}

#     # Get the 2d-array of property values, of shape (n_nodes, n_prop)
#     # where n_nodes is the number of nodes in the graph and 
#     # n_prop is the number of properties (1 if scalar property)
#     n_nodes = G.number_of_nodes()

#     prop_dict = nx.get_node_attributes(G, node_attr)
#     if len(prop_dict) == 0:
#         raise ValueError(f'No value for specified node attribute ({node_attr})')

#     v_first = np.asarray(list(prop_dict.values())[0]) # first value entry as array
#     if v_first.ndim > 1:
#         raise ValueError('Value of the specified attribute should by a scalar or a sequence of scalar (value is an array of dimension greater than 1)')

#     n_prop = np.asarray(v_first).size # size of first value entry
#     n_nodes = G.number_of_nodes() # total number of nodes

#     if len(prop_dict) != n_nodes:
#         # There is missing value for some nodes, set default (np.nan)
#         if v_first.ndim == 0:
#             # attribute is a scalar (n_prop = 1)
#             prop_dict = nx.get_node_attributes(G, node_attr, default=np.nan)
#         else: # v_first.ndim == 1
#             # attribute is a sequence
#             prop_dict = nx.get_node_attributes(G, node_attr, default=np.full(n_prop, np.nan))

#     try:
#         v = np.asarray(list(prop_dict.values())).reshape(G.number_of_nodes(), -1)
#     except:
#         raise ValueError('Unable to set specified attribute (entry for every node should have the same shape)')

#     if np.all(np.isnan(v), axis=0).any():
#         raise ValueError('The specified attribute (or one of its component) has only undefined (`nan`) values')


#     # Initialize random number generator
#     if seed is not None:
#         np.random.seed(seed)

#     # Compute variogram cloud
#     h, gamma = [], []

#     if verbose > 0:
#         progress_old = 0

#     for k in range(nsample_pair):
#         if verbose > 0:
#             progress = int(k/nsample_pair*100.0)
#             if progress > progress_old:
#                 print(f'Progress: {progress:3d}%')
#                 progress_old = progress
                
#         i, j = np.random.choice(n_nodes, 2, replace=False)
        
#         if exclude_nan and (np.isnan(v[i]).any() or np.isnan(v[j]).any()):
#             continue

#         ui, uj = node_index2label[i], node_index2label[j]

#         try:
#             # path = nx.shortest_path(G, source=ui, target=uj, weight=edge_length_attr, method='dijkstra')
#             # hk = float(np.sum(np.asarray([G.edges[(a, b)][edge_length_attr] for a, b in zip(path[:-1], path[1:])])))
#             hk = nx.shortest_path_length(G, source=ui, target=uj, weight=edge_length_attr, method='dijkstra')
#         except:
#             # shortest_path_length may "fail" if there is no path between u and v (e.g. if graph G has several connected components)
#             continue
        
#         if hmax is not None and hk > hmax:
#             continue

#         # gammak = 0.5*(G.nodes[u][node_attr] - G.nodes[v][node_attr])**2
#         gammak = 0.5*(v[i] - v[j])**2

#         h.append(hk)
#         gamma.append(gammak)

#     if verbose > 0:
#         progress = 100
#         print(f'Progress: {progress:3d}%')

#     h = np.asarray(h)
#     gamma = np.asarray(gamma)

#     if n_prop == 1:
#         gamma = gamma.reshape(-1)

#     if make_plot:
#         # default linestyle and marker for plot (if not specified)
#         if 'linestyle' not in kwargs.keys() and 'ls' not in kwargs.keys():
#             kwargs['linestyle'] = ''
#         if 'marker' not in kwargs.keys():
#             kwargs['marker'] = '.'
#         if n_prop == 1:
#             plt.plot(h, gamma, **kwargs)
#         else: 
#             plt.plot(h, gamma, label=[f'index {i}' for i in range(n_prop)], **kwargs)
#             plt.legend()
#         plt.xlabel('h (distance btween pair of nodes)')
#         # plt.xlabel('h')
#         # plt.ylabel(r'$1/2(Z(x)-Z(x+h))^2$ (with $Z$ considered variable)')
#         plt.ylabel(r'$1/2(Z(x)-Z(x+h))^2$')

#     return h, gamma, len(h)
# # ----------------------------------------------------------------------------

# -------------------------------------------------------------------------
def variogram_cloud(
        G, 
        node_attr,
        nsample_start_nodes,
        edge_length_attr=None,
        exclude_nan=True,
        hmax=None,
        seed=None,
        make_plot=True,
        verbose=0,
        **kwargs):
    """
    Computes the variogram cloud for the given node attribute `node_attr` 
    which is assumed to be a scalar (float) or a sequence of floats.

    For each sample pair of nodes i, j, let
    
    - h(i, j) the shortest path length between nodes i and j
    - gamma(i, j) = 0.5*(v[i]-v[j])**2 (gamma values for the attribute) \
    where v[k] is the scalar (or vector) attribute attached to node k

    This method retrieves all points (h(i, j), gamma(i, j)) for which h(i, j) is
    less than or equal to `hmax` (if given), among the `nsample_pair` pairs
    of nodes randomly chosen in the entire graph `G`.

    Parameters
    ----------
    G : networkx.Graph
        input graph

    node_attr : str
        name of the node attribute for which the variogram cloud is computed;
        this node attribute should be a scalar (float) or a sequence of floats

    nsample_start_nodes : int
        number of sample start nodes (randomly chosen) in the original graph; for
        each sample start node i, all nodes j (not equal to i), at distance less 
        than or equal to `hmax` (if specified), are considered for a pair (i, j)
        to be included in the variogram cloud

    edge_length_attr : str, optional
        name of the edge attribute for length (used to compute the distance between
        nodes);
        by default (`None`): the edges have a length of one

    hmax : float, optional
        maximal distance between two nodes to be integrated in the variogram cloud

    exclude_nan : bool, default: `True`
        if `True`, any pair of nodes i, j within a branch, with at least 
        one `nan` value for the considered attribute at these nodes, is 
        excluded from the variogram cloud

    seed : int, optional
        seed to initialize the random number generator

    make_plot : bool, default: `True`
        indicates if the variogram cloud is plotted (in the current
        figure axis)

    verbose : int, default: 0
        verbose mode, larger value implies more info printed

    kwargs : dict
        keyword arguments passed to the function `matplotlib.pyplot.plot`
        (used if `make_plot=True`)

    Returns
    -------
    h : 1d numpy array of floats
        distance between pair of nodes (1st coordinate) of the variogram cloud points

    gamma : 1d or 2d numpy array of floats
        gamma values for the specified attribute of the variogram
        cloud points:
        
        - if the specified attribute (`node_attr`) is a scalar, \
        `gamma` is a 1d array and the points (h, gamma) constitutes the variogram cloud \
        for that attribute (may contains `nan` values if `exclude_nan=False`)
        - if the specified attribute (`node_attr`) is a sequence, \
        `gamma` is a 2d array and the points (h, gamma[:, i]) constitutes the variogram cloud \
        for the component i of that attribute (may contains `nan` values if `exclude_nan=False`)
    
    npair : int
        number of points (pairs of data points considered) in the variogram cloud 
        (length of `h` or `gamma`)

    Examples
    --------
    >>> h, gamma, npair = myKGraph.variogram_cloud('node_prop', nsample_pair=1000, seed=235)
    """

    # Check node attribute
    if node_attr not in kn.utils.get_node_attribute_names(G):
        raise ValueError(f'{node_attr} is not a node attribute')

    # Check edge length attribute
    if edge_length_attr is not None and edge_length_attr not in kn.utils.get_edge_attribute_names(G):
        raise ValueError(f'{edge_length_attr} is not an edge attribute')

    # Set dictionary to convert node label (id) to node index, and vice versa
    node_label2index = {u:i for i, u in enumerate(G.nodes())}
    node_index2label = {i:u for i, u in enumerate(G.nodes())}

    # Get the 2d-array of property values, of shape (n_nodes, n_prop)
    # where n_nodes is the number of nodes in the graph and 
    # n_prop is the number of properties (1 if scalar property)
    n_nodes = G.number_of_nodes()

    prop_dict = nx.get_node_attributes(G, node_attr)
    if len(prop_dict) == 0:
        raise ValueError(f'No value for specified node attribute ({node_attr})')

    v_first = np.asarray(list(prop_dict.values())[0]) # first value entry as array
    if v_first.ndim > 1:
        raise ValueError('Value of the specified attribute should by a scalar or a sequence of scalar (value is an array of dimension greater than 1)')

    n_prop = np.asarray(v_first).size # size of first value entry
    n_nodes = G.number_of_nodes() # total number of nodes

    if len(prop_dict) != n_nodes:
        # There is missing value for some nodes, set default (np.nan)
        if v_first.ndim == 0:
            # attribute is a scalar (n_prop = 1)
            prop_dict = nx.get_node_attributes(G, node_attr, default=np.nan)
        else: # v_first.ndim == 1
            # attribute is a sequence
            prop_dict = nx.get_node_attributes(G, node_attr, default=np.full(n_prop, np.nan))

    try:
        v = np.asarray(list(prop_dict.values())).reshape(G.number_of_nodes(), -1)
    except:
        raise ValueError('Unable to set specified attribute (entry for every node should have the same shape)')

    if np.all(np.isnan(v), axis=0).any():
        raise ValueError('The specified attribute (or one of its component) has only undefined (`nan`) values')


    # Initialize random number generator
    if seed is not None:
        np.random.seed(seed)

    # Compute variogram cloud
    h, gamma = [], []

    if verbose > 0:
        progress_old = 0

    nsample_start_nodes = min(nsample_start_nodes, n_nodes)
    start_node_ind = np.random.choice(n_nodes, nsample_start_nodes, replace=False)

    for k, i in enumerate(start_node_ind):
        if verbose > 0:
            progress = int(k/nsample_start_nodes*100.0)
            if progress > progress_old:
                print(f'Progress: {progress:3d}%')
                progress_old = progress
                
        vi = v[i]
        if exclude_nan and np.isnan(vi).any():
            continue

        ui = node_index2label[i]

        length_dict = nx.single_source_dijkstra_path_length(G, ui, cutoff=hmax, weight=edge_length_attr)
        for uj, hj in length_dict.items():
            vj = v[node_label2index[uj]]
            if exclude_nan and np.isnan(vj).any():
                continue

            h_ij = hj 
            # gamma_ij = 0.5*(G.nodes[ui][node_attr] - G.nodes[uj][node_attr])**2
            gamma_ij = 0.5*(vi - vj)**2

            h.append(h_ij)
            gamma.append(gamma_ij)

    if verbose > 0:
        progress = 100
        print(f'Progress: {progress:3d}%')

    h = np.asarray(h)
    gamma = np.asarray(gamma)

    if n_prop == 1:
        gamma = gamma.reshape(-1)

    if make_plot:
        # default linestyle and marker for plot (if not specified)
        if 'linestyle' not in kwargs.keys() and 'ls' not in kwargs.keys():
            kwargs['linestyle'] = ''
        if 'marker' not in kwargs.keys():
            kwargs['marker'] = '.'
        if n_prop == 1:
            plt.plot(h, gamma, **kwargs)
        else: 
            plt.plot(h, gamma, label=[f'index {i}' for i in range(n_prop)], **kwargs)
            plt.legend()
        plt.xlabel('h (distance btween pair of nodes)')
        # plt.xlabel('h')
        # plt.ylabel(r'$1/2(Z(x)-Z(x+h))^2$ (with $Z$ considered variable)')
        plt.ylabel(r'$1/2(Z(x)-Z(x+h))^2$')

    return h, gamma, len(h)
# ----------------------------------------------------------------------------

# -------------------------------------------------------------------------
def variogram_exp_from_variogram_cloud(
        variogramCloud,
        hmax=None,
        ncla=10,
        cla_center=None,
        cla_length=None,
        make_plot=True,
        show_count=True,
        **kwargs):
    """
    Computes the experimental variogram from a given variogram cloud.

    For each pair of nodes i, j within a branch, let
        
    - h(i, j) the length between nodes i and j
    - gamma(i, j) = 0.5*(v[i]-v[j])**2 (gamma values for the property) \
    where v[k] is the property scalar (or vector) of attached to node k

    The variogram cloud consists of all points (h(i, j), gamma(i, j)), with
    h the distance and gamma the gamma value, for which h(i, j) is less than or 
    equal to `hmax` (if given) (see function `variogram_cloud`).

    The mean point in each class is retrieved from the variogram cloud; the
    i-th class is determined by its center `cla_center[i]` and its length
    `cla_length[i]`, and corresponds to the interval
    `]cla_center[i]-cla_length[i]/2, cla_center[i]+cla_length[i]/2]`
    along h (lag) axis (abscissa).

    Parameters
    ----------
    variogramCloud : 3-tuple
        `variogramCloud`=(h, gamma, npair) is a variogram cloud (already computed,
        e.g. returned by the function `variogram_cloud`)

    hmax : float, optional
        maximal distance between two nodes to be integrated in the 
        experimental variogram

    ncla : int, default: 10
        number of classes, the parameter is used if `cla_center=None`, in that
        situation `ncla` classes are considered and the class centers are set to
        
        - `cla_center[i] = (i+0.5)*l, i=0,...,ncla-1`
        
        with l = H / ncla, H being the max of the distance between two points of
        the considered pairs (in the variogram cloud);
        if `cla_center` is specified (not `None`), the number of classes (`ncla`)
        is set to the length of the sequence `cla_center` (ignoring the value
        passed as argument)
    
    cla_center : 1D array-like of floats, optional
        sequence of floats, center of each class (in abscissa) in the experimental
        variogram; by default (`None`): `cla_center` is defined from `ncla` (see
        above)

    cla_length : 1D array-like of floats, or float, optional
        length of each class centered at `cla_center` (in abscissa) in the
        experimental variogram:

        - if `cla_length` is a sequence, it should be of length `ncla`
        - if `cla_length` is a float, the value is repeated `ncla` times
        - if `cla_length=None` (default), the minimum of difference between two \
        sucessive class centers (`np.inf` if one class) is used and repeated `ncla` \
        times
    

    make_plot : bool, default: `True`
        indicates if the experimental variogram is plotted (in the current
        figure axis)

    show_count : bool, default: `True`
        indicates if counters (`cexp`) are displayed on the plot
        (used if `make_plot=True`)

    kwargs : dict
        keyword arguments passed to the function `matplotlib.pyplot.plot`
        (used if `make_plot=True`)

    Returns
    -------
    hexp : 1d or 2d numpy array of floats
        distance (1st coordinate) of the experimental variogram

        - if experimental variogram is computed for one scalar variable \
        (i.e. if `variogramCloud[1]` is a 1d array), then `hexp` is a 1d array
        - if experimental variogram is computed for more than one scalar variable \
        (i.e. if `variogramCloud[1]` is a 2d array), then `hexp` is a 2d array \
        where each column corresponds to each variable 

    gexp : 1d or 2d numpy array of floats
        gamma values (2nd coordinate) of the experimental variogram

        - the shape of the array `gexp` is the same as the shape of `hexp` \
        (see `hexp` above)

    cexp : 1d or 2d numpy array of ints
        counters, i.e. number of points (pairs of data points considered) 
        in each class in the variogram cloud

        - the shape of the array `cexp` is the same as the shape of `hexp` \
        (see `hexp` above)
    """

    h, gamma, npair = variogramCloud

    if hmax is not None:
        ind = h <= hmax
        h = h[ind]
        gamma = gamma[ind]
        npair = len(h)

    if npair == 0:
        print('No point in the variogram cloud (nothing is done).')
        return None, None, None

    # Set classes
    if cla_center is not None:
        cla_center = np.asarray(cla_center, dtype='float').reshape(-1)
        ncla = len(cla_center)
    else:
        length = np.max(h) / ncla
        cla_center = (np.arange(ncla, dtype='float') + 0.5) * length

    if cla_length is not None:
        cla_length = np.asarray(cla_length, dtype='float').reshape(-1)
        if len(cla_length) == 1:
            cla_length = np.repeat(cla_length, ncla)
        elif len(cla_length) != ncla:
            print(f"ERROR: `cla_length` not valid")
            return None, None, None
    else:
        if ncla == 1:
            cla_length = np.array([np.inf], dtype='float')
        else:
            cla_length = np.repeat(np.min(np.diff(cla_center)), ncla)

    # Set gamma as a 2d array
    gamma = gamma.reshape(npair, -1)

    n_prop = gamma.shape[1]

    # Compute experimental variogram
    hexp = np.full((ncla, n_prop), np.nan)
    gexp = np.full((ncla, n_prop), np.nan)
    cexp = np.zeros((ncla, n_prop), dtype='int')

    for i, (c, l) in enumerate(zip(cla_center, cla_length)):
        d = 0.5*l
        ind = np.all((h > c-d , h <= c+d), axis=0)
        ind_not_nan = ~np.isnan(gamma[ind])
        cexp[i] = np.sum(ind_not_nan, axis=0)
        hexp[i] = np.sum(h[ind].reshape(len(h[ind]), 1)*ind_not_nan, axis=0) / cexp[i]
        gexp[i] = np.nanmean(gamma[ind], axis=0)

    # If a scalar property is considered, reshape hexp, gexp and cexp
    if n_prop == 1:
        hexp = hexp.reshape(-1)
        gexp = gexp.reshape(-1)
        cexp = cexp.reshape(-1)

    if make_plot:
        # default linestyle and marker for plot (if not specified)
        if 'linestyle' not in kwargs.keys() and 'ls' not in kwargs.keys():
            kwargs['linestyle'] = 'solid'
        if 'marker' not in kwargs.keys():
            kwargs['marker'] = '+'

        if n_prop == 1:
            plt.plot(hexp, gexp, **kwargs)
            if show_count:
                for i, c in enumerate(cexp):
                    if c > 0:
                        plt.text(hexp[i], gexp[i], str(c), ha='left', va='top')
        else:
            if 'label' in kwargs.keys():
                if kwargs['label'] is not None:
                    label = f"{kwargs['label']} - "
                else:
                    label = 'None'
                del kwargs['label']
            else:
                label = None

            for j in range(n_prop):
                if label is not None:
                    label = f'{label} - index {j}'
                plt.plot(hexp[:, j], gexp[:, j], label=label, **kwargs)
                if show_count:
                    for i, c in enumerate(cexp[:, j]):
                        if c > 0:
                            plt.text(hexp[i, j], gexp[i, j], str(c), ha='left', va='top')
                plt.legend()
        plt.xlabel('h (distance btween pair of nodes)')
        # plt.xlabel('h')
        # plt.ylabel(r'$1/2(Z(x)-Z(x+h))^2$ (with $Z$ considered variable)')
        plt.ylabel(r'$1/2(Z(x)-Z(x+h))^2$')
    
    return hexp, gexp, cexp
# -------------------------------------------------------------------------

# =============================================================================
# Kriging and Sequential Gaussian Simulation (SGS) on graph
# =============================================================================

# ----------------------------------------------------------------------------
def krige_node_attribute(
        G,
        cov_model,
        node_attr,
        err_std=0.0,
        edge_length_attr=None,
        method='ordinary_kriging',
        mean=None,
        var=None,
        neighbors_on_same_branch_only=False,
        branches=None,
        use_unique_neighborhood=False,
        searchRadius=None,
        searchRadiusRelative=1.2,
        nneighborMax=12,
        update_graph=True,
        kriging_est_node_attr=None,
        kriging_std_node_attr=None,
        i0=None, 
        i1=None,
        pid=None,
        verbose=1):
    """
    Interpolates a node attribute on a graph by kriging.

    Parameters
    ----------
    G : networkx.Graph
        input graph

    cov_model : :class:`geone.covModel.CovModel1D`
        covariance model in 1D

    node_attr : str
        name of the node attribute to be kriged

    err_std : float, or str, default : 0.0
        standard deviation of error (zero-mean Gaussian of given std):
        - if float: the same error std is used for data at any graph node; 
        - if str: name of the node attribute of the error std

    edge_length_attr : str, optional
        name of the edge attribute for length (used to compute the distance between
        nodes);
        by default (`None`): the edges have a length of one

    method : str {'simple_kriging', 'ordinary_kriging'}, default: 'ordinary_kriging'
        type of kriging;
        note: if `method='ordinary_kriging'`, the parameter `mean` is not used

    mean : float, or str, optional
        kriging mean value :
        - if float: the value is used at all graph nodes; 
        - if str: name of the node attribute for the mean 
        
        by default (`None`), the mean of data node values is used for any node
        
        note: if `method=ordinary_kriging`, parameter `mean` is ignored

    var : float, or str, optional
        kriging variance value :
        - if float: the value is used at all graph nodes; 
        - if str: name of the node attribute for the mean 

        note: if `method=ordinary_kriging`, parameter `var` is ignored

    neighbors_on_same_branch_only : bool, default: False
        - if `True`: only neighbors on the branch of the estimated graph node are used
        - if `False`: neighbors everywhere in the graph are used

        note: `neighbors_on_same_branch_only=True` requires `branches` (see below)
        
    branches : list, optional
        list of branches, required if `neighbors_on_same_branch_only=True`;
        each element (branch) is a list of node labels, typically `branches=kg.branches`
        should be used, where `kg` is the parent `KGraph` of `G`
        
    use_unique_neighborhood : bool, default: False
        indicates if a unique neighborhood is used:

        - if True: data at any graph node is taken into account, and the kriging matrix \
        is computed once; the parameters `searchRadius`, `searchRadiusRelative`, \
        `nneighborMax` are not used
        - if False: only data at graph node within a search neighborhood are \
        taken into account according to `searchRadius`, `searchRadiusRelative`, `nneighborMax`

        note: `use_unique_neighborhood=True` is not compatible with `neighbors_on_same_branch_only=True`

    searchRadius : float, optional
        if specified, i.e. not `None`: search radius, i.e. 
        the data at graph node at distance to the estimated graph node greater 
        than `searchRadius` are not taken into account in the kriging system; 
        if `searchRadius` is not `None`, then `searchRadiusRelative` is not used;
        by default (`searchRadius=None`): `searchRadiusRelative` is used;

    searchRadiusRelative : float, default: 1.2
        used only if `searchRadius` is `None`;
        the search radius is set to `searchRadiusRelative` times the range of the 
        covariance model `cov_model`

    nneighborMax : int, default: 12
        maximal number of neighbors (data at graph nodes) taken into account in the
        kriging system; the data at graph nodes the closest to the estimated graph node are
        taken into account;
        note: if `nneighborMax=None` or `nneighborMax<0`, then `nneighborMax` is
        set to the number of graph nodes with informed data

    update_graph : bool, default True
        - if `True`: the kriging estimates is set as a node attribute (named \
        `kriging_est_node_attr` (see below)) and the kriging standard deviation \
        is set as a node attribute (named `kriging_std_node_attr` (see below))
        - if `False`: the graph is not updated (`kriging_est_node_attr` and 
        `kriging_std_node_attr` are not used)

    kriging_est_node_attr : str, optional
        name of the node attribute in output for kriging estimates;
        by default (`None`),  the name `node_attr` + '_krig_est' is used

    kriging_std_node_attr : str, optional
        name of the node attribute in output for kriging standard deviation;
        by default (`None`),  the name `node_attr` + '_krig_std' is used

    i0 : int, optional
        starting index in `list(G.nodes())` of nodes to be kriged (used with multiprocessing)

    i1 : int, optional
        ending index (excluded) in `list(G.nodes())` of nodes to be kriged (used with multiprocessing)

    pid : int, optional
        process id of the caller (used with multiprocessing)

    verbose : int, default: 0
        verbose mode, higher implies more printing (info)

    Returns
    -------
    kriging_est : dict
        dictionary with node labels as keys and kriging estimates as values

    kriging_std : dict
        dictionary with node labels as keys and kriging standard deviation as values
    """
    if verbose > 0:
        if pid is not None:
            pid_str = f'[pid={pid}] '
        else:
            pid_str = ''

    # Check cov_model
    if not isinstance(cov_model, geone.covModel.CovModel1D):
        raise ValueError(f'{pid_str}`cov_model` must be an instance of `geone.covModel.CovModel1D`')

    # Prevent calculation if covariance model is not stationary
    if not cov_model.is_stationary():
        raise ValueError(f'{pid_str}`cov_model` is not stationary')

    # Check node attribute
    if node_attr not in kn.utils.get_node_attribute_names(G):
        raise ValueError(f'{pid_str}{node_attr} is not a node attribute')

    # Check edge length attribute
    if edge_length_attr is not None and edge_length_attr not in kn.utils.get_edge_attribute_names(G):
        raise ValueError(f'{pid_str}{edge_length_attr} is not an edge attribute')

    # Check neighbors_on_same_branch_only
    if neighbors_on_same_branch_only:
        if branches is None:
            raise ValueError(f'{pid_str}`branches` is required with `neighbors_on_same_branch_only=True`')

        if use_unique_neighborhood:
            raise ValueError(f'{pid_str}`use_unique_neighborhood=True` is incompatible with `neighbors_on_same_branch_only=True`')
        
    # # Set dictionary to convert node label (id) to node index, and vice versa
    # node_label2index = {u:i for i, u in enumerate(G.nodes())}
    # #node_index2label = {i:u for i, u in enumerate(G.nodes())}

    if update_graph:
        # Set node attribute names for output
        if kriging_est_node_attr is not None:
            if not isinstance(kriging_est_node_attr, str):
                raise ValueError(f'{pid_str}`kriging_est_node_attr` must be a string (node attribute name)')
        else:
            kriging_est_node_attr = f'{node_attr}_krig_est'

        if kriging_std_node_attr is not None:
            if not isinstance(kriging_std_node_attr, str):
                raise ValueError(f'{pid_str}`kriging_std_node_attr` must be a string (node attribute name)')
        else:
            kriging_std_node_attr = f'{node_attr}_krig_std'

    # Get the dictionary of node attribute (property) to be kriged
    prop_dict = nx.get_node_attributes(G, node_attr)
    if len(prop_dict) == 0:
        raise ValueError(f'{pid_str}No value for specified node attribute ({node_attr})')

    v_first = np.asarray(list(prop_dict.values())[0]) # first value entry as array
    if v_first.ndim > 0:
        raise ValueError(f'{pid_str}Value of the specified attribute should by a scalar')

    # - set nan at uninformed nodes
    prop_dict = nx.get_node_attributes(G, node_attr, default=np.nan)

    # Covariance function and value at 0
    cov_func = cov_model.func() # covariance function
    cov0 = cov_func(0.)[0] # covariance function at origin (lag=0)

    # Default mean value    
    tmp = np.asarray(list(prop_dict.values()))
    if np.isnan(tmp).all():
        mean_default = 0.0
    else:
        # Set mean of data values
        mean_default = np.nanmean(tmp)
    
    # Method and mean, var
    if method == 'simple_kriging':
        ordinary_kriging = False
        # Get the dictionary of mean kriging values
        if mean is not None:
            if isinstance(mean, float) or isinstance(mean, int):
                mean_dict = {u:mean for u in G.nodes()}

            elif isinstance(mean, str):
                mean_dict = nx.get_node_attributes(G, mean, default=np.nan)
                if np.isnan(np.asarray(list(mean_dict.values()))).any():
                    raise ValueError(f'{pid_str}Specified mean is not defined at all graph nodes')
            
            else:
                raise ValueError(f'{pid_str}Specified mean is not valid')

        else:
            mean_dict = {u:mean_default for u in G.nodes()}
    
        # Get the dictionary of "variance update"
        if var is not None:
            if isinstance(var, float) or isinstance(var, int):
                tmp = np.sqrt(var/cov0)
                var_update_dict = {u:tmp for u in G.nodes()}

            elif isinstance(var, str):
                var_update_dict = nx.get_node_attributes(G, var, default=np.nan)
                if np.isnan(np.asarray(list(var_update_dict.values()))).any():
                    raise ValueError(f'{pid_str}Specified var is not defined at all graph nodes')

                var_update_dict = {k: np.sqrt(v/cov0) for k, v in var_update_dict.items()}

            else:
                raise ValueError(f'{pid_str}Specified var is not valid')

    elif method == 'ordinary_kriging':
        ordinary_kriging = True
        if verbose > 0 and mean is not None:
            print(f"{pid_str}WARNING: `mean` is ignored with `method='ordinary_kriging'`")

        mean = None

        if verbose > 0 and var is not None:
            print(f"{pid_str}WARNING: `var` is ignored with `method='ordinary_kriging'`")

        var = None

    else:
        raise ValueError(f'{pid_str}`method` ({method}) unknown')

    # List of data node labels / index
    data_node_labels = [u for u, v in prop_dict.items() if not np.isnan(v)]

    # Number of data nodes 
    n = len(data_node_labels)

    # Number of nodes in the graph
    n_nodes = G.number_of_nodes()
    if i0 is None:
        i0 = 0
    if i1 is None:
        i1 = n_nodes

    node_list = list(G.nodes())[i0:i1]
    node_list_len = i1 - i0

    if n == 0:
        if ordinary_kriging:
            krig_est_dict = {u: 0.0 for u in node_list}

            tmp = float(np.sqrt(cov0))
            krig_std_dict = {u: tmp for u in node_list}

        else: # simple kriging
            krig_est_dict = {u: float(mean_dict[u]) for u in node_list}

            tmp = float(np.sqrt(cov0))
            if var is not None:
                krig_std_dict = {u: tmp*var_update_dict[u] for u in node_list}
            else:
                krig_std_dict = {u: tmp for u in node_list}

        if update_graph:
            # Set node attributes (output)
            nx.set_node_attributes(G, krig_est_dict, kriging_est_node_attr)
            nx.set_node_attributes(G, krig_std_dict, kriging_std_node_attr)

        return krig_est_dict, krig_std_dict
    
    # Get the dictionary of error variance
    if err_std is None:
        err_std = 0.0

    if isinstance(err_std, float) or isinstance(err_std, int):
        err_var = err_std**2
        err_var_dict = {u:err_var for u in data_node_labels}

    elif isinstance(err_std, str):
        err_var_dict = nx.get_node_attributes(G, err_std)
        if np.any(np.asarray([u not in err_var_dict.keys() for u in data_node_labels])):
            raise ValueError(f'{pid_str}Specified error std is must be defined at all data nodes')
        
        err_var_dict = {u:s**2 for u, s in err_var_dict.items()}

    # Do kriging of nodes in node_list
    if use_unique_neighborhood:
        # Initialize 
        # - kriging matrix (mat) of order nmat
        # - right hand side of all kriging systems (b), matrix of dimension nmat x node_list_len
        if ordinary_kriging:
            nmat = n+1
            mat = np.ones((nmat, nmat))
            mat[-1,-1] = 0.0
        else:
            nmat = n
            mat = np.ones((nmat, nmat))

        b = np.ones((nmat, node_list_len))

        # Set kriging matrix (mat) and right hand side of all kriging systems (b)
        for i, ui in enumerate(data_node_labels):
            length_dict = nx.single_source_dijkstra_path_length(G, ui, cutoff=None, weight=edge_length_attr)
            h = np.asarray([length_dict[u] for u in data_node_labels])
            mat[i, :n] = cov_func(h)
            h = np.asarray([length_dict[u] for u in node_list])
            b[i, :] = cov_func(h)
        
        # Add error variance on the diagonal of the kriging matrix
        for i, u in enumerate(data_node_labels):
            mat[i, i] = mat[i, i] + err_var_dict[u]

        # Solve all kriging systems
        w = np.linalg.solve(mat, b) # w: matrix of dimension nmat x n_nodes

        # Kriged values
        if mean is not None:
            # simple kriging
            if var is not None:
                tmp = np.asarray([1.0/var_update_dict[ui] * (prop_dict[ui] - mean_dict[ui]) for ui in data_node_labels]).dot(w)
                krig_est_dict = {u: float(mean_dict[u] + var_update_dict[u]*v) for u, v in zip(node_list, tmp)}
            else:
                tmp = np.asarray([prop_dict[ui] - mean_dict[ui] for ui in data_node_labels]).dot(w)
                krig_est_dict = {u: float(mean_dict[u] + v) for u, v in zip(node_list, tmp)}
        else:
            # ordinary kriging
            tmp = np.asarray([prop_dict[ui] for ui in data_node_labels]).dot(w[:n, :])
            krig_est_dict = {u: float(v) for u, v in zip(node_list, tmp)}

        # Kriged standard deviation
        tmp = np.sqrt(np.maximum(0.0, cov0 - np.array([np.dot(w[:,i], b[:,i]) for i in range(node_list_len)])))
        krig_std_dict = {u: float(v) for u, v in zip(node_list, tmp)}

    else:
        # Limited search neighborhood

        # Set dmax (search radius)
        if searchRadius is not None:
            if searchRadius <= 0.0:
                raise ValueError(f'{pid_str}`searchRadius` not valid (negative)')

            dmax = searchRadius

        else:
            # use searchRadiusRelative
            if searchRadiusRelative <= 0.0:
                raise ValueError(f'{pid_str}`searchRadiusRelative` (factor) not valid (negative)')
            
            dmax = searchRadiusRelative * cov_model.r()

        if nneighborMax is None or nneighborMax > n or nneighborMax < 0:
            nneighborMax = n

        # Initialize kriging matrix, second member and property values at neigbhors
        mat = np.ones((nneighborMax+1, nneighborMax+1))
        b = np.ones(nneighborMax+1) 
        prop_val = np.ones(nneighborMax)

        # Initialize the dictionaries of kriging estimates and standard deviation
        krig_est_dict = {u:np.nan for u in node_list}
        krig_std_dict = {u:np.nan for u in node_list}

        if verbose > 0:
            progress_old = 0

        for k, u in enumerate(node_list):
            if verbose > 0:
                progress = int(k/node_list_len*100.0)
                if progress > progress_old:
                    print(f'{pid_str}Kriging: {progress:3d}%')
                    progress_old = progress
            
            if u in data_node_labels and err_var_dict[u] == 0.0:
                krig_est_dict[u] = prop_dict[u]
                krig_std_dict[u] = 0.0
                continue
            
            length_dict = nx.single_source_dijkstra_path_length(G, u, cutoff=dmax, weight=edge_length_attr)
            length_keys_list = list(length_dict.keys())
            ind = np.argsort(np.asarray(list(length_dict.values())))
            neighbor_labels = []
            nn = 0

            if neighbors_on_same_branch_only:
                # Get branch id
                br_ids = G.nodes[u]['branch_ids_list']

                for i in ind:
                    ui = length_keys_list[i]
                    if not np.isnan(prop_dict[ui]) and np.any(np.asarray([ui in branches[id] for id in br_ids])):
                        prop_val[nn] = prop_dict[ui]
                        b[nn] = cov_model(length_dict[ui])[0]
                        neighbor_labels.append(ui)
                        nn = nn+1
                        if nn >= nneighborMax:
                            break

            else:
                for i in ind:
                    ui = length_keys_list[i]
                    if not np.isnan(prop_dict[ui]):
                        prop_val[nn] = prop_dict[ui]
                        b[nn] = cov_model(length_dict[ui])[0]
                        neighbor_labels.append(ui)
                        nn = nn+1
                        if nn >= nneighborMax:
                            break

            if nn == 0:
                if mean is not None:
                    # simple kriging
                    krig_est_dict[u] = float(mean_dict[u])
                else:
                    # ordinary kriging
                    krig_est_dict[u] = mean_default

                krig_std_dict[u] = float(np.sqrt(cov0))
                continue
                
            for i in range(nn-1):
                ui = neighbor_labels[i]
                for j in range(i+1, nn):
                    uj = neighbor_labels[j]
                    h = nx.dijkstra_path_length(G, source=ui, target=uj, weight=edge_length_attr)
                    # h = nx.shortest_path_length(G, source=ui, target=uj, weight=edge_length_attr, method='dijkstra')
                    cov_h = cov_func(h)[0]
                    mat[i, j] = cov_h
                    mat[j, i] = cov_h
                
                mat[i, i] = cov0 + err_var_dict[ui]
            
            mat[nn-1, nn-1] = cov0 + err_var_dict[neighbor_labels[nn-1]]

            if ordinary_kriging:
                nmat = nn+1
                mat[nn, :] = 1.0
                mat[:, nn] = 1.0
                mat[nn, nn] = 0.0
                b[nn] = 1.0
            
            else:
                nmat = nn

            # Solve kriging system
            w = np.linalg.solve(mat[:nmat, :nmat], b[:nmat]) # w: vector of dimension nmat

            if mean is not None:
                # simple kriging
                if var is not None:
                    krig_est_dict[u] = float(mean_dict[u] + var_update_dict[u] * (np.asarray([1.0/var_update_dict[ui] for ui in neighbor_labels])*(prop_val[:nn] - np.asarray([mean_dict[ui] for ui in neighbor_labels]))).dot(w))
                else:
                    krig_est_dict[u] = float(mean_dict[u] + (prop_val[:nn] - np.asarray([mean_dict[ui] for ui in neighbor_labels])).dot(w))
            else:
                # ordinary kriging
                krig_est_dict[u] = float(prop_val[:nn].dot(w[:nn]))

            # Kriged standard deviation
            krig_std_dict[u] = float(np.sqrt(np.maximum(0.0, cov0 - b[:nmat].dot(w))))
        
        if verbose > 0:
            progress = 100
            if progress > progress_old:
                print(f'{pid_str}Kriging: {progress:3d}%')
                progress_old = progress

    if var is not None:
        krig_std_dict = {u: var_update_dict[u]*v for u, v in krig_std_dict.items()}

    if update_graph:
        # Set node attributes (output)
        nx.set_node_attributes(G, krig_est_dict, kriging_est_node_attr)
        nx.set_node_attributes(G, krig_std_dict, kriging_std_node_attr)
    
    return krig_est_dict, krig_std_dict
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
def krige_node_attribute_mp(
        G,
        cov_model,
        node_attr,
        err_std=0.0,
        edge_length_attr=None,
        method='ordinary_kriging',
        mean=None,
        var=None,
        neighbors_on_same_branch_only=False,
        branches=None,
        use_unique_neighborhood=False,
        searchRadius=None,
        searchRadiusRelative=1.2,
        nneighborMax=12,
        update_graph=True,
        kriging_est_node_attr=None,
        kriging_std_node_attr=None,
        verbose=1,
        nproc=-1):
    """
    Computes the same as the function :func:`krige_node_attribute`, using multiprocessing.

    All the parameters except `nproc` are the same as those of the function
    :func:`krige_node_attribute`.

    This function launches parallel processes [parallel calls of the
    function :func:`krige_node_attribute`]; the set of nodes to be kriged is distributed in a 
    balanced way over the processes.

    The number of processes used (in parallel) is determined by the parameter `nproc` 
    (int, default: -1); a negative number (or zero), -n <= 0, can be specified 
    to use the total number of cpu(s) of the system except n; `nproc` is finally
    at maximum equal to `G.number_of_nodes()` but at least 1 by applying:
        
    - if `nproc >= 1`, then `nproc = max(min(nproc, G.number_of_nodes()), 1)` is used
    - if `nproc = -n <= 0`, then `nproc = max(min(nmax-n, G.number_of_nodes()), 1)` is used, \
    where nmax is the total number of cpu(s) of the system (retrieved by 
    `multiprocessing.cpu_count()`)

    Note: if `nproc=None`, `nproc=-1` is used.

    See function :func:`krige_node_attribute` for details.
    """
    if update_graph:
        # Set node attribute names for output
        if kriging_est_node_attr is not None:
            if not isinstance(kriging_est_node_attr, str):
                raise ValueError('`kriging_est_node_attr` must be a string (node attribute name)')
        else:
            kriging_est_node_attr = f'{node_attr}_krig_est'

        if kriging_std_node_attr is not None:
            if not isinstance(kriging_std_node_attr, str):
                raise ValueError('`kriging_std_node_attr` must be a string (node attribute name)')
        else:
            kriging_std_node_attr = f'{node_attr}_krig_std'

    # Number of graph nodes
    n_nodes = G.number_of_nodes()

    # Set number of process(es): nproc
    if nproc is None:
        nproc = -1
    
    if nproc <= 0:
        nproc = max(min(multiprocessing.cpu_count() + nproc, n_nodes), 1)
    else:
        nproc_tmp = nproc
        nproc = max(min(int(nproc), n_nodes), 1)
        if verbose > 0 and nproc != nproc_tmp:
            print(f'Number of processes has been changed (now: nproc={nproc})')

    # Set index for distributing realizations
    q, r = np.divmod(n_nodes, nproc)
    ids_proc = [i*q + min(i, r) for i in range(nproc+1)]

    if verbose > 0:
        print(f'Running `krige_node_attribute` on {nproc} processes...')

    # Set pool of nproc workers
    pool = multiprocessing.Pool(nproc)
    out_pool = []
    for i in range(nproc):
        # Set i-th process
        args = (G, cov_model, node_attr)
        kwargs = dict(
                    err_std=err_std,
                    edge_length_attr=edge_length_attr, 
                    method=method,
                    mean=mean,
                    var=var,
                    neighbors_on_same_branch_only=neighbors_on_same_branch_only,
                    branches=branches,
                    use_unique_neighborhood=use_unique_neighborhood,
                    searchRadius=searchRadius, 
                    searchRadiusRelative=searchRadiusRelative, 
                    nneighborMax=nneighborMax,
                    update_graph=False, # do not update graph if multiprocessing is enabled
                    kriging_est_node_attr=kriging_est_node_attr,
                    kriging_std_node_attr=kriging_std_node_attr,
                    i0=ids_proc[i],
                    i1=ids_proc[i+1], 
                    pid=i,
                    verbose=verbose*(i==0))
        out_pool.append(pool.apply_async(krige_node_attribute, args=args, kwds=kwargs))

    # Properly end working process
    pool.close() # Prevents any more tasks from being submitted to the pool,
    pool.join()  # then, wait for the worker processes to exit.

    # Get result from each process
    out = [w.get() for w in out_pool]
    if np.any([x is None for x in out]):
        err_msg = f'`krige_node_attribute_mp`: an error occured on a process (worker)'
        raise ValueError(err_msg)

    # Gather dictionaries of all processes
    krig_est_dict = {}
    krig_std_dict = {}
    for krig_est_dict_pid_i, krig_std_dict_pid_i in out:
        krig_est_dict.update(krig_est_dict_pid_i)
        krig_std_dict.update(krig_std_dict_pid_i)

    if update_graph:
        # Set node attributes (output)
        nx.set_node_attributes(G, krig_est_dict, kriging_est_node_attr)
        nx.set_node_attributes(G, krig_std_dict, kriging_std_node_attr)
    
    return krig_est_dict, krig_std_dict
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
def sgs_node_attribute(
        G,
        cov_model,
        node_attr,
        err_std=0.0,
        edge_length_attr=None,
        method='ordinary_kriging',
        mean=None,
        var=None,
        neighbors_on_same_branch_only=False,
        branches=None,
        searchRadius=None,
        searchRadiusRelative=1.2,
        nneighborMax=12,
        nreal=1,
        seed=None,
        update_graph=True,
        sgs_node_attr=None,
        pid=None,
        verbose=1):
    """
    Performs Sequential Gaussian Simulation (SGS) of a node attribute on a graph.

    Parameters
    ----------
    G : networkx.Graph
        input graph

    cov_model : :class:`geone.covModel.CovModel1D`
        covariance model in 1D

    node_attr : str
        name of the node attribute to be simulated

    err_std : float, or str, default : 0.0
        standard deviation of error (zero-mean Gaussian of given std):
        - if float: the same error std is used for data at any graph node; 
        - if str: name of the node attribute of the error std

    edge_length_attr : str, optional
        name of the edge attribute for length (used to compute the distance between
        nodes);
        by default (`None`): the edges have a length of one

    method : str {'simple_kriging', 'ordinary_kriging'}, default: 'ordinary_kriging'
        type of kriging;
        note: if `method='ordinary_kriging'`, the parameter `mean` is not used

    mean : float, or str, default : 0.0
        kriging mean value :
        - if float: the value is used at all graph nodes; 
        - if str: name of the node attribute for the mean 
        
        note: if `method=ordinary_kriging`, parameter `mean` is ignored

    var : float, or str, optional
        kriging variance value :
        - if float: the value is used at all graph nodes; 
        - if str: name of the node attribute for the mean 

        note: if `method=ordinary_kriging`, parameter `var` is ignored
        
    neighbors_on_same_branch_only : bool, default: False
        - if `True`: only neighbors on the branch of the simulated graph node are used
        - if `False`: neighbors everywhere in the graph are used

        note: `neighbors_on_same_branch_only=True` requires `branches` (see below)

    branches : list, optional
        list of branches, required if `neighbors_on_same_branch_only=True`;
        each element (branch) is a list of node labels, typically `branches=kg.branches`
        should be used, where `kg` is the parent `KGraph` of `G`

    searchRadius : float, optional
        if specified, i.e. not `None`: search radius, i.e. 
        the data at graph node at distance to the estimated graph node greater 
        than `searchRadius` are not taken into account in the kriging system; 
        if `searchRadius` is not `None`, then `searchRadiusRelative` is not used;
        by default (`searchRadius=None`): `searchRadiusRelative` is used;

    searchRadiusRelative : float, default: 1.2
        used only if `searchRadius` is `None`;
        the search radius is set to `searchRadiusRelative` times the range of the 
        covariance model `cov_model`

    nneighborMax : int, default: 12
        maximal number of neighbors (data at graph nodes) taken into account in the
        kriging system; the data at graph nodes the closest to the estimated graph node are
        taken into account;
        note: if `nneighborMax=None` or `nneighborMax<0`, then `nneighborMax` is
        set to the number of graph nodes with informed data

    nreal : int, default: 1
        number of realization(s)

    seed : int, optional
        seed for initializing random number generator

    update_graph : bool, default True
        - if `True`: the simulations (realizations) are set as a node attribute (named \
        `sgs_node_attr` (see below))
        - if `False`: the graph is not updated (`sgs_node_attr` is not used)

    sgs_node_attr : str, optional
        name of the node attribute in output for simulation;
        by default (`None`),  the name `node_attr` + '_sgs' is used;
        the values of the this output node attributes are lists of length `nreal`
        containing all the realizations

    pid : int, optional
        process id of the caller (used with multiprocessing)
    
    verbose : int, default: 0
        verbose mode, higher implies more printing (info)

    Returns
    -------
    sgs : dict
        dictionary with node labels as keys and simulations (realizations) as values;
        each value is a list of length `nreal` containing all the realizations
    """
    if verbose > 0:
        if pid is not None:
            pid_str = f'[pid={pid}] '
        else:
            pid_str = ''

    # Check cov_model
    if not isinstance(cov_model, geone.covModel.CovModel1D):
        raise ValueError(f'{pid_str}`cov_model` must be an instance of `geone.covModel.CovModel1D`')

    # Prevent calculation if covariance model is not stationary
    if not cov_model.is_stationary():
        raise ValueError(f'{pid_str}`cov_model` is not stationary')

    # Check node attribute
    if node_attr not in kn.utils.get_node_attribute_names(G):
        raise ValueError(f'{pid_str}{node_attr} is not a node attribute')

    # Check edge length attribute
    if edge_length_attr is not None and edge_length_attr not in kn.utils.get_edge_attribute_names(G):
        raise ValueError(f'{pid_str}{edge_length_attr} is not an edge attribute')

    # Check neighbors_on_same_branch_only
    if neighbors_on_same_branch_only:
        if branches is None:
            raise ValueError(f'{pid_str}`branches` is required with `neighbors_on_same_branch_only=True`')

    # # Set dictionary to convert node label (id) to node index, and vice versa
    # node_label2index = {u:i for i, u in enumerate(G.nodes())}
    # #node_index2label = {i:u for i, u in enumerate(G.nodes())}

    if update_graph:
        # Set node attribute name for output
        if sgs_node_attr is not None:
            if not isinstance(sgs_node_attr, str):
                raise ValueError(f'{pid_str}`sgs_node_attr` must be a string (node attribute name)')
        else:
            sgs_node_attr = f'{node_attr}_sgs'

    # Get the dictionary of node attribute (property) to be kriged
    prop_dict = nx.get_node_attributes(G, node_attr)
    if len(prop_dict) == 0:
        raise ValueError(f'{pid_str}No value for specified node attribute ({node_attr})')

    v_first = np.asarray(list(prop_dict.values())[0]) # first value entry as array
    if v_first.ndim > 0:
        raise ValueError(f'{pid_str}Value of the specified attribute should by a scalar')

    # - set nan at uninformed nodes
    prop_dict = nx.get_node_attributes(G, node_attr, default=np.nan)

    # Covariance function and value at 0
    cov_func = cov_model.func() # covariance function
    cov0 = cov_func(0.)[0] # covariance function at origin (lag=0)

    # Default mean value    
    tmp = np.asarray(list(prop_dict.values()))
    if np.isnan(tmp).all():
        mean_default = 0.0
    else:
        # Set mean of data values
        mean_default = np.nanmean(tmp)
    
    # Method and mean, var
    if method == 'simple_kriging':
        ordinary_kriging = False
        # Get the dictionary of mean kriging values
        if mean is not None:
            if isinstance(mean, float) or isinstance(mean, int):
                mean_dict = {u:mean for u in G.nodes()}

            elif isinstance(mean, str):
                mean_dict = nx.get_node_attributes(G, mean, default=np.nan)
                if np.isnan(np.asarray(list(mean_dict.values()))).any():
                    raise ValueError(f'{pid_str}Specified mean is not defined at all graph nodes')

        else:
            mean_dict = {u:mean_default for u in G.nodes()}
    
        # Get the dictionary of "variance update"
        if var is not None:
            if isinstance(var, float) or isinstance(var, int):
                tmp = np.sqrt(var/cov0)
                var_update_dict = {u:tmp for u in G.nodes()}

            elif isinstance(var, str):
                var_update_dict = nx.get_node_attributes(G, var, default=np.nan)
                if np.isnan(np.asarray(list(var_update_dict.values()))).any():
                    raise ValueError(f'{pid_str}Specified var is not defined at all graph nodes')

                var_update_dict = {k: np.sqrt(v/cov0) for k, v in var_update_dict.items()}

            else:
                raise ValueError(f'{pid_str}Specified var is not valid')

    elif method == 'ordinary_kriging':
        ordinary_kriging = True
        if verbose > 0 and mean is not None:
            print(f"{pid_str}WARNING: `mean` is ignored with `method='ordinary_kriging'`")

        mean = None

        if verbose > 0 and var is not None:
            print(f"{pid_str}WARNING: `var` is ignored with `method='ordinary_kriging'`")

        var = None

    else:
        raise ValueError(f'{pid_str}`method` ({method}) unknown')

    # List of data node labels / index
    data_node_labels = [u for u, v in prop_dict.items() if not np.isnan(v)]

    # Number of data nodes 
    n = len(data_node_labels)

    # Get the dictionary of error variance
    # - initialize
    err_var_dict = {u:0.0 for u in G.nodes()}
    
    if err_std is None:
        err_std = 0.0

    if isinstance(err_std, float) or isinstance(err_std, int):
        err_var = err_std**2
        for u in data_node_labels:
            err_var_dict[u] = err_var 
            
    elif isinstance(err_std, str):
        tmp = nx.get_node_attributes(G, err_std)
        if np.any(np.asarray([u not in tmp.keys() for u in data_node_labels])):
            raise ValueError(f'{pid_str}Specified error std is must be defined at all data nodes')
        
        for u in data_node_labels:
            err_var_dict[u] = tmp[u]**2 

    # Number of nodes in the graph
    n_nodes = G.number_of_nodes()

    # Limited search neighborhood

    # Set dmax (search radius)
    if searchRadius is not None:
        if searchRadius <= 0.0:
            raise ValueError(f'{pid_str}`searchRadius` not valid (negative)')

        dmax = searchRadius

    else:
        # use searchRadiusRelative
        if searchRadiusRelative <= 0.0:
            raise ValueError(f'{pid_str}`searchRadiusRelative` (factor) not valid (negative)')
        
        dmax = searchRadiusRelative * cov_model.r()

    if nneighborMax is None or nneighborMax < 0:
        raise ValueError(f'{pid_str}`nneigbhorMax` not valid')

    # Initialize kriging matrix, second member and property values at neigbhors
    mat = np.ones((nneighborMax+1, nneighborMax+1))
    b = np.ones(nneighborMax+1) 
    prop_val = np.ones(nneighborMax)

    # Initialize the dictionary of sgs
    sgs_all_dict = {u:[] for u in G.nodes()}

    if seed is None:
        seed = np.random.randint(1, 1000000)
    seed = int(seed)

    if verbose > 0:
        progress_old = 0

    for ireal in range(nreal):
        # Initialize random number generator
        np.random.seed(seed+ireal)
        # sel_dict = {u:False for u in G.nodes()}
        # for u in data_node_labels:
        #     sel_dict[u] = True
        err_var_curr_dict = err_var_dict.copy()
        sgs_dict = {u:np.nan for u in G.nodes()}
        for u in data_node_labels:
            sgs_dict[u] = prop_dict[u]

        for k, u in enumerate(G.nodes()):
            if verbose > 0:
                progress = int((k+ireal*n_nodes)/(nreal*n_nodes)*100.0)
                if progress > progress_old:
                    print(f'{pid_str}SGS: {progress:3d}% ({ireal:3d} realizations done of {nreal})')
                    progress_old = progress
            
            # if u in data_node_labels and err_var_curr_dict[u] == 0.0:
            if not np.isnan(sgs_dict[u]) and err_var_curr_dict[u] == 0.0:
                continue

            length_dict = nx.single_source_dijkstra_path_length(G, u, cutoff=dmax, weight=edge_length_attr)
            length_keys_list = list(length_dict.keys())
            ind = np.argsort(np.asarray(list(length_dict.values())))
            neighbor_labels = []
            nn = 0

            if neighbors_on_same_branch_only:
                # Get branch id
                br_ids = G.nodes[u]['branch_ids_list']

                for i in ind:
                    ui = length_keys_list[i]
                    if not np.isnan(prop_dict[ui]) and np.any(np.asarray([ui in branches[id] for id in br_ids])):
                        prop_val[nn] = prop_dict[ui]
                        b[nn] = cov_model(length_dict[ui])[0]
                        neighbor_labels.append(ui)
                        nn = nn+1
                        if nn >= nneighborMax:
                            break

            else:
                for i in ind:
                    ui = length_keys_list[i]
                    if not np.isnan(sgs_dict[ui]):
                        prop_val[nn] = sgs_dict[ui]
                        b[nn] = cov_model(length_dict[ui])[0]
                        neighbor_labels.append(ui)
                        nn = nn+1
                        if nn >= nneighborMax:
                            break

            if nn == 0:
                # Mean and std (by kriging)
                if mean is not None:
                    mu = mean_dict[u]
                else:
                    mu = mean_default

                std = float(np.sqrt(cov0))
                if var is not None:
                    std = var_update_dict[u] * std
            else:                
                for i in range(nn-1):
                    ui = neighbor_labels[i]
                    for j in range(i+1, nn):
                        uj = neighbor_labels[j]
                        h = nx.dijkstra_path_length(G, source=ui, target=uj, weight=edge_length_attr)
                        # h = nx.shortest_path_length(G, source=ui, target=uj, weight=edge_length_attr, method='dijkstra')
                        cov_h = cov_func(h)[0]
                        mat[i, j] = cov_h
                        mat[j, i] = cov_h
                    
                    mat[i, i] = cov0 + err_var_dict[ui]
                
                mat[nn-1, nn-1] = cov0 + err_var_dict[neighbor_labels[nn-1]]

                if ordinary_kriging:
                    nmat = nn+1
                    mat[nn, :] = 1.0
                    mat[:, nn] = 1.0
                    mat[nn, nn] = 0.0
                    b[nn] = 1.0
                
                else:
                    nmat = nn

                # Solve kriging system
                w = np.linalg.solve(mat[:nmat, :nmat], b[:nmat]) # w: vector of dimension nmat

                # Mean and std (by kriging)
                if mean is not None:
                    # simple kriging
                    std = float(np.sqrt(np.maximum(0.0, cov0 - b[:nmat].dot(w))))
                    if var is not None:
                        mu = float(mean_dict[u] + var_update_dict[u] * (np.asarray([1.0/var_update_dict[ui] for ui in neighbor_labels])*(prop_val[:nn] - np.asarray([mean_dict[ui] for ui in neighbor_labels]))).dot(w))
                        std = var_update_dict[u] * std
                    else:
                        mu = float(mean_dict[u] + (prop_val[:nn] - np.asarray([mean_dict[ui] for ui in neighbor_labels])).dot(w))
                else:
                    # ordinary kriging
                    mu = float(prop_val[:nn].dot(w[:nn]))
                    std = float(np.sqrt(np.maximum(0.0, cov0 - b[:nmat].dot(w))))
        
            # Draw value in N(mu, std^2)
            sgs_dict[u] = np.random.normal(loc=mu, scale=std)
            err_var_curr_dict[u] = 0.0 # now the location is simulated, no more error taken into account
            
        # Store k-th realization
        for u in G.nodes():
            sgs_all_dict[u].append(sgs_dict[u])

    if verbose > 0:
        progress = 100
        if progress > progress_old:
            print(f'{pid_str}SGS: {progress:3d}% ({nreal:3d} realizations done of {nreal})')
            progress_old = progress

    if update_graph:
        # Set node attributes (output)
        nx.set_node_attributes(G, sgs_all_dict, sgs_node_attr)
    
    return sgs_all_dict
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
def sgs_node_attribute_mp(
        G,
        cov_model,
        node_attr,
        err_std=0.0,
        edge_length_attr=None,
        method='ordinary_kriging',
        mean=None,
        var=None,
        neighbors_on_same_branch_only=False,
        branches=None,
        searchRadius=None,
        searchRadiusRelative=1.2,
        nneighborMax=12,
        nreal=1,
        seed=None,
        update_graph=True,
        sgs_node_attr=None,
        nproc=-1,
        verbose=1):
    """
    Computes the same as the function :func:`sgs_node_attribute`, using multiprocessing.

    All the parameters except `nproc` are the same as those of the function
    :func:`sgs_node_attribute`.

    This function launches parallel processes [parallel calls of the
    function :func:`sgs_node_attribute`]; the set of realizations (specified by `nreal`) is
    distributed in a balanced way over the processes.

    The number of processes used (in parallel) is determined by the parameter `nproc` 
    (int, default: -1); a negative number (or zero), -n <= 0, can be specified 
    to use the total number of cpu(s) of the system except n; `nproc` is finally
    at maximum equal to `nreal` but at least 1 by applying:
        
    - if `nproc >= 1`, then `nproc = max(min(nproc, nreal), 1)` is used
    - if `nproc = -n <= 0`, then `nproc = max(min(nmax-n, nreal), 1)` is used, \
    where nmax is the total number of cpu(s) of the system (retrieved by 
    `multiprocessing.cpu_count()`)

    Note: if `nproc=None`, `nproc=-1` is used.

    Note: specifying a `seed` guarantees reproducible results whatever the number
    of processes used.

    See function :func:`sgs_node_attribute` for details.
    """
    if update_graph:
        # Set node attribute name for output
        if sgs_node_attr is not None:
            if not isinstance(sgs_node_attr, str):
                raise ValueError(f'`sgs_node_attr` must be a string (node attribute name)')
        else:
            sgs_node_attr = f'{node_attr}_sgs'

    # Set number of process(es): nproc
    if nproc is None:
        nproc = -1
    
    if nproc <= 0:
        nproc = max(min(multiprocessing.cpu_count() + nproc, nreal), 1)
    else:
        nproc_tmp = nproc
        nproc = max(min(int(nproc), nreal), 1)
        if verbose > 0 and nproc != nproc_tmp:
            print(f'Number of processes has been changed (now: nproc={nproc})')
    
    # Set index for distributing realizations
    q, r = np.divmod(nreal, nproc)
    ids_proc = [i*q + min(i, r) for i in range(nproc+1)]

    if verbose > 0:
        print(f'Running `sgs_node_attr` on {nproc} processes...')

    # Set seed (base)
    if seed is None:
        seed = np.random.randint(1, 1000000)
    seed = int(seed)

    # Set pool of nproc workers
    pool = multiprocessing.Pool(nproc)
    out_pool = []
    for i in range(nproc):
        # Set i-th process
        args = (G, cov_model, node_attr)
        kwargs = dict(
                    err_std=err_std,
                    edge_length_attr=edge_length_attr,
                    method=method,
                    mean=mean,
                    var=var,
                    neighbors_on_same_branch_only=neighbors_on_same_branch_only,
                    branches=branches,
                    searchRadius=searchRadius,
                    searchRadiusRelative=searchRadiusRelative,
                    nneighborMax=nneighborMax,
                    nreal=ids_proc[i+1]-ids_proc[i], 
                    seed=seed+ids_proc[i],
                    update_graph=False, # do not update graph if multiprocessing is enabled
                    sgs_node_attr=sgs_node_attr,
                    pid=i,
                    verbose=verbose*(i==0))
        out_pool.append(pool.apply_async(sgs_node_attribute, args=args, kwds=kwargs))

    # Properly end working process
    pool.close() # Prevents any more tasks from being submitted to the pool,
    pool.join()  # then, wait for the worker processes to exit.

    # Get result from each process
    out = [w.get() for w in out_pool]
    if np.any([x is None for x in out]):
        err_msg = f'`sgs_node_attribute_mp`: an error occured on a process (worker)'
        raise ValueError(err_msg)

    # Gather dictionaries of all processes
    sgs_all_dict = {u:[] for u in G.nodes()}
    for sgs_all_dict_pid_i in out:
        for u in G.nodes():
            sgs_all_dict[u].extend(sgs_all_dict_pid_i[u])

    if update_graph:
        # Set node attributes (output)
        nx.set_node_attributes(G, sgs_all_dict, sgs_node_attr)
    
    return sgs_all_dict
# ----------------------------------------------------------------------------

# =============================================================================
# Function for K-means
# =============================================================================

# -----------------------------------------------------------------------------
def kmeans(
        points, 
        nclusters, 
        fixed_cluster_list=None, 
        coord_weight=None, 
        max_iter=100, 
        tol=1e-4, 
        seed=None, 
        verbose=0):
    """
    Performs K-means clustering on a set of points.

    Parameters
    ----------
    points : 2d numpy array of shape (n, d)
        points in dimension d:

        - points[i] : coordinates of i-th points 
    
    nclusters : int (> 0)
        number of clusters

    fixed_cluster_list: list of lists of ints, optional
        if given:
        
        - fixed_cluster_list[i]: list of the indices of points that \
        form one cluster

        by default (`None`): no fixed cluster set as constraint

    coord_weight : 1d numpy array of shape (d, ), optional
        weights for each coordinate dimension to compute the distance between points; 
        by default (`None`): all dimensions are equally weighted;

    max_iter : int, default: 100
        maximum number of iterations to perform
    
    tol : float (> 0), default: 1e-4
        tolerance for convergence, i.e. the algorithm stops if the
        norm of the difference between new and old centers is less than `tol`
    
    seed : int, optional
        random seed for reproducibility; the first set of cluster centers is 
        selected randomly among the input points;
        by default (`None`): no seed is set
    
    verbose : int, default: 0
        verbosity level: higher value means more verbose output

    Returns
    -------
    centers : 2d numpy array of shape (nclusters, d)
        coordinates of the cluster centers, in dimension d:
        
        - centers[i] : coordinates of the i-th cluster center

    cluster_ind : 1d numpy array of shape (n,)
        index of the cluster to which each point belongs, with:

        - cluster_ind[i] : index of the cluster to which the i-th point belongs        
    """
    if not isinstance(points, np.ndarray) or points.ndim != 2:
        raise ValueError("`kmeans`: points must be a 2D array (numpy.ndarray)")

    if nclusters <= 0:
        raise ValueError("`kmeans`: number of clusters must be positive")

    if nclusters > points.shape[0]:
        raise ValueError("`kmeans`: number of clusters must be less than or equal to number of points")
    
    if fixed_cluster_list is not None:
        if not isinstance(fixed_cluster_list, list) and not isinstance(fixed_cluster_list, tuple):
            raise ValueError("`kmeans`: `fixed_cluster_list` should be a list (or tuple) of lists of ints")
        
        nfixed_clusters = len(fixed_cluster_list)
        if nfixed_clusters > nclusters:
            raise ValueError("`kmeans`: number of fixed clusters (length of `fixed_cluster_list`) is greater than `nclusters`")

        all_nodes = list(range(points.shape[0]))
        for node_list in fixed_cluster_list:
            if not isinstance(node_list, list) and not isinstance(node_list, tuple) and not (isinstance(node_list, np.ndarray) and node_list.ndim == 1):
                raise ValueError("`kmeans`: any element of `fixed_cluster_list` should be a sequence of ints")
        
            if len(np.intersect1d(node_list, all_nodes)) != len(node_list):                
                raise ValueError("`kmeans`: an element of `fixed_cluster_list` is invalid")
            
        all_fixed_nodes = [i for node_list in fixed_cluster_list for i in node_list]
        if len(np.unique(all_fixed_nodes)) != len(all_fixed_nodes):
            raise ValueError("`kmeans`: fixed clusters (elements of `fixed_cluster_list`) are not disjoint")
    else:
        nfixed_clusters = 0

    if coord_weight is not None:
        if len(coord_weight) != points.shape[1]:
            raise ValueError("`kmeans`: `coord_weight` must have the same length as the number of dimensions of points")
        
        if np.any(coord_weight <= 0):
            raise ValueError("`kmeans`: `coord_weight` must be positive for all dimensions")
        
        points = points * coord_weight  # Apply coordinate weights to points

    # else:
    #     # If coord_weight is not provided, set it to 1 for each dimension
    #     coord_weight = np.ones(points.shape[1])

    if nfixed_clusters > 0:
        # Compute centers (and cluster indices) for fixed clusters
        fixed_centers = np.vstack([points[node_list].mean(axis=0) for node_list in fixed_cluster_list])
        fixed_cluster_ind = np.zeros(points.shape[0], dtype='int')

        for i, node_list in enumerate(fixed_cluster_list):
            fixed_cluster_ind[node_list] = i

        # Get list of "free" points and number of "free" clusters
        free_nodes = np.setdiff1d(all_nodes, all_fixed_nodes)
        points = np.vstack([points[i] for i in free_nodes])

        nclusters = nclusters - nfixed_clusters
        if nclusters > points.shape[0]:
            raise ValueError("`kmeans`: number of free clusters must be less than or equal to number of free points")
        
        if nclusters == 0:
            if points.shape[0] > 0:
                raise ValueError("`kmeans`: no free cluster, but it remains free points")

            return fixed_centers, fixed_cluster_ind

    if seed is not None:
        np.random.seed(seed)

    # Initialize centers randomly from the points
    ind = np.random.choice(points.shape[0], nclusters, replace=False)
    centers = points[ind]

    # Initialize convergence flag
    converge = False

    for iter in range(max_iter):
        # Assign each point to the nearest center
        cluster_ind = np.asarray([np.argmin(np.sum((point - centers)**2, axis=1)) for point in points])
        # cluster_ind = np.asarray([np.argmin(np.linalg.norm(point - centers, axis=1)) for point in points])
        # cluster_ind = np.asarray([np.argmin(np.linalg.norm(coord_weight * (point - centers), axis=1)) for point in points])

        # Recalculate centers as the mean of each cluster
        new_centers = np.asarray([points[cluster_ind == i].mean(axis=0) for i in range(nclusters)])

        # Check for convergence
        if np.all(np.linalg.norm(new_centers - centers, axis=1) < tol):
            converge = True

        centers = new_centers

        if converge:
            break

    if not converge:
        raise ValueError(f'K-means did not converge after {max_iter} iterations (tol={tol})')

    if verbose > 0:
        print(f'K-means converged after {iter+1} iterations (tol={tol})')

    if nfixed_clusters > 0:
        # Integrate fixed centers and update cluster indices
        # - centers
        centers = np.vstack((fixed_centers, centers))
        # - update cluster indices
        cluster_ind = cluster_ind + nfixed_clusters        
        for i, j in enumerate(free_nodes):
            fixed_cluster_ind[j] = cluster_ind[i]
        cluster_ind = fixed_cluster_ind

    if coord_weight is not None:
        # If `coord_weight` was applied, apply its inverse for the centers
        centers = centers / coord_weight

    return centers, cluster_ind
# -----------------------------------------------------------------------------

# =============================================================================
# Spectral embedding - reducing graph based on K-means
# =============================================================================

# ----------------------------------------------------------------------------
def reduce_graph_spectral_kmeans(
        G, 
        nclusters,
        embedding_space_dim,
        edge_weight=None,
        node_attr_mode='mean', 
        node_attr_list=None,
        fixed_cluster_dict=None,
        eigval=None, 
        eigvec=None, 
        return_new_id_dict=False,
        return_orig_id_dict=False,
        return_medoid_id_dict=False,
        return_cluster_diameter_dict=False,
        seed=None,
        kmeans_kwargs=None,
        verbose=0):
    """
    Reducing (coarsening) graph via spectral embedding and K-means.

    The input graph `G` is assumed to be connected (one connected component); 
    it is reduced (coarsened) by applying the following procedure.
    
    1. Spectral embedding.
    The `n` graph nodes are embedded in a space of dimension `m=embedding_space_dim`
    as follows. Let :math:`u_j` and :math:`\\lambda_j`, :math:`j=0, \\ldots, n-1`, be
    the eigenvectors and eigenvalues of the Laplacian matrix of `G`, with 
    :math:`0 = \\lambda_0 < \\lambda_1 \leq \\ldots \\lambda_{n-1}`. The graph node of index
    `i` (:math:`i = 0, \\ldots, n-1`) is represented by the point of coordinate 
    :math:`(u_1(i), \\ldots, u_m(i))`; i.e. with `U` the :math:`n\times m` matrix
    whose columns are the first eigenvectors with strictly positive eigenvalues
    of the Laplacian matrix, the rows of `U` are the points representing the
    graph nodes

    2. K-means (with coordinate weight) in the embedding space.
    Then the K-means algorithm is applied on the points representing the graph nodes
    in the embedding space, with the coordinate `j` re-scaled (weighted) by 
    :math:`\\lambda_j^{-1/2}`, i.e. the points
    :math:`(\\lambda_1^{-1/2}\\cdot u_1(i), \\ldots, \\lambda_m^{-1/2}\\cdot u_m(i))`
    are considered

    3. Building the reduced (coarse) graph.
    The clusters resulting from K-means are retrieved and form the new 
    nodes (super-nodes) in the output graph. The node attributes (properties) in 
    the output graph (if any) are defined according to the keyword argument 
    `node_attr_mode` from the attributes of the original graph nodes in the cluster 
    
    - `node_attr_mode='mean'`: mean is used
    - `node_attr_mode='medoid'`: medoid is used, in that case each cluster must be \
    connected, i.e. must have one connected component; the medoid of the subgraph of the \
    input graph induced by the nodes in the cluster is computed, where the medoid is the \
    node realizing the minimum of the sum of the shortest path lengths wrt. inverse of edge \
    weight (i.e., wrt. resistances (instead of conductances determining the Laplacian matrix)) \
    to any other node
    
    Finally, an edge is set between two super-nodes if the two corresponding clusters are 
    linked in the original graph.

    An explicit list of node attributes can be specified by the parameters
    `node_attr_list`, and the corresponding operation for the aggregation 
    (in cluster) can be specified in a list by the parameters `node_attr_mode`. 
    By default, all node attributes are considered, and the same operation is 
    applied.

    Moreover, a pre-defined clusters may be set with the parameter 
    `fixed_cluster_dict`: dictionary with input graph node labels as keys and
    integer as values, where the nodes with the same value > 0 form a cluster
    (super-node in the output graph ; the nodes with value <= 0 (or not in the
    dictionary) are treated with the K-means algorithm (no constraint on those
    nodes).
    
    Parameters
    ----------
    G : networkx.Graph
        input graph
    
    ncluster : int (> 0)
        number of clusters, i.e. number of (super-)nodes in output graph,
        including pre-defined fixed cluster (if any, given via `fixed_cluster_dict`)
    
    embedding_space_dim : int (> 0)
        space dimension for spectral embedding (dimension of the space on which
        K-means is applied)

    edge_weight : str, optional
        name of the edge attribute used as weight (for the Laplacian matrix);
        by default (`None`) : all edges have a weight of 1

    node_attr_mode : str {'mean', 'medoid'} (or list of strs), default: 'mean'
        string or list of strings, the strings indicate how the node attributes 
        are computed in the reduced graph (see above)

        - if string: the same operation (mode) is used for all considered \
        node attributes
        - if list of strings: the length of the list must be equal to the length \
        of the list `node_attr_list`, and `node_attr_mode[i]` indicates the \
        operation (mode) used for the node attribute `node_attr_list[i]`

    node_attr_list : list of strs, optional
        list of node attributes of the input graph to be included in 
        the reduced graph;
        by default (`None`): the list of all nodes attributes is considered
    
    fixed_cluster_dict : dict, optional
        dictionary defining pre-defined clusters, the keys are input graph node 
        labels, and the values integer ids; the id 0 is the default id attached to 
        any node that is not listed in the keys; the nodes with the same id > 0 form
        a pre-defined cluster;
        by default (`None`): not applied, i.e. no extra constraint on clusters
      
    eigval : 1d-array, optional
        array of shape (m, ), with `embedding_space_dim+1` <= m <= n, where
        n is the number of input graph nodes, the m first eigen values of the 
        Laplacian matrix, in ascending order, i.e. for a connected graph, 
        `0 = eigval[0] < eigenval[1] <= eigenval[2] ...`;
        by default (`None`): the first m = `embedding_space_dim+1` eigen values 
        are computed

    eigvec : 2d-array, optional
        array of shape (n, m), where `embedding_space_dim`+1 <= m <= n, with n
        the number of input graph nodes, the m first eigen vectors in columns 
        (orthogonal and of norm 1) corresponding to the eigen values `eigval`; 
        for a connected graph, `eigvec[:,0]` is proportional to the 
        vector `[1, 1, ..., 1]`;
        by default (`None`): the first m = `embedding_space_dim+1` eigen vectors 
        are computed

    return_new_id_dict : bool, default: False
        if `True`, the dictionary `new_id_dict` is returned, where
        the keys are the node ids (labels) in the input graph, and the values
        are the corresponding node ids in the output graph (cluster ids)
        
    return_orig_id_dict : bool, default: False
        if `True`, the dictionary `orig_id_dict` is returned, where
        the keys are the node ids (labels) in the output graph (cluster), and the 
        values are the lists of the corresponding node ids (labels) in the input graph
        (original ids)
        
    return_medoid_id_dict : bool, default: False
        if `True`, the dictionary `medoid_id_dict` is returned, where
        the keys are the node ids (labels) in the output graph (cluster), and the 
        values are the node ids (labels) in the input graph of the cluster medoids;
        the medoid of the subgraph of the input graph induced by the nodes in a
        considered cluster is computed, where the medoid is the node realizing the
        minimum of the sum of the shortest path lengths (wrt. inverse of edge weight,
        i.e. wrt. resistances (instead of conductances determining the Laplacian matrix))
        to any other node

    return_cluster_diameter_dict : bool, default: False
        if `True`, the dictionary `cluster_diameter_dict` is returned, where
        the keys are the node ids (labels) in the output graph (cluster), and the 
        values are the diameter of the cluster (box), i.e. the maximal
        distance (in the original graph, wrt. inverse of edge weight) between two nodes 
        in the cluster
        
    seed : int, optional
        seed number for initializing the random number generator

    kmeans_kwargs : dict, optional
        keyword arguments passed to the function `kmeans`, the keys
        `fixed_cluster_list`, `coord_weight`, `verbose` must not be given 
        (set in the current function)

    verbose : int, default: 0
        verbose mode, larger value implies more info printed

    Returns
    -------
    G_red : networkx.Graph
        reduced graph (see above)
    
    new_id_dict : dict, optional
        returned if `return_new_id_dict=True`: dictionary giving the new id
        (label of node in the output graph) (value) for each original id (label of 
        node in the input graph) (key) (see `return_new_id_dict` above);
        note: this property can be attached to nodes of the input graph `G` with:
        `nx.set_node_attributes(G, new_id_dict, 'new_id')`

    orig_id_dict : dict, optional
        returned if `return_orig_id_dict=True`: dictionary giving the list of 
        original ids (labels of nodes in the input graph) (value) for each new id 
        (label of node in the output graph) (key) (see `return_orig_id_dict` above);
        note: this property can be attached to nodes of the output graph `G_red` with:
        `nx.set_node_attributes(G_red, orig_id_dict, 'orig_id')`

    medoid_id_dict : dict, optional
        returned if `return_medoid_id_dict=True`: dictionary giving the medoid
        id of the cluster (box) (label of node in the input graph) (value), for each 
        new id (label of node in the output graph) (key) (see `return_medoid_id_dict` 
        above);
        note: this property can be attached to nodes of the output graph `G_red` with:
        `nx.set_node_attributes(G_red, medoid_id_dict, 'medoid_id')`

    cluster_diameter_dict : dict, optional
        returned if `return_cluster_diameter_dict=True`: dictionary giving the diameter 
        of the cluster (box) (value), for each new id (label of node in the output graph) 
        (key) (see `return_cluster_diameter_dict` above);
        note: this property can be attached to nodes of the output graph `G_red` with:
        `nx.set_node_attributes(G_red, cluster_diameter_dict, 'cluster_diameter')`
    """
    if verbose > 0:
        print(f'Reduce (coarsen) graph via spectral embedding and K-means...')

    # Check if input graph is connected
    if nx.number_connected_components(G) > 1:
        raise ValueError('Input graph is not connected (more than one connected component)')
    
    # Check embedding space dimension
    if not isinstance(embedding_space_dim, int) or embedding_space_dim < 1 or embedding_space_dim > G.number_of_nodes()-1:
        raise ValueError('`embedding_space_dim` is not valid (must be an integer in [1, number of graph nodes - 1])')

    # Check node attributes of the input graph to be considered in the reduced graph
    # and corresponding operation (mode)

    # Get keys (attribute names) in the input graph
    # # - all keys
    # keys = [list(G_copy.nodes[i].keys()) for i in G_copy.nodes()] # list of lists
    # node_attr_all = np.unique([xij for xi in keys for xij in xi])
    # # - from one node
    # node_attr_all = G_copy.nodes[list(G_copy.nodes())[0]].keys()
    node_attr_all = kn.utils.get_node_attribute_names(G)

    # Check given node attributes name
    if node_attr_list is not None:
        if not np.all([k in node_attr_all for k in node_attr_list]):
            raise ValueError('Attribute name does not exist, check `attr_node_list` parameters')
    else:
        node_attr_list = node_attr_all

    # Check given mode
    if isinstance(node_attr_mode, list):
        if len(node_attr_mode) != len(node_attr_list):
            raise ValueError('Length of the list `node_attr_mode` is not valid')

        if np.any([s not in ('mean', 'medoid') for s in node_attr_mode]):
            raise ValueError('Entry of the list `node_attr_mode` is not valid')

    else: # `node_attr_mode` assumed to be a string
        if node_attr_mode not in ('mean', 'medoid'):
            raise ValueError('`node_attr_mode` is not valid')
        node_attr_mode = len(node_attr_list)*[node_attr_mode]

    # Set correspondance between node label and node index
    node_label2index = {u:i for i, u in enumerate(G.nodes())}
    node_index2label = {i:u for i, u in enumerate(G.nodes())}

    if seed is not None:
        np.random.seed(seed)

    # Set list of fixed cluster for K-means
    fixed_cluster_list = None # Initialization
    if fixed_cluster_dict is not None:
        ids = np.unique(list(fixed_cluster_dict.values())) # all ids
        ids = np.asarray([id for id in ids if id != 0])    # remove id 0
        nfixed_clusters = len(ids)
        if nfixed_clusters:
            fixed_cluster_list = [[] for _ in range(nfixed_clusters)]
            for u, id in fixed_cluster_dict.items():
                if id != 0:
                    j = np.where(ids == id)[0][0]
                    fixed_cluster_list[j].append(node_label2index[u])

    # Compute (if needed) the eigenvalue and eigenvectors of the Laplacian matrix
    if verbose > 0:
        print(f'Get or compute eigenvalues and eigenvectors of the Laplacian matrix...')
    
    if eigvec is None or eigval is None:
        eigval, eigvec = laplacian_eigs(G, edge_weight=edge_weight, k=embedding_space_dim+1)
    elif eigvec.shape[1] < embedding_space_dim+1 or len(eigval) < embedding_space_dim+1:
        raise ValueError('dimension of `eigvec` and/or `eigval` is not valid')

    # Do K-means:
    #   Let u[j] and lambda[j] the eigenvectors and eigenvalues of the Laplacian matrix,
    #   with lambda[0] = 0 < lambda[1] <= ... 
    #   K-means is applied on the points representing the graph nodes:
    #   - point i : node i of the graph with coordinates [u[1][i], ..., u[embedding_space_dim][i]]
    #   - coordinate j of all points are rescaled by lambda[j]^{-1/2}
    if verbose > 0:
        print(f'Apply K-means in the spectral embedding space (dim={embedding_space_dim}, nclusters={nclusters})...')

    if kmeans_kwargs is not None:
        if 'fixed_cluster_list' in kmeans_kwargs.keys():
            raise ValueError(f"'fixed_cluster_list' must not be given in `kmeans_kwargs`") 
        if 'coord_weight' in kmeans_kwargs.keys():
            raise ValueError(f"'coord_weight' must not be given in `kmeans_kwargs`") 
        if 'verbose' in kmeans_kwargs.keys():
            raise ValueError(f"'verbose' must not be given in `kmeans_kwargs`") 
    else:
        kmeans_kwargs = {}
    
    points = eigvec[:, 1:embedding_space_dim+1]
    coord_weight = 1./np.sqrt(eigval[1:embedding_space_dim+1])
    try:
        _, cluster_ind = kmeans(
                            points,
                            nclusters,
                            fixed_cluster_list=fixed_cluster_list,
                            coord_weight=coord_weight,
                            **kmeans_kwargs)
    except Exception as exc:
        raise ValueError(f'exc')



    # Set dictionary `nodes_orig_id_dict` according to cluster index
    nodes_new_id_dict = {node_index2label[i]:int(cluster_ind[i]) for i in range(G.number_of_nodes())}

    # Set the list of original ids for each new id (in a dictionary)
    # (returned optionally, and used for computing mean or medoid of node attributes over cluster (see below))
    nodes_orig_id_dict = {i:[] for i in range(nclusters)}
    for id in G.nodes():
        nodes_orig_id_dict[nodes_new_id_dict[id]].append(id)

    # # Check if every cluster is connected
    # ok = True
    # for i in range(nclusters):
    #     if nx.number_connected_components(nx.subgraph(G, nodes_orig_id_dict[i])) > 1:
    #         ok = False
    #         break
    
    # if not ok:
    #     raise ValueError('A cluster is not connected (after K-means)')

    # Check if every cluster is connected
    cluster_ncc = np.asarray([nx.number_connected_components(nx.subgraph(G, nodes_orig_id_dict[i])) for i in range(nclusters)])
    if np.any(cluster_ncc !=1 ):
        ind = np.where(cluster_ncc > 1)[0]
        for i in ind:
            print(f'... error: cluster {i} has {cluster_ncc[i]} connected components')
        
        raise ValueError('Cluster(s) not connected (after K-means)')

    # Build the reduced (coarse) graph (output)
    if verbose > 0:
        print(f'Build the output graph (reduced / coarse graph) ...')
    
    G_red = nx.Graph()
    # Set nodes (clusters in the reduced graph)
    G_red.add_nodes_from(range(nclusters))
    # Set edges (links between clusters in the reduced graph)
    for u, v in G.edges():
        if nodes_new_id_dict[u] != nodes_new_id_dict[v]:
            G_red.add_edge(nodes_new_id_dict[u], nodes_new_id_dict[v])

    # Compute medoid of each cluster (box) and set cluster (box) diameter, if needed
    compute_medoid_dict = return_medoid_id_dict or 'medoid' in node_attr_mode
    if compute_medoid_dict or return_cluster_diameter_dict:
        if compute_medoid_dict:
            medoid_id_dict = {} # initialization of the dictionary
        
        if return_cluster_diameter_dict:
            cluster_diameter_dict = {} # initialization of the dictionary

        for i in G_red.nodes():
            G_cluster = nx.subgraph(G, nodes_orig_id_dict[i])
            # if nx.number_connected_components(G_cluster) > 1: # already checked above
            #     raise ValueError('A cluster is not connected (after K-means)')

            if edge_weight is not None:
                w = nx.get_node_attributes(G_cluster, edge_weight) # get "conductance" in Laplacian matrix 
                w = {u:1/v for u, v in w.items()}                  # set "resistance" (1/conductance) to compute distance
                nx.set_node_attributes(G_cluster, w, 'w')
                dist_edge_weight = 'w'
            else:
                dist_edge_weight = None
        
            if compute_medoid_dict:
                d_sum_min = np.inf
                for (u, dist_from_u) in nx.all_pairs_dijkstra_path_length(G_cluster, weight=dist_edge_weight):
                    d_sum = np.asarray(list(dist_from_u.values())).sum()
                    if d_sum < d_sum_min:
                        medoid = u
                        d_sum_min = d_sum

                medoid_id_dict[i] = medoid

            if return_cluster_diameter_dict:
                cluster_diameter_dict[i] = nx.diameter(G_cluster, weight=dist_edge_weight)

    # Set node attributes in the reduced graph (according to `node_attr_list` and `node_attr_mode`)
    for attr, mode in zip(node_attr_list, node_attr_mode):
        d = nx.get_node_attributes(G, attr, default=np.nan) # dictionary original_id:value_of_attribute

        if mode == 'mean':
            v0 = list(d.values())[0] # value of one node in G
            if hasattr(v0, '__len__'):
                # attr_type = type(v0)[0]     # attribute type
                for i in G_red.nodes():
                    # G_red.nodes[i][attr] = attr_type(np.nanmean(np.asarray([np.atleast_1d(d[j]) for j in nodes_orig_id_dict[i]]), axis=0))
                    G_red.nodes[i][attr] = np.nanmean(np.asarray([np.atleast_1d(d[j]) for j in nodes_orig_id_dict[i]]), axis=0)
            else:
                for i in G_red.nodes():
                    G_red.nodes[i][attr] = np.nanmean(np.asarray([np.atleast_1d(d[j]) for j in nodes_orig_id_dict[i]]), axis=0)[0]
        
        elif mode == 'medoid':
            for i in G_red.nodes():
                G_red.nodes[i][attr] = d[medoid_id_dict[i]]

    out = [G_red]
    if return_new_id_dict:
        out.append(nodes_new_id_dict)
    if return_orig_id_dict:
        out.append(nodes_orig_id_dict)
    if return_medoid_id_dict:
        out.append(medoid_id_dict)
    if return_cluster_diameter_dict:
        out.append(cluster_diameter_dict)
   
    if len(out) == 1:
        out = out[0]
    else:
        out = tuple(out)

    return out
# ----------------------------------------------------------------------------

# =============================================================================
# Maximum Excluded Mass Burning (MEMB) algorithm
# =============================================================================

# ----------------------------------------------------------------------------
def reduce_graph_memb(
        G, 
        r=1, 
        edge_weight=None,
        node_attr_mode='mean', 
        node_attr_list=None,
        node_group_id_dict=None,
        return_new_id_dict=False,
        return_orig_id_dict=False,
        return_center_id_dict=False,
        return_medoid_id_dict=False,
        return_cluster_diameter_dict=False,
        seed=None,
        verbose=0):
    """
    Maximum Excluded Mass Burning (MEMB) - for weighted and unweigted graphs.

    The Maximum-Excluded Mass Burning (MEMB) algorithm, based 
    on a radius `r`, consists of a box covering of the input graph nodes
    where each cluster (box) has a diameter of at most `2*r` and each cluster
    is connected, i.e. for any pair of nodes in the cluster, there exists a path 
    entirely in the cluster between the two nodes.

    A reduced graph is provided by the box covering, where
    - the (new) nodes are the clusters (boxes) (also called super-nodes)
    - two (super-)nodes C1, C2, are linked by an edge if there exists two \
    nodes u1, u2 of the input graph in the clusters corresponding to C1, C2 \
    that are linked with an edge

    The idea of the MEMB algorithm is to used clusters centered on nodes 
    of the input graph: for a node c in the graph, the center of a cluster, 
    all nodes u at a distance from c less than or equal to `r` are can be 
    considered in the cluster of center c.
    
    The algorithm proceeds in 3 steps.

    1. Identifying the centers in the original graph.    
    
    (i) All nodes are marked as uncovered and non-centers.
    
    (ii) For all non-center nodes x (including already covered nodes),
    the excluded mass is computed, i.e. the number of uncovered nodes
    in the cluster centered at x, and among these nodes, select a 
    node c realizing the maximum of the excluded mass, and mark c 
    as a center.

    (iii) Mark all the nodes in the cluster centered at c (i.e. nodes 
    at a distance from c less than or equal to r) as covered

    (iv) Repeat steps (ii) and (iii) until all nodes are covered.

    2. Creating the clusters.

    (i) Assign an id to each center.

    (ii) For all nodes, compute the central distance, i.e. the distance
    to the nearest center.

    (iii) Set a list of all non-center nodes, sorted according to 
    increasing central distance.

    (iv) Take the first node u of the list, select one of its 
    neighbors, v, with a smaller central distance, and assign the id
    of v to u. Remove the node u from the list.

    (v) Repeat (iv) until the list is empty.

    3. Building the reduced graph.

    Finally, the nodes having the same id form one cluster (connected
    in the original graph by construction). The reduced graph is 
    then computed such that:
    
    - one node is one cluster, with node attributes (if any) defined \
    according to the keyword argument `node_attr_mode` from the attributes \
    of the original graph nodes in the cluster (`node_attr_mode='center'`: \
    atttribute from the central node is used, `node_attr_mode='mean'`: mean \
    is used, `node_attr_mode='medoid'`: medoid is used, the medoid of the subgraph \
    of the input graph induced by the nodes in the cluster is computed, where the \
    medoid is the node realizing the minimum of the sum of the shortest path lengths \
    (wrt. edge weight) to any other node
    - an edge is set between two nodes if the two corresponding clusters \
    are linked in the original graph

    An explicit list of node attributes can be specified by the parameters
    `node_attr_list`, and the corresponding operation for the aggregation 
    (in cluster) can be specified in a list by the parameters `node_attr_mode`. 
    By default, all node attributes are considered, and the same operation is 
    applied.

    Moreover, a "group id" can be assigned to the nodes of the original graph
    with the parameter `node_group_id_dict`, with the following meaning: two 
    nodes with a different group id must be in a different cluster (box). By
    default, the group id 0 is used (for node without group id). The algorithm
    starts by "deconnecting" the different groups, which ensures that two 
    nodes with different group ids will never be in a same cluster.
    
    Parameters
    ----------
    G : networkx.Graph
        input graph
    
    r : int or float, default: 1
        radius of a cluster (distance to a central node; every node connected 
        to a central node with at most `r` may belong to the same cluster)
    
    edge_weight : str, optional
        name of the edge attribute used as weight (for computing distance
        between nodes);
        by default (`None`) : all edges have a weight of 1

    node_attr_mode : str {'mean', 'center', 'medoid'} (or list of strs), default: 'mean'
        string or list of strings, the strings indicate how the node attributes 
        are computed in the reduced graph (see above)

        - if string: the same operation (mode) is used for all considered \
        node attributes
        - if list of strings: the length of the list must be equal to the length \
        of the list `node_attr_list`, and `node_attr_mode[i]` indicates the \
        operation (mode) used for the node attribute `node_attr_list[i]`

    node_attr_list : list of strs, optional
        list of node attributes of the input graph to be included in 
        the reduced graph;
        by default (`None`): the list of all nodes attributes is considered
    
    node_group_id_dict : dict, optional
        dictionary of "group id", the keys are the input graph node labels, and
        the values an integer group id; the id 0 is the default id attached to 
        any node that is not listed in the keys; nodes with different group
        ids will be grouped in different clusters;
        by default (`None`): not applied, i.e. no extra constraint on clusters
   
    return_new_id_dict : bool, default: False
        if `True`, the dictionary `new_id_dict` is returned, where
        the keys are the node ids (labels) in the input graph, and the values
        are the corresponding node ids in the output graph (cluster ids)
        
    return_orig_id_dict : bool, default: False
        if `True`, the dictionary `orig_id_dict` is returned, where
        the keys are the node ids (labels) in the output graph (cluster), and the 
        values are the lists of the corresponding node ids (labels) in the input graph
        (original ids)
        
    return_center_id_dict : bool, default: False
        if `True`, the dictionary `center_id_dict` is returned, where
        the keys are the node ids (labels) in the output graph (cluster), and the 
        values are the node ids (labels) in the input graph of the cluster centers
        (according to the MEMB algorithm)

    return_medoid_id_dict : bool, default: False
        if `True`, the dictionary `medoid_id_dict` is returned, where
        the keys are the node ids (labels) in the output graph (cluster), and the 
        values are the node ids (labels) in the input graph of the cluster medoids;
        the medoid of the subgraph of the input graph induced by the nodes in a
        considered cluster is computed, where the medoid is the node realizing the
        minimum of the sum of the shortest path lengths (wrt. edge weight) to any 
        other node

    return_cluster_diameter_dict : bool, default: False
        if `True`, the dictionary `cluster_diameter_dict` is returned, where
        the keys are the node ids (labels) in the output graph (cluster), and the 
        values are the diameter of the cluster (box), i.e. the maximal
        distance (in the original graph, wrt. edge weight) between two nodes in 
        the cluster
    
    seed : int, optional
        seed number for initializing the random number generator

    verbose : int, default: 0
        verbose mode, larger value implies more info printed

    Returns
    -------
    G_red : networkx.Graph
        reduced graph (see above)
    
    new_id_dict : dict, optional
        returned if `return_new_id_dict=True`: dictionary giving the new id
        (label of node in the output graph) (value) for each original id (label of 
        node in the input graph) (key) (see `return_new_id_dict` above);
        note: this property can be attached to nodes of the input graph `G` with:
        `nx.set_node_attributes(G, new_id_dict, 'new_id')`

    orig_id_dict : dict, optional
        returned if `return_orig_id_dict=True`: dictionary giving the list of 
        original ids (labels of nodes in the input graph) (value) for each new id 
        (label of node in the output graph) (key) (see `return_orig_id_dict` above);
        note: this property can be attached to nodes of the output graph `G_red` with:
        `nx.set_node_attributes(G_red, orig_id_dict, 'orig_id')`

    center_id_dict : dict, optional
        returned if `return_center_id_dict=True`: dictionary giving the MEMB center id
        of the cluster (box) (label of node in the input graph) (value) for each new id 
        (label of node in the output graph) (key) (see `return_center_id_dict` above);
        note: this property can be attached to nodes of the output graph `G_red` with:
        `nx.set_node_attributes(G_red, center_id_dict, 'center_id')`

    medoid_id_dict : dict, optional
        returned if `return_medoid_id_dict=True`: dictionary giving the medoid
        id of the cluster (box) (label of node in the input graph) (value), for each 
        new id (label of node in the output graph) (key) (see `return_medoid_id_dict` 
        above);
        note: this property can be attached to nodes of the output graph `G_red` with:
        `nx.set_node_attributes(G_red, medoid_id_dict, 'medoid_id')`

    cluster_diameter_dict : dict, optional
        returned if `return_cluster_diameter_dict=True`: dictionary giving the diameter 
        of the cluster (box) (value), for each new id (label of node in the output graph) 
        (key) (see `return_cluster_diameter_dict` above);
        note: this property can be attached to nodes of the output graph `G_red` with:
        `nx.set_node_attributes(G_red, cluster_diameter_dict, 'cluster_diameter')`

    References
    ----------
    - C. Song, L. K. Gallos, S. Havlin, and H. A. Makse (2007), \
    How to calculate the fractal dimension of a complex network: \
    the box covering algorithm, doi: 10.1088/1742-5468/2007/03/P03006}
    """
    if node_group_id_dict is not None:
        for x in G.nodes():
            if x not in node_group_id_dict.keys():
                node_group_id_dict[x] = 0

        edges_rm_list = [e for e in G.edges() if node_group_id_dict[e[0]] != node_group_id_dict[e[1]]]
        G_copy = G.copy()
        G_copy.remove_edges_from(edges_rm_list)
        # gc.remove_edges_from(edges_list)
        # gc = nx.connected_components(gc)
        # return gc
    else:
        G_copy = G
        edges_rm_list = []

    if seed is not None:
        np.random.seed(seed)

    # Check node attributes of the input graph to be considered in the reduced graph
    # and corresponding operation (mode)

    # Get keys (attribute names) in the input graph
    # # - all keys
    # keys = [list(G_copy.nodes[i].keys()) for i in G_copy.nodes()] # list of lists
    # node_attr_all = np.unique([xij for xi in keys for xij in xi])
    # # - from one node
    # node_attr_all = G_copy.nodes[list(G_copy.nodes())[0]].keys()
    node_attr_all = kn.utils.get_node_attribute_names(G_copy)

    # Check given node attributes name
    if node_attr_list is not None:
        if not np.all([k in node_attr_all for k in node_attr_list]):
            raise ValueError('Attribute name does not exist, check `attr_node_list` parameters')
    else:
        node_attr_list = node_attr_all

    # Check given mode
    if isinstance(node_attr_mode, list):
        if len(node_attr_mode) != len(node_attr_list):
            raise ValueError('Length of the list `node_attr_mode` is not valid')

        if np.any([s not in ('center', 'mean', 'medoid') for s in node_attr_mode]):
            raise ValueError('Entry of the list `node_attr_mode` is not valid')

    else: # `node_attr_mode` assumed to be a string
        if node_attr_mode not in ('center', 'mean', 'medoid'):
            raise ValueError('`node_attr_mode` is not valid')
        node_attr_mode = len(node_attr_list)*[node_attr_mode]


    # Step 1: Identifying the centers in the original graph
    # ------

    # (i) All nodes are marked as uncovered and non-centers.
    nodes_uncovered_dict = {id:True for id in G_copy.nodes()}
    nodes_center_dict = {id:False for id in G_copy.nodes()}

    if edge_weight is None:
        if verbose > 0:
            print(f'MEMB (unweighted graph) - Computing shortest path length for any pair of nodes (with cutoff {r})...')
        cost_from = dict(nx.all_pairs_shortest_path_length(G_copy, cutoff=r))
    else:
        if verbose > 0:
            print(f'MEMB (weighted graph) - Computing shortest path length for any pair of nodes (with cutoff {r})...')
        cost_from = dict(nx.all_pairs_dijkstra_path_length(G_copy, cutoff=r, weight=edge_weight))

    n_nodes_uncovered = G_copy.number_of_nodes() # or any positive integer
    
    if verbose > 0:
        print('MEMB - Selecting node centers...')
    
    # Loop until all nodes are covered
    n_centers = 0
    while n_nodes_uncovered:
        if verbose > 1:
            print(f'... Number of centers: {n_centers:5d}, Number of uncovered nodes:  {n_nodes_uncovered:5d}')

        # (ii) For all non-center nodes x (including already covered nodes),
        # the excluded mass is computed, i.e. the number of uncovered nodes
        # in the cluster centered at x
        nodes_em_dict = {id:0 for id in G_copy.nodes()}

        for id in G_copy.nodes():
            if nodes_center_dict[id]:
                continue

            nodes_em_dict[id] = len([j for j in cost_from[id].keys() if nodes_uncovered_dict[j]])

        # print('...', 'nodes_em_dict', nodes_em_dict)

        # Select a node c realizing the maximum of the excluded mass, and mark c as a center
        em_max = np.asarray(list(nodes_em_dict.values())).max()
        nodes_em_max_list = [k for k, v in nodes_em_dict.items() if v == em_max]
        if len(nodes_em_max_list) > 1:
            id = nodes_em_max_list[np.random.randint(len(nodes_em_max_list))]
        else: # len(nodes_em_max_list) == 1
            id = nodes_em_max_list[0]
        nodes_center_dict[id] = True
        n_centers = n_centers + 1

        # print('...', 'new center', id)
        
        # (iii) Mark all the nodes in the cluster centered at c (i.e. nodes 
        # at a distance from c less than or equal to r) as covered
        for j in cost_from[id].keys():
            nodes_uncovered_dict[j] = False

        # Update the number of uncovered nodes
        n_nodes_uncovered = np.asarray(list(nodes_uncovered_dict.values())).sum()

    if verbose > 0:
        print(f'MEMB - Number of centers: {n_centers:d}')

    # print('n_nodes_uncovered', n_nodes_uncovered)
    # print('nodes_uncovered_dict', nodes_uncovered_dict)
    # print('nodes_center_dict', nodes_center_dict)

    # Step 2. Creating the clusters.
    # -------
        
    if verbose > 0:
        print('MEMB - Assigning new id to node centers...')
    
    # (i) Assign an id to each center.
    nodes_center_list = [k for k, v in nodes_center_dict.items() if v]
    nodes_new_id_dict = {i:-1 for i in G_copy.nodes()}
    for i, k in enumerate(nodes_center_list):
        nodes_new_id_dict[k] = i

    if verbose > 0:
        print('MEMB - Computing distance to node centers for all nodes...')

    # (ii) For all nodes, compute the central distance, i.e. the distance
    # to the nearest center.    
    nodes_central_dist_dict = {id:np.asarray([v for j, v in cost_from[id].items() if nodes_center_dict[j]]).min() for id in G_copy.nodes()}

    # (iii) Sort the nodes according to increasing central distance
    nodes_central_dist_dict = {k:v for k, v in sorted(nodes_central_dist_dict.items(), key=lambda item: item[1])}

    if verbose > 0:
        print('MEMB - Assigning new id for non-center nodes...')

    # (iv) For each non-center node u (visiting in order such that the central 
    # distance is increasing), select one of its neighbors, v, with a smaller 
    # central distance, and assign the id of v to u.
    if edge_weight is None:
        for id, d in nodes_central_dist_dict.items():
            if nodes_center_dict[id]:
                continue
            
            id_smaller_dist_list = [j for j in G_copy.neighbors(id) if nodes_central_dist_dict[j] < d]

            if len(id_smaller_dist_list) > 1:
                id_smaller_dist = id_smaller_dist_list[np.random.randint(len(id_smaller_dist_list))]
            else: # len(id_smaller_dist_list) == 1
                id_smaller_dist = id_smaller_dist_list[0]
            nodes_new_id_dict[id] = nodes_new_id_dict[id_smaller_dist]
    
    else:
        for id, d in nodes_central_dist_dict.items():
            if nodes_center_dict[id]:
                continue
        
            id_smaller_dist_list = [j for j in G_copy.neighbors(id) if nodes_central_dist_dict[j] < d and nodes_central_dist_dict[j] + G_copy.edges[id, j][edge_weight] <= r]

            if len(id_smaller_dist_list) > 1:
                id_smaller_dist = id_smaller_dist_list[np.random.randint(len(id_smaller_dist_list))]
            else: # len(id_smaller_dist_list) == 1
                id_smaller_dist = id_smaller_dist_list[0]
            nodes_new_id_dict[id] = nodes_new_id_dict[id_smaller_dist]
    
    # Set the list of original ids for each new id (in a dictionary)
    # (returned optionally, and used for computing mean or medoid of node attributes over cluster (see below))
    nodes_orig_id_dict = {i:[] for i in range(n_centers)}
    for id in G_copy.nodes():
        nodes_orig_id_dict[nodes_new_id_dict[id]].append(id)

    if verbose > 0:
        print('MEMB - Building the reduced graph...')

    # Step 3. Building the reduced graph.
    # -------
    # WORK ON G FROM HERE (NOT G_copy)!
    G_red = nx.Graph()
    # Set nodes (clusters in the reduced graph)
    G_red.add_nodes_from(range(n_centers))
    # Set edges (links between clusters in the reduced graph)
    for u, v in G.edges():
        if nodes_new_id_dict[u] != nodes_new_id_dict[v]:
            G_red.add_edge(nodes_new_id_dict[u], nodes_new_id_dict[v])

    # # Set cluster (box) diameter if needed
    # if return_cluster_diameter_dict:
    #     cluster_diameter_dict = {i: nx.diameter(nx.subgraph(G, nodes_orig_id_dict[i]), weight=edge_weight) for i in G_red.nodes()}

    # # Set node attributes in the reduced graph (according to `node_attr_list` and `node_attr_mode`)
    # for attr, mode in zip(node_attr_list, node_attr_mode):
    #     d = nx.get_node_attributes(G, attr, default=np.nan) # dictionary original_id:value_of_attribute
    #     if mode == 'center':
    #         for c in nodes_center_list:
    #             G_red.nodes[nodes_new_id_dict[c]][attr] = d[c]
    #     elif mode == 'mean':
    #         v0 = list(d.values())[0] # value of one node in G
    #         if hasattr(v0, '__len__'):
    #             # attr_type = type(v0)[0]     # attribute type
    #             for i in G_red.nodes():
    #                 # G_red.nodes[i][attr] = attr_type(np.nanmean(np.asarray([np.atleast_1d(d[j]) for j in nodes_orig_id_dict[i]]), axis=0))
    #                 G_red.nodes[i][attr] = np.nanmean(np.asarray([np.atleast_1d(d[j]) for j in nodes_orig_id_dict[i]]), axis=0)
    #         else:
    #             for i in G_red.nodes():
    #                 G_red.nodes[i][attr] = np.nanmean(np.asarray([np.atleast_1d(d[j]) for j in nodes_orig_id_dict[i]]), axis=0)[0]

    # Compute medoid of each cluster (box) and set cluster (box) diameter, if needed
    compute_medoid_dict = return_medoid_id_dict or 'medoid' in node_attr_mode
    if compute_medoid_dict or return_cluster_diameter_dict:
        if compute_medoid_dict:
            medoid_id_dict = {} # initialization of the dictionary
        
        if return_cluster_diameter_dict:
            cluster_diameter_dict = {} # initialization of the dictionary

        for i in G_red.nodes():
            G_cluster = nx.subgraph(G, nodes_orig_id_dict[i])
        
            if compute_medoid_dict:
                d_sum_min = np.inf
                for (u, dist_from_u) in nx.all_pairs_dijkstra_path_length(G_cluster, weight=edge_weight):
                    d_sum = np.asarray(list(dist_from_u.values())).sum()
                    if d_sum < d_sum_min:
                        medoid = u
                        d_sum_min = d_sum

                medoid_id_dict[i] = medoid

            if return_cluster_diameter_dict:
                cluster_diameter_dict[i] = nx.diameter(G_cluster, weight=edge_weight)

    # Set node attributes in the reduced graph (according to `node_attr_list` and `node_attr_mode`)
    for attr, mode in zip(node_attr_list, node_attr_mode):
        d = nx.get_node_attributes(G, attr, default=np.nan) # dictionary original_id:value_of_attribute
        if mode == 'center':
            for c in nodes_center_list:
                G_red.nodes[nodes_new_id_dict[c]][attr] = d[c]

        elif mode == 'mean':
            v0 = list(d.values())[0] # value of one node in G
            if hasattr(v0, '__len__'):
                # attr_type = type(v0)[0]     # attribute type
                for i in G_red.nodes():
                    # G_red.nodes[i][attr] = attr_type(np.nanmean(np.asarray([np.atleast_1d(d[j]) for j in nodes_orig_id_dict[i]]), axis=0))
                    G_red.nodes[i][attr] = np.nanmean(np.asarray([np.atleast_1d(d[j]) for j in nodes_orig_id_dict[i]]), axis=0)
            else:
                for i in G_red.nodes():
                    G_red.nodes[i][attr] = np.nanmean(np.asarray([np.atleast_1d(d[j]) for j in nodes_orig_id_dict[i]]), axis=0)[0]
        
        elif mode == 'medoid':
            for i in G_red.nodes():
                G_red.nodes[i][attr] = d[medoid_id_dict[i]]

    out = [G_red]
    if return_new_id_dict:
        out.append(nodes_new_id_dict)
    if return_orig_id_dict:
        out.append(nodes_orig_id_dict)
    if return_center_id_dict:
        center_id_dict = {i:u for u, i in nodes_new_id_dict.items() if nodes_center_dict[u]}
        out.append(center_id_dict)
    if return_medoid_id_dict:
        out.append(medoid_id_dict)
    if return_cluster_diameter_dict:
        out.append(cluster_diameter_dict)
    
    if len(out) == 1:
        out = out[0]
    else:
        out = tuple(out)

    return out
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
def graph_fractal_dimension_with_memb(
        G,
        r_values=np.arange(1, 101),
        min_n_edges_fraction=0.02,
        min_n_edges_lim=3,
        data_points_removed_first_fraction=0.33,
        data_points_removed_last_fraction=0.0,
        confidence_level=0.95,
        pos_attr='pos',
        memb_edge_weight=None,
        memb_node_attr_mode='mean',
        memb_node_attr_list=None,
        memb_node_group_id_dict=None,
        return_poly_fit_params=True,
        return_poly_fit_params_cov=True,
        return_index_used_for_fit=True,
        return_G_red_mean_edge_length_list=True,
        return_G_red_n_nodes_list=True,
        return_G_red_n_edges_list=True,
        return_G_red_r_list=False,
        return_G_red_list=False,
        seed=None,
        verbose=1,
        make_plot=True,
        plot_in_log_scale=True):
    """
    Computes the fractal dimension of a graph by using MEMB algorithm.

    This function computes a sequence of reduced graphs by applying MEMB to 
    the given graph (with increasing radius) (function `reduce_graph_memb`); 
    then, assuming the relationship
    
    .. math::
        N \propto \\lambda^{-d_f}

    where :math:`N` is number of graph nodes and :math:`\\lambda` the mean edge 
    length in reduced graph, the fractal dimension :math:`d_f` is obtained by 
    fitting a line on the log-log plot of :math:`N` as function of :math:`\\lambda`.

    Parameters
    ----------
    G : networkx.Graph
        input graph
    
    r_values : sequence of ints or floats
        sequence of radius values for reducing input graph, the sequence should
        be in ascending order);
        see function `reduce_graph_memb`
        
    min_n_edges_fraction : float, default: 0.02
        see `min_n_edges_lim`

    min_n_edges_lim : int, default: 3
        with `m = max(min_n_edges_lim, int(min_n_edges_fraction * G_n_edges))`, where 
        `G_n_edges` is the number of edges in the input graph:
        as soon as the reduced graph has a number of edges smaller than `m`, this 
        graph is not retained and the sequence of reduced graph is stopped (hence, it 
        is important that the sequence `r_values` is in ascending order)

    data_points_removed_first_fraction : float, default : 0.33
        see `data_points_removed_last_fraction`

    data_points_removed_last_fraction : float, default : 0.0
        `data_points_removed_first_fraction` and `data_points_removed_last_fraction` are
        positive numbers in (0, 1); sequences of mean edge lengths and number of nodes,
        including input graph and successive reduced graphs constitute the data
        points; 
        the fraction `data_points_removed_first_fraction`, resp.
        `data_points_removed_last_fraction`, at the beginning, resp. end, of the
        data points are removed before doing the fit (line on the log-log plot);
        this allows to deal with the "extreme" data points (the sequence  `r_values` 
        should be in ascending order)

    confidence_level : float, default: 0.95
        confidence level, float in the interval (0, 1), to compute the 
        uncertainty for the line fitting (and for the fractal dimension)

    pos_attr : str, default: 'pos'
        name of the node attribute for position

    memb_edge_weight : str, optional
        name of the edge attribute used as weight for the MEMB algorithm;
        parameter `edge_weight` passed to the function `reduce_graph_memb`

    memb_node_attr_mode : str {'mean', 'center'} (or list of strs), default: 'mean'
        string or list of strings, the strings indicate how the node attributes 
        are computed for the MEMB algorithm;
        parameter `node_attr_mode` passed to the function `reduce_graph_memb`

    memb_node_attr_list : list of strs, optional
        list of node attributes of the input graph to be included in 
        the reduced graphs for the MEMB algorithm;
        parameter `node_attr_list` passed to the function `reduce_graph_memb`;
        by default (`None`): `memb_node_attr_list` is set [`pos_attr`]

    memb_node_group_id_dict : dict, optional
        dictionary of "group id", the keys are the input graph node labels, and
        the values an integer id; the id 0 is the default id attached to 
        any node that is not listed in the keys; nodes with different group
        ids will be in grouped in different clusters;
        parameter `node_group_id_dict` passed to the function `reduce_graph_memb`;

    return_poly_fit_params : bool, default: True
        if `True`, the array (shape (2, )) of parameters of the fitted 
        polynom (line) is returned

    return_poly_fit_params_cov : bool, default: True
        if `True`, the array (shape (2, 2)) of parameters covariance of the 
        fitted polynom (line) is returned
    
    return_index_used_for_fit: bool, default: True
        if `True`, the starting (included) and ending (excluded) indices of the sequence 
        of data points (from original graph followed by successive reduced graph) used 
        for fitting is returned
    
    return_G_red_mean_edge_length_list : bool, default: True
        if `True`, the list of mean edge length of all computed reduced graphs 
        is returned; note: the mean edge length of the original graph (`G`) is 
        also returned

    return_G_red_n_nodes_list : bool, default: True
        if `True`, the list of number of nodes of all computed reduced graphs
        is returned; note: the number of nodes of the original graph (`G`) is 
        also returned

    return_G_red_n_edges_list : bool, default: True
        if `True`, the list of number of edges of all computed reduced graphs
        is returned; note: the number of edges of the original graph (`G`) is 
        also returned

    return_G_red_r_list : bool, default: False
        if `True`, the list of radius values for computing reduced graphs 
        is returned

    return_G_red_list : bool, default: False
        if `True`, the list of all computed reduced graphs is returned
    
    seed : int, optional
        seed number for initializing the random number generator

    verbose : int, default: 1
        verbose mode, larger value implies more info printed

    make_plot : bool, default: True
        indicates if the plot of the data points and fitting curve is displayed 
        (in the current figure axis)

    plot_in_log_scale : bool, default: True
        used if `make_plot=True`, indicates if the plot is in log-log scale
        (along x- and y-axes)

    Returns
    -------
    out_dict : dict
        dictionary with the following keys / values:
        
        - df : float
            fractal dimension
    
        - df_delta : float
            uncertainty for fractal dimension, the interval `df +/- df_delta` 
            corresponds to the confidence interval at a confidence level 
            `confidence_level` derived from Gaussian distribution
    
        - confidence_level : float
            confidence level used to compute the uncertainty for the line fitting 
            
        - poly_fit_params : array of shape (2, ), optional
            returned if `return_poly_fit_params=True`:
            array of parameters of the fitted polynom (line); in particular:
            `df = -poly_fit_params[0]`

        - poly_fit_params_cov : array of shape (2, 2), optional
            returned if `return_poly_fit_params_cov=True`:
            array (shape (2, 2)) of parameters covariance of the fitted polynom 
            (line); in particular: 
            `df_delta = scipy.stats.norm.ppf((1+confidence_level)/2) * np.sqrt(poly_fit_params_cov[0, 0])`
        
        - index_used_for_fit : list of two ints, optional
            returned if `return_index_used_for_fit=True`:
            starting (included) index and ending (excluded) index of the sequence of
            data points (from original graph followed by successive reduced graph) used for fitting
        
        - G_mean_edge_length : float, optional
            returned if `return_G_red_mean_edge_length_list=True`: 
            mean edge length of the input graph

        - G_n_nodes : float, optional
            returned if `return_G_red_n_nodes_list=True`: 
            number of nodes of the input graph

        - G_n_edges : float, optional
            returned if `return_G_red_n_edges_list=True`: 
            number of nodes of the input graph

        - G_red_mean_edge_length_list : list of floats, optional
            returned if `return_G_red_mean_edge_length_list=True`: 
            the list of mean edge length of graphs used to compute the 
            fractal dimension (by line fitting on the log-log plot (see above)) is 
            returned; note that if `account_for_original_graph=True`: the entry at
            index 0 corresponds to the input graph, and the next entries to the 
            reduced graphs

        - G_red_n_nodes_list : list of ints, optional
            returned if `return_G_red_n_nodes_list=True`: 
            the list of number of nodes of reduced graphs

        - G_red_n_edges_list : list of ints, optional
            returned if `return_G_red_n_edges_list=True`:
            the list of number of edges of reduced graphs

        - G_red_r_list : bool, list of ints or floats, optional
            returned if `return_G_red_r_list=True`:
            the list of radius values for computing reduced graphs

        - G_red_list : list of networkx.Graph, optional
            returned if `return_G_red_list=True`: 
            the list of reduced graphs
    """
    if seed is not None:
        np.random.seed(seed)

    if memb_node_attr_list is None:
        memb_node_attr_list = [pos_attr]

    G_n_nodes = G.number_of_nodes()
    G_n_edges = G.number_of_edges()
    min_n_edges = max(min_n_edges_lim, int(min_n_edges_fraction * G.number_of_edges()))
    
    if verbose > 0:
        print(f'Compute Reduced graphs (MEMB); original graph: #nodes = {G_n_nodes}, #edges = {G_n_edges}; min #edges for red. graph: {min_n_edges}')

    G_mean_edge_length = kn.utils.mean_edge_length(G, pos_attr=pos_attr)
    # G_mean_edge_length = \
    #     np.mean(np.sqrt(
    #         np.sum(np.asarray(
    #             [np.asarray(G.nodes[u][pos_attr]) - np.asarray(G.nodes[v][pos_attr]) for u, v in G.edges()]
    #         )**2, axis=1)
    #     ))

    G_red_mean_edge_length_list = []
    G_red_n_nodes_list = []

    if return_G_red_n_edges_list:
        G_red_n_edges_list = []

    if return_G_red_r_list:
        G_red_r_list = []
    
    if return_G_red_list:
        G_red_list = []
        
    t1_all = time.time()
    for r in r_values:
        t1 = time.time()
        G_red = reduce_graph_memb(
            G, 
            r=r, 
            edge_weight=memb_edge_weight,
            node_attr_mode=memb_node_attr_mode,
            node_attr_list=memb_node_attr_list,
            node_group_id_dict=memb_node_group_id_dict,
            return_new_id_dict=False, 
            return_orig_id_dict=False,
            return_center_id_dict=False, 
            return_medoid_id_dict=False,
            return_cluster_diameter_dict=False,
            seed=None,
            verbose=verbose-1)
        t2 = time.time()

        G_red_n_nodes = G_red.number_of_nodes()
        G_red_n_edges = G_red.number_of_edges()
        if verbose > 0:
            print(f'MEMB, r = {r:.3g}: #nodes = {G_red_n_nodes}, #edges = {G_red_n_edges}, elapsed time : {t2-t1:.2f} sec')

        if G_red_n_edges < min_n_edges:
            if verbose > 0:
                print(f'Stop MEMB (#edges = {G_red_n_edges} < {min_n_edges})')
            break

        G_red_mean_edge_length_list.append(kn.utils.mean_edge_length(G_red, pos_attr=pos_attr))
        # G_red_mean_edge_length_list.append(
        #     np.mean(np.sqrt(
        #         np.sum(np.asarray(
        #             [np.asarray(G_red.nodes[u][pos_attr]) - np.asarray(G_red.nodes[v][pos_attr]) for u, v in G_red.edges()]
        #         )**2, axis=1)
        #     ))
        # )
        
        G_red_n_nodes_list.append(G_red_n_nodes)
        
        if return_G_red_n_edges_list:
            G_red_n_edges_list.append(G_red_n_edges)
    
        if return_G_red_r_list:
            G_red_r_list.append(r)

        if return_G_red_list:
            G_red_list.append(G_red)
        
    t2_all = time.time()
    if verbose > 0:
        print(f'Total elapsed time (MEMBs): {t2_all-t1_all:.2f} sec')

    if verbose > 0:
        print(f'Compute fractal dimension (line fitting on log-log plot of number of nodes as function of mean edge length)')

    # Compute fractal dimension (line fitting on log-log plot of number of nodes as function of mean edge length)
    # -----------------------------------------------------------------------------------------------------------
    # Data (all)
    x0 = G_mean_edge_length
    y0 = G_n_nodes
    x_all = np.asarray([x0] + G_red_mean_edge_length_list)
    y_all = np.asarray([y0] + G_red_n_nodes_list)

    poly_fit_params, poly_fit_params_cov, x, y, n1, n2 = _loglog_fit(x_all, y_all, data_points_removed_first_fraction, data_points_removed_last_fraction, verbose=verbose)

    # Get fractal dimension df, and df_delta
    if poly_fit_params is not None:
        df = - poly_fit_params[0] # - slope
        fit_ok = True
    else:
        df = np.nan
        fit_ok = False
    
    if poly_fit_params_cov is not None:
        t = scipy.stats.norm.ppf((1+confidence_level)/2)
        df_delta = t * np.sqrt(poly_fit_params_cov[0, 0]) # slope_delta
    else:
        df_delta = np.nan

    # Make plot (if needed)
    if make_plot:
        plot_fit_fractal_dimension(
            x_all, y_all, x, y, x0, y0, df, df_delta, confidence_level, poly_fit_params, 
            plot_fit=fit_ok, plot_in_log_scale=plot_in_log_scale)

    # Set output dictionary
    out_dict = {'df': df, 'df_delta': df_delta, 'confidence_level': confidence_level}

    if return_poly_fit_params:
        out_dict['poly_fit_params'] = poly_fit_params

    if return_poly_fit_params_cov:
        out_dict['poly_fit_params_cov'] = poly_fit_params_cov

    if return_index_used_for_fit:
        out_dict['index_used_for_fit'] = [n1, n2]

    if return_G_red_mean_edge_length_list:
        out_dict['G_red_mean_edge_length_list'] = G_red_mean_edge_length_list
        out_dict['G_mean_edge_length'] = G_mean_edge_length

    if return_G_red_n_nodes_list:
        out_dict['G_red_n_nodes_list'] = G_red_n_nodes_list
        out_dict['G_n_nodes'] = G_n_nodes

    if return_G_red_n_edges_list:
        out_dict['G_red_n_edges_list'] = G_red_n_edges_list
        out_dict['G_n_edges'] = G_n_edges

    if return_G_red_r_list:
        out_dict['G_red_r_list'] = G_red_r_list
    
    if return_G_red_list:
        out_dict['G_red_list'] = G_red_list

    return out_dict
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
def graph_fractal_dimension_with_memb_update(
        out_dict,
        data_points_removed_first_fraction=0.33,
        data_points_removed_last_fraction=0.0,
        confidence_level=0.95,
        verbose=1,
        make_plot=True,
        plot_in_log_scale=True):
    """
    Updates the fit (power law) for the fractal dimension.

    It is assumed that the function `graph_fractal_dimension_with_memb` has 
    already been run, and that the output dictionary `out_dict` contains the keys 
    `G_red_mean_edge_length_list`, `G_red_n_nodes_list`, `G_mean_edge_length`, 
    `G_n_nodes`; 
    the fit is updated by changing the fractions of data points removed at the 
    beginning and end of the sequences of data points.

    See function `graph_fractal_dimension_with_memb` for details.

    Note: this function operates inplace on the `out_dict` dictionary (it is
    updated and returned).

    Parameters
    ----------
    out_dict : dict
        output dictionary of the function `graph_fractal_dimension_with_memb`;
        the dictionary should contain the keys `G_red_mean_edge_length_list`, 
        `G_red_n_nodes_list`, `G_mean_edge_length`, `G_n_nodes`
    
    data_points_removed_first_fraction : float, default : 0.33
        see `data_points_removed_last_fraction`

    data_points_removed_last_fraction : float, default : 0.0
        `data_points_removed_first_fraction` and `data_points_removed_last_fraction` are
        positive numbers in (0, 1); sequences of mean edge lengths and number of nodes,
        including input graph and successive reduced graphs constitute the data
        points; 
        the fraction `data_points_removed_first_fraction`, resp.
        `data_points_removed_last_fraction`, at the beginning, resp. end, of the
        data points are removed before doing the fit (line on the log-log plot);
        this allows to deal with the "extreme" data points (the sequence  `r_values` 
        should be in ascending order)

    confidence_level : float, default: 0.95
        confidence level, float in the interval (0, 1), to compute the 
        uncertainty for the line fitting (and for the fractal dimension)

    verbose : int, default: 1
        verbose mode, larger value implies more info printed

    make_plot : bool, default: True
        indicates if the plot of the data points and fitting curve is displayed 
        (in the current figure axis)

    plot_in_log_scale : bool, default: True
        used if `make_plot=True`, indicates if the plot is in log-log scale
        (along x- and y-axes)

    Returns
    -------
    out_dict : dict
        updated output dictionary
    """
    if 'G_red_mean_edge_length_list' not in out_dict.keys() or \
       'G_red_n_nodes_list' not in out_dict.keys() or \
       'G_mean_edge_length' not in out_dict.keys() or \
       'G_n_nodes' not in out_dict.keys():
        raise ValueError('The input dictionary `out_dict` should contain the keys `G_red_mean_edge_length_list`, `G_red_n_nodes_list`, `G_mean_edge_length`, `G_n_nodes`')

    if verbose > 0:
        print(f'Compute fractal dimension (line fitting on log-log plot of number of nodes as function of mean edge length)')

    # Compute fractal dimension (line fitting on log-log plot of number of nodes as function of mean edge length)
    # -----------------------------------------------------------------------------------------------------------
    # Data (all)
    x0 = out_dict['G_mean_edge_length']
    y0 = out_dict['G_n_nodes']
    x_all = np.asarray([x0] + out_dict['G_red_mean_edge_length_list'])
    y_all = np.asarray([y0] + out_dict['G_red_n_nodes_list'])

    poly_fit_params, poly_fit_params_cov, x, y, n1, n2 = _loglog_fit(x_all, y_all, data_points_removed_first_fraction, data_points_removed_last_fraction, verbose=verbose)

    # Get fractal dimension df, and df_delta
    if poly_fit_params is not None:
        df = - poly_fit_params[0] # - slope
        fit_ok = True
    else:
        df = np.nan
        fit_ok = False
    
    if poly_fit_params_cov is not None:
        t = scipy.stats.norm.ppf((1+confidence_level)/2)
        df_delta = t * np.sqrt(poly_fit_params_cov[0, 0]) # slope_delta
    else:
        df_delta = np.nan

    # Make plot (if needed)
    if make_plot:
        plot_fit_fractal_dimension(
            x_all, y_all, x, y, x0, y0, df, df_delta, confidence_level, poly_fit_params, 
            plot_fit=fit_ok, plot_in_log_scale=plot_in_log_scale)

    # Update output dictionary
    out_dict['df'] = df
    out_dict['df_delta'] = df_delta
    out_dict['confidence_level'] = confidence_level

    if 'poly_fit_params' in out_dict.keys():
        out_dict['poly_fit_params'] = poly_fit_params

    if 'poly_fit_params_cov' in out_dict.keys():
        out_dict['poly_fit_params_cov'] = poly_fit_params_cov

    if 'index_used_for_fit' in out_dict.keys():
        out_dict['index_used_for_fit'] = [n1, n2]

    return out_dict
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
def plot_fit_fractal_dimension(
        x_all,
        y_all,
        x,
        y,
        x0,
        y0,
        df,
        df_delta,
        confidence_level,
        poly_fit_params,
        plot_fit=True,
        plot_in_log_scale=True):
    """
    Makes the plot (in current axis figure) of the data points and fitted curve for fractal dimension (power law).

    Parameters
    ----------
    x_all : array of floats
        x-coordinates of all data points (including original graph)

    y_all : array of floats
        y-coordinates of all data points (including original graph)

    x : array of floats
        x-coordinates of data points used for fitting (a subset of `x_all`)
    
    y : array of floats
        y-coordinates of data points used for fitting (a subset of `y_all`)

    x0 : float
        x-coordinate for the original graph
    
    y0 : float
        y-coordinate for the original graph

    df : float
        fractal dimension
    
    df_delta : float
        uncertainty for fractal dimension, the interval `df +/- df_delta` 
        corresponds to the confidence interval at a confidence level 
        `confidence_level` derived from Gaussian distribution
    
    confidence_level : float
        confidence level, float in the interval (0, 1), used to compute the 
        uncertainty for the line fitting (and for the fractal dimension)

    poly_fit_params : array of shape (2, ) or None
        array of parameters of the fitted polynom (line); in particular:
        `df = -poly_fit_params[0]`; if no fit was done, `poly_fit_params` is `None`

    plot_fit : bool ; default : True
        indicates if the plot of the fitted curve is done (in the current figure axis);
        if `False`, only the data points are plotted

    plot_in_log_scale : bool, default: True
        indicates if the plot is in log-log scale (along x- and y-axes)

    Returns
    -------
    None
    """
    df_lim_min = df - df_delta
    df_lim_max = df + df_delta
    
    # Estimation with the fitted line
    if plot_fit:
        xx = np.linspace(x.min(), x.max(), 200)
        yy = np.power(10, poly_fit_params[0] * np.log10(xx) + poly_fit_params[1])

        xx_all = np.linspace(x_all.min(), x_all.max(), 200)
        yy_all = np.power(10, poly_fit_params[0] * np.log10(xx_all) + poly_fit_params[1])

    # Plot parameters
    xlabel = '$\lambda$ : mean edge length'
    ylabel = '$N$ : number of nodes'

    title = '$N\propto \lambda^{-df}$' + f', df={df:.5f} in [{df_lim_min:.5f}, {df_lim_max:.5f}] [{100*confidence_level:.2g}% - inter.]'

    color_orig      = 'blue'
    marker_orig     = 'D'
    markersize_orig = 5

    color_data_all      = 'tab:blue'
    marker_data_all     = 'o'
    markersize_data_all = 5

    color_data_used      = 'tab:orange'
    marker_data_used     = 'o'
    markersize_data_used = 5

    color_fit_data_used = 'red'
    ls_fit_data_used    = 'solid'
    lw_fit_data_used    = 1.

    color_fit_data_all = 'red'
    ls_fit_data_all    = 'dotted'
    lw_fit_data_all    = 1.

    # Figure
    # ------
    # plt.figure(figsize=(10, 6))
    plt.plot(x_all, y_all, ls='', marker=marker_data_all , markersize=markersize_data_all , color=color_data_all , label='data')
    if len(x):
        plt.plot(x    , y    , ls='', marker=marker_data_used, markersize=markersize_data_used, color=color_data_used, label='data used for fitting')
    if plot_fit:
        plt.plot(xx_all, yy_all, ls=ls_fit_data_all , lw=lw_fit_data_all , color=color_fit_data_all)
        plt.plot(xx    , yy    , ls=ls_fit_data_used, lw=lw_fit_data_used, color=color_fit_data_used, label='fit')
    plt.plot(x0, y0, ls='', marker=marker_orig, markersize=markersize_orig, color=color_orig, label='original graph')
    
    if plot_in_log_scale:
        plt.xscale('log')
        plt.yscale('log')
    
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid()
    plt.legend()
    # plt.show()

    return
# ----------------------------------------------------------------------------

# =============================================================================
# Random walks on graphs
# =============================================================================

# ----------------------------------------------------------------------------
def graph_random_walk(
            G, 
            u0, 
            nsteps, 
            edge_weight=None, 
            seed=None):
    """
    Performs a random walk.

    The random walk starts at node `u0` and `nsteps` steps are done.
    At each step, the walker randomly chooses an edge from the currently
    visited node, the probability distribution is proportional to the 
    edge weights (all the same if not given).

    Parameters
    ----------
    G : networkx.Graph
        graph
    
    u0 : node key
        starting node in `G`
    
    nsteps : int
        number of time steps
    
    edge_weight : str, optional
        name of the edge attribute used as weight (probability to choose 
        an edge is proportional to its weight);
        by default (`None`) : all edges have a weight of 1

    seed : int, optional
        seed number for initializing the random number generator
    
    Returns
    -------
    rw : list of nodes
        list of nodes (labels) of length `nstep+1`, the random walk
    """
    if seed is not None:
        np.random.seed(seed)

    rw = [u0]
    u = u0

    r = np.random.random(nsteps)

    if edge_weight is None:
        for i in range(nsteps):
            neigh = [v for v in G.neighbors(u)]
            u = neigh[int(len(neigh)*r[i])]
            rw.append(u)

    else:
        for i in range(nsteps):
            neigh = [v for v in G.neighbors(u)]
            p = np.asarray([G.edges[u, v][edge_weight] for v in neigh]).cumsum()
            p = p / p[-1]
            u = neigh[np.where(r[i] < p)[0][0]]
            rw.append(u)

    return rw
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
def graph_random_walk_square_distance(
            G, 
            nrw, 
            nsteps, 
            edge_weight=None, 
            pos_attr='pos',
            u0=None, 
            seed=None,
            return_full_array=False,
            verbose=0):
    """
    Computes the mean and square distance after each time step, for several random walks.

    A total of `nrw` random walks of `nsteps` are performed, each one starts at node `u0` 
    if given or at a random node otherwise and the mean of Euclidean square distance from 
    the starting node, after each time step, is computed. At each step, the walker randomly 
    chooses an edge from the currently visited node, the probability distribution is 
    proportional to the edge weights (all the same if not given).

    Parameters
    ----------
    G : networkx.Graph
        graph
    
    nrw : int
        number of random walks

    nsteps : int
        number of time steps
    
    edge_weight : str, optional
        name of the edge attribute used as weight (probability to choose 
        an edge is proportional to its weight);
        by default (`None`) : all edges have a weight of 1

    pos_attr : str, default: 'pos'
        name of the attribute attached to nodes defining the position as a 
        sequence of 2 or 3 floats (used to compute distances)
    
    u0 : node key, optional
        starting node in `G` for all random walks;
        by default (`None`) a new node is randomly selected for every random walk

    seed : int, optional
        seed number for initializing the random number generator
    
    return_full_array : bool, default: False
        if `True`, the full array of shape `(nrw, nsteps)` of the 
        square distances from the starting node, after each time step, 
        of every random walk is returned.

    verbose : int, default: 0
        verbose mode, larger value implies more info printed

    Returns
    -------
    rw_square_dist_mean : 1d-array of floats of shape `(nsteps,)`
        mean square distance from the starting node, computed over all random 
        walks
    
    rw_square_dist : 2d-array of floats of shape `(nrw, nsteps)`, optional
        returned if `return_full_array=True`: square distances from the 
        starting node, after each time step, of every random walk
    """

    node_pos = {k:np.asarray(v) for k, v in nx.get_node_attributes(G, pos_attr).items()}

    if seed is not None:
        np.random.seed(seed)

    if u0 is None:
        node_list = list(node_pos.keys())
        ind = (len(node_list)*np.random.random(nrw)).astype('int')
        start_node_list = [node_list[i] for i in ind]
    else:
        start_node_list = nrw*[u0]

    rw_square_dist = np.zeros((nrw, nsteps))

    if verbose > 0:
        progress_old = 0
        progress = -1

    for k in range(nrw):
        if verbose > 0:
            progress = int(100*k/nrw)
            if progress > progress_old:
                print(f'Progress: {progress:3d}% ({k} random walk(s) of {nrw})')
                progress_old = progress

        u = start_node_list[k]
        pos_start = node_pos[u]
        rw = graph_random_walk(G, u, nsteps, edge_weight=edge_weight, seed=None)
        rw_square_dist[k] = np.sum((np.asarray([node_pos[v] for v in rw[1:]]) - pos_start)**2, axis=1)
    
    if verbose > 0:
        progress = 100
        print(f'Progress: {progress:3d}%')

    rw_square_dist_mean = rw_square_dist.mean(axis=0)

    if return_full_array:
        return rw_square_dist_mean, rw_square_dist

    return rw_square_dist_mean
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
def graph_random_walk_escape_time(
            G, 
            nrw, 
            dist, 
            nsteps_max=100000,
            edge_weight=None, 
            pos_attr='pos',
            u0=None, 
            seed=None,
            check_dist=False,
            verbose=0):
    """
    Computes the exit time for a random walk to escape an area, for several random walks.

    A total of `nrw` random walks are performed, each one starts at node `u0` 
    if given or at a random node otherwise and the walker stops as soon as the Euclidean
    distance between the currently visited node and the starting node is greater 
    than the given distance `dist`. At each step, the walker randomly chooses an edge 
    from the currently visited node, the probability distribution is proportional to the 
    edge weights (all the same if not given).

    Parameters
    ----------
    G : networkx.Graph
        graph
    
    nrw : int
        number of random walks

    dist : float (positive)
        limit distance, the walker stops as soon as the distance 
        from the starting node is greater than `dist`
    
    nsteps_max : int, default: 100000
        maximal number of time steps

    edge_weight : str, optional
        name of the edge attribute used as weight (probability to choose 
        an edge is proportional to its weight);
        by default (`None`) : all edges have a weight of 1

    pos_attr : str, default: 'pos'
        name of the attribute attached to nodes defining the position as a 
        sequence of 2 or 3 floats (used to compute distances)
    
    u0 : node key, optional
        starting node in `G` for all random walks;
        by default (`None`) a new node is randomly selected for every random walk

    seed : int, optional
        seed number for initializing the random number generator
    
    check_dist : bool, default: False
        if `True`, a check is done about the distance `dist`: if it is larger 
        than or equal to the half of the outbox diagonal (including the graph)
        then an error is raised, in order to avoid no escape within the 
        `nsteps_max` steps.

    verbose : int, default: 0
        verbose mode, larger value implies more info printed

    Returns
    -------
    escape_time_step : 1d-array of ints of shape `(nrw,)`
        first time step of escape, for each random walk
    """

    node_pos = {k:np.asarray(v) for k, v in nx.get_node_attributes(G, pos_attr).items()}

    if check_dist:
        # Check that dist < outbox_diagonal / 2
        G_outbox_diagonal = kn.utils.outbox_diagonal(G, pos_attr=pos_attr)
        if dist >= 0.5*G_outbox_diagonal:
            raise ValueError(f'`dist` should be smaller than half of outbox diagonal [outbox diagonal = {G_outbox_diagonal:.3g}]')
        
        # # Check that dist < d/2, where d is the maximal distance between any pair of nodes
        # all_pos = np.asarray(list(node_pos.values()))
        # d = np.asarray([np.sum((all_pos[i+1:,:] - all_pos[i])**2, axis=1).max() for i in range(all_pos.shape[0]-1)]).max()
        # d = np.sqrt(d)
        # if dist >= 0.5*d:
        #     raise ValueError(f"`dist` should be smaller than half of maximal distance between two nodes [max(d(u,v))={d:.3g}]")

    if seed is not None:
        np.random.seed(seed)

    if u0 is None:
        node_list = list(node_pos.keys())
        ind = (len(node_list)*np.random.random(nrw)).astype('int')
        start_node_list = [node_list[i] for i in ind]
    else:
        start_node_list = nrw*[u0]

    dist2 = dist**2
    escape_time_step = np.full((nrw,), np.nan)

    if verbose > 0:
        progress_old = 0
        progress = -1

    if edge_weight is None:
        for k in range(nrw):
            if verbose > 0:
                progress = int(100*k/nrw)
                if progress > progress_old:
                    print(f'Progress: {progress:3d}% ({k} random walk(s) of {nrw})')
                    progress_old = progress

            r = np.random.random(nsteps_max)
            u = start_node_list[k]
            pos_start = node_pos[u]

            for i in range(nsteps_max):
                neigh = [v for v in G.neighbors(u)]
                u = neigh[int(len(neigh)*r[i])]

                if np.sum((node_pos[u] - pos_start)**2) > dist2:
                    escape_time_step[k] = i+1
                    break

    else:
        for k in range(nrw):
            if verbose > 0:
                progress = int(100*k/nrw)
                if progress > progress_old:
                    print(f'Progress: {progress:3d}% ({k} random walk(s) of {nrw})')
                    progress_old = progress

            r = np.random.random(nsteps_max)
            u = start_node_list[k]
            pos_start = node_pos[u]

            for i in range(nsteps_max):
                neigh = [v for v in G.neighbors(u)]
                p = np.asarray([G.edges[u, v][edge_weight] for v in neigh]).cumsum()
                p = p / p[-1]
                u = neigh[np.where(r[i] < p)[0][0]]

                if np.sum((node_pos[u] - pos_start)**2) > dist2:
                    escape_time_step[k] = i+1
                    break

    if verbose > 0:
        progress = 100
        print(f'Progress: {progress:3d}%')
 
    return escape_time_step
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
def graph_random_walk_dimension(
        G,
        nrw,
        nsteps=None,
        nsteps_max=10000,
        edge_weight=None,
        pos_attr='pos',
        u0=None, 
        data_points_removed_first_fraction=0.33,
        data_points_removed_last_fraction=0.0,
        confidence_level=0.95,
        return_poly_fit_params=True,
        return_poly_fit_params_cov=True,
        return_index_used_for_fit=True,
        return_rw_square_dist_mean=True,
        return_rw_square_dist_full=False,
        seed=None,
        verbose=1,
        make_plot=True,
        plot_in_log_scale=True):
    """
    Computes the random walk dimension of a graph.

    This function performs several random walks in a graph, and retrieves
    the mean square distance from the strating node after each time step;
    at each step, the walker randomly chooses an edge from the currently visited 
    node, the probability distribution is proportional to the edge weights 
    (all the same if not given); then, assuming the relationship
    
    .. math::
        R^2 \propto t^{2/d_{rw}}
    
    where :math:`R^2` is the mean Euclidean square distance travelled, and
    :math:`t` the numnber of time steps, the random walk dimension (also 
    called random walk exponent or fractral dimension of the random walk),
    :math:`d_{rw}` is obtained by fitting a line on the log-log plot of 
    :math:`R^2` as function of :math:`t`.
    
    Parameters
    ----------
    G : networkx.Graph
        graph
 
    nrw : int
        number of random walks
        
    nsteps : int, optional
        number of time steps for every random walk;
        by default (`None`): `nsteps` is set to the outbox diagonal 
        (including the graph) divided by mean edge length
            
    nsteps_max : int, default: 100000
        maximal number of time steps `nsteps` is set to `nsteps_max`if
        `nsteps` was greater

    edge_weight : str, optional
        name of the edge attribute used as weight (probability to choose 
        an edge is proportional to its weight);
        by default (`None`) : all edges have a weight of 1

    pos_attr : str, default: 'pos'
        name of the attribute attached to nodes defining the position as a 
        sequence of 2 or 3 floats (used to compute distances)
    
    u0 : node key, optional
        starting node in `G` for all random walks;
        by default (`None`) a new node is randomly selected for every random walk

    data_points_removed_first_fraction : float, default : 0.33
        see `data_points_removed_last_fraction`

    data_points_removed_last_fraction : float, default : 0.0
        `data_points_removed_first_fraction` and `data_points_removed_last_fraction` are
        positive numbers in (0, 1); sequences of time steps and mean square distances
        constitute the data points; 
        the fraction `data_points_removed_first_fraction`, resp.
        `data_points_removed_last_fraction`, at the beginning, resp. end, of the
        data points are removed before doing the fit (line on the log-log plot);
        this allows to deal with the "extreme" data points

    confidence_level : float, default: 0.95
        confidence level, float in the interval (0, 1), to compute the 
        uncertainty for the line fitting (and for the random walk dimension)

    return_poly_fit_params : bool, default: True
        if `True`, the array (shape (2, )) of parameters of the fitted 
        polynom (line) is returned

    return_poly_fit_params_cov : bool, default: True
        if `True`, the array (shape (2, 2)) of parameters covariance of the 
        fitted polynom (line) is returned
    
    return_index_used_for_fit: bool, default: True
        if `True`, the starting (included) and ending (excluded) indices of the sequence 
        of data points used for fitting is returned
    
    return_rw_square_dist_mean : bool, default: True
        if `True`, the array of mean square distance is returned
        (mean over all random walks, for each time step)

    return_rw_square_dist_full : bool, default: False
        if `True`, the full array of square distances is returned
        (array of shape `(nrw, nsteps)`)

    seed : int, optional
        seed number for initializing the random number generator

    verbose : int, default: 1
        verbose mode, larger value implies more info printed

    make_plot : bool, default: True
        indicates if the plot of the data points and fitting curve is displayed 
        (in the current figure axis)

    plot_in_log_scale : bool, default: True
        used if `make_plot=True`, indicates if the plot is in log-log scale
        (along x- and y-axes)

    Returns
    -------
    out_dict : dict
        dictionary with the following keys / values:
        
        - drw : float
            random walk dimension
    
        - drw_lim_min : float
            lower limit for random walk dimension, provided from the 
            uncertainty on the slope of the fitted line (where a confidence 
            interval at a confidence level `confidence_level` derived from 
            Gaussian distribution is used)
    
        - drw_lim_max : float
            upper limit for random walk dimension, provided from the 
            uncertainty on the slope of the fitted line (where a confidence 
            interval at a confidence level `confidence_level` derived from 
            Gaussian distribution is used)

        - confidence_level : float
            confidence level used to compute the uncertainty for the line fitting 
            
        - poly_fit_params : array of shape (2, ), optional
            returned if `return_poly_fit_params=True`:
            array of parameters of the fitted polynom (line); in particular:
            `drw = 2.0 / poly_fit_params[0]`

        - poly_fit_params_cov : array of shape (2, 2), optional
            returned if `return_poly_fit_params_cov=True`:
            array (shape (2, 2)) of parameters covariance of the fitted polynom 
            (line); in particular: with
            `delta = scipy.stats.norm.ppf((1+confidence_level)/2) * np.sqrt(poly_fit_params_cov[0, 0])`
            `2.0 / drw +/- delta` gives the interval around the slope of the
            fitted line
        
        - index_used_for_fit : list of two ints, optional
            returned if `return_index_used_for_fit=True`:
            starting (included) index and ending (excluded) index of the sequence of
            data points used for fitting
        
        - rw_square_dist_mean : array of shape (nsteps,), optional
            returned if `return_rw_square_dist_mean=True`: 
            array of mean square distance over the random walks, 
            for each time steps from 1 to nsteps

        - rw_square_dist_full : array of shape (nrw, nsteps), optional
            returned if `return_rw_square_dist_full=True`: 
            the full array of square distances
    """
    if seed is not None:
        np.random.seed(seed)

    # Set nsteps
    if nsteps is None:
        if verbose > 0:
            print(f'Compute nsteps [outbox diagonal / mean edge length]')
        G_outbox_diagonal = kn.utils.outbox_diagonal(G, pos_attr=pos_attr)
        G_mean_edge_length = kn.utils.mean_edge_length(G, pos_attr=pos_attr)
        nsteps = int(G_outbox_diagonal/G_mean_edge_length)
    
    if nsteps > nsteps_max:
        if verbose > 0:
            print(f'WARNING `nsteps` ({nsteps}) too large, reduced to {nsteps_max}')
        nsteps = nsteps_max

    if verbose > 0:
        print(f'Do random walks ({nrw} random walks of {nsteps} time steps)')
    
    t1_all = time.time()
    out_rw = graph_random_walk_square_distance(
            G, nrw, nsteps,
            edge_weight=edge_weight, 
            pos_attr=pos_attr, 
            u0=u0, 
            seed=None, 
            return_full_array=return_rw_square_dist_full,
            verbose=verbose-1)

    if return_rw_square_dist_full:
        rw_square_dist_mean = out_rw[0]
        rw_square_dist_full = out_rw[1]
    else:
        rw_square_dist_mean = out_rw

    t2_all = time.time()
    if verbose > 0:
        print(f'Total elapsed time (random walks): {t2_all-t1_all:.2f} sec')

    if verbose > 0:
        print(f'Compute random walk dimension (line fitting on log-log plot of mean square distance as function of time step)')

    # Compute random walk dimension (line fitting on log-log plot of mean square distance as function of time step)
    # -------------------------------------------------------------------------------------------------------------
    # Data (all)
    x_all = np.arange(1, nsteps+1)
    y_all = rw_square_dist_mean

    poly_fit_params, poly_fit_params_cov, x, y, n1, n2 = _loglog_fit(x_all, y_all, data_points_removed_first_fraction, data_points_removed_last_fraction, verbose=verbose)

    # Get conductance exponent dOmega, and dOmega_delta
    if poly_fit_params is not None:
        slope = poly_fit_params[0]
        drw = 2.0 / slope
        fit_ok = True
    else:
        drw = np.nan
        fit_ok = False
    
    if poly_fit_params_cov is not None:
        # slope is defined (poly_fit_params is not None)
        t = scipy.stats.norm.ppf((1+confidence_level)/2)
        slope_delta = t * np.sqrt(poly_fit_params_cov[0, 0])
        drw_lim_min = 2.0 / (slope + slope_delta)
        drw_lim_max = 2.0 / (slope - slope_delta)
    else:
        drw_lim_min = np.nan
        drw_lim_max = np.nan

    # Make plot (if needed)
    if make_plot:
        plot_fit_random_walk_dimension(
            x_all, y_all, x, y, drw, drw_lim_min, drw_lim_max, confidence_level, poly_fit_params, 
            plot_fit=fit_ok, plot_in_log_scale=plot_in_log_scale)

    # Set output dictionary
    out_dict = {'drw': drw, 'drw_lim_min': drw_lim_min, 'drw_lim_max': drw_lim_max, 'confidence_level': confidence_level}

    if return_poly_fit_params:
        out_dict['poly_fit_params'] = poly_fit_params

    if return_poly_fit_params_cov:
        out_dict['poly_fit_params_cov'] = poly_fit_params_cov

    if return_index_used_for_fit:
        out_dict['index_used_for_fit'] = [n1, n2]

    if return_rw_square_dist_mean:
        out_dict['rw_square_dist_mean'] = rw_square_dist_mean

    if return_rw_square_dist_full:
        out_dict['rw_square_dist_full'] = rw_square_dist_full

    return out_dict
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
def graph_random_walk_dimension_update(
        out_dict,
        data_points_removed_first_fraction=0.33,
        data_points_removed_last_fraction=0.0,
        confidence_level=0.95,
        verbose=1,
        make_plot=True,
        plot_in_log_scale=True):
    """
    Updates the fit (power law) for the random walk dimension.

    It is assumed that the function `graph_random_walk_dimension` has 
    already been run, and that the output dictionary `out_dict` contains the
    key `rw_square_dist_mean`; 
    the fit is updated by changing the fractions of data points removed at the 
    beginning and end of the sequences of data points.

    See function `graph_random_walk_dimension` for details.

    Note: this function operates inplace on the `out_dict` dictionary (it is
    updated and returned).

    Parameters
    ----------
    out_dict : dict
        output dictionary of the function `graph_random_walk_dimension`;
        the dictionary should contain the key `rw_square_dist_mean`

    data_points_removed_first_fraction : float, default : 0.33
        see `data_points_removed_last_fraction`

    data_points_removed_last_fraction : float, default : 0.0
        `data_points_removed_first_fraction` and `data_points_removed_last_fraction` are
        positive numbers in (0, 1); sequences of time steps and mean square distances
        constitute the data points; 
        the fraction `data_points_removed_first_fraction`, resp.
        `data_points_removed_last_fraction`, at the beginning, resp. end, of the
        data points are removed before doing the fit (line on the log-log plot);
        this allows to deal with the "extreme" data points

    confidence_level : float, default: 0.95
        confidence level, float in the interval (0, 1), to compute the 
        uncertainty for the line fitting (and for the random walk dimension)

    verbose : int, default: 1
        verbose mode, larger value implies more info printed

    make_plot : bool, default: True
        indicates if the plot of the data points and fitting curve is displayed 
        (in the current figure axis)

    plot_in_log_scale : bool, default: True
        used if `make_plot=True`, indicates if the plot is in log-log scale
        (along x- and y-axes)

    Returns
    -------
    out_dict : dict
        updated output dictionary
    """
    if 'rw_square_dist_mean' not in out_dict.keys():
        raise ValueError('The input dictionary `out_dict` should contain the key `rw_square_dist_mean`')

    if verbose > 0:
        print(f'Compute random walk dimension (line fitting on log-log plot of mean square distance as function of time step)')

    # Compute random walk dimension (line fitting on log-log plot of mean square distance as function of time step)
    # -------------------------------------------------------------------------------------------------------------
    # Data (all)
    y_all = out_dict['rw_square_dist_mean']
    x_all = np.arange(1, len(y_all)+1)

    poly_fit_params, poly_fit_params_cov, x, y, n1, n2 = _loglog_fit(x_all, y_all, data_points_removed_first_fraction, data_points_removed_last_fraction, verbose=verbose)

    # Get conductance exponent dOmega, and dOmega_delta
    if poly_fit_params is not None:
        slope = poly_fit_params[0]
        drw = 2.0 / slope
        fit_ok = True
    else:
        drw = np.nan
        fit_ok = False
    
    if poly_fit_params_cov is not None:
        # slope is defined (poly_fit_params is not None)
        t = scipy.stats.norm.ppf((1+confidence_level)/2)
        slope_delta = t * np.sqrt(poly_fit_params_cov[0, 0])
        drw_lim_min = 2.0 / (slope + slope_delta)
        drw_lim_max = 2.0 / (slope - slope_delta)
    else:
        drw_lim_min = np.nan
        drw_lim_max = np.nan

    # Make plot (if needed)
    if make_plot:
        plot_fit_random_walk_dimension(
            x_all, y_all, x, y, drw, drw_lim_min, drw_lim_max, confidence_level, poly_fit_params, 
            plot_fit=fit_ok, plot_in_log_scale=plot_in_log_scale)

    # Update output dictionary
    out_dict['drw'] = drw
    out_dict['drw_lim_min'] = drw_lim_min
    out_dict['drw_lim_max'] = drw_lim_max
    out_dict['confidence_level'] = confidence_level

    if 'poly_fit_params' in out_dict.keys():
        out_dict['poly_fit_params'] = poly_fit_params

    if 'poly_fit_params_cov' in out_dict.keys():
        out_dict['poly_fit_params_cov'] = poly_fit_params_cov

    if 'index_used_for_fit' in out_dict.keys():
        out_dict['index_used_for_fit'] = [n1, n2]

    return out_dict
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
def plot_fit_random_walk_dimension(
        x_all,
        y_all,
        x,
        y,
        drw,
        drw_lim_min,
        drw_lim_max,
        confidence_level,
        poly_fit_params,
        plot_fit=True,
        plot_in_log_scale=True):
    """
    Makes the plot (in current axis figure) of the data points and fitted curve for random walk dimension (power law).

    Parameters
    ----------
    x_all : array of floats
        x-coordinates of all data points

    y_all : array of floats
        y-coordinates of all data points

    x : array of floats
        x-coordinates of data points used for fitting (a subset of `x_all`)
    
    y : array of floats
        y-coordinates of data points used for fitting (a subset of `y_all`)

    drw : float
        random walk dimension
    
    drw_lim_min : float
        lower limit for random walk dimension, provided from the 
        uncertainty on the slope of the fitted line (where a confidence 
        interval at a confidence level `confidence_level` derived from 
        Gaussian distribution is used)

    drw_lim_max : float
        upper limit for random walk dimension, provided from the 
        uncertainty on the slope of the fitted line (where a confidence 
        interval at a confidence level `confidence_level` derived from 
        Gaussian distribution is used)

    confidence_level : float
        confidence level, float in the interval (0, 1), used to compute the 
        uncertainty for the line fitting (and for the fractal dimension)

    poly_fit_params : array of shape (2, ) or None
        array of parameters of the fitted polynom (line); in particular:
        `df = -poly_fit_params[0]`; if no fit was done, `poly_fit_params` is `None`

    plot_fit : bool ; default : True
        indicates if the plot of the fitted curve is done (in the current figure axis);
        if `False`, only the data points are plotted

    plot_in_log_scale : bool, default: True
        indicates if the plot is in log-log scale (along x- and y-axes)

    Returns
    -------
    None
    """
    # Estimation with the fitted line
    if plot_fit:
        xx = np.linspace(x.min(), x.max(), 200)
        yy = np.power(10, poly_fit_params[0] * np.log10(xx) + poly_fit_params[1])

        xx_all = np.linspace(x_all.min(), x_all.max(), 200)
        yy_all = np.power(10, poly_fit_params[0] * np.log10(xx_all) + poly_fit_params[1])

    # Plot parameters
    xlabel = '$t$ : time step'
    ylabel = '$R^2$ : mean square distance'

    title = '$R^2\propto t^{2/d_{rw}}$' + ', $d_{rw}$' + f'={drw:.5f} in [{drw_lim_min:.5f}, {drw_lim_max:.5f}] [{100*confidence_level:.2g}% - inter.]' 

    color_data_all      = 'tab:blue'
    marker_data_all     = '.'
    markersize_data_all = 5

    color_data_used      = 'tab:orange'
    marker_data_used     = '.'
    markersize_data_used = 5

    color_fit_data_used = 'red'
    ls_fit_data_used    = 'solid'
    lw_fit_data_used    = 1.

    color_fit_data_all = 'red'
    ls_fit_data_all    = 'dotted'
    lw_fit_data_all    = 1.

    # Figure
    # ------
    # plt.figure(figsize=(10, 6))
    plt.plot(x_all, y_all, ls='', marker=marker_data_all , markersize=markersize_data_all , color=color_data_all , label='data')
    if len(x):
        plt.plot(x    , y    , ls='', marker=marker_data_used, markersize=markersize_data_used, color=color_data_used, label='data used for fitting')
    if plot_fit:
        plt.plot(xx_all, yy_all, ls=ls_fit_data_all , lw=lw_fit_data_all , color=color_fit_data_all)
        plt.plot(xx    , yy    , ls=ls_fit_data_used, lw=lw_fit_data_used, color=color_fit_data_used, label='fit')
    
    if plot_in_log_scale:
        plt.xscale('log')
        plt.yscale('log')
    
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid()
    plt.legend()
    # plt.show()

    return
# ----------------------------------------------------------------------------

# ============================================================================
# Tools for graph related to laplacian
# ============================================================================

# ----------------------------------------------------------------------------
def laplacian_eigs(
        G, 
        edge_weight=None, 
        k=None,
        normalized_laplacian=False,
        **kwargs):
    """
    Computes the first eigen values and vectors of the Laplacian matrix.

    Let :math:`L` be the (weighted) Laplacian matrix of the graph, i.e.
    :math:`L = D - A` is a :math:`n \times n` symmetric matrix, where 
    :math:`n` is the number of graph nodes, :math:`D` is the diagonal matrix 
    with :math:`D(i)` equal to the  degree of the i-th node (or the sum of the 
    weight of the edges with an extremity at the i-th node) and :math:`A` the 
    adjacency matrix with :math:`A(i,j)` equal -1 (or :math:`-w(e)`) if the an
    edge :math:`e` exists between the i-th and j-th nodes.

    For a connected graph (i.e. with one connected component), 
    :math:`\operatorname{Ker}(L) = <(1, 1, ..., 1)^\\top>`, and the eigen values
    verify :math:`0 = \\lambda_1 < \\lambda_2 \leq \\ldots \\lambda_n`.

    This function computes the first `k` eigen values and eigen vectors, 
    see parameter `k` below.

    Parameters
    ----------
    G : networkx.Graph
        graph

    k : int, optional
        number (max) of eigen values and eigen vectors to compute, 
        the first `k` eigen values in ascending order are considered;
        by default (`None`): all eigen values and eigen vectors are 
        computed, i.e. `k` is set to `G.number_of_nodes()`

    edge_weight : str, optional
        name of the edge attribute used as weight
        by default (`None`) : all edges have a weight of 1
    
    normalized_laplacian : bool, default: False
        if `True`, the first eigenvalues and eigenvectors of the 
        normalized Laplacian matrix, i.e. 
        :math:`D^{-1/2}\\cdot L \\cdot D^{-1/2}`, where :math:`D` is 
        the (weighted) degree diagonal matrix of the graph;
        note that the eigenvalues of the normalized Laplacian matrix
        are all in the interval [0, 2]

    kwargs : dict
        keyword arguments passed to the function 
        `scipy.sparse.linalg.eigsh` to compute the `k` first eigen values
        and eigen vectors (unused if `k >= G.number_of_nodes()`, in this 
        case, the function `scipy.linalg.eigh` is used); 
        note: the parameter `which`, if not given in `kwargs`, is set to 
        'SM' (smallest eigen values are considered)

    Returns
    -------
    eigval : 1d-array of shape (k, )
        array of the first k eighen values of the Laplacian matrix, in 
        ascending order; 
        note: for a connected graph, 
        `0 = eigval[0] < eigenval[1] <= eigenval[2] ...`

    eigvec : 2d-array of shape (n, k)
        array of the first k eighen vectors (in columns) of the Laplacian matrix, 
        (n is the number of nodes in the graph), the eigen vectors are 
        normalized (norm of 1);
        note: for a connected graph, `eigvec[:,0]` is proportional to the 
        vector `[1, 1, ..., 1]`
    """
    if normalized_laplacian:
        # L: normalized Laplacian matrix (sparse matrix in CSR format)
        L = nx.normalized_laplacian_matrix(G, weight=edge_weight)
    else:
        # L: Laplacian matrix (sparse matrix in CSR format)
        L = nx.laplacian_matrix(G, weight=edge_weight)

    if k is None or k >= G.number_of_nodes():
        # Compute all eigen values and eigen vectors
        eigval, eigvec = np.linalg.eigh(L.toarray())
        
        # # Equivalent:
        # eigval, eigvec = scipy.linalg.eigh(L.toarray())

        # # Equivalent (slower, symmetry of L is not assumed), and the eigen values are not sorted
        # eigval, eigvec = np.linalg.eig(L.toarray())
        # # or :
        # eigval, eigvec = scipy.linalg.eig(L.toarray())
        # # Sort eigen values in ascending order and set eigen vectors accordingly
        # ind = np.argsort(eigval)
        # eigval = eigval[ind]
        # eigvec = eigvec[:, ind]
    else:
        if 'which' not in kwargs.keys():
            kwargs['which'] = 'SM' # smallest eigen values
        
        eigval, eigvec = scipy.sparse.linalg.eigsh(L.astype('float'), k=k, **kwargs)

    return eigval, eigvec
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
def laplacian_pseudo_inverse(
        G, 
        edge_weight=None, 
        k=None, 
        eigval=None, 
        eigvec=None, 
        **kwargs):
    """
    Computes the pseudo inverse of the Laplacian matrix, for a connected graph.

    Let :math:`L` be the (weighted) Laplacian matrix of the graph, i.e.
    :math:`L = D - A` is a :math:`n \times n` symmetric matrix, where 
    :math:`n` is the number of graph nodes, :math:`D` is the diagonal matrix 
    with :math:`D(i)` equal to the  degree of the i-th node (or the sum of the 
    weight of the edges with an extremity at the i-th node) and :math:`A` the 
    adjacency matrix with :math:`A(i,j)` equal -1 (or :math:`-w(e)`) if the an
    edge :math:`e` exists between the i-th and j-th nodes.

    For a connected graph (i.e. with one connected component), 
    :math:`\operatorname{Ker}(L) = <(1, 1, ..., 1)^\\top>`, and the eigen values
    verify :math:`0 = \\lambda_1 < \\lambda_2 \leq \\ldots \\lambda_n`.

    Let :math:`L^{+}` be the pseudo inverse of the Laplacian matrix, then
    
    .. math::
        L^{+}=\\sum_{i} \\frac{1}{\\lambda_i} u_i \\cdot u_i^\\top,
    
    where :math:`u_i` is a eigen vector of norm 1 for the eigen value
    :math:`\\lambda_i`, and where the sum goes through the non-zero eigen 
    values.

    Furthermore, 
    
    .. math::
        L^{+} \\cdot L = L \\cdot L^{+} = \\sum_{i} u_i \\cdot u_i^\\top,

    with the same sum indices as above, is the orthogonal projection onto 
    :math:`\operatorname{Im}(L) = \operatorname{Ker}(L)^\\bot`.

    This function computes the approximation of the pseudo inverse of 
    the Laplacian matrix by accounting for the first `k` eigen values and 
    eigen vectors (of the Laplacian matrix), see parameter `k` below.

    Parameters
    ----------
    G : networkx.Graph
        graph, must be connected (i.e. with one connected component)

    edge_weight : str, optional
        name of the edge attribute used as weight; 
        unused if `eigval` and `eigvec` are specified (not `None`);
        by default (`None`) : all edges have a weight of 1

    k : int, optional
        number (max) of eigen values and eigen vectors (of the Laplacian
        matrix) to compute, the first `k` eigen values in ascending order 
        are considered;
        unused if `eigval` and `eigvec` are specified (not `None`);
        by default (`None`): all eigen values and eigen vectors are 
        computed, i.e. `k` is set to `G.number_of_nodes()`

    eigval : 1d-array, optional
        array of shape (m, ), with m less than or equal to the number 
        of nodes in the graph, m first eigen values of the Laplacian 
        matrix, in ascending order, i.e. for a connected graph, 
        `0 = eigval[0] < eigenval[1] <= eigenval[2] ...`;
        by default (`None`): the first `k` eigen values are computed

    eigvec : 2d-array, optional
        array of shape (n, m), where n is the number of nodes in the
        graph and m <= n, m first eigen vectors in columns (orthogonal 
        and of norm 1) corresponding to the eigen values `eigval`; 
        for a connected graph, `eigvec[:,0]` is proportional to the 
        vector `[1, 1, ..., 1]`;
        by default (`None`): the first `k` eigen vectors are computed

    kwargs : dict
        keyword arguments passed to the function 
        `scipy.sparse.linalg.eigsh` to compute the `k` first eigen values
        and eigen vectors (unused if `k >= G.number_of_nodes()`, in this 
        case, the function `scipy.linalg.eigh` is used); 
        note: the parameter `which`, if not given in `kwargs`, is set to 
        'SM' (smallest eigen values are considered)

    Returns
    -------
    Lpinv : 2d-array of shape (n, n)
        array of shape (n, n), (approximation of) pseudo inverse of the 
        Laplacian matrix (see above)
    """
    if eigvec is None or eigval is None:
        eigval, eigvec = laplacian_eigs(G, edge_weight=edge_weight, k=k, **kwargs)

    # Psendo inverse (approximation) of the Laplacian matrix
    Lpinv = (eigvec[:,1:] * 1.0/eigval[1:]) @ eigvec[:,1:].T
    # Lpinv = eigvec[:,1:] @ np.diag(1.0/eigval[1:]) @ eigvec[:,1:].T # equiv
    # Lpinv = np.sum([1/eigval[i] * np.outer(eigvec[:, i], eigvec[:, i]) for i in range(1, len(eigval))], axis=0) # equiv

    # # Projection onto (ker L)-orthogonal
    # G_proj = eigvec[:,1:] @ eigvec[:,1:].T
    # # G_proj = np.sum([np.outer(eigvec[:, i], eigvec[:, i]) for i in range(1, len(eigval))], axis=0) # equiv.
    # np.allclose(Lpinv @ L, G_proj) # should be True

    return Lpinv
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
def solve_flow_without_bc(
        G,
        q_nodes_dict,
        edge_weight=None,
        k=None, 
        eigval=None, 
        eigvec=None, 
        **kwargs):
    """
    Solves the flow problem in a connected graph.

    Computes the potential at each node in the graph, given the source term 
    (entering fluxes) at the graph nodes.

    The source term (entering fluxes) at the graph nodes is given by `q_nodes_dict`.

    The edge weights represent the conductances of the edges; if not specified
    (`edge_weight=None`), all edges have a conductance of 1.

    Denoting :math:`L` is the Laplacian matrix of the graph (with given edge 
    weights), the potential at the inner nodes are given by 
    
    Let :math:`L^{+}` be the pseudo inverse of the Laplacian matrix, i.e.
    
    .. math::
        L^{+}=\\sum_{i} \\frac{1}{\\lambda_i} u_i \\cdot u_i^\\top,
    
    where :math:`u_i` is a eigen vector of norm 1 for the eigen value
    :math:`\\lambda_i`, and where the sum goes through the non-zero eigen 
    values.

    The solution (potential at the graph nodes) is given by

    .. math::
        p = L^{+} \\cdot q,

    where :math:`q` is the source term (vector of entering fluxes at the graph nodes).

    This function uses the approximation of the pseudo inverse of 
    the Laplacian matrix by accounting for the first `k` eigen values and 
    eigen vectors (of the Laplacian matrix), see parameter `k` below.

    Note: it is recommended to pre-compute (first) eigen values and eigen vectors
    (and specify parameters `eigval` and `eigvec`).

    Parameters
    ----------
    G : networkx.Graph
        graph, must be connected (i.e. with one connected component)

    q_nodes_dict : dict
        dictionary of entering flux at the graph nodes, the keys are the node
        labels and the values are the entering flux values
            
    edge_weight : str, optional
        name of the edge attribute used as weight (conductance)

    k : int, optional
        number (max) of eigen values and eigen vectors (of the Laplacian
        matrix) to compute, the first `k` eigen values in ascending order 
        are considered;
        unused if `eigval` and `eigvec` are specified (not `None`);
        by default (`None`): all eigen values and eigen vectors are 
        computed, i.e. `k` is set to `G.number_of_nodes()`

    eigval : 1d-array, optional
        array of shape (m, ), with m less than or equal to the number 
        of nodes in the graph, m first eigen values of the Laplacian 
        matrix, in ascending order, i.e. for a connected graph, 
        `0 = eigval[0] < eigenval[1] <= eigenval[2] ...`;
        by default (`None`): the first `k` eigen values are computed

    eigvec : 2d-array, optional
        array of shape (n, m), where n is the number of nodes in the
        graph and m <= n, m first eigen vectors in columns (orthogonal 
        and of norm 1) corresponding to the eigen values `eigval`; 
        for a connected graph, `eigvec[:,0]` is proportional to the 
        vector `[1, 1, ..., 1]`;
        by default (`None`): the first `k` eigen vectors are computed

    kwargs : dict
        keyword arguments passed to the function 
        `scipy.sparse.linalg.eigsh` to compute the `k` first eigen values
        and eigen vectors (unused if `k >= G.number_of_nodes()`, in this 
        case, the function `scipy.linalg.eigh` is used); 
        note: the parameter `which`, if not given in `kwargs`, is set to 
        'SM' (smallest eigen values are considered) 

    Returns
    -------
    p_nodes_dict : dict
        dictionary of potential at all nodes in the graph, the keys are the node
        labels and the values are the potential values
    """
    if eigvec is None or eigval is None:
        eigval, eigvec = laplacian_eigs(G, edge_weight=edge_weight, k=k, **kwargs)

    # Set q vector (source term at the graph nodes, sorted by node labels)
    q = np.asarray([q_nodes_dict[k] for k in G.nodes()])

    # Solve the system
    p = (eigvec[:, 1:] * 1.0/eigval[1:]) @ eigvec[:, 1:].T @ q
    # # or:
    # p = np.sum([1/val * vec.dot(q) * vec for val, vec in zip(eigval[1:], eigvec[:, 1:].T)], axis=0)
    
    p_nodes_dict = {k: v for k, v in zip(G.nodes(), p)}

    return p_nodes_dict
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
def solve_flow_with_bc_inner_connected(
        G,
        pB_nodes_dict,
        qI_nodes_dict,
        edge_weight=None,
        return_entering_fluxes=True):
    """
    Solves the flow problem in a connected graph, with boundary condition (potential).

    Computes the potential at each node in the graph, given the boundary
    potential and the source term (entering fluxes) at the inner nodes.

    The potential at the boundary nodes is given by `pB_nodes_dict`, and the
    source term (entering fluxes) at the inner nodes is given by `qI_nodes_dict`.
    The boundary nodes B (keys of `pB_nodes_dict`) and the inner nodes I (keys of 
    `qI_nodes_dict`) must form a partition of all graph nodes (i.e. disjoint sets 
    with union equal to the set of all graph nodes). Moreover, the boundary nodes B 
    must be a non empty set, and the subgraph induced by the inner nodes I must be 
    connected (unless I is empty).

    The edge weights represent the conductances of the edges; if not specified
    (`edge_weight=None`), all edges have a conductance of 1.

    Denoting :math:`L` is the Laplacian matrix of the graph (with given edge 
    weights), the potential at the inner nodes are given by 
    
    .. math::
        p_I = {L_{II}}^{-1}\\cdot (q_I - L_{IB} p_B)
    
    where the submatrix :math:`L_{II}` is invertible provided that the subgraph
    induced by the nodes I is connected, and the entering fluxes at the boundary 
    nodes are given by 

    .. math::
        q_B = L_{BI}\\cdot p_I + L_{BB} p_B
    
    Parameters
    ----------
    G : networkx.Graph
        graph, must be connected (i.e. with one connected component)

    pB_nodes_dict : dict
        dictionary of potential at boundary nodes, the keys are the node
        labels and the values are the potential values

    qI_nodes_dict : dict
        dictionary of entering flux at inner nodes, the keys are the node
        labels and the values are the entering flux values
            
    edge_weight : str, optional
        name of the edge attribute used as weight (conductance)
    
    return_entering_fluxes : bool, default: True
        if `True`, the entering flux at all nodes in the graph is returned

    Returns
    -------
    p_nodes_dict : dict
        dictionary of potential at all nodes in the graph, the keys are the node
        labels and the values are the potential values
    
    q_nodes_dict : dict, optional
        if `return_entering_fluxes=True`, dictionary of entering flux at all nodes 
        in the graph, the keys are the node labels and the values are the entering 
        flux values
    """
    # Set dictionary to convert node label to node index
    node_label2index = {u:i for i, u in enumerate(G.nodes())}
    # node_index2label = kn.utils.get_node_label2index(G) # equivalent

    B_nodes_label_list = pB_nodes_dict.keys()
    I_nodes_label_list = qI_nodes_dict.keys()

    B_nodes_index_list = [node_label2index[k] for k in B_nodes_label_list]
    I_nodes_index_list = [node_label2index[k] for k in I_nodes_label_list]

    if len(B_nodes_index_list) == 0:
        raise ValueError('Boundary nodes must be a non empty set.')
    
    if len(np.intersect1d(B_nodes_index_list, I_nodes_index_list)) > 0:
        raise ValueError('Boundary nodes and internal nodes must be disjoint sets.')

    if len(B_nodes_index_list) + len(I_nodes_index_list) != G.number_of_nodes():
        raise ValueError('Boundary nodes and internal nodes must cover all nodes in the graph.')

    if len(B_nodes_index_list) == G.number_of_nodes():
        # All nodes are boundary nodes, no inner nodes
        # "Sort" the dictionary to have the same order as the nodes (labels) in the graph

        # Potential at all nodes
        p_nodes_dict = {k: pB_nodes_dict[k] for k in G.nodes()}        

        if return_entering_fluxes:
            p = np.array(list(pB_nodes_dict.values()))
        
            # Compute entering flux at all nodes
            L = nx.laplacian_matrix(G, weight=edge_weight) # sparse matrix in CSR format
            q = L @ p
            q_nodes_dict = {k: q[k] for k in G.nodes()}

            return p_nodes_dict, q_nodes_dict

        return p_nodes_dict

    if not nx.is_connected(G.subgraph(I_nodes_label_list)):
        raise ValueError('Removing boundary nodes disconnect the graph, flow is not solved.')

    # Solve the flow problem: get pI (potential at inner nodes)
    pB = np.array(list(pB_nodes_dict.values()))
    qI = np.array(list(qI_nodes_dict.values()))

    L = nx.laplacian_matrix(G, weight=edge_weight) # sparse matrix in CSR format

    L_I = L[I_nodes_index_list, :]
    L_II = L_I[:, I_nodes_index_list]
    L_IB = L_I[:, B_nodes_index_list]

    pI = scipy.sparse.linalg.spsolve(L_II, qI - L_IB @ pB)

    # Set dictionary of potential at inner nodes
    pI_nodes_dict = {k: v for k, v in zip(I_nodes_label_list, pI)}

    # Concatenate dictionaries pI_nodes_dict and pB_nodes_dict
    p_nodes_dict = pI_nodes_dict.copy()
    p_nodes_dict.update(pB_nodes_dict)  
    # dict(**pB_nodes_dict, **pI_nodes_dict) # concatenate dictionaries (works if keys are strings)

    # "Sort" the dictionary to have the same order as the nodes (labels) in the graph
    p_nodes_dict = {k: p_nodes_dict[k] for k in G.nodes()}

    if return_entering_fluxes:
        # Solve the flow problem: get qB (entering flux at boundary nodes)
        L_BB = L[B_nodes_index_list, :][:, B_nodes_index_list]

        qB = L_IB.transpose() @ pI + L_BB @ pB

        # Set dictionary of entering flux at boundary nodes
        qB_nodes_dict = {k: v for k, v in zip(B_nodes_label_list, qB)}


        # Concatenate dictionaries qI_nodes_dict and qB_nodes_dict
        q_nodes_dict = qI_nodes_dict.copy()
        q_nodes_dict.update(qB_nodes_dict)  

        # "Sort" the dictionary to have the same order as the nodes (labels) in the graph
        q_nodes_dict = {k: q_nodes_dict[k] for k in G.nodes()}

        return p_nodes_dict, q_nodes_dict

    return p_nodes_dict
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
def solve_flow_with_bc(
        G,
        pB_nodes_dict,
        qI_nodes_dict,
        edge_weight=None,
        return_entering_fluxes=True):
    """
    Solves the flow problem in a connected graph, with boundary condition (potential).

    Computes the potential at each node in the graph, given the boundary
    potential and the source term (entering fluxes) at the inner nodes.

    The potential at the boundary nodes is given by `pB_nodes_dict`, and the
    source term (entering fluxes) at the inner nodes is given by `qI_nodes_dict`.
    The boundary nodes B (keys of `pB_nodes_dict`) and the inner nodes I (keys of 
    `qI_nodes_dict`) must form a partition of all graph nodes (i.e. disjoint sets 
    with union equal to the set of all graph nodes). Moreover, the boundary nodes B 
    must be a non empty set.
    
    Let G(I) the subgraph induced by the set of inner nodes I, and let
    I_1, ..., I_m the partition of the set I, such that G(I_1), ..., G(I_m) are the 
    connected components of G(I). For each subset I_k of I, let G_k be the connected 
    component of the subgraph induced by the union of I_k and B that contains the 
    nodes I_k. Let B_k be the nodes of B in G_k (some nodes in B may have been removed
    to have a connected graph, that is B_k may be not equal to B). Then, G_k = G(I_k, B_k) 
    is the subgraph induced by the inner nodes I_k and the boundary nodes B_k. 
    The flow problem is solved for each G_k, using the function 
    `solve_flow_with_bc_inner_connected`, providing the potentials p_I_k at the nodes I_k
    and (if desired) the fluxes q_B_k entering at the nodes B_k. The solution of the flow
    problem in the entire graph is then obtained by gathering the vectors p_I_k to get 
    the potentials p_I at the nodes I, and (if desired) by gathering the vectors q_B_k
    and summing the values at the common nodes, to get the flux q_B entering at the nodes B.

    The edge weights represent the conductances of the edges; if not specified
    (`edge_weight=None`), all edges have a conductance of 1.

    Parameters
    ----------
    G : networkx.Graph
        graph, must be connected (i.e. with one connected component)

    pB_nodes_dict : dict
        dictionary of potential at boundary nodes, the keys are the node
        labels and the values are the potential values

    qI_nodes_dict : dict
        dictionary of entering flux at inner nodes, the keys are the node
        labels and the values are the entering flux values
            
    edge_weight : str, optional
        name of the edge attribute used as weight (conductance)
    
    return_entering_fluxes : bool, default: True
        if `True`, the entering flux at all nodes in the graph is returned

    Returns
    -------
    p_nodes_dict : dict
        dictionary of potential at all nodes in the graph, the keys are the node
        labels and the values are the potential values
    
    q_nodes_dict : dict, optional
        if `return_entering_fluxes=True`, dictionary of entering flux at all nodes 
        in the graph, the keys are the node labels and the values are the entering 
        flux values
    """
    # Set dictionary to convert node label to node index
    node_label2index = {u:i for i, u in enumerate(G.nodes())}
    # node_index2label = kn.utils.get_node_label2index(G) # equivalent

    B_nodes_label_list = pB_nodes_dict.keys()
    I_nodes_label_list = qI_nodes_dict.keys()

    B_nodes_index_list = [node_label2index[k] for k in B_nodes_label_list]
    I_nodes_index_list = [node_label2index[k] for k in I_nodes_label_list]

    if len(B_nodes_index_list) == 0:
        raise ValueError('Boundary nodes must be a non empty set.')
    
    if len(np.intersect1d(B_nodes_index_list, I_nodes_index_list)) > 0:
        raise ValueError('Boundary nodes and internal nodes must be disjoint sets.')

    if len(B_nodes_index_list) + len(I_nodes_index_list) != G.number_of_nodes():
        raise ValueError('Boundary nodes and internal nodes must cover all nodes in the graph.')

    if len(B_nodes_index_list) == G.number_of_nodes():
        # All nodes are boundary nodes, no inner nodes
        # "Sort" the dictionary to have the same order as the nodes (labels) in the graph

        # Potential at all nodes
        p_nodes_dict = {k: pB_nodes_dict[k] for k in G.nodes()}        

        if return_entering_fluxes:
            p = np.array(list(pB_nodes_dict.values()))
        
            # Compute entering flux at all nodes
            L = nx.laplacian_matrix(G, weight=edge_weight) # sparse matrix in CSR format
            q = L @ p
            q_nodes_dict = {k: q[k] for k in G.nodes()}

            return p_nodes_dict, q_nodes_dict

        return p_nodes_dict

    # Sub-graph induced by inner nodes I
    G_sub = G.subgraph(I_nodes_label_list)
    
    if nx.is_connected(G_sub) == 1:
        # Solve flow on the entire graph
        return solve_flow_with_bc_inner_connected(
            G,
            pB_nodes_dict,
            qI_nodes_dict,
            edge_weight=edge_weight,
            return_entering_fluxes=return_entering_fluxes)

    # Initialize solution of node potentials
    p_nodes_dict = {k: 0.0 for k in G.nodes()}
    p_nodes_dict.update(pB_nodes_dict)

    if return_entering_fluxes:
        # Initialize solution of entering fluxes
        q_nodes_dict = {k: 0.0 for k in G.nodes()}
    
    for G_sub_nodes in nx.connected_components(G_sub):
        # -> G_sub_nodes : set of nodes I_i (forming a connected component of G_sub)
        # Get the entering fluxes at I_i
        qI_nodes_dict_i = {k: qI_nodes_dict[k] for k in G_sub_nodes}
        
        # Get one (first) node of I_i
        u0 = next(iter(G_sub_nodes))

        # Get the the subgraph G(I_i, B) induced by the union of I_i and B
        G_sub_nodes.update(set(B_nodes_label_list)) # add nodes from B
        G_i = G.subgraph(G_sub_nodes) # get the induced subgraph

        if nx.is_connected(G_i):
            pB_nodes_dict_i = pB_nodes_dict
        else:
            # if G_i is not connected, extract its connected component containing the nodes I_i 
            # (it suffices to get the component containing u0)
            for G_i_nodes in nx.connected_components(G_i):
                if u0 in G_i_nodes:
                    G_i = G_i.subgraph(G_i_nodes)
                    pB_nodes_dict_i = {k: v for k, v in pB_nodes_dict.items() if k in G_i_nodes}
                    break
        
        # Solve flow on G_i
        out = solve_flow_with_bc_inner_connected(
                G_i,
                pB_nodes_dict_i,
                qI_nodes_dict_i,
                edge_weight=edge_weight,
                return_entering_fluxes=return_entering_fluxes)
        
        if return_entering_fluxes:
            p_nodes_dict_i, q_nodes_dict_i = out
            
            # Update solution of entering fluxes
            for k, v in q_nodes_dict_i.items():
                q_nodes_dict[k] = q_nodes_dict[k] + v

        else:
            p_nodes_dict_i = out

        # Update solution of node potentials
        p_nodes_dict.update(p_nodes_dict_i)

    if return_entering_fluxes:
        return p_nodes_dict, q_nodes_dict

    return p_nodes_dict
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
def compute_abs_flow_along_edges(
        G,
        edge_weight=None,
        node_potential='p'):
    """
    Computes the absolute flows along edges in a graph.
    
    This function computes the absolute flows along each edge in a graph `G`, 
    given the edge conductances `w` (edge attribute) and the node potential `p` 
    (node attribute), based on the relation `q(u, v) = w(u, v) * (p(u) - p(v))`
    for an edge `(u, v)`

    Parameters
    ----------
    G : networkx.Graph
        graph
    
    edge_weight : str, optional
        name of the edge attribute used as weight (conductance); 
        by default (`None`) : all edges have a weight of 1

    node_potential : str, default: 'p'
        name of the node attribute used as potential

    Returns
    -------
    edge_abs_flow_dict : dict
        dictionary of absolute flows along edges, the keys are the edges and 
        the values are the absolute flows
    """
    # Compute absolute flows along edges
    if edge_weight is None:
        edge_abs_flow_dict = {(u, v): np.abs(G.nodes[u][node_potential] - G.nodes[v][node_potential]) for u, v in G.edges()}
    else:
        edge_abs_flow_dict = {(u, v): np.abs(G.edges[u, v][edge_weight] * (G.nodes[u][node_potential] - G.nodes[v][node_potential])) for u, v in G.edges()}

    return edge_abs_flow_dict    
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
def compute_directed_graph_of_flow(
        G,
        edge_weight=None,
        node_potential='p',
        edge_flow='q'):
    """
    Computes a directed graph of flows.
    
    Given a graph `G` with edge conductances `w` (edge attribute) and node potential `p` 
    (node attribute), this function returns a directed graph `G_dir`, with the same edges 
    as in `G` but directed according to flow direction, with edge flows `q` (edge attribute), 
    based on the relation `q(u, v) = - w(u, v) * (p(u) - p(v))` for an edge `(u, v)`.

    Parameters
    ----------
    G : networkx.Graph
        graph
    
    edge_weight : str, optional
        name of the edge attribute used as weight (conductance); 
        by default (`None`) : all edges have a weight of 1

    node_potential : str, default: 'p'
        name of the node attribute used as potential

    edge_flow : str, default: 'q'
        name of the edge attribute used to store the (positive) flow on edges
        of the output directed graph

    Returns
    -------
    G_dir : networkx.DiGraph
        directed graph, obtained by conversion of `G` (with the method `to_directed`), and where
        only the edges with positive (or zero) flow have been kept, these flows are stored as an
        edge attribute named `edge_flow`
    """
    # Convert G to a directed graph (copy)
    G_dir = G.to_directed()
    # Note : for an edge `(u, v)` in `G`: the two edges `(u, v)` and `(v, u)` are in `G_dir`

    # Compute flows along edges
    if edge_weight is None:
        edge_flow_dict = {(u, v): G_dir.nodes[u][node_potential] - G_dir.nodes[v][node_potential] for u, v in G_dir.edges()}
    else:
        edge_flow_dict = {(u, v): G_dir.edges[u, v][edge_weight] * (G_dir.nodes[u][node_potential] - G_dir.nodes[v][node_potential]) for u, v in G_dir.edges()}

    # Set edge flow as edge attribute in `G_dir`
    nx.set_edge_attributes(G_dir, edge_flow_dict, edge_flow)

    # Remove edges with negative flows
    edgelist_to_remove = {k for k, v in edge_flow_dict.items() if v < 0}
    # edgelist_to_remove = {k for k, v in nx.get_edge_attributes(G_dir, edge_flow).items() if v < 0}
    G_dir.remove_edges_from(edgelist_to_remove)

    return G_dir
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
def equivalent_conductance(
        G, 
        pair_nodes_list, 
        edge_weight=None, 
        k=None, 
        eigval=None, 
        eigvec=None, 
        **kwargs):
    """
    Computes the equivalent conductance between pairs of nodes in a connected graph.

    Let `w(e)` the conductance of the edge `e` where `w` is given by 
    `edge_weight` (`w(e)=1` for any edge if not specified). This function
    computes the equivalent conductance between the two nodes given by 
    `pair_nodes_list[i]`, for all i, by using (an approximation of) the pseudo 
    inverse of the Laplacian matrix of the graph, which is determined by the 
    (first) eigen values and eigen vectors of the Laplacian matrix.

    Let :math:`L^{+}` be the pseudo inverse of the Laplacian matrix, i.e.
    
    .. math::
        L^{+}=\\sum_{i} \\frac{1}{\\lambda_i} u_i \\cdot u_i^\\top,
    
    where :math:`u_i` is a eigen vector of norm 1 for the eigen value
    :math:`\\lambda_i`, and where the sum goes through the non-zero eigen 
    values.

    Then, for a pair of nodes :math:`(a, b)`, with :math:`q=\delta_a-\delta_b`
    where :math:`\delta_a(u) = 1` if :math:`u=a` and :math:`\delta_a(u) = 0`
    otherwise, the equivalent conductance :math:`w_{eq}(a, b)` between 
    :math:`a` and :math:`b` satisfies
    
    .. math::
        w_{eq}(a, b)^{-1} = q^\\top \\cdot L^{+}\\cdot q

    i.e.

    .. math::
        w_{eq}(a, b)^{-1} = \\sum_{i} \\frac{1}{\\lambda_i} (u_i^\\top \\cdot q)^2
        
    This function uses the approximation of the pseudo inverse of 
    the Laplacian matrix by accounting for the first `k` eigen values and 
    eigen vectors (of the Laplacian matrix), see parameter `k` below.

    Note: it is recommended to pre-compute (first) eigen values and eigen vectors
    (and specify parameters `eigval` and `eigvec`).

    Parameters
    ----------
    G : networkx.Graph
        graph, must be connected (i.e. with one connected component)

    pair_nodes_list : list
        list of pair(s) of nodes (labels) in the graph `G`:

        - `pair_nodes_list[i] = [u1, u2]`, where `u1`, `u2` are two \
        nodes (labels) in the graph
    
    edge_weight : str, optional
        name of the edge attribute used as weight (conductance); 
        unused if `eigval` and `eigvec` are specified (not `None`)
        by default (`None`) : all edges have a weight of 1

    k : int, optional
        number (max) of eigen values and eigen vectors (of the Laplacian
        matrix) to compute, the first `k` eigen values in ascending order 
        are considered;
        unused if `eigval` and `eigvec` are specified (not `None`);
        by default (`None`): all eigen values and eigen vectors are 
        computed, i.e. `k` is set to `G.number_of_nodes()`

    eigval : 1d-array, optional
        array of shape (m, ), with m less than or equal to the number 
        of nodes in the graph, m first eigen values of the Laplacian 
        matrix, in ascending order, i.e. for a connected graph, 
        `0 = eigval[0] < eigenval[1] <= eigenval[2] ...`;
        by default (`None`): the first `k` eigen values are computed

    eigvec : 2d-array, optional
        array of shape (n, m), where n is the number of nodes in the
        graph and m <= n, m first eigen vectors in columns (orthogonal 
        and of norm 1) corresponding to the eigen values `eigval`; 
        for a connected graph, `eigvec[:,0]` is proportional to the 
        vector `[1, 1, ..., 1]`;
        by default (`None`): the first `k` eigen vectors are computed

    kwargs : dict
        keyword arguments passed to the function 
        `scipy.sparse.linalg.eigsh` to compute the `k` first eigen values
        and eigen vectors (unused if `k >= G.number_of_nodes()`, in this 
        case, the function `scipy.linalg.eigh` is used); 
        note: the parameter `which`, if not given in `kwargs`, is set to 
        'SM' (smallest eigen values are considered) 

    Returns
    -------
    weq : 1d-array of floats
        1d-array of same length as the length of the list `pair_nodes_list`:

        - `weq[i]` : equivalent conductance between the two nodes \
        `pair_nodes_list[i]`
    """
    if eigvec is None or eigval is None:
        eigval, eigvec = laplacian_eigs(G, edge_weight=edge_weight, k=k, **kwargs)

    # Set dictionary to convert node label to node index
    node_label2index = {k:i for i, k in enumerate(G.nodes())}
    # node_label2index = kn.utils.get_node_label2index_dict(G) # equivalent

    t = 1.0 / eigval[1:]
    
    weq = np.full((len(pair_nodes_list),), np.nan)
    for i, (u_in, u_out) in enumerate(pair_nodes_list):
        i_in = node_label2index[u_in]  # index of inlet node
        i_out = node_label2index[u_out] # index of outlet node
        weq[i] = np.sum(t*(eigvec[i_in, 1:] - eigvec[i_out, 1:])**2)
    
    weq = 1.0 / weq

    # # Equiv.
    # # ------
    # # Psendo inverse (approximation) of the Laplacian matrix
    # Lpinv = (eigvec[:,1:] * 1.0/eigval[1:]) @ eigvec[:,1:].T
    # # Lpinv = eigvec[:,1:] @ np.diag(1.0/eigval[1:]) @ eigvec[:,1:].T # equiv
    # # Lpinv = np.sum([1/eigval[i] * np.outer(eigvec[:, i], eigvec[:, i]) for i in range(1, len(eigval))], axis=0) # equiv

    # n = G.number_of_nodes()

    # weq = np.full((len(pair_nodes_list),), np.nan)
    # for i, (u_in, u_out) in enumerate(pair_nodes_list):
    #     # Set source term correpsonding to flux of 1 at inlet node (u_in) and flux of -1 at outlet node (u_out)  
    #     qext = np.zeros(n)
    #     qext[node_label2index[u_in]] = 1
    #     qext[node_label2index[u_out]] = -1
    #     # Compute equivalent conductance between inlet and outlet nodes
    #     weq[i] = 1.0 / (qext @ Lpinv @ qext)

    #     # p = Lpinv @ qext # solve L p = qext, p potential at nodes
    #     # pdiff = qext @ p # get potential difference between inlet and outlet nodes
    #     # weq[i] = 1.0/pdiff # equivalent conductance
    #     # # req[i] = 1.0/weq[i]  # equivalent resistance

    return weq
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
def equivalent_conductance_sets(
        G, 
        pair_node_sets_list, 
        edge_weight=None, 
        k=None, 
        eigval=None, 
        eigvec=None, 
        **kwargs):
    """
    Computes the equivalent conductance between pairs of disjoint node sets in a connected graph.

    The definition of the equivlalent conductance between two disjoint sets of nodes
    is based on the constraint that the difference between the average node potential
    over the first set and the average node potential over the second set is fixed.

    Let `w(e)` the conductance of the edge `e` where `w` is given by 
    `edge_weight` (`w(e)=1` for any edge if not specified). This function
    computes the equivalent conductance between the two disjoint sets of nodes given by 
    `pair_node_sets_list[i]`, for all i, by using (an approximation of) the pseudo 
    inverse of the Laplacian matrix of the graph, which is determined by the 
    (first) eigen values and eigen vectors of the Laplacian matrix.

    Let :math:`L^{+}` be the pseudo inverse of the Laplacian matrix, i.e.
    
    .. math::
        L^{+}=\\sum_{i} \\frac{1}{\\lambda_i} u_i \\cdot u_i^\\top,
    
    where :math:`u_i` is a eigen vector of norm 1 for the eigen value
    :math:`\\lambda_i`, and where the sum goes through the non-zero eigen 
    values.

    Then, for a pair of disjoint node sets :math:`(A, B)`, 
    with :math:`q=1/|A|\delta_A-1/|B|\delta_B`
    where :math:`\delta_A(u) = 1` if :math:`u\in a` and :math:`\delta_A(u) = 0`
    otherwise, the equivalent conductance :math:`w_{eq}(A, B)` between 
    :math:`A` and :math:`B` is defined as
    
    .. math::
        w_{eq}(A, B)^{-1} = q^\\top \\cdot L^{+}\\cdot q

    i.e.

    .. math::
        w_{eq}(A, B)^{-1} = \\sum_{i} \\frac{1}{\\lambda_i} (u_i^\\top \\cdot q)^2
        
    This function uses the approximation of the pseudo inverse of 
    the Laplacian matrix by accounting for the first `k` eigen values and 
    eigen vectors (of the Laplacian matrix), see parameter `k` below.
    
    Note: it is recommended to pre-compute (first) eigen values and eigen vectors
    (and specify parameters `eigval` and `eigvec`).

    Parameters
    ----------
    G : networkx.Graph
        graph, must be connected (i.e. with one connected component)

    pair_node_sets_list : list
        list of pair(s) of disjoint sets of nodes (labels) in the graph `G`:

        - `pair_nodes_list[i] = [A1, A2]`, where `A1`, `A2` are two \
        disjoint sets of nodes (labels) in the graph
    
    edge_weight : str, optional
        name of the edge attribute used as weight (conductance); 
        unused if `eigval` and `eigvec` are specified (not `None`)
        by default (`None`) : all edges have a weight of 1

    k : int, optional
        number (max) of eigen values and eigen vectors (of the Laplacian
        matrix) to compute, the first `k` eigen values in ascending order 
        are considered;
        unused if `eigval` and `eigvec` are specified (not `None`);
        by default (`None`): all eigen values and eigen vectors are 
        computed, i.e. `k` is set to `G.number_of_nodes()`

    eigval : 1d-array, optional
        array of shape (m, ), with m less than or equal to the number 
        of nodes in the graph, m first eigen values of the Laplacian 
        matrix, in ascending order, i.e. for a connected graph, 
        `0 = eigval[0] < eigenval[1] <= eigenval[2] ...`;
        by default (`None`): the first `k` eigen values are computed

    eigvec : 2d-array, optional
        array of shape (n, m), where n is the number of nodes in the
        graph and m <= n, m first eigen vectors in columns (orthogonal 
        and of norm 1) corresponding to the eigen values `eigval`; 
        for a connected graph, `eigvec[:,0]` is proportional to the 
        vector `[1, 1, ..., 1]`;
        by default (`None`): the first `k` eigen vectors are computed

    kwargs : dict
        keyword arguments passed to the function 
        `scipy.sparse.linalg.eigsh` to compute the `k` first eigen values
        and eigen vectors (unused if `k >= G.number_of_nodes()`, in this 
        case, the function `scipy.linalg.eigh` is used); 
        note: the parameter `which`, if not given in `kwargs`, is set to 
        'SM' (smallest eigen values are considered) 

    Returns
    -------
    weq : 1d-array of floats
        1d-array of same length as the length of the list `pair_node_sets_list`:

        - `weq[i]` : equivalent conductance between the two disjoint set of nodes \
        `pair_node_sets_list[i]`
    """
    # Set dictionary to convert node label to node index
    node_label2index = {k:i for i, k in enumerate(G.nodes())}
    # node_label2index = kn.utils.get_node_label2index_dict(G) # equivalent

    # Convert node label to node index
    pair_node_index_sets_list = [[[node_label2index[k] for k in A], [node_label2index[k] for k in B]] for A, B in pair_node_sets_list]

    # Check
    for A_index, B_index in pair_node_index_sets_list:
        if len(np.unique(A_index)) != len(A_index) or len(np.unique(B_index)) != len(B_index):
            raise ValueError('A set of nodes contains duplicated nodes')
        if len(np.intersect1d(A_index, B_index)):
            raise ValueError('A pair of sets of nodes is not disjoint (non empty intersection)')

    if eigvec is None or eigval is None:
        eigval, eigvec = laplacian_eigs(G, edge_weight=edge_weight, k=k, **kwargs)
    
    t = 1.0 / eigval[1:]
    
    weq = np.full((len(pair_node_index_sets_list),), np.nan)
    for i, (A_index, B_index) in enumerate(pair_node_index_sets_list):
        weq[i] = np.sum(t * (eigvec[A_index, 1:].mean(axis=0) - eigvec[B_index, 1:].mean(axis=0))**2)

    # # Equiv.    
    # weq = np.full((len(pair_node_index_sets_list),), np.nan)
    # for i, (A_index, B_index) in enumerate(pair_node_index_sets_list):
    #     n_A = len(A_index)
    #     n_B = len(B_index)
    #     qext = np.zeros(G.number_of_nodes())
    #     for k in A_index:
    #         qext[k] = 1.0/n_A
    #     for k in B_index:
    #         qext[k] = -1.0/n_B
    #     weq[i] = qext @ eigvec[:, 1:] @ (1.0/eigval[1:] * (eigvec[:, 1:].T @ qext))

    weq = 1.0 / weq

    return weq
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
def graph_conductance_exponent(
        G, 
        npairs,
        nclasses=21,
        class_in_log_scale=True,
        data_points_removed_first_fraction=0.66,
        data_points_removed_last_fraction=0.0,
        confidence_level=0.95,
        edge_weight=None, 
        k=None, 
        eigval=None, 
        eigvec=None,
        return_poly_fit_params=True,
        return_poly_fit_params_cov=True,
        return_index_used_for_fit=True,
        return_dist_mean=True,
        return_weq_mean=True,
        return_nb_points_per_class=True,
        return_class_lim=True,
        return_dist_full=False,
        return_weq_full=False,
        return_class_id_full=False,
        return_pair_nodes_list=False,
        seed=None,
        verbose=1,
        make_plot=True,
        plot_in_log_scale=True,
        plot_full_data_points=True,
        **kwargs):
    """
    Computes the conductance exponent (also called resistance exponent) of a graph.

    This function estimates the equivalent conductance (or resistance) for several 
    pairs of nodes, randomly selected; then, assuming the relationship
    
    .. math::
        \omega_{eq} \propto r^{-d_{\Omega}}
    
    or

    .. math::
        R_{eq} \propto r^{d_{\Omega}}

    where :math:`\omega_{eq}` (resp. :math:`R_{eq}`) is the equivalent conductance
    (resp. equivalent resistance) between two nodes, and :math:`r` the Euclidean 
    distance between the two nodes, the exponent (called resistance exponent or 
    conductance exponent) :math:`d_{\Omega}` is obtained by fitting 
    a line on the log-log plot of :math:`\omega_{eq}` as function of :math:`r`.

    This function uses the approximation of the pseudo inverse of 
    the Laplacian matrix by accounting for the first `k` eigen values and 
    eigen vectors (of the Laplacian matrix), see parameter `k` below.

    Note: it is recommended to pre-compute (first) eigen values and eigen vectors
    (and specify parameters `eigval` and `eigvec`).

    Parameters
    ----------
    G : networkx.Graph
        graph, must be connected (i.e. with one connected component)

    npairs : int
        number of pairs of nodes to select
    
    nclasses : int, default: 21
        number of classes for distances; let r(u, v) the Euclidean distance 
        between the nodes u and v; let rmin and rmax the minimum and maximum
        respectively of r(u,v) over all selected pairs of nodes (u, v); the 
        interval [rmin, rmax] is divided into `nclasses` sub-intervals (classes)
        of same length (in log scale if `class_in_log_scale=True`); then
        each pair of nodes is assigned to one class according to the distance
        between the two nodes; mean distance and mean equivalent conductance 
        on each class will then be computed

    class_in_log_scale : bool, default: True
        - if `True`: the classes of distances are set in log scale
        - if `False`: the classes of distances are set in ususal scale
        
    data_points_removed_first_fraction : float, default : 0.66
        see `data_points_removed_last_fraction`

    data_points_removed_last_fraction : float, default : 0.0
        `data_points_removed_first_fraction` and `data_points_removed_last_fraction` are
        positive numbers in (0, 1); sequences of mean distances and mean 
        equivalent conductances (according to the classes determined by 
        `nclasses`), constitute the data points; 
        the fraction `data_points_removed_first_fraction`, resp.
        `data_points_removed_last_fraction`, at the beginning, resp. end, of the
        data points are removed before doing the fit (line on the log-log plot);
        this allows to deal with the "extreme" data points

    confidence_level : float, default: 0.95
        confidence level, float in the interval (0, 1), to compute the 
        uncertainty for the line fitting (and for the conductance exponent)

    edge_weight : str, optional
        name of the edge attribute used as weight (conductance); 
        unused if `eigval` and `eigvec` are specified (not `None`)
        by default (`None`) : all edges have a weight of 1

    k : int, optional
        number (max) of eigen values and eigen vectors (of the Laplacian
        matrix) to compute, the first `k` eigen values in ascending order 
        are considered;
        unused if `eigval` and `eigvec` are specified (not `None`);
        by default (`None`): all eigen values and eigen vectors are 
        computed, i.e. `k` is set to `G.number_of_nodes()`

    eigval : 1d-array, optional
        array of shape (m, ), with m less than or equal to the number 
        of nodes in the graph, m first eigen values of the Laplacian 
        matrix, in ascending order, i.e. for a connected graph, 
        `0 = eigval[0] < eigenval[1] <= eigenval[2] ...`;
        by default (`None`): the first `k` eigen values are computed

    eigvec : 2d-array, optional
        array of shape (n, m), where n is the number of nodes in the
        graph and m <= n, m first eigen vectors in columns (orthogonal 
        and of norm 1) corresponding to the eigen values `eigval`; 
        for a connected graph, `eigvec[:,0]` is proportional to the 
        vector `[1, 1, ..., 1]`;
        by default (`None`): the first `k` eigen vectors are computed

    return_poly_fit_params : bool, default: True
        if `True`, the array (shape (2, )) of parameters of the fitted 
        polynom (line) is returned

    return_poly_fit_params_cov : bool, default: True
        if `True`, the array (shape (2, 2)) of parameters covariance of the 
        fitted polynom (line) is returned
    
    return_index_used_for_fit: bool, default: True
        if `True`, the starting (included) and ending (excluded) indices of the sequence 
        of data points used for fitting is returned

    return_dist_mean : bool, default: True
        if `True`, the array of mean distances (in each class) is returned
        (mean over all pair of nodes in one class, for each class)

    return_weq_mean : bool, default: True
        if `True`, the array of mean equivalent conductances (in each class) 
        is returned (mean over all pair of nodes in one class, for each class)

    return_nb_points_per_class : bool, default: True
        if `True`, number of points (pairs) per class is returned

    return_class_lim : bool, default: True
        if `True`, array of class bounds is returned

    return_dist_full : bool, default: False
        if `True`, the array of all distances is returned

    return_weq_full : bool, default: False
        if `True`, the array of all equivalent conductances is returned

    return_class_id_full : bool, default: False
        if `True`, the array of class id of all points (pairs) is returned

    return_pair_nodes_list : bool, default: False
        if `True`, the list of selected pairs of nodes (labels) is returned

    seed : int, optional
        seed number for initializing the random number generator

    verbose : int, default: 1
        verbose mode, larger value implies more info printed

    make_plot : bool, default: True
        indicates if the plot of the data points and fitting curve is displayed 
        (in the current figure axis)

    plot_in_log_scale : bool, default: True
        used if `make_plot=True`, indicates if the plot is in log-log scale
        (along x- and y-axes)

    plot_full_data_points : bool, default: True
        used if `make_plot=True`, indicates if full data points (pairs) are
        displayed on the plot

    kwargs : dict
        keyword arguments passed to the function 
        `scipy.sparse.linalg.eigsh` to compute the `k` first eigen values
        and eigen vectors (unused if `k >= G.number_of_nodes()`, in this 
        case, the function `scipy.linalg.eigh` is used); 
        note: the parameter `which`, if not given in `kwargs`, is set to 
        'SM' (smallest eigen values are considered) 

    Returns
    -------
    out_dict : dict
        dictionary with the following keys / values:
        
        - dOmega : float
            conductance exponent (or resistance exponent)
    
        - dOmega_delta : float
            uncertainty for conductance exponent, the interval 
            `dOmega +/- dOmega_delta` corresponds to the confidence interval 
            at a confidence level `confidence_level` derived from Gaussian 
            distribution

        - confidence_level : float
            confidence level used to compute the uncertainty for the line fitting 
            
        - poly_fit_params : array of shape (2, ), optional
            returned if `return_poly_fit_params=True`:
            array of parameters of the fitted polynom (line); in particular:
            `dOmega = -poly_fit_params[0]`

        - poly_fit_params_cov : array of shape (2, 2), optional
            returned if `return_poly_fit_params_cov=True`:
            array (shape (2, 2)) of parameters covariance of the fitted polynom 
            (line); in particular: with
            `dOmega_delta = scipy.stats.norm.ppf((1+confidence_level)/2) * np.sqrt(poly_fit_params_cov[0, 0])`
        
        - index_used_for_fit : list of two ints, optional
            returned if `return_index_used_for_fit=True`:
            starting (included) index and ending (excluded) index of the sequence of
            data points used for fitting

        - dist_mean : 1d-array of shape (nclasses,), optional
            returned if `return_dist_mean=True`:
            array of mean on each class of Euclidean distances between the 
            two nodes of the selected pair of nodes

        - weq_mean : 1d-array of shape (nclasses,), optional
            returned if `return_weq_mean=True`:
            array of mean on each class of equivalent conductance for the
            selected pair of nodes

        - nb_points_per_class : 1d-array of shape (nclasses,), optional
            returned if `return_nb_points_per_class=True`:
            array of number of points (pairs) in each class

        - class_lim : 1d-array of shape (nclasses+1,), optional
            returned if `return_class_lim=True`:
            array of class limits: increasing float numbers determining the
            classes

        - dist_full : 1d-array of shape (npairs,), optional
            returned if `return_dist_full=True`:
            array of Euclidean distances between the two nodes of all selected 
            pair of nodes

        - weq_full : 1d-array of shape (npairs,), optional
            returned if `return_weq_full=True`:
            array of equivalent conductances between the two nodes of all selected 
            pair of nodes (corresponding to distances in `dist_full`)

        - class_id_full : 1d-array of shape (npairs, ), optional
            returned if `return_class_id_full=True`:
            array of class id (int in (0, ..., nclass-1)) of all points (pairs)

        - pair_nodes_list : list of length npairs, optional
            returned if `return_pair_nodes_list=True`:
            list of selected pairs of nodes (labels) in the graph `G`:

            - `pair_nodes_list[i] = [u1, u2]`, where `u1`, `u2` are two \
            nodes (labels) in the graph
    """
    if seed is not None:
        np.random.seed(seed)
    
    if nclasses < 1:
        print('No class!')
        return
    
    if npairs < nclasses:
        print(f'Two few pairs ({npairs}, less than number of classes ({nclasses}))')
        return
    
    # Set dictionary to convert node index to node label
    node_index2label = {i:u for i, u in enumerate(G.nodes())}
    # node_index2label = kn.utils.get_node_index2label(G) # equivalent

    n_nodes = G.number_of_nodes()
    pair_index_list = [np.random.choice(n_nodes, size=2, replace=False) for _ in range(npairs)]
    pair_nodes_list = [(node_index2label[i], node_index2label[j]) for i, j in pair_index_list]
 
    # Compute the Euclidean distance between the two nodes for each pair
    dist_full = np.asarray([np.asarray(G.nodes[x]['pos']) - np.asarray(G.nodes[y]['pos']) for x, y in pair_nodes_list])
    dist_full = np.sqrt(np.sum(dist_full**2, axis=1))
    
    # Compute the equivalent conductance between the two nodes for each pair
    if verbose > 0:
        if eigvec is None or eigval is None:
            print(f'Compute Laplacian pseudo inverse and equivalent conductance ({npairs} pair of nodes)')
        else:
            print(f'Compute equivalent conductance ({npairs} pair of nodes) (Laplacian pseudo inverse pre-computed)')
    
    t1_all = time.time()
    weq_full = kn.tools.equivalent_conductance(G, pair_nodes_list, edge_weight=edge_weight, k=k, eigval=eigval, eigvec=eigvec, **kwargs)
    t2_all = time.time()
    if verbose > 0:
        print(f'Total elapsed time (w_eq for pair of nodes): {t2_all-t1_all:.2f} sec')

    if verbose > 0:
        print(f'Divide distances into classes, and compute mean on each class (of distance and of equivalent conductance)')
    
    rmin = dist_full.min()
    rmax = dist_full.max()
    if class_in_log_scale:
        class_lim = np.linspace(np.log10(rmin), np.log10(rmax), nclasses+1)
        class_id_full = np.digitize(np.log10(dist_full), class_lim) - 1
    else:
        class_lim = np.linspace(rmin, rmax, nclasses+1)
        class_id_full = np.digitize(dist_full, class_lim) - 1

    nb_points_per_class = np.asarray([np.sum(class_id_full == i) for i in range(nclasses)])

    # dist_mean = np.asarray([np.mean(dist_full[class_id_full == i]) for i in range(nclasses)])
    # weq_mean = np.asarray([np.mean(weq_full[class_id_full == i]) for i in range(nclasses)])

    dist_mean = np.full(nclasses, np.nan)
    weq_mean = np.full(nclasses, np.nan)
    for i in range(nclasses):
        if nb_points_per_class[i] > 0:
            dist_mean[i] = np.mean(dist_full[class_id_full==i])    
            weq_mean[i] = np.mean(weq_full[class_id_full==i])    
    
    if verbose > 0:
        print(f'Compute conductance exponent (line fitting on log-log plot of equivalent conductance as function of distance)')

    # Compute conductance exponent (line fitting on log-log plot of equivalent conductance as function of distance)
    # -------------------------------------------------------------------------------------------------------------
    # Data (all)
    x_all = dist_mean
    y_all = weq_mean

    poly_fit_params, poly_fit_params_cov, x, y, n1, n2 = _loglog_fit(x_all, y_all, data_points_removed_first_fraction, data_points_removed_last_fraction, verbose=verbose)

    # Get conductance exponent dOmega, and dOmega_delta
    if poly_fit_params is not None:
        dOmega = - poly_fit_params[0] # - slope
        fit_ok = True
    else:
        dOmega = np.nan
        fit_ok = False
    
    if poly_fit_params_cov is not None:
        t = scipy.stats.norm.ppf((1+confidence_level)/2)
        dOmega_delta = t * np.sqrt(poly_fit_params_cov[0, 0]) # slope_delta
    else:
        dOmega_delta = np.nan

    # Make plot (if needed)
    if make_plot:
        plot_fit_conductance_exponent(
            x_all, y_all, x, y, dist_full, weq_full, class_id_full, nb_points_per_class, dOmega, dOmega_delta, confidence_level, poly_fit_params, 
            plot_fit=fit_ok, plot_in_log_scale=plot_in_log_scale, plot_full_data_points=plot_full_data_points)

    # Set output dictionary
    out_dict = {'dOmega': dOmega, 'dOmega_delta': dOmega_delta, 'confidence_level': confidence_level}

    if return_poly_fit_params:
        out_dict['poly_fit_params'] = poly_fit_params

    if return_poly_fit_params_cov:
        out_dict['poly_fit_params_cov'] = poly_fit_params_cov

    if return_index_used_for_fit:
        out_dict['index_used_for_fit'] = [n1, n2]

    if return_dist_mean:
        out_dict['dist_mean'] = dist_mean

    if return_weq_mean:
        out_dict['weq_mean'] = weq_mean

    if return_nb_points_per_class:
        out_dict['nb_points_per_class'] = nb_points_per_class

    if return_class_lim:
        out_dict['class_lim'] = class_lim

    if return_dist_full:
        out_dict['dist_full'] = dist_full

    if return_weq_full:
        out_dict['weq_full'] = weq_full

    if return_class_id_full:
        out_dict['class_id_full'] = class_id_full

    if return_pair_nodes_list:
        out_dict['pair_nodes_list'] = pair_nodes_list

    return out_dict
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
def graph_conductance_exponent_update(
        out_dict,
        nclasses=21,
        class_in_log_scale=True,
        data_points_removed_first_fraction=0.33,
        data_points_removed_last_fraction=0.0,
        confidence_level=0.95,
        verbose=1,
        make_plot=True,
        plot_in_log_scale=True,
        plot_full_data_points=True):
    """
    Updates the fit (power law) for the conductance exponent.

    It is assumed that the function `graph_conductance_exponent` has 
    already been run, and that the output dictionary `out_dict` contains the
    keys `dist_full`, `weq_full`; 
    the fit is updated by changing the classes and / or the fractions of 
    data points removed at the beginning and end of the sequences of data points.

    See function `graph_conductance_exponent` for details.

    Note: this function operates inplace on the `out_dict` dictionary (it is
    updated and returned).

    Parameters
    ----------
    out_dict : dict
        output dictionary of the function `graph_conductance_exponent`;
        the dictionary should contain the key `xxx`

    nclasses : int, default: 21
        number of classes for distances; let r(u, v) the Euclidean distance 
        between the nodes u and v; let rmin and rmax the minimum and maximum
        respectively of r(u,v) over all selected pairs of nodes (u, v); the 
        interval [rmin, rmax] is divided into `nclasses` sub-intervals (classes)
        of same length (in log scale if `class_in_log_scale=True`); then
        each pair of nodes is assigned to one class according to the distance
        between the two nodes; mean distance and mean equivalent conductance 
        on each class will then be computed

    class_in_log_scale : bool, default: True
        - if `True`: the classes of distances are set in log scale
        - if `False`: the classes of distances are set in ususal scale

    data_points_removed_first_fraction : float, default : 0.33
        see `data_points_removed_last_fraction`

    data_points_removed_last_fraction : float, default : 0.0
        `data_points_removed_first_fraction` and `data_points_removed_last_fraction` are
        positive numbers in (0, 1); sequences of time steps and mean square distances
        constitute the data points; 
        the fraction `data_points_removed_first_fraction`, resp.
        `data_points_removed_last_fraction`, at the beginning, resp. end, of the
        data points are removed before doing the fit (line on the log-log plot);
        this allows to deal with the "extreme" data points

    confidence_level : float, default: 0.95
        confidence level, float in the interval (0, 1), to compute the 
        uncertainty for the line fitting (and for the random walk dimension)

    verbose : int, default: 1
        verbose mode, larger value implies more info printed

    make_plot : bool, default: True
        indicates if the plot of the data points and fitting curve is displayed 
        (in the current figure axis)

    plot_in_log_scale : bool, default: True
        used if `make_plot=True`, indicates if the plot is in log-log scale
        (along x- and y-axes)

    plot_full_data_points : bool, default: True
        used if `make_plot=True`, indicates if full data points (pairs) are
        displayed on the plot

    Returns
    -------
    out_dict : dict
        updated output dictionary
    """
    if 'dist_full' not in out_dict.keys() or \
        'weq_full' not in out_dict.keys():
        raise ValueError('The input dictionary `out_dict` should contain the keys `dist_full` and `weq_full`')

    if nclasses < 1:
        print('No class!')
        return
    
    dist_full = out_dict['dist_full']
    weq_full = out_dict['weq_full']

    npairs = len(dist_full)
    if npairs < nclasses:
        print(f'Two few pairs ({npairs}, less than number of classes ({nclasses}))')
        return
    
    if verbose > 0:
        print(f'Divide distances into classes, and compute mean on each class (of distance and of equivalent conductance)')
    
    rmin = dist_full.min()
    rmax = dist_full.max()
    if class_in_log_scale:
        class_lim = np.linspace(np.log10(rmin), np.log10(rmax), nclasses+1)
        class_id_full = np.digitize(np.log10(dist_full), class_lim) - 1
    else:
        class_lim = np.linspace(rmin, rmax, nclasses+1)
        class_id_full = np.digitize(dist_full, class_lim) - 1

    nb_points_per_class = np.asarray([np.sum(class_id_full == i) for i in range(nclasses)])

    # dist_mean = np.asarray([np.mean(dist_full[class_id_full == i]) for i in range(nclasses)])
    # weq_mean = np.asarray([np.mean(weq_full[class_id_full == i]) for i in range(nclasses)])

    dist_mean = np.full(nclasses, np.nan)
    weq_mean = np.full(nclasses, np.nan)
    for i in range(nclasses):
        if nb_points_per_class[i] > 0:
            dist_mean[i] = np.mean(dist_full[class_id_full==i])    
            weq_mean[i] = np.mean(weq_full[class_id_full==i])    
    
    if verbose > 0:
        print(f'Compute conductance exponent (line fitting on log-log plot of equivalent conductance as function of distance)')

    # Compute conductance exponent (line fitting on log-log plot of equivalent conductance as function of distance)
    # -------------------------------------------------------------------------------------------------------------
    # Data (all)
    x_all = dist_mean
    y_all = weq_mean

    poly_fit_params, poly_fit_params_cov, x, y, n1, n2 = _loglog_fit(x_all, y_all, data_points_removed_first_fraction, data_points_removed_last_fraction, verbose=verbose)

    # Get conductance exponent dOmega, and dOmega_delta
    if poly_fit_params is not None:
        dOmega = - poly_fit_params[0] # - slope
        fit_ok = True
    else:
        dOmega = np.nan
        fit_ok = False
    
    if poly_fit_params_cov is not None:
        t = scipy.stats.norm.ppf((1+confidence_level)/2)
        dOmega_delta = t * np.sqrt(poly_fit_params_cov[0, 0]) # slope_delta
    else:
        dOmega_delta = np.nan

    # Make plot (if needed)
    if make_plot:
        plot_fit_conductance_exponent(
            x_all, y_all, x, y, dist_full, weq_full, class_id_full, nb_points_per_class, dOmega, dOmega_delta, confidence_level, poly_fit_params, 
            plot_fit=fit_ok, plot_in_log_scale=plot_in_log_scale, plot_full_data_points=plot_full_data_points)

    # Update output dictionary
    out_dict['dOmega'] = dOmega
    out_dict['dOmega_delta'] = dOmega_delta
    out_dict['confidence_level'] = confidence_level

    if 'poly_fit_params' in out_dict.keys():
        out_dict['poly_fit_params'] = poly_fit_params

    if 'poly_fit_params_cov' in out_dict.keys():
        out_dict['poly_fit_params_cov'] = poly_fit_params_cov

    if 'index_used_for_fit' in out_dict.keys():
        out_dict['index_used_for_fit'] = [n1, n2]

    if 'dist_mean' in out_dict.keys():
        out_dict['dist_mean'] = dist_mean

    if 'weq_mean' in out_dict.keys():
        out_dict['weq_mean'] = weq_mean

    if 'nb_points_per_class' in out_dict.keys():
        out_dict['nb_points_per_class'] = nb_points_per_class

    if 'class_lim' in out_dict.keys():
        out_dict['class_lim'] = class_lim

    if 'class_id_full' in out_dict.keys():
        out_dict['class_id_full'] = class_id_full

    return out_dict
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
def plot_fit_conductance_exponent(
        x_all,
        y_all,
        x,
        y,
        x_full,
        y_full,
        class_id_full,
        nb_points_per_class,
        dOmega,
        dOmega_delta,
        confidence_level,
        poly_fit_params,
        plot_fit=True,
        plot_in_log_scale=True,
        plot_full_data_points=True):
    """
    Makes the plot (in current axis figure) of the data points and fitted curve for conductance exponent (power law).

    Parameters
    ----------
    x_all : array of floats
        x-coordinates of all data points (classes)

    y_all : array of floats
        y-coordinates of all data points (classes)

    x : array of floats
        x-coordinates of data points used for fitting (a subset of `x_all`)
    
    y : array of floats
        y-coordinates of data points used for fitting (a subset of `y_all`)

    x_full : array of floats
        x-coordinates of single data points ("full", before making the classes)
    
    y_full : array of floats
        y-coordinates of single data points ("full", before making the classes)

    class_id_full : array of ints
        class id of single data points ("full", id of the class to which each 
        data point belongs)
    
    nb_points_per_class : array of ints
        number of data points per classes, could be calculated by
        `nb_points_per_class = np.asarray([np.sum(class_id_full == i) for i in range(nclasses)])`

    drw : float
        random walk dimension
    
    drw_lim_min : float
        lower limit for random walk dimension, provided from the 
        uncertainty on the slope of the fitted line (where a confidence 
        interval at a confidence level `confidence_level` derived from 
        Gaussian distribution is used)

    drw_lim_max : float
        upper limit for random walk dimension, provided from the 
        uncertainty on the slope of the fitted line (where a confidence 
        interval at a confidence level `confidence_level` derived from 
        Gaussian distribution is used)

    confidence_level : float
        confidence level, float in the interval (0, 1), used to compute the 
        uncertainty for the line fitting (and for the fractal dimension)

    poly_fit_params : array of shape (2, ) or None
        array of parameters of the fitted polynom (line); in particular:
        `df = -poly_fit_params[0]`; if no fit was done, `poly_fit_params` is `None`

    plot_fit : bool ; default : True
        indicates if the plot of the fitted curve is done (in the current figure axis);
        if `False`, only the data points are plotted

    plot_in_log_scale : bool, default: True
        indicates if the plot is in log-log scale (along x- and y-axes)

    plot_full_data_points : bool, default: True
        indicates if full data points (pairs) are displayed on the plot

    Returns
    -------
    None
    """
    dOmega_lim_min = dOmega - dOmega_delta
    dOmega_lim_max = dOmega + dOmega_delta
    
    nclasses = len(x_all)

    # Estimation with the fitted line
    if plot_fit:
        xx = np.linspace(x.min(), x.max(), 200)
        yy = np.power(10, poly_fit_params[0] * np.log10(xx) + poly_fit_params[1])

        xx_all = np.linspace(np.nanmin(x_all), np.nanmax(x_all), 200)
        yy_all = np.power(10, poly_fit_params[0] * np.log10(xx_all) + poly_fit_params[1])

    # Plot parameters
    xlabel = '$r$ : distance'
    ylabel = '$\omega_{eq}$ : equivalent conductance'

    title = '$\omega_{eq}\propto r^{-d_{\Omega}}$' + ', $d_{\Omega}$' + f'={dOmega:.5f} in [{dOmega_lim_min:.5f}, {dOmega_lim_max:.5f}] [{100*confidence_level:.2g}% - inter.]'

    color_data_all      = 'tab:blue'
    marker_data_all     = 'o'
    markersize_data_all = 5

    color_data_used      = 'tab:orange'
    marker_data_used     = 'o'
    markersize_data_used = 5

    color_fit_data_used = 'red'
    ls_fit_data_used    = 'solid'
    lw_fit_data_used    = 1.

    color_fit_data_all = 'red'
    ls_fit_data_all    = 'dotted'
    lw_fit_data_all    = 1.

    # Figure
    # ------
    # plt.figure(figsize=(10, 6))
    if plot_full_data_points:
        # Add single data points (pair) (before computing mean on classes)
        marker_single     = '.'
        markersize_single = 1

        cmap = plt.get_cmap('viridis')

        label_single = 'full data point (pair)'
        for i in range(nclasses):
            plt.plot(x_full[class_id_full==i], y_full[class_id_full==i], alpha=.3, ls='', marker=marker_single, markersize=markersize_single, color=cmap(i/(nclasses-1)), label=label_single)
            label_single = None
    
    plt.plot(x_all, y_all, ls='', marker=marker_data_all , markersize=markersize_data_all , color=color_data_all , label='data')
    if len(x) > 0:
        plt.plot(x    , y    , ls='', marker=marker_data_used, markersize=markersize_data_used, color=color_data_used, label='data used for fitting')
    if plot_fit:
        plt.plot(xx_all, yy_all, ls=ls_fit_data_all , lw=lw_fit_data_all , color=color_fit_data_all)
        plt.plot(xx    , yy    , ls=ls_fit_data_used, lw=lw_fit_data_used, color=color_fit_data_used, label='fit')
    
    # Add number of points in each class (text)
    for xi, yi, ni in zip(x_all, y_all, nb_points_per_class):
        if not np.isnan(xi) and not np.isnan(yi):
            plt.text(xi, yi, f'{ni}', ha='left', va='bottom')
    
    if plot_in_log_scale:
        plt.xscale('log')
        plt.yscale('log')
    
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid()
    plt.legend()
    # plt.show()

    return
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
def invert_Req(
        G,
        Req_target_dict,
        T_estimate_attr='T_estimate',
        T_estimate_start_dict=None,
        k=None,
        nit_max=100,
        tol=1.e-10,
        tol_step=1.e-10,
        T_clip_min=1.e-20,
        rescale_T=True,
        Req_err_type='SSE',
        return_Req_estimate_dict=True,
        return_Req_err=True,
        return_nit=False,
        return_Req_estimate_arr_list=False,
        return_Req_err_list=False,
        return_T_estimate_arr_list=False,
        verbose=0):
    """
    Inverts edges equivalent resistances.

    This function inverts equivalent resistances (:math:`R_{eq}`) given on every edge 
    of a connected graph by a Gauss-Newton optimisation scheme, i.e. it estimates the 
    conductances (`T`) on each edge such that the resulting edge equivalent 
    resistances are (as close as possible to) the target ones.

    This function uses the approximation of the pseudo inverse of 
    the Laplacian matrix by accounting for the first `k` eigen values and 
    eigen vectors (of the Laplacian matrix), see parameter `k` below.
    Note that edge weights for (the Laplacian matrix) is set to `T` estimates, 
    which change at each iteration.

    This function operates inplace on the graph `G`: the edge attribute
    `T_estimate_attr` contains the (final) edge conductance estimates.
    
    Parameters
    ----------
    G : networkx.Graph
        graph, must be connected (i.e. with one connected component)

    Req_target_dict : dict
        dictionary of target equivalent resistances (Req) on edges, where 
        the keys are the edges of `G` and the values are the target Req
    
    T_estimate_attr : str, default: 'T_estimate'
        name of the edge attribute used to store the (final) edge conductance estimates

    T_estimates_start_dict : dict, optional
        dictionary of starting edge conductance estimates, where 
        the keys are the edges of `G` and the values are the conductances;
        by default (`None`): starting conductance for edge uv, `T(uv)`, is set proportional
        to `1.0 / 0.5*(G.degree[u] + G.degree[v]) * Req_target_dict[(u, v)])` with
        the sum of `T(uv)*Req_target_dict[(u, v)]` over all edges `uv` is equal to the 
        number of graph nodes minus 1 (according to leverage scores property)

    k : int, optional
        number (max) of eigen values and eigen vectors (of the Laplacian
        matrix) to compute, the first `k` eigen values in ascending order 
        are considered;
        by default (`None`): all eigen values and eigen vectors are 
        computed, i.e. `k` is set to `G.number_of_nodes()`
        
    nit_max : int, default: 100
        maximum number of iteration(s) (in Gauss-Newton algorithm), or maximum
        number of updates of the edge conductance estimates

    tol : float, default: 1.e-10
        positive float number, the tolerance for the error between target Req and 
        Req estimates, stop criterion for the Gauss-Newton algorithm
    
    tol_sep : float, default: 1.e-10
        positive float number, the tolerance for the relative absolute difference
        between two successive errors (between target Req and Req estimates), 
        stop criterion for the Gauss-Newton algorithm
    
    T_clip_min : float, default: 1.e-20
        the new estimate of `T` are clipped with min value `T_clip_min` at each 
        iteration of the Gauss-Newton algorithm (useful to avoid negative value)
    
    rescale_T : bool, default: True
        if `True`, the new estimate of `T` (after clipping) is rescaled such that
        the sum of `T(uv)*Req_target_dict[(u, v)]` over all edges `uv` is equal to the 
        number of graph nodes minus 1 (according to leverage scores property)

    Req_err_type : str {'MaxAE'|'MAE'|'MSE'|'RMSE'|'SAE'|'SSE'}, default:'SSE'
        type of error computed (between target Req and Req estimates):

        - 'MaxAE' : Maximum Absolute Error
        - 'MAE' : Mean Absolute Error
        - 'MSE' : Mean Square Error
        - 'RMSE' : Root Mean Square Error
        - 'SAE' : Sum of Absolute Error
        - 'SSE' : Sum of Square Error
    
    return_Req_estimate_dict : bool, default: True
        if `True`, the (final) edge equivalent resistance estimates corresponding to 
        the (final) edge conductance estimates are returned

    return_Req_err : bool, default: True
        if `True`, the error between the (final) edge equivalent resistances and
        the target ones is returned

    return_nit : bool, default: False
        if `True`, the number of iteration(s) done (in Gaussian-Newton algorithm) 
        is returned

    return_Req_estimate_arr_list : bool, default: False
        if `True`, the list of arrays of edge equivalent resistance estimates, 
        at each iteration, is returned

    return_Req_err_list : bool, default: False
        if `True`, the list of errors between the edge equivalent resistance estimates and
        the target ones, at each iteration, is returned

    return_T_estimate_arr_list : bool, default: False
        if `True`, the list of arrays of edge conductance estimates, at each iteration, 
        is returned

    Returns
    -------
    Req_estimate_dict : dict, optional
        if `return_Req_estimate_dict=True`, the dictionary of the (final) edge equivalent 
        resistances (corresponding to the (final) edge conductance estimates (see parameter 
        `T_estimate_attr`)): the keys are the edges, and the values are the (final) 
        equivalent resistance estimates

    Req_err : float, optional
        if `return_Req_err=True`, the error between the (final) edge equivalent resistance
        estimates and the target ones (type of error given by `Req_err_type`)
    
    nit : int
        if `return_nit=True`, the number of iteration(s) done (in Gaussian-Newton algorithm),
        or the number of updates of the edge conductance estimates

    Req_estimate_arr_list : list of 1D-arrays of shape (n_edges, )
        if `return_Req_estimate_arr_list=True`, the list of arrays of edge 
        equivalent resistance estimates, at each iteration: `Req_estimate_arr_list[i][j]` 
        is the value at the `i`-th iteration for the `j`-th edge; list of length `nit + 1`

    Req_err_list : list of floats
        if `return_Req_err_list=True`, the list of errors between the edge equivalent 
        resistance estimates and the target ones, at each iteration: `Req_err_list[i]`
        is the value for the `i`-th iteration; list of length `nit + 1`

    T_estimate_arr_list : list of 1D-arrays of shape (n_edges, )
        if `return_T_estimate_arr_list=True`, the list of arrays of edge
        conductance estimates, at each iteration: `T_estimate_arr_list[i][j]` is the
        value at the `i`-th iteration for the `j`-th edge; list of length `nit + 1`

    Notes
    -----
    - the iteration `i=0` corresponds to initial (starting) values \
    (see parameter `T_estimates_start_dict`)
    - the edge `j` refers to the `j`-th edge returned by `G.edges()` (index starting from 0)
    - `n_edges = G.number_of_edges()` is the number of edges in the graph `G`
    """
    if Req_err_type not in ('MaxAE', 'MAE', 'MSE', 'RMSE', 'SAE', 'SSE'):
        raise ValueError("`Req_err_type` not valid should be 'MaxAE', 'MAE', 'MSE', 'RMSE', 'SAE' or 'SSE'")

    # Set dictionary to convert node label to node index
    node_label2index = {u:i for i, u in enumerate(G.nodes())} 
    # node_label2index = kn.utils.get_node_label2index_dict(G) # equiv.

    # Get number of nodes and number of edges
    n_nodes = G.number_of_nodes()
    n_edges = G.number_of_edges()

    # Initialize lists to be returned (if needed)
    if return_Req_estimate_arr_list:
        Req_estimate_arr_list = []
    if return_Req_err_list:
        Req_err_list = []
    if return_T_estimate_arr_list:
        T_estimate_arr_list = []

    # Set array of Req_target
    Req_target_arr = np.asarray([Req_target_dict[(u, v)] for u, v in G.edges()])
    
    # Set array and dictionary of T_estimate (start)
    if T_estimate_start_dict is not None:
        T_estimate_arr = np.asarray([T_estimate_start_dict[(u, v)] for u, v in G.edges()])
        T_estimate_dict = T_estimate_start_dict
    else:
        # Default values
        T_estimate_arr = 1.0 / np.asarray([0.5*(G.degree[u] + G.degree[v]) for u, v in G.edges()])
        T_estimate_arr = (n_nodes - 1) / T_estimate_arr.sum() * T_estimate_arr * 1.0 / Req_target_arr
        T_estimate_dict = {e: v for e, v in zip(G.edges(), T_estimate_arr)}
    
    # Set T_estimates as graph edge attribute `T_estimate_attr`
    nx.set_edge_attributes(G, T_estimate_dict, T_estimate_attr)
    
    if return_T_estimate_arr_list:
        T_estimate_arr_list.append(T_estimate_arr)

    for nit in range(nit_max+1):
        # if verbose > 0:
        #     print(f'Iteration {nit:3d} ...')

        # Compute (first) eigen values and eigen vectors of the Laplacian matrix
        # (using the current T_estimate) 
        if verbose > 1:
            if k is None:
                print(f'Iteration {nit:3d} : compute all eigen values and eigen vectors of the Laplacian matrix ...')
            else:
                print(f'Iteration {nit:3d} : compute the first {k} (at max.) eigen values and eigen vectors of the Laplacian matrix ...')
        
        eigval, eigvec = kn.tools.laplacian_eigs(G, edge_weight=T_estimate_attr, k=k)

        # Compute Req_estimate and its Jacobian matrix
        if verbose > 1:
            print(f'Iteration {nit:3d} : compute Req estimate and its Jacobian matrix ... ')
        
        Req_estimate_arr = np.zeros(n_edges)
        J_Req_estimate_arr = np.zeros((n_edges, n_edges))
        for i, e in enumerate(G.edges()):
            if verbose > 2:
                print(f'Iteration {nit:3d} : compute Req estimate and its Jacobian matrix: edge {i+1:3d} of {n_edges:3d}')
            
            i0, i1 = node_label2index[e[0]], node_label2index[e[1]]
            p = eigvec[:, 1:] @ (1.0/eigval[1:] * (eigvec[i0, 1:] - eigvec[i1, 1:]))

            Req_estimate_arr[i] = p[i0] - p[i1]
            J_Req_estimate_arr[i, :] = [- (p[node_label2index[f[0]]] - p[node_label2index[f[1]]])**2 for f in G.edges()]

        # # Alternative (equivalent)
        # # -----------
        # Req_estimate_arr = np.zeros(n_edges)
        # J_Req_estimate_arr = np.zeros((n_edges, n_edges))
        # for i, e in enumerate(G.edges()):
        #     if verbose > 2:
        #         print(f'Iteration {nit:3d} : compute Req estimate and its Jacobian matrix: edge {i+1:3d} of {n_edges:3d}')
            
        #     q = np.zeros(n_nodes)
        #     i0, i1 = node_label2index[e[0]], node_label2index[e[1]]
        #     q[i0] = 1.0
        #     q[i1] = -1.0
        #     p = (eigvec[:, 1:] * 1.0/eigval[1:]) @ eigvec[:, 1:].T @ q

        #     Req_estimate_arr[i] = p[i0] - p[i1]
        #     J_Req_estimate_arr[i, :] = [- (p[node_label2index[f[0]]] - p[node_label2index[f[1]]])**2 for f in G.edges()]

        # # Alternative (equivalent)
        # # -----------
        # Req_estimate_arr = np.zeros(n_edges)
        # J_Req_estimate_arr = np.zeros((n_edges, n_edges))
        # for i, e in enumerate(G.edges()):
        #     if verbose > 2:
        #         print(f'Iteration {nit:3d} : compute Req estimate and its Jacobian matrix: edge {i+1:3d} of {n_edges:3d}')
        #
        #     q_nodes_dict = {u: 0.0 for u in G.nodes()}
        #     q_nodes_dict[e[0]] = 1.0
        #     q_nodes_dict[e[1]] = -1.0
        #     p_nodes_dict = solve_flow_without_bc(G, q_nodes_dict=q_nodes_dict, eigval=eigval, eigvec=eigvec)

        #     Req_estimate_arr[i] = p_nodes_dict[e[0]] - p_nodes_dict[e[1]]
        #     J_Req_estimate_arr[i, :] = [- (p_nodes_dict[f[0]] - p_nodes_dict[f[1]])**2 for f in G.edges()]
        # # -----------

        # Compute error for Req
        Req_residual_arr = Req_estimate_arr - Req_target_arr
        if Req_err_type == 'MaxAE':
            Req_err = np.max(np.abs(Req_residual_arr))
        elif Req_err_type == 'SAE':
            Req_err = np.sum(np.abs(Req_residual_arr))
        elif Req_err_type == 'MAE':
            Req_err = np.mean(np.abs(Req_residual_arr))
        elif Req_err_type == 'SSE':
            Req_err = np.sum(Req_residual_arr**2)
        elif Req_err_type == 'MSE':
            Req_err = np.mean(Req_residual_arr**2)
        elif Req_err_type == 'RMSE':
            Req_err = np.sqrt(np.mean(Req_residual_arr**2))

        if return_Req_estimate_arr_list:
            Req_estimate_arr_list.append(Req_estimate_arr)

        if return_Req_err_list:
            Req_err_list.append(Req_err)

        if verbose > 0:
            print(f'Iteration {nit:3d} : Req error ({Req_err_type}) = {Req_err:4e}')

        # Check convergence
        if Req_err < tol:
            if verbose > 0:
                print(f'Iteration {nit:3d} : convergence OK (Req error less than {tol})')
            break

        if nit > 0 and np.abs(Req_err - Req_err_prev) < Req_err_prev*tol_step:
            if verbose > 0:
                print(f'Iteration {nit:3d} : convergence OK (rel. abs. diff. between two last Req err less than {tol_step})')
            break

        if nit < nit_max:
            # Compute delta_T for next estimate
            if verbose > 1:
                print(f'Iteration {nit:3d} : compute update (delta) for T estimate')
           
            delta_T_arr = np.linalg.solve(J_Req_estimate_arr.T @ J_Req_estimate_arr, - J_Req_estimate_arr.T @ Req_residual_arr)

            # Set next estimate for T
            T_estimate_arr = T_estimate_arr + delta_T_arr
            T_estimate_arr = np.clip(T_estimate_arr, a_min=T_clip_min, a_max=None)  # avoid T value <= 0
            if rescale_T: 
                T_estimate_arr = T_estimate_arr / np.sum(T_estimate_arr*Req_target_arr) * (n_nodes - 1) # rescale T

            # Set T_estimates as graph edge attribute `T_estimate_attr`
            T_estimate_dict = {e: v for e, v in zip(G.edges(), T_estimate_arr)}
            nx.set_edge_attributes(G, T_estimate_dict, T_estimate_attr)

            if return_T_estimate_arr_list:
                T_estimate_arr_list.append(T_estimate_arr)

            # Save Req err for next iteration
            Req_err_prev = Req_err

    out = []
    if return_Req_estimate_dict:
        Req_estimate_dict = {e: v for e, v in zip(G.edges(), Req_estimate_arr)}
        out.append(Req_estimate_dict)

    if return_Req_err:
        out.append(Req_err)

    if return_nit:
        out.append(nit)

    if return_Req_estimate_arr_list:
        out.append(Req_estimate_arr_list)

    if return_Req_err_list:
        out.append(Req_err_list)

    if return_T_estimate_arr_list:
        out.append(T_estimate_arr_list)

    if len(out) == 0:
        return None

    if len(out) == 1:
        out = out[0]
    else:
        out = tuple(out)
    
    return out
# ----------------------------------------------------------------------------

# =============================================================================
# "Private functions"
# =============================================================================

# ----------------------------------------------------------------------------
def _loglog_fit(
        x_all, 
        y_all, 
        data_points_removed_first_fraction,
        data_points_removed_last_fraction,
        verbose=0):
    """Internal (private) function to fit line on a log log plot.
    """
    n_all = len(x_all)
    
    # Coordinates used for fitting
    n1 = int(data_points_removed_first_fraction * n_all)
    n2 = n_all - int(data_points_removed_last_fraction * n_all)
    # n1 = int(np.round(data_points_removed_first_fraction * n_all))
    # n2 = n_all - int(np.round(data_points_removed_last_fraction * n_all))

    x = x_all[n1:n2]
    y = y_all[n1:n2]

    # ... extract non nan number
    ind = ~np.any((np.isnan(x), np.isnan(y)), axis=0)
    if not ind.all():
        if verbose > 0:
            print(f'WARNING: undefined point (nan) encountered (removed)')

        x = x[ind]
        y = y[ind]

    n_fit = len(x) # number of data points used for fitting

    if n_fit < 2:
        # No fit
        if verbose > 0:
            print(f'WARNING: not enough points ({n_fit}) for line fitting in log-log plot ; try to increase `fit_last_data_points_fraction`"')

        poly_fit_params, poly_fit_params_cov = None, None

    elif n_fit == 2:
        # Fit without uncertainty  
        if verbose > 0:
            print(f'WARNING: fitting on 2 points in log-log plot ; try to increase `fit_last_data_points_fraction`"')
        
        try:
            poly_fit_params = np.polyfit(np.log10(x), np.log10(y), 1, cov=False) # cov=False, i.e. no uncertainty
            poly_fit_params_cov = None

        except:
            if verbose > 0:
                print(f'WARNING: fitting failed [on {n_fit} points (over {n_all})]')
            
            poly_fit_params, poly_fit_params_cov = None, None

    else:
        # Poly-fit
        if verbose > 0:
            print(f'Fitting on {n_fit} points (over {n_all})')

        try:
            poly_fit_params, poly_fit_params_cov = np.polyfit(np.log10(x), np.log10(y), 1, cov=True)

        except:
            if verbose > 0:
                print(f'WARNING: fitting failed [on {n_fit} points (over {n_all})]')
            
            poly_fit_params, poly_fit_params_cov = None, None

    return poly_fit_params, poly_fit_params_cov, x, y, n1, n2
# ----------------------------------------------------------------------------

# #################### OLD BELOW ##############################################

# # ----------------------------------------------------------------------------
# def krige_node_attribute(
#         G,
#         cov_model,
#         node_attr,
#         err_std=0.0,
#         edge_length_attr=None,
#         method='ordinary_kriging',
#         mean=None,
#         var=None,
#         use_unique_neighborhood=False,
#         searchRadius=None,
#         searchRadiusRelative=1.2,
#         nneighborMax=12,
#         update_graph=True,
#         kriging_est_node_attr=None,
#         kriging_std_node_attr=None,
#         i0=None, 
#         i1=None,
#         pid=None,
#         verbose=1):
#     """
#     Interpolates a node attribute on a graph by kriging.

#     Parameters
#     ----------
#     G : networkx.Graph
#         input graph

#     cov_model : :class:`geone.covModel.CovModel1D`
#         covariance model in 1D

#     node_attr : str
#         name of the node attribute to be kriged

#     err_std : float, or str, default : 0.0
#         standard deviation of error (zero-mean Gaussian of given std):
#         - if float: the same error std is used for data at any graph node; 
#         - if str: name of the node attribute of the error std

#     edge_length_attr : str, optional
#         name of the edge attribute for length (used to compute the distance between
#         nodes);
#         by default (`None`): the edges have a length of one

#     method : str {'simple_kriging', 'ordinary_kriging'}, default: 'ordinary_kriging'
#         type of kriging;
#         note: if `method='ordinary_kriging'`, the parameter `mean` is not used

#     mean : float, or str, optional
#         kriging mean value :
#         - if float: the value is used at all graph nodes; 
#         - if str: name of the node attribute for the mean 
        
#         by default (`None`), the mean of data node values is used for any node
        
#         note: if `method=ordinary_kriging`, parameter `mean` is ignored

#     var : float, or str, optional
#         kriging variance value :
#         - if float: the value is used at all graph nodes; 
#         - if str: name of the node attribute for the mean 

#         note: if `method=ordinary_kriging`, parameter `var` is ignored

#     use_unique_neighborhood : bool, default: False
#         indicates if a unique neighborhood is used:

#         - if True: data at any graph node is taken into account, and the kriging matrix \
#         is computed once; the parameters `searchRadius`, `searchRadiusRelative`, \
#         `nneighborMax` are not used
#         - if False: only data at graph node within a search neighborhood are \
#         taken into account according to `searchRadius`, `searchRadiusRelative`, `nneighborMax`

#     searchRadius : float, optional
#         if specified, i.e. not `None`: search radius, i.e. 
#         the data at graph node at distance to the estimated graph node greater 
#         than `searchRadius` are not taken into account in the kriging system; 
#         if `searchRadius` is not `None`, then `searchRadiusRelative` is not used;
#         by default (`searchRadius=None`): `searchRadiusRelative` is used;

#     searchRadiusRelative : float, default: 1.2
#         used only if `searchRadius` is `None`;
#         the search radius is set to `searchRadiusRelative` times the range of the 
#         covariance model `cov_model`

#     nneighborMax : int, default: 12
#         maximal number of neighbors (data at graph nodes) taken into account in the
#         kriging system; the data at graph nodes the closest to the estimated graph node are
#         taken into account;
#         note: if `nneighborMax=None` or `nneighborMax<0`, then `nneighborMax` is
#         set to the number of graph nodes with informed data

#     update_graph : bool, default True
#         - if `True`: the kriging estimates is set as a node attribute (named \
#         `kriging_est_node_attr` (see below)) and the kriging standard deviation \
#         is set as a node attribute (named `kriging_std_node_attr` (see below))
#         - if `False`: the graph is not updated (`kriging_est_node_attr` and 
#         `kriging_std_node_attr` are not used)

#     kriging_est_node_attr : str, optional
#         name of the node attribute in output for kriging estimates;
#         by default (`None`),  the name `node_attr` + '_krig_est' is used

#     kriging_std_node_attr : str, optional
#         name of the node attribute in output for kriging standard deviation;
#         by default (`None`),  the name `node_attr` + '_krig_std' is used

#     i0 : int, optional
#         starting index in `list(G.nodes())` of nodes to be kriged (used with multiprocessing)

#     i1 : int, optional
#         ending index (excluded) in `list(G.nodes())` of nodes to be kriged (used with multiprocessing)

#     pid : int, optional
#         process id of the caller (used with multiprocessing)

#     verbose : int, default: 0
#         verbose mode, higher implies more printing (info)

#     Returns
#     -------
#     kriging_est : dict
#         dictionary with node labels as keys and kriging estimates as values

#     kriging_std : dict
#         dictionary with node labels as keys and kriging standard deviation as values
#     """
#     if verbose > 0:
#         if pid is not None:
#             pid_str = f'[pid={pid}] '
#         else:
#             pid_str = ''

#     # Check cov_model
#     if not isinstance(cov_model, geone.covModel.CovModel1D):
#         raise ValueError(f'{pid_str}`cov_model` must be an instance of `geone.covModel.CovModel1D`')

#     # Prevent calculation if covariance model is not stationary
#     if not cov_model.is_stationary():
#         raise ValueError(f'{pid_str}`cov_model` is not stationary')

#     # Check node attribute
#     if node_attr not in kn.utils.get_node_attribute_names(G):
#         raise ValueError(f'{pid_str}{node_attr} is not a node attribute')

#     # Check edge length attribute
#     if edge_length_attr is not None and edge_length_attr not in kn.utils.get_edge_attribute_names(G):
#         raise ValueError(f'{pid_str}{edge_length_attr} is not an edge attribute')

#     # # Set dictionary to convert node label (id) to node index, and vice versa
#     # node_label2index = {u:i for i, u in enumerate(G.nodes())}
#     # #node_index2label = {i:u for i, u in enumerate(G.nodes())}

#     if update_graph:
#         # Set node attribute names for output
#         if kriging_est_node_attr is not None:
#             if not isinstance(kriging_est_node_attr, str):
#                 raise ValueError(f'{pid_str}`kriging_est_node_attr` must be a string (node attribute name)')
#         else:
#             kriging_est_node_attr = f'{node_attr}_krig_est'

#         if kriging_std_node_attr is not None:
#             if not isinstance(kriging_std_node_attr, str):
#                 raise ValueError(f'{pid_str}`kriging_std_node_attr` must be a string (node attribute name)')
#         else:
#             kriging_std_node_attr = f'{node_attr}_krig_std'

#     # Get the dictionary of node attribute (property) to be kriged
#     prop_dict = nx.get_node_attributes(G, node_attr)
#     if len(prop_dict) == 0:
#         raise ValueError(f'{pid_str}No value for specified node attribute ({node_attr})')

#     v_first = np.asarray(list(prop_dict.values())[0]) # first value entry as array
#     if v_first.ndim > 0:
#         raise ValueError(f'{pid_str}Value of the specified attribute should by a scalar')

#     # - set nan at uninformed nodes
#     prop_dict = nx.get_node_attributes(G, node_attr, default=np.nan)

#     # Covariance function and value at 0
#     cov_func = cov_model.func() # covariance function
#     cov0 = cov_func(0.)[0] # covariance function at origin (lag=0)

#     # Default mean value    
#     tmp = np.asarray(list(prop_dict.values()))
#     if np.isnan(tmp).all():
#         mean_default = 0.0
#     else:
#         # Set mean of data values
#         mean_default = np.nanmean(tmp)
    
#     # Method and mean, var
#     if method == 'simple_kriging':
#         ordinary_kriging = False
#         # Get the dictionary of mean kriging values
#         if mean is not None:
#             if isinstance(mean, float) or isinstance(mean, int):
#                 mean_dict = {u:mean for u in G.nodes()}

#             elif isinstance(mean, str):
#                 mean_dict = nx.get_node_attributes(G, mean, default=np.nan)
#                 if np.isnan(np.asarray(list(mean_dict.values()))).any():
#                     raise ValueError(f'{pid_str}Specified mean is not defined at all graph nodes')
            
#             else:
#                 raise ValueError(f'{pid_str}Specified mean is not valid')

#         else:
#             mean_dict = {u:mean_default for u in G.nodes()}
    
#         # Get the dictionary of "variance update"
#         if var is not None:
#             if isinstance(var, float) or isinstance(var, int):
#                 tmp = np.sqrt(var/cov0)
#                 var_update_dict = {u:tmp for u in G.nodes()}

#             elif isinstance(var, str):
#                 var_update_dict = nx.get_node_attributes(G, var, default=np.nan)
#                 if np.isnan(np.asarray(list(var_update_dict.values()))).any():
#                     raise ValueError(f'{pid_str}Specified var is not defined at all graph nodes')

#                 var_update_dict = {k: np.sqrt(v/cov0) for k, v in var_update_dict.items()}

#             else:
#                 raise ValueError(f'{pid_str}Specified var is not valid')

#     elif method == 'ordinary_kriging':
#         ordinary_kriging = True
#         if verbose > 0 and mean is not None:
#             print(f"{pid_str}WARNING: `mean` is ignored with `method='ordinary_kriging'`")

#         mean = None

#         if verbose > 0 and var is not None:
#             print(f"{pid_str}WARNING: `var` is ignored with `method='ordinary_kriging'`")

#         var = None

#     else:
#         raise ValueError(f'{pid_str}`method` ({method}) unknown')

#     # List of data node labels / index
#     data_node_labels = [u for u, v in prop_dict.items() if not np.isnan(v)]

#     # Number of data nodes 
#     n = len(data_node_labels)

#     # Number of nodes in the graph
#     n_nodes = G.number_of_nodes()
#     if i0 is None:
#         i0 = 0
#     if i1 is None:
#         i1 = n_nodes

#     node_list = list(G.nodes())[i0:i1]
#     node_list_len = i1 - i0

#     if n == 0:
#         if ordinary_kriging:
#             krig_est_dict = {u: 0.0 for u in node_list}

#             tmp = float(np.sqrt(cov0))
#             krig_std_dict = {u: tmp for u in node_list}

#         else: # simple kriging
#             krig_est_dict = {u: float(mean_dict[u]) for u in node_list}

#             tmp = float(np.sqrt(cov0))
#             if var is not None:
#                 krig_std_dict = {u: tmp*var_update_dict[u] for u in node_list}
#             else:
#                 krig_std_dict = {u: tmp for u in node_list}

#         if update_graph:
#             # Set node attributes (output)
#             nx.set_node_attributes(G, krig_est_dict, kriging_est_node_attr)
#             nx.set_node_attributes(G, krig_std_dict, kriging_std_node_attr)

#         return krig_est_dict, krig_std_dict
    
#     # Get the dictionary of error variance
#     if err_std is None:
#         err_std = 0.0

#     if isinstance(err_std, float) or isinstance(err_std, int):
#         err_var = err_std**2
#         err_var_dict = {u:err_var for u in data_node_labels}

#     elif isinstance(err_std, str):
#         err_var_dict = nx.get_node_attributes(G, err_std)
#         if np.any(np.asarray([u not in err_var_dict.keys() for u in data_node_labels])):
#             raise ValueError(f'{pid_str}Specified error std is must be defined at all data nodes')
        
#         err_var_dict = {u:s**2 for u, s in err_var_dict.items()}

#     # Do kriging of nodes in node_list
#     if use_unique_neighborhood:
#         # Initialize 
#         # - kriging matrix (mat) of order nmat
#         # - right hand side of all kriging systems (b), matrix of dimension nmat x node_list_len
#         if ordinary_kriging:
#             nmat = n+1
#             mat = np.ones((nmat, nmat))
#             mat[-1,-1] = 0.0
#         else:
#             nmat = n
#             mat = np.ones((nmat, nmat))

#         b = np.ones((nmat, node_list_len))

#         # Set kriging matrix (mat) and right hand side of all kriging systems (b)
#         for i, ui in enumerate(data_node_labels):
#             length_dict = nx.single_source_dijkstra_path_length(G, ui, cutoff=None, weight=edge_length_attr)
#             h = np.asarray([length_dict[u] for u in data_node_labels])
#             mat[i, :n] = cov_func(h)
#             h = np.asarray([length_dict[u] for u in node_list])
#             b[i, :] = cov_func(h)
        
#         # Add error variance on the diagonal of the kriging matrix
#         for i, u in enumerate(data_node_labels):
#             mat[i, i] = mat[i, i] + err_var_dict[u]

#         # Solve all kriging systems
#         w = np.linalg.solve(mat, b) # w: matrix of dimension nmat x n_nodes

#         # Kriged values
#         if mean is not None:
#             # simple kriging
#             if var is not None:
#                 tmp = np.asarray([1.0/var_update_dict[ui] * (prop_dict[ui] - mean_dict[ui]) for ui in data_node_labels]).dot(w)
#                 krig_est_dict = {u: float(mean_dict[u] + var_update_dict[u]*v) for u, v in zip(node_list, tmp)}
#             else:
#                 tmp = np.asarray([prop_dict[ui] - mean_dict[ui] for ui in data_node_labels]).dot(w)
#                 krig_est_dict = {u: float(mean_dict[u] + v) for u, v in zip(node_list, tmp)}
#         else:
#             # ordinary kriging
#             tmp = np.asarray([prop_dict[ui] for ui in data_node_labels]).dot(w[:n, :])
#             krig_est_dict = {u: float(v) for u, v in zip(node_list, tmp)}

#         # Kriged standard deviation
#         tmp = np.sqrt(np.maximum(0.0, cov0 - np.array([np.dot(w[:,i], b[:,i]) for i in range(node_list_len)])))
#         krig_std_dict = {u: float(v) for u, v in zip(node_list, tmp)}

#     else:
#         # Limited search neighborhood

#         # Set dmax (search radius)
#         if searchRadius is not None:
#             if searchRadius <= 0.0:
#                 raise ValueError(f'{pid_str}`searchRadius` not valid (negative)')

#             dmax = searchRadius

#         else:
#             # use searchRadiusRelative
#             if searchRadiusRelative <= 0.0:
#                 raise ValueError(f'{pid_str}`searchRadiusRelative` (factor) not valid (negative)')
            
#             dmax = searchRadiusRelative * cov_model.r()

#         if nneighborMax is None or nneighborMax > n or nneighborMax < 0:
#             nneighborMax = n

#         # Initialize kriging matrix, second member and property values at neigbhors
#         mat = np.ones((nneighborMax+1, nneighborMax+1))
#         b = np.ones(nneighborMax+1) 
#         prop_val = np.ones(nneighborMax)

#         # Initialize the dictionaries of kriging estimates and standard deviation
#         krig_est_dict = {u:np.nan for u in node_list}
#         krig_std_dict = {u:np.nan for u in node_list}

#         if verbose > 0:
#             progress_old = 0

#         for k, u in enumerate(node_list):
#             if verbose > 0:
#                 progress = int(k/node_list_len*100.0)
#                 if progress > progress_old:
#                     print(f'{pid_str}Kriging: {progress:3d}%')
#                     progress_old = progress
            
#             if u in data_node_labels and err_var_dict[u] == 0.0:
#                 krig_est_dict[u] = prop_dict[u]
#                 krig_std_dict[u] = 0.0
#                 continue

#             length_dict = nx.single_source_dijkstra_path_length(G, u, cutoff=dmax, weight=edge_length_attr)
#             length_keys_list = list(length_dict.keys())
#             ind = np.argsort(np.asarray(list(length_dict.values())))
#             neighbor_labels = []
#             nn = 0
#             for i in ind:
#                 ui = length_keys_list[i]
#                 if not np.isnan(prop_dict[ui]):
#                     prop_val[nn] = prop_dict[ui]
#                     b[nn] = cov_model(length_dict[ui])[0]
#                     neighbor_labels.append(ui)
#                     nn = nn+1
#                     if nn >= nneighborMax:
#                         break

#             if nn == 0:
#                 if mean is not None:
#                     # simple kriging
#                     krig_est_dict[u] = float(mean_dict[u])
#                 else:
#                     # ordinary kriging
#                     krig_est_dict[u] = mean_default

#                 krig_std_dict[u] = float(np.sqrt(cov0))
#                 continue
                
#             for i in range(nn-1):
#                 ui = neighbor_labels[i]
#                 for j in range(i+1, nn):
#                     uj = neighbor_labels[j]
#                     h = nx.dijkstra_path_length(G, source=ui, target=uj, weight=edge_length_attr)
#                     # h = nx.shortest_path_length(G, source=ui, target=uj, weight=edge_length_attr, method='dijkstra')
#                     cov_h = cov_func(h)[0]
#                     mat[i, j] = cov_h
#                     mat[j, i] = cov_h
                
#                 mat[i, i] = cov0 + err_var_dict[ui]
            
#             mat[nn-1, nn-1] = cov0 + err_var_dict[neighbor_labels[nn-1]]

#             if ordinary_kriging:
#                 nmat = nn+1
#                 mat[nn, :] = 1.0
#                 mat[:, nn] = 1.0
#                 mat[nn, nn] = 0.0
#                 b[nn] = 1.0
            
#             else:
#                 nmat = nn

#             # Solve kriging system
#             w = np.linalg.solve(mat[:nmat, :nmat], b[:nmat]) # w: vector of dimension nmat

#             if mean is not None:
#                 # simple kriging
#                 if var is not None:
#                     krig_est_dict[u] = float(mean_dict[u] + var_update_dict[u] * (np.asarray([1.0/var_update_dict[ui] for ui in neighbor_labels])*(prop_val[:nn] - np.asarray([mean_dict[ui] for ui in neighbor_labels]))).dot(w))
#                 else:
#                     krig_est_dict[u] = float(mean_dict[u] + (prop_val[:nn] - np.asarray([mean_dict[ui] for ui in neighbor_labels])).dot(w))
#             else:
#                 # ordinary kriging
#                 krig_est_dict[u] = float(prop_val[:nn].dot(w[:nn]))

#             # Kriged standard deviation
#             krig_std_dict[u] = float(np.sqrt(np.maximum(0.0, cov0 - b[:nmat].dot(w))))
        
#         if verbose > 0:
#             progress = 100
#             if progress > progress_old:
#                 print(f'{pid_str}Kriging: {progress:3d}%')
#                 progress_old = progress

#     if var is not None:
#         krig_std_dict = {u: var_update_dict[u]*v for u, v in krig_std_dict.items()}

#     if update_graph:
#         # Set node attributes (output)
#         nx.set_node_attributes(G, krig_est_dict, kriging_est_node_attr)
#         nx.set_node_attributes(G, krig_std_dict, kriging_std_node_attr)
    
#     return krig_est_dict, krig_std_dict
# # ----------------------------------------------------------------------------

# # ----------------------------------------------------------------------------
# def krige_node_attribute_mp(
#         G,
#         cov_model,
#         node_attr,
#         err_std=0.0,
#         edge_length_attr=None,
#         method='ordinary_kriging',
#         mean=None,
#         var=None,
#         use_unique_neighborhood=False,
#         searchRadius=None,
#         searchRadiusRelative=1.2,
#         nneighborMax=12,
#         update_graph=True,
#         kriging_est_node_attr=None,
#         kriging_std_node_attr=None,
#         verbose=1,
#         nproc=-1):
#     """
#     Computes the same as the function :func:`krige_node_attribute`, using multiprocessing.

#     All the parameters except `nproc` are the same as those of the function
#     :func:`krige_node_attribute`.

#     This function launches parallel processes [parallel calls of the
#     function :func:`krige_node_attribute`]; the set of nodes to be kriged is distributed in a 
#     balanced way over the processes.

#     The number of processes used (in parallel) is determined by the parameter `nproc` 
#     (int, default: -1); a negative number (or zero), -n <= 0, can be specified 
#     to use the total number of cpu(s) of the system except n; `nproc` is finally
#     at maximum equal to `G.number_of_nodes()` but at least 1 by applying:
        
#     - if `nproc >= 1`, then `nproc = max(min(nproc, G.number_of_nodes()), 1)` is used
#     - if `nproc = -n <= 0`, then `nproc = max(min(nmax-n, G.number_of_nodes()), 1)` is used, \
#     where nmax is the total number of cpu(s) of the system (retrieved by 
#     `multiprocessing.cpu_count()`)

#     Note: if `nproc=None`, `nproc=-1` is used.

#     See function :func:`krige_node_attribute` for details.
#     """
#     if update_graph:
#         # Set node attribute names for output
#         if kriging_est_node_attr is not None:
#             if not isinstance(kriging_est_node_attr, str):
#                 raise ValueError('`kriging_est_node_attr` must be a string (node attribute name)')
#         else:
#             kriging_est_node_attr = f'{node_attr}_krig_est'

#         if kriging_std_node_attr is not None:
#             if not isinstance(kriging_std_node_attr, str):
#                 raise ValueError('`kriging_std_node_attr` must be a string (node attribute name)')
#         else:
#             kriging_std_node_attr = f'{node_attr}_krig_std'

#     # Number of graph nodes
#     n_nodes = G.number_of_nodes()

#     # Set number of process(es): nproc
#     if nproc is None:
#         nproc = -1
    
#     if nproc <= 0:
#         nproc = max(min(multiprocessing.cpu_count() + nproc, n_nodes), 1)
#     else:
#         nproc_tmp = nproc
#         nproc = max(min(int(nproc), n_nodes), 1)
#         if verbose > 0 and nproc != nproc_tmp:
#             print(f'Number of processes has been changed (now: nproc={nproc})')

#     # Set index for distributing realizations
#     q, r = np.divmod(n_nodes, nproc)
#     ids_proc = [i*q + min(i, r) for i in range(nproc+1)]

#     if verbose > 0:
#         print(f'Running `krige_node_attribute` on {nproc} processes...')

#     # Set pool of nproc workers
#     pool = multiprocessing.Pool(nproc)
#     out_pool = []
#     for i in range(nproc):
#         # Set i-th process
#         args = (G, cov_model, node_attr)
#         kwargs = dict(
#                     err_std=err_std,
#                     edge_length_attr=edge_length_attr, 
#                     method=method,
#                     mean=mean,
#                     var=var,
#                     use_unique_neighborhood=use_unique_neighborhood,
#                     searchRadius=searchRadius, 
#                     searchRadiusRelative=searchRadiusRelative, 
#                     nneighborMax=nneighborMax,
#                     update_graph=False, # do not update graph if multiprocessing is enabled
#                     kriging_est_node_attr=kriging_est_node_attr,
#                     kriging_std_node_attr=kriging_std_node_attr,
#                     i0=ids_proc[i],
#                     i1=ids_proc[i+1], 
#                     pid=i,
#                     verbose=verbose*(i==0))
#         out_pool.append(pool.apply_async(krige_node_attribute, args=args, kwds=kwargs))

#     # Properly end working process
#     pool.close() # Prevents any more tasks from being submitted to the pool,
#     pool.join()  # then, wait for the worker processes to exit.

#     # Get result from each process
#     out = [w.get() for w in out_pool]
#     if np.any([x is None for x in out]):
#         err_msg = f'`krige_node_attribute_mp`: an error occured on a process (worker)'
#         raise ValueError(err_msg)

#     # Gather dictionaries of all processes
#     krig_est_dict = {}
#     krig_std_dict = {}
#     for krig_est_dict_pid_i, krig_std_dict_pid_i in out:
#         krig_est_dict.update(krig_est_dict_pid_i)
#         krig_std_dict.update(krig_std_dict_pid_i)

#     if update_graph:
#         # Set node attributes (output)
#         nx.set_node_attributes(G, krig_est_dict, kriging_est_node_attr)
#         nx.set_node_attributes(G, krig_std_dict, kriging_std_node_attr)
    
#     return krig_est_dict, krig_std_dict
# # ----------------------------------------------------------------------------

# # ----------------------------------------------------------------------------
# def sgs_node_attribute(
#         G,
#         cov_model,
#         node_attr,
#         err_std=0.0,
#         edge_length_attr=None,
#         method='ordinary_kriging',
#         mean=None,
#         var=None,
#         searchRadius=None,
#         searchRadiusRelative=1.2,
#         nneighborMax=12,
#         nreal=1,
#         seed=None,
#         update_graph=True,
#         sgs_node_attr=None,
#         pid=None,
#         verbose=1):
#     """
#     Performs Sequential Gaussian Simulation (SGS) of a node attribute on a graph.

#     Parameters
#     ----------
#     G : networkx.Graph
#         input graph

#     cov_model : :class:`geone.covModel.CovModel1D`
#         covariance model in 1D

#     node_attr : str
#         name of the node attribute to be simulated

#     err_std : float, or str, default : 0.0
#         standard deviation of error (zero-mean Gaussian of given std):
#         - if float: the same error std is used for data at any graph node; 
#         - if str: name of the node attribute of the error std

#     edge_length_attr : str, optional
#         name of the edge attribute for length (used to compute the distance between
#         nodes);
#         by default (`None`): the edges have a length of one

#     method : str {'simple_kriging', 'ordinary_kriging'}, default: 'ordinary_kriging'
#         type of kriging;
#         note: if `method='ordinary_kriging'`, the parameter `mean` is not used

#     mean : float, or str, default : 0.0
#         kriging mean value :
#         - if float: the value is used at all graph nodes; 
#         - if str: name of the node attribute for the mean 
        
#         note: if `method=ordinary_kriging`, parameter `mean` is ignored

#     var : float, or str, optional
#         kriging variance value :
#         - if float: the value is used at all graph nodes; 
#         - if str: name of the node attribute for the mean 

#         note: if `method=ordinary_kriging`, parameter `var` is ignored
        
#     searchRadius : float, optional
#         if specified, i.e. not `None`: search radius, i.e. 
#         the data at graph node at distance to the estimated graph node greater 
#         than `searchRadius` are not taken into account in the kriging system; 
#         if `searchRadius` is not `None`, then `searchRadiusRelative` is not used;
#         by default (`searchRadius=None`): `searchRadiusRelative` is used;

#     searchRadiusRelative : float, default: 1.2
#         used only if `searchRadius` is `None`;
#         the search radius is set to `searchRadiusRelative` times the range of the 
#         covariance model `cov_model`

#     nneighborMax : int, default: 12
#         maximal number of neighbors (data at graph nodes) taken into account in the
#         kriging system; the data at graph nodes the closest to the estimated graph node are
#         taken into account;
#         note: if `nneighborMax=None` or `nneighborMax<0`, then `nneighborMax` is
#         set to the number of graph nodes with informed data

#     nreal : int, default: 1
#         number of realization(s)

#     seed : int, optional
#         seed for initializing random number generator

#     update_graph : bool, default True
#         - if `True`: the simulations (realizations) are set as a node attribute (named \
#         `sgs_node_attr` (see below))
#         - if `False`: the graph is not updated (`sgs_node_attr` is not used)

#     sgs_node_attr : str, optional
#         name of the node attribute in output for simulation;
#         by default (`None`),  the name `node_attr` + '_sgs' is used;
#         the values of the this output node attributes are lists of length `nreal`
#         containing all the realizations

#     pid : int, optional
#         process id of the caller (used with multiprocessing)
    
#     verbose : int, default: 0
#         verbose mode, higher implies more printing (info)

#     Returns
#     -------
#     sgs : dict
#         dictionary with node labels as keys and simulations (realizations) as values;
#         each value is a list of length `nreal` containing all the realizations
#     """
#     if verbose > 0:
#         if pid is not None:
#             pid_str = f'[pid={pid}] '
#         else:
#             pid_str = ''

#     # Check cov_model
#     if not isinstance(cov_model, geone.covModel.CovModel1D):
#         raise ValueError(f'{pid_str}`cov_model` must be an instance of `geone.covModel.CovModel1D`')

#     # Prevent calculation if covariance model is not stationary
#     if not cov_model.is_stationary():
#         raise ValueError(f'{pid_str}`cov_model` is not stationary')

#     # Check node attribute
#     if node_attr not in kn.utils.get_node_attribute_names(G):
#         raise ValueError(f'{pid_str}{node_attr} is not a node attribute')

#     # Check edge length attribute
#     if edge_length_attr is not None and edge_length_attr not in kn.utils.get_edge_attribute_names(G):
#         raise ValueError(f'{pid_str}{edge_length_attr} is not an edge attribute')

#     # # Set dictionary to convert node label (id) to node index, and vice versa
#     # node_label2index = {u:i for i, u in enumerate(G.nodes())}
#     # #node_index2label = {i:u for i, u in enumerate(G.nodes())}

#     if update_graph:
#         # Set node attribute name for output
#         if sgs_node_attr is not None:
#             if not isinstance(sgs_node_attr, str):
#                 raise ValueError(f'{pid_str}`sgs_node_attr` must be a string (node attribute name)')
#         else:
#             sgs_node_attr = f'{node_attr}_sgs'

#     # Get the dictionary of node attribute (property) to be kriged
#     prop_dict = nx.get_node_attributes(G, node_attr)
#     if len(prop_dict) == 0:
#         raise ValueError(f'{pid_str}No value for specified node attribute ({node_attr})')

#     v_first = np.asarray(list(prop_dict.values())[0]) # first value entry as array
#     if v_first.ndim > 0:
#         raise ValueError(f'{pid_str}Value of the specified attribute should by a scalar')

#     # - set nan at uninformed nodes
#     prop_dict = nx.get_node_attributes(G, node_attr, default=np.nan)

#     # Covariance function and value at 0
#     cov_func = cov_model.func() # covariance function
#     cov0 = cov_func(0.)[0] # covariance function at origin (lag=0)

#     # Default mean value    
#     tmp = np.asarray(list(prop_dict.values()))
#     if np.isnan(tmp).all():
#         mean_default = 0.0
#     else:
#         # Set mean of data values
#         mean_default = np.nanmean(tmp)
    
#     # Method and mean, var
#     if method == 'simple_kriging':
#         ordinary_kriging = False
#         # Get the dictionary of mean kriging values
#         if mean is not None:
#             if isinstance(mean, float) or isinstance(mean, int):
#                 mean_dict = {u:mean for u in G.nodes()}

#             elif isinstance(mean, str):
#                 mean_dict = nx.get_node_attributes(G, mean, default=np.nan)
#                 if np.isnan(np.asarray(list(mean_dict.values()))).any():
#                     raise ValueError(f'{pid_str}Specified mean is not defined at all graph nodes')

#         else:
#             mean_dict = {u:mean_default for u in G.nodes()}
    
#         # Get the dictionary of "variance update"
#         if var is not None:
#             if isinstance(var, float) or isinstance(var, int):
#                 tmp = np.sqrt(var/cov0)
#                 var_update_dict = {u:tmp for u in G.nodes()}

#             elif isinstance(var, str):
#                 var_update_dict = nx.get_node_attributes(G, var, default=np.nan)
#                 if np.isnan(np.asarray(list(var_update_dict.values()))).any():
#                     raise ValueError(f'{pid_str}Specified var is not defined at all graph nodes')

#                 var_update_dict = {k: np.sqrt(v/cov0) for k, v in var_update_dict.items()}

#             else:
#                 raise ValueError(f'{pid_str}Specified var is not valid')

#     elif method == 'ordinary_kriging':
#         ordinary_kriging = True
#         if verbose > 0 and mean is not None:
#             print(f"{pid_str}WARNING: `mean` is ignored with `method='ordinary_kriging'`")

#         mean = None

#         if verbose > 0 and var is not None:
#             print(f"{pid_str}WARNING: `var` is ignored with `method='ordinary_kriging'`")

#         var = None

#     else:
#         raise ValueError(f'{pid_str}`method` ({method}) unknown')

#     # List of data node labels / index
#     data_node_labels = [u for u, v in prop_dict.items() if not np.isnan(v)]

#     # Number of data nodes 
#     n = len(data_node_labels)

#     # Get the dictionary of error variance
#     # - initialize
#     err_var_dict = {u:0.0 for u in G.nodes()}
    
#     if err_std is None:
#         err_std = 0.0

#     if isinstance(err_std, float) or isinstance(err_std, int):
#         err_var = err_std**2
#         for u in data_node_labels:
#             err_var_dict[u] = err_var 
            
#     elif isinstance(err_std, str):
#         tmp = nx.get_node_attributes(G, err_std)
#         if np.any(np.asarray([u not in tmp.keys() for u in data_node_labels])):
#             raise ValueError(f'{pid_str}Specified error std is must be defined at all data nodes')
        
#         for u in data_node_labels:
#             err_var_dict[u] = tmp[u]**2 

#     # Number of nodes in the graph
#     n_nodes = G.number_of_nodes()

#     # Limited search neighborhood

#     # Set dmax (search radius)
#     if searchRadius is not None:
#         if searchRadius <= 0.0:
#             raise ValueError(f'{pid_str}`searchRadius` not valid (negative)')

#         dmax = searchRadius

#     else:
#         # use searchRadiusRelative
#         if searchRadiusRelative <= 0.0:
#             raise ValueError(f'{pid_str}`searchRadiusRelative` (factor) not valid (negative)')
        
#         dmax = searchRadiusRelative * cov_model.r()

#     if nneighborMax is None or nneighborMax > n or nneighborMax < 0:
#         nneighborMax = n

#     # Initialize kriging matrix, second member and property values at neigbhors
#     mat = np.ones((nneighborMax+1, nneighborMax+1))
#     b = np.ones(nneighborMax+1) 
#     prop_val = np.ones(nneighborMax)

#     # Initialize the dictionary of sgs
#     sgs_all_dict = {u:[] for u in G.nodes()}

#     if seed is None:
#         seed = np.random.randint(1, 1000000)
#     seed = int(seed)

#     if verbose > 0:
#         progress_old = 0

#     for ireal in range(nreal):
#         # Initialize random number generator
#         np.random.seed(seed+ireal)
#         # sel_dict = {u:False for u in G.nodes()}
#         # for u in data_node_labels:
#         #     sel_dict[u] = True
#         err_var_curr_dict = err_var_dict.copy()
#         sgs_dict = {u:np.nan for u in G.nodes()}
#         for u in data_node_labels:
#             sgs_dict[u] = prop_dict[u]

#         for k, u in enumerate(G.nodes()):
#             if verbose > 0:
#                 progress = int((k+ireal*n_nodes)/(nreal*n_nodes)*100.0)
#                 if progress > progress_old:
#                     print(f'{pid_str}SGS: {progress:3d}% ({ireal:3d} realizations done of {nreal})')
#                     progress_old = progress
            
#             # if u in data_node_labels and err_var_curr_dict[u] == 0.0:
#             if not np.isnan(sgs_dict[u]) and err_var_curr_dict[u] == 0.0:
#                 continue

#             length_dict = nx.single_source_dijkstra_path_length(G, u, cutoff=dmax, weight=edge_length_attr)
#             length_keys_list = list(length_dict.keys())
#             ind = np.argsort(np.asarray(list(length_dict.values())))
#             neighbor_labels = []
#             nn = 0
#             for i in ind:
#                 ui = length_keys_list[i]
#                 if not np.isnan(sgs_dict[ui]):
#                     prop_val[nn] = sgs_dict[ui]
#                     b[nn] = cov_model(length_dict[ui])[0]
#                     neighbor_labels.append(ui)
#                     nn = nn+1
#                     if nn >= nneighborMax:
#                         break

#             if nn == 0:
#                 # Mean and std (by kriging)
#                 if mean is not None:
#                     mu = mean_dict[u]
#                 else:
#                     mu = mean_default

#                 std = float(np.sqrt(cov0))
#                 if var is not None:
#                     std = var_update_dict[u] * std
#             else:                
#                 for i in range(nn-1):
#                     ui = neighbor_labels[i]
#                     for j in range(i+1, nn):
#                         uj = neighbor_labels[j]
#                         h = nx.dijkstra_path_length(G, source=ui, target=uj, weight=edge_length_attr)
#                         # h = nx.shortest_path_length(G, source=ui, target=uj, weight=edge_length_attr, method='dijkstra')
#                         cov_h = cov_func(h)[0]
#                         mat[i, j] = cov_h
#                         mat[j, i] = cov_h
                    
#                     mat[i, i] = cov0 + err_var_dict[ui]
                
#                 mat[nn-1, nn-1] = cov0 + err_var_dict[neighbor_labels[nn-1]]

#                 if ordinary_kriging:
#                     nmat = nn+1
#                     mat[nn, :] = 1.0
#                     mat[:, nn] = 1.0
#                     mat[nn, nn] = 0.0
#                     b[nn] = 1.0
                
#                 else:
#                     nmat = nn

#                 # Solve kriging system
#                 w = np.linalg.solve(mat[:nmat, :nmat], b[:nmat]) # w: vector of dimension nmat

#                 # Mean and std (by kriging)
#                 if mean is not None:
#                     # simple kriging
#                     std = float(np.sqrt(np.maximum(0.0, cov0 - b[:nmat].dot(w))))
#                     if var is not None:
#                         mu = float(mean_dict[u] + var_update_dict[u] * (np.asarray([1.0/var_update_dict[ui] for ui in neighbor_labels])*(prop_val[:nn] - np.asarray([mean_dict[ui] for ui in neighbor_labels]))).dot(w))
#                         std = var_update_dict[u] * std
#                     else:
#                         mu = float(mean_dict[u] + (prop_val[:nn] - np.asarray([mean_dict[ui] for ui in neighbor_labels])).dot(w))
#                 else:
#                     # ordinary kriging
#                     mu = float(prop_val[:nn].dot(w[:nn]))
#                     std = float(np.sqrt(np.maximum(0.0, cov0 - b[:nmat].dot(w))))
        
#             # Draw value in N(mu, std^2)
#             sgs_dict[u] = np.random.normal(loc=mu, scale=std)
#             err_var_curr_dict[u] = 0.0 # now the location is simulated, no more error taken into account
            
#         # Store k-th realization
#         for u in G.nodes():
#             sgs_all_dict[u].append(sgs_dict[u])

#     if verbose > 0:
#         progress = 100
#         if progress > progress_old:
#             print(f'{pid_str}SGS: {progress:3d}% ({nreal:3d} realizations done of {nreal})')
#             progress_old = progress

#     if update_graph:
#         # Set node attributes (output)
#         nx.set_node_attributes(G, sgs_all_dict, sgs_node_attr)
    
#     return sgs_all_dict
# # ----------------------------------------------------------------------------

# # ----------------------------------------------------------------------------
# def sgs_node_attribute_mp(
#         G,
#         cov_model,
#         node_attr,
#         err_std=0.0,
#         edge_length_attr=None,
#         method='ordinary_kriging',
#         mean=None,
#         searchRadius=None,
#         searchRadiusRelative=1.2,
#         nneighborMax=12,
#         nreal=1,
#         seed=None,
#         update_graph=True,
#         sgs_node_attr=None,
#         nproc=-1,
#         verbose=1):
#     """
#     Computes the same as the function :func:`sgs_node_attribute`, using multiprocessing.

#     All the parameters except `nproc` are the same as those of the function
#     :func:`sgs_node_attribute`.

#     This function launches parallel processes [parallel calls of the
#     function :func:`sgs_node_attribute`]; the set of realizations (specified by `nreal`) is
#     distributed in a balanced way over the processes.

#     The number of processes used (in parallel) is determined by the parameter `nproc` 
#     (int, default: -1); a negative number (or zero), -n <= 0, can be specified 
#     to use the total number of cpu(s) of the system except n; `nproc` is finally
#     at maximum equal to `nreal` but at least 1 by applying:
        
#     - if `nproc >= 1`, then `nproc = max(min(nproc, nreal), 1)` is used
#     - if `nproc = -n <= 0`, then `nproc = max(min(nmax-n, nreal), 1)` is used, \
#     where nmax is the total number of cpu(s) of the system (retrieved by 
#     `multiprocessing.cpu_count()`)

#     Note: if `nproc=None`, `nproc=-1` is used.

#     Note: specifying a `seed` guarantees reproducible results whatever the number
#     of processes used.

#     See function :func:`sgs_node_attribute` for details.
#     """
#     if update_graph:
#         # Set node attribute name for output
#         if sgs_node_attr is not None:
#             if not isinstance(sgs_node_attr, str):
#                 raise ValueError(f'`sgs_node_attr` must be a string (node attribute name)')
#         else:
#             sgs_node_attr = f'{node_attr}_sgs'

#     # Set number of process(es): nproc
#     if nproc is None:
#         nproc = -1
    
#     if nproc <= 0:
#         nproc = max(min(multiprocessing.cpu_count() + nproc, nreal), 1)
#     else:
#         nproc_tmp = nproc
#         nproc = max(min(int(nproc), nreal), 1)
#         if verbose > 0 and nproc != nproc_tmp:
#             print(f'Number of processes has been changed (now: nproc={nproc})')
    
#     # Set index for distributing realizations
#     q, r = np.divmod(nreal, nproc)
#     ids_proc = [i*q + min(i, r) for i in range(nproc+1)]

#     if verbose > 0:
#         print(f'Running `sgs_node_attr` on {nproc} processes...')

#     # Set seed (base)
#     if seed is None:
#         seed = np.random.randint(1, 1000000)
#     seed = int(seed)

#     # Set pool of nproc workers
#     pool = multiprocessing.Pool(nproc)
#     out_pool = []
#     for i in range(nproc):
#         # Set i-th process
#         args = (G, cov_model, node_attr)
#         kwargs = dict(
#                     err_std=err_std,
#                     edge_length_attr=edge_length_attr,
#                     method=method,
#                     mean=mean,
#                     searchRadius=searchRadius,
#                     searchRadiusRelative=searchRadiusRelative,
#                     nneighborMax=nneighborMax,
#                     nreal=ids_proc[i+1]-ids_proc[i], 
#                     seed=seed+ids_proc[i],
#                     update_graph=False, # do not update graph if multiprocessing is enabled
#                     sgs_node_attr=sgs_node_attr,
#                     pid=i,
#                     verbose=verbose*(i==0))
#         out_pool.append(pool.apply_async(sgs_node_attribute, args=args, kwds=kwargs))

#     # Properly end working process
#     pool.close() # Prevents any more tasks from being submitted to the pool,
#     pool.join()  # then, wait for the worker processes to exit.

#     # Get result from each process
#     out = [w.get() for w in out_pool]
#     if np.any([x is None for x in out]):
#         err_msg = f'`sgs_node_attribute_mp`: an error occured on a process (worker)'
#         raise ValueError(err_msg)

#     # Gather dictionaries of all processes
#     sgs_all_dict = {u:[] for u in G.nodes()}
#     for sgs_all_dict_pid_i in out:
#         for u in G.nodes():
#             sgs_all_dict[u].extend(sgs_all_dict_pid_i[u])

#     if update_graph:
#         # Set node attributes (output)
#         nx.set_node_attributes(G, sgs_all_dict, sgs_node_attr)
    
#     return sgs_all_dict
# # ----------------------------------------------------------------------------

# # ----------------------------------------------------------------------------
# def krige_node_attribute(
#         G,
#         cov_model,
#         node_attr,
#         err_std=0.0,
#         edge_length_attr=None,
#         method='ordinary_kriging',
#         mean=None,
#         use_unique_neighborhood=False,
#         searchRadius=None,
#         searchRadiusRelative=1.2,
#         nneighborMax=12,
#         update_graph=True,
#         kriging_est_node_attr=None,
#         kriging_std_node_attr=None,
#         verbose=1):
#     """
#     Interpolates a node attribute on a graph by kriging.

#     Parameters
#     ----------
#     G : networkx.Graph
#         input graph

#     cov_model : :class:`geone.covModel.CovModel1D`
#         covariance model in 1D

#     node_attr : str
#         name of the node attribute to be kriged

#     err_std : float, or str, default : 0.0
#         standard deviation of error (zero-mean Gaussian of given std):
#         - if float: the same error std is used for data at any graph node; 
#         - if str: name of the node attribute of the error std

#     edge_length_attr : str, optional
#         name of the edge attribute for length (used to compute the distance between
#         nodes);
#         by default (`None`): the edges have a length of one

#     method : str {'simple_kriging', 'ordinary_kriging'}, default: 'ordinary_kriging'
#         type of kriging;
#         note: if `method='ordinary_kriging'`, the parameter `mean` is not used

#     mean : float, or str, default : 0.0
#         kriging mean value :
#         - if float: the value is used at all graph nodes; 
#         - if str: name of the node attribute for the mean 

#         note: if `method=ordinary_kriging`, parameter `mean` is ignored

#     use_unique_neighborhood : bool, default: False
#         indicates if a unique neighborhood is used:

#         - if True: data at any graph node is taken into account, and the kriging matrix \
#         is computed once; the parameters `searchRadius`, `searchRadiusRelative`, \
#         `nneighborMax` are not used
#         - if False: only data at graph node within a search neighborhood are \
#         taken into account according to `searchRadius`, `searchRadiusRelative`, `nneighborMax`

#     searchRadius : float, optional
#         if specified, i.e. not `None`: search radius, i.e. 
#         the data at graph node at distance to the estimated graph node greater 
#         than `searchRadius` are not taken into account in the kriging system; 
#         if `searchRadius` is not `None`, then `searchRadiusRelative` is not used;
#         by default (`searchRadius=None`): `searchRadiusRelative` is used;

#     searchRadiusRelative : float, default: 1.2
#         used only if `searchRadius` is `None`;
#         the search radius is set to `searchRadiusRelative` times the range of the 
#         covariance model `cov_model`

#     nneighborMax : int, default: 12
#         maximal number of neighbors (data at graph nodes) taken into account in the
#         kriging system; the data at graph nodes the closest to the estimated graph node are
#         taken into account;
#         note: if `nneighborMax=None` or `nneighborMax<0`, then `nneighborMax` is
#         set to the number of graph nodes with informed data

#     update_graph : bool, default True
#         - if `True`: the kriging estimates is set as a node attribute (named \
#         `kriging_est_node_attr` (see below)) and the kriging standard deviation \
#         is set as a node attribute (named `kriging_std_node_attr` (see below))
#         - if `False`: the graph is not updated (`kriging_est_node_attr` and 
#         `kriging_std_node_attr` are not used)

#     kriging_est_node_attr : str, optional
#         name of the node attribute in output for kriging estimates;
#         by default (`None`),  the name `node_attr` + '_krig_est' is used

#     kriging_std_node_attr : str, optional
#         name of the node attribute in output for kriging standard deviation;
#         by default (`None`),  the name `node_attr` + '_krig_std' is used

#     verbose : int, default: 0
#         verbose mode, higher implies more printing (info)

#     Returns
#     -------
#     kriging_est : dict
#         dictionary with node labels as keys and kriging estimates as values

#     kriging_std : dict
#         dictionary with node labels as keys and kriging standard deviation as values
#     """
#     # Check cov_model
#     if not isinstance(cov_model, geone.covModel.CovModel1D):
#         raise ValueError('`cov_model` must be an instance of `geone.covModel.CovModel1D`')

#     # Prevent calculation if covariance model is not stationary
#     if not cov_model.is_stationary():
#         raise ValueError('`cov_model` is not stationary')

#     # Check node attribute
#     if node_attr not in kn.utils.get_node_attribute_names(G):
#         raise ValueError(f'{node_attr} is not a node attribute')

#     # Check edge length attribute
#     if edge_length_attr is not None and edge_length_attr not in kn.utils.get_edge_attribute_names(G):
#         raise ValueError(f'{edge_length_attr} is not an edge attribute')

#     # Set dictionary to convert node label (id) to node index, and vice versa
#     node_label2index = {u:i for i, u in enumerate(G.nodes())}
#     #node_index2label = {i:u for i, u in enumerate(G.nodes())}

#     # Set node attribute names for output
#     if update_graph:
#         if kriging_est_node_attr is not None:
#             if not isinstance(kriging_est_node_attr, str):
#                 raise ValueError('`kriging_est_node_attr` must be a string (node attribute name)')
#         else:
#             kriging_est_node_attr = f'{node_attr}_krig_est'

#         if kriging_std_node_attr is not None:
#             if not isinstance(kriging_std_node_attr, str):
#                 raise ValueError('`kriging_std_node_attr` must be a string (node attribute name)')
#         else:
#             kriging_std_node_attr = f'{node_attr}_krig_std'

#     # Get the dictionary of node attribute (property) to be kriged
#     prop_dict = nx.get_node_attributes(G, node_attr)
#     if len(prop_dict) == 0:
#         raise ValueError(f'No value for specified node attribute ({node_attr})')

#     v_first = np.asarray(list(prop_dict.values())[0]) # first value entry as array
#     if v_first.ndim > 0:
#         raise ValueError('Value of the specified attribute should by a scalar')

#     # - set nan at uninformed nodes
#     prop_dict = nx.get_node_attributes(G, node_attr, default=np.nan)

#     # Default mean value    
#     tmp = np.asarray(list(prop_dict.values()))
#     if np.isnan(tmp).all():
#         mean_default = 0.0
#     else:
#         # Set mean of data values
#         mean_default = np.nanmean(tmp)
    
#     # Method and mean
#     if method == 'simple_kriging':
#         ordinary_kriging = False
#         # Get the dictionary of mean kriging values
#         if mean is not None:
#             if isinstance(mean, float) or isinstance(mean, int):
#                 mean_dict = {u:mean for u in G.nodes()}

#             elif isinstance(mean, str):
#                 mean_dict = nx.get_node_attributes(G, mean, default=np.nan)
#                 if np.isnan(np.asarray(list(mean_dict.values()))).any():
#                     raise ValueError('Specified mean is not defined at all graph nodes')

#         else:
#             mean_dict = {u:mean_default for u in G.nodes()}
    
#     elif method == 'ordinary_kriging':
#         ordinary_kriging = True
#         if verbose > 0 and mean is not None:
#             print("WARNING: `mean` is ignored with `method='ordinary_kriging'`")

#         mean = None

#     else:
#         raise ValueError(f'`method` ({method}) unknown')

#     # Covariance function and value at 0
#     cov_func = cov_model.func() # covariance function
#     cov0 = cov_func(0.)[0] # covariance function at origin (lag=0)

#     # List of data node labels / index
#     data_node_labels = [u for u, v in prop_dict.items() if not np.isnan(v)]

#     # Number of data nodes 
#     n = len(data_node_labels)

#     if n == 0:
#         if ordinary_kriging:
#             krig_est_dict = {u: 0.0 for u in G.nodes()}

#         else: # simple kriging
#             krig_est_dict = mean_dict.copy()
#             # krig_est_dict = {u: float(mean_dict[u]) for u in G.nodes()}

#         tmp = float(np.sqrt(cov0))
#         krig_std_dict = {u: tmp for u in G.nodes()}

#         if update_graph:
#             # Set node attributes (output)
#             nx.set_node_attributes(G, krig_est_dict, kriging_est_node_attr)
#             nx.set_node_attributes(G, krig_std_dict, kriging_std_node_attr)

#         return krig_est_dict, krig_std_dict
    
#     # Get the dictionary of error variance
#     if err_std is None:
#         err_std = 0.0

#     if isinstance(err_std, float) or isinstance(err_std, int):
#         err_var = err_std**2
#         err_var_dict = {u:err_var for u in data_node_labels}

#     elif isinstance(err_std, str):
#         err_var_dict = nx.get_node_attributes(G, err_std)
#         if np.any(np.asarray([u not in err_var_dict.keys() for u in data_node_labels])):
#             raise ValueError('Specified error std is must be defined at all data nodes')
        
#         err_var_dict = {u:s**2 for u, s in err_var_dict.items()}

#     # Number of nodes in the graph
#     n_nodes = G.number_of_nodes()

#     # Do kriging
#     if use_unique_neighborhood:
#         # Initialize 
#         # - kriging matrix (mat) of order nmat
#         # - right hand side of all kriging systems (b), matrix of dimension nmat x n_nodes
#         if ordinary_kriging:
#             nmat = n+1
#             mat = np.ones((nmat, nmat))
#             mat[-1,-1] = 0.0
#         else:
#             nmat = n
#             mat = np.ones((nmat, nmat))

#         b = np.ones((nmat, n_nodes))

#         # Set kriging matrix (mat) and right hand side of all kriging systems (b)
#         for i, ui in enumerate(data_node_labels):
#             length_dict = nx.single_source_dijkstra_path_length(G, ui, cutoff=None, weight=edge_length_attr)
#             h = np.asarray([length_dict[u] for u in data_node_labels])
#             mat[i, :n] = cov_func(h)
#             h = np.asarray([length_dict[u] for u in G.nodes()])
#             b[i, :] = cov_func(h)
        
#         # Add error variance on the diagonal of the kriging matrix
#         for i, u in enumerate(data_node_labels):
#             mat[i, i] = mat[i, i] + err_var_dict[u]

#         # Solve all kriging systems
#         w = np.linalg.solve(mat, b) # w: matrix of dimension nmat x n_nodes

#         # Kriged values
#         if mean is not None:
#             # simple kriging
#             tmp = np.asarray([prop_dict[ui] - mean_dict[ui] for ui in data_node_labels]).dot(w)
#             krig_est_dict = {u: float(mean_dict[u] + v) for u, v in zip(G.nodes(), tmp)}
#         else:
#             # ordinary kriging
#             tmp = np.asarray([prop_dict[ui] for ui in data_node_labels]).dot(w[:n, :])
#             krig_est_dict = {u: float(v) for u, v in zip(G.nodes(), tmp)}

#         # Kriged standard deviation
#         ind = [node_label2index[u] for u in G.nodes()]
#         tmp = np.sqrt(np.maximum(0.0, cov0 - np.array([np.dot(w[:,i], b[:,i]) for i in ind])))
#         krig_std_dict = {u: float(v) for u, v in zip(G.nodes(), tmp)}

#     else:
#         # Limited search neighborhood

#         # Set dmax (search radius)
#         if searchRadius is not None:
#             if searchRadius <= 0.0:
#                 raise ValueError('`searchRadius` not valid (negative)')

#             dmax = searchRadius

#         else:
#             # use searchRadiusRelative
#             if searchRadiusRelative <= 0.0:
#                 raise ValueError('`searchRadiusRelative` (factor) not valid (negative)')
            
#             dmax = searchRadiusRelative * cov_model.r()

#         if nneighborMax is None or nneighborMax > n or nneighborMax < 0:
#             nneighborMax = n

#         # Initialize kriging matrix, second member and property values at neigbhors
#         mat = np.ones((nneighborMax+1, nneighborMax+1))
#         b = np.ones(nneighborMax+1) 
#         prop_val = np.ones(nneighborMax)

#         # Initialize the dictionaries of kriging estimates and standard deviation
#         krig_est_dict = {u:np.nan for u in G.nodes()}
#         krig_std_dict = {u:np.nan for u in G.nodes()}

#         if verbose > 0:
#             progress_old = 0

#         for k, u in enumerate(G.nodes()):
#             if verbose > 0:
#                 progress = int(k/n_nodes*100.0)
#                 if progress > progress_old:
#                     print(f'Kriging: {progress:3d}%')
#                     progress_old = progress
            
#             if u in data_node_labels and err_var_dict[u] == 0.0:
#                 krig_est_dict[u] = prop_dict[u]
#                 krig_std_dict[u] = 0.0
#                 continue

#             length_dict = nx.single_source_dijkstra_path_length(G, u, cutoff=dmax, weight=edge_length_attr)
#             length_keys_list = list(length_dict.keys())
#             ind = np.argsort(np.asarray(list(length_dict.values())))
#             neighbor_labels = []
#             nn = 0
#             for i in ind:
#                 ui = length_keys_list[i]
#                 if not np.isnan(prop_dict[ui]):
#                     prop_val[nn] = prop_dict[ui]
#                     b[nn] = cov_model(length_dict[ui])[0]
#                     neighbor_labels.append(ui)
#                     nn = nn+1
#                     if nn >= nneighborMax:
#                         break

#             if nn == 0:
#                 if mean is not None:
#                     # simple kriging
#                     krig_est_dict[u] = float(mean_dict[u])
#                 else:
#                     # ordinary kriging
#                     krig_est_dict[u] = mean_default

#                 krig_std_dict[u] = float(np.sqrt(cov0))
#                 continue
                
#             for i in range(nn-1):
#                 ui = neighbor_labels[i]
#                 for j in range(i+1, nn):
#                     uj = neighbor_labels[j]
#                     h = nx.dijkstra_path_length(G, source=ui, target=uj, weight=edge_length_attr)
#                     # h = nx.shortest_path_length(G, source=ui, target=uj, weight=edge_length_attr, method='dijkstra')
#                     cov_h = cov_func(h)[0]
#                     mat[i, j] = cov_h
#                     mat[j, i] = cov_h
                
#                 mat[i, i] = cov0 + err_var_dict[ui]
            
#             mat[nn-1, nn-1] = cov0 + err_var_dict[neighbor_labels[nn-1]]

#             if ordinary_kriging:
#                 nmat = nn+1
#                 mat[nn, :] = 1.0
#                 mat[:, nn] = 1.0
#                 mat[nn, nn] = 0.0
#                 b[nn] = 1.0
            
#             else:
#                 nmat = nn

#             # Solve kriging system
#             w = np.linalg.solve(mat[:nmat, :nmat], b[:nmat]) # w: vector of dimension nmat

#             # Kriged values
#             if mean is not None:
#                 # simple kriging
#                 krig_est_dict[u] = float(mean_dict[u] + (prop_val[:nn] - np.asarray([mean_dict[ui] for ui in neighbor_labels])).dot(w))
#             else:
#                 # ordinary kriging
#                 krig_est_dict[u] = float(prop_val[:nn].dot(w[:nn]))

#             # Kriged standard deviation
#             krig_std_dict[u] = float(np.sqrt(np.maximum(0.0, cov0 - b[:nmat].dot(w))))
        
#         if verbose > 0:
#             progress = 100
#             if progress > progress_old:
#                 print(f'Kriging: {progress:3d}%')
#                 progress_old = progress

#     if update_graph:
#         # Set node attributes (output)
#         nx.set_node_attributes(G, krig_est_dict, kriging_est_node_attr)
#         nx.set_node_attributes(G, krig_std_dict, kriging_std_node_attr)
    
#     return krig_est_dict, krig_std_dict
# # ----------------------------------------------------------------------------

# # ----------------------------------------------------------------------------
# def krige_node_attribute_by_branches(
#         kg,
#         cov_model,
#         node_attr,
#         err_std=0.0,
#         edge_length_attr=None,
#         method='simple_kriging',
#         mean_on_branches=None,
#         var_on_branches=None,
#         use_unique_neighborhood=False,
#         searchRadius=None,
#         searchRadiusRelative=1.2,
#         nneighborMax=12,
#         update_graph=True,
#         kriging_est_node_attr=None,
#         kriging_std_node_attr=None,
#         i0=None, 
#         i1=None,
#         pid=None,
#         verbose=1):
#     """
#     Interpolates a node attribute on a graph by kriging.

#     Parameters
#     ----------
#     G : networkx.Graph
#         input graph

#     cov_model : :class:`geone.covModel.CovModel1D`
#         covariance model in 1D

#     node_attr : str
#         name of the node attribute to be kriged

#     err_std : float, or str, default : 0.0
#         standard deviation of error (zero-mean Gaussian of given std):
#         - if float: the same error std is used for data at any graph node; 
#         - if str: name of the node attribute of the error std

#     edge_length_attr : str, optional
#         name of the edge attribute for length (used to compute the distance between
#         nodes);
#         by default (`None`): the edges have a length of one

#     method : str {'simple_kriging', 'ordinary_kriging'}, default: 'ordinary_kriging'
#         type of kriging;
#         note: if `method='ordinary_kriging'`, the parameter `mean` is not used

#     mean : float, or str, default : 0.0
#         kriging mean value :
#         - if float: the value is used at all graph nodes; 
#         - if str: name of the node attribute for the mean 

#         note: if `method=ordinary_kriging`, parameter `mean` is ignored

#     use_unique_neighborhood : bool, default: False
#         indicates if a unique neighborhood is used:

#         - if True: data at any graph node is taken into account, and the kriging matrix \
#         is computed once; the parameters `searchRadius`, `searchRadiusRelative`, \
#         `nneighborMax` are not used
#         - if False: only data at graph node within a search neighborhood are \
#         taken into account according to `searchRadius`, `searchRadiusRelative`, `nneighborMax`

#     searchRadius : float, optional
#         if specified, i.e. not `None`: search radius, i.e. 
#         the data at graph node at distance to the estimated graph node greater 
#         than `searchRadius` are not taken into account in the kriging system; 
#         if `searchRadius` is not `None`, then `searchRadiusRelative` is not used;
#         by default (`searchRadius=None`): `searchRadiusRelative` is used;

#     searchRadiusRelative : float, default: 1.2
#         used only if `searchRadius` is `None`;
#         the search radius is set to `searchRadiusRelative` times the range of the 
#         covariance model `cov_model`

#     nneighborMax : int, default: 12
#         maximal number of neighbors (data at graph nodes) taken into account in the
#         kriging system; the data at graph nodes the closest to the estimated graph node are
#         taken into account;
#         note: if `nneighborMax=None` or `nneighborMax<0`, then `nneighborMax` is
#         set to the number of graph nodes with informed data

#     update_graph : bool, default True
#         - if `True`: the kriging estimates is set as a node attribute (named \
#         `kriging_est_node_attr` (see below)) and the kriging standard deviation \
#         is set as a node attribute (named `kriging_std_node_attr` (see below))
#         - if `False`: the graph is not updated (`kriging_est_node_attr` and 
#         `kriging_std_node_attr` are not used)

#     kriging_est_node_attr : str, optional
#         name of the node attribute in output for kriging estimates;
#         by default (`None`),  the name `node_attr` + '_krig_est' is used

#     kriging_std_node_attr : str, optional
#         name of the node attribute in output for kriging standard deviation;
#         by default (`None`),  the name `node_attr` + '_krig_std' is used

#     i0 : int, optional
#         starting index in `list(G.nodes())` of nodes to be kriged (used with multiprocessing)

#     i1 : int, optional
#         ending index (excluded) in `list(G.nodes())` of nodes to be kriged (used with multiprocessing)

#     pid : int, optional
#         process id of the caller (used with multiprocessing)

#     verbose : int, default: 0
#         verbose mode, higher implies more printing (info)

#     Returns
#     -------
#     kriging_est : dict
#         dictionary with node labels as keys and kriging estimates as values

#     kriging_std : dict
#         dictionary with node labels as keys and kriging standard deviation as values
#     """
#     if verbose > 0:
#         if pid is not None:
#             pid_str = f'[pid={pid}] '
#         else:
#             pid_str = ''

#     G = kg.graph
#     Gb = kg.graph_of_branches

#     # Check cov_model
#     if not isinstance(cov_model, geone.covModel.CovModel1D):
#         raise ValueError(f'{pid_str}`cov_model_base` must be an instance of `geone.covModel.CovModel1D`')

#     # Prevent calculation if covariance model is not stationary
#     if not cov_model.is_stationary():
#         raise ValueError(f'{pid_str}`cov_model` is not stationary')

#     # Check node attribute
#     if node_attr not in kn.utils.get_node_attribute_names(G):
#         raise ValueError(f'{pid_str}{node_attr} is not a node attribute')

#     # Check edge length attribute
#     if edge_length_attr is not None and edge_length_attr not in kn.utils.get_edge_attribute_names(G):
#         raise ValueError(f'{pid_str}{edge_length_attr} is not an edge attribute')

#     # # Set dictionary to convert node label (id) to node index, and vice versa
#     # node_label2index = {u:i for i, u in enumerate(G.nodes())}
#     # #node_index2label = {i:u for i, u in enumerate(G.nodes())}

#     if update_graph:
#         # Set node attribute names for output
#         if kriging_est_node_attr is not None:
#             if not isinstance(kriging_est_node_attr, str):
#                 raise ValueError(f'{pid_str}`kriging_est_node_attr` must be a string (node attribute name)')
#         else:
#             kriging_est_node_attr = f'{node_attr}_krig_est'

#         if kriging_std_node_attr is not None:
#             if not isinstance(kriging_std_node_attr, str):
#                 raise ValueError(f'{pid_str}`kriging_std_node_attr` must be a string (node attribute name)')
#         else:
#             kriging_std_node_attr = f'{node_attr}_krig_std'

#     # Get the dictionary of node attribute (property) to be kriged
#     prop_dict = nx.get_node_attributes(G, node_attr)
#     if len(prop_dict) == 0:
#         raise ValueError(f'{pid_str}No value for specified node attribute ({node_attr})')

#     v_first = np.asarray(list(prop_dict.values())[0]) # first value entry as array
#     if v_first.ndim > 0:
#         raise ValueError(f'{pid_str}Value of the specified attribute should by a scalar')

#     # - set nan at uninformed nodes
#     prop_dict = nx.get_node_attributes(G, node_attr, default=np.nan)

#     # Covariance function and value at 0
#     cov_func = cov_model.func() # covariance function
#     cov0 = cov_func(0.)[0] # covariance function at origin (lag=0)

#     # Default mean value on branches
#     mean_default_on_branches_dict = {}
#     for i, br in enumerate(kg.branches):
#         tmp = np.asarray([prop_dict[u] for u in br])
#         if np.isnan(tmp).all():
#             mean_default_on_branches_dict[i] = 0.0
#         else:
#             mean_default_on_branches_dict[i] = np.nanmean(tmp)

#     # Method and mean, var
#     if method == 'simple_kriging':
#         ordinary_kriging = False
#         # Get the dictionary of mean kriging values on branches
#         if mean_on_branches is not None:
#             if isinstance(mean_on_branches, float) or isinstance(mean_on_branches, int):
#                 mean_on_branches_dict = {u:mean_on_branches for u in Gb.nodes()}

#             elif isinstance(mean_on_branches, str):
#                 mean_on_branches_dict = nx.get_node_attributes(Gb, mean_on_branches, default=np.nan)
#                 if np.isnan(np.asarray(list(mean_on_branches_dict.values()))).any():
#                     raise ValueError(f'{pid_str}Specified mean on branches is not defined on all graph branches')
                
#             else:
#                 raise ValueError(f'{pid_str}Specified mean on branches is not valid')

#         else:
#             mean_on_branches_dict = mean_default_on_branches_dict
    
#         if var_on_branches is not None:
#             if isinstance(var_on_branches, float) or isinstance(var_on_branches, int):
#                 tmp = np.sqrt(var_on_branches/cov0)
#                 var_update_on_branches_dict = {u:tmp for u in Gb.nodes()}

#             elif isinstance(var_on_branches, str):
#                 var_update_on_branches_dict = nx.get_node_attributes(Gb, var_on_branches, default=np.nan)
#                 if np.isnan(np.asarray(list(var_update_on_branches_dict.values()))).any():
#                     raise ValueError(f'{pid_str}Specified variance on branches is not defined on all graph branches')
                
#                 var_update_on_branches_dict = {u: np.sqrt(v/cov0) for u, v in var_update_on_branches_dict.items()}

#             else:
#                 raise ValueError(f'{pid_str}Specified variance on branches is not valid')

#         # else:
#         #     var_update_on_branches_dict = None

#     elif method == 'ordinary_kriging':
#         ordinary_kriging = True
#         if verbose > 0 and mean_on_branches is not None:
#             print(f"{pid_str}WARNING: `mean_on_branches` is ignored with `method='ordinary_kriging'`")

#         mean_on_branches = None

#         if verbose > 0 and var_on_branches is not None:
#             print(f"{pid_str}WARNING: `var_on_branches` is ignored with `method='ordinary_kriging'`")

#         var_on_branches = None

#     else:
#         raise ValueError(f'{pid_str}`method` ({method}) unknown')

#     # List of data node labels / index
#     data_node_labels = [u for u, v in prop_dict.items() if not np.isnan(v)]

#     # Number of data nodes 
#     n = len(data_node_labels)

#     # Number of nodes in the graph
#     n_nodes = G.number_of_nodes()
#     if i0 is None:
#         i0 = 0
#     if i1 is None:
#         i1 = n_nodes

#     node_list = list(G.nodes())[i0:i1]
#     node_list_len = i1 - i0

#     if n == 0:
#         if ordinary_kriging:
#             krig_est_dict = {u: 0.0 for u in node_list}

#             tmp = float(np.sqrt(cov0))
#             krig_std_dict = {u: tmp for u in node_list}

#         else: # simple kriging
#             krig_est_dict = {u: float(mean_on_branches_dict[G.nodes[u]['branch_ids_list'][0]]) for u in node_list}

#             tmp = float(np.sqrt(cov0))
#             if var_on_branches is not None:
#                 krig_std_dict = {u: tmp*var_update_on_branches_dict[G.nodes[u]['branch_ids_list'][0]] for u in node_list}
#             else:
#                 krig_std_dict = {u: tmp for u in node_list}

#         if update_graph:
#             # Set node attributes (output)
#             nx.set_node_attributes(G, krig_est_dict, kriging_est_node_attr)
#             nx.set_node_attributes(G, krig_std_dict, kriging_std_node_attr)

#         return krig_est_dict, krig_std_dict
    
#     # Get the dictionary of error variance
#     if err_std is None:
#         err_std = 0.0

#     if isinstance(err_std, float) or isinstance(err_std, int):
#         err_var = err_std**2
#         err_var_dict = {u:err_var for u in data_node_labels}

#     elif isinstance(err_std, str):
#         err_var_dict = nx.get_node_attributes(G, err_std)
#         if np.any(np.asarray([u not in err_var_dict.keys() for u in data_node_labels])):
#             raise ValueError(f'{pid_str}Specified error std is must be defined at all data nodes')
        
#         err_var_dict = {u:s**2 for u, s in err_var_dict.items()}

#     # Do kriging of nodes in node_list
#     if not use_unique_neighborhood:
#         # Limited search neighborhood

#         # Set dmax (search radius)
#         if searchRadius is not None:
#             if searchRadius <= 0.0:
#                 raise ValueError(f'{pid_str}`searchRadius` not valid (negative)')

#             dmax = searchRadius

#         else:
#             # use searchRadiusRelative
#             if searchRadiusRelative <= 0.0:
#                 raise ValueError(f'{pid_str}`searchRadiusRelative` (factor) not valid (negative)')
            
#             dmax = searchRadiusRelative * cov_model.r()

#         if nneighborMax is None or nneighborMax > n or nneighborMax < 0:
#             nneighborMax = n

#     else:
#         nneighborMax = n

#     # Initialize the dictionaries of kriging estimates and standard deviation
#     krig_est_dict = {u:np.nan for u in node_list}
#     krig_std_dict = {u:np.nan for u in node_list}

#     if verbose > 0:
#         progress_old = 0

#     for k, u in enumerate(node_list):
#         if verbose > 0:
#             progress = int(k/node_list_len*100.0)
#             if progress > progress_old:
#                 print(f'{pid_str}Kriging: {progress:3d}%')
#                 progress_old = progress
        
#         if u in data_node_labels and err_var_dict[u] == 0.0:
#             krig_est_dict[u] = prop_dict[u]
#             krig_std_dict[u] = 0.0
#             continue
        
#         # Get branch id
#         br_id = G.nodes[u]['branch_ids_list'][0] # first one if in more than one branch
#         # Get branch
#         br = kg.branches[br_id]
#         # Get index in the branch
#         ind_in_br = int(np.where(np.asarray(br)==u)[0][0])
        
#         # Compute length from u to informed nodes along the branch
#         br_len_dict = {}
#         if not np.isnan(prop_dict[u]):
#             br_len_dict[u] = 0.0

#         if edge_length_attr is None:
#             br_len = 0.0
#             for j in range(ind_in_br-1, -1, -1):
#                 br_len = br_len+1
#                 if not np.isnan(prop_dict[br[j]]):
#                     br_len_dict[br[j]] = br_len

#             br_len = 0.0
#             for j in range(ind_in_br+1, len(br)):
#                 br_len = br_len+1
#                 if not np.isnan(prop_dict[br[j]]):
#                     br_len_dict[br[j]] = br_len

#         else:
#             br_len = 0.0
#             for j in range(ind_in_br-1, -1, -1):
#                 br_len = br_len + G.edges[(br[j], br[j+1])][edge_length_attr]
#                 if not np.isnan(prop_dict[br[j]]):
#                     br_len_dict[br[j]] = br_len

#             br_len = 0.0
#             for j in range(ind_in_br+1, len(br)):
#                 br_len = br_len + G.edges[(br[j-1], br[j])][edge_length_attr]
#                 if not np.isnan(prop_dict[br[j]]):
#                     br_len_dict[br[j]] = br_len

#         br_len_keys_list = list(br_len_dict.keys())
#         ind = np.argsort(np.asarray(list(br_len_dict.values())))
#         nn = min(nneighborMax, len(ind))

#         if nn == 0:
#             if mean_on_branches is not None:
#                 # simple kriging
#                 krig_est_dict[u] = float(mean_on_branches_dict[br_id])
#             else:
#                 # ordinary kriging
#                 krig_est_dict[u] = float(mean_default_on_branches_dict[br_id])

#             krig_std_dict[u] = float(np.sqrt(cov0))
#             if var_on_branches is not None:
#                 krig_std_dict[u] = var_update_on_branches_dict[br_id] * krig_std_dict[u]
#             continue
            
#         neighbor_labels = br_len_keys_list[:nn]
#         prop_val = np.asarray([prop_dict[ui] for ui in br_len_keys_list[:nn]])
#         b = cov_model(np.asarray([br_len_dict[ui] for ui in br_len_keys_list[:nn]]))

#         if ordinary_kriging:
#             nmat = nn+1
#             mat = np.ones((nmat, nmat))
#             mat[nn, nn] = 0.0
#             b = np.hstack((b, [1.0]))
        
#         else:
#             nmat = nn
#             mat = np.zeros((nmat, nmat))

#         # could be optimized without calling nk.dijkstra_path_length below ...
#         for i in range(nn-1):
#             ui = neighbor_labels[i]
#             for j in range(i+1, nn):
#                 uj = neighbor_labels[j]
#                 h = nx.dijkstra_path_length(G, source=ui, target=uj, weight=edge_length_attr)
#                 # h = nx.shortest_path_length(G, source=ui, target=uj, weight=edge_length_attr, method='dijkstra')
#                 cov_h = cov_func(h)[0]
#                 mat[i, j] = cov_h
#                 mat[j, i] = cov_h
            
#             mat[i, i] = cov0 + err_var_dict[ui]
        
#         mat[nn-1, nn-1] = cov0 + err_var_dict[neighbor_labels[nn-1]]

#         # Solve kriging system
#         w = np.linalg.solve(mat, b) # w: vector of dimension nmat

#         if mean_on_branches is not None:
#             # simple kriging
#             # (no difference whatever is var_on_branches)
#             krig_est_dict[u] = float(mean_on_branches_dict[br_id] + (prop_val - mean_on_branches_dict[br_id]*np.ones(nn)).dot(w))
#         else:
#             # ordinary kriging
#             krig_est_dict[u] = float(prop_val.dot(w[:nn]))

#         # Kriged standard deviation
#         krig_std_dict[u] = float(np.sqrt(np.maximum(0.0, cov0 - b.dot(w))))
#         if var_on_branches is not None:
#             krig_std_dict[u] = var_update_on_branches_dict[br_id] * krig_std_dict[u]

#     if verbose > 0:
#         progress = 100
#         if progress > progress_old:
#             print(f'{pid_str}Kriging: {progress:3d}%')
#             progress_old = progress

#     if update_graph:
#         # Set node attributes (output)
#         nx.set_node_attributes(G, krig_est_dict, kriging_est_node_attr)
#         nx.set_node_attributes(G, krig_std_dict, kriging_std_node_attr)
    
#     return krig_est_dict, krig_std_dict
# # ----------------------------------------------------------------------------

# # ----------------------------------------------------------------------------
# def krige_node_attribute_by_branches_mp(
#         kg,
#         cov_model,
#         node_attr,
#         err_std=0.0,
#         edge_length_attr=None,
#         method='simple_kriging',
#         mean_on_branches=None,
#         var_on_branches=None,
#         use_unique_neighborhood=False,
#         searchRadius=None,
#         searchRadiusRelative=1.2,
#         nneighborMax=12,
#         update_graph=True,
#         kriging_est_node_attr=None,
#         kriging_std_node_attr=None,
#         verbose=1,
#         nproc=-1):
#     """
#     Computes the same as the function :func:`krige_node_attribute`, using multiprocessing.

#     All the parameters except `nproc` are the same as those of the function
#     :func:`krige_node_attribute`.

#     This function launches parallel processes [parallel calls of the
#     function :func:`krige_node_attribute`]; the set of nodes to be kriged is distributed in a 
#     balanced way over the processes.

#     The number of processes used (in parallel) is determined by the parameter `nproc` 
#     (int, default: -1); a negative number (or zero), -n <= 0, can be specified 
#     to use the total number of cpu(s) of the system except n; `nproc` is finally
#     at maximum equal to `G.number_of_nodes()` but at least 1 by applying:
        
#     - if `nproc >= 1`, then `nproc = max(min(nproc, G.number_of_nodes()), 1)` is used
#     - if `nproc = -n <= 0`, then `nproc = max(min(nmax-n, G.number_of_nodes()), 1)` is used, \
#     where nmax is the total number of cpu(s) of the system (retrieved by 
#     `multiprocessing.cpu_count()`)

#     Note: if `nproc=None`, `nproc=-1` is used.

#     See function :func:`krige_node_attribute` for details.
#     """
#     if update_graph:
#         # Set node attribute names for output
#         if kriging_est_node_attr is not None:
#             if not isinstance(kriging_est_node_attr, str):
#                 raise ValueError('`kriging_est_node_attr` must be a string (node attribute name)')
#         else:
#             kriging_est_node_attr = f'{node_attr}_krig_est'

#         if kriging_std_node_attr is not None:
#             if not isinstance(kriging_std_node_attr, str):
#                 raise ValueError('`kriging_std_node_attr` must be a string (node attribute name)')
#         else:
#             kriging_std_node_attr = f'{node_attr}_krig_std'

#     G = kg.graph

#     # Number of graph nodes
#     n_nodes = G.number_of_nodes()

#     # Set number of process(es): nproc
#     if nproc is None:
#         nproc = -1
    
#     if nproc <= 0:
#         nproc = max(min(multiprocessing.cpu_count() + nproc, n_nodes), 1)
#     else:
#         nproc_tmp = nproc
#         nproc = max(min(int(nproc), n_nodes), 1)
#         if verbose > 0 and nproc != nproc_tmp:
#             print(f'Number of processes has been changed (now: nproc={nproc})')

#     # Set index for distributing realizations
#     q, r = np.divmod(n_nodes, nproc)
#     ids_proc = [i*q + min(i, r) for i in range(nproc+1)]

#     if verbose > 0:
#         print(f'Running `krige_node_attribute_by_branches` on {nproc} processes...')

#     # Set pool of nproc workers
#     pool = multiprocessing.Pool(nproc)
#     out_pool = []
#     for i in range(nproc):
#         # Set i-th process
#         args = (kg, cov_model, node_attr)
#         kwargs = dict(
#                     err_std=err_std,
#                     edge_length_attr=edge_length_attr, 
#                     method=method,
#                     mean_on_branches=mean_on_branches,
#                     var_on_branches=var_on_branches,
#                     use_unique_neighborhood=use_unique_neighborhood,
#                     searchRadius=searchRadius, 
#                     searchRadiusRelative=searchRadiusRelative, 
#                     nneighborMax=nneighborMax,
#                     update_graph=False, # do not update graph if multiprocessing is enabled
#                     kriging_est_node_attr=kriging_est_node_attr,
#                     kriging_std_node_attr=kriging_std_node_attr,
#                     i0=ids_proc[i],
#                     i1=ids_proc[i+1], 
#                     pid=i,
#                     verbose=verbose*(i==0))
#         out_pool.append(pool.apply_async(krige_node_attribute_by_branches, args=args, kwds=kwargs))

#     # Properly end working process
#     pool.close() # Prevents any more tasks from being submitted to the pool,
#     pool.join()  # then, wait for the worker processes to exit.

#     # Get result from each process
#     out = [w.get() for w in out_pool]
#     if np.any([x is None for x in out]):
#         err_msg = f'`krige_node_attribute_mp`: an error occured on a process (worker)'
#         raise ValueError(err_msg)

#     # Gather dictionaries of all processes
#     krig_est_dict = {}
#     krig_std_dict = {}
#     for krig_est_dict_pid_i, krig_std_dict_pid_i in out:
#         krig_est_dict.update(krig_est_dict_pid_i)
#         krig_std_dict.update(krig_std_dict_pid_i)

#     if update_graph:
#         # Set node attributes (output)
#         nx.set_node_attributes(G, krig_est_dict, kriging_est_node_attr)
#         nx.set_node_attributes(G, krig_std_dict, kriging_std_node_attr)
    
#     return krig_est_dict, krig_std_dict
# # ----------------------------------------------------------------------------

# # ----------------------------------------------------------------------------
# def fractal_dimension_with_memb(
#         g, 
#         r_values,
#         nrepeat=1,
#         min_n_edges_factor=0.01,
#         account_for_original_graph=True,
#         confidence_level=0.95,
#         pos_attr='pos',
#         memb_edge_weight=None,
#         memb_node_attr_mode='mean',
#         memb_node_attr_list=None,
#         return_poly_fit_params=True,
#         return_poly_fit_params_cov=True,
#         return_mean_edge_length_list=True,
#         return_n_nodes_list=True,
#         return_n_edges_list=True,
#         return_g_red_r_list=False,
#         return_g_red_list=False,
#         seed=None,
#         verbose=0):
#     """
#     Computes fractal dimension using MEMB algorithm.

#     This function computes a sequence of reduced graphs by applying MEMB to 
#     the given graph (with increasing radius) (function `reduce_graph_memb`); 
#     then assuming the relationship
    
#     .. math::
#         N \propto \lambda^{-d_f}

#     where :math:`N` is number of graph nodes and :math:`\lambda` the mean edge 
#     length in reduced graph, the fractal dimension :math:`d_f` is obtained by 
#     fitting a line on the log-log plot of :math:`N` as function of :math:`\lambda`.

#     Parameters
#     ----------
#     g : networkx.Graph
#         input graph
    
#     r_values : sequence of ints or floats
#         sequence of radius values for reducing input graph, the sequence should
#         be in ascending order);
#         see function `reduce_graph_memb`
    
#     n_repeat : int, default: 1
#         number of times each radius value in `r_values` are used to compute
#         a reduced graph
    
#     min_n_edges_factor : float, default: 0.01
#         factor used to compute the minimal number of edges for a reduced graph;
#         it is defined as the greatest integer less than or equal to  
#         `min_n_edges_factor` times the number of edges in the input graph;
#         as soon as the reduced graph has a smaller number of edges, this graph is
#         not retained and the sequence of reduced graph is stopped (hence, it 
#         is important that the sequence `r_values` is in ascending order)

#     account_for_original_graph : bool, default: True
#         indicates if the input graph is taken into account for fitting the 
#         line on the log-log plot to compute the fractal dimension (see above)

#     confidence_level : float, default: 0.95
#         confidence level, float in the interval (0, 1), to compute the 
#         uncertainty for the fractal dimension

#     pos_attr : str, default: 'pos'
#         name of the node attribute for position

#     memb_edge_weight : str, optional
#         name of the edge attribute used as weight for the MEMB algorithm;
#         parameter `edge_weight` passed to the function `reduce_graph_memb`

#     memb_node_attr_mode : str {'mean', 'center'} (or list of strs), default: 'mean'
#         string or list of strings, the strings indicate how the node attributes 
#         are computed for the MEMB algorithm;
#         parameter `node_attr_mode` passed to the function `reduce_graph_memb`

#     memb_node_attr_list : list of strs, optional
#         list of node attributes of the input graph to be included in 
#         the reduced graphs for the MEMB algorithm;
#         parameter `node_attr_list` passed to the function `reduce_graph_memb`
#         node_attr_list;
#         by default (`None`): `memb_node_attr_list` is set [`pos_attr`]

#     return_poly_fit_params : bool, default: True
#         if `True`, the array (shape (2, )) of parameters of the fitted 
#         polynom (line) is returned

#     return_poly_fit_params_cov : bool, default: True
#         if `True`, the array (shape (2, 2)) of parameters covariance of the 
#         fitted polynom (line) is returned

#     return_mean_edge_length_list : bool, default: True
#         if `True`, the list of mean edge length of graphs used to compute the 
#         fractal dimension (by line fitting on the log-log plot (see above)) is 
#         returned; note that if `account_for_original_graph=True`: the entry at
#         index 0 corresponds to the input graph, and the next entries to the 
#         reduced graphs

#     return_n_nodes_list : bool, default: True
#         if `True`, the list of number of nodes of graphs used to compute the 
#         fractal dimension (by line fitting on the log-log plot (see above)) is 
#         returned; note that if `account_for_original_graph=True`: the entry at
#         index 0 corresponds to the input graph, and the next entries to the 
#         reduced graphs

#     return_n_edges_list : bool, default: True
#         if `True`, the list of number of edges of graphs used to compute the 
#         fractal dimension (by line fitting on the log-log plot (see above)) is 
#         returned; note that if `account_for_original_graph=True`: the entry at
#         index 0 corresponds to the input graph, and the next entries to the 
#         reduced graphs

#     return_g_red_r_list : bool, default: False
#         if `True`, the list of radius values for computing reduced graphs used 
#         to compute the fractal dimension (by line fitting on the log-log plot 
#         (see above)) is returned

#     return_g_red_list : bool, default: False
#         if `True`, the list of reduced graphs used to compute the fractal 
#         dimension (by line fitting on the log-log plot (see above)) is returned
    
#     seed : int, optional
#         seed number for initializing the random number generator

#     verbose : int, default: 0
#         verbose mode, larger value implies more info printed

#     Returns
#     -------
#     df : float
#         fractal dimension
    
#     df_delta : float
#         uncertainty for fractal dimension, the interval `df +/- df_delta` 
#         corresponds to the confidence interval at a confidence level 
#         `confidence_level` derived from Gaussian distribution
    
#     poly_fit_params : array of shape (2, ), optional
#         returned if `return_poly_fit_params=True`:
#         array of parameters of the fitted polynom (line); in particular:
#         `df = -poly_fit_params[0]`

#     poly_fit_params_cov : array of shape (2, 2), optional
#         returned if `return_poly_fit_params_cov=True`:
#         array (shape (2, 2)) of parameters covariance of the fitted polynom 
#         (line); in particular: 
#         `df_delta = scipy.stats.norm.ppf((1+confidence_level)/2) * np.sqrt(poly_fit_params_cov[0, 0])`
    
#     mean_edge_length_list : list of floats, optional
#         returned if `return_mean_edge_length_list=True`: 
#         the list of mean edge length of graphs used to compute the 
#         fractal dimension (by line fitting on the log-log plot (see above)) is 
#         returned; note that if `account_for_original_graph=True`: the entry at
#         index 0 corresponds to the input graph, and the next entries to the 
#         reduced graphs

#     n_nodes_list : list of ints, optional
#         returned if `return_n_nodes_list=True`:
#         if `True`, the list of number of nodes of graphs used to compute the 
#         fractal dimension (by line fitting on the log-log plot (see above)) is 
#         returned; note that if `account_for_original_graph=True`: the entry at
#         index 0 corresponds to the input graph, and the next entries to the 
#         reduced graphs

#     n_edges_list : list of ints, optional
#         returned if `return_n_edges_list=True`:
#         if `True`, the list of number of edges of graphs used to compute the 
#         fractal dimension (by line fitting on the log-log plot (see above)) is 
#         returned; note that if `account_for_original_graph=True`: the entry at
#         index 0 corresponds to the input graph, and the next entries to the 
#         reduced graphs

#     g_red_r_list : bool, list of ints or floats, optional
#         returned if `return_g_red_r_list=True`:
#         if `True`, the list of radius values for computing reduced graphs used 
#         to compute the fractal dimension (by line fitting on the log-log plot 
#         (see above)) is returned

#     g_red_list : list of networkx.Graph, optional
#         returned if `return_g_red_list=True`:
#         if `True`, the list of reduced graphs used to compute the fractal 
#         dimension (by line fitting on the log-log plot (see above)) is returned
#     """
#     if seed is not None:
#         np.random.seed(seed)

#     if memb_node_attr_list is None:
#         memb_node_attr_list = [pos_attr]

#     g_n_nodes = g.number_of_nodes()
#     g_n_edges = g.number_of_edges()
#     min_n_edges = int(min_n_edges_factor * g.number_of_edges())
    
#     if verbose > 0:
#         print(f'Original graph: #nodes = {g_n_nodes}, #edges = {g_n_edges}; min #edges for red. graph: {min_n_edges}')

#     mean_edge_length_list = []
#     n_nodes_list = []

#     if return_n_edges_list:
#         n_edges_list = []

#     if return_g_red_r_list:
#         g_red_r_list = []
    
#     if return_g_red_list:
#         g_red_list = []

#     if account_for_original_graph:
#         mean_edge_length_list.append(
#             np.mean(np.sqrt(
#                 np.sum(np.asarray(
#                     [np.asarray(g.nodes[u][pos_attr]) - np.asarray(g.nodes[v][pos_attr]) for u, v in g.edges()]
#                 )**2, axis=1)
#             ))
#         )
        
#         n_nodes_list.append(g_n_nodes)
        
#         if return_n_edges_list:
#             n_edges_list.append(g_n_edges)
        
#     ok = True
#     t1_all = time.time()
#     for r in r_values:
#         if not ok:
#             break
#         for _ in range(nrepeat):
#             t1 = time.time()
#             g_red = reduce_graph_memb(
#                 g, 
#                 r=r, 
#                 edge_weight=memb_edge_weight,
#                 node_attr_mode=memb_node_attr_mode,
#                 node_attr_list=memb_node_attr_list,
#                 return_center_dict=False, 
#                 return_new_id_dict=False, 
#                 return_orig_id_dict=False,
#                 return_box_diameter_dict=False,
#                 seed=None,
#                 verbose=verbose-1)
#             t2 = time.time()

#             g_red_n_nodes = g_red.number_of_nodes()
#             g_red_n_edges = g_red.number_of_edges()
#             if verbose > 0:
#                 print(f'MEMB, r = {r:.3g}: #nodes = {g_red_n_nodes}, #edges = {g_red_n_edges}, elapsed time : {t2-t1:.2g} sec')

#             if g_red_n_edges < min_n_edges:
#                 if verbose > 0:
#                     print(f'Reduced graph has too few edges (less than {min_n_edges}, not taken into account, stop reducing graph!')
#                 ok = False
#                 break

#             mean_edge_length_list.append(
#                 np.mean(np.sqrt(
#                     np.sum(np.asarray(
#                         [np.asarray(g_red.nodes[u][pos_attr]) - np.asarray(g_red.nodes[v][pos_attr]) for u, v in g_red.edges()]
#                     )**2, axis=1)
#                 ))
#             )
            
#             n_nodes_list.append(g_red_n_nodes)
            
#             if return_n_edges_list:
#                 n_edges_list.append(g_red_n_edges)
        
#             if return_g_red_r_list:
#                 g_red_r_list.append(r)

#             if return_g_red_list:
#                 g_red_list.append(g_red)
        
#     t2_all = time.time()
#     if verbose > 0:
#         print(f'Total elapsed time (MEMBs): {t2_all-t1_all:.2g} sec')

#     if verbose > 0:
#         print(f'Compute fractal dimension (line fitting on log-log plot of number of nodes as function of mean edge length)')

#     # Compute fractal dimension (line fitting on log-log plot of number of nodes as function of mean edge length)
#     poly_fit_params, poly_fit_params_cov = np.polyfit(np.log10(mean_edge_length_list), np.log10(n_nodes_list), 1, cov=True)

#     # Get fractal dimension df, and df_delta for given confidence level
#     t = scipy.stats.norm.ppf((1+confidence_level)/2)
#     df = -poly_fit_params[0]
#     df_delta = t * np.sqrt(poly_fit_params_cov[0, 0])

#     out = [df, df_delta]

#     if return_poly_fit_params:
#         out.append(poly_fit_params)

#     if return_poly_fit_params_cov:
#         out.append(poly_fit_params_cov)

#     if return_mean_edge_length_list:
#         out.append(mean_edge_length_list)

#     if return_n_nodes_list:
#         out.append(n_nodes_list)

#     if return_n_edges_list:
#         out.append(n_edges_list)

#     if return_g_red_r_list:
#         out.append(g_red_r_list)

#     if return_g_red_list:
#         out.append(g_red_list)

#     out = tuple(out)

#     return out
# # ----------------------------------------------------------------------------

# # ----------------------------------------------------------------------------
# # Multiprocessing : not efficient...
# # ---------------------------------------------
# # Set local function to use in pool below
# # ---------------------------------------------
# def f_loc(g, nsteps, start_node_list, seed_list, edge_weight=None, pos_attr='pos'):
#     node_pos = {k:np.asarray(v) for k, v in nx.get_node_attributes(g, pos_attr).items()}

#     nrw = len(start_node_list)
#     rw_square_dist_loc = np.zeros((nrw, nsteps))

#     for k in range(nrw):
#         u = start_node_list[k]
#         pos_start = node_pos[u]
#         seed_rw = seed_list[k]
#         rw = graph_random_walk(g, u, nsteps, edge_weight=edge_weight, seed=seed_rw)
#         rw_square_dist_loc[k] = np.sum((np.asarray([node_pos[v] for v in rw[1:]]) - pos_start)**2, axis=1)
    
#     return rw_square_dist_loc
# # ---------------------------------------------
# def graph_random_walk_square_distance(
#             g, 
#             nrw, 
#             nsteps, 
#             edge_weight=None, 
#             pos_attr='pos',
#             u0=None, 
#             seed=None,
#             return_full_array=False,
#             use_multiprocessing=True,
#             nproc=-1):
#     """
#     Computes the mean and std. dev. of square distance after each time step, for several random walks.

#     A total of `nrw` random walks of `nsteps` are performed, each one starts at node `u0` 
#     if given or at a random node otherwise and the mean and standard deviation of Euclidean 
#     square distance from the starting node, after each time step, is computed. At each step, 
#     the walker randomly chooses an edge from the currently visited node, the probability 
#     distribution is proportional to the edge weights (all the same if not given).

#     Parameters
#     ----------
#     g : networkx.Graph
#         graph
    
#     nrw : int
#         number of random walks

#     nsteps : int
#         number of time steps
    
#     edge_weight : str, optional
#         name of the edge attribute used as weight (probability to choose 
#         an edge is proportional to its weight);
#         by default (`None`) : all edges have a weight of 1

#     pos_attr : str, default: 'pos'
#         name of the attribute attached to nodes defining the position as a 
#         sequence of 2 or 3 floats (used to compute distances)
    
#     u0 : node key, optional
#         starting node in `g` for all random walks;
#         by default (`None`) a new node is randomly selected for every random walk

#     seed : int, optional
#         seed number for initializing the random number generator
    
#     return_full_array : bool, default: False
#         if `True`, the full array of shape `(nrw, nsteps)` of the 
#         square distances from the starting node, after each time step, 
#         of every random walk is returned.

#     use_multi_processing : bool, default: True
#         if `True`: multiprocessing is used (see `nproc`)

#     nproc : int, default: -1
#         if `use_multi_processing=True`: the number of processes used 
#         (in parallel) is n, and determined by the parameter `nproc` as 
#         follows:

#         - if `nproc>0`: n = `nproc`;
#         - if `nproc<= 0`: n = max(nmax+`nproc`, 1), where nmax is the total
#         number of cpu(s) of the system (retrieved by `multiprocessing.cpu_count()`), 
#         i.e. all cpus except `nproc` is used (but at least one).

#     Returns
#     -------
#     rw_square_dist_mean : 1d-array of floats of shape `(nsteps,)`
#         mean square distance from the starting node, computed over all random 
#         walks
    
#     rw_square_dist : 2d-array of floats of shape `(nrw, nsteps)`, optional
#         returned if `return_full_array=True`: square distances from the 
#         starting node, after each time step, of every random walk
#     """

#     node_pos = {k:np.asarray(v) for k, v in nx.get_node_attributes(g, pos_attr).items()}

#     if seed is not None:
#         np.random.seed(seed)

#     if u0 is None:
#         node_list = list(node_pos.keys())
#         ind = (len(node_list)*np.random.random(nrw)).astype('int')
#         start_node_list = [node_list[i] for i in ind]
#     else:
#         start_node_list = nrw*[u0]
    
#     seed_start = np.random.randint(2**32 - nrw + 1)
#     seed_list = list(range(seed_start, seed_start+nrw))

#     if use_multiprocessing:

#         import multiprocessing

#         # Set number of processes (n)
#         if nproc > 0:
#             n = nproc
#         else:
#             n = min(multiprocessing.cpu_count()+nproc, 1)

#         # Set index for sharing task
#         q, r = np.divmod(nrw, n)
#         rw_ids_proc = [i*q + min(i, r) for i in range(n+1)]

#         # Set pool of n workers
#         pool = multiprocessing.Pool(n)
#         out_pool = []
#         for i in range(n):
#             # Set i-th process
#             out_pool.append(pool.apply_async(
#                         f_loc, 
#                         args=(g, nsteps, start_node_list[rw_ids_proc[i]:rw_ids_proc[i+1]], seed_list[rw_ids_proc[i]:rw_ids_proc[i+1]]),
#                         kwds=dict(edge_weight=edge_weight, pos_attr=pos_attr)
#                     ))

#         # Properly end working process
#         pool.close() # Prevents any more tasks from being submitted to the pool,
#         pool.join()  # then, wait for the worker processes to exit.

#         # Get result from each process
#         out = [p.get() for p in out_pool]
#         rw_square_dist = np.vstack(out)

#     else:
#         rw_square_dist = np.zeros((nrw, nsteps))

#         for k in range(nrw):
#             u = start_node_list[k]
#             pos_start = node_pos[u]
#             seed_rw = seed_list[k]
#             rw = graph_random_walk(g, u, nsteps, edge_weight=edge_weight, seed=seed_rw)
#             rw_square_dist[k] = np.sum((np.asarray([node_pos[v] for v in rw[1:]]) - pos_start)**2, axis=1)
        
#     rw_square_dist_mean = rw_square_dist.mean(axis=0)

#     if return_full_array:
#         return rw_square_dist_mean, rw_square_dist

#     return rw_square_dist_mean
# # ----------------------------------------------------------------------------

# # ----------------------------------------------------------------------------
# def remove_simple_cycles(
#         g, 
#         max_cycle_len=None, 
#         pos_attr='pos'):
#     """
#     Simplifies a graph by removing simple cycles.

#     The function `networkx.simple_cyles` is used to identify the cycles.
#     The nodes of any cycle whose length (in number of nodes) not exceeding 
#     `max_cylcle_len` (if specified) are removed and replaced by one node, 
#     and the edges issued from the removed nodes are removed and replaced 
#     by new edges issued from the new node. The position of the new nodes is 
#     set to the mean of the position of the nodes that it replaces.

#     Notes: 
    
#     - edge attributes and node attributes (other than the position) \
#     are not considered for new edges and new nodes.
#     - the node labels in the output graph are converted to integers \
#     starting from 0
        
#     Parameters
#     ----------
#     g : networkx.Graph
#         input graph

#     max_cycle_len : int, optional
#         maximum cycle length (in number of nodes): cycles with more nodes
#         are kept;
#         by default (`None`): all cycles are removed    
    
#     pos_attr : str, default: 'pos'
#         name of the node attribute for position

#     Returns
#     -------
#     g_out : networkx.Graph
#         output graph (simplifed as described above)
#     """

#     # Copy graph to initialize output graph (g_out)
#     g_out = g.copy()
    
#     # Number of nodes
#     n_nodes = g_out.number_of_nodes()

#     # Convert node labels to integers from 0 to n_nodes
#     g_out = nx.convert_node_labels_to_integers(g_out)

#     # Get simple cycle to be removed
#     if max_cycle_len is not None:
#         cycle_nodes_list = [c for c in nx.simple_cycles(g_out) if len(c) <= max_cycle_len]
#     else:
#         cycle_nodes_list = [c for c in nx.simple_cycles(g_out)]

#     # Number of new nodes
#     n_new_nodes = len(cycle_nodes_list)

#     if n_new_nodes == 0:
#         return g_out

#     # Set list of node to remove
#     nodes_list_to_remove = [u for c in cycle_nodes_list  for u in c]

#     # Set new nodes (id): cycle "i" becomes new node "n_nodes+i"
#     new_node_id = range(n_nodes, n_nodes + n_new_nodes)

#     # Get node position of new nodes (mean of nodes in the cycle)
#     new_node_pos = np.asarray([np.mean(np.asarray([np.asarray(g_out.nodes[u][pos_attr]) for u in c]), axis=0) for c in cycle_nodes_list])

#     # Set dictionary with
#     # - key: current node id
#     # - value: list of the cycle ids containing the node
#     node2cycle = {k:[] for k in range(n_nodes)}
#     for i, c in enumerate(cycle_nodes_list):
#         for u in c:
#             node2cycle[u].append(i)
#     # -> node u s.t. node2cycle[u] is a non empty list will be removed

#     # Get new edges according to the edges involving nodes that will be removed
#     new_edges = np.empty((0, 2), dtype='int')
#     for i, c in enumerate(cycle_nodes_list):
#         for u in c:
#             edges = []
#             for v in g.neighbors(u):
#                 if len(node2cycle[v]) == 0:
#                     edges.append((n_nodes + i, v))
#                 else:
#                     for j in node2cycle[v]:
#                         if j != i:
#                             edges.append((n_nodes + min(i, j), n_nodes + max(i, j)))
#             if len(edges):
#                 edges = np.unique(edges, axis=0)
#                 new_edges = np.unique(np.vstack((new_edges, edges)), axis=0)

#             # edges = [(n_nodes + min(i, j), n_nodes + max(i, j)) for j in node2cycle[u] if j != i]
#             # if len(edges):
#             #     edges = np.unique(edges, axis=0)
#             #     new_edges = np.unique(np.vstack((new_edges, edges)), axis=0)
#             # edges = [(n_nodes + i, v) for v in g.neighbors(u) if len(node2cycle[v]) == 0]
#             # if len(edges):
#             #     new_edges = np.vstack((new_edges, edges))

#     # # Get new edges according to the edges involving nodes that will be removed
#     # new_edges = []
#     # for i, c in enumerate(cycle_nodes_list):
#     #     for u in c:
#     #         for v in g.neighbors(u):
#     #             if len(node2cycle[v]) == 0:
#     #                 new_edges.append((n_nodes + i, v))
#     #                 print('b', new_edges)
#     #             else:
#     #                 for j in node2cycle[v]:
#     #                     if j != i:
#     #                         new_edges.append((n_nodes + i, n_nodes + j))
#     #                         print('c', new_edges)
    
#     # Remove nodes in cycle
#     g_out.remove_nodes_from(nodes_list_to_remove)

#     # Add new nodes with position 
#     g_out.add_nodes_from(new_node_id)
#     nx.set_node_attributes(g_out, {k:v for k, v in zip(new_node_id, new_node_pos)}, pos_attr)

#     # Add new edges
#     g_out.add_edges_from(new_edges)
    
#     # Finally: convert node labels to integers
#     g_out = nx.convert_node_labels_to_integers(g_out)

#     return g_out    
# # ----------------------------------------------------------------------------


#####
# # ----------------------------------------------------------------------------
# def reduce_graph_memb(
#         g, 
#         r=1, 
#         edge_weight=None,
#         node_attr_mode='mean', 
#         node_attr_list=None,
#         return_center_dict=False,
#         return_new_id_dict=False,
#         return_orig_id_dict=False,
#         return_box_diameter_dict=False,
#         seed=None,
#         verbose=0):
#     """
#     Maximum Excluded Mass Burning (MEMB) - for weighted and unweigted graphs.

#     The Maximum-Excluded Mass Burning (MEMB) algorithm, based 
#     on a radius `r`, consists of a box covering of the input graph nodes
#     where each box (cluster) has a diameter of at most `2*r` and each box
#     is connected, i.e. for any pair of nodes in the box, there exists a path 
#     entirely in the box between the two nodes.

#     A reduced graph is provided by the box covering, where
#     - the (new) nodes are the boxes (clusters)
#     - two nodes C1, C2, are linked by an edge if there exists two nodes \
#     u1, u2 of the input graph in the clusters corresponding to C1, C2 that \
#     are linked with an edge

#     The idea of the MEMB algorithm is to used clusters centered on nodes 
#     of the input graph: for a node c in the graph, the center of a cluster, 
#     all nodes u at a distance from c less than or equal to `r` are can be 
#     considered in the cluster of center c.
    
#     The algorithm proceeds in 3 steps.

#     1. Identifying the centers in the original graph.    
    
#     (i) All nodes are marked as uncovered and non-centers.
    
#     (ii) For all non-center nodes x (including already covered nodes),
#     the excluded mass is computed, i.e. the number of uncovered nodes
#     in the cluster centered at x, and among these nodes, select a 
#     node c realizing the maximum of the excluded mass, and mark c 
#     as a center.

#     (iii) Mark all the nodes in the cluster centered at c (i.e. nodes 
#     at a distance from c less than or equal to r) as covered

#     (iv) Repeat steps (ii) and (iii) until all nodes are covered.

#     2. Creating the clusters.

#     (i) Assign an id to each center.

#     (ii) For all nodes, compute the central distance, i.e. the distance
#     to the nearest center.

#     (iii) Set a list of all non-center nodes, sorted according to 
#     increasing central distance.

#     (iv) Take the first node u of the list, select one of its 
#     neighbors, v, with a smaller central distance, and assign the id
#     of v to u. Remove the node u from the list.

#     (v) Repeat (iv) until the list is empty.

#     3. Building the reduced graph.

#     Finally, the nodes having the same id form one cluster (connected
#     in the original graph by construction). The reduced graph is 
#     then computed such that:
    
#     - one node is one cluster, with node attributes (if any) defined \
#     according to the keyword argument `node_attr_mode` from the attributes \
#     of the original graph nodes in the cluster (`node_attr_mode='center'`: \
#     atttribute from the central node is used, `node_attr_mode='mean'`: mean \
#     is used)
#     - an edge is set between two nodes if the two corresponding clusters \
#     are linked in the original graph

#     An explicit list of node attributes can be specified by the parameters
#     `node_attr_list`, and the corresponding operation for the aggregation 
#     (in cluster) can be specified in a list by the parameters `node_attr_mode`. 
#     By default, all node attributes are considered, and the same operation is 
#     applied.
    
#     Parameters
#     ----------
#     g : networkx.Graph
#         input graph
    
#     r : int or float, default: 1
#         radius of a cluster (distance to a central node; every node connected 
#         to a central node with at most `r` may belong to the same cluster)
    
#     edge_weight : str, optional
#         name of the edge attribute used as weight (for computing distance
#         between nodes);
#         by default (`None`) : all edges have a weight of 1

#     node_attr_mode : str {'mean', 'center'} (or list of strs), default: 'mean'
#         string or list of strings, the strings indicate how the node attributes 
#         are computed in the reduced graph (see above)

#         - if string: the same operation (mode) is used for all considered \
#         node attributes
#         - if list of strings: the length of the list must be equal to the \
#         the list `node_attr_list`, and `node_attr_mode[i]` indicates the \
#         operation (mode) used for the node attribute `node_attr_list[i]`

#     node_attr_list : list of strs, optional
#         list of node attributes of the input graph to be included in 
#         the reduced graph;
#         by default (`None`): the list of all nodes attributes is considered

#     return_center_dict : bool, default: False
#         if `True`, the dictionary `center_dict` is returned, where
#         the keys are the node ids in the input graph, and the values
#         are booleans: `True` if the node is a center is the output graph
#         and `False` otherwise
    
#     return_new_id_dict : bool, default: False
#         if `True`, the dictionary `new_id_dict` is returned, where
#         the keys are the node ids in the input graph, and the values
#         are the corresponding node ids in the output graph (cluster ids)
        
#     return_orig_id_dict : bool, default: False
#         if `True`, the dictionary `orig_id_dict` is returned, where
#         the keys are the node ids in the output graph (cluster), and the 
#         values are the lists of the corresponding node id in the input graph
#         (original id)
        
#     return_box_diameter_dict : bool, default: False
#         if `True`, the dictionary `box_diameter_dict` is returned, where
#         the keys are the node ids in the output graph (cluster), and the 
#         values are the diameter of the box (cluster), i.e. the maximal
#         distance (in the original graph) between two nodes in the cluster
    
#     seed : int, optional
#         seed number for initializing the random number generator

#     verbose : int, default: 0
#         verbose mode, larger value implies more info printed

#     Returns
#     -------
#     g_red : networkx.Graph
#         reduced graph (see above)
    
#     center_dict : dict, optional
#         returned if `return_center_dict=True`: dictionary indicating if 
#         a node in the input graph is a center (see `return_center_dict` above);
#         note: this property can be attached to nodes of the input graph `g` with:
#         `nx.set_node_attributes(g, center_dict, 'center')`

#     new_id_dict : dict, optional
#         returned if `return_new_id_dict=True`: dictionary giving the new id
#         (of node in the output graph) (value) for each original id (of node in 
#         the input graph) (key) (see `return_new_id_dict` above);
#         note: this property can be attached to nodes of the input graph `g` with:
#         `nx.set_node_attributes(g, new_id_dict, 'new_id')`

#     orig_id_dict : dict, optional
#         returned if `return_orig_id_dict=True`: dictionary giving the list of 
#         original ids (of nodes in the input graph) (value) for each new id (of 
#         node in the output graph) (key) (see `return_orig_id_dict` above);
#         note: this property can be attached to nodes of the output graph `g_red` with:
#         `nx.set_node_attributes(g_red, orig_id_dict, 'orig_id')`

#     box_diameter_dict : dict, optional
#         returned if `return_box_diameter_dict=True`: dictionary giving the diameter 
#         of the box (cluster) (value), for each new id (of node in the output graph) 
#         (key) (see `return_box_diameter_dict` above);
#         note: this property can be attached to nodes of the output graph `g_red` with:
#         `nx.set_node_attributes(g_red, box_diameter_dict, 'box_diameter')`

#     References
#     ----------
#     - C. Song, L. K. Gallos, S. Havlin, and H. A. Makse (2007), \
#     How to calculate the fractal dimension of a complex network: \
#     the box covering algorithm, doi: 10.1088/1742-5468/2007/03/P03006}
#     """
#     if seed is not None:
#         np.random.seed(seed)

#     # Check node attributes of the input graph to be considered in the reduced graph
#     # and corresponding operation (mode)

#     # Get keys (attribute names) in the input graph
#     # # - all keys
#     # keys = [list(g.nodes[i].keys()) for i in g.nodes()] # list of lists
#     # node_attr_all = np.unique([xij for xi in keys for xij in xi])
#     # # - from one node
#     # node_attr_all = g.nodes[list(g.nodes())[0]].keys()
#     node_attr_all = kn.utils.get_node_attribute_names(g)

#     # Check given node attributes name
#     if node_attr_list is not None:
#         if not np.all([k in node_attr_all for k in node_attr_list]):
#             raise ValueError('Attribute name does not exist, check `attr_node_list` parameters')
#     else:
#         node_attr_list = node_attr_all

#     # Check given mode
#     if isinstance(node_attr_mode, list):
#         if len(node_attr_mode) != len(node_attr_list):
#             raise ValueError('Length of the list `node_attr_mode` is not valid')

#         if np.any([s not in ('center', 'mean') for s in node_attr_mode]):
#             raise ValueError('Entry of the list `node_attr_mode` is not valid')

#     else: # `node_attr_mode` assumed to be a string
#         if node_attr_mode not in ('center', 'mean'):
#             raise ValueError('`node_attr_mode` is not valid')
#         node_attr_mode = len(node_attr_list)*[node_attr_mode]


#     # Step 1: Identifying the centers in the original graph
#     # ------

#     # (i) All nodes are marked as uncovered and non-centers.
#     nodes_uncovered_dict = {id:True for id in g.nodes()}
#     nodes_center_dict = {id:False for id in g.nodes()}

#     if edge_weight is None:
#         if verbose > 0:
#             print(f'MEMB (unweighted graph) - Computing shortest path length for any pair of nodes (with cutoff {r})...')
#         cost_from = dict(nx.all_pairs_shortest_path_length(g, cutoff=r))
#     else:
#         if verbose > 0:
#             print(f'MEMB (weighted graph) - Computing shortest path length for any pair of nodes (with cutoff {r})...')
#         cost_from = dict(nx.all_pairs_dijkstra_path_length(g, cutoff=r, weight=edge_weight))

#     n_nodes_uncovered = g.number_of_nodes() # or any positive integer
    
#     if verbose > 0:
#         print('MEMB - Selecting node centers...')
    
#     # Loop until all nodes are covered
#     n_centers = 0
#     while n_nodes_uncovered:
#         if verbose > 1:
#             print(f'... Number of centers: {n_centers:5d}, Number of uncovered nodes:  {n_nodes_uncovered:5d}')

#         # (ii) For all non-center nodes x (including already covered nodes),
#         # the excluded mass is computed, i.e. the number of uncovered nodes
#         # in the cluster centered at x
#         nodes_em_dict = {id:0 for id in g.nodes()}

#         for id in g.nodes():
#             if nodes_center_dict[id]:
#                 continue

#             nodes_em_dict[id] = len([j for j in cost_from[id].keys() if nodes_uncovered_dict[j]])

#         # print('...', 'nodes_em_dict', nodes_em_dict)

#         # Select a node c realizing the maximum of the excluded mass, and mark c as a center
#         em_max = np.asarray(list(nodes_em_dict.values())).max()
#         nodes_em_max_list = [k for k, v in nodes_em_dict.items() if v == em_max]
#         if len(nodes_em_max_list) > 1:
#             id = nodes_em_max_list[np.random.randint(len(nodes_em_max_list))]
#         else: # len(nodes_em_max_list) == 1
#             id = nodes_em_max_list[0]
#         nodes_center_dict[id] = True
#         n_centers = n_centers + 1

#         # print('...', 'new center', id)
        
#         # (iii) Mark all the nodes in the cluster centered at c (i.e. nodes 
#         # at a distance from c less than or equal to r) as covered
#         for j in cost_from[id].keys():
#             nodes_uncovered_dict[j] = False

#         # Update the number of uncovered nodes
#         n_nodes_uncovered = np.asarray(list(nodes_uncovered_dict.values())).sum()

#     if verbose > 0:
#         print(f'MEMB - Number of centers: {n_centers:d}')

#     # print('n_nodes_uncovered', n_nodes_uncovered)
#     # print('nodes_uncovered_dict', nodes_uncovered_dict)
#     # print('nodes_center_dict', nodes_center_dict)

#     # Step 2. Creating the clusters.
#     # -------
        
#     if verbose > 0:
#         print('MEMB - Assigning new id to node centers...')
    
#     # (i) Assign an id to each center.
#     nodes_center_list = [k for k, v in nodes_center_dict.items() if v]
#     nodes_new_id_dict = {i:-1 for i in g.nodes()}
#     for i, k in enumerate(nodes_center_list):
#         nodes_new_id_dict[k] = i

#     if verbose > 0:
#         print('MEMB - Computing distance to node centers for all nodes...')

#     # (ii) For all nodes, compute the central distance, i.e. the distance
#     # to the nearest center.    
#     nodes_central_dist_dict = {id:np.asarray([v for j, v in cost_from[id].items() if nodes_center_dict[j]]).min() for id in g.nodes()}

#     # (iii) Sort the nodes according to increasing central distance
#     nodes_central_dist_dict = {k:v for k, v in sorted(nodes_central_dist_dict.items(), key=lambda item: item[1])}

#     if verbose > 0:
#         print('MEMB - Assigning new id for non-center nodes...')

#     # (iv) For each non-center node u (visiting in order such that the central 
#     # distance is increasing), select one of its neighbors, v, with a smaller 
#     # central distance, and assign the id of v to u.
#     if edge_weight is None:
#         for id, d in nodes_central_dist_dict.items():
#             if nodes_center_dict[id]:
#                 continue
            
#             id_smaller_dist_list = [j for j in g.neighbors(id) if nodes_central_dist_dict[j] < d]

#             if len(id_smaller_dist_list) > 1:
#                 id_smaller_dist = id_smaller_dist_list[np.random.randint(len(id_smaller_dist_list))]
#             else: # len(id_smaller_dist_list) == 1
#                 id_smaller_dist = id_smaller_dist_list[0]
#             nodes_new_id_dict[id] = nodes_new_id_dict[id_smaller_dist]
    
#     else:
#         for id, d in nodes_central_dist_dict.items():
#             if nodes_center_dict[id]:
#                 continue
        
#             id_smaller_dist_list = [j for j in g.neighbors(id) if nodes_central_dist_dict[j] < d and nodes_central_dist_dict[j] + g.edges[id, j][edge_weight] <= r]

#             if len(id_smaller_dist_list) > 1:
#                 id_smaller_dist = id_smaller_dist_list[np.random.randint(len(id_smaller_dist_list))]
#             else: # len(id_smaller_dist_list) == 1
#                 id_smaller_dist = id_smaller_dist_list[0]
#             nodes_new_id_dict[id] = nodes_new_id_dict[id_smaller_dist]
    
#     # Set the list of original ids for each new id (in a dictionary)
#     # (returned optionally, and used for computing mean of node attributes over cluster (see below))
#     nodes_orig_id_dict = {i:[] for i in range(n_centers)}
#     for id in g.nodes():
#         nodes_orig_id_dict[nodes_new_id_dict[id]].append(id)

#     if verbose > 0:
#         print('MEMB - Building the reduced graph...')

#     # Step 3. Building the reduced graph.
#     # -------

#     g_red = nx.Graph()
#     # Set nodes (clusters in the original graph)
#     g_red.add_nodes_from(range(n_centers))
#     # Set edges (links between clusters in the original graph)
#     for u, v in g.edges():
#         if nodes_new_id_dict[u] != nodes_new_id_dict[v]:
#             g_red.add_edge(nodes_new_id_dict[u], nodes_new_id_dict[v])

#     # Set box (cluster) diameter if needed
#     if return_box_diameter_dict:
#        box_diameter_dict = {i: nx.diameter(nx.subgraph(g, nodes_orig_id_dict[i]), weight=edge_weight) for i in g_red.nodes()}

#     # Set node attributes in the reduced graph (according to `node_attr_list` and `node_attr_mode`)
#     for attr, mode in zip(node_attr_list, node_attr_mode):
#         d = nx.get_node_attributes(g, attr, default=np.nan) # dictionary original_id:value_of_attribute
#         if mode == 'center':
#             for c in nodes_center_list:
#                 g_red.nodes[nodes_new_id_dict[c]][attr] = d[c]
#         elif mode == 'mean':
#             v0 = list(d.values())[0] # value of one node in g
#             if hasattr(v0, '__len__'):
#                 attr_type = type(v0)[0]     # attribute type
#                 for i in g_red.nodes():
#                     # g_red.nodes[i][attr] = attr_type(np.nanmean(np.asarray([np.atleast_1d(d[j]) for j in nodes_orig_id_dict[i]]), axis=0))
#                     g_red.nodes[i][attr] = np.nanmean(np.asarray([np.atleast_1d(d[j]) for j in nodes_orig_id_dict[i]]), axis=0)
#             else:
#                 for i in g_red.nodes():
#                     g_red.nodes[i][attr] = np.nanmean(np.asarray([np.atleast_1d(d[j]) for j in nodes_orig_id_dict[i]]), axis=0)[0]

#     out = [g_red]
#     if return_center_dict:
#         out.append(nodes_center_dict)
#     if return_new_id_dict:
#         out.append(nodes_new_id_dict)
#     if return_orig_id_dict:
#         out.append(nodes_orig_id_dict)
#     if return_box_diameter_dict:
#         out.append(box_diameter_dict)
    
#     if len(out) == 1:
#         out = out[0]
#     else:
#         out = tuple(out)

#     return out
# # ----------------------------------------------------------------------------

# # ----------------------------------------------------------------------------
# def graph_conductance_exponent(
#         g, 
#         npairs,
#         nclasses=21,
#         class_in_log_scale=True,
#         fit_last_data_points_fraction=0.66,
#         confidence_level=0.95,
#         edge_weight=None, 
#         k=None, 
#         eigval=None, 
#         eigvec=None,
#         return_poly_fit_params=True,
#         return_poly_fit_params_cov=True,
#         return_starting_index_used_for_fit=True,
#         return_dist_mean=True,
#         return_weq_mean=True,
#         return_nb_points_per_class=True,
#         return_class_lim=True,
#         return_dist_full=False,
#         return_weq_full=False,
#         return_class_id_full=False,
#         return_pair_nodes_list=False,
#         seed=None,
#         verbose=1,
#         make_plot=True,
#         plot_in_log_scale=True,
#         plot_single_data_points=True,
#         **kwargs):
#     """
#     Computes the conductance exponent (also called resistance exponent) of a graph.

#     This function estimates the equivalent conductance (or resistance) for several 
#     pairs of nodes, randomly selected; then, assuming the relationship
    
#     .. math::
#         \omega_{eq} \propto r^{-d_{\Omega}}
    
#     or

#     .. math::
#         R_{eq} \propto r^{d_{\Omega}}

#     where :math:`\omega_{eq}` (resp. :math:`R_{eq}`) is the equivalent conductance
#     (resp. equivalent resistance) between two nodes, and :math:`r` the Euclidean 
#     distance between the two nodes, the exponent (called resistance exponent or 
#     conductance exponent) :math:`d_{\Omega}` is obtained by fitting 
#     a line on the log-log plot of :math:`\omega_{eq}` as function of :math:`r`.

#     Note: it is recommended to pre-compute (first) eigen values and eigen vectors
#     (and specify parameters `eigval` and `eigvec`).

#     Parameters
#     ----------
#     g : networkx.Graph
#         graph, must be connected (i.e. with one connected component)

#     npairs : int
#         number of pairs of nodes to select
    
#     nclasses : int, default: 21
#         number of classes for distances; let r(u, v) the Euclidean distance 
#         between the nodes u and v; let rmin and rmax the minimum and maximum
#         respectively of r(u,v) over all selected pairs of nodes (u, v); the 
#         interval [rmin, rmax] is divided into `nclasses` sub-intervals (classes)
#         of same length (in log scale if `class_in_log_scale=True`); then
#         each pair of nodes is assigned to one class according to the distance
#         between the two nodes; mean distance and mean equivalent conductance 
#         on each class will then be computed

#     class_in_log_scale : bool, default: True
#         - if `True`: the classes of distances are set in log scale
#         - if `False`: the classes of distances are set in ususal scale
        
#     fit_last_data_points_fraction : float, default : 0.66
#         positive number in (0, 1); sequences of mean distances and mean 
#         equivalent conductances (according to the classes determined by 
#         `nclasses`), constitute the data points; 
#         only the ending fraction `fit_last_data_points_fraction` 
#         of the data points are used to do the fit (line on the log-log plot);
#         this allows to only account for the "tail" of the data points

#     confidence_level : float, default: 0.95
#         confidence level, float in the interval (0, 1), to compute the 
#         uncertainty for the line fitting (and for the conductance exponent)

#     edge_weight : str, optional
#         name of the edge attribute used as weight, for conductance; 
#         unused if `eigval` and `eigvec` are specified (not `None`)
#         by default (`None`) : all edges have a weight of 1

#     k : int, optional
#         number of eigen values and eigen vectors to use, the first
#         `k` eigen values in ascending order are considered;
#         unused if `eigval` and `eigvec` are specified (not `None`);
#         by default (`None`): all eigen values and eigen vectors are 
#         computed, i.e. `k` is set to `g.number_of_nodes()`

#     eigval : 1d-array, optional
#         array of shape (m, ), with m less than or equal to the number 
#         of nodes in the graph, m first eigen values of the laplacian 
#         matrix, in ascending order, i.e. for a connected graph, 
#         `0 = eigval[0] < eigenval[1] <= eigenval[2] ...`;
#         by default (`None`): the first `k` eigen values are computed

#     eigvec : 2d-array, optional
#         array of shape (n, m), where n is the number of nodes in the
#         graph and m <= n, m first eigen vectors in columns (orthogonal 
#         and of norm 1) corresponding to the eigen values `eigval`; 
#         for a connected graph, `eigvec[:,0]` is proportional to the 
#         vector `[1, 1, ..., 1]`;
#         by default (`None`): the first `k` eigen vectors are computed

#     return_poly_fit_params : bool, default: True
#         if `True`, the array (shape (2, )) of parameters of the fitted 
#         polynom (line) is returned

#     return_poly_fit_params_cov : bool, default: True
#         if `True`, the array (shape (2, 2)) of parameters covariance of the 
#         fitted polynom (line) is returned
    
#     return_starting_index_used_for_fit: bool, default: True
#         if `True`, the starting index of the sequence of data points
#         used for fitting is returned

#     return_dist_mean : bool, default: True
#         if `True`, the array of mean distances (in each class) is returned
#         (mean over all pair of nodes in one class, for each class)

#     return_weq_mean : bool, default: True
#         if `True`, the array of mean equivalent conductances (in each class) 
#         is returned (mean over all pair of nodes in one class, for each class)

#     return_nb_points_per_class : bool, default: True
#         if `True`, number of points (pairs) per class is returned

#     return_class_lim : bool, default: True
#         if `True`, array of class bounds is returned

#     return_dist_full : bool, default: False
#         if `True`, the array of all distances is returned

#     return_weq_full : bool, default: False
#         if `True`, the array of all equivalent conductances is returned

#     return_class_id_full : bool, default: False
#         if `True`, the array of class id of all points (pairs) is returned

#     return_pair_nodes_list : bool, default: False
#         if `True`, the list of selected pairs of nodes (labels) is returned

#     seed : int, optional
#         seed number for initializing the random number generator

#     verbose : int, default: 1
#         verbose mode, larger value implies more info printed

#     make_plot : bool, default: True
#         indicates if the plot of the data points and fitting curve is displayed 
#         (in the current figure axis)

#     plot_in_log_scale : bool, default: True
#         used if `make_plot=True`, indicates if the plot is in log-log scale
#         (along x- and y-axes)

#     plot_single_data_points : bool, default: True
#         used if `make_plot=True`, indicates if single data points (pairs) are
#         displayed on the plot

#     kwargs : dict
#         keyword arguments passed to the function 
#         `scipy.sparse.linalg.eigsh` to compute the `k` first eigen values
#         and eigen vectors (unused if `k >= g.number_of_nodes()`, in this 
#         case, the function `scipy.linalg.eigh` is used); 
#         note: the parameter `which`, if not given in `kwargs`, is set to 
#         'SM' (smallest eigen values are considered) 

#     Returns
#     -------
#     out_dict : dict
#         dictionary with the following keys / values:
        
#         - dOmega : float
#             conductance exponent (or resistance exponent)
    
#         - dOmega_delta : float
#             uncertainty for conductance exponent, the interval 
#             `dOmega +/- dOmega_delta` corresponds to the confidence interval 
#             at a confidence level `confidence_level` derived from Gaussian 
#             distribution

#         - poly_fit_params : array of shape (2, ), optional
#             returned if `return_poly_fit_params=True`:
#             array of parameters of the fitted polynom (line); in particular:
#             `dOmega = -poly_fit_params[0]`

#         - poly_fit_params_cov : array of shape (2, 2), optional
#             returned if `return_poly_fit_params_cov=True`:
#             array (shape (2, 2)) of parameters covariance of the fitted polynom 
#             (line); in particular: with
#             `dOmega_delta = scipy.stats.norm.ppf((1+confidence_level)/2) * np.sqrt(poly_fit_params_cov[0, 0])`
        
#         - dist_mean : 1d-array of shape (nclasses,), optional
#             returned if `return_dist_mean=True`:
#             array of mean on each class of Euclidean distances between the 
#             two nodes of the selected pair of nodes

#         - weq_mean : 1d-array of shape (nclasses,), optional
#             returned if `return_weq_mean=True`:
#             array of mean on each class of equivalent conductance for the
#             selected pair of nodes

#         - nb_points_per_class : 1d-array of shape (nclasses,), optional
#             returned if `return_nb_points_per_class=True`:
#             array of number of points (pairs) in each class

#         - class_lim : 1d-array of shape (nclasses+1,), optional
#             returned if `return_class_lim=True`:
#             array of class limits: increasing float numbers determining the
#             classes

#         - dist_full : 1d-array of shape (npairs,), optional
#             returned if `return_dist_full=True`:
#             array of Euclidean distances between the two nodes of all selected 
#             pair of nodes

#         - weq_full : 1d-array of shape (npairs,), optional
#             returned if `return_weq_full=True`:
#             array of equivalent conductances between the two nodes of all selected 
#             pair of nodes (corresponding to distances in `dist_full`)

#         - class_id_full : 1d-array of shape (npairs, ), optional
#             returned if `return_class_id_full=True`:
#             array of class id (int in (0, ..., nclass-1)) of all points (pairs)

#         - pair_nodes_list : list of length npairs, optional
#             returned if `return_pair_nodes_list=True`:
#             list of selected pairs of nodes (labels) in the graph `g`:

#             - `pair_nodes_list[i] = [u1, u2]`, where `u1`, `u2` are two \
#             nodes (labels) in the graph
#     """
#     if seed is not None:
#         np.random.seed(seed)
    
#     if nclasses < 1:
#         print('No class!')
#         return
    
#     if npairs < nclasses:
#         print(f'Two few pairs ({npairs}, less than number of classes ({nclasses}))')
#         return
    
#     # Set dictionary to convert node index to node label
#     node_index2label = {i:u for i, u in enumerate(g.nodes())}
#     # node_index2label = kn.utils.get_node_index2label(g) # equivalent

#     n_nodes = g.number_of_nodes()
#     pair_index_list = [np.random.choice(n_nodes, size=2, replace=False) for _ in range(npairs)]
#     pair_nodes_list = [(node_index2label[i], node_index2label[j]) for i, j in pair_index_list]
 
#     # Compute the Euclidean distance between the two nodes for each pair
#     dist_full = np.asarray([np.asarray(g.nodes[x]['pos']) - np.asarray(g.nodes[y]['pos']) for x, y in pair_nodes_list])
#     dist_full = np.sqrt(np.sum(dist_full**2, axis=1))
    
#     # Compute the equivalent conductance between the two nodes for each pair
#     if verbose > 0:
#         if eigvec is None or eigval is None:
#             print(f'Compute Laplacian pseudo inverse and equivalent conductance ({npairs} pair of nodes)')
#         else:
#             print(f'Compute equivalent conductance ({npairs} pair of nodes) (Laplacian pseudo inverse pre-computed)')
    
#     t1_all = time.time()
#     weq_full = equivalent_conductance(g, pair_nodes_list, edge_weight=edge_weight, k=k, eigval=eigval, eigvec=eigvec, **kwargs)
#     t2_all = time.time()
#     if verbose > 0:
#         print(f'Total elapsed time (w_eq for pair of nodes): {t2_all-t1_all:.2f} sec')

#     if verbose > 0:
#         print(f'Divide distances into classes, and compute mean on each class (of distance and of equivalent conductance)')
    
#     rmin = dist_full.min()
#     rmax = dist_full.max()
#     if class_in_log_scale:
#         class_lim = np.linspace(np.log10(rmin), np.log10(rmax), nclasses+1)
#         class_id_full = np.digitize(np.log10(dist_full), class_lim) - 1
#     else:
#         class_lim = np.linspace(rmin, rmax, nclasses+1)
#         class_id_full = np.digitize(dist_full, class_lim) - 1

#     nb_points_per_class = np.asarray([np.sum(class_id_full == i) for i in range(nclasses)])

#     # dist_mean = np.asarray([np.mean(dist_full[class_id_full == i]) for i in range(nclasses)])
#     # weq_mean = np.asarray([np.mean(weq_full[class_id_full == i]) for i in range(nclasses)])

#     dist_mean = np.full(nclasses, np.nan)
#     weq_mean = np.full(nclasses, np.nan)
#     for i in range(nclasses):
#         if nb_points_per_class[i] > 0:
#             dist_mean[i] = np.mean(dist_full[class_id_full==i])    
#             weq_mean[i] = np.mean(weq_full[class_id_full==i])    
    
#     if verbose > 0:
#         print(f'Compute conductance exponent (line fitting on log-log plot of equivalent conductance as function of distance)')

#     # Compute conductance exponent (line fitting on log-log plot of equivalent conductance as function of distance)
#     # -------------------------------------------------------------------------------------------------------------
#     # Data (all)
#     x_all = dist_mean
#     y_all = weq_mean

#     n_all = len(x_all)
#     if n_all < 2:
#         raise ValueError(f"Not enough points ({n_all}) for line fitting in log-log plot ...")

#     # Coordinates used for fitting
#     n = min(n_all - int(fit_last_data_points_fraction * n_all), n_all - 3) # starting index (at least 3 points)
#     if verbose > 0:
#         print(f'Fitting on last {n_all - n} data points over {n_all}')
#     x = x_all[n:]
#     y = y_all[n:]

#     ind = ~np.any((np.isnan(x), np.isnan(y)), axis=0)
#     x = x[ind]
#     y = y[ind]
                 
#     # Poly-fit
#     poly_fit_params, poly_fit_params_cov = np.polyfit(np.log10(x), np.log10(y), 1, cov=True)

#     # Get slope and slope_delta for given confidence level
#     t = scipy.stats.norm.ppf((1+confidence_level)/2)
#     slope = poly_fit_params[0]
#     slope_delta = t * np.sqrt(poly_fit_params_cov[0, 0])

#     # Get conductance exponent dOmega, and dOmega_delta
#     dOmega = - slope
#     dOmega_delta = slope_delta

#     # Make plot (if needed)
#     if make_plot:
#         dOmega_lim_min = dOmega - dOmega_delta
#         dOmega_lim_max = dOmega + dOmega_delta
        
#         # Estimation with the fitted line
#         xx = np.linspace(x.min(), x.max(), 200)
#         yy = np.power(10, poly_fit_params[0] * np.log10(xx) + poly_fit_params[1])

#         xx_all = np.linspace(np.nanmin(x_all), np.nanmax(x_all), 200)
#         yy_all = np.power(10, poly_fit_params[0] * np.log10(xx_all) + poly_fit_params[1])

#         # Plot parameters
#         xlabel = '$r$ : distance'
#         ylabel = '$\omega_{eq}$ : equivalent conductance'

#         title = '$\omega_{eq}\propto r^{-d_{\Omega}}$' + ', $d_{\Omega}$' + f'={dOmega:.5f} in [{dOmega_lim_min:.5f}, {dOmega_lim_max:.5f}] [{100*confidence_level:.2g}% - inter.]'

#         color_data_all      = 'tab:blue'
#         marker_data_all     = 'o'
#         markersize_data_all = 5

#         color_data_used      = 'tab:orange'
#         marker_data_used     = 'o'
#         markersize_data_used = 5

#         color_fit_data_used = 'red'
#         ls_fit_data_used    = 'solid'
#         lw_fit_data_used    = 1.

#         color_fit_data_all = 'red'
#         ls_fit_data_all    = 'dotted'
#         lw_fit_data_all    = 1.

#         # Figure
#         # ------
#         # plt.figure(figsize=(10, 6))
#         if plot_single_data_points:
#             # Add single data points (pair) (before computing mean on classes)
#             marker_single     = '.'
#             markersize_single = 1

#             cmap = plt.get_cmap('viridis')

#             label_single = 'single data point (pair)'
#             for i in range(nclasses):
#                 plt.plot(dist_full[class_id_full==i], weq_full[class_id_full==i], alpha=.3, ls='', marker=marker_single, markersize=markersize_single, color=cmap(i/(nclasses-1)), label=label_single)
#                 label_single = None

#             # # or:
#             # color_single_all      = 'tab:blue'
#             # marker_single_all     = '.'
#             # markersize_single_all = 1

#             # color_single_used      = 'tab:orange'
#             # marker_single_used     = '.'
#             # markersize_single_used = 1

#             # label_single_used = 'single pair (used)'
#             # label_single_all = 'single pair'
#             # for i in range(nclasses):
#             #     if i >= n1 and i < n2:
#             #         plt.plot(dist_full[class_id_full==i], weq_full[class_id_full==i], alpha=.5, ls='', marker=marker_single_used, markersize=markersize_single_used, color=color_single_used, label=label_single_used)
#             #         label_single_used = None
#             #     else:
#             #         plt.plot(dist_full[class_id_full==i], weq_full[class_id_full==i], alpha=.5, ls='', marker=marker_single_all , markersize=markersize_single_all , color=color_single_all , label=label_single_all)
#             #         label_single_all = None
        
#         plt.plot(x_all, y_all, ls='', marker=marker_data_all , markersize=markersize_data_all , color=color_data_all , label='data')
#         plt.plot(x    , y    , ls='', marker=marker_data_used, markersize=markersize_data_used, color=color_data_used, label='data used for fitting')
#         plt.plot(xx_all, yy_all, ls=ls_fit_data_all , lw=lw_fit_data_all , color=color_fit_data_all)
#         plt.plot(xx    , yy    , ls=ls_fit_data_used, lw=lw_fit_data_used, color=color_fit_data_used, label='fit')
        
#         # Add number of points in each class (text)
#         for xi, yi, ni in zip(x_all, y_all, nb_points_per_class):
#             if not np.isnan(xi) and not np.isnan(yi):
#                 plt.text(xi, yi, f'{ni}', ha='left', va='bottom')
        
#         if plot_in_log_scale:
#             plt.xscale('log')
#             plt.yscale('log')
        
#         plt.xlabel(xlabel)
#         plt.ylabel(ylabel)
#         plt.title(title)
#         plt.grid()
#         plt.legend()
#         # plt.show()

#     # Set output dictionary
#     out_dict = {'dOmega': dOmega, 'dOmega_delta': dOmega_delta}

#     if return_poly_fit_params:
#         out_dict['poly_fit_params'] = poly_fit_params

#     if return_poly_fit_params_cov:
#         out_dict['poly_fit_params_cov'] = poly_fit_params_cov

#     if return_starting_index_used_for_fit:
#         out_dict['starting_index_used_for_fit'] = n

#     if return_dist_mean:
#         out_dict['dist_mean'] = dist_mean

#     if return_weq_mean:
#         out_dict['weq_mean'] = weq_mean

#     if return_nb_points_per_class:
#         out_dict['nb_points_per_class'] = nb_points_per_class

#     if return_class_lim:
#         out_dict['class_lim'] = class_lim

#     if return_dist_full:
#         out_dict['dist_full'] = dist_full

#     if return_weq_full:
#         out_dict['weq_full'] = weq_full

#     if return_class_id_full:
#         out_dict['class_id_full'] = class_id_full

#     if return_pair_nodes_list:
#         out_dict['pair_nodes_list'] = pair_nodes_list

#     return out_dict
# # ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
# def _loglog_fit(x_all, y_all, fit_last_data_points_fraction, verbose=0):
#     """Internal (private) function to fit line on a log log plot.
#     """
#     n_all = len(x_all)
    
#     # Coordinates used for fitting
#     n = int(np.round((1.0 - fit_last_data_points_fraction) * n_all))

#     x = x_all[n:]
#     y = y_all[n:]

#     # ... extract non nan number
#     ind = ~np.any((np.isnan(x), np.isnan(y)), axis=0)
#     if not ind.all():
#         if verbose > 0:
#             print(f'WARNING: undefined point (nan) encountered (removed)')

#         x = x[ind]
#         y = y[ind]

#     n_fit = len(x) # number of data points used for fitting

#     if n_fit < 2:
#         # No fit
#         if verbose > 0:
#             print(f'WARNING: not enough points ({n_fit}) for line fitting in log-log plot ; try to increase `fit_last_data_points_fraction`"')

#         poly_fit_params, poly_fit_params_cov = None, None

#     elif n_fit == 2:
#         # Fit without uncertainty  
#         if verbose > 0:
#             print(f'WARNING: fitting on 2 points in log-log plot ; try to increase `fit_last_data_points_fraction`"')
        
#         try:
#             poly_fit_params = np.polyfit(np.log10(x), np.log10(y), 1, cov=False) # cov=False, i.e. no uncertainty
#             poly_fit_params_cov = None

#         except:
#             if verbose > 0:
#                 print(f'WARNING: fitting failed [on {n_fit} points (over {n_all})]')
            
#             poly_fit_params, poly_fit_params_cov = None, None

#     else:
#         # Poly-fit
#         if verbose > 0:
#             print(f'Fitting on {n_fit} points (over {n_all})')

#         try:
#             poly_fit_params, poly_fit_params_cov = np.polyfit(np.log10(x), np.log10(y), 1, cov=True)

#         except:
#             if verbose > 0:
#                 print(f'WARNING: fitting failed [on {n_fit} points (over {n_all})]')
            
#             poly_fit_params, poly_fit_params_cov = None, None

#     return poly_fit_params, poly_fit_params_cov, x, y
# # ----------------------------------------------------------------------------

# # ----------------------------------------------------------------------------
# def remove_simple_cycles(
#         g, 
#         max_cycle_len=None,
#         n_iter_max = None,
#         pos_attr='pos',
#         verbose=0):
#     """
#     Simplifies a graph by removing simple cycles.

#     The function `networkx.simple_cyles` is used to identify the cycles.
#     The nodes of any cycle whose length (in number of nodes) not exceeding 
#     `max_cylcle_len` (if specified) are removed and replaced by one node, 
#     and the edges issued from the removed nodes are removed and replaced 
#     by new edges issued from the new node. The position of the new nodes is 
#     set to the mean of the position of the nodes that it replaces.

#     The procedure is iterative: at each iteration, only cycles that
#     don't share common nodes are treated.

#     Notes: 
    
#     - edge attributes and node attributes (other than the position) \
#     are not considered for new edges and new nodes.
#     - the node labels in the output graph are converted to integers \
#     starting from 0
        
#     Parameters
#     ----------
#     g : networkx.Graph
#         input graph

#     max_cycle_len : int, optional
#         maximum cycle length (in number of nodes): cycles with more nodes
#         are kept;
#         by default (`None`): cycle length not limited
    
#     n_iter_max : int, optional
#         maximal number of iteration(s)
#         by default (`None`): as many iterations as needed to remove all 
#         cycles (not exceeding a length of `max_cycle_len`)

#     pos_attr : str, default: 'pos'
#         name of the node attribute for position

#     verbose : int, default: 0
#         verbose mode, larger value implies more info printed

#     Returns
#     -------
#     g_out : networkx.Graph
#         output graph (simplifed as described above)
#     """

#     # Copy graph to initialize output graph (g_out)
#     g_out = g.copy()
    
#     # Convert node labels to integers from 0 to n_nodes
#     g_out = nx.convert_node_labels_to_integers(g_out)

#     if n_iter_max is None:
#         n_iter_max = g_out.number_of_nodes()

#     for n_iter in range(n_iter_max):
#         # Number of nodes
#         n_nodes = g_out.number_of_nodes()

#         # Get simple cycles to be removed
#         if max_cycle_len is not None:
#             cycle_nodes_list = [c for c in nx.simple_cycles(g_out) if len(c) <= max_cycle_len]
#         else:
#             cycle_nodes_list = [c for c in nx.simple_cycles(g_out)]

#         # Remove cycles so that two cycles do not have common nodes
#         i0 = 0
#         while i0 < len(cycle_nodes_list) - 1:
#             ind = np.ones(len(cycle_nodes_list), dtype='bool')
#             for i in range(i0+1, len(cycle_nodes_list)):
#                 if len(np.intersect1d(cycle_nodes_list[i0], cycle_nodes_list[i])) > 0:
#                     ind[i] = False
#             cycle_nodes_list = [cycle_nodes_list[i] for i in np.where(ind)[0]]
#             i0 = i0+1
            
#         # Number of new nodes
#         n_new_nodes = len(cycle_nodes_list)

#         if verbose > 0:
#             print(f'Iteration {n_iter+1} : number of cycles to simplify = {n_new_nodes}')

#         if n_new_nodes == 0:
#             return g_out

#         # Set list of node to remove
#         nodes_list_to_remove = [u for c in cycle_nodes_list  for u in c]

#         # Set position of new nodes: cycle "i" becomes new node "n_nodes+i"
#         new_node_id = range(n_nodes, n_nodes + n_new_nodes)

#         # Get node position of new nodes (mean of nodes in the cycle)
#         new_node_pos = np.asarray([np.mean(np.asarray([np.asarray(g_out.nodes[u][pos_attr]) for u in c]), axis=0) for c in cycle_nodes_list])

#         # Set dictionary with
#         # - key: current node id
#         # - value: list of the cycle ids containing the node
#         node2cycle = {k:[] for k in range(n_nodes)}
#         for i, c in enumerate(cycle_nodes_list):
#             for u in c:
#                 node2cycle[u].append(i)
#         # -> node u s.t. node2cycle[u] is a non empty list will be removed

#         # Get new edges according to the edges involving nodes that will be removed
#         new_edges = np.empty((0, 2), dtype='int')
#         for i, c in enumerate(cycle_nodes_list):
#             for u in c:
#                 edges = []
#                 for v in g_out.neighbors(u):
#                     if len(node2cycle[v]) == 0:
#                         edges.append((n_nodes + i, v))
#                     else:
#                         for j in node2cycle[v]:
#                             if j != i:
#                                 edges.append((n_nodes + min(i, j), n_nodes + max(i, j)))
#                 if len(edges):
#                     edges = np.unique(edges, axis=0)
#                     new_edges = np.unique(np.vstack((new_edges, edges)), axis=0)
    
#         # Remove nodes in cycle
#         g_out.remove_nodes_from(nodes_list_to_remove)

#         # Add new nodes with position 
#         g_out.add_nodes_from(new_node_id)
#         nx.set_node_attributes(g_out, {k:v for k, v in zip(new_node_id, new_node_pos)}, pos_attr)

#         # Add new edges
#         g_out.add_edges_from(new_edges)
        
#         # Convert node labels to integers
#         g_out = nx.convert_node_labels_to_integers(g_out)

#     return g_out    
# # ----------------------------------------------------------------------------

# # ----------------------------------------------------------------------------
# def remove_simple_cycles(
#         g, 
#         max_cycle_len=None,
#         n_iter_max = None,
#         node_attr_mode='mean',
#         node_attr_list=None,
#         return_node_merging_count_dict=False,
#         return_node_merging_list_dict=False,
#         verbose=0):
#     """
#     Simplifies a graph by removing simple cycles.

#     The function `networkx.simple_cyles` is used to identify the cycles.
#     The nodes of any cycle whose length (in number of nodes) not exceeding 
#     `max_cylcle_len` (if specified) are removed and replaced by one node, 
#     and the edges issued from the removed nodes are removed and replaced 
#     by new edges issued from the new node.

#     The procedure is iterative: at each iteration, only cycles that
#     don't share common nodes are treated.

#     At each iteration, a new node has the label of the first node in the 
#     cycle it comes from, and its attributes are computed according to the 
#     attributes of the nodes in the cycle it comes from, using the mode(s) 
#     specified in `node_attr_mode`. In particular, the set of output nodes
#     (labels) is a subset of the set of input nodes (labels).

#     The edge attributes are not considered for new edges.

#     The function returns an output graph (simplified as described above)
#     and an output dictionary where the keys are the node labels in the 
#     input graph, and the values the number of nodes in the input graph that
#     have been merged to form this node in the output graph. (A value of 0 
#     indicates that the node has been removed in the output graph.)

#     Note: this function should not be used if too many cycles are present
#     (function `nx.simple_cycles` can be very time consuming in this case).

#     Parameters
#     ----------
#     g : networkx.Graph
#         input graph

#     max_cycle_len : int, optional
#         maximum cycle length (in number of nodes): cycles with more nodes
#         are kept;
#         by default (`None`): cycle length not limited
    
#     n_iter_max : int, optional
#         maximal number of iteration(s)
#         by default (`None`): as many iterations as needed to remove all 
#         cycles (not exceeding a length of `max_cycle_len`)

#     node_attr_mode : str {'first', 'mean'} (or list of strs), default: 'mean'
#         string or list of strings, the strings indicate how the node attributes
#         (for nodes that are part of a cycle) are computed at each iteration
#          of the  in the output graph (see above)

#         - if string: the same operation (mode) is used for all considered \
#         node attributes
#         - if list of strings: the length of the list must be equal to the \
#         the list `node_attr_list`, and `node_attr_mode[i]` indicates the \
#         operation (mode) used for the node attribute `node_attr_list[i]`

#     node_attr_list : list of strs, optional
#         list of node attributes of the input graph to be included in 
#         the output graph;
#         by default (`None`): the list of all nodes attributes is considered

#     return_node_merging_count_dict : bool, default: False
#         if `True`, the dictionary `node_merging_count_dict` is returned
#         (see below)

#     return_node_merging_list_dict : bool, default: False
#         if `True`, the dictionary `node_merging_list_dict` is returned
#         (see below)

#     verbose : int, default: 0
#         verbose mode, larger value implies more info printed

#     Returns
#     -------
#     g_out : networkx.Graph
#         output graph (simplifed as described above)

#     node_merging_count_dict : dict
#         return if `return_node_merging_count_dict=True`; dictionary with

#         - key: node label in the input graph
#         - value:
#             - 0 if the node has been removed in the output graph
#             - k > 0 if the node is in the output graph, k being the number of \
#             nodes in the input graph that have been merged to form this node \
#             in the output graph

#     node_merging_list_dict : dict, optional
#         return if `return_node_merging_list_dict=True`; dictionary with

#         - key: node label in the input graph
#         - value: the list of nodes in the input graph that have been merged to \
#         form this node in the output graph
#     """

#     # Check node attributes of the input graph to be considered in the output graph
#     # and corresponding operation (mode)

#     # Get keys (attribute names) in the input graph
#     # # - all keys
#     # keys = [list(gg.nodes[i].keys()) for i in gg.nodes()] # list of lists
#     # node_attr_all = np.unique([xij for xi in keys for xij in xi])
#     # # - from one node
#     # node_attr_all = gg.nodes[list(gg.nodes())[0]].keys()
#     node_attr_all = kn.utils.get_node_attribute_names(g)

#     # Check given node attributes name
#     if node_attr_list is not None:
#         if not np.all([k in node_attr_all for k in node_attr_list]):
#             raise ValueError('Attribute name does not exist, check `attr_node_list` parameters')
#     else:
#         node_attr_list = node_attr_all

#     # Check given mode
#     if isinstance(node_attr_mode, list):
#         if len(node_attr_mode) != len(node_attr_list):
#             raise ValueError('Length of the list `node_attr_mode` is not valid')

#         if np.any([s not in ('first', 'mean') for s in node_attr_mode]):
#             raise ValueError('Entry of the list `node_attr_mode` is not valid')

#     else: # `node_attr_mode` assumed to be a string
#         if node_attr_mode not in ('first', 'mean'):
#             raise ValueError('`node_attr_mode` is not valid')
#         node_attr_mode = len(node_attr_list)*[node_attr_mode]

#     # Initialize `node_merging_count_dict` and `node_merging_list_dict`
#     node_merging_count_dict = {u:1 for u in g.nodes()} # used even if not returned
    
#     if return_node_merging_list_dict:
#         node_merging_list_dict = {u:[u] for u in g.nodes()}

#     # Copy graph to initialize output graph (g_out)
#     g_out = g.copy()
    
#     # # Convert node labels to integers from 0 to n_nodes
#     # g_out = nx.convert_node_labels_to_integers(g_out)

#     if n_iter_max is None:
#         n_iter_max = g_out.number_of_nodes()

#     for n_iter in range(n_iter_max):
#         # Get simple cycles to be removed
#         if max_cycle_len is not None:
#             cycle_nodes_list = [c for c in nx.simple_cycles(g_out) if len(c) <= max_cycle_len]
#         else:
#             cycle_nodes_list = [c for c in nx.simple_cycles(g_out)]

#         # Remove cycles so that two cycles do not have common nodes
#         i0 = 0
#         while i0 < len(cycle_nodes_list) - 1:
#             ind = np.ones(len(cycle_nodes_list), dtype='bool')
#             for i in range(i0+1, len(cycle_nodes_list)):
#                 for x in cycle_nodes_list[i]:
#                     if x in cycle_nodes_list[i0]:
#                         ind[i] = False
#                         break
#             cycle_nodes_list = [cycle_nodes_list[i] for i in np.where(ind)[0]]
#             i0 = i0+1
            
#         # Number of new nodes
#         n_new_nodes = len(cycle_nodes_list)

#         if verbose > 0:
#             print(f'Iteration {n_iter+1} : number of cycles to simplify = {n_new_nodes}')

#         if n_new_nodes == 0:
#             break # exit for loop

#         # Set list of nodes to remove (all except first node in each cycle)
#         nodes_list_to_remove = [u for c in cycle_nodes_list for u in c[1:]]

#         # Set new (or kept) node labels (first node in each cycle)
#         new_node_labels = [c[0] for c in cycle_nodes_list]

#         # Update `node_merging_list_dict` if needed
#         if return_node_merging_list_dict:
#             for c in cycle_nodes_list:
#                 for u in c[1:]:
#                     node_merging_list_dict[c[0]].extend(node_merging_list_dict[u])
#                     node_merging_list_dict[u] = []

#         # Set node attributes for new nodes (according to `node_attr_list` and `node_attr_mode`)
#         for attr, mode in zip(node_attr_list, node_attr_mode):
#             d = nx.get_node_attributes(g_out, attr, default=np.nan) # dictionary node_label:value_of_attribute
#             if mode == 'first':
#                 pass
#                 # for u in new_node_labels:
#                 #     g_out.nodes[u][attr] = d[u]
#             elif mode == 'mean':
#                 v0 = list(d.values())[0] # value of one node in g
#                 if hasattr(v0, '__len__'):
#                     for u, c in zip(new_node_labels, cycle_nodes_list):
#                         values = np.asarray([np.atleast_1d(d[v]) for v in c])
#                         w = np.asarray([node_merging_count_dict[v] for v in c]).reshape(-1, 1)
#                         ind = np.all(~np.isnan(values), axis=1)
#                         g_out.nodes[u][attr] = np.sum(values[ind]*w[ind], axis=0)/np.sum(w[ind])
#                 else:
#                     for u, c in zip(new_node_labels, cycle_nodes_list):
#                         values = np.asarray([np.atleast_1d(d[v]) for v in c])
#                         w = np.asarray([node_merging_count_dict[v] for v in c]).reshape(-1, 1)
#                         ind = np.all(~np.isnan(values), axis=1)
#                         g_out.nodes[u][attr] = (np.sum(values[ind]*w[ind], axis=0)/np.sum(w[ind]))[0]
#                 # if hasattr(v0, '__len__'):
#                 #     for u, c in zip(new_node_labels, cycle_nodes_list):
#                 #         g_out.nodes[u][attr] = np.nanmean(np.asarray([np.atleast_1d(d[v]) for v in c]), axis=0)
#                 # else:
#                 #     for u, c in zip(new_node_labels, cycle_nodes_list):
#                 #         g_out.nodes[u][attr] = np.nanmean(np.asarray([np.atleast_1d(d[v]) for v in c]), axis=0)[0]

#         # Update `node_merging_count_dict`
#         for c in cycle_nodes_list:
#             node_merging_count_dict[c[0]] = np.asarray([node_merging_count_dict[u] for u in c]).sum()
#             for u in c[1:]:
#                 node_merging_count_dict[u] = 0

#         # Set dictionary with
#         # - key: current node label
#         # - value: list of the cycle index containing the node
#         node2cycle = {u:[] for u in g_out.nodes()}
#         for i, c in enumerate(cycle_nodes_list):
#             for u in c:
#                 node2cycle[u].append(i)
#         # -> node u s.t. node2cycle[u] is a non empty list will be removed

#         # Get new edges according to the edges involving nodes that will be removed
#         new_edges = []
#         for i, c in enumerate(cycle_nodes_list):
#             for u in c[1:]:
#                 edges = []
#                 for v in g_out.neighbors(u):
#                     if len(node2cycle[v]) == 0:
#                         edges.append((new_node_labels[i], v))
#                     else:
#                         for j in node2cycle[v]:
#                             if j != i:
#                                 edges.append((new_node_labels[i], new_node_labels[j]))
#                 for e in edges:
#                     if [e[0], e[1]] not in new_edges and [e[1], e[0]] not in new_edges:
#                         new_edges.append(e)

#         # Remove nodes in cycle
#         g_out.remove_nodes_from(nodes_list_to_remove)

#         # Add new edges
#         g_out.add_edges_from(new_edges)
                
#     out = [g_out]
#     if return_node_merging_count_dict:
#         out.append(node_merging_count_dict)
#     if return_node_merging_list_dict:
#         out.append(node_merging_list_dict)
    
#     if len(out) == 1:
#         out = out[0]
#     else:
#         out = tuple(out)
    
#     return out
# # ----------------------------------------------------------------------------

# # ----------------------------------------------------------------------------
# # ===== OLD =====
# def solve_flow_with_bc(
#         G,
#         pB_nodes_dict,
#         qI_nodes_dict,
#         edge_weight=None,
#         return_entering_fluxes=True):
#     """
#     Solves the flow problem in a connected graph, with boundary condition (potential).

#     Computes the potential at each node in the graph, given the boundary
#     potential and the source term (entering fluxes) at the inner nodes.

#     The potential at the boundary nodes is given by `pB_nodes_dict`, and the
#     source term (entering fluxes) at the inner nodes is given by `qI_nodes_dict`.
#     The boundary nodes B (keys of `pB_nodes_dict`) and the inner nodes I (keys of 
#     `qI_nodes_dict`) must form a partition of all graph nodes (i.e. disjoint sets 
#     with union equal to the set of all graph nodes). Moreover, the boundary nodes B 
#     must be a non empty set, and the subgraph induced by the inner nodes I must be 
#     connected (unless I is empty).

#     The edge weights represent the conductances of the edges; if not specified
#     (`edge_weight=None`), all edges have a conductance of 1.

#     Denoting :math:`L` is the Laplacian matrix of the graph (with given edge 
#     weights), the potential at the inner nodes are given by 
    
#     .. math::
#         p_I = {L_{II}}^{-1}\\cdot (q_I - L_{IB} p_B)
    
#     where the submatrix :math:`L_{II}` is invertible provided that the subgraph
#     induced by the nodes I is connected, and the entering fluxes at the boundary 
#     nodes are given by 

#     .. math::
#         q_B = L_{BI}\\cdot p_I + L_{BB} p_B
    
#     Parameters
#     ----------
#     G : networkx.Graph
#         graph, must be connected (i.e. with one connected component)

#     pB_nodes_dict : dict
#         dictionary of potential at boundary nodes, the keys are the node
#         labels and the values are the potential values

#     qI_nodes_dict : dict
#         dictionary of entering flux at inner nodes, the keys are the node
#         labels and the values are the entering flux values
            
#     edge_weight : str, optional
#         name of the edge attribute used as weight (conductance)
    
#     return_entering_fluxes : bool, default: True
#         if `True`, the entering flux at all nodes in the graph is returned

#     Returns
#     -------
#     p_nodes_dict : dict
#         dictionary of potential at all nodes in the graph, the keys are the node
#         labels and the values are the potential values
    
#     q_nodes_dict : dict, optional
#         if `return_entering_fluxes=True`, dictionary of entering flux at all nodes 
#         in the graph, the keys are the node labels and the values are the entering 
#         flux values
#     """
#     # Set dictionary to convert node label to node index
#     node_label2index = {u:i for i, u in enumerate(G.nodes())}
#     # node_index2label = kn.utils.get_node_label2index(G) # equivalent

#     B_nodes_label_list = pB_nodes_dict.keys()
#     I_nodes_label_list = qI_nodes_dict.keys()

#     B_nodes_index_list = [node_label2index[k] for k in B_nodes_label_list]
#     I_nodes_index_list = [node_label2index[k] for k in I_nodes_label_list]

#     if len(B_nodes_index_list) == 0:
#         raise ValueError('Boundary nodes must be a non empty set.')
    
#     if len(np.intersect1d(B_nodes_index_list, I_nodes_index_list)) > 0:
#         raise ValueError('Boundary nodes and internal nodes must be disjoint sets.')

#     if len(B_nodes_index_list) + len(I_nodes_index_list) != G.number_of_nodes():
#         raise ValueError('Boundary nodes and internal nodes must cover all nodes in the graph.')

#     if len(B_nodes_index_list) == G.number_of_nodes():
#         # All nodes are boundary nodes, no inner nodes
#         # "Sort" the dictionary to have the same order as the nodes (labels) in the graph

#         # Potential at all nodes
#         p_nodes_dict = {k: pB_nodes_dict[k] for k in G.nodes()}        

#         if return_entering_fluxes:
#             p = np.array(list(pB_nodes_dict.values()))
        
#             # Compute entering flux at all nodes
#             L = nx.laplacian_matrix(G, weight=edge_weight) # sparse matrix in CSR format
#             q = L @ p
#             q_nodes_dict = {k: q[k] for k in G.nodes()}

#             return p_nodes_dict, q_nodes_dict

#         return p_nodes_dict

#     if nx.number_connected_components(G.subgraph(I_nodes_label_list)) > 1:
#         raise ValueError('Removing boundary nodes disconnect the graph.')

#     # Solve the flow problem: get pI (potential at inner nodes)
#     pB = np.array(list(pB_nodes_dict.values()))
#     qI = np.array(list(qI_nodes_dict.values()))

#     L = nx.laplacian_matrix(G, weight=edge_weight) # sparse matrix in CSR format

#     L_I = L[I_nodes_index_list, :]
#     L_II = L_I[:, I_nodes_index_list]
#     L_IB = L_I[:, B_nodes_index_list]

#     pI = scipy.sparse.linalg.spsolve(L_II, qI - L_IB @ pB)

#     # Set dictionary of potential at inner nodes
#     pI_nodes_dict = {k: v for k, v in zip(I_nodes_label_list, pI)}

#     # Concatenate dictionaries pI_nodes_dict and pB_nodes_dict
#     p_nodes_dict = pI_nodes_dict.copy()
#     p_nodes_dict.update(pB_nodes_dict)  
#     # dict(**pB_nodes_dict, **pI_nodes_dict) # concatenate dictionaries (works if keys are strings)

#     # "Sort" the dictionary to have the same order as the nodes (labels) in the graph
#     p_nodes_dict = {k: p_nodes_dict[k] for k in G.nodes()}

#     if return_entering_fluxes:
#         # Solve the flow problem: get qB (entering flux at boundary nodes)
#         L_BB = L[B_nodes_index_list, :][:, B_nodes_index_list]

#         qB = L_IB.transpose() @ pI + L_BB @ pB

#         # Set dictionary of entering flux at boundary nodes
#         qB_nodes_dict = {k: v for k, v in zip(B_nodes_label_list, qB)}


#         # Concatenate dictionaries qI_nodes_dict and qB_nodes_dict
#         q_nodes_dict = qI_nodes_dict.copy()
#         q_nodes_dict.update(qB_nodes_dict)  

#         # "Sort" the dictionary to have the same order as the nodes (labels) in the graph
#         q_nodes_dict = {k: q_nodes_dict[k] for k in G.nodes()}

#         return p_nodes_dict, q_nodes_dict

#     return p_nodes_dict
# # ----------------------------------------------------------------------------

