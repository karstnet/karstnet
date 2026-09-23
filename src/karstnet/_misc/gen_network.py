"""
Module `_misc.gen_network`
--------------------------

The module `_misc.gen_network` is accessible as `misc`. It contains functions
to generate miscellaneous network.
"""

import numpy as np
import networkx as nx
import scipy

import karstnet as kn

# ===== Functions for generating grid graph ===================================
# Julien Straubhaar

# -----------------------------------------------------------------------------
def graph_grid(
        dimension, 
        spacing=None, 
        origin=None, 
        pos_attr='pos', 
        integer_labels=True):
    """
    Generates a grid in 1D, 2D or 3D as graph.

    Node position are sequences of 2 (resp. 3) floats for 1D and 2D (resp. 3D).
    
    Parameters
    ----------
    dimension : [sequence of] int(s)
        number of nodes along each axis

    spacing : [sequence of] float(s), optional
        distance between adjacent grid nodes along each axis

        by default (`None`): 1.0 along each axis

    origin : [sequence of] float(s), optional
        origin of the grid (position of the first node)

        by default (`None`): 0.0 along each axis

    pos_attr : str, default: 'pos'
        name of the node attribute for position

    integer_labels : bool, default: True
        - if `True`: the nodes labels are integers starting from 0
        - if `False`: the nodes labels are 2-tuples
  
    Returns
    -------
    g : networkx.Graph
        graph of the grid
    """
    dimension = np.atleast_1d(dimension).reshape(-1).astype('int')
    
    # Set space dimension
    d = dimension.shape[0]

    if spacing is None:
        spacing = np.ones(d)
    else:
        spacing = np.atleast_1d(spacing).reshape(-1).astype('float')
    
    if origin is None:
        origin = np.zeros(d)
    else:
        origin = np.atleast_1d(origin).reshape(-1).astype('float')

    if d == 1:
        g = nx.grid_2d_graph(dimension[0], 1, periodic=False)
        pos = {(i, j): np.array([origin[0] + i*spacing[0], 0.0]) for i, j in g.nodes()} # dictionary of node positions
    
    elif d == 2:
        g = nx.grid_2d_graph(dimension[0], dimension[1], periodic=False)
        pos = {(i, j): np.array([origin[0] + i*spacing[0], origin[1] + j*spacing[1]]) for i, j in g.nodes()} # dictionary of node positions

    elif d == 3:
        g = nx.grid_graph(dim=tuple(dimension[::-1]), periodic=False)
        pos = {(i, j, k): np.array([origin[0] + i*spacing[0], origin[1] + j*spacing[1], origin[2] + k*spacing[2]]) for i, j, k in g.nodes()} # dictionary of node positions
    
    else:
        raise ValueError('Not a 1D, 2D or 3D grid!')
    
    nx.set_node_attributes(g, pos, pos_attr)
    
    if integer_labels:
        g = nx.convert_node_labels_to_integers(g)
    
    return g
# -----------------------------------------------------------------------------

# ===== Functions for generating "honeycomb" graph (2D) =======================
# Julien Straubhaar

# -----------------------------------------------------------------------------
def graph_honeycomb(
        nrow, 
        ncol,
        edge_length=1.0,
        origin=None,
        pos_attr='pos', 
        integer_labels=True):
    """
    Generates a "honeycomb" grid in 2D (hexagonal mesh).

    This function uses the function `hexagonal_lattice_graph` from `networkx`.

    Parameters
    ----------
    nrow : int
        number of rows of hexagons

    ncol : int
        number of columns of hexagons

    edge_length : float, default: 1.0
        length of the edge of the hexagon

    origin : [sequence of] float(s), optional
        "origin", position of the bottom left node of the 
        bottom left hexagon is `origin + (edge_length, 0)`

        by default (`None`): origin = (0.0, 0.0)

    pos_attr : str, default: 'pos'
        name of the node attribute for position

    integer_labels : bool, default: True
        - if `True`: the nodes labels are integers starting from 0
        - if `False`: the nodes labels are 2-tuples

    Returns
    -------
    g : networkx.Graph
        graph of the "honeycomb" grid
    """   

    g = nx.hexagonal_lattice_graph(nrow, ncol, with_positions=True, periodic=False)

    # Get position
    pos = nx.get_node_attributes(g, 'pos')

    if origin is None:
        origin = np.zeros(2)
    else:
        origin = np.asarray(origin)

    # Set position (according to edge_length (scaling) and origin (shift)) as numpy array
    pos = {u: origin + edge_length * np.asarray(p) for u, p in pos.items()}

    if pos_attr != 'pos':
        # remove 'pos' attribute
        for u in g.nodes():
            del(g.nodes[u]['pos'])

    # Set position attribute
    nx.set_node_attributes(g, pos, pos_attr)
    
    if integer_labels:
        g = nx.convert_node_labels_to_integers(g)
    
    return g
# -----------------------------------------------------------------------------

# ===== Functions for generating graph on torus (periodic grid) ===============
# Julien Straubhaar

# -----------------------------------------------------------------------------
def graph_torus(
        rh, rs, nh, ns, 
        pos_attr='pos', 
        integer_labels=True):
    """
    Generates a 2D periodic grid with node positions on a torus.

    Parameters
    ----------
    rh : float
        positive number, radius of the "hole circle" of the torus 
        (distance from the origin to the center of "section circle")
    
    rs : float
        positive number, radius of the "circle section" of the torus

    nh : int
        number of nodes along "hole circle"
    
    ns : int
        number of nodes along "section circle"

    pos_attr : str, default: 'pos'
        name of the node attribute for position

    integer_labels : bool, default: True
        - if `True`: the nodes labels are integers starting from 0
        - if `False`: the nodes labels are 2-tu
  
    Returns
    -------
    g : networkx.Graph
        graph of the grid
    """

    g = nx.grid_2d_graph(ns, nh, periodic=True)

    angle_h = np.linspace(0, 2*np.pi, nh+1) # angles along the "hole circle"
    angle_s = np.linspace(0, 2*np.pi, ns+1) # angles along the "section circle"
    pos = {(i, j): np.array([
                        (rh + rs*np.cos(angle_s[i]))*np.cos(angle_h[j]), 
                        (rh + rs*np.cos(angle_s[i]))*np.sin(angle_h[j]),
                        rs*np.sin(angle_s[i])
        ]) for i, j in g.nodes()} # dictionary of node positions
    
    nx.set_node_attributes(g, pos, pos_attr)

    if integer_labels:
        g = nx.convert_node_labels_to_integers(g)

    return g
# -----------------------------------------------------------------------------

# ===== Functions for generating graph from Delaunay triangulation ============
# Julien Straubhaar

# -----------------------------------------------------------------------------
def graph_from_delaunay_tri(
        points, 
        pos_attr='pos', 
        min_edge_len=None, 
        max_edge_len=None):
    """
    Generates a graph from a Delaunay triangulation of points.

    Edges whose length is less than `min_edge_len` (if specified)
    or greater than `max_edge_len` (if specified), are removed.

    The node labels are integers from 0 and corresponds to the indices
    of the points.

    Parameters
    ----------
    points : 2d array
        points, each row is the coordinates of one point

    pos_attr : str, default: 'pos'
        name of the node attribute for position

    min_edge_len : float, optional
        minimum edge length

    max_edge_len : float, optional
        maximum edge length
  
    Returns
    -------
    g : networkx.Graph
        graph of the grid
    """
    # Get Delaunay triangulation (scipy)
    tri = scipy.spatial.Delaunay(points)

    # Set networkx graph of the triangulation
    g = nx.Graph()
    for simp in tri.simplices:
        nx.add_path(g, np.hstack((simp, simp[0:1])).tolist())
    
    nx.set_node_attributes(g, {i: p for i, p in enumerate(points)}, pos_attr)

    if min_edge_len is not None or max_edge_len is not None:
        # Compute edge square length (dictionary)
        edge_len_2 = {e:np.sum((np.asarray(g.nodes[e[0]][pos_attr]) - np.asarray(g.nodes[e[1]][pos_attr]))**2) for e in g.edges()}
        
        # Get list of edges to be removed
        if min_edge_len is None:
            min_edge_len_2 = 0.0
        else:
            min_edge_len_2 = min_edge_len**2

        if max_edge_len is None:
            max_edge_len_2 = np.inf
        else:
            max_edge_len_2 = max_edge_len**2

        edges_to_remove = [e for e in g.edges() if edge_len_2[e] < min_edge_len_2 or edge_len_2[e] > max_edge_len_2]

        # Remove edges
        g.remove_edges_from(edges_to_remove)

    return g
# -----------------------------------------------------------------------------

# ===== Functions for generating classical fractals ===========================
# Julien Straubhaar

# -----------------------------------------------------------------------------
def von_koch(
        n, 
        first_edge=[(0.0, 0.0), (1.0, 0.0)], 
        pos_attr='pos', 
        close=False):
    """
    Generates von Koch snow flakes fractal as graph.
    
    Parameters
    ----------
    n : int
        number of iterations
    
    first_edge : sequence of 2 sequences of 2 floats, default: [(0.0, 0.0), (1.0, 0.0)]
        position of the first two nodes, at iteration 0

    pos_attr : str, default: 'pos'
        name of the node attribute for position

    close : bool, default: False
        - if True : a "closed" curve (loop) is generated, \
        starting  with a triangle at iteration 0
        - if True : a curve is generated, \
        starting with a segment at iteration 0
    
    Returns
    -------
    g : networkx.Graph
        graph of the "curve"
    """
    sq3 = np.sqrt(3)

    g = nx.Graph()
    g.add_edge(0, 1)

    p0 = np.asarray(first_edge[0])
    p1 = np.asarray(first_edge[1])
    nx.set_node_attributes(g, {0:p0, 1:p1}, pos_attr)

    if close:
        p = p1 - p0
        pn = np.array([-p[1], p[0]])
        p2 = p0 + 0.5 * (p + sq3*pn)
        g.add_node(2, pos=p2)
        g.add_edges_from([(1, 2), (2, 0)])

    nn = g.number_of_nodes()
    for _ in range(n):
        ne = g.number_of_edges()
        u = 0
        for _ in range(ne):
            v = list(g.edges(u))[0][-1]          
            g.remove_edge(u, v)
            pu = g.nodes[u][pos_attr]
            pv = g.nodes[v][pos_attr]            
            p1 = (2*pu + pv)/3
            p3 = (pu + 2*pv)/3
            p = p3 - p1
            pn = np.array([-p[1], p[0]])
            if close:
                pn = - pn
            p2 = p1 + 0.5 * (p + sq3*pn)
            g.add_node(nn, pos=p1)
            g.add_node(nn+1, pos=p2)
            g.add_node(nn+2, pos=p3)
            g.add_edges_from([(u, nn), (nn, nn+1), (nn+1, nn+2), (nn+2, v)])
            nn = nn+3
            u = v

    return g
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
def sierpinski(
        n, 
        first_edge=[(0.0, 0.0), (1.0, 0.0)], 
        pos_attr='pos'):
    """
    Generates Sierpinski triangle fractal as graph.
    
    Parameters
    ----------
    n : int
        number of iterations

    first_edge : sequence of 2 sequences of 2 floats, default: [(0.0, 0.0), (1.0, 0.0)]
        position of the first two nodes, at iteration 0

    pos_attr : str, default: 'pos'
        name of the node attribute for position

    Returns
    -------
    g : networkx.Graph
        graph of the "curve"
    """
    sq3 = np.sqrt(3)

    p0 = np.asarray(first_edge[0])
    p1 = np.asarray(first_edge[1])
    p = p1 - p0
    pn = np.array([-p[1], p[0]])
    p2 = p0 + 0.5 * (p + sq3*pn)

    points = np.vstack((p0, p1, p2))
    simplices = np.array([[0, 1, 2]])
    
    for _ in range(n):
        nn = len(points)
        new_points = np.empty((0, 2), dtype='float')
        new_simplices = np.empty((0, 3), dtype='int')
        for s in simplices:
            p0, p1, p2 = points[s]
            new_points = np.vstack((new_points, 0.5*np.array([p0+p1, p1+p2, p2+p0])))
            new_simplices = np.vstack((new_simplices, np.array([[s[0], nn, nn+2], [s[1], nn+1, nn], [s[2], nn+2, nn+1]])))
            nn = nn+3

        points = np.vstack((points, new_points))
        simplices = new_simplices

    g = nx.Graph()
    for simp in simplices:
        nx.add_path(g, np.hstack((simp, simp[0:1])).tolist())

    nx.set_node_attributes(g, {i:p for i, p in enumerate(points)}, pos_attr)

    return g
# ----------------------------------------------------------------------------

# ============================================================================
# Functions for generating random clusters in regular grid (percolation network)
# ============================================================================
# Julien Straubhaar

# ----------------------------------------------------------------------------
def graph_from_binary_image_grid_2d(
        arr,
        sx=1.0, 
        sy=1.0,
        pos_attr='pos'):
    """
    Returns the graph corresponding to a binary image in a 2D grid.

    The graph of activated cells, with edges linking adjacent cells, is computed.

    Parameters
    ----------
    arr : 2d array
        binary image, array of shape (my, mx), where mx (resp. my) is the number
        of cells along x-axis (resp. y-axis); cell with non zero values are
        considered as activated

    sx : float, default: 1.0
        cell size along x-axis

    sy : float, default: 1.0
        cell size along y-axis

    pos_attr : str, default: 'pos'
        name of the node attribute for position

    Returns
    -------
    g : networkx.Graph
        graph corresponding to the input array `arr`, with position
        (sequence of 2 floats) as node attributes
    """
    # Position
    idy, idx = np.where(arr)
    pos = np.vstack((idx*sx, idy*sy)).T

    # Array of node id in the array (-1 for no node)
    arr_id = -1 * np.ones_like(arr, dtype='int')
    arr_id[arr!=0] = np.arange((arr!=0).sum())

    # Initialize graph
    g = nx.Graph()

    # Add nodes
    g.add_nodes_from(np.arange(len(pos)))
    
    # Set node position
    nx.set_node_attributes(g, {i:p for i, p in enumerate(pos)}, pos_attr)

    # Add edges along x-axis
    idy, idx = np.where(arr[:, :-1] * arr[:, 1:])
    g.add_edges_from(zip(arr_id[idy, idx], arr_id[idy, idx+1]))

    # Add edges along y-axis
    idy, idx = np.where(arr[:-1, :] * arr[1:, :])
    g.add_edges_from(zip(arr_id[idy, idx], arr_id[idy+1, idx]))
        
    return g
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
def graph_from_binary_image_grid_3d(
        arr,
        sx=1.0, 
        sy=1.0,
        sz=1.0,
        pos_attr='pos'):
    """
    Returns the graph corresponding to a binary image in a 3D grid.

    The graph of activated cells, with edges linking adjacent cells, is computed.

    Parameters
    ----------
    arr : 3d array
        binary image, array of shape (mz, my, mx), where mx (resp. my, mz) is the 
        number of cells along x-axis (resp. y-axis, z-axis); cell with non zero 
        values are considered as activated

    sx : float, default: 1.0
        cell size along x-axis

    sy : float, default: 1.0
        cell size along y-axis

    sz : float, default: 1.0
        cell size along z-axis

    pos_attr : str, default: 'pos'
        name of the node attribute for position

    Returns
    -------
    g : networkx.Graph
        graph corresponding to the input array `arr`, with position
        (sequence of 3 floats) as node attributes 
    """
    # Position
    idz, idy, idx = np.where(arr)
    pos = np.vstack((idx*sx, idy*sy, idz*sz)).T

    # Array of node id in the array (-1 for no node)
    arr_id = -1 * np.ones_like(arr, dtype='int')
    arr_id[arr!=0] = np.arange((arr!=0).sum())

    # Initialize graph
    g = nx.Graph()

    # Add nodes
    g.add_nodes_from(np.arange(len(pos)))
    
    # Set node position
    nx.set_node_attributes(g, {i:p for i, p in enumerate(pos)}, pos_attr)

    # Add edges along x-axis
    idz, idy, idx = np.where(arr[:, :, :-1] * arr[:, :, 1:])
    g.add_edges_from(zip(arr_id[idz, idy, idx], arr_id[idz, idy, idx+1]))

    # Add edges along y-axis
    idz, idy, idx = np.where(arr[:, :-1, :] * arr[:, 1:, :])
    g.add_edges_from(zip(arr_id[idz, idy, idx], arr_id[idz, idy+1, idx]))
    
    # Add edges along z-axis
    idz, idy, idx = np.where(arr[:-1, :, :] * arr[1:, :, :])
    g.add_edges_from(zip(arr_id[idz, idy, idx], arr_id[idz+1, idy, idx]))
        
    return g
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
def random_cluster_in_grid_2d(
        mx, my, probability, 
        sx=1.0, 
        sy=1.0,
        pos_attr='pos', 
        extract_connected_components=True,
        renumbering_nodes=True,
        seed=None, 
        return_array=False):
    """
    Generates graph(s) from a binary medium in a 2D grid.

    Each cell in a 2D grid (array) is activated with a given
    probability `probability` (Bernoulli distribution). 
    Then the graph of activated cells, with edges linking 
    adjacent cells, is computed.

    Parameters
    ----------
    mx : int
        number of cells along x-axis

    my : int
        number of cells along y-axis
    
    probability : float
        float in (0, 1), probability of cell (node) activation

    sx : float, default: 1.0
        cell size along x-axis

    sy : float, default: 1.0
        cell size along y-axis

    pos_attr : str, default: 'pos'
        name of the node attribute for position

    extract_connected_components : bool, default: True
        - if `True`: all the connected components of the graph \
        are extracted as subgraph, and returned in a list, \
        with the number of nodes in decreasing order
        - if `False`: The whole graph is returned
    
    renumbering_nodes : bool, default: True
        used if `extract_connected_components=True`:
        
        - if `True`: nodes of each connected components are \
        renumbered from 0
        - if `False`: original node numbers are kept in each \
        connected components

    seed : int, optional
        seed number for initializing the random number generator
    
    return_array : bool, default: False
        indicates if the binary array is returned

    Returns
    -------
    out : networkx.Graph, or list of network.Graph
        resulting graph(s):

        - if `extract_connected_components=True`, `out` \
        is a list of graphs (one graph for each connected \
        component)
        - if `extract_connected_component=False`, `out` \
        is a graph

    arr : 2d array, optional
        returned if `return_arrray=True`, 2d array of 
        shape `(my, mx)` of 0 and 1 for inactive and 
        active cells respectively
    """
    if seed is not None:
        np.random.seed(seed)

    # Generate binary array
    arr = 1 * (np.random.random(size=(my, mx)) < probability)

    # Array of node id in the array (-1 for no node)
    arr_id = arr - 1
    arr_id[arr==1] = np.arange((arr==1).sum())

    # Position
    # my, mx = arr.shape
    mmy, mmx = np.meshgrid(np.arange(my), np.arange(mx), indexing='ij')
    id_xy = np.vstack((mmx[arr==1].reshape(-1), mmy[arr==1].reshape(-1))).T
    pos = id_xy*np.array([sx, sy])

    # Initialize graph
    g = nx.Graph()

    # Add nodes
    g.add_nodes_from(np.arange(len(pos)).tolist())
    
    # Set node position
    nx.set_node_attributes(g, {i:p for i, p in enumerate(pos)}, pos_attr)

    # Add edges along x-axis
    idy, idx = np.where(arr[:, :-1] * arr[:, 1:])
    g.add_edges_from(zip(arr_id[idy, idx].tolist(), arr_id[idy, idx+1].tolist()))

    # Add edges along y-axis
    idy, idx = np.where(arr[:-1, :] * arr[1:, :])
    g.add_edges_from(zip(arr_id[idy, idx].tolist(), arr_id[idy+1, idx].tolist()))
    
    if extract_connected_components:
        g_list = kn.utils.get_list_of_connected_components(g, order='descending', renumbering_nodes=renumbering_nodes)
        out = g_list
    else:
        out = g

    if return_array:
        return out, arr
    
    return out
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
def random_cluster_in_grid_3d(
        mx, my, mz, probability, 
        sx=1.0, 
        sy=1.0, 
        sz=1.0,
        pos_attr='pos',
        extract_connected_components=True,
        renumbering_nodes=True,
        seed=None, 
        return_array=False):
    """
    Generates graph(s) from a binary medium in a 3D grid.

    Each cell in a 3D grid (array) is activated with a given
    probability `probability` (Bernoulli distribution). 
    Then the graph of activated cells, with edges linking 
    adjacent cells, is computed.

    Parameters
    ----------
    mx : int
        number of cells along x-axis

    my : int
        number of cells along y-axis
    
    mz : int
        number of cells along z-axis
    
    probability : float
        float in (0, 1), probability of activation
        (of each cell)

    sx : float, default: 1.0
        cell size along x-axis

    sy : float, default: 1.0
        cell size along y-axis

    sz : float, default: 1.0
        cell size along z-axis

    pos_attr : str, default: 'pos'
        name of the node attribute for position

    extract_connected_components : bool, default: True
        - if `True`: all the connected components of the graph \
        are extracted as subgraph, and returned in a list, \
        with the number of nodes in decreasing order
        - if `False`: The whole graph is returned
    
    renumbering_nodes : bool, default: True
        used if `extract_connected_components=True`:
        
        - if `True`: nodes of each connected components are \
        renumbered from 0
        - if `False`: original node numbers are kept in each \
        connected components

    seed : int, optional
        seed number for initializing the random number generator
    
    return_array : bool, default: False
        indicates if the binary array is returned

    Returns
    -------
    out : networkx.Graph, or list of network.Graph
        resulting graph(s):

        - if `extract_connected_components=True`, `out` \
        is a list of graphs (one graph for each connected \
        component)
        - if `extract_connected_component=False`, `out` \
        is a graph

    arr : 3d array, optional
        returned if `return_arrray=True`, 3d array of 
        shape `(mz, my, mx)` of 0 and 1 for inactive and 
        active cells respectively
    """
    if seed is not None:
        np.random.seed(seed)

    # Generate binary array
    arr = 1 * (np.random.random(size=(mz, my, mx)) < probability)

    # Array of node id in the array (-1 for no node)
    arr_id = arr - 1
    arr_id[arr==1] = np.arange((arr==1).sum())

    # Position
    # mz, my, mx = arr.shape
    mmz, mmy, mmx = np.meshgrid(np.arange(mz), np.arange(my), np.arange(mx), indexing='ij')
    id_xyz = np.vstack((mmx[arr==1].reshape(-1), mmy[arr==1].reshape(-1), mmz[arr==1].reshape(-1))).T
    pos = id_xyz*np.array([sx, sy, sz])

    # Initialize graph
    g = nx.Graph()

    # Add nodes
    g.add_nodes_from(np.arange(len(pos)))
    
    # Set node position
    nx.set_node_attributes(g, {i:p for i, p in enumerate(pos)}, pos_attr)

    # Add edges along x-axis
    idz, idy, idx = np.where(arr[:, :, :-1] * arr[:, :, 1:])
    g.add_edges_from(zip(arr_id[idz, idy, idx], arr_id[idz, idy, idx+1]))

    # Add edges along y-axis
    idz, idy, idx = np.where(arr[:, :-1, :] * arr[:, 1:, :])
    g.add_edges_from(zip(arr_id[idz, idy, idx], arr_id[idz, idy+1, idx]))
    
    # Add edges along z-axis
    idz, idy, idx = np.where(arr[:-1, :, :] * arr[1:, :, :])
    g.add_edges_from(zip(arr_id[idz, idy, idx], arr_id[idz+1, idy, idx]))

    if extract_connected_components:
        g_list = kn.utils.get_list_of_connected_components(g, order='descending', renumbering_nodes=renumbering_nodes)
        out = g_list
    else:
        out = g

    if return_array:
        return out, arr
    
    return out
# ----------------------------------------------------------------------------

# ===== Functions for Diffusion-Limited Aggregation (DLA) - Brownian tree ====
# Julien Straubhaar

# ----------------------------------------------------------------------------
def dla_in_grid_2d(
        arr, 
        nparticles,
        nparticles_try_max=None, 
        n_random_walk_step_max=1000000,
        stop_at_percolation=False, 
        link_to_all_neighbors=False,
        sx=1.0, 
        sy=1.0,
        pos_attr='pos',
        dla_step_attr='dla_step',
        seed=None,
        inplace=False,
        return_array=False,
        verbose=1):
    """
    DLA (Diffusion-Limited Aggregation).

    A simulation consists in generating random walk of particles in a regular 
    grid from a starting area to an ending (attracting) area, retaining the 
    last position (beside the ending area) and adding this position to the 
    ending area for the next particle.

    In this way, a DLA simulation array is generated, with values:
        
    - 1 : in starting area
    - n > 1: at ending point of the (n-1)-th particle
    - 0 : elsewhere

    Then, a graph is computed from the DLA simulation array; each cell in 
    the array with a value greater than 0 is a node of the graph, and each 
    node is linked (with edge(s)) to: one neighbor with a value (from the array) 
    that is smaller than (or equal to) its own value (if 
    `link_to_all_neighbor=False`), or to all all its neighbors (adjacent cell in 
    the array) (if `link_to_all_neighbor=True`).
       
    Parameters
    ----------
    arr: 2d array
        starting array, with value:
        
        - -1 : in starting area of the particles
        - 1 : in ending area
        - 0 : elsewhere

        Note: point in the border of the array will be removed from 
        the starting area and ending area
    
    nparticles : int
        number of particles for DLA simulation
    
    nparticles_try_max : int, optional
        maximal number of trials (particles);
        by default: set to 100*nparticles
    
    n_random_walk_step_max : int, default: 1000000
        maximal number of steps in the random walk of each particle

    stop_at_percolation : bool, default: False
        if `True`: the simulation is stopped, i.e. no new particle is considered,
        when percolation is reached, i.e. when one cell in the starting area is 
        in the ending area
        
    link_to_all_neighbors : bool, default : False
        - if `True`: each graph node is linked with edges to all its neighbors
        - if `False`: each graph node is linked with an edge to one of its \
        neighbor (selected randomly)

    sx : float, default: 1.0
        cell size along x-axis

    sy : float, default: 1.0
        cell size along y-axis

    pos_attr : str, default: 'pos'
        name of the node attribute for position

    dla_step_attr : str, default: 'dla_step'
        name of the node attribute for DLA simulation step (particule number)

    seed : int, optional
        seed number for initializing the random number generator

    inplace : bool, default: False
        - if `True`: input array `arr` is modified (inplace operation)
        - if `False`: input array `mat` is not modified

    return_array : bool, default: False
        indicates if the DLA simulation arrray is returned

    verbose : int, default: 1
        verbose mode, larger value implies more info printed
    
    Returns
    -------
    g : networkx.Graph
        resulting graph, with the following node attributes:
        
        - 'pos' : sequence of 2 floats, position of the node
        - 'dla_step' : int, corresponding value in the DLA simulation \
        array (see `dla_arr` below)

    dla_arr: 2d array, optional
        returned if `return_array=True`, resulting DLA simulation array, 
        with values:
        
        - 1 : in starting area
        - n > 1: at ending point of the (n-1)-th particle
        - 0 : elsewhere

        note: if `inplace=True`, this output array is the 
        input array (modified)
    """
    # Compute DLA simulation array
    # ----------------------------
    if verbose > 0:
        print(f'Compute DLA simulation array')

    if inplace:
        dla_arr = arr
    else:
        dla_arr = np.copy(arr)

    # Set border
    dla_arr[np.array([0, -1]), :] = -2
    dla_arr[:, np.array([0, -1])] = -2
    
    ind_start_list = np.vstack(np.where(dla_arr==-1)).T
    nstart = len(ind_start_list)

    if nstart <= 0:
        raise ValueError('No point in starting area (in `arr`)')

    if not np.any(dla_arr==1):
        raise ValueError('No point in ending area (in `arr`)')

    if seed is not None:
        np.random.seed(seed)

    if nparticles_try_max is None:
        nparticles_try_max = 100*nparticles

    nn = 1 # current particle
    for k in range(nparticles_try_max):
        if verbose > 0:
            print(f'Current particle {nn} over {nparticles} ; attempt {k+1} over {nparticles_try_max}')

        # Select a starting point
        ind = ind_start_list[np.random.randint(nstart)].copy()

        if dla_arr[tuple(ind)] > 0:
            # node already selected in the ending area (percolation)
            if stop_at_percolation:
                if verbose > 0:
                    print('  Particle at start already in ending area (percolation reached): stop DLA simulation')
                break
            else:
                if verbose > 0:
                    print('  Particle at start already in ending area (percolation reached): start with a new particle')
                continue
        
        for j in range(n_random_walk_step_max):
            if dla_arr[tuple(ind)] == -2:
                # node on the border
                if verbose > 1:
                    print(f'   Particle on the border (step = {j}): start with a new particle')
                break
            
            if         dla_arr[ind[0]-1, ind[1]] > 0 \
                    or dla_arr[ind[0]+1, ind[1]] > 0 \
                    or dla_arr[ind[0], ind[1]-1] > 0 \
                    or dla_arr[ind[0], ind[1]+1] > 0:
                # beside ending area
                if verbose > 2:
                    print(f'   Particle aggregated (step = {j})')
                nn = nn+1
                dla_arr[tuple(ind)] = nn
                break
            
            ii = np.random.randint(4)
            if ii == 0:
                ind[0] = ind[0]-1
            elif ii == 1:
                ind[0] = ind[0]+1
            elif ii == 2:
                ind[1] = ind[1]-1
            else: # ii == 3:
                ind[1] = ind[1]+1

        if nn > nparticles:
            break

    # Reset starting and border area
    dla_arr[dla_arr<0] = 0
    
    # Set graph
    # ---------
    if verbose > 0:
        print(f'Set graph')

    # Index of DLA simulation array corresponding to the graph nodes
    iy, ix = np.where(dla_arr>0)

    # Number of graph nodes
    nn = len(ix)

    # Index (x, y) of matrix cells
    my, mx = dla_arr.shape
    yy, xx = np.meshgrid(np.arange(my), np.arange(mx), indexing='ij')

    # Attributes of nodes
    pos = np.vstack((xx[iy, ix].reshape(-1)*sx, yy[iy, ix].reshape(-1)*sy)).T # position (accounting for cell size)
    step = dla_arr[iy, ix]

    # Initialize graph
    g = nx.Graph()
    
    if link_to_all_neighbors:
        for i in range(nn-1):
            neigh = np.where(np.abs(iy[i] - iy[i+1:]) + np.abs(ix[i] - ix[i+1:]) == 1)[0].tolist()
            g.add_edges_from([(i, i+1+ne) for ne in neigh])
    else:
        for i in range(nn):
            neigh = np.where(np.abs(iy[i] - iy[:]) + np.abs(ix[i] - ix[:]) == 1)[0].tolist()
            k = np.where(step[i] >= step[neigh])[0]
            if len(k) == 1:
                ne = neigh[k[0]]
                g.add_edge(i, ne)
            elif len(k) > 1:
                ne = neigh[k[np.random.randint(len(k))]]
                g.add_edge(i, ne)

    # Set node attributes
    nx.set_node_attributes(g, {i:p for i, p in enumerate(pos)}, pos_attr)
    nx.set_node_attributes(g, {i:s for i, s in enumerate(step)}, dla_step_attr)

    # Return
    # ------
    if return_array:
        return g, dla_arr
    
    return g
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
def dla_in_grid_3d(
        arr, 
        nparticles,
        nparticles_try_max=None, 
        n_random_walk_step_max=1000000, 
        stop_at_percolation=False,
        link_to_all_neighbors=False,
        sx=1.0, 
        sy=1.0,
        sz=1.0,
        pos_attr='pos',
        dla_step_attr='dla_step',
        seed=None,
        inplace=False,
        return_array=False,
        verbose=1):
    """
    DLA (Diffusion-Limited Aggregation).

    A simulation consists in generating random walk of particles in a regular 
    grid from a starting area to an ending (attracting) area, retaining the 
    last position (beside the ending area) and adding this position to the 
    ending area for the next particle.

    In this way, a DLA simulation array is generated, with values:
        
    - 1 : in starting area
    - n > 1: at ending point of the (n-1)-th particle
    - 0 : elsewhere

    Then, a graph is computed from the DLA simulation array; each cell in 
    the array with a value greater than 0 is a node of the graph, and each 
    node is linked (with edge(s)) to: one neighbor with a value (from the array) 
    that is smaller than (or equal to) its own value (if 
    `link_to_all_neighbor=False`), or to all all its neighbors (adjacent cell in 
    the array) (if `link_to_all_neighbor=True`).
       
    Parameters
    ----------
    arr: 3d array
        starting array, with value:
        
        - -1 : in starting area of the particles
        - 1 : in ending area
        - 0 : elsewhere

        Note: point in the border of the array will be removed from 
        the starting area and ending area
    
    nparticles : int
        number of particles for DLA simulation
    
    nparticles_try_max : int, optional
        maximal number of trials (particles);
        by default: set to 100*nparticles
    
    n_random_walk_step_max : int, default: 1000000
        maximal number of steps in the random walk of each particle

    stop_at_percolation : bool, default: False
        if `True`: the simulation is stopped, i.e. no new particle is considered,
        when percolation is reached, i.e. when one cell in the starting area is 
        in the ending area

    link_to_all_neighbors : bool, default : False
        - if `True`: each graph node is linked with edges to all its neighbors
        - if `False`: each graph node is linked with an edge to one of its \
        neighbor (selected randomly)

    sx : float, default: 1.0
        cell size along x-axis

    sy : float, default: 1.0
        cell size along y-axis

    sz : float, default: 1.0
        cell size along z-axis

    pos_attr : str, default: 'pos'
        name of the node attribute for position

    dla_step_attr : str, default: 'dla_step'
        name of the node attribute for DLA simulation step (particule number)

    seed : int, optional
        seed number for initializing the random number generator

    inplace : bool, default: False
        - if `True`: input array `arr` is modified (inplace operation)
        - if `False`: input array `mat` is not modified

    return_array : bool, default: False
        indicates if the DLA simulation arrray is returned

    verbose : int, default: 1
        verbose mode, larger value implies more info printed
    
    Returns
    -------
    g : networkx.Graph
        resulting graph, with the following node attributes:
        
        - 'pos' : sequence of 2 floats, position of the node
        - 'dla_step' : int, corresponding value in the DLA simulation \
        array (see `dla_arr` below)

    dla_arr: 3d array, optional
        returned if `return_array=True`, resulting DLA simulation array, 
        with values:
        
        - 1 : in starting area
        - n > 1: at ending point of the (n-1)-th particle
        - 0 : elsewhere

        note: if `inplace=True`, this output array is the 
        input array (modified)
    """
    # Compute DLA simulation array
    # ----------------------------
    if verbose > 0:
        print(f'Compute DLA simulation array')

    if inplace:
        dla_arr = arr
    else:
        dla_arr = np.copy(arr)

    # Set border
    dla_arr[np.array([0, -1]), :, :] = -2
    dla_arr[:, np.array([0, -1]), :] = -2
    dla_arr[:, :, np.array([0, -1])] = -2
    
    ind_start_list = np.vstack(np.where(dla_arr==-1)).T
    nstart = len(ind_start_list)

    if nstart <= 0:
        raise ValueError('No point in starting area (in `arr`)')

    if not np.any(dla_arr==1):
        raise ValueError('No point in ending area (in `arr`)')

    if seed is not None:
        np.random.seed(seed)

    if nparticles_try_max is None:
        nparticles_try_max = 100*nparticles

    nn = 1 # current particle
    for k in range(nparticles_try_max):
        if verbose > 0:
            print(f'Current particle {nn} over {nparticles} ; attempt {k+1} over {nparticles_try_max}')

        # Select a starting point
        ind = ind_start_list[np.random.randint(nstart)].copy()

        if dla_arr[tuple(ind)] > 0:
            # node already selected in the ending area (percolation)
            if stop_at_percolation:
                if verbose > 0:
                    print('  Particle at start already in ending area (percolation reached): stop DLA simulation')
                break
            else:
                if verbose > 0:
                    print('  Particle at start already in ending area (percolation reached): start with a new particle')
                continue
        
        for j in range(n_random_walk_step_max):
            if dla_arr[tuple(ind)] == -2:
                # node on the border
                if verbose > 1:
                    print(f'   Particle on the border (random walk step = {j}): start with a new particle')
                break
            
            if         dla_arr[ind[0]-1, ind[1], ind[2]] > 0 \
                    or dla_arr[ind[0]+1, ind[1], ind[2]] > 0 \
                    or dla_arr[ind[0], ind[1]-1, ind[2]] > 0 \
                    or dla_arr[ind[0], ind[1]+1, ind[2]] > 0 \
                    or dla_arr[ind[0], ind[1], ind[2]-1] > 0 \
                    or dla_arr[ind[0], ind[1], ind[2]+1] > 0:
                # beside ending area
                if verbose > 2:
                    print(f'   Particle aggregated (random walk step = {j})')
                nn = nn+1
                dla_arr[tuple(ind)] = nn
                break
            
            ii = np.random.randint(6)
            if ii == 0:
                ind[0] = ind[0]-1
            elif ii == 1:
                ind[0] = ind[0]+1
            elif ii == 2:
                ind[1] = ind[1]-1
            elif ii == 3:
                ind[1] = ind[1]+1
            elif ii == 4:
                ind[2] = ind[2]-1
            else: # ii == 5:
                ind[2] = ind[2]+1

        if nn > nparticles:
            break

    # Reset starting and border area
    dla_arr[dla_arr<0] = 0
    
    # Set graph
    # ---------
    if verbose > 0:
        print(f'Set graph')

    # Index of DLA simulation array corresponding to the graph nodes
    iz, iy, ix = np.where(dla_arr>0)

    # Number of graph nodes
    nn = len(ix)

    # Index (x, y) of matrix cells
    mz, my, mx = dla_arr.shape
    zz, yy, xx = np.meshgrid(np.arange(mz), np.arange(my), np.arange(mx), indexing='ij')

    # Attributes of nodes
    pos = np.vstack((xx[iz, iy, ix].reshape(-1)*sx, yy[iz, iy, ix].reshape(-1)*sy, zz[iz, iy, ix].reshape(-1)*sz)).T # position (accounting for cell size)
    step = dla_arr[iz, iy, ix]

    # Initialize graph
    g = nx.Graph()
    
    if link_to_all_neighbors:
        for i in range(nn-1):
            neigh = np.where(np.abs(iz[i] - iz[i+1:]) + np.abs(iy[i] - iy[i+1:]) + np.abs(ix[i] - ix[i+1:]) == 1)[0].tolist()
            g.add_edges_from([(i, i+1+ne) for ne in neigh])
    else:
        for i in range(nn):
            neigh = np.where(np.abs(iz[i] - iz[:]) + np.abs(iy[i] - iy[:]) + np.abs(ix[i] - ix[:]) == 1)[0].tolist()
            k = np.where(step[i] >= step[neigh])[0]
            if len(k) == 1:
                ne = neigh[k[0]]
                g.add_edge(i, ne)
            elif len(k) > 1:
                ne = neigh[k[np.random.randint(len(k))]]
                g.add_edge(i, ne)

    # Set node attributes
    nx.set_node_attributes(g, {i:p for i, p in enumerate(pos)}, pos_attr)
    nx.set_node_attributes(g, {i:s for i, s in enumerate(step)}, dla_step_attr)

    # Return
    # ------
    if return_array:
        return g, dla_arr
    
    return g
# ----------------------------------------------------------------------------
