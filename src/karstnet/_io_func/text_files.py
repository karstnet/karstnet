"""
Module `_io_func.text_files`
----------------------------

The module `_io_func.text_files` is accessible as `io_func`. It contains functions 
to import / export networkx graphs from / to text files.
"""


import pandas as pd 
import numpy as np
import networkx as nx
import karstnet as kn


# ===== Import functions from csv file(s) for a networkx graph ===============
# Julien Straubhaar

# ----------------------------------------------------------------------------
def networkx_from_csv(
        basename,
        suffix_nodes='_nodes.csv', 
        suffix_edges='_edges.csv', 
        delimiter_nodes=';',  
        delimiter_edges=';',
        pos_attr='pos',
        import_node_scalar_properties=True,
        import_edge_scalar_properties=True,
        default_start_id=0,
        verbose=True):
    """
    Creates a networkx graph from two csv files (nodes, and edges (links)).

    The file for nodes has the columns: 
        - ['id',] 'x', 'y'[, 'z', '<node_prop0>', '<node_prop1>', ...]

    Default id (if column 'id' is not given) are integer starting from 0.
    The columns with a name other than 'id', 'x', 'y', 'z', if any, are 
    additional scalar properties (attributes) attached to the nodes.
    The column 'z' can be omitted for 2D case.

    The file for edges (links) has the columns:
        - '<from>', '<to>'[, '<edge_prop0>', '<edge_prop1>', ...]

    where the first two columns, '<from>', '<to>', are the node ids 
    of the two extremities of an edge (link), the next columns (if any) are
    additional scalar properties (attributes) attached to the edges.
    
    Parameters
    ----------
    basename : str
        the base name (prefix) used for the input files

    suffix_nodes : str, default: '_nodes.csv'
        the input file containing the nodes is 
        `basename``suffix_nodes`
    
    suffix_edges : str, default: '_edges.csv'
        the input file containing the edges (links) is 
        `basename``suffix_edges`
    
    delimiter_nodes : str, default: ';'
        delimiter used in the file for nodes
    
    delimiter_edges : str, default: ';'
        delimiter used in the file for edges
    
    pos_attr : str, default: 'pos'
        name of the node attribute for position

    import_node_scalar_properties : bool, default: True
        indicates if the scalar properties (other than the position) 
        associated to the nodes (if present) are imported; 
        all properties are assumed to be scalar (corresponding to column in 
        the input file for nodes with a name other than 'id', 'x', 'y', 'z')

    import_edge_scalar_properties : bool, default: True
        indicates if the scalar properties associated to the edges 
        (if present) are imported; 
        all properties are assumed to be scalar (corresponding to any column
        except the two first ones)

    default_start_id : int, default: 0
        default starting id for nodes if no 'id' column in node file

    verbose : bool, default: True
        indicates if info is printed

    Returns
    -------
    G : networkx.Graph
        networkx graph with nodes and edges given by the input files,
        with:

        - node attributes :
            - `pos_attr` : sequence of 2 or 3 floats, the position
                in 2D or 3D given by the columns 'x', 'y'[, 'z']
                in input file for nodes

            - optional scalar properties given by their name
                in input file for nodes

        - edge attributes :
            - optional scalar properties given by their name
                in input file for edges

    Examples
    --------
        >>> G = kn.io_func.networkx_from_csv('my_graph')
    """
    # Files
    filename_nodes = f'{basename}{suffix_nodes}'
    filename_edges = f'{basename}{suffix_edges}'

    # Read node file (set data frame)
    try:
        nodes_df = pd.read_csv(filename_nodes, delimiter=delimiter_nodes)
    except OSError:
        print(f'IMPORT ERROR: Could not import nodes (file: {filename_nodes})')
        return

    # Read edge file (set data frame)
    try:
        edges_df = pd.read_csv(filename_edges, delimiter=delimiter_edges)
    except OSError:
        print(f'IMPORT ERROR: Could not import edges (file: {filename_edges})')
        return

    # Initialize networkx graph
    G = nx.Graph()

    # Set node id (list)
    if 'id' not in nodes_df.columns:
        nodes_list = list(range(default_start_id, default_start_id + nodes_df.shape[0]))
    else:
        nodes_list = list(nodes_df['id'])

    G.add_nodes_from(nodes_list)

    # Set node position
    if 'x' not in nodes_df.columns:
        raise ValueError(f"Column 'x' not present")
    if 'y' not in nodes_df.columns:
        raise ValueError(f"Column 'y' not present")
    if 'z' in nodes_df.columns:
        pos = dict(zip(nodes_list, nodes_df[['x', 'y', 'z']].values.tolist()))
    else:
        pos = dict(zip(nodes_list, nodes_df[['x', 'y']].values.tolist()))

    nx.set_node_attributes(G, pos, pos_attr)

    if import_node_scalar_properties:
        # Set optional node properties
        prop_name_list = [name for name in nodes_df.columns if name not in ('id', 'x', 'y', 'z')]
        for name in prop_name_list:
            print("Setting node property:", name)
            nx.set_node_attributes(G, dict(zip(nodes_list, nodes_df[name].values.tolist())), name)

    # Set edges (list)
    if edges_df.columns.size < 2:
        raise ValueError("Edge file must have at least 2 columns")
    
    edges_arr = edges_df[edges_df.columns[:2]].values
    if not np.all([id in nodes_list for id in edges_arr.flatten()]):
        raise ValueError("Node id in an edge not existing")
    
    edges_list = list(map(tuple, edges_arr.tolist()))
    # edges_list = [(int(i), int(j)) for i, j in edges_arr] # equivalent
    # edges_list = [(i, j) for i, j in edges_arr.tolist()] # equivalent

    G.add_edges_from(edges_list)

    if import_edge_scalar_properties:
        # Set optional edge properties
        prop_name_list = [name for name in edges_df.columns[2:]]
        for name in prop_name_list:
            nx.set_edge_attributes(G, dict(zip(edges_list, edges_df[name].values.tolist())), name)

    if verbose:
        print("Networkx graph successfully created from csv files !\n")

    return G
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
def networkx_node_scalar_properties_from_csv(
        G, 
        filename,
        properties_name_list=None,
        delimiter=';',  
        verbose=True):
    """
    Imports node scalar properties from a csv file and add them to a networkx graph `G`.

    The file has the columns: 
        - ['id', 'x', 'y', 'z',] '<node_prop0>', ['<node_prop1>', ...]

    Default id (if column 'id' is not given) is integer starting from 0.
    The columns with a name other than 'id' and 'x', 'y'[, 'z'] (giving the 
    position of the nodes), are additional scalar properties (attributes)
    attached to the nodes.

    Note that a line in the file with an id corresponding to none of the nodes
    in the graph `G` is ommitted. 

    Note also that inplace operations are done, i.e. the graph `G` is modified.
    
    Parameters
    ----------
    G : networkx.Graph
        networkx graph

    filename : str
        the name of the input file
    
    properties_name_list : list of strs, optional
        list of property names that are imported;
        by default (`None`): all properties are (except the position 
        ('x', 'y'[, 'z'])) are imported

    delimiter : str, default: ';'
        delimiter used in the file
    
    verbose : bool, default: True
        indicates if info is printed

    Returns
    -------
    G : networkx.Graph
        updated networkx graph

    Examples
    --------
        >>> G = kn.io_func.networkx_node_scalar_properties_from_csv(G, 'my_graph_node_properties.csv')
    """
    # Read node properties file (set data frame)
    try:
        nodes_df = pd.read_csv(filename, delimiter=delimiter)
    except OSError:
        print(f'IMPORT ERROR: Could not import node scalar properties (file: {filename})')
        return

    # Set node id (list)
    if 'id' not in nodes_df.columns:
        nodes_list = list(range(nodes_df.shape[0]))
    else:
        nodes_list = list(nodes_df['id'])

    # Set node properties
    if properties_name_list is None:
        properties_name_list = [name for name in nodes_df.columns if name not in ('id', 'x', 'y', 'z')]

    for name in properties_name_list:
        nx.set_node_attributes(G, dict(zip(nodes_list, nodes_df[name].values.tolist())), name)

    if verbose:
        print("Networkx graph: node scalar properties successfully imported from csv file !\n")

    return G
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
def networkx_edge_scalar_properties_from_csv(
        G, 
        filename,
        properties_name_list=None,
        delimiter=';',  
        verbose=True):
    """
    Imports edge scalar properties from a csv file and add them to a networkx graph `G`.

    The file has the columns: 
        - '<from>', '<to>', '<edge_prop0>'[, '<edge_prop1>', ...]

    where the first two columns, '<from>', '<to>', are the node ids 
    of the two extremities of an edge (link), the next columns are
    additional scalar properties (attributes) attached to the edges.

    Note that a line in the file with a pair of ids ('<from>', '<to>') 
    corresponding to none of the edges in the graph `G` is ommitted. 

    Note also that inplace operations are done, i.e. the graph `G` is modified.
    
    Parameters
    ----------
    G : networkx.Graph
        networkx graph

    filename : str
        the name of the input file
    
    properties_name_list : list of strs, optional
        list of property names that are imported;
        by default (`None`): all properties are imported
    
    delimiter : str, default: ';'
        delimiter used in the file
    
    verbose : bool, default: True
        indicates if info is printed

    Returns
    -------
    G : networkx.Graph
        updated networkx graph

    Examples
    --------
        >>> G = kn.io_func.networkx_edge_scalar_properties_from_csv(G, 'my_graph_edge_properties.csv')
    """
    # Read edge file (set data frame)
    try:
        edges_df = pd.read_csv(filename, delimiter=delimiter)
    except OSError:
        print(f'IMPORT ERROR: Could not import edge scalar properties (file: {filename})')
        return

    # Set edges (list)
    if edges_df.columns.size < 2:
        raise ValueError(f"Edge file must have at least 2 columns")
    
    edges_arr = edges_df[edges_df.columns[:2]].values

    edges_list = list(map(tuple, edges_arr.tolist()))
    # edges_list = [(int(i), int(j)) for i, j in edges_arr] # equivalent
    # edges_list = [(i, j) for i, j in edges_arr.tolist()] # equivalent
    
    # Set edge properties
    if properties_name_list is None:
        properties_name_list = [name for name in edges_df.columns[2:]]

    for name in properties_name_list:
        nx.set_edge_attributes(G, dict(zip(edges_list, edges_df[name].values.tolist())), name)

    if verbose:
        print("Networkx graph: edge scalar properties successfully imported from csv file !\n")

    return G
# ----------------------------------------------------------------------------

# ===== Export functions to csv file(s) for a networkx graph =================
# Added by Julien Straubhaar

# ----------------------------------------------------------------------------
def networkx_to_csv(
        G, 
        basename,
        suffix_nodes='_nodes.csv',
        suffix_edges='_edges.csv',
        delimiter_nodes=';',
        delimiter_edges=';',
        pos_attr='pos',
        edge_name=['from', 'to'],
        export_node_scalar_properties=True,
        export_edge_scalar_properties=True,
        round_node_properties=None,
        round_edge_properties=None,
        verbose=True):
    """
    Exports a networkx graph to two csv files (nodes, and edges (links)).

    The file for nodes has the columns:
        - 'id', 'x', 'y'[, 'z', '<node_prop0>', '<node_prop1>', ...]

    where the column 'id' is the id of the nodes, the columns 'x', 'y'[, 'z']
    give the position of the nodes and the next columns are additional
    scalar properties (attributes) attached to the nodes.
    The column 'z' is omitted for 2D case.

    The file for edges (links) has the columns:
        - '<from>', '<to>'[, '<edge_prop0>', '<edge_prop1>', ...]

    where the first two columns, '<from>', '<to>' (names given by `edge_name`), are the 
    node ids of the two extremities of an edge (link), the next columns (if any) 
    are additional scalar properties (attributes) attached to the edges.

    Note: if a property (attribute) is an array (non-scalar), it is splitted into
    several columns with names '<prop_name>_0', '<prop_name>_1', ..., where '<prop_name>' 
    is the name of the property.

    Parameters
    ----------
    G : networkx.Graph
        networkx graph (that contains at least one node property 
        corresponding to the position of the nodes)

    basename : str
        the base name (prefix) used for the output files

    suffix_nodes : str, default: '_nodes.csv'
        the output file containing the nodes is
        `basename``suffix_nodes`

    suffix_edges : str, default: '_edges.csv'
        the output file containing the edges (links) is
        `basename``suffix_edges`

    delimiter_nodes : str, default:';'
        delimiter used in file for nodes

    delimiter_edges : str, default:';'
        delimiter used in file for edges

    pos_attr : str, default: 'pos'
        name of the node attribute for position

    edge_name : list of 2 str, default: ['from', 'to']
        names of the columns corresponding to the extremities of 
        the edges in the output file for edges

    export_node_scalar_properties : bool, default: True
        indicates if the properties (other than 'pos') associated to 
        the nodes (if present) are exported

    export_edge_scalar_properties : bool, default: True
        indicates if the properties associated to 
        the edges (if present) are exported

    round_node_properties : int or dict-like, optional
        if `int`, number of decimals for rounding the node attribute values;
        if `dict-like`, keys are the node attribute names, and values are
        the number of decimals for rounding the corresponding node attribute values;
        if `None` (default), no rounding is performed

    round_edge_properties : int or dict-like, optional
        if `int`, number of decimals for rounding the edge attribute values;
        if `dict-like`, keys are the edge attribute names, and values are
        the number of decimals for rounding the corresponding edge attribute values;
        if `None` (default), no rounding is performed

    verbose : bool, default: True
        indicates if info is printed

    Examples
    --------
        >>> kn.io_func.networkx_to_csv(G, "MyKarst")
    """

    # Files
    filename_nodes = f'{basename}{suffix_nodes}'
    filename_edges = f'{basename}{suffix_edges}'

    # Data frame for nodes
    if export_node_scalar_properties:
        attr_list = None
    else:
        attr_list = []

    # Get data frame for nodes
    try:
        nodes_df = kn.utils.get_node_data_frame(
                        G, 
                        pos_attr=pos_attr, 
                        include_pos=True, 
                        attr_list=attr_list,
                        round=round_node_properties)
    except:
        raise ValueError(f'Export error : can not set data frame for nodes')

    # Get data frame for edges
    if export_edge_scalar_properties:
        attr_list = None
    else:
        attr_list = []

    try:
        edges_df = kn.utils.get_edge_data_frame(
                        G, 
                        edge_name=edge_name, 
                        attr_list=attr_list,
                        round=round_edge_properties)
    except:
        raise ValueError(f'Export error : can not set data frame for edges')

    # Export
    try:
        nodes_df.to_csv(filename_nodes, index=False, sep=delimiter_nodes)
    except:
        raise ValueError(f'Export error : can not export csv file for nodes ({filename_nodes})')

    try:
        edges_df.to_csv(filename_edges, index=False, sep=delimiter_edges)
    except:
        raise ValueError(f'Export error : can not export csv file for edges ({filename_edges})')

    if verbose:
        print(f"Networkx graph successfully exported to csv files ({filename_nodes}, {filename_edges}) !\n")

    return
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
def networkx_node_scalar_properties_to_csv(
        G, 
        filename,
        properties_name_list=None,
        delimiter=';',
        pos_attr='pos',
        round=None,
        verbose=True):
    """
    Exports node scalar properties from a networkx graph to a csv file.

    The file has the columns:
        - 'id', '<node_prop0>'[, '<node_prop1>', ...]

    where the column 'id' is the id of the nodes and the next columns are 
    additional scalar properties (attributes) attached to the nodes.

    Note: if a property (attribute) is an array (non-scalar), it is splitted into
    several columns with names '<prop_name>_0', '<prop_name>_1', ..., where '<prop_name>' 
    is the name of the property.

    Parameters
    ----------
    G : networkx.Graph
        networkx graph

    filename : str
        the name of the output file
    
    properties_name_list : list of strs, optional
        list of property names that are exported;
        by default (`None`): all properties (except the position)
        are exported and assumed to be scalar

    delimiter : str, default: ';'
        delimiter used in the file
    
    pos_attr : str, default: 'pos'
        name of the node attribute for position
    
    round : int or dict-like, optional
        if `int`, number of decimals for rounding the node attribute values;
        if `dict-like`, keys are the node attribute names, and values are
        the number of decimals for rounding the corresponding node attribute values;
        if `None` (default), no rounding is performed
        
    verbose : bool, default: True
        indicates if info is printed

    Examples
    --------
        >>> kn.io_func.networkx_node_scalar_properties_to_csv(G, "MyKarst_node_properties.csv")
    """
    # Get data frame
    try:
        df = kn.utils.get_node_data_frame(
                        G, 
                        pos_attr=pos_attr, 
                        include_pos=False, 
                        attr_list=properties_name_list,
                        round=round)
    except:
        raise ValueError(f'Export error : can not set data frame for node properties')

    # Rounding
    if round is not None:
        df = df.round(round)

    # Export
    try:
        df.to_csv(filename, index=False, sep=delimiter)
    except:
        raise ValueError(f'Export error : can not export csv file for node properties ({filename})')

    if verbose:
        print(f"Networkx graph: node scalar properties successfully exported to csv file ({filename}) !\n")

    return
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
def networkx_edge_scalar_properties_to_csv(
        G, 
        filename,
        properties_name_list=None,
        delimiter=';',
        edge_name=['from', 'to'],
        round=None,
        verbose=True):
    """
    Exports edge scalar properties from a networkx graph to a csv file.

    The file has the columns:
        - '<from>', '<to>', '<edge_prop0>'[, '<edge_prop1>', ...]

    where the first two columns, '<from>', '<to>' (names given by `edge_name`), are the 
    node ids of the two extremities of an edge (link), the next columns (if any) 
    are additional scalar properties (attributes) attached to the edges.

    Note: if a property (attribute) is an array (non-scalar), it is splitted into
    several columns with names '<prop_name>_0', '<prop_name>_1', ..., where '<prop_name>' 
    is the name of the property.
    
    Parameters
    ----------
    G : networkx.Graph
        networkx graph

    filename : str
        the name of the output file
    
    properties_name_list : list of strs, optional
        list of property names that are exported;
        by default (`None`): all properties are exported and assumed to 
        be scalar
    
    delimiter : str, default: ';'
        delimiter used in the file
    
    edge_name : list of 2 str, default: ['from', 'to']
        names of the columns corresponding to the extremities of 
        the edges in the output file for edges

    round_edge_properties : int or dict-like, optional
        if `int`, number of decimals for rounding the edge attribute values;
        if `dict-like`, keys are the edge attribute names, and values are
        the number of decimals for rounding the corresponding edge attribute values;
        if `None` (default), no rounding is performed

    verbose : bool, default: True
        indicates if info is printed

    Examples
    --------
        >>> kn.io_func.networkx_edge_scalar_properties_to_csv(G, "MyKarst_edge_properties.csv")
    """
    # Get data frame
    try:
        df = kn.utils.get_edge_data_frame(
                        G, 
                        edge_name=edge_name, 
                        attr_list=properties_name_list,
                        round=round)
    except:
        raise ValueError(f'Export error : can not set data frame for edge properties')

    # Rounding
    if round is not None:
        df = df.round(round)

    # Export
    try:
        df.to_csv(filename, index=False, sep=delimiter)
    except:
        raise ValueError(f'Export error : can not export csv file for edge properties ({filename})')

    if verbose:
        print(f"Networkx graph: edge scalar properties successfully exported to csv file ({filename}) !\n")

    return
# ----------------------------------------------------------------------------
