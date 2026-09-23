"""
Module `_view.plot`
-------------------

The module `_view.plot` is accessible as `view`. It contains functions
to plot networkx graphs in 2D and 3D.
"""

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.colors

import karstnet as kn


# ===== Functions for plot in 2D (networkx graph) ============================
# Julien Straubhaar

# ----------------------------------------------------------------------------
def plot_graph2d(
        G,
        ax=None,
        node_attr=None,
        edge_attr=None,
        pos_attr='pos',
        proj_mode='xy',
        centering=False,
        plot_nodes=True,
        plot_edges=True,
        add_options_nodes=None,
        add_options_edges=None,
        show_node_colorbar=True,
        show_edge_colorbar=True,
        options_node_colorbar=None,
        options_edge_colorbar=None,
        show_axis_labels=True,
        xlabel=None,
        ylabel=None,
        options_xlabel=None,
        options_ylabel=None,
        show_ticks=False,
        options_xaxis_ticks=None,
        options_yaxis_ticks=None,
        axis_equal=False,
        **kwargs):
    """
    Plot a graph in 2D.

    The graph `G` in input is assumed to have nodes localized in 2D or 3D space.

    Parameters
    ----------
    G : networkx.Graph object
        graph, with 2D or 3D position as node attribute
    
    ax : matplotlib.axes._axes.Axes, optional
        axis in which the plot is done;
        by default (`None`): the current axis, retrieved by `ax = plt.gca()` is used

    node_attr : str, optional
        node attribute to be plotted; this should be a scalar attribute

    edge_attr : str, optional
        edge attribute to be plotted; this should be a scalar attribute

    pos_attr : str, default: 'pos'
        name of the attribute attached to nodes defining the position as a 
        sequence of 2 or 3 floats; if position are in 2D, a coordinate 
        of zero is set in the 3rd dimension
    
    proj_mode : str {'xy', 'xz', 'yz', 'pca_xy', 'pca_xz', 'pca_yz'}, default: 'xy'
        string defining the projection mode:
        
        - 'xy', 'xz', 'yz': orthognal projection onto a 2D section spanned \
        by the two given axes of original 3D axes system (xy, or xz, or yz)
        - 'pca_12', 'pca_13', 'pca_23': orthognal projection onto a 2D section spanned \
        by the two given axes (12, 13, or 23) of the 3D axes system retrieved from the \
        principal component analysis (pca); this axes system is orthonormal and the \
        principal axes 1 2 and 3 are ordered such that the variance of the coordinates \
        of the node position projected onto these axes are decreasing
    
    centering: bool, default: False
        - if `True`: the coordinates of the position are centered (mean equals to zero)
        - if `False`: the origin of the system is not changed

    plot_nodes : bool, default: True
        - if `True`: the nodes are plotted
        - if `False`: the nodes are not plotted (other parameters about \
        plotting nodes will be ignored)

    plot_edges : bool, default: True
        - if `True`: the edges are plotted
        - if `False`: the edges are not plotted (other parameters about \
        plotting edges will be ignored)

    add_options_nodes : dict, optional
        dictionary of additional options for plotting nodes when a node attribute
        is specified (`node_attr` is not `None`), using the function 
        `nx.draw_networkx_nodes` ; note that valid keyword arguments from `kwargs`
        are also used

    add_options_edges : dict, optional
        dictionary of additional options for plotting edges when an edge attribute
        is specified (`edge_attr` is not `None`), using the function 
        `nx.draw_networkx_edges` ; note that valid keyword arguments from `kwargs`
        are also used

    show_node_colorbar : bool, default: True
        indicates if the colorbar for the node attribute is displayed;
        used if `node_attr` is not None
        
    show_edge_colorbar : bool, default: True
        indicates if the colorbar for the edge attribute is displayed;
        used if `edge_attr` is not None
        
    options_node_colorbar : dict, optional
        options (keyword arguments) passed to the function `plt.colorbar`
        for plotting the colorbar for the node attribute
        (e.g. possible keys: 'orientation', 'location', 'shrink', etc.);
        used if `node_attr` is not `None` and `show_node_colorbar=True`

    options_edge_colorbar : dict, optional
        options (keyword arguments) passed to the function `plt.colorbar`
        for plotting the colorbar for the edge attribute
        (e.g. possible keys: 'orientation', 'location', 'shrink', etc.);
        used if `edge_attr` is not `None` and `show_edge_colorbar=True`

    show_axis_labels : bool, default: True
        if `True`, axis labels are displayed
    
    xlabel : str, optional
        name for x axis label; by default (`None`): automatically set;
        note: if specified, `ylabel` must also be specified
    
    ylabel : str, optional
        name for y axis label; by default (`None`): automatically set
        note: if specified, `xlabel` must also be specified

    options_xlabel : dict, optional
        options (keyword arguments) passed to the function `plt.xlabel`
        (e.g. possible keys: 'fontsize', 'weight', etc.)

    options_ylabel : dict, optional
        options (keyword arguments) passed to the function `plt.ylabel`
        (e.g. possible keys: 'fontsize', 'weight', etc.)

    show_ticks : bool, default: False
        if `True`, ticks are displayed along axes

    options_xaxis_ticks : dict, optional
        options (keyword arguments) passed to the function 
        `plt.gca().xaxis.set_tick_params`
        (e.g. possible keys: 'rotation', 'labelsize', 'color', 
        'labelcolor', etc.)

    options_yaxis_ticks : dict, optional
        options (keyword arguments) passed to the function 
        `plt.gca().yaxis.set_tick_params`
        (e.g. possible keys: 'rotation', 'labelsize', 'color', 
        'labelcolor', etc.)

    axis_equal : bool, default: False
        if `True`, same scale is used for both axes

    kwargs : dict
        keyword arguments used in function `nx.draw_networkx`, 
        `nx.draw_networkx_nodes` or `nx.draw_networkx_edges`
        (relevant / valid keywords arguments are considered);
        first, `nx.draw_networkx` is called (`with_labels` set to
        `False` if not specified), then 
        `nx.draw_networkx_nodes` if `node_attr` is not `None`, 
        and finally
        `nx.draw_networkx_edges` if `edge_attr` is not `None`
        
        Node size:
        
        - the keyword argument `node_size` can be set to a float (or int), \
        or to a sequence of n_nodes values (floats or ints), where n_nodes \
        is the number of nodes in the entire graph `G` and where the values \
        are the desired size for nodes given by `G.nodes()`, even if the keywork \
        argument `nodelist` (specifying the node to be plotted) is a list \
        of length < n_nodes (the node sizes used are then retrieved from the \
        sequence of node sizes defined for all graph nodes)

        Colors for nodes:

        - if `node_attr` is not given (`None`): `node_color` is used
        - if `node_attr` is given (not `None`): `cmap` is used and \
        and for nan values, `cmap.get_cbad()` or 'node_color' \
        (if specified, prevails) is used for nan values; moreover \
        `vmin` and / or `vmax` may be specified to limit the range of \
        colorbar, in this case `cmap.get_under()`, resp. `cmap.get_over()`, \
        are used for values out of the range

        Colors for edges:

        - if `edge_attr` is not given (`None`): `edge_color` is used
        - if `edge_attr` is given (not `None`): `edge_cmap` is used and \
        and for nan values, `edge_cmap.get_cbad()` or 'edge_color' \
        (if specified, prevails) is used for nan values; moreover \
        `edge_vmin` and / or `edge_vmax` may be specified to limit the range of \
        colorbar, in this case `edge_cmap.get_under()`, resp. `edge_cmap.get_over()`, \
        are used for values out of the range

    Examples
    --------
        >>> kn.view.plot_graph2d(G)
    """
    import inspect

    if proj_mode == 'xy':
        pca_mode = False
        axes_ind = np.array([0, 1])
    elif proj_mode == 'xz':
        pca_mode = False
        axes_ind = np.array([0, 2])
    elif proj_mode == 'yz':
        pca_mode = False
        axes_ind = np.array([1, 2])
    elif proj_mode == 'pca_12':
        pca_mode = True
        axes_ind = np.array([0, 1])
    elif proj_mode == 'pca_13':
        pca_mode = True
        axes_ind = np.array([0, 2])
    elif proj_mode == 'pca_23':
        pca_mode = True
        axes_ind = np.array([1, 2])
    else:
        raise ValueError('Parameter `proj_mode` not valid')

    # Get dictionary of 3d node position
    pos = kn.utils.get_pos3d(G, pos_attr=pos_attr)

    if len(pos) == 0:
        print('No node in graph!')
        return
    
    # Get the 2d array of node positions:
    # - pos_arr[i] : node position of the node list(pos.keys())[i]
    pos_arr = np.asarray(list(pos.values()))

    if centering:
        # Substact the mean position
        pos_arr = pos_arr - pos_arr.mean(axis=0)

    if pca_mode:
        # Do the pca, get the pca axes and pca variances
        pca_axes, pca_var = kn.utils.pca_of_node_position(G)
        # Set node position in the pca axes system
        pos_arr = pos_arr.dot(pca_axes)
        
    # Retain specified axes
    pos_arr = pos_arr[:, axes_ind]    

    # Set dictionary of 2d position to be plotted
    pos = {k:pos_arr[i] for i, k in enumerate(pos.keys())}

    # Draw graph (not accounting for node / edge attribute)
    kwargs_base = kwargs.copy()
    if node_attr is not None and plot_nodes:
        if 'cmap' in kwargs_base.keys():
            if 'node_color' not in kwargs_base.keys():
                # Get bad color from cmap
                cmap = plt.get_cmap(kwargs_base['cmap'])
                cbad = cmap.get_bad()
                # cbad = matplotlib.colors.to_rgba(cbad)
                cbad = matplotlib.colors.to_hex(cbad)
                kwargs_base['node_color'] = cbad
        
            # Remove 'cmap' key
            del(kwargs_base['cmap'])

        # Remove also related option 'vmin' / 'vmax' if specified
        if 'vmin' in kwargs_base.keys():
            del(kwargs_base['vmin'])

        if 'vmax' in kwargs_base.keys():
            del(kwargs_base['vmax'])

        # Set list of nodes with nan value for the specified node attribute
        node_attr_dict = nx.get_node_attributes(G, node_attr, default=np.nan)
        nodelist = [k for k, v in node_attr_dict.items() if np.isnan(v)]
        if 'nodelist' in kwargs_base.keys():
            nodelist = [u for u in kwargs_base['nodelist'] if u in nodelist]

        kwargs_base['nodelist'] = nodelist

        if 'node_size' in kwargs_base.keys():
            if hasattr(kwargs_base['node_size'], '__len__'):
                kwargs_base['node_size'] = [kwargs_base['node_size'][i] for i, u in enumerate(G.nodes()) if u in nodelist]

    elif node_attr is None and plot_nodes:
        if 'node_size' in kwargs_base.keys() and hasattr(kwargs_base['node_size'], '__len__') and 'nodelist' in kwargs_base.keys():
            nodelist = kwargs_base['nodelist']
            kwargs_base['node_size'] = [kwargs_base['node_size'][i] for i, u in enumerate(G.nodes()) if u in nodelist]
                
    if edge_attr is not None and plot_edges:
        if 'edge_cmap' in kwargs_base.keys():
            if 'edge_color' not in kwargs_base.keys():
                # Get bad color from edge_cmap
                edge_cmap = plt.get_cmap(kwargs_base['edge_cmap'])
                edge_cbad = edge_cmap.get_bad()
                # edge_cbad = matplotlib.colors.to_rgba(edge_cbad)
                edge_cbad = matplotlib.colors.to_hex(edge_cbad)
                kwargs_base['edge_color'] = edge_cbad

            # Remove 'edge_cmap' key
            del(kwargs_base['edge_cmap'])

        # Set list of edges with nan value for the specified edge attribute
        edge_attr_dict = nx.get_edge_attributes(G, edge_attr, default=np.nan)
        edgelist = [k for k, v in edge_attr_dict.items() if np.isnan(v)]
        kwargs_base['edgelist'] = edgelist

    if not plot_nodes:
        # Set node size to zero
        kwargs_base['node_size'] = 0

    if not plot_edges:
        # Set edge width to zero
        kwargs_base['width'] = 0

    if 'with_labels' not in kwargs_base.keys():
        # Set to False
        kwargs_base['with_labels'] = False

    nx.draw_networkx(G, pos=pos, **kwargs_base)

    if node_attr is not None and plot_nodes:
        if add_options_nodes is None:
            add_options_nodes = {}
        
        # Get valid kwargs for nx.draw_networkx_nodes
        all_args_name = set([val.name for val in inspect.signature(nx.draw_networkx_nodes).parameters.values()])
        kwargs_nodes_keys = all_args_name.intersection(kwargs.keys())
        kwargs_nodes = {k:kwargs[k] for k in kwargs_nodes_keys}
        
        # node_attr_dict = nx.get_node_attributes(G, node_attr, default=np.nan) # already defined above

        # Set (update) node_color option
        kwargs_nodes['node_color'] = np.asarray(list(node_attr_dict.values()))

        # Set (update) node_size option
        if 'nodelist' in kwargs_nodes.keys():
            node_attr_dict = {k:v for k, v in node_attr_dict.items() if k in kwargs_nodes['nodelist']}
            if 'node_size' in kwargs_nodes.keys():
                if hasattr(kwargs_nodes['node_size'], '__len__'):
                    kwargs_nodes['node_size'] = [kwargs_nodes['node_size'][i] for i, u in enumerate(G.nodes()) if u in node_attr_dict.keys()]
        
        # Plot nodes
        im_nodes = nx.draw_networkx_nodes(G, pos=pos, **add_options_nodes, **kwargs_nodes)
        if show_node_colorbar:
            if options_node_colorbar is None:
                options_node_colorbar = {}

            if 'label' not in options_node_colorbar.keys():
                options_node_colorbar['label'] = f'node: {node_attr}'

            plt.colorbar(im_nodes, **options_node_colorbar)

    if edge_attr is not None and plot_edges:
        if add_options_edges is None:
            add_options_edges = {}
        
        # Get valid kwargs for nx.draw_networkx_edges
        all_args_name = set([val.name for val in inspect.signature(nx.draw_networkx_edges).parameters.values()])
        kwargs_edges_keys = all_args_name.intersection(kwargs.keys())
        kwargs_edges = {k:kwargs[k] for k in kwargs_edges_keys}

        # Set edge_cmap        
        if 'edge_cmap' not in kwargs_edges.keys():
            edge_cmap = 'viridis' # default color map
        else:
            edge_cmap = kwargs_edges['edge_cmap']

        edge_cmap = plt.get_cmap(edge_cmap)
        if 'edge_color' in kwargs_edges.keys():
            # Set edge_color as bad color
            edge_cbad = kwargs_edges['edge_color']
            edge_cmap.set_bad(edge_cbad)

        # Note: string for color map is not accepted in `nx.draw_networkx_edges`
        kwargs_edges['edge_cmap'] = edge_cmap

        # Set (update) edge_color option
        # edge_attr_dict = nx.get_edge_attributes(G, edge_attr, default=np.nan) # already defined above
        if 'edgelist' in kwargs_edges.keys():
            edge_attr_dict = {k:v for k, v in edge_attr_dict.items() if k in kwargs_edges['edgelist']}
        
        kwargs_edges['edge_color'] = np.asarray(list(edge_attr_dict.values()))

        # Plot edges
        im_edges = nx.draw_networkx_edges(G, pos=pos, **add_options_edges, **kwargs_edges)
        if show_edge_colorbar:
            if options_edge_colorbar is None:
                options_edge_colorbar = {}

            if 'label' not in options_edge_colorbar.keys():
                options_edge_colorbar['label'] = f'edge: {edge_attr}'

            plt.colorbar(im_edges, **options_edge_colorbar)

    if show_axis_labels:
        if options_xlabel is None:
            options_xlabel = {}
        
        if options_ylabel is None:
            options_ylabel = {}

        if xlabel is None or ylabel is None:
            if pca_mode:
                axes_label_name = ['pca-1st', 'pca-2nd', 'pca-3rd']
                pca_var_percent = pca_var/pca_var.sum() * 100.0
                # xlabel = axes_label_name[axes_ind[0]] + f', {pca_var_percent[axes_ind[0]]:5.1f}% of tot. var.' # + ', vec = (' + ', '.join([f'{v:.2g}' for v in pca_axes[:, axes_ind[0]]]) + ')' 
                # ylabel = axes_label_name[axes_ind[1]] + f', {pca_var_percent[axes_ind[1]]:5.1f}% of tot. var.' # + ', vec = (' + ', '.join([f'{v:.2g}' for v in pca_axes[:, axes_ind[1]]]) + ')' 
                xlabel = axes_label_name[axes_ind[0]] + ': (' + ', '.join([f'{v:.2g}' for v in pca_axes[:, axes_ind[0]]]) + ')' + f'\n{pca_var_percent[axes_ind[0]]:5.1f}% of tot. var.'
                ylabel = axes_label_name[axes_ind[1]] + ': (' + ', '.join([f'{v:.2g}' for v in pca_axes[:, axes_ind[1]]]) + ')' + f'\n{pca_var_percent[axes_ind[1]]:5.1f}% of tot. var.'
            else:
                axes_label_name = ['x', 'y', 'z']
                xlabel = axes_label_name[axes_ind[0]]
                ylabel = axes_label_name[axes_ind[1]]
        plt.xlabel(xlabel, **options_xlabel)
        plt.ylabel(ylabel, **options_ylabel)

    if show_ticks:
        if options_xaxis_ticks is None:
            options_xaxis_ticks = {}

        if options_yaxis_ticks is None:
            options_yaxis_ticks = {}

        if ax is None:
            ax = plt.gca()
        ax.tick_params(left=True, bottom=True, labelleft=True, labelbottom=True)
        ax.xaxis.set_tick_params(**options_xaxis_ticks)
        ax.yaxis.set_tick_params(**options_yaxis_ticks)

    if axis_equal:
        plt.axis('equal')
# ----------------------------------------------------------------------------

# ===== Functions for plot in 3D (networkx graph) based on matplotlib ========
# Julien Straubhaar

# ----------------------------------------------------------------------------
def plot_graph3d(
        G,
        ax=None,
        node_attr=None,
        edge_attr=None,
        pos_attr='pos',
        mode='xyz',
        centering=False,
        plot_nodes=True,
        plot_edges=True,
        options_nodes=None,
        options_edges=None,
        show_node_colorbar=True,
        show_edge_colorbar=True,
        options_node_colorbar=None,
        options_edge_colorbar=None,
        show_axis_labels=True,
        xlabel=None,
        ylabel=None,
        zlabel=None,
        options_xlabel=None,
        options_ylabel=None,
        options_zlabel=None,
        options_xaxis_ticks=None,
        options_yaxis_ticks=None,
        options_zaxis_ticks=None,
        zrotation=30.0, 
        xyrotation=60.0,
        box_aspect=None,
        aspect_equal=False):
    """
    Plot a graph in 3D, based on `mpl_toolkits.mplot3d.axes3d.Axes3D`.

    The graph `G` in input is assumed to have nodes localized in 2D or 3D space.

    Parameters
    ----------
    G : networkx.Graph object
        graph, with 2D or 3D position as node attribute
    
    ax : mpl_toolkits.mplot3d.axes3d.Axes3D, optional
        axis in which the plot is done;
        by default (`None`): the current figure is retrieved or created with
        `fig = plt.gcf()`, then axis are created with 
        `ax = fig.add_subplot(projection='3d')`

    node_attr : str, optional
        node attribute to be plotted; this should be a scalar attribute

    edge_attr : str, optional
        edge attribute to be plotted; this should be a scalar attribute

    pos_attr : str, default: 'pos'
        name of the attribute attached to nodes defining the position as a 
        sequence of 2 or 3 floats; if position are in 2D, a coordinate 
        of zero is set in the 3rd dimension
    
    mode : str {'xyz', 'pca'}, default: 'xyz'
        - if `mode='xyz'`: original 3D axes system is used
        - if `mode='pca'`: the 3D axes system retrieved from the \
        principal component analysis (pca) is used; this axes system is \
        orthonormal and the principal axes are ordered such that the variance \
        of the coordinates of the node position projected onto these axes are \
        decreasing
            
    centering: bool, default: False
        - if `True`: the coordinates of the position are centered (mean equals to zero)
        - if `False`: the origin of the system is not changed

    plot_nodes : bool, default: True
        - if `True`: the nodes are plotted
        - if `False`: the nodes are not plotted (other parameters about \
        plotting nodes will be ignored)

    plot_edges : bool, default: True
        - if `True`: the edges are plotted
        - if `False`: the edges are not plotted (other parameters about \
        plotting edges will be ignored)

    options_nodes : dict, optional
        dictionary of options that is passed to `ax.scatter` to draw 
        the nodes; e.g. possible keys: 's' (size), etc.
        
        for size:
        
        - options `s` can be set to a float (or int), \
        or to a sequence of n_nodes values (floats or ints), where n_nodes \
        is the number of nodes in the entire graph `G` and where the values \
        are the desired size for nodes given by `G.nodes()`
        
        for colors:

        - if `node_attr` is not given (`None`): 'c' or 'color' is used
        - if `node_attr` is given (not `None`): `cmap` is used and \
        and for nan values, `cmap.get_cbad()` or 'c' or 'color' \
        (if specified, prevails) is used for nan values; moreover \
        'vmin' and / or 'vmax' may be specified to limit the range of \
        colorbar, in this case `cmap.get_under()`, resp. `cmap.get_over()`, \
        are used for values out of the range

    options_edges : dict, optional
        dictionary of options that is passed to `ax.plot` to draw each individual
        edge; e.g. possible keys: 'linestyle' or 'ls', 'linewidth' or 'lw', etc.

        for colors:
        
        - if `edge_attr` is not given (`None`): 'c' or 'color' is used
        - if `edge_attr` is given (not `None`): `cmap` is used and \
        and for nan values, `cmap.get_cbad()` or 'c' or 'color' \
        (if specified, prevails) is used for nan values; moreover \
        'vmin' and / or 'vmax' may be specified to limit the range of \
        colorbar, in this case `cmap.get_under()`, resp. `cmap.get_over()`, \
        are used for values out of the range

    show_node_colorbar : bool, default: True
        used only when `node_attr` is not None: if `True`, colorbar for node 
        attribute is displayed
        
    show_edge_colorbar : bool, default: True
        used only when `edge_attr` is not None: if `True`, colorbar for edge 
        attribute is displayed

    options_node_colorbar : dict, optional
        options (keyword arguments) passed to the function `plt.colorbar`
        for plotting the colorbar for the node attribute
        (e.g. possible keys: 'orientation', 'location', 'shrink', etc.);
        used if `node_attr` is not `None` and `show_node_colorbar=True`
        
    options_edge_colorbar : dict, optional
        options (keyword arguments) passed to the function `plt.colorbar`
        for plotting the colorbar for the edge attribute
        (e.g. possible keys: 'orientation', 'location', 'shrink', etc.);
        used if `edge_attr` is not `None` and `show_edge_colorbar=True`

    show_axis_labels : bool, default: True
        if `True`, axis labels are displayed

    xlabel : str, optional
        name for x axis label; by default (`None`): automatically set;
        note: if specified, `ylabel` and `zlabel` must also be specified
    
    ylabel : str, optional
        name for y axis label; by default (`None`): automatically set
        note: if specified, `xlabel` and `zlabel` must also be specified

    zlabel : str, optional
        name for z axis label; by default (`None`): automatically set
        note: if specified, `xlabel` and `ylabel` must also be specified

    options_xlabel : dict, optional
        options (keyword arguments) passed to the function `ax.set_xlabel`
        (e.g. possible keys: 'fontsize', 'weight', etc.)

    options_ylabel : dict, optional
        options (keyword arguments) passed to the function `ax.set_ylabel`
        (e.g. possible keys: 'fontsize', 'weight', etc.)

    options_zlabel : dict, optional
        options (keyword arguments) passed to the function `ax.set_zlabel`
        (e.g. possible keys: 'fontsize', 'weight', etc.)

    options_xaxis_ticks : dict, optional
        options (keyword arguments) passed to the function 
        `ax.xaxis.set_tick_params`
        (e.g. possible keys: 'rotation', 'labelsize', 'color', 
        'labelcolor', etc.)

    options_yaxis_ticks : dict, optional
        options (keyword arguments) passed to the function 
        `ax.yaxis.set_tick_params`
        (e.g. possible keys: 'rotation', 'labelsize', 'color', 
        'labelcolor', etc.)

    options_zaxis_ticks : dict, optional
        options (keyword arguments) passed to the function 
        `ax.zaxis.set_tick_params`
        (e.g. possible keys: 'rotation', 'labelsize', 'color', 
        'labelcolor', etc.)

    box_aspect : sequence of 3 floats, optional
        if given, specify the box aspect ratios

    aspect_equal : bool, default: False
        if `True`, same scale is used for the three axes (this prevails over
        `box_aspect`)

    zrotation : float, default: 30.0
        angle in degrees between horizontal plane and viewpoint

    xyrotation : float, default: 60.0
        angle in degree for the horizontal rotation of the viewpoint;
        if `xyrotation=0`, the view is from the South toward North
      
    Examples
    --------
       >>> kn.view.plot_graph3d(G)
    """
    from mpl_toolkits.mplot3d.axes3d import Axes3D

    if mode == 'xyz':
        pca_mode = False
    elif mode == 'pca':
        pca_mode = True
    else:
        raise ValueError('Parameter `mode` not valid')

    # Get dictionary of 3d node position
    pos = kn.utils.get_pos3d(G, pos_attr=pos_attr)

    if len(pos) == 0:
        print('No node in graph!')
        return

    # Get the 2d array of node positions:
    # - pos_arr[i] : node position of the node list(pos.keys())[i]
    pos_arr = np.asarray(list(pos.values()))

    if centering:
        # Substact the mean position
        pos_arr = pos_arr - pos_arr.mean(axis=0)

    if pca_mode:
        # Do the pca, get the pca axes and pca variances
        pca_axes, pca_var = kn.utils.pca_of_node_position(G)
        # Set node position in the pca axes system
        pos_arr = pos_arr.dot(pca_axes)
        
    # Set dictionary of 3d position to be plotted
    pos = {k:pos_arr[i] for i, k in enumerate(pos.keys())}

    if options_edges is None:
        options_edges = {}
    
    if options_nodes is None:
        options_nodes = {}

    # Get / set current figure and axis for plotting
    fig = plt.gcf()
    if ax is not None:
        if not isinstance(ax, Axes3D):
            raise TypeError('`ax` should be an instance of `mpl_toolkits.mplot3d.axes3d.Axes3D`')
    else:
        ax = fig.add_subplot(projection='3d')

    # Plot nodes
    if plot_nodes:
        if node_attr is not None:
            v = np.asarray(list(nx.get_node_attributes(G, node_attr, default=np.nan).values()))

            # Get bad color from cmap
            if 'cmap' in options_nodes.keys():
                cmap = plt.get_cmap(options_nodes['cmap'])
            else:
                cmap = plt.get_cmap('viridis') # default cmap
            
            cbad = cmap.get_bad()

            if 'c' in options_nodes.keys():
                # cbad = matplotlib.colors.to_rgba(options_nodes['c'])
                cbad = matplotlib.colors.to_hex(options_nodes['c'])
                del(options_nodes['c'])

            if 'color' in options_nodes.keys():
                # cbad = matplotlib.colors.to_rgba(options_nodes['color'])
                cbad = matplotlib.colors.to_hex(options_nodes['color'])
                del(options_nodes['color'])

            # Set (update) s (node size) option
            if 's' in options_nodes.keys() and hasattr(options_nodes['s'], '__len__'):
                s_not_nan = [options_nodes['s'][i] for i, vi in enumerate(v) if ~np.isnan(vi)]
                s_nan = [options_nodes['s'][i] for i, vi in enumerate(v) if np.isnan(vi)]
                options_nodes['s'] = s_not_nan
                s_has_been_updated = True
            else:
                s_has_been_updated = False

            im_nodes =  ax.scatter(*pos_arr[~np.isnan(v)].T, c=v[~np.isnan(v)], **options_nodes)

            if np.isnan(v).any():
                if 'cmap' in options_nodes.keys():
                    del(options_nodes['cmap'])
                
                if 'vmin' in options_nodes.keys():
                    del(options_nodes['vmin'])

                if 'vmax' in options_nodes.keys():
                    del(options_nodes['vmax'])
                
                if s_has_been_updated:
                    options_nodes['s'] = s_nan

                ax.scatter(*pos_arr[np.isnan(v), :].T, color=cbad, **options_nodes)
            
            if show_node_colorbar:
                if options_node_colorbar is None:
                    options_node_colorbar = {}

                if 'label' not in options_node_colorbar.keys():
                    options_node_colorbar['label'] = f'node: {node_attr}'

                fig.colorbar(im_nodes, ax=ax, **options_node_colorbar)
        else:
            ax.scatter(*pos_arr.T, **options_nodes)
        
    # Plot edges
    if plot_edges:
        if edge_attr is not None:
            v = np.asarray(list(nx.get_edge_attributes(G, edge_attr, default=np.nan).values()))

            if 'cmap' not in options_edges.keys():
                cmap = 'viridis' # default cmap
            else:
                cmap = options_edges['cmap']
                del(options_edges['cmap'])

            if 'vmin' not in options_edges.keys():
                vmin = np.nanmin(v)
            else:
                vmin = options_edges['vmin']
                del(options_edges['vmin'])

            if 'vmax' not in options_edges.keys():
                vmax =  np.nanmax(v)
            else:
                vmax = options_edges['vmax']
                del(options_edges['vmax'])

            cmap = plt.get_cmap(cmap)
            # cbad = matplotlib.colors.to_rgba(cmap.get_bad())
            cbad = matplotlib.colors.to_hex(cmap.get_bad())
            if 'c' in options_edges.keys():
                # cbad = matplotlib.colors.to_rgba(options_edges['c'])
                cbad = matplotlib.colors.to_hex(options_edges['c'])
                del(options_edges['c'])

            if 'color' in options_edges.keys():
                # cbad = matplotlib.colors.to_rgba(options_edges['color'])
                cbad = matplotlib.colors.to_hex(options_edges['color'])
                del(options_edges['color'])

            t = 1.0/(vmax - vmin)
            edge_col = []
            for vi in v:
                if ~ np.isnan(vi):
                    edge_col.append(cmap(t*(vi-vmin)))
                else:
                    edge_col.append(cbad)

            for i, e in enumerate(G.edges()):
                x = np.array((pos[e[0]][0], pos[e[1]][0]))
                y = np.array((pos[e[0]][1], pos[e[1]][1]))
                z = np.array((pos[e[0]][2], pos[e[1]][2]))
                
                # Plot the connecting lines
                ax.plot(x, y, z, color=edge_col[i], **options_edges)
        
            if show_edge_colorbar:
                # Fake scatter plot (with point size s=0) to retrieve mappable image to draw colorbar
                tmpx = len(v)*[pos_arr[0, 0]]
                tmpy = len(v)*[pos_arr[0, 1]]
                tmpz = len(v)*[pos_arr[0, 2]]
                im_tmp = ax.scatter(tmpx, tmpy, tmpz, s=0, c=v, cmap=cmap, vmin=vmin, vmax=vmax)

                if options_edge_colorbar is None:
                    options_edge_colorbar = {}

                if 'label' not in options_edge_colorbar.keys():
                    options_edge_colorbar['label'] = f'edge: {edge_attr}'

                fig.colorbar(im_tmp, ax=ax, **options_edge_colorbar)

        else:
            if 'c' not in options_edges.keys() and 'color' not in options_edges.keys():
                # Set a default color: black
                options_edges['color'] = 'black'

            for e in G.edges():
                x = np.array((pos[e[0]][0], pos[e[1]][0]))
                y = np.array((pos[e[0]][1], pos[e[1]][1]))
                z = np.array((pos[e[0]][2], pos[e[1]][2]))

                # Plot the connecting lines
                ax.plot(x, y, z, **options_edges)

    # Set the view
    ax.view_init(zrotation, -xyrotation - 90)

    if show_axis_labels:
        if options_xlabel is None:
            options_xlabel = {}
        
        if options_ylabel is None:
            options_ylabel = {}

        if options_zlabel is None:
            options_zlabel = {}

        if xlabel is None or ylabel is None or zlabel is None:
            if pca_mode:
                axes_label_name = ['pca-1st', 'pca-2nd', 'pca-3rd']
                pca_var_percent = pca_var/pca_var.sum() * 100.0
                # xlabel = axes_label_name[0] + f', {pca_var_percent[0]:5.1f}% of tot. var.' # + ', vec = (' + ', '.join([f'{v:.2g}' for v in pca_axes[:, 0]]) + ')' 
                # ylabel = axes_label_name[1] + f', {pca_var_percent[1]:5.1f}% of tot. var.' # + ', vec = (' + ', '.join([f'{v:.2g}' for v in pca_axes[:, 1]]) + ')' 
                # zlabel = axes_label_name[2] + f', {pca_var_percent[2]:5.1f}% of tot. var.' # + ', vec = (' + ', '.join([f'{v:.2g}' for v in pca_axes[:, 2]]) + ')' 
                xlabel = axes_label_name[0] + ': (' + ', '.join([f'{v:.2g}' for v in pca_axes[:, 0]]) + ')' + f'\n{pca_var_percent[0]:5.1f}% of tot. var.' 
                ylabel = axes_label_name[1] + ': (' + ', '.join([f'{v:.2g}' for v in pca_axes[:, 1]]) + ')' + f'\n{pca_var_percent[1]:5.1f}% of tot. var.' 
                zlabel = axes_label_name[2] + ': (' + ', '.join([f'{v:.2g}' for v in pca_axes[:, 2]]) + ')' + f'\n{pca_var_percent[2]:5.1f}% of tot. var.' 
            else:
                axes_label_name = ['x', 'y', 'z']
                xlabel = axes_label_name[0]
                ylabel = axes_label_name[1]
                zlabel = axes_label_name[2]

        ax.set_xlabel(xlabel, **options_xlabel)
        ax.set_ylabel(ylabel, **options_ylabel)
        ax.set_zlabel(zlabel, **options_zlabel)

    if options_xaxis_ticks is None:
        options_xaxis_ticks = {}

    if options_yaxis_ticks is None:
        options_yaxis_ticks = {}

    if options_zaxis_ticks is None:
        options_zaxis_ticks = {}

    ax.xaxis.set_tick_params(**options_xaxis_ticks)
    ax.yaxis.set_tick_params(**options_yaxis_ticks)
    ax.zaxis.set_tick_params(**options_zaxis_ticks)
    
    if aspect_equal:
        ax.set_aspect('equal')
    elif box_aspect is not None:
        ax.set_box_aspect(box_aspect)
# ----------------------------------------------------------------------------

# ===== Functions for plot in 3D (networkx graph) based on pyvista ===========
# Note: `import pyvista as pv` is set in the functions, so that the karsnet 
# package can still be used if pyvista is not installed (only some specific 
# functions will not be available).
#
# Julien Straubhaar

# ----------------------------------------------------------------------------
def get_mesh_from_graph(G, pos_attr='pos', return_node_label_list=True):
    """
    Get mesh of a graph (for plotting with pyvista).

    The graph `G` in input is assumed to have nodes localized in 2D or 3D space.

    Parameters
    ----------
    G : networkx.Graph object
        graph, with 2D or 3D position as node attribute
    
    pos_attr : str, default: 'pos'
        name of the attribute attached to nodes defining the position as a 
        sequence of 2 or 3 floats; if position are in 2D, a coordinate 
        of zero is set in the 3rd dimension
    
    return_node_label_list : bool, default: True
        if `True`, the list of node label is returned

    Returns
    -------
    mesh : pyvista.PolyData object
        mesh with points (node of the graph `G`) and lines (edges of 
        the graph `G`), that can be plotted with pyvista
    
    node_label_list : list, optional
        returned if `return_node_label_list=True`: list of node labels,
        the node i (integer id) has the label (in input graph `G`) `node_label_list[i]`

    Examples
    --------
    >>> G_mesh, G_node_label_list = kn.view.get_mesh_from_graph(G)
    >>> 
    >>> # Plot (pyvista)
    >>> # ----
    >>> import pyvista as pv
    >>> # pp = pv.Plotter(notebook=False) # pop-up window
    >>> pp = pv.Plotter()
    >>> pp.add_mesh(G_mesh, line_width=3, color='tab:red', opacity=1) # draw edges
    >>> pp.add_mesh(G_mesh.points, render_points_as_spheres=True, point_size=15, color='tab:blue') # draw nodes
    >>> pp.add_point_labels(G_mesh.points, G_node_label_list, point_size=15, font_size=36) # add node labels
    >>> pp.show_grid() # show grid
    >>> # pp.show_bounds() # show bounds
    >>> # pp.screenshot('file.png', transparent_background=False) # save file (png)
    >>> pp.show()
    """
    import pyvista as pv

    # Get dictionary of 3d node position
    pos = kn.utils.get_pos3d(G, pos_attr=pos_attr)

    # Set dictionary to convert node label (id) to node index
    node_label2index = {u:i for i, u in enumerate(G.nodes())}
    
    # Set 2d-array of points: of shape (n_nodes, 3), row i is the position in 3D of the node i (index)
    points = np.asarray(list(pos.values()))
    # if points.shape[1] == 2: # 2D
    #     points = np.insert(points, 2, 0.0, axis=1) # add 0.0 as 3-rd coordinate (z)

    # Set lines: 2d-array of shape (n_edges, 3), row i is [2, j0, j1] where (j0, j1) is an edge btw node j0 and node j1 (integer ids)
    lines = np.asarray([[node_label2index[u], node_label2index[v]] for u, v in G.edges()])
    if len(lines):
        lines = np.insert(lines, 0, 2, axis=1)
        lines = lines.ravel()
    else:
        lines = None

    # Set mesh
    mesh = pv.PolyData(points, lines=lines)

    if return_node_label_list:
        node_label_list = list(node_label2index.keys())
        return mesh, node_label_list

    return mesh
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
def pv_plot_graph(
        G,         
        plotter=None,
        node_attr=None,
        edge_attr=None,
        pos_attr='pos',
        mode='xyz',
        centering=False,
        plot_nodes=True,
        plot_edges=True,
        options_nodes=None,
        options_edges=None, 
        show_node_colorbar=True,
        show_edge_colorbar=True,
        options_node_colorbar=None,
        options_edge_colorbar=None,
        with_labels=False, 
        options_labels=None,
        show_grid=True,
        options_show_grid=None,
        show_outline=False,
        options_show_outline=None,
        cpos=None,
        return_cpos=False):
    """
    Plot a graph in 3D using pyvista.

    The graph `G` in input is assumed to have nodes localized in 2D or 3D space.

    Parameters
    ----------
    G : networkx.Graph object
        graph, with 2D or 3D position as node attribute
    
    plotter : :class:`pyvista.Plotter`, optional
        - if given (not `None`), add element to the plotter, a further call to \
        `plotter.show()` will be required to show the plot
        - if not given (`None`, default): a plotter is created and the plot \
        is shown

    node_attr : str, optional
        node attribute to be plotted; this should be a scalar attribute

    edge_attr : str, optional
        edge attribute to be plotted; this should be a scalar attribute

    pos_attr : str, default: 'pos'
        name of the attribute attached to nodes defining the position as a 
        sequence of 2 or 3 floats; if position are in 2D, a coordinate 
        of zero is set in the 3rd dimension
    
    mode : str {'xyz', 'pca'}, default: 'xyz'
        - if `mode='xyz'`: original 3D axes system is used
        - if `mode='pca'`: the 3D axes system retrieved from the \
        principal component analysis (pca) is used; this axes system is \
        orthonormal and the principal axes are ordered such that the variance \
        of the coordinates of the node position projected onto these axes are \
        decreasing

    centering : bool, default: False
        - if `True`: the coordinates of the position are centered (mean equals to zero)
        - if `False`: the origin of the system is not changed

    plot_nodes : bool, default: True
        - if `True`: the nodes are plotted
        - if `False`: the nodes are not plotted (other parameters about \
        plotting nodes will be ignored)

    plot_edges : bool, default: True
        - if `True`: the edges are plotted
        - if `False`: the edges are not plotted (other parameters about \
        plotting edges will be ignored)
        
    options_nodes : dict, optional
        dictionary of options that is passed to `pyvista.add_mesh` to draw 
        the points, e.g. possible keys: 'point_size', 'opacity', etc.; 

        for size:
        
        - options `point_size` can be set to a float (or int), \
        or to a sequence of n_nodes values (floats or ints), where n_nodes \
        is the number of nodes in the entire graph `G` and where the values \
        are the desired size for nodes given by `G.nodes()`

        for colors:

        - if `node_attr` is not given (`None`): 'color' is used
        - if `node_attr` is given (not `None`): `cmap` is used and \
        and for nan values, `cmap.get_cbad()` or 'color' (if specified, \
        prevails), or 'nan_color' (if specified, prevails) is used for \
        nan values; moreover \
        'vmin' and / or 'vmax' (or 'clim', prevails if specified) may be \
        specified to limit the range of colorbar, in this case \
        `cmap.get_under()` (or 'below_color', prevails if specified)) resp. \
        `cmap.get_over()` (or 'above_color', prevails if specified)) \
        are used for values out of the range

    options_edges : dict, optional
        dictionary of options that is passed to `pyvista.add_mesh` to draw 
        the edges (e.g. possible keys: 'line_width', 'opacity', etc.;

        for colors:

        - if `edge_attr` is not given (`None`): 'color' is used
        - if `edge_attr` is given (not `None`): `cmap` is used and \
        and for nan values, `cmap.get_cbad()` or 'color' (if specified, \
        prevails), or 'nan_color' (if specified, prevails) is used for \
        nan values; moreover \
        'vmin' and / or 'vmax' (or 'clim', prevails if specified) may be \
        specified to limit the range of colorbar, in this case \
        `cmap.get_under()` (or 'below_color', prevails if specified)) resp. \
        `cmap.get_over()` (or 'above_color', prevails if specified)) \
        are used for values out of the range

    show_node_colorbar : bool, default: True
        used only when `node_attr` is not None: if `True`, colorbar (scalar bar)
        for node attribute is displayed
        
    show_edge_colorbar : bool, default: True
        used only when `edge_attr` is not None: if `True`, colorbar (scalar bar)
        for edge attribute is displayed

    options_node_colorbar : dict, optional
        options (keyword arguments) passed to the function `pyvista.add_scalar_bar`
        for plotting the colorbar (scalar bar) for the node attribute
        (e.g. possible keys: 'orientation', 'location', 'shrink', etc.);
        used if `node_attr` is not `None` and `show_node_colorbar=True`

    options_edge_colorbar : dict, optional
        options (keyword arguments) passed to the function `pyvista.add_scalar_bar`
        for plotting the colorbar (scalar bar) for the edge attribute
        (e.g. possible keys: 'orientation', 'location', 'shrink', etc.);
        used if `edge_attr` is not `None` and `show_edge_colorbar=True`

    with_labels : bool, default: False
        if `True`, add node labels

    options_labels : dict, optional
        dictionary of options that is passed to `pyvista.add_point_labels` to 
        draw the nodes labels (e.g. possible key: 'font_size');
        ignored if `with_labels=False`
            
    show_grid : bool, default: True
        indicates if the grid with axe labels are shown
        
    options_show_grid : dict, optional
        options (keyword arguments) passed to the function `pyvista.show_grid`
        (e.g. possible keys: 'xtitle', 'ytitle', 'ztitle' (for x, y, z labels), etc.)

    show_outline : bool, default: False
        indicates if the outline (box) is displayed
        
    options_show_outline : dict, optional
        options (keyword arguments) passed to the function `pyvista.add_mesh`
        (e.g. possible keys: 'line_width', 'color', etc.)

    cpos : sequence[sequence[float]], optional
        camera position (unused if `plotter=None`);
        `cpos` = [camera_location, focus_point, viewup_vector], with

        - camera_location: (tuple of length 3) camera location ("eye")
        - focus_point    : (tuple of length 3) focus point
        - viewup_vector  : (tuple of length 3) viewup vector (vector \
        attached to the "head" and pointed to the "sky")

        note: in principle, (focus_point - camera_location) is orthogonal to
        viewup_vector

    return_pos : bool, default: False
        if `True`, the camera position `cpos` is returned (useful in interactive plot
        for retrieving the camera position which can be used for further plots)
    
    Returns
    -------
    cpos : sequence
        returned if `return_cpos=True`

    Examples
    --------
        >>> kn.view.pv_plot_graph(G)
    """
    import pyvista as pv

    if mode == 'xyz':
        pca_mode = False
    elif mode == 'pca':
        pca_mode = True
    else:
        raise ValueError('Parameter `mode` not valid')

    mesh, node_label_list = get_mesh_from_graph(G, pos_attr=pos_attr, return_node_label_list=True)

    if len(mesh.points) == 0:
        print('No node in graph!')
        return

    if centering:
        # Substact the mean position
        mesh.points = mesh.points - mesh.points.mean(axis=0)

    if pca_mode:
        # Do the pca, get the pca axes and pca variances
        pca_axes, pca_var = kn.utils.pca_of_node_position(G)
        # Set node position in the pca axes system
        mesh.points = mesh.points.dot(pca_axes)

    if options_edges is None:
        options_edges = {}
    
    if options_nodes is None:
        options_nodes = {}
    
    if options_labels is None:
        options_labels = {}

    if plotter is not None:
        pp = plotter
    else:
        pp = pv.Plotter()

            # Set (update) s (node size) option
    if 'point_size' in options_nodes.keys() and hasattr(options_nodes['point_size'], '__len__'):
        point_size = options_nodes['point_size']
        del(options_nodes['point_size'])
        scale_point = True
    else:
        scale_point = False

    if node_attr is not None or scale_point:
        meshp = pv.PolyData(mesh.points)
    else:
        meshp = mesh.points

    if scale_point:
        meshp['__point_size__'] = point_size

    if node_attr is not None:
        meshp[node_attr] = np.asarray(list(nx.get_node_attributes(G, node_attr, default=np.nan).values()))
        meshp.set_active_scalars(node_attr)

    if edge_attr is not None:
        scalars = np.asarray(list(nx.get_edge_attributes(G, edge_attr).values()))
    else:
        scalars = None

    if plot_nodes:
        # Draw nodes
        if node_attr is not None:
            # Get bad color from cmap
            if 'cmap' in options_nodes.keys():
                cmap = plt.get_cmap(options_nodes['cmap'])
            else:
                cmap = plt.get_cmap('viridis') # default cmap
            
            cbad = cmap.get_bad()

            if 'color' in options_nodes.keys():
                # cbad = matplotlib.colors.to_rgba(options_nodes['color'])
                cbad = matplotlib.colors.to_hex(options_nodes['color'])
                del(options_nodes['color'])

            if 'nan_color' not in options_nodes.keys():
                options_nodes['nan_color'] = cbad

            # Check if limit (clim, vmin, vmax are given), set color for below (under) and above (over), 
            # and plot accordingly
            below_color = None
            above_color = None

            if 'clim' in options_nodes.keys():
                vmin, vmax = clim
                below_color = cmap.get_under()
                above_color = cmap.get_over()

            elif 'vmin' or 'vmax' in options_nodes.keys():
                if 'vmin' in options_nodes.keys():
                    vmin = options_nodes['vmin']
                    below_color = cmap.get_under()
                    del(options_nodes['vmin'])
                else:
                    vmin = np.nanmin(meshp[node_attr])

                if 'vmax' in options_nodes.keys():
                    vmax = options_nodes['vmax']
                    above_color = cmap.get_over()
                    del(options_nodes['vmax'])
                else:
                    vmax = np.nanmax(meshp[node_attr])

                clim = [vmin, vmax]
                options_nodes['clim'] = clim

            if below_color is not None and 'below_color' not in options_nodes.keys():
                # below_color = matplotlib.colors.to_rgba(below_color)
                below_color = matplotlib.colors.to_hex(below_color)
                options_nodes['below_color'] = below_color

            if above_color is not None and 'above_color' not in options_nodes.keys():
                # above_color = matplotlib.colors.to_rgba(above_color)
                above_color = matplotlib.colors.to_hex(above_color)
                options_nodes['above_color'] = above_color

            if scale_point:
                actor = pp.add_mesh(meshp, style='points_gaussian', emissive=False, render_points_as_spheres=True, **options_nodes, show_scalar_bar=False)
                actor.mapper.scale_array = '__point_size__'
            else:
                pp.add_mesh(meshp, render_points_as_spheres=True, **options_nodes, show_scalar_bar=False)

            if show_node_colorbar:
                if options_node_colorbar is None:
                    options_node_colorbar = {}

                if 'title' not in options_node_colorbar.keys():
                    options_node_colorbar['title'] = f'node: {node_attr}'
        
                pp.add_scalar_bar(**options_node_colorbar)

        else:
            if scale_point:
                actor = pp.add_mesh(meshp, style='points_gaussian', emissive=False, render_points_as_spheres=True, **options_nodes, show_scalar_bar=False)
                actor.mapper.scale_array = '__point_size__'
            else:
                pp.add_mesh(meshp, render_points_as_spheres=True, **options_nodes, show_scalar_bar=False)

        # Draw node labels (if needed)
        if with_labels:
            pp.add_point_labels(meshp, node_label_list, **options_labels)

    if plot_edges:
        # Draw edges
        if edge_attr is not None:
            # Get bad color from cmap
            if 'cmap' in options_edges.keys():
                cmap = plt.get_cmap(options_edges['cmap'])
            else:
                cmap = plt.get_cmap('viridis') # default cmap
            
            cbad = cmap.get_bad()

            if 'color' in options_edges.keys():
                # cbad = matplotlib.colors.to_rgba(options_edges['color'])
                cbad = matplotlib.colors.to_hex(options_edges['color'])
                del(options_edges['color'])

            if 'nan_color' not in options_edges.keys():
                options_edges['nan_color'] = cbad

            # Check if limit (clim, vmin, vmax are given), set color for below (under) and above (over), 
            # and plot accordingly
            below_color = None
            above_color = None

            if 'clim' in options_edges.keys():
                vmin, vmax = clim
                below_color = cmap.get_under()
                above_color = cmap.get_over()

            elif 'vmin' or 'vmax' in options_edges.keys():
                if 'vmin' in options_edges.keys():
                    vmin = options_edges['vmin']
                    below_color = cmap.get_under()
                    del(options_edges['vmin'])
                else:
                    vmin = np.nanmin(scalars)

                if 'vmax' in options_edges.keys():
                    vmax = options_edges['vmax']
                    above_color = cmap.get_over()
                    del(options_edges['vmax'])
                else:
                    vmax = np.nanmax(scalars)

                clim = [vmin, vmax]
                options_edges['clim'] = clim

            if below_color is not None and 'below_color' not in options_edges.keys():
                # below_color = matplotlib.colors.to_rgba(below_color)
                below_color = matplotlib.colors.to_hex(below_color)
                options_edges['below_color'] = below_color

            if above_color is not None and 'above_color' not in options_edges.keys():
                # above_color = matplotlib.colors.to_rgba(above_color)
                above_color = matplotlib.colors.to_hex(above_color)
                options_edges['above_color'] = above_color

            pp.add_mesh(mesh, scalars=scalars, **options_edges, show_scalar_bar=False)

            if show_edge_colorbar:
                if options_edge_colorbar is None:
                    options_edge_colorbar = {}

                if 'title' not in options_edge_colorbar.keys():
                    options_edge_colorbar['title'] = f'edge: {edge_attr}'
        
                pp.add_scalar_bar(**options_edge_colorbar)

        else:
            pp.add_mesh(mesh, scalars=scalars, **options_edges, show_scalar_bar=False)

    if show_grid:
        if options_show_grid is None:
            options_show_grid = {}
    
        if 'xtitle' not in options_show_grid.keys() or 'ytitle' not in options_show_grid.keys() or 'ztitle' not in options_show_grid.keys():
            if pca_mode:
                axes_label_name = ['pca-1st', 'pca-2nd', 'pca-3rd']
                pca_var_percent = pca_var/pca_var.sum() * 100.0
                # xlabel = axes_label_name[0] + f', {pca_var_percent[0]:5.1f}% of tot. var.' # + ', vec = (' + ', '.join([f'{v:.2g}' for v in pca_axes[:, 0]]) + ')' 
                # ylabel = axes_label_name[1] + f', {pca_var_percent[1]:5.1f}% of tot. var.' # + ', vec = (' + ', '.join([f'{v:.2g}' for v in pca_axes[:, 1]]) + ')' 
                # zlabel = axes_label_name[2] + f', {pca_var_percent[2]:5.1f}% of tot. var.' # + ', vec = (' + ', '.join([f'{v:.2g}' for v in pca_axes[:, 2]]) + ')' 
                xlabel = axes_label_name[0] + ': (' + ', '.join([f'{v:.2g}' for v in pca_axes[:, 0]]) + ')' + f'\n{pca_var_percent[0]:5.1f}% of tot. var.' 
                ylabel = axes_label_name[1] + ': (' + ', '.join([f'{v:.2g}' for v in pca_axes[:, 1]]) + ')' + f'\n{pca_var_percent[1]:5.1f}% of tot. var.' 
                zlabel = axes_label_name[2] + ': (' + ', '.join([f'{v:.2g}' for v in pca_axes[:, 2]]) + ')' + f'\n{pca_var_percent[2]:5.1f}% of tot. var.'
            else:
                axes_label_name = ['x', 'y', 'z']
                xlabel = axes_label_name[0]
                ylabel = axes_label_name[1]
                zlabel = axes_label_name[2]

            options_show_grid['xtitle'] = xlabel
            options_show_grid['ytitle'] = ylabel
            options_show_grid['ztitle'] = zlabel
        
        pp.show_grid(**options_show_grid)

    if show_outline:
        if options_show_outline is None:
            options_show_outline = {}

        pp.add_mesh(mesh.outline(), **options_show_outline)

    if plotter is None:
        cpos = pp.show(cpos=cpos, return_cpos=return_cpos)
    else:
        cpos = None

    if return_cpos:
        return cpos
# ----------------------------------------------------------------------------

# # ----------------------------------------------------------------------------
# def pv_plot_graph(
#         G,         
#         plotter=None,
#         node_attr=None,
#         edge_attr=None,
#         pos_attr='pos',
#         mode='xyz',
#         centering=False,
#         plot_nodes=True,
#         plot_edges=True,
#         options_nodes=None,
#         options_edges=None, 
#         show_node_colorbar=True,
#         show_edge_colorbar=True,
#         options_node_colorbar=None,
#         options_edge_colorbar=None,
#         with_labels=False, 
#         options_labels=None,
#         show_grid=True,
#         options_show_grid=None,
#         show_outline=False,
#         options_show_outline=None,
#         cpos=None,
#         return_cpos=False):
#     """
#     Plot a graph in 3D using pyvista.

#     The graph `G` in input is assumed to have nodes localized in 2D or 3D space.

#     Parameters
#     ----------
#     G : networkx.Graph object
#         graph, with 2D or 3D position as node attribute
    
#     plotter : :class:`pyvista.Plotter`, optional
#         - if given (not `None`), add element to the plotter, a further call to \
#         `plotter.show()` will be required to show the plot
#         - if not given (`None`, default): a plotter is created and the plot \
#         is shown

#     node_attr : str, optional
#         node attribute to be plotted; this should be a scalar attribute

#     edge_attr : str, optional
#         edge attribute to be plotted; this should be a scalar attribute

#     pos_attr : str, default: 'pos'
#         name of the attribute attached to nodes defining the position as a 
#         sequence of 2 or 3 floats; if position are in 2D, a coordinate 
#         of zero is set in the 3rd dimension
    
#     mode : str {'xyz', 'pca'}, default: 'xyz'
#         - if `mode='xyz'`: original 3D axes system is used
#         - if `mode='pca'`: the 3D axes system retrieved from the \
#         principal component analysis (pca) is used; this axes system is \
#         orthonormal and the principal axes are ordered such that the variance \
#         of the coordinates of the node position projected onto these axes are \
#         decreasing

#     centering : bool, default: False
#         - if `True`: the coordinates of the position are centered (mean equals to zero)
#         - if `False`: the origin of the system is not changed

#     plot_nodes : bool, default: True
#         - if `True`: the nodes are plotted
#         - if `False`: the nodes are not plotted (other parameters about \
#         plotting nodes will be ignored)

#     plot_edges : bool, default: True
#         - if `True`: the edges are plotted
#         - if `False`: the edges are not plotted (other parameters about \
#         plotting edges will be ignored)
        
#     options_nodes : dict, optional
#         dictionary of options that is passed to `pyvista.add_mesh` to draw 
#         the points, e.g. possible keys: 'point_size', 'opacity', etc.; 

#         for colors:

#         - if `node_attr` is not given (`None`): 'color' is used
#         - if `node_attr` is given (not `None`): `cmap` is used and \
#         and for nan values, `cmap.get_cbad()` or 'color' (if specified, \
#         prevails), or 'nan_color' (if specified, prevails) is used for \
#         nan values; moreover \
#         'vmin' and / or 'vmax' (or 'clim', prevails if specified) may be \
#         specified to limit the range of colorbar, in this case \
#         `cmap.get_under()` (or 'below_color', prevails if specified)) resp. \
#         `cmap.get_over()` (or 'above_color', prevails if specified)) \
#         are used for values out of the range

#     options_edges : dict, optional
#         dictionary of options that is passed to `pyvista.add_mesh` to draw 
#         the edges (e.g. possible keys: 'line_width', 'opacity', etc.;

#         for colors:

#         - if `edge_attr` is not given (`None`): 'color' is used
#         - if `edge_attr` is given (not `None`): `cmap` is used and \
#         and for nan values, `cmap.get_cbad()` or 'color' (if specified, \
#         prevails), or 'nan_color' (if specified, prevails) is used for \
#         nan values; moreover \
#         'vmin' and / or 'vmax' (or 'clim', prevails if specified) may be \
#         specified to limit the range of colorbar, in this case \
#         `cmap.get_under()` (or 'below_color', prevails if specified)) resp. \
#         `cmap.get_over()` (or 'above_color', prevails if specified)) \
#         are used for values out of the range

#     show_node_colorbar : bool, default: True
#         used only when `node_attr` is not None: if `True`, colorbar (scalar bar)
#         for node attribute is displayed
        
#     show_edge_colorbar : bool, default: True
#         used only when `edge_attr` is not None: if `True`, colorbar (scalar bar)
#         for edge attribute is displayed

#     options_node_colorbar : dict, optional
#         options (keyword arguments) passed to the function `pyvista.add_scalar_bar`
#         for plotting the colorbar (scalar bar) for the node attribute
#         (e.g. possible keys: 'orientation', 'location', 'shrink', etc.);
#         used if `node_attr` is not `None` and `show_node_colorbar=True`

#     options_edge_colorbar : dict, optional
#         options (keyword arguments) passed to the function `pyvista.add_scalar_bar`
#         for plotting the colorbar (scalar bar) for the edge attribute
#         (e.g. possible keys: 'orientation', 'location', 'shrink', etc.);
#         used if `edge_attr` is not `None` and `show_edge_colorbar=True`

#     with_labels : bool, default: False
#         if `True`, add node labels

#     options_labels : dict, optional
#         dictionary of options that is passed to `pyvista.add_point_labels` to 
#         draw the nodes labels (e.g. possible key: 'font_size');
#         ignored if `with_labels=False`
            
#     show_grid : bool, default: True
#         indicates if the grid with axe labels are shown
        
#     options_show_grid : dict, optional
#         options (keyword arguments) passed to the function `pyvista.show_grid`
#         (e.g. possible keys: 'xtitle', 'ytitle', 'ztitle' (for x, y, z labels), etc.)

#     show_outline : bool, default: False
#         indicates if the outline (box) is displayed
        
#     options_show_outline : dict, optional
#         options (keyword arguments) passed to the function `pyvista.add_mesh`
#         (e.g. possible keys: 'line_width', 'color', etc.)

#     cpos : sequence[sequence[float]], optional
#         camera position (unused if `plotter=None`);
#         `cpos` = [camera_location, focus_point, viewup_vector], with

#         - camera_location: (tuple of length 3) camera location ("eye")
#         - focus_point    : (tuple of length 3) focus point
#         - viewup_vector  : (tuple of length 3) viewup vector (vector \
#         attached to the "head" and pointed to the "sky")

#         note: in principle, (focus_point - camera_location) is orthogonal to
#         viewup_vector

#     return_pos : bool, default: False
#         if `True`, the camera position `cpos` is returned (useful in interactive plot
#         for retrieving the camera position which can be used for further plots)
    
#     Returns
#     -------
#     cpos : sequence
#         returned if `return_cpos=True`

#     Examples
#     --------
#         >>> kn.view.pv_plot_graph(G)
#     """
#     import pyvista as pv

#     if mode == 'xyz':
#         pca_mode = False
#     elif mode == 'pca':
#         pca_mode = True
#     else:
#         raise ValueError('Parameter `mode` not valid')

#     mesh, node_label_list = get_mesh_from_graph(G, pos_attr=pos_attr, return_node_label_list=True)

#     if len(mesh.points) == 0:
#         print('No node in graph!')
#         return

#     if centering:
#         # Substact the mean position
#         mesh.points = mesh.points - mesh.points.mean(axis=0)

#     if pca_mode:
#         # Do the pca, get the pca axes and pca variances
#         pca_axes, pca_var = kn.utils.pca_of_node_position(G)
#         # Set node position in the pca axes system
#         mesh.points = mesh.points.dot(pca_axes)

#     if options_edges is None:
#         options_edges = {}
    
#     if options_nodes is None:
#         options_nodes = {}
    
#     if options_labels is None:
#         options_labels = {}

#     if plotter is not None:
#         pp = plotter
#     else:
#         pp = pv.Plotter()

#     if node_attr is not None:
#         meshp = pv.PolyData(mesh.points)
#         meshp[node_attr] = np.asarray(list(nx.get_node_attributes(G, node_attr, default=np.nan).values()))
#         meshp.set_active_scalars(node_attr) # not mandatory
#     else:
#         meshp = mesh.points

#     if edge_attr is not None:
#         scalars = np.asarray(list(nx.get_edge_attributes(G, edge_attr).values()))
#     else:
#         scalars = None

#     if plot_nodes:
#         # Draw nodes
#         if node_attr is not None:
#             # Get bad color from cmap
#             if 'cmap' in options_nodes.keys():
#                 cmap = plt.get_cmap(options_nodes['cmap'])
#             else:
#                 cmap = plt.get_cmap('viridis') # default cmap
            
#             cbad = cmap.get_bad()

#             if 'color' in options_nodes.keys():
#                 cbad = matplotlib.colors.to_hex(options_nodes['color'])
#                 del(options_nodes['color'])

#             if 'nan_color' not in options_nodes.keys():
#                 options_nodes['nan_color'] = cbad

#             # Check if limit (clim, vmin, vmax are given), set color for below (under) and above (over), 
#             # and plot accordingly
#             below_color = None
#             above_color = None

#             if 'clim' in options_nodes.keys():
#                 vmin, vmax = clim
#                 below_color = cmap.get_under()
#                 above_color = cmap.get_over()

#             elif 'vmin' or 'vmax' in options_nodes.keys():
#                 if 'vmin' in options_nodes.keys():
#                     vmin = options_nodes['vmin']
#                     below_color = cmap.get_under()
#                     del(options_nodes['vmin'])
#                 else:
#                     vmin = np.nanmin(meshp[node_attr])

#                 if 'vmax' in options_nodes.keys():
#                     vmax = options_nodes['vmax']
#                     above_color = cmap.get_over()
#                     del(options_nodes['vmax'])
#                 else:
#                     vmax = np.nanmax(meshp[node_attr])

#                 clim = [vmin, vmax]
#                 options_nodes['clim'] = clim

#             if below_color is not None and 'below_color' not in options_nodes.keys():
#                 below_color = matplotlib.colors.to_hex(below_color)
#                 options_nodes['below_color'] = below_color

#             if above_color is not None and 'above_color' not in options_nodes.keys():
#                 above_color = matplotlib.colors.to_hex(above_color)
#                 options_nodes['above_color'] = above_color

#             pp.add_mesh(meshp, render_points_as_spheres=True, **options_nodes, show_scalar_bar=False)

#             if show_node_colorbar:
#                 if options_node_colorbar is None:
#                     options_node_colorbar = {}

#                 if 'title' not in options_node_colorbar.keys():
#                     options_node_colorbar['title'] = f'node: {node_attr}'
        
#                 pp.add_scalar_bar(**options_node_colorbar)

#         else:
#             pp.add_mesh(meshp, render_points_as_spheres=True, **options_nodes, show_scalar_bar=False)

#         # Draw node labels (if needed)
#         if with_labels:
#             pp.add_point_labels(meshp, node_label_list, **options_labels)

#     if plot_edges:
#         # Draw edges
#         if edge_attr is not None:
#             # Get bad color from cmap
#             if 'cmap' in options_edges.keys():
#                 cmap = plt.get_cmap(options_edges['cmap'])
#             else:
#                 cmap = plt.get_cmap('viridis') # default cmap
            
#             cbad = cmap.get_bad()

#             if 'color' in options_edges.keys():
#                 cbad = matplotlib.colors.to_hex(options_edges['color'])
#                 del(options_edges['color'])

#             if 'nan_color' not in options_edges.keys():
#                 options_edges['nan_color'] = cbad

#             # Check if limit (clim, vmin, vmax are given), set color for below (under) and above (over), 
#             # and plot accordingly
#             below_color = None
#             above_color = None

#             if 'clim' in options_edges.keys():
#                 vmin, vmax = clim
#                 below_color = cmap.get_under()
#                 above_color = cmap.get_over()

#             elif 'vmin' or 'vmax' in options_edges.keys():
#                 if 'vmin' in options_edges.keys():
#                     vmin = options_edges['vmin']
#                     below_color = cmap.get_under()
#                     del(options_edges['vmin'])
#                 else:
#                     vmin = np.nanmin(scalars)

#                 if 'vmax' in options_edges.keys():
#                     vmax = options_edges['vmax']
#                     above_color = cmap.get_over()
#                     del(options_edges['vmax'])
#                 else:
#                     vmax = np.nanmax(scalars)

#                 clim = [vmin, vmax]
#                 options_edges['clim'] = clim

#             if below_color is not None and 'below_color' not in options_edges.keys():
#                 below_color = matplotlib.colors.to_hex(below_color)
#                 options_edges['below_color'] = below_color

#             if above_color is not None and 'above_color' not in options_edges.keys():
#                 above_color = matplotlib.colors.to_hex(above_color)
#                 options_edges['above_color'] = above_color

#             pp.add_mesh(mesh, scalars=scalars, **options_edges, show_scalar_bar=False)

#             if show_edge_colorbar:
#                 if options_edge_colorbar is None:
#                     options_edge_colorbar = {}

#                 if 'title' not in options_edge_colorbar.keys():
#                     options_edge_colorbar['title'] = f'edge: {edge_attr}'
        
#                 pp.add_scalar_bar(**options_edge_colorbar)

#         else:
#             pp.add_mesh(mesh, scalars=scalars, **options_edges, show_scalar_bar=False)

#     if show_grid:
#         if options_show_grid is None:
#             options_show_grid = {}
    
#         if 'xtitle' not in options_show_grid.keys() or 'ytitle' not in options_show_grid.keys() or 'ztitle' not in options_show_grid.keys():
#             if pca_mode:
#                 axes_label_name = ['pca-1st', 'pca-2nd', 'pca-3rd']
#                 pca_var_percent = pca_var/pca_var.sum() * 100.0
#                 # xlabel = axes_label_name[0] + f', {pca_var_percent[0]:5.1f}% of tot. var.' # + ', vec = (' + ', '.join([f'{v:.2g}' for v in pca_axes[:, 0]]) + ')' 
#                 # ylabel = axes_label_name[1] + f', {pca_var_percent[1]:5.1f}% of tot. var.' # + ', vec = (' + ', '.join([f'{v:.2g}' for v in pca_axes[:, 1]]) + ')' 
#                 # zlabel = axes_label_name[2] + f', {pca_var_percent[2]:5.1f}% of tot. var.' # + ', vec = (' + ', '.join([f'{v:.2g}' for v in pca_axes[:, 2]]) + ')' 
#                 xlabel = axes_label_name[0] + ': (' + ', '.join([f'{v:.2g}' for v in pca_axes[:, 0]]) + ')' + f'\n{pca_var_percent[0]:5.1f}% of tot. var.' 
#                 ylabel = axes_label_name[1] + ': (' + ', '.join([f'{v:.2g}' for v in pca_axes[:, 1]]) + ')' + f'\n{pca_var_percent[1]:5.1f}% of tot. var.' 
#                 zlabel = axes_label_name[2] + ': (' + ', '.join([f'{v:.2g}' for v in pca_axes[:, 2]]) + ')' + f'\n{pca_var_percent[2]:5.1f}% of tot. var.'
#             else:
#                 axes_label_name = ['x', 'y', 'z']
#                 xlabel = axes_label_name[0]
#                 ylabel = axes_label_name[1]
#                 zlabel = axes_label_name[2]

#             options_show_grid['xtitle'] = xlabel
#             options_show_grid['ytitle'] = ylabel
#             options_show_grid['ztitle'] = zlabel
        
#         pp.show_grid(**options_show_grid)

#     if show_outline:
#         if options_show_outline is None:
#             options_show_outline = {}

#         pp.add_mesh(mesh.outline(), **options_show_outline)

#     if plotter is None:
#         cpos = pp.show(cpos=cpos, return_cpos=return_cpos)
#     else:
#         cpos = None

#     if return_cpos:
#         return cpos
# # ----------------------------------------------------------------------------
