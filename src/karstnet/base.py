#    Copyright (C) 2018-2025 by
#    Philippe Renard <philippe.renard@unine.ch>
#    Pauline Collon <pauline.collon@univ-lorraine.fr>
#    All rights reserved.
#    MIT license.
#
"""
The module `base` contains the class `KGraph` for karst network object.
"""

# External libraries importations
import numpy as np
import networkx as nx
import scipy.stats as st
import matplotlib.pyplot as plt
import mplstereonet # for projection='stereo'
# noinspection PyUnresolvedReferences

# Internal libraries importations
import karstnet as kn #from karstnet import utils


# *****************************************************************************
# Test function
# *****************************************************************************
def test_kn():
    print("test ok")
    print("relance ok")


# *****************************************************************************
# KGraph class
# *****************************************************************************
class KGraph:
    """
    Class dedicated to the construction and manipulation of graphs
    representing Karstic network.

    Attributes
    ----------
    - graph : networkx.Graph
        the original graph of the karstic network;
        each station is a node, each line-of-sight is an edge;
        `graph` is a networkx.Graph object, with :

        - node attributes :
            - 'pos' : sequence of 2 or 3 floats
                position in 2D or 3D
            
            - 'branch_ids_list' : list
                list of ids of branches containing the node
                (or, equivalently, list of node labels of graph_of_branches,
                corresponding to branches containing the node)

            - other optional properties; this allows to
                store additional information

        - edge attributes :
            - 'length' : float
                real length of the edge

            - 'length2d': float
                length of the projction of the edge in the xy plane

            - 'azimuth' : float
                azimuth in degree in [0, 360), from North towards East
                (counter-clockwise); it provide "bearings" for stereonet;
                for vertical edge, the azimuth is set to Nan

            - 'dip' : float
                dip in degree in [0, 90], "plunging" direction, towards
                under the xy plane; it provide "plunges" for stereonet

            - other optional properties; this allows to
                store additional information

    - graph_simpl : networkx.Graph
        the simplified graph of the karstic network,
        i.e. all nodes of degree 2 are removed (from the original graph)
        except in loops (2 types of loops, cf Collon et al., 2017,
        Geomorphology), so that the topology is preserved;
        `graph_simple` is a networkx.Graph object, with :

        - node attributes :
            - 'pos' : sequence of 2 or 3 floats
                position in 2D or 3D (of retained nodes)

        - edge attributes :
            - 'length' : float
                real length of the corresponding branch in the
                original graph

    - graph_of_branches : networkx.Graph
        the graph of branches of the karstic network, where a node represents 
        a branch in the original graph (labeled by integers from 0) and an edge
        exists between two nodes if the two branches represented by the two nodes 
        have a common extremity, with

        - node attributes :
            - 'pos' : sequence of 2 or 3 floats
                mean position of the original graph nodes in the branch
                represented by the node

            - 'length' : flaot
                real length of the branch represented by the node
        
    - branches : list
        list of branches, each entry is a list of nodes forming a
        path in the original graph corresponding to a branch;
        note: `branches[i]` is the branch corresponding to the node `i`
        in `graph_of_branches`

    - br_lengths : list of floats
        list of branch lengths, each entry is the real length of
        the path forming the corresponding branch in `branches`

    - br_tort : list of floats
        list of branch tortuosities, i.e. real length of the branches
        divided by the distance between the two extremities

    - list_simpl_edges : list of tuple
        list of simple edges, i.e. the edges in the simplified graph,
        necessary to export graph to plines 
    
    Methods
    -------
    See below.
    """
    def __init__(self, edges, coordinates,
                 node_property=None,
                 node_property_name=None,
                 edge_property=None,
                 edge_property_name=None,
                 networkx_graph=None,
                 verbose=True):
        """
        Creates a Kgraph from nodes and edges.

        Parameters
        ----------
        edges : list
            a list of edges

        coordinates : dict
            coordinates of the nodes, dictionary where keys are node names, and
            values are sequences of length 2 or 3 (dimension)

        node_property : dict or list of dicts, optional
            additional node property(ies), keys are node names, values are the 
            property values

        node_property_name : str or list of strs, optional
            name of the additional node property(ies)
            if `node_property` is a dictionary, `node_property_name` must be a string;
            if `node_property` is a list of dictionaries, `node_property_name` must
            be a list of strings

        edge_property : dict or list of dicts, optional
            additional edge property(ies), keys are edge names, values are the 
            property values

        edge_property_name : str or list of strs, optional
            name of the additional edge property(ies)
            if `edge_property` is a dictionary, `edge_property_name` must be a string;
            if `edge_property` is a list of dictionaries, `edge_property_name` must
            be a list of strings

        networkx_graph : networkx.Graph, optional
            if given, it is assumed that all information is given in `networkx_graph`,
            and the attribute `graph` is set to a copy of `networkx_graph`; in this
            case, all parameters above are ignored

        verbose : bool, default: True
            indicates if some information is printed (`True`) during the creation
            of the object

        Examples
        --------
            >>> myKGraph = KGraph([],{}) # create an empty karstic network
        """

        self.verbose = verbose

        if networkx_graph is not None:
            # This assumes that a netorwkx graph is already defined, with
            # node attribute 'pos' for the position, and optional additional
            # node attributes and optional additional edge attributes
            self.graph = networkx_graph.copy()
        
        else:
            # Initialization of the graph (networkx)
            self.graph = nx.Graph()

            # Set nodes and their position
            self.graph.add_nodes_from(list(coordinates.keys()))
            nx.set_node_attributes(self.graph, coordinates, 'pos')

            # Set edges
            self.graph.add_edges_from(edges)

            # Set additional node properties
            if node_property is not None:
                if isinstance(node_property, list):
                    # A list of properties (dictionaries) is given
                    if not isinstance(node_property_name, list):
                        print(f'ERROR: `node_property_name` must be a list when `node_property` is a list of dict')
                        return
                    for node_prop, node_prop_name in zip(node_property, node_property_name):
                        nx.set_node_attributes(self.graph, node_prop, node_prop_name)
                else:
                    # One property (dictionary) is given
                    nx.set_node_attributes(self.graph, node_property, node_property_name)

            # Set additional edge properties
            if edge_property is not None:
                if isinstance(edge_property, list):
                    # A list of properties (dictionaries) is given
                    if not isinstance(edge_property_name, list):
                        print(f'ERROR: `edge_property_name` must be a list when `edge_property` is a list of dict')
                        return
                    for edge_prop, edge_prop_name in zip(edge_property, edge_property_name):
                        nx.set_edge_attributes(self.graph, edge_prop, edge_prop_name)
                else:
                    # One property (dictionary) is given
                    nx.set_edge_attributes(self.graph, edge_property, edge_property_name)

        # ----- Notes -----
        # Position in 2d and in 3d, can be retrieved by
        # using the methods pos2d and pos3d respectively.
        # # self.pos2d, self.pos3d = _pos_initialization(coordinates) # OLD
        # -----

        if self.verbose:
            print(
                "\nThis network contains",
                nx.number_connected_components(self.graph),
                "connected components")

        # Compute and set edge attributes
        self._set_graph_lengths()
        self._set_graph_orientations()

        # Compute the branches of the graph
        # self.branches is necessary to export graph to plines
        self.branches, self.br_lengths, self.br_tort = self._getallbranches()

        # Construct the simplified graph
        # self.list_simpl_edges is necessary to export graph to plines
        self.list_simpl_edges, self.graph_simpl = self._simplify_graph()

        # Construct the graph of branches
        self.graph_of_branches = self._get_graph_of_branches()

        # Set the node attribute 'branch_ids_list' on original graph
        branch_ids_list = {u:[] for u in self.graph.nodes()}
        for i, node_list in enumerate(self.branches):
            for ui in node_list:
                branch_ids_list[ui].append(i)
            
        nx.set_node_attributes(self.graph, branch_ids_list, 'branch_ids_list')

    # *************************************************************************
    # Methods for getting node position
    # *************************************************************************

    # -------------------------------------------------------------------------
    def pos2d(self):
        """
        Gets node position in 2D (ignoring z-coordinate if it exists).
        
        The original graph is considered.

        Returns
        -------
        pos2d: dict
            keys are node names, values are the 2D position of the nodes
        """
        pos2d = kn.utils.get_pos2d(self.graph, pos_attr='pos')
        return pos2d
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    def pos3d(self):
        """
        Gets node position in 3D (setting z-coordinate to zero if it does not exist).

        The original graph is considered.

        Returns
        -------
        pos3d: dict
            keys are node names, values are the 3D position of the nodes
        """
        pos3d = kn.utils.get_pos3d(self.graph, pos_attr='pos')
        return pos3d
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    def posz(self):
        """
        Gets z coordinates of node position of a networkx graph.

        The original graph is considered.

        Returns
        -------
        posz: dict
            keys are node names, values are the z position of the nodes
        """
        posz = kn.utils.get_posz(self.graph, pos_attr='pos')
        return posz
    # -------------------------------------------------------------------------

    # *************************************************************************
    # Methods for plots
    # *************************************************************************

    # -------------------------------------------------------------------------
    def plot2d(self, graph_type=0, **kwargs):
        """
        Plots a 2D view of the karstic network.

        See also function `kn.view.plot_graph2d`.

        Parameters
        ----------
        graph_type : int or str {'original', 'simplified'}
            - if 0 or 'original', displays the original graph,
            - if 1 or 'simplified' displays the simplified graph

        kwargs : dict, optional
            keyword arguments passed to the function `kn.view.plot_graph2d`

        Examples
        --------
            >>> myKGraph = KGraph([],{})
            >>> myKGraph.plot2d()
            >>> myKGraph.plot2d(1)
        """
        if isinstance(graph_type, str):
            if graph_type == 'original':
                graph_type = 0
            elif graph_type == 'simplified':
                graph_type = 1
            else:
                raise ValueError('Parameter `graph_type` not valid')
        
        if graph_type == 0:
            kn.view.plot_graph2d(self.graph, **kwargs)
            plt.title('original')
        elif graph_type == 1:
            kn.view.plot_graph2d(self.graph_simpl, **kwargs)
            plt.title('simplified')
        else:
            raise ValueError('Parameter `graph_type` not valid')
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    def plot3d(self, graph_type=0, **kwargs):
        """
        Plots a 3D view of the karstic network.

        See also function `kn.view.plot_graph3d`.

        Parameters
        ----------
        graph_type : int or str {'original', 'simplified'}
            - if 0 or 'original', displays the original graph,
            - if 1 or 'simplified' displays the simplified graph

        kwargs : dict, optional
            keyword arguments passed to the function `kn.view.plot_graph3d`
            
        Examples
        --------
            >>> myKGraph = KGraph([],{})
            >>> myKGraph.plot3d()
            >>> myKGraph.plot3d(1)
        """
        # 2D  plot
        if isinstance(graph_type, str):
            if graph_type == 'original':
                graph_type = 0
            elif graph_type == 'simplified':
                graph_type = 1
            else:
                raise ValueError('Parameter `graph_type` not valid')
        
        if graph_type == 0:
            kn.view.plot_graph3d(self.graph, **kwargs)
            plt.title('original')
        elif graph_type == 1:
            kn.view.plot_graph3d(self.graph_simpl, **kwargs)
            plt.title('simplified')
        else:
            raise ValueError('Parameter `graph_type` not valid')
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    def pv_plot(self, graph_type=0, **kwargs):
        """
        Plots a 3D view of the karstic network, using `pyvista`.

        See also function `kn.view.pv_plot_graph`.

        Parameters
        ----------
        graph_type : int or str {'original', 'simplified'}
            - if 0 or 'original', displays the original graph,
            - if 1 or 'simplified' displays the simplified graph

        kwargs : dict, optional
            keyword arguments passed to the function `kn.view.pv_plot_graph`
            
        Examples
        --------
            >>> myKGraph = KGraph([],{})
            >>> myKGraph.plot3d()
            >>> myKGraph.plot3d(1)
        """
        # 2D  plot
        if isinstance(graph_type, str):
            if graph_type == 'original':
                graph_type = 0
            elif graph_type == 'simplified':
                graph_type = 1
            else:
                raise ValueError('Parameter `graph_type` not valid')
        
        if graph_type == 0:
            kn.view.pv_plot_graph(self.graph, **kwargs)
        elif graph_type == 1:
            kn.view.pv_plot_graph(self.graph_simpl, **kwargs)
        else:
            raise ValueError('Parameter `graph_type` not valid')
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    def plot2(self,
              graph_type=0,
              figsize=(6, 3),
              return_fig=False,
              **kwargs):
        """
        Plots a 2D view of the karstic network.

        A new figure is created, and can be returned; the method
        `plot2d` is used for the plot.

        Parameters
        ----------
        graph_type : int or str {'original', 'simplified'}
            - if 0 or 'original', displays the original graph,
            - if 1 or 'simplified' displays the simplified graph

        figsize : tuple, default: (6, 3)
            contains the (width, height) dimension of the figure

        return_fig : bool, default: False
            indicates if the figure object is returned

        kwargs : dict, optional
            keyword arguments passed to the method `plot2d`

        Notes
        -----
        Obsolete method, used `plot2d` instead'.

        Returns
        -------
        fig: Figure object (matplotlib), optional        
            figure, returned if `return_fig=True`

        Examples
        --------
            >>> myKGraph = KGraph([],{})
            >>> myKGraph.plot2()
            >>> myKGraph.plot2(1)
        """
        print('Note: method `plot2` is obsolete, use method `plot2d` instead')

        fig = plt.figure(figsize=figsize)

        self.plot2d(graph_type=graph_type, **kwargs)

        plt.show()

        if return_fig:
            return fig
        return
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    def plot3(self,
              graph_type=0,
              figsize=(6, 3),
              return_fig=False,
              **kwargs):
        """
        Plots a 3D view of the karstic network.

        A new figure is created, and can be returned; the method
        `plot3d` is used for the plot.

        Parameters
        ----------
        graph_type : int or str {'original', 'simplified'}
            - if 0 or 'original', displays the original graph,
            - if 1 or 'simplified' displays the simplified graph

        figsize : tuple, default: (6, 3)
            contains the (width, height) dimension of the figure

        return_fig : bool, default: False
            indicates if the figure object is returned

        kwargs : dict, optional
            keyword arguments passed to the method `plot3d`

        Notes
        -----
        Obsolete method, used `plot3d` instead'.

        Returns
        -------
        fig: Figure object (matplotlib), optional        
            figure, returned if `return_fig=True`

        Examples
        --------
            >>> myKGraph = KGraph([],{})
            >>> myKGraph.plot3()
            >>> myKGraph.plot3(1)
        """
        print('Note: method `plot3` is obsolete, use method `plot3d` instead')

        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(projection='3d')

        self.plot3d(graph_type=graph_type, ax=ax, **kwargs)

        plt.show()

        if return_fig:
            return fig
        return
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    def plot(self, figsize=(12, 5), return_fig=False, **kwargs):
        """
        Simple 2D map of the original and simplified karstic network.

        The two maps are ploted side by side. This function allows
        to check rapidly the data after an import for example.

        A new figure is created with two subplots, and can be returned; 
        the method `plot2d` is used for the two plots.

        Parameters
        ----------
        figsize : tuple, default: (12, 5)
            contains the (width, height) dimension of the figure

        return_fig : bool, default: False
            indicates if the figure object is returned

        kwargs : dict, optional
            keyword arguments passed to the method `plot2d` (for 
            the two plots of original and simplified graphs)

        Returns
        -------
        fig: Figure object (matplotlib), optional        
            figure, returned if `return_fig=True`

        Examples
        --------
            >>> myKGraph.plot()
        """
        fig = plt.figure(figsize=figsize)

        plt.subplot(121)
        self.plot2d(graph_type=0, **kwargs)

        plt.subplot(122)
        self.plot2d(graph_type=1, **kwargs)

        plt.show()

        if return_fig:
            return fig
        return
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    def plotxz(self, figsize=(12, 5), return_fig=False):
        """
        Calls method `plot` with keyword arguments xlabel='x', ylabel='z'.
        """
        print('Note: method `plotxz` is obsolete, use method `plot` instead')

        return self.plot(xlabel='x', ylabel='z', figsize=figsize, return_fig=return_fig)
    # -------------------------------------------------------------------------

    # ------
    # Member function written by Philippe Vernant 2019/11/25
    # Modified by Pauline Collon (aug. 2020) to weight density map by lengths

    # Some additional modif by Julien Straubhaar 2025/01
    
    def stereo(self, weighted=True, figsize=(16, 8), return_fig=False):
        """
        Density map of orientations and rose diagram of the karstic network.

        The two maps are ploted side by side.

        The density of orientations of the edges of the original graph is
        ploted according to Schmidt\'s projection on lower hemisphere (stereo)
        on the left map ("3D orientations"), and as a rose diagram accounting only
        for the azimuth, i.e. orientation of the projection of the edges in the
        xy plane ("2D orientations") on the right map.

        By default, stereo and rose diagram are weighted by the length of the
        edges: for stereo, the length in 3D (real length), and for rose diagram,
        the length in 2D (length of the projection in the xy plane).

        Parameters
        ----------
        weighted : bool, default: True
            indicates if the maps are weighted by length (True), or if
            each edge has the same weight (False)

        figsize : tuple, default: (16, 8)
            contains the (width, height) dimension of the figure

        return_fig : bool, default: False
            indicates if the figure object is returned

        Returns
        -------
        fig: Figure object (matplotlib), optional        
            figure, returned if `return_fig=True`

        Examples
        --------
            >>> myKGraph.stereo()
            >>> myKGraph.stereo(weighted = False)
        """

        # Create an np.array of azimuths and dips
        # and lengths (projected(2d) and real (3d))
        azim = np.asarray(list((nx.get_edge_attributes(self.graph, 'azimuth')).values()))
        azim_not_Nan = azim[~np.isnan(azim)] # for rose diagram (exclude vertical edges)
        bearing_dc = np.nan_to_num(azim) # convert nan to zero (and inf to very large number)
        plunge_dc = np.array(list((nx.get_edge_attributes(self.graph, 'dip')).values()))
        if weighted:
            l2d = np.array(list((nx.get_edge_attributes(self.graph, 'length2d')).values()))
            l2d_not_Nan = l2d[~np.isnan(azim)] # for rose diagram
            l3d = np.array(list((nx.get_edge_attributes(self.graph, 'length')).values()))
        else:
            l2d_not_Nan = None
            l3d = None

        # Making colormap, based on Collon et al.(2017)
        # we saturate the colormap at 40%
        from matplotlib import colormaps
        from matplotlib.colors import ListedColormap
        from matplotlib.gridspec import GridSpec

        nbint = 15
        levels = np.linspace(0, 1, nbint)
        rainbow = colormaps['rainbow']
        newcolors = rainbow(levels)
        white = np.ones(4) # np.array([256 / 256, 256 / 256, 256 / 256, 1])
        newcolors[:1, :] = white
        newcmp = ListedColormap(newcolors)

        # Define the grid for plotting maps and the figure
        gs = GridSpec(nrows=20, ncols=2)
        fig = plt.figure(figsize=figsize)

        # ----- Stereo -----
        # Density map - Allows to consider almost vertical conduits
        # The data are weigthted by the real length of the segments (l3d)
        # Use the traditional "Schmidt" method : 1% count
        dc = fig.add_subplot(gs[:-2, 0], projection='stereonet')
        cdc = dc.density_contourf(plunge_dc,
                                  bearing_dc,
                                  measurement='lines',
                                  method='schmidt',
                                  levels=np.arange(0, nbint * 2 + 1, 2),
                                  extend='both',
                                  cmap=newcmp,
                                  weights=l3d)
        dc.set_title('Density map of orientations [Schmidt\'s projection]',
                     y=1.10,
                     fontsize=15)
        dc.grid()
        dc._polar.set_position(dc.get_position())
        dc.set_azimuth_ticks(np.arange(0, 351, 10))

        dc_cb = fig.add_subplot(gs[-1:, 0]) # axe for the colorbar, made invisible
        for spine in ["top", "bottom", "left", "right"]:
            dc_cb.spines[spine].set_visible(False)
            dc_cb.set_xticks([])
            dc_cb.set_yticks([])

        # Colorbar of the density map
        cbar = plt.colorbar(cdc, ax=dc_cb,
                            fraction=0.95,
                            pad=0.01,
                            orientation='horizontal')
        cbar.set_label('[%]')

        # ----- Rose diagram -----
        # The azimuth data are weighted by the projected length (l2d)
        bin_edges = np.arange(-5, 366, 10)
        number_of_strikes, bin_edges = np.histogram(azim_not_Nan,
                                                    bin_edges,
                                                    weights=l2d_not_Nan)
        number_of_strikes[0] += number_of_strikes[-1]
        half = np.sum(np.split(number_of_strikes[:-1], 2), 0)
        two_halves = np.concatenate([half, half])

        # Rose diagram on right hand side of picture
        rs = fig.add_subplot(gs[:-2,1], projection='polar')
        rs.bar(np.deg2rad(np.arange(0, 360, 10)),
               two_halves,
               width=np.deg2rad(10),
               bottom=0.0,
               color='.8',
               edgecolor='k')
        rs.set_theta_zero_location('N')
        rs.set_theta_direction(-1)
        rs.set_thetagrids(np.arange(0, 360, 10), labels=np.arange(0, 360, 10))
        rs.set_title('Rose Diagram of the cave survey segments',
                     y=1.10,
                     fontsize=15)

        # fig.tight_layout()
        plt.show()

        if return_fig:
            return fig
        return
    # -------------------------------------------------------------------------

    # end modif PV 2019/11/25

    # *************************************************************************
    # Methods for import scalar node or edge properties from csv files
    # *************************************************************************
    # Julien Straubhaar

    # -------------------------------------------------------------------------
    def add_node_scalar_properties_from_csv(
            self, 
            filename,
            properties_name_list=None,
            delimiter=';',
            verbose=True):
        """
        Imports node scalar properties from a csv file and add them to the original graph.

        The file has the columns:
            - ['id', 'x', 'y', 'z',] '<node_prop0>', ['<node_prop1>', ...]
   
        Default id (if column 'id' is not given) is integer starting from 0.
        The columns with a name other than 'id' and 'x', 'y'[, 'z'] (giving the 
        position of the nodes), are additional scalar properties (attributes)
        attached to the nodes.

        The original graph is considered.

        Parameters
        ----------
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

        Examples
        --------
            >>> Kg.add_node_scalar_properties_from_csv("MyKarst_node_properties.csv")
        """
        G = kn.io_func.networkx_node_scalar_properties_from_csv(
            self.graph,
            filename,
            properties_name_list=properties_name_list,
            delimiter=delimiter,
            verbose=False)

        if G is not None:
            # update self.graph
            self.graph = G

            if verbose:
                print(f"Karstnet graph: node scalar properties successfully imported from csv file ({filename}) !\n")

        return
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    def add_edge_scalar_properties_from_csv(
            self, 
            filename,
            properties_name_list=None,
            delimiter=';',
            verbose=True):
        """
        Imports edge scalar properties from a csv file and add them to the original graph.

        The file has the columns:
            - '<from_id>', '<to_id>', '<edge_prop0>'[, '<edge_prop1>', ...]

        where the first two columns, '<from_id>', '<to_id>', are the node ids 
        of the two extremities of an edge (link), the next columns are
        additional scalar properties (attributes) attached to the edges.

        The original graph is considered.

        Parameters
        ----------
        filename : str
            the name of the input file
        
        properties_name_list : list of strs, optional
            list of property names that are exported;
            by default (`None`): all properties are imported
        
        delimiter : str, default: ';'
            delimiter used in the file
        
        verbose : bool, default: True
            indicates if info is printed

        Examples
        --------
            >>> Kg.add_edge_scalar_properties_from_csv("MyKarst_edge_properties.csv")
        """
        G = kn.io_func.networkx_edge_scalar_properties_from_csv(
            self.graph,
            filename,
            properties_name_list=properties_name_list,
            delimiter=delimiter,
            verbose=False)

        if G is not None:
            # update self.graph
            self.graph = G

            if verbose:
                print(f"Karstnet graph: edge scalar properties successfully imported from csv file ({filename}) !\n")

        return
    # -------------------------------------------------------------------------

    # *************************************************************************
    # Methods for export to csv files
    # *************************************************************************
    # Julien Straubhaar

    # -------------------------------------------------------------------------
    def to_csv(
            self,
            basename,
            suffix_nodes='_nodes.csv',
            suffix_edges='_edges.csv',
            delimiter_nodes=';',
            delimiter_edges=';',
            export_node_scalar_properties=True,
            export_edge_scalar_properties=True,
            verbose=True):
        """
        Exports the karstnet graph (original graph) to two csv files (nodes, and edges (links)).

        The file for nodes has the columns:
            - 'id', 'x', 'y'[, 'z', '<node_prop0>', '<node_prop1>', ...]

        where the column 'id' is the id of the nodes, the columns 'x', 'y'[, 'z']
        give the position of the nodes and the next columns are additional
        scalar properties (attributes) attached to the nodes.
        The column 'z' is omitted for 2D case.

        The file for edges (links) has the columns:
            - 'from_id', 'to_id'[, '<edge_prop0>', '<edge_prop1>', ...]

        where the first two columns, 'from_id', 'to_id', are the node ids 
        of the two extremities of an edge (link), the next columns (if any) are
        additional scalar properties (attributes) attached to the edges.

        Parameters
        ----------
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

        export_node_scalar_properties : bool, default: True
            indicates if the properties (other than 'pos') associated to 
            the nodes (if present) are exported; all properties are assumed to
            be scalar

        export_edge_scalar_properties : bool, default: True
            indicates if the properties associated to 
            the edges (if present) are exported; all properties are assumed to
            be scalar

        verbose : bool, default: True
            indicates if info is printed

        Examples
        --------
            >>> Kg.to_csv("MyKarst")
        """
        kn.io_func.networkx_to_csv(
            self.graph, 
            basename,
            suffix_nodes=suffix_nodes,
            suffix_edges=suffix_edges,
            delimiter_nodes=delimiter_nodes,
            delimiter_edges=delimiter_edges,
            pos_attr='pos',
            export_node_scalar_properties=export_node_scalar_properties,
            export_edge_scalar_properties=export_edge_scalar_properties,
            verbose=False)

        if verbose:
            print(f"Karstnet graph successfully exported to csv files !\n")

        return
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    def node_scalar_properties_to_csv(
            self, 
            filename,
            properties_name_list=None,
            delimiter=';',
            verbose=True):
        """
        Exports node scalar properties from the karstnet graph (original graph) to a csv file.

        The file has the columns:
            - 'id', '<node_prop0>'[, '<node_prop1>', ...]

        where the column 'id' is the id of the nodes and the next columns are 
        additional scalar properties (attributes) attached to the nodes.

        Parameters
        ----------
        filename : str
            the name of the output file
        
        properties_name_list : list of strs, optional
            list of property names that are exported;
            by default (`None`): all properties (except the position)
            are exported and assumed to be scalar
        
        delimiter : str, default: ';'
            delimiter used in the file
        
        verbose : bool, default: True
            indicates if info is printed

        Examples
        --------
            >>> Kg.node_scalar_properties_to_csv("MyKarst_node_properties.csv")
        """
        kn.io_func.networkx_node_scalar_properties_to_csv(
            self.graph,
            filename,
            properties_name_list=properties_name_list,
            delimiter=delimiter,
            verbose=False)

        if verbose:
            print(f"Karstnet graph: node scalar properties successfully exported to csv file ({filename}) !\n")

        return
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    def edge_scalar_properties_to_csv(
            self, 
            filename,
            properties_name_list=None,
            delimiter=';',
            verbose=True):
        """
        Exports edge scalar properties from the karstnet graph (original graph) to a csv file.

        The file has the columns:
            - 'from_id', 'to_id', '<edge_prop0>'[, '<edge_prop1>', ...]

        where the first two columns, 'from_id', 'to_id', are the node ids 
        of the two extremities of an edge (link), the next columns are
        additional scalar properties (attributes) attached to the edges.

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
        
        verbose : bool, default: True
            indicates if info is printed

        Examples
        --------
            >>> Kg.edge_scalar_properties_to_csv("MyKarst_edge_properties.csv")
        """
        kn.io_func.networkx_edge_scalar_properties_to_csv(
            self.graph,
            filename,
            properties_name_list=properties_name_list,
            delimiter=delimiter,
            verbose=False)

        if verbose:
            print(f"Karstnet graph edge scalar properties successfully exported to csv file ({filename}) !\n")

        return
    # -------------------------------------------------------------------------

    # *************************************************************************
    # Methods for export to pline (GOCAD ASCII)
    # *************************************************************************

    # -------------------------------------------------------------------------
    def to_pline(self, basename):
        """
        Exports the original graph to Pline (GOCAD ASCII object).

        Manages the colocated vertices indicated by the mention "ATOM" in
        the ASCII file.

        `Warning`: this version does not export the properties on the nodes.

        Parameters
        ----------
        basename : str
            the base name of the file name used for the Pline file,
            the name contains no extension, it will be added by the function

        Examples
        --------
        The following command saves the file "MyKarst_exported.pl" :

            >>> myKGraph.to_pline("MyKarst")
        """
        # For a original graph the list of Ilines corresponds to self.branches

        self._ilines_to_pline(self.branches, basename)

        return
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    def simpleGraph_to_pline(self, basename):
        """
        Exports the simplified graph to pline (GOCAD ASCII object).

        Manages the colocated vertices indicated by the mention "ATOM"
        in the ASCII file.

        `Warning`: this version does not export the properties on the nodes.

        Parameters
        ----------
        basename : str
            the base name of the file name used for the Pline file,
            the name contains no extension, it will be added by the function

        Examples
        --------
        The following command saves the file "MyKarst_simpl_exported.pl" :

            >>> myKGraph.to_pline("MyKarst")
        """
        # For a simplified graph the list of Ilines will corresponds to
        # the edges of the simple graph
        # one iline is created for each edge, which is not exactly a branch
        # but does not prevent exportation

        # to clearly explicit it is the simplified graph
        basename = basename + "_simpl_"
        self._ilines_to_pline(self.list_simpl_edges, basename)

        return
    # -------------------------------------------------------------------------

    # *************************************************************************
    # Methods for statistical analysis
    # *************************************************************************

    # -------------------------------------------------------------------------
    def basic_analysis(self):
        """
        Prints the basic statistics of a karstic network graph analysis.

        Examples
        --------
            >>> t = myKGraph.basic_analysis()
        """

        # On the original graph
        nb_nodes_comp = nx.number_of_nodes(self.graph)
        nb_edges_comp = nx.number_of_edges(self.graph)

        # On the simplified graph
        nb_nodes = nx.number_of_nodes(self.graph_simpl)
        nb_edges = nx.number_of_edges(self.graph_simpl)
        nb_connected_components = nx.number_connected_components(
            self.graph_simpl)

        nb_cycles = nb_edges - nb_nodes + nb_connected_components

        # Compute all extremities and junction nodes (on the simple graph)
        nb_extremity_nodes = 0
        nb_junction_nodes = 0
        for i in self.graph_simpl.nodes():
            if self.graph_simpl.degree(i) == 1:
                nb_extremity_nodes += 1
            elif self.graph_simpl.degree(i) > 2:
                nb_junction_nodes += 1

        # Print these basics
        print(
            "\nThis network (original graph) contains :\n",
            "   ", nb_nodes_comp, "nodes (stations) and ", nb_edges_comp, "edges.\n",
            "On the simplified graph, there are :\n",
            "   ", nb_nodes, "nodes (stations) and ", nb_edges, "edges,\n",
            "   ", nb_extremity_nodes, "are extremity nodes (entries or exits) and\n",
            "   ", nb_junction_nodes, "are junction nodes.\n"
            "There is/are", nb_connected_components, "connected component.s and", 
            nb_cycles, "cycle.s.\n")

        # Howard's parameters
        # (Howard, A. D., Keetch, M. E., & Vincent, C. L. (1970).
        # Topological and geometrical properties of braided patterns.
        # Water Resources Research, 6(6), 1674–1688.)
        # Rmq: All nodes of the simplified network have to be considered,
        # even those of degree 2 in the cycles or, it is not becoming
        # consistent with the exemples of braided rivers given by Howard.
        # This is indeed consistent with what has been done in
        # Collon, P., Bernasconi, D., Vuilleumier, C., & Renard, P. (2017).
        # Statistical metrics for the characterization of karst network
        # geometry and topology. Geomorphology, 283, 122–142.
        alpha = nb_cycles / (2 * (nb_nodes) - 5)
        beta = nb_edges / (nb_nodes)
        gamma = nb_edges / (3 * (nb_nodes - 2))
        print("\nHoward's parameter are (Howard, 1970) :",
              "\n   alpha :", alpha,
              "\n   beta :", beta,
              "\n   gamma :", gamma)
        print("\nNote that this computation considers the node of degree 2",
              "necessary to loop preservations as Seed Nodes, in order to",
              "stay consistent with Howard's illustrations.")
        return {'alpha': alpha, 'beta': beta, 'gamma': gamma}
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    def mean_tortuosity(self):
        """
        Computes the mean tortuosity of a karstic network.

        Returns
        -------
        t : float
            mean tortuosity of the branches

        Examples
        --------
            >>> t = myKGraph.mean_tortuosity()
        """
        nb_of_Nan = np.isnan(self.br_tort).sum()

        if self.verbose:
            if nb_of_Nan != 0:
                print(
                    "\nWARNING: This network contains",
                    nb_of_Nan,
                    "looping branche.s, which is.are not considered for the",
                    "mean tortuosity computation.")

        return np.nanmean(self.br_tort)
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    def mean_length(self):
        """
        Computes the mean length of the branches of a karstic networkx

        Returns
        -------
        l : float
            mean length of the branches

        Examples
        --------
            >>> l = myKGraph.mean_length()
        """
        return np.mean(self.br_lengths)
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    def coef_variation_length(self):
        """
        Computes the coefficient of variation of length of the branches of a
        karstic network.

        Returns
        -------
        cvl : float
            coefficient of variation of the length of the branches

        Examples
        --------
            >>> cvl = myKGraph.coef_variation_length()
        """

        # std is used with ddof=1 to use the estimate of std (divides by N-1)
        # The test below is to avoid an error with np.std computation
        # for networks having a single branch.
        if len(self.br_lengths) > 1:
            return np.std(self.br_lengths, ddof=1) / np.mean(self.br_lengths)
        else:
            return 0 # No
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    def length_entropy(self, mode="default"):
        r"""
        Computes the entropy of lengths of the branches of a karstic network.

        The entropy of the branch lengths, normalized in [0, 100] is computed
        as

        .. math::
            H=-\sum_i p_i \log_n(p_i)

        where :math:`p_i` is the probability mass in the `i`-th bin, and the
        base of logarithm, `n`, is the number of bins (see parameter `mode`
        below).

        Parameters
        ----------
        mode : str
            the mode style should be equal to "default" or "sturges";
            the property is normalized as compared to the
            maximum observed value

            - If mode="default", the entropy is computed according to the \
            implementation proposed in Collon et al. 2017 : \
            using a fixed bin number of 10 and ranging from 0 to 100
            - If mode="sturges" : alternative calculation using Sturges \
            rules to define the number of bins

        Returns
        -------
        entropy : float
            entropy of the lengths of the branches

        Examples
        --------
            >>> l_entrop = myKGraph.length_entropy()
        """

        v = self.br_lengths
        # In the paper of 2017, we normalize the length to get comparable
        # results between networks
        v = v[np.nonzero(v)]
        v_normalized = v / np.amax(v) * 100

        if (len(v) > 1):
            if mode == "sturges":
                # Sturges rule to define the number of bins fron nb of samples
                nbins = int(np.ceil(1 + np.log2(len(v_normalized))))
            else:
                # Fixed nb of bins to facilitate comparison,
                # as done in Collon et al 2017 and other papers on roads
                # orientation entropy
                nbins = 10

            # Pauline : range should be 0 and 100 or we get very different
            # results from what we should have
            # (I verify that 0 and 100 are included in the counts)
            counts, _ = np.histogram(v_normalized,
                                     bins=nbins,
                                     range=(0, 100))
            freq = counts / np.sum(counts)  # Computes the frequencies
            entropy = st.entropy(freq, base=nbins)
        else:
            entropy = 0  # v contains a single value - no uncertainty

        return entropy
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    def orientation_entropy(self, mode="default", weighted=True):
        r"""
        Computes the entropy of orientation (2d) of the segments of a
        karstic network.

        The 2D orientation is considered, i.e. the orientation of the
        projection of the edges (links) in xy plane; the entropy of the
        distribution of the azimuth angles in degrees in [0, 180) is
        computed as

        .. math::
            H=-\sum_i p_i \log_n(p_i)

        where :math:`p_i` is the probability mass in the `i`-th bin, and the
        base of logarithm, `n`, is the number of bins (see parameter `mode`
        below).

        By default, the distribution of azimuth is weighted by the length of the
        projection of the edges in the xy plane.

        Parameters
        ----------
        mode : string
            the mode style should be equal to "default" or "sturges"

            - If mode="default", the entropy is computed according to \
            the implementation proposed in Collon et al. 2017 : \
            using a fixed bin number of 18 and ranging from 0 to 180
            - If mode="sturges" : alternative calculation using Sturges \
            rules to define the number of bins

        weighted : bool, default: True
            indicates if the distribution are weighted by length (True), or if
            each edge has the same weight (False)

        Returns
        -------
        entropy : float
            entropy of the segment orientation

        Examples
        --------
            >>> or_entropy = myKGraph.orientation_entropy()
        """

        # Get azimuths
        azim = np.array(list((nx.get_edge_attributes(self.graph, 'azimuth')).values()))
        # Removing NAN Azimuth values (vertical edges, length 2d is zero)
        ind = ~np.isnan(azim)
        azim = azim[ind]
        azim = azim % 180 # identify angle azimuth angles that differ from 180 (Julien Straubhaar)

        if weighted:
            # Get length 2d for weights
            weights = np.array(list((nx.get_edge_attributes(self.graph, 'length2d')).values()))
            weights = weights[ind]
        else:
            weights = None

        if len(azim) > 1:
            if mode == "sturges":
                # Sturges rule to define the number of bins fron nb of samples
                nbins = int(np.ceil(1 + np.log2(len(azim))))
            else:
                # Fixed nb of bins to facilitate comparison,
                # as done in Collon et al 2017 and other papers on roads
                # orientation entropy
                nbins = 18

            # Pauline : range should be 0 and 180 or we get very different
            # results from what we had (change the position of bin borders).
            # Also, I verified that it is consistent :
            # 0 and 180 are counted, not excluded
            counts, _ = np.histogram(azim,
                                     bins=nbins,
                                     range=(0, 180),
                                     weights=weights)
            freq = counts / np.sum(counts)  # Computes the frequencies
            return st.entropy(freq, base=nbins)
        else:
            return 0
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    # Julien Straubhaar
    def orientation3d_entropy(self, nbins_azim=36, nbins_dip=6, weighted=True):
        r"""
        Computes the entropy of orientation of the segments of a
        karstic network.

        The 3D orientation is considered, i.e. orientation of an edge
        (link) is defined by the dip angle in degrees in [0, 90], and
        the azimuth angle in degrees in [0, 360).

        A pair of angles :math:`(\alpha, \theta)`, where :math:`\alpha`
        is the azimuth and :math:`\theta` the dip, defines a point on
        the lower hemisphere

        .. math::
            (\sin\alpha \cos\theta, \cos\alpha \cos\theta, -\sin\theta)

        The entropy of the joint distribution of the azimuth, dip angles
        is computed as

        .. math::
            H=-\sum_i p_i \log_n(p_i)

        where :math:`p_i` is the probability mass in the `i`-th bin, and the
        base of logarithm, `n`, is the number of bins.

        The bins cover the hemisphere, and every bin has the same area.
        Denoting :math:`n_\alpha` and :math:`n_\theta` the number of bins
        for azimuth and dip respectively (given as parameters), the total
        number of bins for the joint distribution is

        .. math::
            n = 1 + n_\alpha \cdot (n_\theta -1)

        where, with
        :math:`0=\alpha_0 < \alpha_1 <\cdots < \alpha_{n_\alpha}=360` and
        :math:`0=\theta_0 < \theta_1 <\cdots < \theta_{n_\theta}=90` the limits
        angles for azimuth and dip respectively, the 2D bins are

        .. math::
            \{(\alpha, \theta)\ :\ \theta_{n_{\theta-1}} \leqslant \theta\} \text{ "sphere cap"}\\
            \{(\alpha, \theta)\ :\ \alpha_{i} \leqslant \alpha < \alpha_{i+1}, \ \theta_{j} \leqslant \theta < \theta_{j+1}\}

        for :math:`i=0, \ldots, n_{\alpha-1}`, :math:`j=0, \ldots, n_{\theta-1}`.

        By default, the distribution of angles is weighted by the length of the
        of the edges.

        Parameters
        ----------
        nbins_azim : int, default: 36
            number of bins for azimuth angles (see above), at least 1

        nbins_dip : int, default: 6
            number of bins for dip angles (see above), at least 2

        weighted : bool, default: True
            indicates if the distribution are weighted by length (True), or if
            each edge has the same weight (False)

        Returns
        -------
        entropy : float
            entropy of the segment orientation

        Examples
        --------
            >>> or3d_entropy = myKGraph.orientation3d_entropy()
        """

        # Determine limits angles so that area of each bin is equal.
        #
        # With na = nbins_azim: alpha_i = i/(2pi), i=0, ..., na.
        #
        # With nt = nbins_dip: to determine theta_j, j=0, ..., nt.
        # Let
        #   theta_0=0, theta_{nt}=90
        # The area of a crown (on the hemisphere) defined
        # by an interval theta in [t1, t2], is equal to
        # sin(t2) - sin(t1).
        #
        # Then, with b = 1/nbins, and nbins = 1+na*(nt-1) the total
        # number of 2d bins, let
        #   theta_{nt-1} = arcsin(1-b)
        # such that the area of the sphere cap, for theta in
        # [theta_{nt_-1}, 90], is equal to b=1/nbins.
        #
        # Let
        #   theta_j = arcsin(j*na*b), j=0, ..., nt-1, (*)
        # such that the crown of the sphere, for theta in
        # [theta_{j}, theta_{j+1}], is equal to na*b.
        # Hence, as the crown will be divided in na parts (according
        # to bins for azimuth), each part will have an area of b.
        #
        # Note that (*) for j=0 and j=nt-1, give the same value
        # as the values f or theta_0 and theta_{nt-1} given before.
        #

        if nbins_azim < 1:
            print('ERROR: invalid bins for azimuth')
            return None
        if nbins_dip < 2:
            print('ERROR: invalid bins for dip')
            return None

        # Define bins limits
        nbins = 1 + nbins_azim*(nbins_dip-1)
        alphas = np.arange(nbins_azim+1)*2*np.pi/nbins_azim
        thetas = np.arcsin(np.arange(nbins_dip)*nbins_azim/nbins) # theta_{nbins_dip} = 90 not stored
        # thetas = np.hstack((np.arcsin(np.arange(nbins_dip)*nbins_azim/nbins), np.array([np.pi/2])))

        alphas = np.rad2deg(alphas)
        thetas = np.rad2deg(thetas)

        # Get azimuths and dips
        azim = np.array(list((nx.get_edge_attributes(self.graph, 'azimuth')).values()))
        dip = np.array(list((nx.get_edge_attributes(self.graph, 'dip')).values()))
        azim = np.nan_to_num(azim) # convert nan to zero (and inf to very large number

        if weighted:
            # Get length
            length = np.array(list((nx.get_edge_attributes(self.graph, 'length')).values()))

        # Compute probability mass in each bins (probability weighted by length of edges)
        counts = []
        for j in range(nbins_dip-1):
            ind = np.all((dip >= thetas[j], dip < thetas[j+1]), axis=0)
            if weighted:
                weights = length[ind]
            else:
                weights = None
            c, _ = np.histogram(azim[ind], bins=nbins_azim, range=(0, 360), weights=weights)
            counts.append(c)

        ind = dip >= thetas[nbins_dip-1]
        if weighted:
            c = np.array([np.sum(length[ind])])
        else:
            c = np.array([np.sum(ind)])
        counts.append(c)

        counts = np.hstack(counts)

        freq = counts / np.sum(counts)  # Computes the frequencies
        return st.entropy(freq, base=nbins)
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    def mean_degree_and_CV(self):
        """
        Computes the average and the coefficient of variation of the degree.

        The computation is done on the simplified graph.

        Returns
        -------
        tuple
            meandeg, cvdeg : the mean and coefficient of variation
            of the node degrees

        Examples
        --------
            >>> meandeg, cvde = myKGraph.coef_variation_degree()
        """
        # Vector of degrees
        d = np.asarray(list(dict(self.graph_simpl.degree()).values()))
        #d = np.asarray(self.graph_simpl.degree())[:, 1]
        
        # Mean degree
        meandeg = np.mean(d)

        # Coefficient of variation of the degrees : std is used with ddof=1 to
        # use the estimate of std (divides by N-1)
        cvde = np.std(d, ddof=1) / np.mean(d)

        return meandeg, cvde
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    def correlation_vertex_degree(self, cvde=None):
        """
        Computes the correlation of vertex degree.

        The computation is done on the simplified graph.

        Parameters
        ----------
        cvde : float (or bool: False), optional
            Optional input: coefficient of variation of the degree.
            If not provided (`None` or `False`), it is computed automatically internally.

        Returns
        -------
        cvd : float
            Correlation of Vertex Degree

        Examples
        --------
           >>> cvd = myKGraph.correlation_vertex_degree()
        """
        if cvde is None or not cvde:
            _, cvde = self.mean_degree_and_CV()

        # To avoid division by 0 when computing correlation coef
        if cvde != 0:
            cvd = nx.degree_pearson_correlation_coefficient(self.graph_simpl)
        else:
            cvd = 1

        return cvd
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    def central_point_dominance(self):
        """
        Computes central point dominance.

        The computation is done on the simplified graph.

        Returns
        -------
        cpd : float
            central point dominance

        Examples
        --------
            >>> cpd = myKGraph.central_point_dominance()
        """
        bet_cen = nx.betweenness_centrality(self.graph_simpl)
        bet_cen = list(bet_cen.values())
        cpd = np.sum(max(bet_cen) - np.array(bet_cen)) / (len(bet_cen) - 1)

        return cpd
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    def average_SPL(self, dist_weight=False):
        """
        Computes the average shortest path length.

        Notes
        -----

        The computation is done on the simplified graph.

        The function handles the case of several connected components
        which is not the case for the Networkx function
        "average_shortest_path_length".

        In case of several connected components, the average_SPL
        is the average of each SPL weighted by the number of nodes of each
        connected component.

        Parameters
        ----------
        dist_weight : bool, default: False
            if True, shortest path lengths are computed by weighting the
            edges by their length

        Returns
        -------
        aspl : float
            average shortest path length

        Examples
        --------
            >>> aspl = myKGraph.average_SPL()
        """

        sum_aspl = 0  # initialize the sum
        # Compute average spl on each connected component with Networkx
        for c in (self.graph_simpl.subgraph(c).copy()
                  for c in nx.connected_components(self.graph_simpl)):
            if not dist_weight:
                sum_aspl += nx.average_shortest_path_length(
                    c) * nx.number_of_nodes(c)
            else:
                sum_aspl += nx.average_shortest_path_length(
                    c, weight="length") * nx.number_of_nodes(c)

        av_SPL = sum_aspl / nx.number_of_nodes(self.graph_simpl)

        return av_SPL
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    def characterize_graph(self, compute_orientation_3d_entropy=False, verbose=False):
        """
        Computes the set of metrics used to characterize a graph.

        Parameters
        ----------
        compute_orientation_3d_entropy : bool, default : False
            indicates if the orientation 3d entropy is computed;
            by default (`False`): not computed (for back-compatibility)

        verbose : bool, default: False
            If True, the function displays information about the
            progress of the computation, and the results.

        Returns
        -------
        results : dict
            All the statistical metrics are stored in a dictionary. The
            keys of the dictionary are provided below with
            corresponding explanation.

            - `mean length` : mean length of the branches
            - `cv length` : coefficient of variation of length of branches
            - `length entropy` : entropy of the length of the branches
            - `mean tortuosity` : mean tortuosity of the branches
            - `orientation entropy` : entropy of the orientation (2d, azimuth) of the conduits
            - `orientation 3d entropy` : entropy of the orientation (3d, azimuth and dip) of the conduits
            - `aspl` : average shortest path length
            - `cpd` : central point dominance
            - `mean degree` : mean of the vertex degrees
            - `cv degrees` : coefficient of variation of vertex degrees
            - `correlation vertex degree` :  correlation of vertex degrees

        Examples
        --------
            >>> results = myKGraph.characterize_graph()
        """

        results = {}

        if verbose:
            print('Computing:')
            print(' - mean length', end='', flush=True)

        results["mean length"] = self.mean_length()

        if verbose:
            print(', cv length', end='', flush=True)

        results["cv length"] = self.coef_variation_length()

        if verbose:
            print(', length entropy', end='', flush=True)

        results["length entropy"] = self.length_entropy()

        if verbose:
            print(', mean tortuosity', end='', flush=True)

        results["tortuosity"] = self.mean_tortuosity()

        if verbose:
            print('', end='\n', flush=True)
            print(' - orientation entropy', end='', flush=True)

        results["orientation entropy"] = self.orientation_entropy()

        if compute_orientation_3d_entropy:
            if verbose:
                print(', orientation 3d entropy', end='', flush=True)

            results["orientation 3d entropy"] = self.orientation3d_entropy()

        if verbose:
            print('', end='\n', flush=True)
            print(' - aspl', end='', flush=True)

        results["aspl"] = self.average_SPL()

        if verbose:
            print(', cpd', end='', flush=True)

        results["cpd"] = self.central_point_dominance()

        if verbose:
            print(', md, cv degree', end='', flush=True)

        md, cvde = self.mean_degree_and_CV()
        results["mean degree"] = md
        results["cv degree"] = cvde

        if verbose:
            print(', cvd', end='', flush=True)

        cvd = self.correlation_vertex_degree(cvde=cvde)
        results["correlation vertex degree"] = cvd

        if verbose:
            print('', end='\n', flush=True)
            print('\nResults:')
            print("--------------------------------------")
            for key in results.keys():
                print(f' {key:25s} = {results[key]:6.3f}')
                # print(" %25s = %5.3f" % (key, results[key]))
            print("--------------------------------------")

        return results
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    def upscale_edge_resistances(self, edge_property_name):
        """
        Upscales the edge resistances from the original graph to the simplified graph.
        
        As an edge in the simplified graph replaces a (part of a) branch (path) in 
        the original graph, the edge resistance in the simplified graph is simply the 
        sum of the edge resistances along the corresponding path in the original graph.
        
        Parameters
        ----------
        edge_property_name : str
            name of the edge property for resistance:
            
            - in the original graph : the property must exist
            - in the simplified graph : the property is created (or erased if already \
            existing) with the upscaled resistance
        """
        # Initialize dictionary of upscaled resistances
        r_upsc = {}

        for e in self.graph_simpl.edges():
            # Set path: part of the branch from e[0] to e[1]
            for b in self.branches:
                if e[0] in b and e[1] in b:
                    # Edge e is in branch b
                    for i, bi in enumerate(b):
                        if e[0] == bi:
                            i0 = i
                        if e[1] == bi:
                            i1 = i
                    if i1 > i0:
                        path = b[i0:i1+1]
                    else:
                        path = b[i1:i0+1]

            # Compute upscaled resistance for the edge e
            r_upsc[e] = np.sum(np.asarray([self.graph.edges[(a0, a1)][edge_property_name] for a0, a1 in zip(path[:-1], path[1:])]))
        
        # Set upscaled resistances as edge property in simlified graph
        nx.set_edge_attributes(self.graph_simpl, r_upsc, edge_property_name)
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    def upscale_edge_conductances(self, edge_property_name):
        """
        Upscales the edge conductances from the original graph to the simplified graph.
        
        As an edge in the simplified graph replaces a (part of a) branch (path) in 
        the original graph, the edge conductance in the simplified graph is simply the 
        inverse of the sum of the inverse of edge conductances along the corresponding 
        path in the original graph.
        
        Parameters
        ----------
        edge_property_name : str
            name of the edge property for conductance:
            
            - in the original graph : the property must exist
            - in the simplified graph : the property is created (or erased if already \
            existing) with the upscaled conductances
        """
        # Initialize dictionary of upscaled conductances
        w_upsc = {}

        for e in self.graph_simpl.edges():
            # Set path: part of the branch from e[0] to e[1]
            for b in self.branches:
                if e[0] in b and e[1] in b:
                    # Edge e is in branch b
                    for i, bi in enumerate(b):
                        if e[0] == bi:
                            i0 = i
                        if e[1] == bi:
                            i1 = i
                    if i1 > i0:
                        path = b[i0:i1+1]
                    else:
                        path = b[i1:i0+1]

            # Compute upscaled conductance for the edge e
            w_upsc[e] = 1.0 / np.sum(np.asarray([1.0 / self.graph.edges[(a0, a1)][edge_property_name] for a0, a1 in zip(path[:-1], path[1:])]))
        
        # Set upscaled conductances as edge property in simlified graph
        nx.set_edge_attributes(self.graph_simpl, w_upsc, edge_property_name)
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    def variogram_cloud_intra_branches(
            self, 
            node_attr,
            exclude_branch_extremities=True,
            exclude_nan=True,
            hmax=None, 
            make_plot=True, 
            **kwargs):
        """
        Computes the variogram cloud intra-branches for the given node attribute
        `node_attr` which is assumed to be a scalar (float) or a sequence 
        of floats.

        For each pair of nodes i, j within a branch (see `exclude_branch_extremities`
        and `exclude_nan` below), let
        
        - h(i, j) the length between nodes i and j
        - gamma(i, j) = 0.5*(v[i]-v[j])**2 (gamma values for the property) \
        where v[k] is the property scalar (or vector) of attached to node k

        This method retrieves all points (h(i, j), gamma(i, j)) for which h(i, j) is
        less than or equal to `hmax` (if given).

        Parameters
        ----------
        node_attr : str
            name of the node attribute for which the variogram cloud is computed;
            this node attribute should be a scalar (float) or a sequence of floats

        hmax : float, optional
            maximal distance between two nodes to be integrated in the variogram cloud

        exclude_branch_extremities : bool, default: `True`
            if `True`, any pair of nodes i, j within a branch, with i or j at an 
            extremity of the branch is excluded from the variogram cloud

        exclude_nan : bool, default: `True`
            if `True`, any pair of nodes i, j within a branch, with at least 
            one `nan` value for the considered attribute at these nodes, is 
            excluded from the variogram cloud

        make_plot : bool, default: `True`
            indicates if the variogram cloud is plotted (in the current
            figure axis)

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
        >>> h, gamma, npair = myKGraph.variogram_cloud_intra_branches('node_prop')
        """
        # Set dictionary to convert node label (id) to node index, and vice versa
        node_label2index = {u:i for i, u in enumerate(self.graph.nodes())}
        #node_index2label = {i:u for i, u in enumerate(self.graph.nodes())}

        # Get the 2d-array of property values, of shape (n_nodes, n_prop)
        # where n_nodes is the number of nodes in the graph and 
        # n_prop is the number of properties (1 if scalar property)
        n_nodes = self.graph.number_of_nodes()
        
        prop_dict = nx.get_node_attributes(self.graph, node_attr)
        if len(prop_dict) == 0:
            raise ValueError('No value for specified property')
        
        v_first = np.asarray(list(prop_dict.values())[0]) # first value entry as array
        if v_first.ndim > 1:
            raise ValueError('Value of the property should by a scalar or a sequence of scalar (value is an array of dimension greater than 1)')
        
        n_prop = np.asarray(v_first).size # size of first value entry

        if len(prop_dict) != n_nodes:
            # There is missing value for some nodes, set default (np.nan)
            if v_first.ndim == 0:
                # attribute is a scalar (n_prop = 1)
                prop_dict = nx.get_node_attributes(self.graph, node_attr, default=np.nan)
            else: # v_first.ndim == 1
                # attribute is a sequence
                prop_dict = nx.get_node_attributes(self.graph, node_attr, default=np.full(n_prop, np.nan))

        try:
            v = np.asarray(list(prop_dict.values())).reshape(self.graph.number_of_nodes(), -1)
        except:
            raise ValueError('Unable to set property (entry for every node should have the same shape)')

        if np.all(np.isnan(v), axis=0).any():
            raise ValueError('The property (or one of its component) has only undefined (`nan`) values')

        h, gamma = np.array([]), np.empty(shape=(0, n_prop))

        branches = self.branches
        edge_lengths = nx.get_edge_attributes(self.graph, 'length')
        for br in branches:
            if exclude_branch_extremities:
                br = br[1:-1]

            if len(br) <= 1:
                continue

            br_local_lengths = []
            for i in range(len(br)-1):
                if (br[i], br[i+1]) in edge_lengths.keys():
                    br_l = edge_lengths[(br[i], br[i+1])]
                else:
                    br_l = edge_lengths[(br[i+1], br[i])]

                br_local_lengths.append(br_l)

            br_prop = np.asarray([v[node_label2index[i]] for i in br])
            for i in range(len(br)-1):
                br_h = np.cumsum(br_local_lengths[i:])
                br_gamma = 0.5*(br_prop[i:i+1] - br_prop[i+1:])**2
                
                if exclude_nan:
                    ind = ~np.isnan(br_gamma).any(axis=1)
                
                br_h = br_h[ind]
                br_gamma = br_gamma[ind]
                
                if hmax is not None:
                    ind = br_h <= hmax
                    br_h = br_h[ind]
                    br_gamma = br_gamma[ind]
                
                h = np.hstack((h, br_h))
                gamma = np.vstack((gamma, br_gamma))

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
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    def variogram_exp_intra_branches(
            self,
            node_attr=None,
            exclude_branch_extremities=True,
            exclude_nan=True,
            hmax=None, 
            ncla=10,
            cla_center=None,
            cla_length=None,
            variogramCloud=None,
            make_plot=True,
            show_count=True,
            **kwargs):
        """
        Computes the experimental variogram intra-branches for the given node 
        attribute `node_attr` which is assumed to be a scalar (float) or a sequence 
        of floats.

        First, the variogram cloud intra-branches is computed with the method
        `variogram_cloud_intra_branches`, or directly retrieved from the 
        parameter `variogramCloud` if specified.

        The variogram cloud intra-branches consists of all points 
        (h(i, j), gamma(i, j)), with h the distance and gamma the gamma value, 
        for which the pair of nodes i, j is on the same branch and for which h(i, j) 
        is less than or equal to `hmax` (if given) (see method
        `variogram_cloud_intra_branches`).

        The mean point in each class is retrieved from the variogram cloud; the
        i-th class is determined by its center `cla_center[i]` and its length
        `cla_length[i]`, and corresponds to the interval
        `]cla_center[i]-cla_length[i]/2, cla_center[i]+cla_length[i]/2]`
        along h (lag) axis (abscissa).

        Parameters
        ----------
        node_attr : str
            name of the node attribute for which the variogram cloud is computed;
            this node attribute should be a scalar (float) or a sequence of floats

        hmax : float, optional
            maximal distance between two nodes to be integrated in the variogram cloud

        exclude_branch_extremities : bool, default: `True`
            if `True`, any pair of nodes i, j within a branch, with i or j at an 
            extremity of the branch is excluded from the variogram cloud

        exclude_nan : bool, default: `True`
            if `True`, any pair of nodes i, j within a branch, with at least 
            one `nan` value for the considered attribute at these nodes, is 
            excluded from the variogram cloud

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
        
        variogramCloud : 3-tuple, optional
            if given (not `None`), `variogramCloud`=(h, gamma, npair) is a variogram
            cloud (already computed and returned by the method
            `variogram_cloud_intra_branches`; in this case, `node_attr` is not used;
            by default (`None`): the variogram cloud is computed by using the
            method `variogram_cloud_intra_branches`

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
            (i.e. if the specified attribute (`node_attr`) is a scalar, \
            or if `variogramCloud[1]` is a 1d array), then `hexp` is a 1d array
            - if experimental variogram is computed for more than one scalar variable \
            (i.e. if the specified attribute (`node_attr`) is a vector, \
            or if `variogramCloud[1]` is a 2d array), then `hexp` is a 2d array \
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

        Examples
        --------
        >>> # Example for a scalar property 'node_prop' on KGraph `myKGraph`
        >>> # 1. Plot experimental variogram only
        >>> hexp, gexp, cexp = myKGraph.variogram_exp_intra_branches('node_prop')
        >>> #
        >>> # 2. Plot variogram cloud and experimental variogram
        >>> plt.figure(figsize=(12,6))
        >>> h, gamma, npair = myKGraph.variogram_cloud_intra_branches('node_prop', color='lightblue')
        >>> hexp, gexp, cexp = myKGraph.variogram_exp_intra_branches(variogramCloud=(h, gamma, npair), color='tab:blue')
        >>> plt.grid()
        >>> #plt.xlabel('h')
        >>> #plt.ylabel(r'$1/2(Z(x)-Z(x+h))^2$')
        >>> plt.show()
        """

        if variogramCloud is None:
            try:
                h, gamma, npair = self.variogram_cloud_intra_branches(
                                        node_attr,
                                        exclude_branch_extremities=exclude_branch_extremities,
                                        exclude_nan=exclude_nan,
                                        hmax=hmax, 
                                        make_plot=False)
            except Exception as exc:
                raise ValueError('variogram cloud can not be computed') from exc


        else:
            try:
                h, gamma, npair = variogramCloud
            except:
                raise ValueError('`variogramCloud` can not be used')

        hexp, gexp, cexp = kn.tools.variogram_exp_from_variogram_cloud(
                                (h, gamma, npair),
                                hmax=hmax,
                                ncla=ncla,
                                cla_center=cla_center,
                                cla_length=cla_length,
                                make_plot=make_plot,
                                show_count=show_count,
                                **kwargs)
        
        return hexp, gexp, cexp
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    def variogram_cloud(
            self, 
            node_attr,
            nsample_start_nodes,
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
        of nodes randomly chosen in the entire original graph.

        Parameters
        ----------
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
            if `True`, any pair of nodes i, j, with at least 
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
        h, gamma, npair = kn.tools.variogram_cloud(
                            self.graph,
                            node_attr,
                            nsample_start_nodes,
                            edge_length_attr='length',
                            exclude_nan=exclude_nan,
                            hmax=hmax,
                            seed=seed,
                            make_plot=make_plot,
                            verbose=verbose,
                            **kwargs)

        return h, gamma, npair
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    def variogram_exp(
            self,
            node_attr=None,
            nsample_start_nodes=None,
            exclude_nan=True,
            hmax=None,
            seed=None,
            ncla=10,
            cla_center=None,
            cla_length=None,
            variogramCloud=None,
            make_plot=True,
            show_count=True,
            **kwargs):
        """
        Computes the experimental variogram for the given node attribute `node_attr` 
        which is assumed to be a scalar (float) or a sequence of floats.

        First, the variogram cloud in the entire original graph is computed with the 
        method `variogram_cloud`, or directly retrieved from the parameter `variogramCloud` 
        if specified.

        The variogram cloud consists of all points (h(i, j), gamma(i, j)), with h the distance 
        and gamma the gamma value, for any pair of nodes i, j in the original graph for 
        which h(i, j) is less than or equal to `hmax` (if given) (see method
        `variogram_cloud`).

        The mean point in each class is retrieved from the variogram cloud; the
        i-th class is determined by its center `cla_center[i]` and its length
        `cla_length[i]`, and corresponds to the interval
        `]cla_center[i]-cla_length[i]/2, cla_center[i]+cla_length[i]/2]`
        along h (lag) axis (abscissa).

        Parameters
        ----------
        node_attr : str, optional
            name of the node attribute for which the variogram cloud is computed;
            this node attribute should be a scalar (float) or a sequence of floats;
            if not specified (`None`), `variogramCloud` is used (must be specified)

        nsample_start_nodes : int, optional
            number of sample start nodes (randomly chosen) in the original graph; for
            each sample start node i, all nodes j (not equal to i), at distance less 
            than or equal to `hmax` (if specified), are considered for a pair (i, j)
            to be included in the variogram cloud;
            if not specified (`None`), `variogramCloud` is used (must be specified)

        hmax : float, optional
            maximal distance between two nodes to be integrated in the variogram cloud

        exclude_nan : bool, default: `True`
            if `True`, any pair of nodes i, j, with at least 
            one `nan` value for the considered attribute at these nodes, is 
            excluded from the variogram cloud

        seed : int, optional
            seed to initialize the random number generator; used if
            `variogramCloud` is not specified

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
        
        variogramCloud : 3-tuple, optional
            if given (not `None`), `variogramCloud`=(h, gamma, npair) is a variogram
            cloud (already computed and returned by the method
            `variogram_cloud_intra_branches`; in this case, `node_attr` is not used;
            by default (`None`): the variogram cloud is computed by using the
            method `variogram_cloud_intra_branches`

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
            (i.e. if the specified attribute (`node_attr`) is a scalar, \
            or if `variogramCloud[1]` is a 1d array), then `hexp` is a 1d array
            - if experimental variogram is computed for more than one scalar variable \
            (i.e. if the specified attribute (`node_attr`) is a vector, \
            or if `variogramCloud[1]` is a 2d array), then `hexp` is a 2d array \
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

        Examples
        --------
        >>> # Example for a scalar property 'node_prop' on KGraph `myKGraph`
        >>> # 1. Plot experimental variogram only
        >>> hexp, gexp, cexp = myKGraph.variogram_exp('node_prop', nsample_pair=2000, seed=924)
        >>> #
        >>> # 2. Plot variogram cloud and experimental variogram
        >>> plt.figure(figsize=(12,6))
        >>> h, gamma, npair = myKGraph.variogram_cloud('node_prop', nsample_pair=2000, seed=924, color='lightblue')
        >>> hexp, gexp, cexp = myKGraph.variogram_exp(variogramCloud=(h, gamma, npair), color='tab:blue')
        >>> plt.grid()
        >>> #plt.xlabel('h')
        >>> #plt.ylabel(r'$1/2(Z(x)-Z(x+h))^2$')
        >>> plt.show()
        """

        if variogramCloud is None:
            h, gamma, npair = self.variogram_cloud(                        
                                    node_attr,
                                    nsample_start_nodes,
                                    exclude_nan=exclude_nan,
                                    hmax=hmax,
                                    seed=seed,
                                    make_plot=False)
        else:
            try:
                h, gamma, npair = variogramCloud
            except:
                raise ValueError('`variogramCloud` can not be used')

        hexp, gexp, cexp = kn.tools.variogram_exp_from_variogram_cloud(
                                (h, gamma, npair),
                                hmax=hmax,
                                ncla=ncla,
                                cla_center=cla_center,
                                cla_length=cla_length,
                                make_plot=make_plot,
                                show_count=show_count,
                                **kwargs)
        
        return hexp, gexp, cexp
    # -------------------------------------------------------------------------

    # # -------------------------------------------------------------------------
    # def variogram_cloud(
    #         self, 
    #         node_attr,
    #         nsample_pair,
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
    #     of nodes randomly chosen in the entire original graph.

    #     Parameters
    #     ----------
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
    #         if `True`, any pair of nodes i, j, with at least 
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
    #     h, gamma, npair = kn.tools.variogram_cloud(
    #                         self.graph,
    #                         node_attr,
    #                         nsample_pair,
    #                         edge_length_attr='length',
    #                         exclude_nan=exclude_nan,
    #                         hmax=hmax,
    #                         seed=seed,
    #                         make_plot=make_plot,
    #                         verbose=verbose,
    #                         **kwargs)

    #     return h, gamma, npair
    # # -------------------------------------------------------------------------

    # # -------------------------------------------------------------------------
    # def variogram_exp(
    #         self,
    #         node_attr=None,
    #         nsample_pair=None,
    #         exclude_nan=True,
    #         hmax=None,
    #         seed=None,
    #         ncla=10,
    #         cla_center=None,
    #         cla_length=None,
    #         variogramCloud=None,
    #         make_plot=True,
    #         show_count=True,
    #         **kwargs):
    #     """
    #     Computes the experimental variogram for the given node attribute `node_attr` 
    #     which is assumed to be a scalar (float) or a sequence of floats.

    #     First, the variogram cloud in the entire original graph is computed with the 
    #     method `variogram_cloud`, or directly retrieved from the parameter `variogramCloud` 
    #     if specified.

    #     The variogram cloud consists of all points (h(i, j), gamma(i, j)), with h the distance 
    #     and gamma the gamma value, for any pair of nodes i, j in the original graph for 
    #     which h(i, j) is less than or equal to `hmax` (if given) (see method
    #     `variogram_cloud`).

    #     The mean point in each class is retrieved from the variogram cloud; the
    #     i-th class is determined by its center `cla_center[i]` and its length
    #     `cla_length[i]`, and corresponds to the interval
    #     `]cla_center[i]-cla_length[i]/2, cla_center[i]+cla_length[i]/2]`
    #     along h (lag) axis (abscissa).

    #     Parameters
    #     ----------
    #     node_attr : str, optional
    #         name of the node attribute for which the variogram cloud is computed;
    #         this node attribute should be a scalar (float) or a sequence of floats;
    #         if not specified, `variogramCloud` is used (must be specified)

    #     nsample_pair : int, optional
    #         number of sample pairs of nodes (randomly chosen) in the original graph
    #         if not specified, `variogramCloud` is used (must be specified)

    #     hmax : float, optional
    #         maximal distance between two nodes to be integrated in the variogram cloud

    #     exclude_nan : bool, default: `True`
    #         if `True`, any pair of nodes i, j, with at least 
    #         one `nan` value for the considered attribute at these nodes, is 
    #         excluded from the variogram cloud

    #     seed : int, optional
    #         seed to initialize the random number generator; used if
    #         `variogramCloud` is not specified

    #     ncla : int, default: 10
    #         number of classes, the parameter is used if `cla_center=None`, in that
    #         situation `ncla` classes are considered and the class centers are set to
            
    #         - `cla_center[i] = (i+0.5)*l, i=0,...,ncla-1`
            
    #         with l = H / ncla, H being the max of the distance between two points of
    #         the considered pairs (in the variogram cloud);
    #         if `cla_center` is specified (not `None`), the number of classes (`ncla`)
    #         is set to the length of the sequence `cla_center` (ignoring the value
    #         passed as argument)
       
    #     cla_center : 1D array-like of floats, optional
    #         sequence of floats, center of each class (in abscissa) in the experimental
    #         variogram; by default (`None`): `cla_center` is defined from `ncla` (see
    #         above)

    #     cla_length : 1D array-like of floats, or float, optional
    #         length of each class centered at `cla_center` (in abscissa) in the
    #         experimental variogram:

    #         - if `cla_length` is a sequence, it should be of length `ncla`
    #         - if `cla_length` is a float, the value is repeated `ncla` times
    #         - if `cla_length=None` (default), the minimum of difference between two \
    #         sucessive class centers (`np.inf` if one class) is used and repeated `ncla` \
    #         times
        
    #     variogramCloud : 3-tuple, optional
    #         if given (not `None`), `variogramCloud`=(h, gamma, npair) is a variogram
    #         cloud (already computed and returned by the method
    #         `variogram_cloud_intra_branches`; in this case, `node_attr` is not used;
    #         by default (`None`): the variogram cloud is computed by using the
    #         method `variogram_cloud_intra_branches`

    #     make_plot : bool, default: `True`
    #         indicates if the experimental variogram is plotted (in the current
    #         figure axis)

    #     show_count : bool, default: `True`
    #         indicates if counters (`cexp`) are displayed on the plot
    #         (used if `make_plot=True`)

    #     kwargs : dict
    #         keyword arguments passed to the function `matplotlib.pyplot.plot`
    #         (used if `make_plot=True`)

    #     Returns
    #     -------
    #     hexp : 1d or 2d numpy array of floats
    #         distance (1st coordinate) of the experimental variogram

    #         - if experimental variogram is computed for one scalar variable \
    #         (i.e. if the specified attribute (`node_attr`) is a scalar, \
    #         or if `variogramCloud[1]` is a 1d array), then `hexp` is a 1d array
    #         - if experimental variogram is computed for more than one scalar variable \
    #         (i.e. if the specified attribute (`node_attr`) is a vector, \
    #         or if `variogramCloud[1]` is a 2d array), then `hexp` is a 2d array \
    #         where each column corresponds to each variable 

    #     gexp : 1d or 2d numpy array of floats
    #         gamma values (2nd coordinate) of the experimental variogram

    #         - the shape of the array `gexp` is the same as the shape of `hexp` \
    #         (see `hexp` above)

    #     cexp : 1d or 2d numpy array of ints
    #         counters, i.e. number of points (pairs of data points considered) 
    #         in each class in the variogram cloud

    #         - the shape of the array `cexp` is the same as the shape of `hexp` \
    #         (see `hexp` above)

    #     Examples
    #     --------
    #     >>> # Example for a scalar property 'node_prop' on KGraph `myKGraph`
    #     >>> # 1. Plot experimental variogram only
    #     >>> hexp, gexp, cexp = myKGraph.variogram_exp('node_prop', nsample_pair=2000, seed=924)
    #     >>> #
    #     >>> # 2. Plot variogram cloud and experimental variogram
    #     >>> plt.figure(figsize=(12,6))
    #     >>> h, gamma, npair = myKGraph.variogram_cloud('node_prop', nsample_pair=2000, seed=924, color='lightblue')
    #     >>> hexp, gexp, cexp = myKGraph.variogram_exp(variogramCloud=(h, gamma, npair), color='tab:blue')
    #     >>> plt.grid()
    #     >>> #plt.xlabel('h')
    #     >>> #plt.ylabel(r'$1/2(Z(x)-Z(x+h))^2$')
    #     >>> plt.show()
    #     """

    #     if variogramCloud is None:
    #         h, gamma, npair = self.variogram_cloud(                        
    #                                 node_attr,
    #                                 nsample_pair,
    #                                 exclude_nan=exclude_nan,
    #                                 hmax=hmax,
    #                                 seed=seed,
    #                                 make_plot=False)
    #     else:
    #         try:
    #             h, gamma, npair = variogramCloud
    #         except:
    #             raise ValueError('`variogramCloud` can not be used')

    #     hexp, gexp, cexp = kn.tools.variogram_exp_from_variogram_cloud(
    #                             (h, gamma, npair),
    #                             hmax=hmax,
    #                             ncla=ncla,
    #                             cla_center=cla_center,
    #                             cla_length=cla_length,
    #                             make_plot=make_plot,
    #                             show_count=show_count,
    #                             **kwargs)
        
    #     return hexp, gexp, cexp
    # # -------------------------------------------------------------------------

    # *************************************************************************
    # Non Public member functions of KGraph class
    # *************************************************************************

    # # =========================================================================
    # # Private functions for plots - OLD
    # # =========================================================================
    # # -------------------------------------------------------------------------
    # def _plot2(self,
    #            G,
    #            with_labels=True,
    #            node_size=300,
    #            node_color='lightblue',
    #            show_ticks=False,
    #            axis_equal=False,
    #            figsize=(6, 3)):
    #     """
    #     NOT PUBLIC
    #     Plots a 2D view of a graph G that could be the simplified of the
    #     original one.
    #     Requires self.pos2d() member function.
    #     Called by the plot2() public function.
    #     """

    #     # 2D  plot
    #     fig = plt.figure(figsize=figsize)

    #     nx.draw_networkx(G,
    #                      with_labels=with_labels,
    #                      pos=self.pos2d(),
    #                      node_size=node_size,
    #                      node_color=node_color)
    #     if show_ticks:
    #         plt.gca().tick_params(left=True, bottom=True, labelleft=True, labelbottom=True)
    #     if axis_equal:
    #         plt.axis('equal')

    #     return fig
    # # -------------------------------------------------------------------------

    # # -------------------------------------------------------------------------
    # def _plot3(self, G, zrotation=30, xyrotation=0, show_nodes=False, figsize=(8, 5)):
    #     """
    #     NOT PUBLIC
    #     Plots a 3D view of a graph G that could be the simplified or the
    #     original one.
    #     Requires self.pos3d() member function.
    #     Called by the plot3() public function.
    #     """

    #     # 3D  plot
    #     # try:
    #     #     from mpl_toolkits.mplot3d import Axes3D
    #     # except ImportError:
    #     #     raise ImportError("karstnet.plot3 requires mpl_toolkits.mplot3d ")

    #     # fig = plt.figure(figsize=figsize)
    #     # ax = Axes3D(fig)

    #     fig = plt.figure(figsize=figsize)
    #     ax = fig.add_subplot(projection='3d')

    #     pos3d = self.pos3d()
    #     for e in G.edges():
    #         x = np.array((pos3d[e[0]][0], pos3d[e[1]][0]))
    #         y = np.array((pos3d[e[0]][1], pos3d[e[1]][1]))
    #         z = np.array((pos3d[e[0]][2], pos3d[e[1]][2]))

    #         # Plot the connecting lines
    #         ax.plot(x, y, z, c='black', alpha=0.5)

    #     if show_nodes:
    #         ax.scatter(*np.asarray(list(pos3d.values())).T, c='blue', alpha=0.5)

    #     # Set the view
    #     ax.view_init(zrotation, -xyrotation - 90)

    #     ax.set_xlabel('X')
    #     ax.set_ylabel('Y')
    #     ax.set_zlabel('Z')

    #     return fig
    # # -------------------------------------------------------------------------
    
    # =========================================================================
    # Private function for export
    # =========================================================================
    # -------------------------------------------------------------------------
    def _ilines_to_pline(self, list_Iline, basename):
        """
        Writes a Pline (Gocad ascii object) from a list of ilines in a file.

        Used to export either original or simplified graph.

        Manages the colocated vertices indicated by the mention "ATOM"
        in the ASCII file.

        WARNING: this version does not export (for the moment)
        the properties on the nodes.

        Parameters
        ----------
        list_Iline : list
            list of the ilines to write

        basename: str
            string containing the base name of output file;
            the output file name is `basename`_exported.pl

        Examples
        --------
            >>> myKGraph.to_Pline("MyKarst" )
        """

        # Pline file creation
        output_file_name = basename + '_exported.pl'
        f_pline = open(output_file_name, 'w')

        # Header writing
        f_pline.write('GOCAD PLine 1\n')
        f_pline.write('HEADER {\n')
        f_pline.write('name:' + output_file_name + '\n')
        f_pline.write('}\n')
        f_pline.write('GOCAD_ORIGINAL_COORDINATE_SYSTEM\n')
        f_pline.write('NAME Default\nAXIS_NAME "U" "V" "W"\n')
        f_pline.write('AXIS_UNIT "m" "m" "m"\n')
        f_pline.write('ZPOSITIVE Elevation\n')
        f_pline.write('END_ORIGINAL_COORDINATE_SYSTEM\n')
        f_pline.write('PROPERTY_CLASS_HEADER Z {\n')
        f_pline.write('is_z:on\n}\n')

        # Create an empty dictionary of nodes already written in one iline
        # key is the same as in dico_nodes, the node index
        # value is the corresponding number of vtrx in the iline file,
        # due to the specific numbering of this file format,
        # it is different of the node index
        dico_added_nodes = {}

        pos3d = self.pos3d()

        # To count vertices: in plines,
        # vertices are virtually duplicated in the counting
        cpt_vrtx = 1
        # Each branch would be an iline
        for iline in list_Iline:
            f_pline.write('ILINE\n')

            # Memorize counting state to write correctly the segments
            cpt_vrtx_deb = cpt_vrtx

            # Each node of a iline is written as a vertex or atom
            for node in iline:
                # First, verify that this node has not already been
                # added to choose between vrtx or atom
                if node not in dico_added_nodes:
                    f_pline.write('VRTX ' + str(cpt_vrtx) + ' ' +
                                  str(pos3d[node][0]) + ' ' +
                                  str(pos3d[node][1]) + ' ' +
                                  str(pos3d[node][2]) + '\n')
                    # Update dico_added_nodes to indicate that the node
                    # has already been declared
                    # and store the correct index in the pline domain
                    dico_added_nodes[node] = cpt_vrtx
                # if node is in dico_added_nodes, we must build an atom
                # refering to the vrtx number in the pline
                else:
                    f_pline.write('ATOM ' + str(cpt_vrtx) + ' ' +
                                  str(dico_added_nodes[node]) + '\n')
                # Update vrtx counting to treat the next node of the iline
                cpt_vrtx += 1
            # When all nodes of a branch have been written, write the list
            # of segments using new numbers
            for i in range(len(iline) - 1):
                f_pline.write('SEG ' + str(cpt_vrtx_deb + i) + ' ' +
                              str(cpt_vrtx_deb + i + 1) + '\n')
            # One Iline has been written, go to next one

        # All ilines have been written
        f_pline.write('END\n')

        if self.verbose:
            print('File created')

        # Close the file
        f_pline.close()
        return
    # -------------------------------------------------------------------------

    # =========================================================================
    # Private functions used by constructors
    # =========================================================================
    # -------------------------------------------------------------------------
    def _set_graph_lengths(self):
        """NON PUBLIC.
        Computes edge length at the creation of KGraph object.
        This function is called by all constructors.
        It updates graph.
        """

        pos = nx.get_node_attributes(self.graph, 'pos')

        # Creation of a dictionary to store the length of each edge
        length = {e:float(np.sqrt(np.sum((np.asarray(pos[e[1]])-np.asarray(pos[e[0]]))**2)))
                  for e in self.graph.edges()}
        # length = {}
        # for e in self.graph.edges():
        #     dx = self.pos3d[e[0]][0] - self.pos3d[e[1]][0]
        #     dy = self.pos3d[e[0]][1] - self.pos3d[e[1]][1]
        #     dz = self.pos3d[e[0]][2] - self.pos3d[e[1]][2]
        #     length[e] = np.sqrt(dx ** 2 + dy ** 2 + dz ** 2)
        # Storing the length as an edge attribute
        nx.set_edge_attributes(self.graph, length, 'length')

        return

    def _simplify_graph(self):
        """
        Constructs a simplified graph by removing nodes of degree 2, while
        preserving the topology.

        Member function:
          Use self.graph (with its "length" attribute on edges) and self.pos3d()
          Use self.branches produced by self._get_allbranches()

        Returns:
        --------
        list_simpl_edges : list
            list of the simple edges (necessary for export to pline)

        Gs : networkx.Graph
            the simplified output graph
        """

        # Deals with cycles and loops to ensure that topology is not changed
        simpl_edges = _split_branches(self.branches)

        # Creates a new empty graph
        Gs = nx.Graph()

        # list of simpl_edges for export
        list_simpl_edges = []

        # Fills the graph with simpl_edges
        for i in simpl_edges:
            list_simpl_edges.append([i[0], i[-1]])

            # Compute the length of the current edge
            l_edge = float(np.sum(np.asarray([self.graph.edges[(node0, node1)]['length']
                                              for node0, node1 in zip(i[:-1], i[1:])])))
            # Notes:
            #  - `self.graph.edges[(node0, node1)]` gets the edge (node0, node1) or (node1, node0),
            #    i.e. order of endpoints does not matter
            #  - if the dictionary is first extracted: `length = nx.get_edge_attributes(self.graph, 'length')`,
            #    then order in the key (node0, node1) or (node1, node0) matters, i.e. if the first key is not found,
            #    then the second key must be used !
            #
            # ----- using the extracted dictionary `length` -----
            # l_edge = 0
            # for node0, node1 in zip(i[:-1], i[1:]):
            #     local_edge = (node0, node1)
            #     if length.__contains__(local_edge):
            #         br_len += length[local_edge]
            #     else:
            #         local_edge = (node1, node0)
            #         if length.__contains__(local_edge):
            #             br_len += length[local_edge]
            #         else:
            #             print("Warning: could not find ",
            #                   "1 edge when computing length")
            # -----

            # Add the edge corresponding to the simpl_edges, with attribute 'length'
            Gs.add_edge(i[0], i[-1], length=l_edge)

        # # ----- Alternative -----
        # # Creates the dictionary for length of edges
        # edges_length = {}
        # # Fills the graph with simpl_edges

        # for i in simpl_edges:
        #     list_simpl_edges.append((i[0], i[-1]))

        #     # Compute the length of the current edge
        #     l_edge = np.sum(np.asarray([self.graph.edges[(node0, node1)]['length']
        #                                 for node0, node1 in zip(i[:-1], i[1:])]))

        #     edges_length[list_simpl_edges[-1]] = l_edge

        # # Stores the results
        # Gs.add_edges_from(list_simpl_edges)
        # nx.set_edge_attributes(Gs, edges_length, 'length')
        # # -----
        
        # Set node position from the original graph (removed nodes are ignored)
        nx.set_node_attributes(Gs, nx.get_node_attributes(self.graph, 'pos'), 'pos')

        return list_simpl_edges, Gs
    # -------------------------------------------------------------------------

    # # ----- Alternative -----
    # -------------------------------------------------------------------------
    # def _simplify_graph(self):
    #     """
    #     Constructs a simplified graph by removing nodes of degree 2, while
    #     preserving the topology.

    #     Member function:
    #       Use self.graph (with its "length" attribute on edges) and self.pos3d()
    #       Use self.branches produced by self._get_allbranches()

    #     Returns:
    #     --------
    #     list_simpl_edges : list
    #         list of the simple edges (necessary for export to pline)

    #     Gs : networkx.Graph
    #         the simplified output graph
    #     """

    #     # Deals with cycles and loops to ensure that topology is not changed
    #     simpl_edges = _split_branches(self.branches)

    #     # Creates a new empty graph
    #     Gs = nx.Graph()

    #     # list of simpl_edges for export
    #     list_simpl_edges = []

    #     # Creates the dictionary for length of edges
    #     edges_length = {}

    #     # Fills the graph with simpl_edges

    #     for i in simpl_edges:
    #         list_simpl_edges.append((i[0], i[-1]))

    #         # Compute the length of the current edge
    #         l_edge = np.sum(np.asarray([self.graph.edges[(node0, node1)]['length']
    #                                     for node0, node1 in zip(i[:-1], i[1:])]))

    #         edges_length[list_simpl_edges[-1]] = l_edge

    #     # Stores the results
    #     Gs.add_edges_from(list_simpl_edges)
    #     nx.set_edge_attributes(Gs, edges_length, 'length')

    #     return list_simpl_edges, Gs
    # -------------------------------------------------------------------------
    # # ----- Alternative -----

    # -------------------------------------------------------------------------
    def _get_graph_of_branches(self):
        """
        Constructs the graph of branches.

        The graph of branches is a graph such that:
        - a node represents a branch, labeled by integers from 0
        - an edge exists between two nodes if the two branches represented \
        by the two nodes have a common extremity

        with

        - node attributes :
            - 'pos' : sequence of 2 or 3 floats
                mean position of the original graph nodes in the branch
                represented by the node

            - 'length' : flaot
                real length of the branch represented by the node
            
        Member function:
          Use self.graph (with its "pos" attribute on node)
          Use self.branches and self.br_lengths produced by self._get_allbranches()
        
        Returns:
        --------
        Gb : networkx.Graph
            the output graph of branch
        """

        # Number of branches in the original graph: number of nodes in the graph of branches
        nbranches = len(self.branches)

        # Set the list of nodes (label) in the graph of branches
        branch_nodes = list(range(nbranches))

        # Set the list of edges in the graph of branches
        branch_edges = []
        for i in range(nbranches-1):
            for j in range(i+1, nbranches):
                if     self.branches[i][0] == self.branches[j][0] \
                    or self.branches[i][0] == self.branches[j][-1] \
                    or self.branches[i][-1] == self.branches[j][0] \
                    or self.branches[i][-1] == self.branches[j][-1]:
                    branch_edges.append((i, j))

        # Creates the graph of branches
        Gb = nx.Graph()
        Gb.add_nodes_from(branch_nodes)
        Gb.add_edges_from(branch_edges)

        # Set position (mean position of original graph nodes in the branch)
        pos = {i: list(map(float, np.asarray([self.graph.nodes[u]['pos'] for u in br]).mean(axis=0))) 
               for i, br in enumerate(self.branches)}
        nx.set_node_attributes(Gb, pos, 'pos')

        # Set length (length of the branch)
        length = {i: float(l) for i, l in enumerate(self.br_lengths)}
        nx.set_node_attributes(Gb, length, 'length')

        return Gb
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    def _getallbranches(self):
        """
        Constructs the list of all branches of the karstic graph self_graph.
        Computes lengths and tortuosities.
        """        
        # Initialisations
        target = []
        degreeTarget = []

        # Create one subgraph per connected components
        list_sub_gr = [self.graph.subgraph(c).copy()
                       for c in nx.connected_components(self.graph)]

        # Get list of nodes that are at the extremities of a branch, and their degree;
        # for loop with all nodes of degree 2, keep one node
        for sub_gr in list_sub_gr:
            #node_id = np.asarray(list(d.keys()))
            # Set dictionary to convert node index to node label
            node_index2label = {i:u for i, u in enumerate(sub_gr.nodes())}
            # Dictionary of nodes degree
            d = dict(sub_gr.degree())
            node_degree = np.asarray(list(d.values()))
            ind = node_degree != 2

            if ind.sum() == 0:
                # all nodes of degree 2: add last one
                ind[-1] = True

            # target.append(node_id[ind])
            target.append([node_index2label[i] for i in np.where(ind)[0]])
            degreeTarget.append(node_degree[ind])

        if len(target) == 0:
            # no target node
            if self.verbose:
                print("Warning:", 
                      "\n   This network contains zero branch.\n")

            branches = []
            br_lengths = []
            br_tort = []
            return branches, np.array(br_lengths), np.array(br_tort)

        # target = np.hstack(target)
        target = [u for list_t in target for u in list_t]
        degreeTarget = np.hstack(degreeTarget)

        # Identifies all the neighbors of those nodes,
        # to create all the initial paths
        listStartBranches = [[i, n] for i in target for n in self.graph.neighbors(i)]

        # Follow all these initial paths to get all the branches
        branches = []
        for path in listStartBranches:
            go = True
            # Check all existing branches to avoid adding a branch twice
            # if starting from other extremity
            for knownbranch in branches:
                if ((path[0] == knownbranch[-1]) & (path[1] == knownbranch[-2])):
                    go = False
                    break
            if go:
                new_branch = self._getbranch(path)
                # Be sure that starting and ending nodes of the new branch
                # do not correspond to ending and starting nodes respectively of
                # an existing branch: if it is the case, the order of nodes in
                # the new branch are reversed
                # [this could be removed, because this should not happen due to
                # the way `listStartBranches` is created]
                for knownbranch in branches:
                    if new_branch[0] == knownbranch[-1] and new_branch[-1] == knownbranch[0]:
                        new_branch = new_branch[::-1]
                        break
                # Append new branch
                branches.append(new_branch)

        # Compute the list of branch lengths and tortuosities
        br_lengths = []
        br_tort = []

        # length = nx.get_edge_attributes(self.graph, 'length') # not used see note in the for loop below
        for br in branches:

            # Computes the distance between extremities
            dist = np.sqrt(np.sum((np.asarray(self.graph.nodes[br[0]]['pos'])
                                   - np.asarray(self.graph.nodes[br[-1]]['pos']))**2))

            # Computes the length of the current branch
            br_len = np.sum(np.asarray([self.graph.edges[(node0, node1)]['length']
                                        for node0, node1 in zip(br[:-1], br[1:])]))
            # Notes:
            #  - `self.graph.edges[(node0, node1)]` gets the edge (node0, node1) or (node1, node0),
            #    i.e. order of endpoints does not matter
            #  - if the dictionary is extracted first: `length = nx.get_edge_attributes(self.graph, 'length')`,
            #    then order in the key (node0, node1) or (node1, node0) matters, i.e. if the first key is not found,
            #    then the second key must be used !
            #
            # ----- using the extracted dictionary `length` -----
            # br_len = 0
            # for node0, node1 in zip(br[:-1], br[1:]):
            #     local_edge = (node0, node1)
            #     if length.__contains__(local_edge):
            #         br_len += length[local_edge]
            #     else:
            #         local_edge = (node1, node0)
            #         if length.__contains__(local_edge):
            #             br_len += length[local_edge]
            #         else:
            #             print("Warning: could not find ",
            #                   "1 edge when computing length")
            # -----

            br_lengths.append(br_len)

            # Set the toruosity of the current branch
            if np.isclose(dist, 0.0):
                tort = np.nan
            else:
                tort = br_len / dist

            br_tort.append(tort)

        br_lengths = np.asarray(br_lengths)
        br_tort = np.asarray(br_tort)

        if self.verbose:
            nb_of_Nan = np.isnan(br_tort).sum()
            if nb_of_Nan:
                print(
                    "Warning:"
                    "\n   This network contains", nb_of_Nan, "looping branch.es.",
                    "\n   Tortuosity is infinite on a looping branch.",
                    "\n   The looping branches are not considered",
                    "for the mean tortuosity computation.\n")
                # print(
                #     "Warning: This network contains",
                #     nb_of_Nan,
                #     "looping branch.es.",
                #     "Tortuosity is infinite on a looping branch.",
                #     "The looping branches are not considered",
                #     "for the mean tortuosity computation.\n")

        return branches, np.array(br_lengths), np.array(br_tort)
    # -------------------------------------------------------------------------

    # =========================================================================
    # Functions related to branches of graphs
    # - a branch is defined between two nodes of degree != 2
    # - exception: loop with all nodes of degree 2, one node is
    #   retained for both extremities of the branch
    # =========================================================================
    # -------------------------------------------------------------------------
    def _nextstep(self, path):
        """
        Works on self.graph
        Adds the next node to a path of self_graph along a branch.
        Stops when reaches a node of degree different from 2.
        """

        current = path[-1]
        # Checks first if the end of the path is already on an end
        if self.graph.degree(current) != 2:
            stopc = False
            return path, stopc

        # This is a security / it may be removed
        if len(path) > 1:
            old = path[-2]
        else:
            old = current

        # Among the neighbors search for the next one
        for nextn in self.graph.neighbors(current):
            if old != nextn:
                break

        # Add the next node to the path and check stopping criteria
        path.append(nextn)

        # Test for a closed loop / even if start node has degree = 2
        testloop = path[0] == path[-1]

        if (self.graph.degree(nextn) != 2) or testloop:
            stopc = False
        else:
            stopc = True

        return path, stopc
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    def _getbranch(self, path):
        """
        Works on self.graph
        Constructs a branch from a starting node.
        """

        path, stopc = self._nextstep(path)
        while stopc:
            path, stopc = self._nextstep(path)
        return path
    # -------------------------------------------------------------------------

    # =========================================================================
    # Private functions used for orientations
    # =========================================================================
    # -------------------------------------------------------------------------
    def _set_graph_orientations(self):
        """NON PUBLIC.
        Computes edge length at the creation of KGraph object.
        This function is called by all constructors.
        It updates graph.
        """

        pos3d = self.pos3d()

        # Creation of a dictionary to store the projected length of each edge,
        # the dip and the azimuth
        length2d = {}
        dip = {}
        azimuth = {}
        for e in self.graph.edges():
            dx = pos3d[e[0]][0] - pos3d[e[1]][0]
            dy = pos3d[e[0]][1] - pos3d[e[1]][1]
            dz = pos3d[e[0]][2] - pos3d[e[1]][2]
            length2d[e] = np.sqrt(dx ** 2 + dy ** 2)

            if length2d[e] != 0:
                dip[e] = np.arctan(dz / length2d[e])  # returns in radians
                dip[e] = np.degrees(dip[e])
                azimuth[e] = np.pi / 2 - np.arctan2(dy, dx) # returns in radians
                azimuth[e] = np.degrees(azimuth[e]) % 360

                if (dz < 0):
                    # negative cave gradients yield positive plunges in a stereonet
                    dip[e] = -dip[e]
                else:
                    # positive cave gradients have bearings 180 apart for the stereonet
                    azimuth[e] = (azimuth[e] + 180) % 360

            else:  # case of nearly pure vertical segments
                azimuth[e] = np.nan
                # azimuth[e] = 0.0 #Convention
                dip[e] = 90  # degrees

            length2d[e] = float(length2d[e])
            dip[e] = float(dip[e])
            azimuth[e] = float(azimuth[e]) 
            
        # Storing the length as an edge attribute
        nx.set_edge_attributes(self.graph, length2d, 'length2d')
        nx.set_edge_attributes(self.graph, azimuth, 'azimuth')
        nx.set_edge_attributes(self.graph, dip, 'dip')

        # azimuth: in [0, 360), from North, towards East (counter-clockwise)        -> bearings for stereonet
        # dip: in [0, 90], "plunging" direction, towards under the horizontal plane -> plunges for stereonet
        return
    # -------------------------------------------------------------------------

# ----- END of KGraph class --------------------------------------------------

# ****************************************************************************
# NON public functions used by KGraph
# (can not be placed in another file)
# ****************************************************************************

# ===== Functions used for graph simplification ==============================

# ----------------------------------------------------------------------------
def _split2(list_):
    """
    Splits a list in 2 sublists.
    """
    list_length = len(list_)
    if (list_length == 2):
        return (list_, [])
    else:
        midpoint = int(list_length / 2)
        return (list_[0:midpoint + 1], list_[midpoint:list_length])
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
def _split3(list_):
    """
    Splits a list in 3 sublists.
    """
    list_length = len(list_)
    if (list_length == 2):
        return (list_, [], [])
    elif (list_length == 3):
        l1, l2 = _split2(list_)
        return (l1, l2, [])
    else:
        k = int(list_length / 3.)
        return (list_[0:k + 1], list_[k:2 * k + 1], list_[2 * k:list_length])
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
def _split_branches(branches):
    """
    Split branches in cases of loop or cycles.
    """
    # Note: if two branches have the same endpoints, then
    # their starting nodes must be the same ones (and then also the
    # ending nodes); hence, any cycle with two junction nodes (of degree
    # greater than 2), is composed of two branches with the same
    # starting node (one of the junction node) and the same
    # ending node (the other junction node); this is guaranteed by
    # construction in method `_get_allbranches` of class:`KGraph`

    # Creates a dictionary to accelerate the search
    # of branches having the same extremities
    list_branches = dict()
    for i, b in enumerate(branches):
        key = (b[0], b[-1])
        if list_branches.__contains__(key):
            list_branches[key].append(i)
        else:
            list_branches[key] = [i]

    # Loop over the branches with same extremities and split them when required
    simpl_edges = []

    for key in list_branches:
        nbb = len(list_branches[key])
        # We test first if this is a loop (same start and end point)
        isloop = key[0] == key[1]

        # Simple case - no cycle but potential loop
        if nbb == 1:
            tmp = branches[list_branches[key][0]]
            if isloop:  # In the loop case, we need to split by 3 the branch
                tmp1, tmp2, tmp3 = _split3(tmp)
                simpl_edges.append(tmp1)
                if (len(tmp2) > 0):
                    simpl_edges.append(tmp2)
                    if (len(tmp3) > 0):
                        simpl_edges.append(tmp3)
            else:
                simpl_edges.append(tmp)

        # Several branches with same extremities - cycles and/or multiple loops
        else:
            if isloop:  # Case with multiple loops, we need to split by 3 each branch
                for i in range(nbb):
                    tmp = branches[list_branches[key][i]]
                    tmp1, tmp2, tmp3 = _split3(tmp)
                    simpl_edges.append(tmp1)
                    if (len(tmp2) > 0):
                        simpl_edges.append(tmp2)
                        if (len(tmp3) > 0):
                            simpl_edges.append(tmp3)
            else:  # Regular branches belonging to cycles, we need to split in 2 each branch
                for i in range(nbb):
                    tmp = branches[list_branches[key][i]]
                    tmp1, tmp2 = _split2(tmp)
                    simpl_edges.append(tmp1)
                    if (len(tmp2) > 0):
                        simpl_edges.append(tmp2)
                # # -----
                # # Note that the last branch could not be splitted, while preserving the topology
                # for i in range(nbb-1):
                #     tmp = branches[list_branches[key][i]]
                #     tmp1, tmp2 = _split2(tmp)
                #     simpl_edges.append(tmp1)
                #     if (len(tmp2) > 0):
                #         simpl_edges.append(tmp2)
                # # last branch
                # tmp = branches[list_branches[key][-1]]
                # simpl_edges.append(tmp)
                # # -----

    return simpl_edges
# ----------------------------------------------------------------------------

# ****************************************************************************
# Public functions used by KGraph
# (can not be placed in another file)
# ****************************************************************************

# ===== Import function from networkx graph ==================================
# Julien Straubhaar

# ----------------------------------------------------------------------------
def from_networkx(G, pos_attr='pos', verbose=True):
    """
    Creates a karstnet graph from a networkx graph.

    It is assumed that the given networkx graph contains at least one 
    node property corresponding to the position of the nodes.
    
    Parameters
    ----------
    G : networkx.Graph
        networkx graph (that contains at least one node property 
        corresponding to the position of the nodes)
    
    pos_attr : str, default: 'pos'
        name of the node attribute for position

    verbose : bool, default: True
        indicates if info is printed

    Returns
    -------
    Kg : karstnet.KGraph
        a KGraph object

    Examples
    --------
        >>> myKGraph = kn.from_networkx(G)
    """
    if pos_attr != 'pos':
        # rename this node property (position of the nodes) as 'pos'
        for u in G.nodes():
            G.nodes[u]['pos'] = G.nodes[u][pos_attr]
            del(G.nodes[u][pos_attr])

    Kg = KGraph(None, None, None, networkx_graph=G, verbose=verbose)
    
    if verbose:
        print("Karstnet graph successfully created from networkx graph !\n")

    return Kg
# ----------------------------------------------------------------------------

# ===== Import functions from csv file(s) ====================================
# Julien Straubhaar
# - function(s) calling the functions in karstnet/_io_func/text_files.py to import 
# the networkx graph, then KGraph object is instanciated

# ----------------------------------------------------------------------------
def from_csv(
        basename,
        suffix_nodes='_nodes.csv', 
        suffix_edges='_edges.csv', 
        delimiter_nodes=';',  
        delimiter_edges=';',
        import_node_scalar_properties=True,
        import_edge_scalar_properties=True,
        default_start_id=0,
        verbose=True):
    """
    Creates a karstnet graph from two csv files (nodes, and edges (links)).

    The file for nodes has the columns: 
        - ['id',] 'x', 'y'[, 'z', '<node_prop0>', '<node_prop1>', ...]

    Default id (if column 'id' is not given) are integer starting from 0.
    The columns with a name other than 'id', 'x', 'y', 'z', if any, are 
    additional scalar properties (attributes) attached to the nodes.
    The column 'z' can be omitted for 2D case.

    The file for edges (links) has the columns:
        - '<from_id>', '<to_id>'[, '<edge_prop0>', '<edge_prop1>', ...]

    where the first two columns, '<from_id>', '<to_id>', are the node ids 
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
    Kg : karstnet.KGraph
        a KGraph object

    Examples
    --------
        >>> Kg = kn.from_csv('my_graph')
    """

    G = kn.io_func.networkx_from_csv(
            basename,
            suffix_nodes=suffix_nodes,
            suffix_edges=suffix_edges, 
            delimiter_nodes=delimiter_nodes,  
            delimiter_edges=delimiter_edges,
            pos_attr='pos',
            import_node_scalar_properties=import_node_scalar_properties,
            import_edge_scalar_properties=import_edge_scalar_properties,
            default_start_id=default_start_id,
            verbose=False)

    Kg = KGraph(None, None, None, networkx_graph=G, verbose=verbose)
    del(G)

    if verbose:
        print("Karstnet graph successfully created from csv files !\n")

    return Kg
# ----------------------------------------------------------------------------

# ===== Import functions from yaml file ======================================
# - function calling the functions in karstnet/_io_func/yamlfile.py to import 
# the networkx graph, then KGraph object is instanciated

# ----------------------------------------------------------------------------
def from_yaml(filepath, add_edge_attr=True, add_node_attr=True, verbose=True):
    """
    Creates a karstnet graph from yaml file.
   
    Parameters
    ----------
    filepath : str
        path to the yaml file with the name and extension 

    add_edge_attr : str, optional
        Load and attach edge attributes, by default True.
        False: does not import any edge attributes.
        True: import all the existing edge attributes.
        ['attr1','attr2',...]: import only the attributes in the list.
        
    add_node_attr : str, optional
        Load and attach node attributes, by default True.
        False: does not import any node attributes.
        True: import all the existing node attributes.
        ['attr1','attr2',...]: import only the attributes in the list.
        'attr': import only one attrivute

    verbose : bool, default: True
        indicates if info is printed
    
    Returns
    -------
    Kg : karstnet.KGraph
        a KGraph object

    Examples
    --------
        >>> Kg = kn.from_yaml('my_graph')
    """

    G = kn.io_func.networkx_from_yaml(
            filepath,
            add_edge_attr=add_edge_attr,
            add_node_attr=add_node_attr)

    Kg = KGraph(None, None, None, networkx_graph=G, verbose=verbose)
    del(G)

    if verbose:
        print("Karstnet graph successfully created from yaml file !\n")

    return Kg
# ----------------------------------------------------------------------------

# ===== Import functions from sql therion ====================================
# - function(s) calling the functions in karstnet/_io_func/therion.py to import 
# the networkx graph, then KGraph object is instanciated

# ----------------------------------------------------------------------------
def from_therion_sql(
        basename,
        pos_attr='pos',
        remove_flagged_edges=False,
        verbose=True):
    """
    Creates a karstnet graph from SQL file exported from a Therion survey file.

    Parameters
    ----------
    basename : str
        the base name (prefix) used for the input file;
        the input file is named using the following convention:

        - `basename`.sql: the containing all the needed informations

    pos_attr : str, default: 'pos'
        name of the node attribute for position

    remove_flagged_edges : bool
        If set to True, all the edges flagged with 'dpl' or 'srf' are removed from the dataset
        by default: False

    verbose : bool, default: True
        indicates if info is printed

    Returns
    -------
    Kg : karstnet.KGraph
        a KGraph object

    Examples
    --------
        >>> myKGraph = kn.from_therion_sql("MyKarst")
    """

    G = kn.io_func.networkx_from_therion_sql(
            basename,
            pos_attr=pos_attr,
            remove_flagged_edges=remove_flagged_edges,
            verbose=False)

    Kg = KGraph(None, None, None, networkx_graph=G, verbose=verbose)
    del(G)

    if verbose:
        print("Karstnet graph successfully created from sql file!\n")

    return Kg
# ----------------------------------------------------------------------------

# # ----------------------------------------------------------------------------
# def from_therion_sql_enhanced(
#         inputfile, 
#         cavename=None, 
#         crs=None, 
#         rights=None,
#         citation=None,
#         verbose=True):
#     """
#     Creates a karstnet graph from SQL file exported from a Therion survey file.

#     See function `karstnet.io_func.networkx_from_therion_sql_enhanced`.

#     Parameters
#     ----------
#     inputfile : str
#         path to the SQL file exported with Therion

#     cavename : str, optional
#         Name of the cave. Will be attached to the graph as metadata. By default None

#     crs : str, optional
#         Coordinate reference. Will be attached to the graph as metadata. By default None

#     rights : str, optional
#         Information about the rights related to the dataset. Example: CC-BY-NC-SA, ... . Will be attached to the graph as metadata. By default None

#     citation : str, optional
#         Description on how to cite the dataset. Will be attached to the graph as metadata. By default None

#     verbose : bool, default: True
#         indicates if info is printed

#     Returns
#     -------
#     Kg : karstnet.KGraph
#         a KGraph object

#     Examples
#     --------
#         >>> myKGraph = kn.from_therion_sql_enhanced('inputfilepath.sql')
#     """

#     G = kn.io_func.networkx_from_therion_sql_enhanced(
#             inputfile, 
#             cavename=cavename, 
#             crs=crs, 
#             rights=rights,
#             citation=citation)

#     Kg = KGraph(None, None, None, networkx_graph=G, verbose=verbose)
#     del(G)

#     if verbose:
#         print("Karstnet graph successfully created from sql file (enhanced)!\n")

#     return Kg
# ----------------------------------------------------------------------------

# ===== Import functions from pline (GOCAD ASCII) ============================
# - function(s) calling the functions in karstnet/_io_func/therion.py to import 
# the networkx graph, then KGraph object is instanciated

# ----------------------------------------------------------------------------
def from_pline(
        filename, 
        verbose=True):
    """
    Creates a karstnet graph from a Pline (Gocad ascii object).

    Parameters
    ----------
    filename : str
        the name of the GOCAD Pline ASCII file

    verbose : bool, default: True
        indicates if info is printed

    Returns
    -------
    Kg : karstnet.KGraph
        a KGraph object

    Examples
    --------
        >>> myKGraph = kn.from_pline("MyKarst.pl")
    """

    G = kn.io_func.networkx_from_pline(
            filename,
            verbose=False)

    Kg = KGraph(None, None, None, networkx_graph=G, verbose=verbose)
    del(G)

    if verbose:
        print("Karstnet graph successfully created from (pline) file !\n")

    return Kg
# ----------------------------------------------------------------------------

# ===== Import functions other formats... ====================================
# * same principle... (if possible)


# ############################################################################


# ****************************************************************************
# Obsolete functions for import
# (kept for compatibility)
# ****************************************************************************

# Import from networkx graph
# ----------------------------------------------------------------------------
def from_nxGraph(nxGraph, coordinates, properties=None, verbose=True):
    """
    Creates a Karst graph from a Networkx graph.

    Takes a graph in the Networkx format and the coordinates of the
    nodes to generate a karstic network.

    Parameters
    ----------
    nxGraph : networkx graph
        the input graph

    coordinates : dictionary
        the coordinates of the node, keys are node names, values are
        sequence of 2 or 3 floats (for 2d or 3d case)

    properties : dictionary, optional
        optional argument containing properties associated with the nodes

    verbose : bool, default: True
        indicates if info is printed

    Returns
    -------
    KGraph
        A KGraph object

    Notes
    -----
    Obsolete function, used `from_networkx` instead'.

    Examples
    --------
        >>> myKGraph = kn.from_nxGraph(G, coord)
        >>> myKGraph = kn.from_nxGraph(G, coord, prop)
    """
    # # Initialization of the original graph
    # edges = nx.to_edgelist(nxGraph)
    # Kg = KGraph(edges, coordinates, properties=properties, verbose=verbose)

    print('Note: `from_nxGraph` is obsolete, use `from_networkx` instead')

    nx.set_node_attributes(nxGraph, coordinates, 'pos')
    if properties is not None:
        nx.set_node_attributes(nxGraph, properties, 'additional_property')
        
    Kg = from_networkx(nxGraph, verbose=verbose)

    return Kg
# ----------------------------------------------------------------------------

# Import from ascii files
# ----------------------------------------------------------------------------
def from_nodlink_dat(
        basename,
        suffix_nodes='_nodes.dat', 
        suffix_links='_links.dat',
        delimiter_nodes=' ',  
        delimiter_links=' ',
        verbose=True):
    """
    Creates the Kgraph from two ascii files (nodes, and links).

    The file for nodes contains the node coordinates and optional 
    properties: one line per node, the first two columns are the 
    x- and y- coordinates, the third column (optional) is the
    z-coordinate and the next columns (if any) are the properties

    The file for links contains the list of edges:
    one edge per line: the node ids of the two extremities (2 columns);
    the id corresponds to the line number in the node file.

    Note: the final node ids (node names) starts from 0 in the output
    object.

    Parameters
    ----------
    basename : str
        the base name (prefix) used for the input files

    suffix_nodes : str, default: '_nodes.dat'
        the input file containing the nodes is 
        `basename``suffix_nodes`
    
    suffix_links : str, default: '_links.dat'
        the input file containing the links (edges) is 
        `basename``suffix_links`
    
    delimiter_nodes : str, default:' '
        delimiter used in the file for nodes
    
    delimiter_links : str, default:' '
        delimiter used in the file for links
 
    verbose : bool, default: True
        indicates if info is printed

    Returns
    -------
    KGraph
        A KGraph object

    Notes
    -----
    Obsolete function, used `from_networkx` instead'.

    Examples
    --------
        >>> myKGraph = kn.from_nodlink_dat("MyKarst")
    """
    print('Note: `from_nodlink_dat` is obsolete, use `from_csv` instead')

    # Files
    filename_links = f'{basename}{suffix_links}'
    filename_nodes = f'{basename}{suffix_nodes}'

    # Read data files if exist - otherwise return empty graph
    try:
        links = np.loadtxt(filename_links, delimiter=delimiter_links).astype(int) - 1
    except OSError:
        print("IMPORT ERROR: Could not import {}".format(filename_links))
        return

    try:
        nodes = np.loadtxt(filename_nodes, delimiter=delimiter_nodes)
    except OSError:
        print("IMPORT ERROR: Could not import {}".format(filename_nodes))
        return
    # Create the dictionary of coordinates
    coord = dict(enumerate(nodes[:, :3].tolist()))

    if len(nodes[0] > 3):
        node_property = dict(enumerate(nodes[:, 3:].tolist()))
        node_property_name = 'additional_properties'
    else:
        node_property = None
        node_property_name = None

    Kg = KGraph(links, coord, node_property=node_property, node_property_name=node_property_name, verbose=verbose)

    if verbose:
        print("Karstnet graph successfully created from ascii files !\n")

    return Kg
# ----------------------------------------------------------------------------
