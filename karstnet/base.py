#    Copyright (C) 2018-2023 by
#    Philippe Renard <philippe.renard@unine.ch>
#    Pauline Collon <pauline.collon@univ-lorraine.fr>
#    All rights reserved.
#    MIT license.
#
"""
Karstnet
========

Karstnet is a Python package for the analysis of karstic networks.

License
-------
Released under the MIT license:
   Copyright (C) 2018-2023 Karstnet Developers
   Philippe Renard <philippe.renard@unine.ch>
   Pauline Collon <pauline.collon@univ-lorraine.fr>
"""

# ----External libraries importations
import numpy as np
import networkx as nx
import scipy.stats as st
import matplotlib.pyplot as plt
import matplotlib as mtb
from typing import List, Dict, Tuple, Optional
# noinspection PyUnresolvedReferences
import mplstereonet


# mplstereonet is used in the notebook - can not be delete


# *************************************************************
# -------------------Test function--------------------------
# *************************************************************
def test_kn():
    print("test ok")
    print("recovery ok")


# **************************************************************
#
# -------------------KGraph class-------------------------------
#
# **************************************************************


class KGraph:
    """
    Class dedicated to the construction and manipulation of graphs
    representing Karstic network.

    Attributes:
    -----------
        - graph : the complete graph of the karstic network.
                Each station is a node, each line-of sight is an edge.
                It is a Networkx graph object, with length of edges
                as attributes.
        - graph_simpl: the simplified graph of the karstic network,
                i.e. all nodes of degree 2 are removed except in loops
                (2 types of loops)(cf Collon et al, 2017, Geomorphology)
                graph_simple is a Networkx graph object, with length of
                edges as attributes.
        - pos2d : a dictionary of the nodes with their 2d coordinates
                as a values [x, y]
        - pos3d : a dictionary of the nodes with their 3d coordinates
                as a values [x, y, z]
        - properties : a dictionary of nodes with their properties in
                case of importing karst data with additional information
        - branches : the list of branches
        - br_lengths : the list of branch lengths
        - br_tort : the list of branches tortuosities
        - list_simpl_edges : the list of simple edges, necessary to export
                graph to plines
        - graph_simpl : the simplified version (without nodes of degree 2)
                of a graph

    """

    def __init__(self, edges: List[int], coordinates: Dict[int, List[float]],
                 properties: Optional[Dict[int, List[float]]] = None,
                 verbose=True):
        """
        Creates a Kgraph from nodes and edges.

        Parameters
        ----------
            edges : list
                a list of edges

            coordinates : dictionary
                coordinates of the nodes, keys are node names

        properties : dictionnary
            optional properties associated to the nodes

        Examples
        --------
           >>> myKGraph = KGraph([],{})
        """

        # Initialization of the complete graph - use functions of kgraph_fc
        # module
        self.graph = nx.Graph()
        self.graph.add_edges_from(edges)
        self.verbose = verbose

        self.pos2d, self.pos3d, self.nodes_xz, self.nodes_yz = _pos_initialization(coordinates)
        if self.verbose:
            print(
                "\n This network contains ",
                nx.number_connected_components(
                    self.graph),
                " connected components")

        # Compute graph length from the pos3d_ to initialize properly the graph
        self._set_graph_lengths()
        self._set_graph_orientations()
        if properties is None:
            self.properties = dict()
        else:
            self.properties = properties
        self._set_graph_properties()

        # Compute branches of the graph
        # self.branches is necessary to export graph to plines
        self.branches, self.br_lengths, self.br_tort = self._getallbranches()

        # Construct the simplified graph
        # self.list_simpl_edges is necessary to export graph to plines
        self.list_simpl_edges, self.graph_simpl = self._simplify_graph()

    # **********************************
    #    Plots
    # **********************************

    def plot_properties2(self, prop_type: Optional[str] = None,
                         figsize: Tuple[float] = (6, 3)):
        """
        Plot a 2D view of the karstic network with properties.
        Called by the plot2 public function

        Parameters
        ----------
            prop_type : string
                plot nodes that get properties data in green "chartreuse"

            figsize : tuple
                contains the (x,y) dimension of the figure

        Examples
        --------
            >>> my_kgraph = KGraph([],{},{})
            >>> my_kgraph.plot_properties2("radius")
        """

        min_v, max_v = self._plot_properties2(self.graph, prop_type, figsize)
        _plot_legend(min_v, max_v)
        plt.title(prop_type)
        plt.xlabel('x')
        plt.ylabel('y')
        plt.show()

    def plot2(self, graph_type: int = 0, prop_type: Optional[str] = None,
              figsize: Tuple[float] = (6, 3)):
        """
        Plot a 2D view of the karstic network

        Parameters
        ----------
            graph_type : int
                if 0 displays the complete graph,
                if 1 displays the simplified graph

            prop_type : string
                if None there is no properties
                if it is a string, there are properties to describe the network

            figsize : tuple
                contains the (x,y) dimension of the figure

        Examples
        --------
           >>> my_kgraph = KGraph([],{})
           >>> my_kgraph.plot2()
           >>> my_kgraph.plot2(1)
           >>> my_kgraph.plot2(graph_type=1, prop_type="radius")
        """
        if prop_type is not None:
            self.plot_properties2(prop_type, figsize)

        else:
            if graph_type == 0:
                self._plot2(self.graph, figsize)
                plt.title('original')

            else:
                self._plot2(self.graph_simpl, figsize)
                plt.title('simplified')

            plt.xlabel('x')
            plt.ylabel('y')
            plt.show()

        return

    def plot_properties3(self, prop_type: Optional[str] = None,
                         zrotation: float = 30, xyrotation: float = 0,
                         figsize: Tuple[float] = (4, 3)):
        """
        Plot a 3D view of the karstic network with properties.
        Called by the plot3 public function

        Parameters
        ----------
            prop_type : string

            zrotation : float
                angle in degrees between horizontal plane and viewpoint

            xyrotation : float
                angle in degree for the horizontal rotation of the viewpoint
                If xyrotation=0, the view is from the South toward North.

            figsize : tuple
                contains the (x,y) dimension of the figure

        Examples
        --------
            >>> my_kgraph = KGraph([],{},{})
            >>> my_kgraph.plot_properties3("radius")
        """
        min_v, max_v = self._plot_properties3(self.graph, prop_type, zrotation,
                                              xyrotation, figsize)
        _plot_legend(min_v, max_v)
        plt.title(prop_type)
        plt.show()

    def plot3(self, graph_type: int = 0, prop_type: Optional[str] = None,
              zrotation: float = 30, xyrotation: float = 0,
              figsize: Tuple[float] = (4, 3)):
        """
        Plot a 3D view of the karstic network.

        Parameters
        ----------
            graph_type : int
                if 0 displays the complete graph,
                if 1 displays the simplified graph

            prop_type : string
                if None there is no properties
                if it is a string, there are properties to describe the network

            zrotation : float
                angle in degrees between horizontal plane and viewpoint

            xyrotation : float
                angle in degree for the horizontal rotation of the viewpoint
                If xyrotation=0, the view is from the South toward North.

            figsize : tuple
                contains the (x,y) dimension of the figure

        Examples
        --------
            >>> my_kgraph = KGraph([],{})
            >>> my_kgraph.plot3()
            >>> my_kgraph.plot3(1, zrotation=20, xyrotation=-30)
            >>> my_kgraph.plot3(prop_type="radius")
        """
        # If there are one or more properties
        if prop_type is not None:
            self.plot_properties3(prop_type, zrotation, xyrotation, figsize)

        # If there are any properties
        else:
            if graph_type == 0:
                self._plot3(self.graph, zrotation, xyrotation, figsize)
                plt.title('original')

            else:
                self._plot3(self.graph_simpl, zrotation, xyrotation, figsize)
                plt.title('simplified')

            plt.show()
        return

    def plot(self):
        """
        Simple 2D map of the original and simplified karstic network.

        The two maps are ploted side by side. This function allows
        to check rapidly the data after an import for example.

        Examples
        ---------
           >>> my_kgraph = KGraph([],{})
           >>> my_kgraph.plot()
        """
        plt.figure(figsize=(12, 5))
        plt.subplot(121)
        nx.draw_networkx(self.graph,
                         pos=self.pos2d,
                         with_labels=False,
                         node_size=0.1)
        plt.xlabel('x')
        plt.ylabel('y')
        plt.subplot(122)
        nx.draw_networkx(self.graph_simpl,
                         pos=self.pos2d,
                         with_labels=False,
                         node_size=0.1)
        plt.xlabel('x')
        plt.ylabel('y')
        plt.show()

    def plotxz(self):
        """
        Simple 2D map of the original and simplified karstic network,
         slice view of the karstic network.

        The two maps are ploted side by side. This function allows
        to check rapidly the data after an import for example.

        Examples
        ---------
           >>> my_kgraph = KGraph([],{})
           >>> my_kgraph.plot()
        """
        plt.figure(figsize=(12, 5))
        plt.subplot(121)
        nx.draw_networkx(self.graph,
                         pos=self.nodes_xz,
                         with_labels=False,
                         node_size=0.1)
        plt.xlabel('x')
        plt.ylabel('z')
        plt.subplot(122)
        nx.draw_networkx(self.graph_simpl,
                         pos=self.nodes_xz,
                         with_labels=False,
                         node_size=0.1)
        plt.xlabel('x')
        plt.ylabel('z')
        plt.show()

    # Member function written by Philippe Vernant 2019/11/25
    # Modified by Pauline Collon (aug. 2020) to weight density map by length

    def stereo(self, weighted: bool = True):
        """
        Density map of orientations and rose diagram of the karstic network.
        The two stereo are ploted side by side.
        By default, stereo and Rose diagram are weighted by lengths.

        Parameters
        ----------
            weighted : boolean
                By default, stereo and Rose diagram are weighted by lengths

        Examples
        ---------
           >>> my_kgraph = KGraph([],{})
           >>> my_kgraph.stereo()
           >>> my_kgraph.stereo(weighted = False)
        """

        # create an np.array of azimuths and dips
        # and lengths (projected(2d) and real (3d))
        azim = np.array(
            list((nx.get_edge_attributes(self.graph, 'azimuth')).values()))
        azim_not_nan = azim[~np.isnan(azim)]
        bearing_dc = np.nan_to_num(azim)
        plunge_dc = np.array(
            list((nx.get_edge_attributes(self.graph, 'dip')).values()))
        if (weighted):
            l2d = np.array(
                list((nx.get_edge_attributes(self.graph,
                                             'length2d')).values()))

            l3d = np.array(
                list((nx.get_edge_attributes(self.graph, 'length')).values()))

            l2d_not_nan = l2d[~np.isnan(azim)]
        else:
            l2d_not_nan = None
            l3d = None
        # Pauline: not sure it is required (?)
        # + not sure it is normal that isnan is parameterised by azim for l2d

        # import matplotlib as mpl

        # Making colormap, based on Collon et al.(2017) \
        # we saturate the colormap at 40%
        from matplotlib import cm
        from matplotlib.colors import ListedColormap
        nbint = 15
        levels = np.linspace(0, 1, nbint)
        rainbow = cm.get_cmap('rainbow')
        newcolors = rainbow(levels)
        white = np.array([256 / 256, 256 / 256, 256 / 256, 1])
        newcolors[:1, :] = white
        newcmp = ListedColormap(newcolors)

        # Density map - Allows to consider almost vertical conduits
        # The data are weighted by the real length of the segments (l3d)
        # Use the traditional "Schmidt" method : 1% count
        fig = plt.figure(figsize=(16, 8))
        dc = fig.add_subplot(121, projection='stereonet')
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

        # colorbar of the density map
        cbar = plt.colorbar(cdc,
                            fraction=0.046,
                            pad=0.04,
                            orientation='horizontal')
        cbar.set_label('[%]')

        # Rose diagram
        # The azimuth data are weighted by the projected length (l2d)
        bin_edges = np.arange(-5, 366, 10)
        number_of_strikes, bin_edges = np.histogram(azim_not_nan,
                                                    bin_edges,
                                                    weights=l2d_not_nan)
        number_of_strikes[0] += number_of_strikes[-1]
        half = np.sum(np.split(number_of_strikes[:-1], 2), 0)
        two_halves = np.concatenate([half, half])

        rs = fig.add_subplot(122, projection='polar')
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

        fig.tight_layout()
        plt.show()

    # end modified PV 2019/11/25

    def variogram_1D(self, variable: str, round_h: int = 2,
                     tolerance_h: int = 0.01):
        """
        Plot the variogram of a properties

        Parameters
        ----------
            variable : string
                Variable of which one whishes to make the 1D variogram

            round_h : int
                Number of digits after the decimal point desired for
                calculating the length of the paths

            tolerance_h : int
                Number of digits after the decimal point desiredfor the values
                of h

        Examples
        --------
            >>> my_kgraph = KGraph([],{},{})
            >>> my_kgraph.variogram_1D("radius")
           """
        if nx.get_node_attributes(self.graph, variable):
            dict_prop = nx.get_node_attributes(self.graph, variable)

            # In case of no lake of data for nodes, we remove the node
            # of dictionary
            length = nx.get_edge_attributes(self.graph, "length")
            shortest_path = nx.shortest_path(self.graph)
            shortest_path_del = _del_node(shortest_path, dict_prop)
            shortest_path_conv = _convert(shortest_path_del)
            shortest_path_h = _calc_length(shortest_path_conv, length, round_h)
            list_h = self._list_h(shortest_path_h, tolerance_h)
            shortest_path_inv = self._inverse(dict_prop, list_h,
                                              shortest_path_h, tolerance_h)

            gamma = self._calc_gamma(list_h, shortest_path_inv)

            plt.plot(list_h, gamma, "o")
            plt.title("1D variogram ")
            plt.xlabel("h")
            plt.ylabel("gamma(h)")
            plt.show()

        else:
            print("No data properties")
        return

    # *************************************************************
    # ----------------------- Export ------------------------------
    # *************************************************************

    def to_pline(self, basename: str):
        """
        Export the complete graph to Pline (GOCAD ASCII object)

        Manages the colocated vertices indicated by the mention "ATOM" in
        the ASCII file.

        `Warning`: this version does not export (for the time) the
        properties on the nodesself.

        Parameters
        ----------
            basename : string
                the base name of the file name used for the Pline file,
                the name contains no extension, it will be added by the
                function

        Examples
        --------
        The following command saves the file "MyKarst_exported.pl" :

            >>> my_kgraph = KGraph([],{})
            >>> my_kgraph.to_pline("MyKarst")
        """
        # For a complete graph the list of Ilines corresponds to self.branches

        self._ilines_to_pline(self.branches, basename)

        return

    def simpleGraph_to_pline(self, basename: str):
        """
        Export the simplified graph to pline (GOCAD ASCII object)

        Arguments:
        ----------
            basename: A string containing the base name used for output file.
            Manage the colocated vertices indicated by the mention "ATOM"
            in the ASCII file
            Warning: this version does not export (for the time) the
            properties on the nodes

        Parameters
        ----------
            basename : string
                the base name of the file name used for the Pline file,
                the name contains no extension, it will be added by the
                function

        Examples
        --------
        The following command saves the file "MyKarst_simpl_exported.pl" :

            >>> my_kgraph = KGraph([],{})
            >>> my_kgraph.to_pline("MyKarst")
        """
        # For a simplified graph the list of Ilines will corresponds to
        # the edges of the simple graph
        # one iline is created for each edge, which is not exactly a branch
        # but does not prevent exportation

        # to clearly explicit it is the simplified graph
        basename = basename + "_simpl_"
        self._ilines_to_pline(self.list_simpl_edges, basename)

        return

    def to_dat(self, basename: str):
        """
        Export the complete graph to dat file

        `Warning`: this version does not export (for the time) the
        properties on the nodesself.

        Parameters
        ----------
            basename : string
                the base name of the file name used for the dat file,
                the name contains no extension, it will be added by the
                function

        Examples
        --------
        The following command saves the files :
            basename_exported_links.dat
            basename_exported_nodes.dat

         Example:
         --------
            >>> my_kgraph = KGraph([],{})
            >>> my_kgraph.to_dat("MyKarst" )

        """
        # For a complete graph the list of Ilines corresponds to self.branches
        print('basename : ', basename)
        print(type(basename))
        self._ilines_to_dat(self.branches, basename)

        return

    # *************************************************************
    # -------------------Computation for Analysis------------------
    # *************************************************************
    def basic_analysis(self):
        """
        Print the basics of a karstic network graph analysis

        Examples
        --------
            >>> my_kgraph = KGraph([],{})
            >>> my_kgraph.basic_analysis()
        """

        # On the complete graph
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
            if (self.graph_simpl.degree(i) == 1):
                nb_extremity_nodes += 1
            elif (self.graph_simpl.degree(i) > 2):
                nb_junction_nodes += 1

        # Print these basics
        print(
            "\n This network contains :\n",
            nb_nodes_comp,
            " nodes (stations) and ",
            nb_edges_comp,
            " edges.\n",
            " On the simplified graph, there are : ",
            nb_nodes,
            " nodes (stations) and ",
            nb_edges,
            " edges,\n",
            nb_extremity_nodes,
            " are extremity nodes (entries or exits) and ",
            nb_junction_nodes,
            " are junction nodes.\nThere is/are ",
            nb_connected_components,
            " connected component.s and ",
            nb_cycles,
            " cycle.s.\n")

        # Howard's parameters
        # (Howard, A. D., Keetch, M. E., & Vincent, C. L. (1970).
        # Topological and geometrical properties of braided patterns.
        # Water Resources Research, 6(6), 1674–1688.)
        # Rmq: All nodes of the simplified network have to be considered,
        # even those of degree 2 in the cycles or, it is not becoming
        # consistent with the examples of braided rivers given by Howard.
        # This is indeed consistent with what has been done in
        # Collon, P., Bernasconi, D., Vuilleumier, C., & Renard, P. (2017).
        # Statistical metrics for the characterization of karst network
        # geometry and topology. Geomorphology, 283, 122–142.
        alpha = nb_cycles / (2 * (nb_nodes) - 5)
        beta = nb_edges / (nb_nodes)
        gamma = nb_edges / (3 * (nb_nodes - 2))
        print("\nHoward's parameter are (Howard, 1970) :",
              " \n alpha: ", alpha,
              "\n beta", beta,
              "\n gamma", gamma)
        print("\nNote that this computation considers the node of degree 2",
              " necessary to loop preservations as Seed Nodes, in order to",
              " stay consistent with Howard's illustrations.")

    def mean_tortuosity(self):
        """
        Compute the mean tortuosity of a karstic network

        Returns
        -------
        float
            Mean tortuosity of the branches

        Examples
        --------
           >>> my_kgraph = KGraph([],{})
           >>> t = my_kgraph.mean_tortuosity()
        """
        nb_of_nan = np.isnan(self.br_tort).sum()
        
        if self.verbose:
            if nb_of_nan != 0:
                print(
                    "\n WARNING: This network contains ",
                    nb_of_nan,
                    " looping branche.s, which is.are not considered for the ",
                    "mean tortuosity computation")
        return (np.nanmean(self.br_tort))

    def mean_length(self):
        """
        Compute the mean length of the branches of a karstic networkx

        Returns
        -------
            l :  Mean length of the branches

        Example:
        --------
           >>> my_kgraph = KGraph([],{})
           >>> l = my_kgraph.mean_length()
        """
        return (np.nanmean(self.br_lengths))

    def coef_variation_length(self) -> float:
        """
        Compute the coefficient of variation of length of the branches of a
        karstic networkx

        Returns
        -------
        float
            cvl :  Coefficient of variation of the length of the branches

        Example:
        --------
           >>> my_kgraph = KGraph([],{})
           >>> cvl = my_kgraph.coef_variation_length()
        """

        # std is used with ddof =1 to use the estimate of std (divides by N-1)
        return (np.nanstd(self.br_lengths, ddof=1) /
                np.nanmean(self.br_lengths))

    def length_entropy(self, mode="default") -> float:
        """
        Compute the entropy of lengths of the branches of a karstic network

        Parameters
        ----------
        mode : string
            the mode style should be equal to "default" or "sturges"
            The property is normalized as compared to the
            maximum observed value
            If mode= "default", the entropy is computed according to the
            implementation proposed in Collon et al. 2017 :
            using a fixed bin number of 10 and ranging from 0 to 100
            If mode = "Sturges" : alternative calculation using Sturges
            rules to define the number of bins

        Returns
        -------
        float
            entropy :  Entropy of the lengths of the branches

        Example:
        --------
           >>> my_kgraph = KGraph([],{})
           >>> l_entrop = my_kgraph.length_entropy()
        """

        v = self.br_lengths
        # In the paper of 2017, we normalize the length to get comparble
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
                                     range=(np.nanmin(v) * 0.97,
                                            np.nanmax(v) * 1.08)) #modif by Louise to be verified
                                    # previous version :   range=(0, 100))
            freq = counts / sum(counts)  # Computes the frequencies
            entropy = st.entropy(freq, base=nbins)
        else:
            entropy = 0  # v contains a single value - no uncertainty

        return entropy

    def orientation_entropy(self, mode="default") -> float:
        """
        Computes the entropy of orientation of the segments of a
        karstic network.

        Parameters
        ----------
        mode : string
            the mode style should be equal to "default" or "sturges"
            If mode= "default", the entropy is computed according to
            the implementation proposed in Collon et al. 2017 :
            using a fixed bin number of 18 and ranging from 0 to 180
            If mode = "Sturges" : alternative calculation using Sturges
            rules to define the number of bins

        Returns
        -------
        float
            entropy:  Entropy of the segment orientation

        Example:
        --------
           >>> my_kgraph = KGraph([],{})
           >>> or_entropy = my_kgraph.orientation_entropy()
        """

        # create an np.array of azimuths and projected lengths
        azim = np.array(
            list((nx.get_edge_attributes(self.graph, 'azimuth')).values()))
        l2d = np.array(
            list((nx.get_edge_attributes(self.graph, 'length2d')).values()))

        # Removing NAN Azimuth values that correspond to length2d=0
        azim_not_Nan = azim[~np.isnan(azim)]
        l2d_not_zero = l2d[np.nonzero(l2d)]

        if len(azim_not_Nan) > 1:

            if mode == "sturges":
                # Sturges rule to define the number of bins fron nb of samples
                nbins = int(np.ceil(1 + np.log2(len(azim_not_Nan))))
            else:
                # Fixed nb of bins to facilitate comparison,
                # as done in Collon et al 2017 and other papers on roads
                # orientation entropy
                nbins = 18

            # Pauline : range should be 0 and 180 or we get very different
            # results from what we had (change the position of bin borders).
            # Also, I verified that it is consistent :
            # 0 and 180 are counted, not excluded
            counts, _ = np.histogram(azim_not_Nan,
                                     bins=nbins,
                                     range=(0, 180),
                                     weights=l2d_not_zero)
            freq = counts / sum(counts)  # Computes the frequencies
            return st.entropy(freq, base=nbins)
        else:
            return 0

    def mean_degree_and_cv(self) -> tuple:
        """
        Computes the average and the coefficient of variation of the degree.

        The computation is done on the simplified graph.

        Returns
        -------
        tuple
            meandeg, cvdeg :  the mean and coefficient of variation

        Examples
        --------
           >>> my_kgraph = KGraph([],{})
           >>> mean_deg, cv_deg = my_kgraph.mean_degree_and_cv()
        """
        # Vector of degrees
        d = np.array(self.graph_simpl.degree())[:, 1]

        # Mean degree
        meandeg = np.nanmean(d)

        # Coefficient of variation of the degrees : std is used with ddof =1 to
        # use the estimate of std (divides by N-1)
        cvde = np.nanstd(d, ddof=1) / np.nanmean(d)

        return meandeg, cvde

    def correlation_vertex_degree(self, cvde: Optional[float] = False) -> float:
        """
        Computes the correlation of vertex degree.

        The computation is done on the simplified graph.

        Parameters
        ----------
            cvde : float
                Optional input: coefficient of variation of the degree.
                If not provided, it is computed automatically internally.

        Returns
        -------
        float
            Correlation of Vertex Degree

        Examples
        --------
            >>> my_kgraph = KGraph([],{})
            >>> _cvd = my_kgraph.correlation_vertex_degree()
        """
        if not cvde:
            _, cvde = self.mean_degree_and_cv()

        # To avoid division by 0 when computing correlation coef
        if cvde != 0:
            cvd = nx.degree_pearson_correlation_coefficient(self.graph_simpl)
        else:
            cvd = 1

        return cvd

    def central_point_dominance(self) -> float:
        """
        Computes central point dominance.

        The computation is done on the simplified graph.

        Returns
        -------
        float
            Central point dominance

        Examples
        --------
            >>> my_kgraph = KGraph([],{})
            >>> cpd = my_kgraph.central_point_dominance()
        """
        bet_cen = nx.betweenness_centrality(self.graph_simpl)
        bet_cen = list(bet_cen.values())
        return sum(max(bet_cen) - np.array(bet_cen)) / (len(bet_cen) - 1)

    def average_spl(self, dist_weight: Optional[float] = False) -> float:
        """
        Computes average shortest path lengths.

        The computation is done on the simplified graph.
        The function handles the case of several connected components
        which is not the case for the Networkx function
        "average_shortest_path_length".
        In case of several connected components, the average_SPL
        is the average of each SPL weighted by the number of nodes of each
        connected component

        Returns
        -------
        float
            average shortest path lengths

        Examples
        --------
            >>> my_kgraph = KGraph([],{})
            >>> aspl = my_kgraph.average_spl()
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

    def mean_std_radius(self) -> Optional[tuple]:
        """
        Returns
        -------
        tuple
            The mean radius and standard deviation value.

        Examples
        --------
            >>> my_kgraph = KGraph([],{},{})
            >>> _mean,_std = my_kgraph.mean_std_radius()
        """
        if nx.get_node_attributes(self.graph, "radius"):
            radius = np.array(
                list((nx.get_node_attributes(self.graph, 'radius')).values()))
            mean = np.nanmean(radius)
            std = np.nanstd(radius)

            return mean, std

        print("There is no radius data")
        return None

    def mean_std_wh_ratio(self) -> Optional[tuple]:
        # attention au cas où il n'y a pas de wh_ratio
        """
        Returns
        -------
        tuple
            The mean wh_ratio and standart deviation value.

        Examples
        --------
            >>> my_kgraph = KGraph([],{})
            >>> _mean, _std = my_kgraph.mean_std_wh_ratio()
        """

        if nx.get_node_attributes(self.graph, "wh_ratio"):
            wh_ratio = np.array(
                list((nx.get_node_attributes(self.graph, 'wh_ratio')).values()))

            mean = np.nanmean(wh_ratio)
            std = np.nanstd(wh_ratio)

            return mean, std

        print("There is no wh_ratio data")
        return None

    def characterize_graph(self, verbose: bool = False) -> dict:
        """
        Computes the set of metrics used to characterize a graph.

        Parameters
        ----------
            verbose : boolean
                If True, the function displays information about the
                progress of the computation, and the results.

        Returns
        -------
        dictionary
            All the statistical metrics are stored in a dictionnary. The
            keys of the dictionnary are provided below with
            corresponding explanation.

            `mean length` : mean length of the branches,
            `cv length` : coefficient of variation of length of branches,
            `length entropy` : entropy of the length of the branches,
            `mean tortuosity` : mean tortuosity of the branches,
            `orientation entropy` : entropy of the orientation of the conduits,
            `aspl` : average shortest path length,
            `cpd` : central point dominance,
            `mean degree` : mean of the vertex degrees,
            `cv degrees` : coefficient of variation of vertex degrees,
            `correlation vertex degree` :  correlation of vertex degrees,

        Examples
        --------
            >>> my_kgraph = KGraph([],{})
            >>> results = my_kgraph.characterize_graph()
        """

        results = {}

        if verbose:
            print('Computing:')
            print(' - mean length', end='', flush=True)

        results["mean length"] = self.mean_length()

        if verbose:
            print(',cv length', end='', flush=True)

        results["cv length"] = self.coef_variation_length()

        if verbose:
            print(',length entropy', end='', flush=True)

        results["length entropy"] = self.length_entropy()

        if verbose:
            print(',mean tortuosity', end='', flush=True)

        results["tortuosity"] = self.mean_tortuosity()

        if verbose:
            print('', end='\n', flush=True)
            print(' - orientation entropy', end='', flush=True)

        results["orientation entropy"] = self.orientation_entropy()

        if verbose:
            print(',aspl', end='', flush=True)

        results["aspl"] = self.average_spl()

        if verbose:
            print(',cpd', end='', flush=True)

        results["cpd"] = self.central_point_dominance()

        if verbose:
            print(',md,cv degree', end='', flush=True)

        md, cvde = self.mean_degree_and_cv()
        results["mean degree"] = md
        results["cv degree"] = cvde
        
        if verbose:
            print(',cvd', end='', flush=True)

        cvd = self.correlation_vertex_degree(cvde=cvde)
        results["correlation vertex degree"] = cvd

        if verbose:
            print(', mean radius, std radius', end='', flush=True)

        if self.mean_std_wh_ratio() is not None:
            mean_radius, std_radius = self.mean_std_radius()
            results["mean radius"] = mean_radius
            results["std radius"] = std_radius
        else:
            print(', There is no radius data', end='', flush=True)


        if verbose:
            print(',mean wh_ratio, std wh_ratio', end='', flush=True)

        if self.mean_std_wh_ratio() is not None:
            mean_wh_ratio, std_wh_ratio = self.mean_std_wh_ratio()
            results["mean wh_ratio"] = mean_wh_ratio
            results["std wh_ratio"] = std_wh_ratio
        else:
            print(', There is no wh_ratio data', end='', flush=True)

        if verbose:
                print('', end='\n', flush=True)
                print("--------------------------------------")
                for key in results.keys():
                    print(" %25s = %5.3f" % (key, results[key]))
                print("--------------------------------------")

        return (results)

    # *************************************************************************
    # Non Public member functions of KGraph class
    # *************************************************************************

    # *******************************
    # Private functions for variogram
    # *******************************

    @staticmethod
    def _list_h(shortest_path: Dict[Tuple[int], int],
                tolerance_h: int) -> List[float]:
        """
        NOT PUBLIC

        Parameters
        ----------
            shortest_path : dict
                a dictionary of the shortest path between each nodes of the
                network with the length of the calculated paths

            tolerance_h : int
                Number of digits after the decimal point desired for the values
                of h

        Returns
        -------
        list
            list of path lengths without repetition

        """
        list_h = list(shortest_path.values())
        list_h.sort()
        list_h_tri = []
        for i_h in range(len(list_h) - 1):
            if not (list_h[i_h] - tolerance_h <= list_h[i_h + 1] <=
                    list_h[i_h] + tolerance_h):
                list_h_tri.append(list_h[i_h + 1])
                if list_h[i_h] not in list_h_tri:
                    list_h_tri.append(list_h[i_h])
        return list_h_tri

    @staticmethod
    def _inverse(dict_prop: Dict[int, float], list_h: List[float],
                 shortest_path: Dict[Tuple[int], int], tolerance_h: int
                 ) -> Dict[float, List[Tuple[float, float]]]:
        # optimisable ?
        """
        NOT PUBLIC

        Parameters
        ----------
            dict_prop : dict
               Dictionary of the proprieties like {nodes: radius_properties}

            list_h : list
                List of path lengths without repetition

            shortest_path : dict
                Dictionary of the shortest path between each nodes of the
                network with the length of the calculated paths

            tolerance_h : int
                Number of digits after the decimal point desired for the values
                of h

        Returns
        -------
        dict
            Dictionary like {length: list of source node - target node
            associated with this length
        """
        dict_gamma = {}
        for h in list_h:
            list_node = []
            for i_h in range(len(shortest_path.values())):
                if (h - tolerance_h <= list(shortest_path.values())[i_h] <=
                        h + tolerance_h):
                    list_node.append(
                        (dict_prop[list(shortest_path.keys())[i_h][0]],
                         dict_prop[list(shortest_path.keys())[i_h][1]]))
            dict_gamma[h] = list_node

        return dict_gamma

    @staticmethod
    def _calc_gamma(list_h: List[float],
                    dict_length: Dict[float, List[Tuple[float, float]]]
                    ) -> List[float]:
        """
        NOT PUBLIC

        Parameters
        ----------
            list_h : list
                List of path lengths without repetition

            dict_length : dict
                Dictionary like {length: list of source node - target node
                associated with this length, returned by _inverse

        Returns
        -------
        list
            List of gamma values, in the same order as in list_h
        """
        gamma = []
        for h in list_h:
            n = len(dict_length[h])
            diff = 0
            for alpha in range(n):
                diff += (dict_length[h][alpha][0] -
                         dict_length[h][alpha][1]) ** 2
            gamma.append((1 / (2 * n)) * diff)
        return gamma

    # *******************************
    # Private functions for plots
    # *******************************

    def _plot_properties2(self, graph: nx.Graph,
                          prop_type: Optional[str] = None,
                          figsize: Tuple[float] = (4, 3)
                          ) -> Tuple[float, float]:
        """
        NOT PUBLIC

        Plot a 2D view of a graph G with colored node according to the value of
        the property.
        Requires self.pos2d => member function.
        Called by the plot_properties2 public function

        Parameters
        ----------
            graph : nx.Graph
                A networkx graph object

            prop_type : string
                if None there is no properties
                if it is a string, there are properties to describe the network

            figsize : tuple
                contains the (x,y) dimension of the figure

        Returns
        -------
        tuple
            The value minimum and maximum of the properties
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            raise ImportError("karstnet.plot3 requires matplotlib.pyplot")

        if prop_type is None:
            print("There is no properties data")
            self._plot2(graph, figsize)

        plt.figure(figsize=figsize)
        nx.draw_networkx(graph, with_labels=True, pos=self.pos2d,
                         node_color='lightblue')

        dict_prop = nx.get_node_attributes(graph, prop_type)
        val = list(dict_prop.values())
        max_v = np.nanmax(val)
        min_v = np.nanmin(val)
        diff = max_v - min_v

        list1, list2, list3, list4, nodata = [], [], [], [], []

        for k, v in dict_prop.items():
            if np.isnan(v):
                nodata.append(k)
            if v <= min_v + diff / 2:
                if v <= min_v + diff / 4:
                    list1.append(k)
                else:
                    list2.append(k)
            elif v >= max_v - diff / 4:
                list4.append(k)
            else:
                list3.append(k)

        nx.draw_networkx_nodes(graph, pos=self.pos2d,
                               nodelist=list1,
                               node_color='chartreuse')
        nx.draw_networkx_nodes(graph,
                               pos=self.pos2d,
                               nodelist=list2,
                               node_color='yellow')
        nx.draw_networkx_nodes(graph,
                               pos=self.pos2d,
                               nodelist=list3,
                               node_color='orange')
        nx.draw_networkx_nodes(graph,
                               pos=self.pos2d,
                               nodelist=list4,
                               node_color='red')
        nx.draw_networkx_nodes(graph,
                               pos=self.pos2d,
                               nodelist=nodata,
                               node_color='lightblue')

        # noinspection PyTypeChecker
        return min_v, max_v

    def _plot2(self, graph: nx.Graph, figsize: Tuple[float] = (4, 3)):
        """
        NOT PUBLIC

        Plot a 2D view of a graph G that could be the simplified of the
        complete one.
        Requires self.pos3d => member function.
        Called by the plot2 public function

        Parameters
        ----------
            graph : nx.Graph
                A networkx graph object

            figsize : tuple
                contains the (x,y) dimension of the figure
        """

        # 2D  plot
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            raise ImportError("karstnet.plot3 requires matplotlib.pyplot")

        plt.figure(figsize=figsize)

        nx.draw_networkx(graph,
                         with_labels=True,
                         pos=self.pos2d,
                         node_color='lightblue')

        return

    def _plot_properties3(self, graph: nx.Graph,
                          prop_type: Optional[str] = None,
                          zrotation: float = 30,
                          xyrotation: float = 10,
                          figsize: Tuple[float] = (4, 3)
                          ) -> Tuple[float, float]:
        """
        NOT PUBLIC

        Plot a 3D view of a graph with colored node according to the value of
        the properties.
        Requires self.pos2d => member function.
        Called by the plot_properties3 public function

        Parameters
        ----------
            graph : nx.Graph
                A networkx graph object

            prop_type : string
                if None there is no properties
                if it is a string, there are properties to describe the network

            zrotation : float
                angle in degrees between horizontal plane and viewpoint

            xyrotation : float
                angle in degree for the horizontal rotation of the viewpoint
                If xyrotation=0, the view is from the South toward North.

            figsize : tuple
                contains the (x,y) dimension of the figure

        Returns
        -------
        tuple
            The value minimum and maximum of the properties
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            raise ImportError("karstnet.plot3 requires matplotlib.pyplot")
        try:
            from mpl_toolkits.mplot3d import Axes3D
        except ImportError:
            raise ImportError("karstnet.plot3 requires mpl_toolkits.mplot3d "
                              )

        if prop_type is None:
            print("There is no properties data")
            self._plot3(self.graph, figsize=figsize)

        dict_prop = nx.get_node_attributes(graph, prop_type)
        val = list(dict_prop.values())
        max_v = np.nanmax(val)
        min_v = np.nanmin(val)
        diff = max_v - min_v

        list1, list2, list3, list4, nd = ([], [], []), ([], [], []), (
            [], [], []), ([], [], []), []

        for k, v in dict_prop.items():
            if np.isnan(v):
                nd.append(k)
            if v <= min_v + diff / 2:
                if v <= min_v + diff / 4:
                    for i in range(3):
                        list1[i].append(self.pos3d[k][i])
                else:
                    for i in range(3):
                        list2[i].append(self.pos3d[k][i])
            else:
                if v >= max_v - diff / 4:
                    for i in range(3):
                        list4[i].append(self.pos3d[k][i])
                else:
                    for i in range(3):
                        list3[i].append(self.pos3d[k][i])
        nodata = (
            [self.pos3d[k][0] for k in nd], [self.pos3d[k][1] for k in nd],
            [self.pos3d[k][2] for k in nd])
        datas = (list1, list2, list3, list4, nodata)
        colors = ('chartreuse', 'yellow', 'orange', 'red', 'lightblue')
        groups = ('list1', 'list2', 'list3', 'list4', 'nodata')

        fig = plt.figure(figsize=figsize)
        # ax = fig.gca(projection='3d')
        ax = Axes3D(fig)

        for d, c, g in zip(datas, colors, groups):
            x, y, z = d
            ax.scatter(x, y, z, alpha=0.8, c=c, s=30, label=g)

        for i, j in enumerate(self.graph.edges()):
            x = np.array((self.pos3d[j[0]][0], self.pos3d[j[1]][0]))
            y = np.array((self.pos3d[j[0]][1], self.pos3d[j[1]][1]))
            z = np.array((self.pos3d[j[0]][2], self.pos3d[j[1]][2]))

            # Plot the connecting lines
            ax.plot(x, y, z, c='black', alpha=0.5)

        # # Set the view
        ax.view_init(elev=zrotation, azim=-xyrotation - 90)

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')

        # noinspection PyTypeChecker
        return min_v, max_v

    def _plot3(self, graph: nx.Graph, zrotation: float = 30,
               xyrotation: float = 0, figsize: Tuple[float] = (4, 3)):
        """
        NOT PUBLIC

        Plot a 3D view of a graph G that could be the simplified or the
        complete one.
        Requires self.pos3d => member function
        Called by the plot3 public function

        Parameters
        ----------
            graph : nx.Graph
                A networkx graph object

            zrotation : float
                angle in degrees between horizontal plane and viewpoint

            xyrotation : float
                angle in degree for the horizontal rotation of the viewpoint
                If xyrotation=0, the view is from the South toward North.

            figsize : tuple
                contains the (x,y) dimension of the figure
        """

        # 3D  plot
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            raise ImportError("karstnet.plot3 requires matplotlib.pyplot")
        try:
            from mpl_toolkits.mplot3d import Axes3D
        except ImportError:
            raise ImportError("karstnet.plot3 requires mpl_toolkits.mplot3d ")

        fig = plt.figure(figsize=figsize)
        ax = Axes3D(fig)

        for i, j in enumerate(graph.edges()):
            x = np.array((self.pos3d[j[0]][0], self.pos3d[j[1]][0]))
            y = np.array((self.pos3d[j[0]][1], self.pos3d[j[1]][1]))
            z = np.array((self.pos3d[j[0]][2], self.pos3d[j[1]][2]))

            # Plot the connecting lines
            ax.plot(x, y, z, c='black', alpha=0.5)

        # Set the view
        ax.view_init(zrotation, -xyrotation - 90)

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')

        return

    # *******************************
    # Private function for export
    # *******************************

    def _ilines_to_pline(self, list_iline: list, basename: str):
        """
        NOT PUBLIC

        Creates a Pline (Gocad ascii object) from a list of ilines
        Used to export either complete or simplified graph

        Arguments:
        ----------
            G: the name of the graph to export
            (allows having a single function for complete or simplified graph)

            list_iline : list of the ilines to write
            basename: A string containing the base name of output file.

            Manages the colocated vertices indicated by the mention "ATOM"
            in the ASCII file.

            WARNING: this version does not export (for the moment)
            the properties on the nodes.

        Returns:
        --------
            Write a output pline file called: MyKarst_exported.pl

        Example:
        --------
            >>> my_kgraph = KGraph([],{})
            >>> my_kgraph.to_pline("MyKarst" )
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
        # key is the same as in dict_nodes, the node index
        # value is the corresponding number of vrtx in the iline file,
        # due to the specific numbering of this file format,
        # it is different of the node index
        dict_added_nodes = {}

        # To count vertices: in plines,
        # vertices are virtually duplicated in the counting
        cpt_vrtx = 1
        # Each branch would be an iline
        print('list_iline : ', list_iline)
        for iline in list_iline:
            print('iline : ', iline)
            f_pline.write('ILINE\n')

            # Memorize counting state to write correctly the segments
            cpt_vrtx_deb = cpt_vrtx

            # Each node of a iline is written as a vertex or atom
            for node in iline:
                # First, verify that this node has not already been
                # added to choose between vrtx or atom
                if node not in dict_added_nodes:
                    f_pline.write('VRTX ' + str(cpt_vrtx) + ' ' +
                                  str(self.pos3d[node][0]) + ' ' +
                                  str(self.pos3d[node][1]) + ' ' +
                                  str(self.pos3d[node][2]) + '\n')
                    # Update dict_added_nodes to indicate that the node
                    # has already been declared
                    # and store the correct index in the pline domain
                    dict_added_nodes[node] = cpt_vrtx
                # if node is in dict_added_nodes, we must build an atom
                # refering to the vrtx number in the pline
                else:
                    f_pline.write('ATOM ' + str(cpt_vrtx) + ' ' +
                                  str(dict_added_nodes[node]) + '\n')
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

    def _ilines_to_dat(self, list_iline, basename):
        """
        NOT PUBLIC

        Creates a dat file from a pline file
        Used to export either complete or simplified graph

        Arguments:
        ----------
        G: the name of the graph to export
        (allows having a single function for complete or simplified graph)

        list_iline : list of the ilines to write
        basename: A string containing the base name of output file.

        WARNING: this version does not export (for the moment)
        the properties on the nodes.

        Returns:
         --------
         Write two output dat files called:
            basename_exported_links.dat
            basename_exported_nodes.dat


         Example:
         --------
            >>> my_kgraph = KGraph([],{})
            >>> my_kgraph.to_dat("MyKarst" )
         """

        output_file_nodes = basename + '_export_nodes.dat'
        output_file_links = basename + '_export_links.dat'
        print('output_file_links : ', output_file_links)
        print('output_file_nodes : ', output_file_nodes)

        f_dat_nodes = open(output_file_nodes, 'w')
        f_dat_links = open(output_file_links, 'w')

        dict_to_add = dict()

        for iline in list_iline:
            f_dat_links.write(str(iline[0] + 1) + ' ' +
                              str(iline[1] + 1) + '\n')
            for node in iline:
                dict_to_add[node + 1] = self.pos3d[node]

        for i in range(len(dict_to_add.keys())):
            f_dat_nodes.write(
                str(self.pos3d[i+1][0]) + ' ' +
                str(self.pos3d[i+1][1]) + ' ' +
                str(self.pos3d[i+1][2]) + '\n')

        f_dat_links.close()
        f_dat_nodes.close()
        print('File created')

    # *******************************
    # Private functions used by constructors
    # *******************************

    def _set_graph_lengths(self):
        """
        NOT PUBLIC

        Compute edge length at the creation of KGraph object.
        This function is called by all constructors.
        It updates graph_.
        """

        # Creation of a dictionary to store the length of each edge
        length = {}
        for e in self.graph.edges():
            dx = self.pos3d[e[0]][0] - self.pos3d[e[1]][0]
            dy = self.pos3d[e[0]][1] - self.pos3d[e[1]][1]
            dz = self.pos3d[e[0]][2] - self.pos3d[e[1]][2]
            length[e] = np.sqrt(dx ** 2 + dy ** 2 + dz ** 2)
        # Storing the length as an edge attribute
        nx.set_edge_attributes(self.graph, length, 'length')

        return

    @staticmethod
    def _calc_radius(height: float, width: float) -> float:
        """
        NOT PUBLIC

        Compute the mean radius of a section.

        Parameters
        ----------
            height : float
                The height of the network in a specific point, or equal to the
                sum of the top and the bottom of the network.

            width : float
                The width of the network in a specific point, or equal to the
                sum of right and left of the network.

        Returns:
        --------
        float
            The mean radius of the section in a specific node.
        """
        return (height + width) / 4

    @staticmethod
    def _calc_wh_ratio(height: float, width: float) -> float:
        """
         NOT PUBLIC

        Compute the Wh ratio of a section.

        Parameters
        ----------
            height : float
                The height of the network in a specific point, or equal to the
                sum of the top and the bottom of the network.

            width : float
                The width of the network in a specific point, or equal to the
                sum of right and left of the network.

        Returns:
        --------
        float
            The wh ratio of the section in a specific node.
        """
        return width / height

    def _create_dict(self) -> List[Dict[int, float]]:
        """
         NOT PUBLIC

        Create a list of dictionary of each property

        Returns:
        --------
        list
            List of dictionary of each property
        """

        dict_prop = []
        if self.properties != {}:
            component_nbr = len(self.properties[1])
            node = len(self.properties.keys())
            dict_prop = []
            radius, wh_ratio = {}, {}

            if 1 <= component_nbr <= 2:
                dict_prop = [
                    {k: v[i] for k, v in self.properties.items() if v[i] != 0}
                    for i in range(component_nbr)]
            for key in range(1, node + 1):
                if component_nbr == 3:
                    radius[key] = self._calc_radius(self.properties[key][1],
                                                    self.properties[key][2])
                    wh_ratio[key] = self._calc_wh_ratio(self.properties[key][1],
                                                        self.properties[key][2])

                if component_nbr == 4:
                    radius[key] = self._calc_radius(
                        self.properties[key][0] + self.properties[key][1],
                        self.properties[key][2] + self.properties[key][3])
                    wh_ratio[key] = self._calc_wh_ratio(
                        self.properties[key][0] + self.properties[key][1],
                        self.properties[key][2] + self.properties[key][3])

            if wh_ratio != {}:
                dict_prop.append(wh_ratio)
            if radius != {}:
                dict_prop.append(radius)
        return dict_prop

    def _set_graph_properties(self):
        """
        NOT PUBLIC

        Compute nodes properties at the creation of KGraph object.
        This function is called by all constructors.
        It updates graph_.
        """
        list_dict = self._create_dict()
        if len(list_dict) > 0:
            nx.set_node_attributes(self.graph, list_dict[0], "radius")
            # Warning : if their are 2 components or more, the user need to
            # input wh_ratio and radius in this order

            if len(list_dict) >= 2:
                nx.set_node_attributes(self.graph, list_dict[1], "wh_ratio")

            if len(list_dict) == 3:
                nx.set_node_attributes(self.graph, list_dict[0], "height")
                nx.set_node_attributes(self.graph, list_dict[1], "width")

            if len(list_dict) == 4:
                nx.set_node_attributes(self.graph, list_dict[0], "up")
                nx.set_node_attributes(self.graph, list_dict[1], "down")
                nx.set_node_attributes(self.graph, list_dict[2], "right")
                nx.set_node_attributes(self.graph, list_dict[3], "left")
        else:
            print("There is no properties data")

        return

    def _simplify_graph(self) -> Tuple[List[List[float]], nx.Graph]:
        """
        NOT PUBLIC

        Constructs a simplified graph by removing nodes of degree 2.
        Member function:
          Use self.graph (with its "length" attribute on edges) and self.pos3d
          Use self.branches produced by get_allbranches_

        Returns:
        --------
        list
            list_simple_edges : the list of simple edges (necessary for export
                to pline)
        KGraph object
            gs: the simplified output graph object

        """

        # Deals with cycles and loops to ensure that topology is not changed
        simple_edges = _split_branches(self.branches)

        # Creates a new empty graph
        gs = nx.Graph()

        # Read the dictionary of length for each edge
        length = nx.get_edge_attributes(self.graph, 'length')

        # Creates the dictionary for length  of edges
        edges_length = {}

        # list of simple_edges for export
        list_simple_edges = []
        # Fills the graph with simple_edges
        for i in simple_edges:
            # Add the edge corresponding to the simple_edges
            gs.add_edge(i[0], i[-1])
            list_simple_edges.append([i[0], i[-1]])

            # Compute the length of the current edge
            l_edge = 0
            for k in range(len(i) - 1):
                local_edge = (i[k], i[k + 1])
                if length.__contains__(local_edge):
                    l_edge += length[local_edge]
                else:
                    local_edge = (i[k + 1], i[k])
                    if length.__contains__(local_edge):
                        l_edge += length[local_edge]
                    else:

                        if self.verbose:
                            print("Warning: could not find ",
                              "1 edge when computing length")

            edges_length[(i[0], i[-1])] = l_edge

        # Stores the results
        nx.set_edge_attributes(gs, edges_length, 'length')
        return list_simple_edges, gs

    def _getallbranches(self) -> (Tuple[List[List[float]], List[float],
                                        List[float]]):
        """
        NOT PUBLIC

        Constructs the list of all branches of the karstic graph self_graph.
        Compute lengths and tortuosities
        """
        # Initialisations
        target = []
        degree_target = []

        # Create one subgraph per connected components
        # to get isolated loops as branches
        # Return a list of connected graphs
        list_sub_gr = [self.graph.subgraph(c).copy()
                       for c in nx.connected_components(self.graph)]

        for sub_gr in list_sub_gr:
            local_counter = 0
            last_node_index = 0
            # Identifies all the extremeties of the branches (nodes of degree
            # != 2)
            for i in sub_gr.nodes():
                if (sub_gr.degree(i) != 2):
                    target.append(i)
                    degree_target.append(nx.degree(sub_gr, i))
                    local_counter += 1
                last_node_index = i
            # to manage cases where a subgraph is only composed of nodes of
            # degree 2
            if (local_counter == 0):
                target.append(last_node_index)
                degree_target.append(nx.degree(sub_gr, last_node_index))

        # Identifies all the neighbors of those nodes,
        # to create all the initial paths
        list_start_branches = []
        for i in target:
            for n in self.graph.neighbors(i):
                list_start_branches.append([i, n])

        # Follow all these initial paths to get all the branches
        branches = []
        for path in list_start_branches:
            go = True
            # Check all existing branches to avoid adding a branch twice
            # if starting from other extremity
            for knownbranch in branches:
                if ((path[0] == knownbranch[-1]) &
                        (path[1] == knownbranch[-2])):
                    go = False
                    break
            if go:
                branches.append(self._getbranch(path))

        # Compute the list of branch lengths and tortuosities
        br_lengths = []
        br_tort = []

        # Read the dictionary of length for each edge
        length = nx.get_edge_attributes(self.graph, 'length')

        # To count the number of looping branches, for which tortuosity is Nan
        nb_of_Nan = 0
        for br in branches:

            # Computes the distance between extremities IF they are different
            if self.pos3d[br[0]] != self.pos3d[br[-1]]:
                dx = self.pos3d[br[0]][0] - self.pos3d[br[-1]][0]
                dy = self.pos3d[br[0]][1] - self.pos3d[br[-1]][1]
                dz = self.pos3d[br[0]][2] - self.pos3d[br[-1]][2]
                dist = np.sqrt(dx ** 2 + dy ** 2 + dz ** 2)
            else:
                dist = 0

            # Compute the length of the current branch
            br_len = 0
            # tort = 0
            # we don't want to treat the last node of the branch
            for k in range(0, len(br) - 1):
                local_edge = (br[k], br[k + 1])
                if length.__contains__(local_edge):
                    # ??????? local_edge in length
                    br_len += length[local_edge]
                else:
                    local_edge = (br[k + 1], br[k])
                    if length.__contains__(local_edge):
                        br_len += length[local_edge]
                    else:
                        print("Warning: could not find ",
                              "1 edge when computing length")

            br_lengths.append(br_len)
            # dist = 0 when positions are not defined
            # or when we have a loop
            if dist != 0:
                tort = br_len / dist
                br_tort.append(tort)
            else:
                # print("Warning: tortuosity is infinite on a looping branch.",
                # "It is set to NAN to avoid further errors.")
                # Pauline : On real systems, this message appears too many
                # times and let the user thinks something goes wrong
                br_tort.append(np.nan)
                nb_of_Nan += 1
        
        if self.verbose:        
            print(
                "Warning: This network contains ",
                nb_of_Nan,
                "looping branche.s",
                "Tortuosity is infinite on a looping branch.",
                "The looping branches are not considered for the mean tortuosity",
                "computation\n")
            
        return branches, np.array(br_lengths), np.array(br_tort)

    # ***********Functions relating to branches of graphs.
    #     A branch is defined between two nodes of degree < > 2

    def _nextstep(self, path: List[int]) -> Tuple[List[int], bool]:
        """
        NOT PUBLIC

        Work on self_graph
        Adds the next node to a path of self_graph along a branch.
        Stops when reaches a node of degree different from 2.

        Parameters
        ----------
            path : list
                A list of nodes to explain a path

        Returns:
        --------
        list
            path : A list of nodes to explain a path

        bool
            stopc
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
        # noinspection PyUnboundLocalVariable
        path.append(nextn)

        # Test for a closed loop / even if start node has degree = 2
        test_loop = path[0] == path[-1]

        if (self.graph.degree(nextn) != 2) or test_loop:
            stopc = False
        else:
            stopc = True

        return path, stopc

    def _getbranch(self, path: List[int]) -> List[int]:
        """
        NOT PUBLIC

        Work on self_graph
        Construct a branch from a starting node.

        Parameters
        ----------
            path : list
                A list of nodes to explain a path

        Returns:
        --------
        tuple
            A tuple included the list of the path from a starting node.
        """

        path, stopc = self._nextstep(path)
        while stopc:
            path, stopc = self._nextstep(path)
        return path

    # *******************************
    # Private functions used for orientations
    # *******************************

    def _set_graph_orientations(self):
        """
        NOT PUBLIC

        Compute edge length at the creation of KGraph object.
        This function is called by all constructors.
        It updates graph_.
        """

        # Creation of a dictionary to store the projected length of each edge,
        # the dip and the azimuth
        length2d = {}
        dip = {}
        azimuth = {}
        for e in self.graph.edges():
            dx = self.pos3d[e[0]][0] - self.pos3d[e[1]][0]
            dy = self.pos3d[e[0]][1] - self.pos3d[e[1]][1]
            dz = self.pos3d[e[0]][2] - self.pos3d[e[1]][2]
            length2d[e] = np.sqrt(dx ** 2 + dy ** 2)

            if length2d[e] != 0:
                dip[e] = np.arctan(abs(dz) / length2d[e])  # returns in radians
                dip[e] = np.degrees(dip[e])
                if (dz < 0):
                    dip[e] = -dip[e]
                if (dx * dy > 0):  # azimuth is comprised between 0 and 90°
                    # returns in radians
                    azimuth[e] = np.arcsin(abs(dx) / length2d[e])
                    # converts in degrees
                    azimuth[e] = np.degrees(azimuth[e])
                else:  # azimuth is comprised between 90° and 180°
                    azimuth[e] = 90 + \
                                 np.degrees(np.arccos(abs(dx) / length2d[e]))

                # to group 0 and 180 inside same bins
                azimuth[e] = np.fmod(azimuth[e], 180)
            else:  # would arrive for pure vertical segments
                azimuth[e] = np.nan
                # azimuth[e] = 0.0 #Convention
                dip[e] = 90  # degrees

        # Storing the length as an edge attribute
        nx.set_edge_attributes(self.graph, length2d, 'length2d')
        nx.set_edge_attributes(self.graph, azimuth, 'azimuth')
        nx.set_edge_attributes(self.graph, dip, 'dip')

        return


# -------------------END of KGraph class-------------------------------

# **************************************************************
#
# -------------------NON public functions used by KGraph
# ----- (can not be placed in another file)       --------------
#
# **************************************************************


def _pos_initialization(coordinates: Dict[int, List[float]]) -> (
        Tuple[Dict[int, List[float]]]):
    """
    NOT PUBLIC

    Create a dictionary of 3d coordinates from 2d or 3d input coordinates.
    Create a dictionary of x and z coordinates for each node.
    If only x, y are provided, z is set to 0

    Parameters
    ----------
        coordinates : dict
            A dictionary of 2d or 3d coordinates

    Returns:
    --------
    dict
        pos2d : a dictionary of 2d coordinates from 2d or 3d input coordinates.
    dict
        pos3d : a dictionary of 2d coordinates from 2d or 3d input coordinates.
                If only x and y are provide, z is set to 0
    dict
        nodes_xz : a dictionary of x and z coordinates for each node.
                    If only x and y are provide, z is set to 0
    dict
        nodes_yz : a dictionary of y and z coordinates for each node.
                    If only x and y are provide, z is set to 0
    """

    coord_are_3d = True
    for key in coordinates.keys():
        if ((len(coordinates[key]) == 3 and np.isnan(coordinates[key][2]))
                or len(coordinates[key]) == 2):
            coord_are_3d = False
            break

    pos3d = {}
    pos2d = {}
    nodes_xz = {}
    nodes_yz = {}

    #  if coordinates are 3d
    if coord_are_3d:
        pos3d = coordinates
        for key, coord in coordinates.items():
            pos2d[key] = [coord[0], coord[1]]
            nodes_xz[key] = [coord[0], coord[2]]
            nodes_yz[key] = [coord[1], coord[2]]

    # if only x and y are provided, set a z value = 0
    else:
        for key, coord in coordinates.items():
            pos2d[key] = [coord[0], coord[1]]
            pos3d[key] = [coord[0], coord[1], 0.]
            nodes_xz[key] = [coord[0], 0.]
            nodes_yz[key] = [coord[1], 0.]
    # noinspection PyTypeChecker
    return pos2d, pos3d, nodes_xz, nodes_yz


# ******Functions used for graph simplification


def _split2(list_: List[float]) -> Tuple[List[float], List[float]]:
    """
    NOT PUBLIC

    Splits a list in 2 sublists.

    Parameters
    ----------
        list : list

    Returns:
    --------
    tuple
        the list input split in 2 sublists included in the tuple.

    """
    list_length = len(list_)
    if (list_length == 2):
        return (list_, [])
    else:
        midpoint = int(list_length / 2)
        return (list_[0:midpoint + 1], list_[midpoint:list_length])


def _split3(list_: List[float]) -> Tuple[List[float], List[float], List[float]]:
    """
    NOT PUBLIC

    Splits a list in 3 sublists.

    Parameters
    ----------
        list : list

    Returns:
    --------
    tuple
        the list input split in 3 sublists included in the tuple.
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


def _split_branches(branches: List[List[float]]) -> List[List[float]]:
    """
    NOT PUBLIC

    Split branches in cases of loop or cycles.

    Parameters
    ----------
        branches : list

    Returns:
    --------
    list
        the branches split in case of loop

    """

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
    simple_edges = []

    for key in list_branches:
        nbb = len(list_branches[key])
        # We test first if this is a loop (same start and end point)
        if key[0] == key[1]:
            is_loop = True
        else:
            is_loop = False

        # Simple case - no cycle but potential loop
        if (nbb == 1):
            tmp = branches[list_branches[key][0]]
            if is_loop:  # In the loop case, we need to split by 3 the branch
                tmp1, tmp2, tmp3 = _split3(tmp)
                simple_edges.append(tmp1)
                if (len(tmp2) > 0):
                    simple_edges.append(tmp2)
                    if (len(tmp3) > 0):
                        simple_edges.append(tmp3)
            else:
                simple_edges.append(tmp)

        # Several branches with same extremities - cycles
        else:
            for i in range(nbb):
                tmp = branches[list_branches[key][i]]
                if is_loop:  # Case with multiple loops
                    tmp1, tmp2, tmp3 = _split3(tmp)
                    simple_edges.append(tmp1)
                    if (len(tmp2) > 0):
                        simple_edges.append(tmp2)
                        if (len(tmp3) > 0):
                            simple_edges.append(tmp3)
                else:
                    # A regular branch belonging to a cycle, need to split in 2
                    tmp1, tmp2 = _split2(tmp)
                    simple_edges.append(tmp1)
                    if len(tmp2) > 0:
                        simple_edges.append(tmp2)
    return simple_edges


def _calc_length(shortest_path: Dict[Tuple[int], List[int]],
                 length: Dict[Tuple[int, int], List[int]],
                 round_h: int = 2) -> Dict[Tuple[int], int]:
    """
    NOT PUBLIC

    Compute the length of a path between two nodes, a source node and a target
    node.

    Parameters
    ----------
        shortest_path : dict
            Dictionary of the shortest path between each nodes of the
            network with the length of the calculated paths

        length : dict
            Dictionary of the length between nodes

        round_h : int
            Number of digits after the decimal point desired for
            calculating the length of the paths

    Returns:
    --------
    dict
        Dictionary like {(source node, target node) : length of the path}
    """
    for key, path in shortest_path.items():
        h = 0
        for i_node in range(len(path) - 1):
            road = path[i_node], path[i_node + 1]
            if path[i_node] == path[i_node + 1]:
                h = 0
            if road in length.keys():
                h += length[road]
            else:
                road = path[i_node + 1], path[i_node]
                h += length[road]
        shortest_path.update({key: round(h, round_h)})
    # noinspection PyTypeChecker
    return shortest_path


def _convert(shortest_path: Dict[int, Dict[int, List[int]]]) -> (
        Dict[Tuple[int, int], List[int]]):
    """
    NOT PUBLIC

    Parameters
    ----------
        shortest_path : dict_length
            Dictionary of the shortest path between each nodes of the
            network with the length of the calculated paths

    Returns:
    --------
    dict_length
         Dictionary like {(source node, target node) : length of the path}
    """
    dict_length = {}
    for key, value in shortest_path.items():
        for target in value.keys():
            node = (key, target)
            node_inv = (target, key)
            if node_inv not in dict_length.keys():
                # We remove the duplicates tuples as (0,1) and (1,0)
                dict_length[node] = value[target]

    # dict_length = {(k, v_k): v_v for k, v in shortest_path.items() for v_k,
    # v_v in v.items()}
    # /!\ We can not remove the duplicates !

    return dict_length


def _del_node(shortest_path: Dict[int, Dict[int, List[int]]],
              prop_dict: Dict[int, float]) -> (
        Dict[int, Dict[int, List[int]]]):
    """
    NOT PUBLIC

    Parameters
    ----------
        shortest_path : dict
            Dictionary of the shortest path between each nodes of the
            network with the length of the calculated paths,
            {source node:{target node: [nodes path], ...}, ...}

        prop_dict : dict
            Dictionary of the properties of each node
            {node: properties, ...}

    Returns:
    --------
    dict
         Shortest_path without nodes with no data properties
    """
    list_key1 = list(shortest_path.copy().keys())

    for key1 in list_key1:
        if np.isnan(prop_dict[key1]):
            del shortest_path[key1]

        else:
            list_key2 = list(shortest_path[key1].copy().keys())
            for key2 in list_key2:
                if np.isnan(prop_dict[key2]):
                    del shortest_path[key1][key2]

    return shortest_path


# ******Functions used for plot

def _plot_legend(min_v: float, max_v: float):
    """
    NOT PUBLIC

    Parameters
    ----------
        min_v : float
            The minimum value of the property plot by of the function
            plot_properties2 or plot_properties3

        max_v : float
            The maximum value of the property plot by of the function
            plot_properties2 or plot_properties3
    """
    green_patch = mtb.patches.Rectangle((0, 0), 0, 0, color='chartreuse')
    yellow_patch = mtb.patches.Rectangle((0, 0), 0, 0, color='yellow')
    orange_patch = mtb.patches.Rectangle((0, 0), 0, 0, color='orange')
    red_patch = mtb.patches.Rectangle((0, 0), 0, 0, color='red')
    blue_patch = mtb.patches.Rectangle((0, 0), 0, 0, color='lightblue')
    plt.legend([red_patch, orange_patch, yellow_patch, green_patch, blue_patch],
               ['[ ' + str(round(max_v - (max_v - min_v) / 4, 2)) + ' , ' + str(
                   round(max_v, 2)) + ' ]',
                '[ ' + str(round(min_v + (max_v - min_v) / 2, 2)) + ' , ' + str(
                    round(max_v - (max_v - min_v) / 4, 2)) + ' ]',
                '[ ' + str(round(min_v + (max_v - min_v) / 4, 2)) + ' , ' + str(
                    round(min_v + (max_v - min_v) / 2, 2)) + ' ]',
                '[ ' + str(round(min_v, 2)) + ' , ' + str(
                    round(min_v + (max_v - min_v) / 4, 2)) + ' ]',
                'without data'], loc='best', frameon=False)
