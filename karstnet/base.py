#    Copyright (C) 2018 by
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
   Copyright (C) 2018 Karstnet Developers
   Philippe Renard <philippe.renard@unine.ch>
   Pauline Collon <pauline.collon@univ-lorraine.fr>
"""

# ----External libraries importations
import numpy as np
import networkx as nx
import scipy.stats as st
import matplotlib.pyplot as plt
import sqlite3
# noinspection PyUnresolvedReferences
import mplstereonet


# *************************************************************
# -------------------Test function--------------------------
# *************************************************************
def test_kn():
    print("test ok")
    print("relance ok")


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
                Each station is a node, each line-of sigth is an edge.
                It is a Networkx graph object, with length of edges
                as attributes.
        - graph_simpl: the simplified graph of the karstic network,
                i.e. all nodes of degree 2 are removed except in loops
                (2 types of loops)(cf Collon et al, 2017, Geomorphology)
                graph_simple is a Networkx graph object, with length of
                edges as attributes.
        - pos2d : a dictionnary of the nodes with their 2d coordinates
                as a values [x, y]
        - pos3d : a dictionnary of the nodes with their 3d coordinates
                as a values [x, y, z]
        - properties : a dictionnary of nodes with their properties in
                case of importing karst data with additional information
        - branches : the list of branches
        - br_lengths : the list of branch lengths
        - br_tort : the list of branche tortuosities
        - list_simpl_edges : the list of simple edges, necessary to export
                graph to plines
        - graph_simpl : the simplified version (without nodes of degree 2)
                of a graph

    """

    def __init__(self, edges, coordinates, properties=None):
        """
        Creates a Kgraph from nodes and edges.

        Parameters
        ----------
        edges : list
            a list of edges

        coord : dictionnary
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

        self.pos2d, self.pos3d = _pos_initialization(coordinates)
        print(
            "\n This network contains ",
            nx.number_connected_components(
                self.graph),
            " connected components")

        # Compute graph length from the pos3d_ to initialize properly the graph
        self._set_graph_lengths()
        self._set_graph_orientations()
        self.properties = properties

        # Compute branches of the graph
        # self.branches is necessary to export graph to plines
        self.branches, self.br_lengths, self.br_tort = self._getallbranches()

        # Construct the simplified graph
        # self.list_simpl_edges is necessary to export graph to plines
        self.list_simpl_edges, self.graph_simpl = self._simplify_graph()

    # **********************************
    #    Plots
    # **********************************

    def plot2(self, graph_type=0, figsize=(6, 3)):
        """
        Plot a 2D view of the karstic network

        Parameters
        ----------
        graph_type : int
            if 0 displays the complete graph,
            if 1 displays the simplified graph

        figsize : tuple
            contains the (x,y) dimension of the figure

        Examples
        --------
           >>> myKGraph = KGraph([],{})
           >>> myKGraph.plot2()
           >>> myKGraph.plot2(1, zrotation=20, xyrotation=-30)
        """

        if (graph_type == 0):
            self._plot2(self.graph, figsize)
            plt.title('original')
            plt.show()
        else:
            self._plot2(self.graph_simpl, figsize)
            plt.title('simplified')
            plt.show()

        return

    def plot3(self, graph_type=0, zrotation=30, xyrotation=0, figsize=(6, 3)):
        """
        Plot a 3D view of the karstic network.

        Parameters
        ----------
        graph_type : int
            if 0 displays the complete graph,
            if 1 displays the simplified graph

        zrotation : float
            angle in degrees between horizontal plane and viewpoint

        xyrotation : float
            angle in degree for the horizontal rotation of the viewpoint.
            If xyrotation=0, the view is from the South toward North.

        figsize : tuple
             contains the (x,y) dimension of the figure

        Examples
        --------
           >>> myKGraph.plot3()
           >>> myKGraph.plot3(1, zrotation=20, xyrotation=-30)

        """
        # 3D  plot

        if (graph_type == 0):
            self._plot3(self.graph, zrotation, xyrotation, figsize)
            plt.title('original')
            plt.show()
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
           >>> myKGraph.plot()

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
        Simple 2D map of the original and simplified karstic network.

        The two maps are ploted side by side. This function allows
        to check rapidly the data after an import for example.

        Examples
        ---------
           >>> myKGraph.plot()

        """
        plt.figure(figsize=(12, 5))
        plt.subplot(121)
        nx.draw_networkx(self.graph,
                         pos=self.pos2d,
                         with_labels=False,
                         node_size=0.1)
        plt.xlabel('x')
        plt.ylabel('z')
        plt.subplot(122)
        nx.draw_networkx(self.graph_simpl,
                         pos=self.pos2d,
                         with_labels=False,
                         node_size=0.1)
        plt.xlabel('x')
        plt.ylabel('z')
        plt.show()

    # Member function written by Philippe Vernant 2019/11/25
    # Modified by Pauline Collon (aug. 2020) to weight density map by lenghts

    def stereo(self, weighted=True):
        """
        Density map of orientations and rose diagram of the karstic network.

        The two stereo are ploted side by side.

        By default, stereo and Rose diagram are weighted by lengths.

        Examples
        ---------
           >>> myKGraph.stereo()
           >>> myKGraph.stereo(weighted = False)

        """

        # create an np.array of azimuths and dips
        # and lengths (projected(2d) and real (3d))
        azim = np.array(
            list((nx.get_edge_attributes(self.graph, 'azimuth')).values()))
        azim_not_Nan = azim[~np.isnan(azim)]
        bearing_dc = np.nan_to_num(azim)
        plunge_dc = np.array(
            list((nx.get_edge_attributes(self.graph, 'dip')).values()))
        if (weighted):
            l2d = np.array(
                list((nx.get_edge_attributes(self.graph,
                                             'length2d')).values()))

            l3d = np.array(
                list((nx.get_edge_attributes(self.graph, 'length')).values()))

            l2d_not_Nan = l2d[~np.isnan(azim)]
        else:
            l2d_not_Nan = None
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
        # The data are weigthted by the real length of the segments (l3d)
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
        number_of_strikes, bin_edges = np.histogram(azim_not_Nan,
                                                    bin_edges,
                                                    weights=l2d_not_Nan)
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

    # end modif PV 2019/11/25

    # *************************************************************
    # ----------------------- Export ------------------------------
    # *************************************************************

    def to_pline(self, basename):
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
            the name contains no extension, it will be added by the function

        Examples
        --------
        The following command saves the file "MyKarst_exported.pl" :

        >>> myKGraph.to_pline("MyKarst")
        """
        # For a complete graph the list of Ilines corresponds to self.branches

        self._ilines_to_pline(self.branches, basename)

        return

    def simpleGraph_to_pline(self, basename):
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

    # *************************************************************
    # -------------------Computation for Analysis------------------
    # *************************************************************
    def basic_analysis(self):
        """
        Print the basics of a karstic network graph analysis

        Examples
        --------
        >>> t = myKGraph.basic_analysis()

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
        # consistent with the exemples of braided rivers given by Howard.
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
           >>> t = myKGraph.mean_tortuosity()
        """
        nb_of_Nan = np.isnan(self.br_tort).sum()
        if nb_of_Nan != 0:
            print(
                "\n WARNING: This network contains ",
                nb_of_Nan,
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
           >>> l = myKGraph.mean_length()
        """
        return (np.mean(self.br_lengths))

    def coef_variation_length(self):
        """
        Compute the coefficient of variation of length of the branches of a
        karstic networkx

        Returns
        -------
            cvl :  Coefficient of variation of the length of the branches

        Example:
        --------
           >>> cvl = myKGraph.coef_variation_length()
        """

        # std is used with ddof =1 to use the estimate of std (divides by N-1)
        return (np.std(self.br_lengths, ddof=1) / np.mean(self.br_lengths))

    def length_entropy(self, mode="default"):
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
            entropy :  Entropy of the lengths of the branches

        Example:
        --------
           >>> l_entrop = myKGraph.length_entropy()
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
                                     range=(0, 100))
            freq = counts / sum(counts)  # Computes the frequencies
            entropy = st.entropy(freq, base=nbins)
        else:
            entropy = 0  # v contains a single value - no uncertainty

        return entropy

    def orientation_entropy(self, mode="default"):
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
            entropy:  Entropy of the segment orientation

        Example:
        --------
           >>> or_entropy = myKGraph.orientation_entropy()
        """

        # create an np.array of azimuths and projected lengths
        azim = np.array(
            list((nx.get_edge_attributes(self.graph, 'azimuth')).values()))
        l2d = np.array(
            list((nx.get_edge_attributes(self.graph, 'length2d')).values()))

        # Removing NAN Azimuth values that correspond to length2d=0
        azim_not_Nan = azim[~np.isnan(azim)]
        l2d_not_zero = l2d[np.nonzero(l2d)]

        if (len(azim_not_Nan) > 1):

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

    def mean_degree_and_CV(self):
        """
        Computes the average and the coefficient of variation of the degree.

        The computation is done on the simplified graph.

        Returns
        -------
        tuple
            meandeg, cvdeg :  the mean and coefficient of variation

        Examples
        --------
           >>> meandeg, cvde = myKGraph.coef_variation_degree()
        """
        # Vector of degrees
        d = np.array(self.graph_simpl.degree())[:, 1]

        # Mean degree
        meandeg = np.mean(d)

        # Coefficient of variation of the degrees : std is used with ddof =1 to
        # use the estimate of std (divides by N-1)
        cvde = np.std(d, ddof=1) / np.mean(d)

        return meandeg, cvde

    def correlation_vertex_degree(self, cvde=False):
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
           >>> cvd = correlation_vertex_degree()
        """
        if not cvde:
            _, cvde = self.mean_degree_and_CV()

        # To avoid division by 0 when computing correlation coef
        if cvde != 0:
            cvd = nx.degree_pearson_correlation_coefficient(self.graph_simpl)
        else:
            cvd = 1

        return cvd

    def central_point_dominance(self):
        """
        Computes central point dominance.

        The computation is done on the simplified graph.

        Returns
        -------
        float
            Central point dominance

        Examples
        --------
           >>> cpd = myKGraph.central_point_dominance()
        """
        bet_cen = nx.betweenness_centrality(self.graph_simpl)
        bet_cen = list(bet_cen.values())
        cpd = sum(max(bet_cen) - np.array(bet_cen)) / (len(bet_cen) - 1)

        return cpd

    def average_SPL(self, dist_weight=False):
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

    def characterize_graph(self, verbose=False):
        """
        Computes the set of metrics used to characterize a graph.

        Parameters
        ----------
        verbose : boolean
            If True, the function displays information about the
            progress of the computation, and the results.

        Returns
        -------
        dictionnary
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
        >>> results = myKGraph.characterize_graph()
        """

        results = {}

        print('Computing:')
        print(' - mean length', end='', flush=True)
        results["mean length"] = self.mean_length()

        print(',cv length', end='', flush=True)
        results["cv length"] = self.coef_variation_length()

        print(',length entropy', end='', flush=True)
        results["length entropy"] = self.length_entropy()

        print(',mean tortuosity', end='', flush=True)
        results["tortuosity"] = self.mean_tortuosity()

        print('', end='\n', flush=True)

        print(' - orientation entropy', end='', flush=True)
        results["orientation entropy"] = self.orientation_entropy()

        print(',aspl', end='', flush=True)
        results["aspl"] = self.average_SPL()

        print(',cpd', end='', flush=True)
        results["cpd"] = self.central_point_dominance()

        print(',md,cv degree', end='', flush=True)
        md, cvde = self.mean_degree_and_CV()
        results["mean degree"] = md
        results["cv degree"] = cvde

        print(',cvd', end='', flush=True)
        cvd = self.correlation_vertex_degree(cvde=cvde)
        results["correlation vertex degree"] = cvd

        print('', end='\n', flush=True)

        if verbose:
            print("--------------------------------------")
            for key in results.keys():
                print(" %25s = %5.3f" % (key, results[key]))
            print("--------------------------------------")

        return (results)

    # *************************************************************************
    # Non Public member functions of KGraph class
    # *************************************************************************

    # *******************************
    # Private functions for plots
    # *******************************

    def _plot2(self, G, figsize=(6, 3)):
        """
        NOT PUBLIC
        Plot a 2D view of a graph G that could be the simplified of the
        complete one.  Requires self.pos3d => member function.
        Called by the plot3 public function

        """

        # 2D  plot
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            raise ImportError("karstnet.plot3 requires matpllotlib.pyplot")

        plt.figure(figsize=figsize)

        nx.draw_networkx(G,
                         with_labels=True,
                         pos=self.pos2d,
                         node_color='lightblue')

        return

    def _plot3(self, G, zrotation=30, xyrotation=0, figsize=(6, 3)):
        """
        NOT PUBLIC
        Plot a 3D view of a graph G that could be the simplified or the
        complete one.
        Requires self.pos3d => member function
        Called by the plot3 public function

        """

        # 3D  plot
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            raise ImportError("karstnet.plot3 requires matpllotlib.pyplot")
        try:
            from mpl_toolkits.mplot3d import Axes3D
        except ImportError:
            raise ImportError("karstnet.plot3 requires mpl_toolkits.mplot3d ")

        fig = plt.figure(figsize=figsize)
        ax = Axes3D(fig)

        for i, j in enumerate(G.edges()):
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

    def _ilines_to_pline(self, list_Iline, basename):
        """
        Creates a Pline (Gocad ascii object) from a list of ilines
        Used to export either complete or simplified graph

        Arguments:
        ----------
            G: the name of the graph to export
            (allows having a single function for complete or simplified graph)

            list_Iline : list of the ilines to write
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

        # Create an empty dictionnary of nodes already written in one iline
        # key is the same as in dico_nodes, the node index
        # value is the corresponding number of vtrx in the iline file,
        # due to the specific numbering of this file format,
        # it is different of the node index
        dico_added_nodes = {}

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
                                  str(self.pos3d[node][0]) + ' ' +
                                  str(self.pos3d[node][1]) + ' ' +
                                  str(self.pos3d[node][2]) + '\n')
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

        print('File created')

        # Close the file
        f_pline.close()
        return

    # *******************************
    # Private functions used by constructors
    # *******************************

    def _set_graph_lengths(self):
        """NON PUBLIC.
        Compute edge length at the creation of KGraph object.
        This function is called by all constructors.
        It updates graph_.
        """

        # Creation of a dictionnary to store the length of each edge
        length = {}
        for e in self.graph.edges():
            dx = self.pos3d[e[0]][0] - self.pos3d[e[1]][0]
            dy = self.pos3d[e[0]][1] - self.pos3d[e[1]][1]
            dz = self.pos3d[e[0]][2] - self.pos3d[e[1]][2]
            length[e] = np.sqrt(dx ** 2 + dy ** 2 + dz ** 2)
        # Storing the length as an edge attribute
        nx.set_edge_attributes(self.graph, length, 'length')

        return

    def _simplify_graph(self):
        """
        Constructs a simplified graph by removing nodes of degree 2.
        Member function:
          Use self.graph (with its "length" attribute on edges) and self.pos3d
          Use self.branches produced by get_allbranches_

        Returns:
        --------
           - list_simpl_edges : the list of simple edges (necessary for export
                to pline)
           - Gs: the simplified output graph object

        """

        # Deals with cycles and loops to ensure that topology is not changed
        simpl_edges = _split_branches(self.branches)

        # Creates a new empty graph
        Gs = nx.Graph()

        # Read the dictionnary of length for each edge
        length = nx.get_edge_attributes(self.graph, 'length')

        # Creates the dictionnary for length  of edges
        edges_length = {}

        # list of simpl_edges for export
        list_simpl_edges = []
        # Fills the graph with simpl_edges
        for i in simpl_edges:
            # Add the edge corresponding to the simpl_edges
            Gs.add_edge(i[0], i[-1])
            list_simpl_edges.append([i[0], i[-1]])

            # Compute the length of the current edge
            l_edge = 0
            for k in range(0, len(i) - 1):
                local_edge = (i[k], i[k + 1])
                if length.__contains__(local_edge):
                    l_edge += length[local_edge]
                else:
                    local_edge = (i[k + 1], i[k])
                    if length.__contains__(local_edge):
                        l_edge += length[local_edge]
                    else:
                        print("Warning: could not find ",
                              "1 edge when computing length")

            edges_length[(i[0], i[-1])] = l_edge

        # Stores the results
        nx.set_edge_attributes(Gs, edges_length, 'length')

        return list_simpl_edges, Gs

    def _getallbranches(self):
        """
        Constructs the list of all branches of the karstic graph self_graph.
        Compute lengths and tortuosities
        """

        # Initialisations
        target = []
        degreeTarget = []

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
                    degreeTarget.append(nx.degree(sub_gr, i))
                    local_counter += 1
                last_node_index = i
            # to manage cases where a subgraph is only composed of nodes of
            # degree 2
            if (local_counter == 0):
                target.append(last_node_index)
                degreeTarget.append(nx.degree(sub_gr, last_node_index))

        # Identifies all the neighbors of those nodes,
        # to create all the initial paths
        listStartBranches = []
        for i in target:
            for n in self.graph.neighbors(i):
                listStartBranches.append([i, n])

        # Follow all these initial paths to get all the branches
        branches = []
        for path in listStartBranches:
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

        # Read the dictionnary of length for each edge
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
            tort = 0
            # we don't want to treat the last node of the branch
            for k in range(0, len(br) - 1):
                local_edge = (br[k], br[k + 1])
                if length.__contains__(local_edge):
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

    def _nextstep(self, path):
        """
        Work on self_graph
        Adds the next node to a path of self_graph along a branch.
        Stops when reaches a node of degree different from 2.
        """

        current = path[-1]
        # Checks first if the end of the path is already on an end
        if (self.graph.degree(current) != 2):
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

        if ((self.graph.degree(nextn) != 2) or testloop):
            stopc = False
        else:
            stopc = True

        return path, stopc

    def _getbranch(self, path):
        """
        Work on self_graph
        Construct a branch from a starting node.
        """

        path, stopc = self._nextstep(path)
        while stopc:
            path, stopc = self._nextstep(path)
        return (path)

    # *******************************
    # Private functions used for orientations
    # *******************************

    def _set_graph_orientations(self):
        """NON PUBLIC.
        Compute edge length at the creation of KGraph object.
        This function is called by all constructors.
        It updates graph_.
        """

        # Creation of a dictionnary to store the projected length of each edge,
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


def _pos_initialization(coordinates):
    '''NON PUBLIC.
    Create a dictionnary of 3d coordinates from 2d or 3d input coordinates.
    If only x, y are provided, z is set to 0
    '''

    coord_are_3d = True
    for key in coordinates.keys():
        if len(coordinates[key]) == 3:
            break  # we only check one value and let coord_are_3d to True
        else:
            coord_are_3d = False
            break
    pos3d = {}
    pos2d = {}
    #  if coordinates are 3d
    if coord_are_3d:
        pos3d = coordinates
        for key, coord in coordinates.items():
            pos2d[key] = [coord[0], coord[1]]

    # if only x and y are provided, set a z value = 0
    else:
        pos2d = coordinates
        for key, coord in coordinates.items():
            pos3d[key] = [coord[0], coord[1], 0]

    return pos2d, pos3d


# ******Functions used for graph simplification


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


def _split_branches(branches):
    """
    Split branches in cases of loop or cycles.
    """

    # Creates a dictionnary to accelerate the search
    # of branches having the same extremities
    list_branches = dict()
    for i, b in enumerate(branches):
        key = (b[0], b[-1])
        if list_branches.__contains__(key):
            list_branches[key].append(i)
        else:
            list_branches[key] = [i]

    # Loop over the branches with same extremitis and split them when required
    simpl_edges = []

    for key in list_branches:
        nbb = len(list_branches[key])
        # We test first if this is a loop (same start and end point)
        if key[0] == key[1]:
            isloop = True
        else:
            isloop = False

        # Simple case - no cycle but potential loop
        if (nbb == 1):
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

        # Several branches with same extremities - cycles
        else:
            for i in range(nbb):
                tmp = branches[list_branches[key][i]]
                if isloop:  # Case with multiple loops
                    tmp1, tmp2, tmp3 = _split3(tmp)
                    simpl_edges.append(tmp1)
                    if (len(tmp2) > 0):
                        simpl_edges.append(tmp2)
                        if (len(tmp3) > 0):
                            simpl_edges.append(tmp3)
                else:
                    # A regular branch belonging to a cycle, need to split in 2
                    tmp1, tmp2 = _split2(tmp)
                    simpl_edges.append(tmp1)
                    if (len(tmp2) > 0):
                        simpl_edges.append(tmp2)

    return simpl_edges
