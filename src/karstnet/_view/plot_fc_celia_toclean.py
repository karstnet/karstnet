# import numpy as np
# import networkx as nx
# import pyvista as pv
# import kntools as kt


# # =============================================================================
# # Utils for plotting graph in 3D with py vista from Julien Straubaar, 
# # github.com/ERC-Karst/karst_networks_gen/utils/graph_utils.py
# # =============================================================================
# # -----------------------------------------------------------------------------
# def get_mesh_from_graph(G, start_id=0, pos_attr='pos'):
#     """
#     Get mesh of a graph (for plotting with pyvista).

#     The graph `G` in input is assumed to have nodes localized in 3D space.

#     Parameters
#     ----------
#     G : networkx.Graph object
#         graph, with nodes with integer id, and 3D position as attribute
#     start_id : int, default: 0
#         start node id
#     pos_attr : str, default: 'pos'
#         name of the attribute attached to nodes given the position as a 
#         sequence of 2 or 3 floats; if position are in 2D, a coordinate 
#         of zero is set in the 3rd dimension

#     Returns
#     -------
#     mesh : pyvista.PolyData object
#         mesh with points (node of the graph `G`) and lines (edges of 
#         the graph `G`), that can be plotted with pyvista
#     """
#     # points: 2d-array of shape (n_nodes, 3), row i is the position in 3D of the node id i
#     points = np.array([G.nodes[ni][pos_attr] for ni in G.nodes()])
#     if points.shape[1] == 2: # 2D
#         points = np.insert(points, 2, 0.0, axis=1) # add 0.0 as 3-rd coordinate (z)

#     # lines: 2d-array of shape (n_edges, 3), row i is [2, j0, j1] where (j0, j1) is an edge btw node id j0 and node id j1
#     lines = np.insert(np.asarray(G.edges())-start_id, 0, 2, axis=1)

#     mesh = pv.PolyData(points, lines=lines.ravel())

#     return mesh
# # -----------------------------------------------------------------------------




# def plot3d_identify_disconnection(G, H, keys_manual=[] , filename='graph3d', show='True'):
#     """This function does not work yet. needs testing and possibly updating

#     Parameters
#     ----------
#     G : networkx.Graph object
#         original cave network graph before removal
#     H : networkx.Graph object
#         clean cave network graph
#     keys_manual : list, optional
#         _description_, by default []
#     filename : str, optional
#         _description_, by default 'graph3d'
#     show : str, optional
#         _description_, by default 'True'
#     """


#     pl = pv.Plotter()
    
#     lines = np.array([G.nodes('pos')[i] for i in np.hstack(G.edges())])
#     pl.add_lines(lines, color='grey', width=4,connected=False)
    
#     points = np.array([value for value in iter(dict(H.nodes('pos')).values())])
#     points_raw = np.array([value for value in iter(dict(G.nodes('pos')).values())])
    
#     pl.add_points(points, render_points_as_spheres=True, point_size=7, color='black')
#     #name the points based on the real index
    
#     if nx.number_connected_components(H)>0:
#         #colors = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)]) for i in range(nx.number_connected_components(G))]
#         colors = [ "mediumblue" , "darkturquoise", "steelblue","dodgerblue","deepskyblue","aquamarine","springgreen","seagreen","orange",'gold',"olivedrab","darkgoldenrod","tan",
#                   "mediumblue" , "darkturquoise", "steelblue","dodgerblue","deepskyblue","aquamarine","springgreen","seagreen","orange",'gold',"olivedrab","darkgoldenrod","tan",
#                   "mediumblue" , "darkturquoise", "steelblue","dodgerblue","deepskyblue","aquamarine","springgreen","seagreen","orange",'gold',"olivedrab","darkgoldenrod","tan",
#                   "mediumblue" , "darkturquoise", "steelblue","dodgerblue","deepskyblue","aquamarine","springgreen","seagreen","orange",'gold',"olivedrab","darkgoldenrod","tan",
#                   "mediumblue" , "darkturquoise", "steelblue","dodgerblue","deepskyblue","aquamarine","springgreen","seagreen","orange",'gold',"olivedrab","darkgoldenrod","tan",]
        
#         for i, subgraph_index in enumerate(nx.connected_components(H)):
#             subgraph = nx.subgraph(H, subgraph_index)
#             lines = np.array([subgraph.nodes('pos')[i] for i in np.hstack(subgraph.edges())])
#             pl.add_lines(lines, color=colors[i], width=7, connected=False)  
        
#         if keys_manual:
#             lines_reco_manual = np.array([H.nodes('pos')[i] for i in np.hstack(keys_manual)])    
#             pl.add_lines(lines_reco_manual, color='red', width=6, connected=False)
            
#         disco_keys = kt.find_disconnected_node(G, H) 
#         if disco_keys:       
#             points_disco = np.array([H.nodes('pos')[i] for i in disco_keys])
#             pl.add_points(points_disco, render_points_as_spheres=True, point_size=15, color='red')
    
#     else:
#         pl.add_lines(lines, color='grey', width=15, connected=False)
        
        
#     pl.add_point_labels(
#             points,
#             H.nodes(),
#             always_visible=True,
#             fill_shape=False,
#             margin=100,
#             shape_opacity=0.0,
#             font_size=20)
    
#     pl.add_point_labels(
#             points_raw,
#             G.nodes(),
#             always_visible=True,
#             fill_shape=False,
#             margin=100,
#             shape_opacity=0.0,
#             font_size=20)



# #---------------------------------------------------
# #plot for the big figure

# def stereo(Kg, dc, rs, weighted=True):
    
#     """
#     Density map of orientations and rose diagram of the karstic network.

#     The two stereo are ploted side by side.

#     By default, stereo and Rose diagram are weighted by lengths.

#     Examples
#     --------
#         >>> myKGraph.stereo()
#         >>> myKGraph.stereo(weighted = False)

#     """
#     import matplotlib.pyplot as plt
#     # create an np.array of azimuths and dips
#     # and lengths (projected(2d) and real (3d))
#     azim = np.array(
#         list((nx.get_edge_attributes(Kg.graph, 'azimuth')).values()))
#     azim_not_Nan = azim[~np.isnan(azim)]
#     bearing_dc = np.nan_to_num(azim)
#     plunge_dc = np.array(
#         list((nx.get_edge_attributes(Kg.graph, 'dip')).values()))
#     if (weighted):
#         l2d = np.array(
#             list((nx.get_edge_attributes(Kg.graph, 'length2d')).values()))

#         l3d = np.array(
#             list((nx.get_edge_attributes(Kg.graph, 'length')).values()))

#         l2d_not_Nan = l2d[~np.isnan(azim)]
#     else:
#         l2d_not_Nan = None
#         l3d = None
#     # Pauline: not sure it is required (?)
#     # + not sure it is normal that isnan is parameterised by azim for l2d

#     # import matplotlib as mpl

#     # Making colormap, based on Collon et al.(2017) \
#     # we saturate the colormap at 40%
#     from matplotlib import cm
#     from matplotlib.colors import ListedColormap
#     nbint = 15
#     levels = np.linspace(0, 1, nbint)
#     rainbow = cm.get_cmap('rainbow')
#     newcolors = rainbow(levels)
#     white = np.array([256 / 256, 256 / 256, 256 / 256, 1])
#     newcolors[:1, :] = white
#     newcmp = ListedColormap(newcolors)

#     # Density map - Allows to consider almost vertical conduits
#     # The data are weigthted by the real length of the segments (l3d)
#     # Use the traditional "Schmidt" method : 1% count
#     # fig = plt.figure(figsize=(16, 8))
#     # dc = fig.add_subplot(121, projection='stereonet')
#     cdc = dc.density_contourf(plunge_dc,
#                                 bearing_dc,
#                                 measurement='lines',
#                                 method='schmidt',
#                                 levels=np.arange(0, nbint * 2 + 1, 2),
#                                 extend='both',
#                                 cmap=newcmp,
#                                 weights=l3d)
#     dc.set_title('Density map of orientations \n [Schmidt\'s projection]',y=1.10)
#                     # y=1.10,
#                     # fontsize=15)
#     dc.grid()

#     # colorbar of the density map
#     cbar = plt.colorbar(cdc,
#                         fraction=0.046,
#                         pad=0.04,
#                         orientation='horizontal')
#     cbar.set_label('[%]')

#     # Rose diagram
#     # The azimuth data are weighted by the projected length (l2d)
#     bin_edges = np.arange(-5, 366, 10)
#     number_of_strikes, bin_edges = np.histogram(azim_not_Nan,
#                                                 bin_edges,
#                                                 weights=l2d_not_Nan)
#     number_of_strikes[0] += number_of_strikes[-1]
#     half = np.sum(np.split(number_of_strikes[:-1], 2), 0)
#     two_halves = np.concatenate([half, half])

#     # rs = fig.add_subplot(122, projection='polar')
#     rs.bar(np.deg2rad(np.arange(0, 360, 10)),
#             two_halves,
#             width=np.deg2rad(10),
#             bottom=0.0,
#             color='.8',
#             edgecolor='k')
#     rs.set_theta_zero_location('N')
#     rs.set_theta_direction(-1)
#     rs.set_thetagrids(np.arange(0, 360, 10), labels=np.arange(0, 360, 10))
#     rs.set_title('Rose Diagram of the \n cave survey segments',y=1.10)
#                     # fontsize=15)


#     # return fig



#     #--------------------------------
#     #unsure if we want to keep this one...


#     def specs_figure(G, 
#                  Kg,
#                  df_info,
#                  cavename = '', 
#                  number = '', 
#                  fullname = '',
#                  section_name = '',
#                  country = '',
#                  original_format = '',
#                  idealized = False
#                  ):

#     import matplotlib.pyplot as plt
#     import seaborn as sns

#     # COMPLETE FIGURE
#     ##################################################
#     ##################################################

#     title = f'{number}: $\\bf{{{cavename}}}{section_name}$'

#     if idealized == False:
#         text = [f'$\\bf{{Full name:}}$ {fullname}\n' +
#         f'$\\bf{{Country:}}$ {country}\n'+
#         f'$\\bf{{original}}$ $\\bf{{format:}}$ {original_format}\n'+
#         f'-------\n'+
#         f'$\\bf{{{df_info.loc[cavename]['number_of_nodes']}}}$ nodes\n'+
#         f'$\\bf{{{int(df_info.loc[cavename]['total_length'])}}}$m total length\n' +
#         f'$\\bf{{{int(df_info.loc[cavename]['percentage_csidm'])}}}$% geometry data (in red)\n'+
#         f'$\\bf{{{df_info.loc[cavename]['number_of_connected_components']}}}$ connected components\n'+
#         f'-------\n'+
#         f'$\\bf{{Status:}}$ {df_info.loc[cavename]['status']}\n'+
#         f'$\\bf{{Rights:}}$ {df_info.loc[cavename]['rights']}\n'+
#         f'$\\bf{{Original}}$ $\\bf{{data}}$ $\\bf{{citation:}}$ {df_info.loc[cavename]['citation']}\n'
#         f'$\\bf{{Thanks}}$ $\\bf{{to}}$ $\\bf{{include:}}$ {df_info.loc[cavename]['thanks']}\n'][0]

#     else:
#         text = [f'$\\bf{{Full name:}}$ {fullname}\n' +
#         f'-------\n'+
#         f'$\\bf{{{df_info.loc[cavename]['number_of_nodes']}}}$ nodes\n'+
#         f'$\\bf{{{int(df_info.loc[cavename]['total_length'])}}}$m total length\n' +
#         f'$\\bf{{{df_info.loc[cavename]['number_of_connected_components']}}}$ connected components\n'+
#         f'-------\n'+
#         f'This is an idealized network shape\n'][0]


#     node_sizes = []
#     labels = {}
#     for n in G.nodes:
#         if G.nodes('csdim')[n] is not None:
#             node_sizes.append( nx.get_node_attributes(G,'csdim')[n][0]*1 )
#             # print(nx.get_node_attributes(G,'csdim')[n][0])
#             # print(type(nx.get_node_attributes(G,'csdim')[n][0]))
#         else:
#             node_sizes.append( 0.0 )


#     # fig2 = Kg.plot3(zrotation=60, xyrotation=10)
#     # fig2.savefig(f'{clean_output_folder}/karstnet_analysis/{cavename}_kn_plot3d.pdf')

#     fig = plt.figure(figsize=(30,10))

#     tx = fig.add_subplot(261)  
#     # tx.text(0, 0, text, fontsize = 12, wrap=True)
#     txt = tx.text(0,0.9, title, wrap=True, fontsize=20)
#     txt = tx.text(0,0.85, text, ha='left', va='top', wrap=True, fontsize=10, transform=tx.transAxes)
#     txt._get_wrap_line_width = lambda : 1200.
#     tx.set_axis_off()

#     #Simplified graph
#     ax2 = fig.add_subplot(262)
#     nx.draw(Kg.graph_simpl,
#                 pos=get_pos2d(G),
#                 with_labels=False,
#                 node_size=0, ax=ax2)
#     # ax2.set_aspect('equal', adjustable='box')
#     ax2.set_title('Simplified graph')#,y=0.8, x=0.8)
#     ax2.set_aspect('equal')


#     #plot Densitiy and direction diagram
#     dc = fig.add_subplot(267, projection='stereonet')
#     rs = fig.add_subplot(268, projection='polar')
#     stereo(Kg, dc, rs, weighted=True)

#     #Full graph
#     ax1 = fig.add_subplot(132)
#     nx.draw(G,get_pos2d(G),node_size=node_sizes ,with_labels=False, node_color='red', ax=ax1)
#     ax1.set_aspect('equal')
#     # ax1.set_aspect('equal', adjustable='box')



#     # PLOT VIOLIN PLOT
#     # df_kn_analysis
#     if idealized == False:
#         list_stat_param =  [#['Total\nlength (m)','total_length', 'total_length'],
#                             ['Average\nwidth (m)', 'average_cs_width', 'average_cs_width'],
#                             ['Average\nheight (m)' , 'average_cs_height', 'average_cs_height'],
#                             ['Cycles' , 'number_of_cycles_(graph_simplified)', 'number_of_cycles_(graph_simplified)'],
#                             ['CV length' , 'cv_length', 'cv_length'],
#                             ['Length\nentropy' , 'length_entropy', 'length_entropy'],
#                             ['Tortuosity'  ,'tortuosity', 'tortuosity'],
#                             ['Orientation\nentropy' , 'orientation_entropy', 'orientation_entropy'],
#                             ['Average\nshortest\npath length ()' , 'aspl', 'aspl'],
#                             ['Central\nPoint\nDominance' , 'cpd', 'cpd'],
#                             ['Mean\nlength (m)' , 'mean_length', 'mean_length'],
#                             ['Mean\ndegree' , 'mean_degree', 'mean_degree'],
#                             ['Branch lengths\nvariation\ncoefficient(%)' , 'cv_degree', 'cv_degree'],
#                             ['Correlation \nVertex Degree', 'correlation_vertex_degree', 'correlation_vertex_degree']]
#     else:    
#         list_stat_param = [ ['Cycles' , 'number_of_cycles_(graph_simplified)', 'number_of_cycles_(graph_simplified)'],
#                             ['CV length' , 'cv_length', 'cv_length'],
#                             ['Length\nentropy' , 'length_entropy', 'length_entropy'],
#                             ['Tortuosity'  ,'tortuosity', 'tortuosity'],
#                             ['Orientation\nentropy' , 'orientation_entropy', 'orientation_entropy'],
#                             ['Average\nshortest\npath length ()' , 'aspl', 'aspl'],
#                             ['Central\nPoint\nDominance' , 'cpd', 'cpd'],
#                             ['Mean\nlength (m)' , 'mean_length', 'mean_length'],
#                             ['Mean\ndegree' , 'mean_degree', 'mean_degree'],
#                             ['Branch lengths\nvariation\ncoefficient(%)' , 'cv_degree', 'cv_degree'],
#                             ['Correlation \nVertex Degree', 'correlation_vertex_degree', 'correlation_vertex_degree']]

#     # list_stat_param = ['tortuosity']

#     # fig,ax = plt.subplots(3,len(list_stat_param),figsize=(10,20))

#     for i,stat_param in enumerate(list_stat_param):
#         # y, x, _ = hist(df_kn_analysis[stat_param[1]])

#         ax = fig.add_subplot(len(list_stat_param)*2,3,(i+1)*6)
#         sns.violinplot(data=df_info, x=stat_param[1],ax=ax,cut=0)
#         ax.plot(df_info.loc[cavename][stat_param[2]],0,color='white', marker='*',markersize=20, markeredgecolor='red')

#         sns.despine(ax=ax,bottom=False,left=True, trim=True)
#         ax.set_yticks([])
#         ax.set_ylabel(stat_param[0],rotation='horizontal', ha='right', va="center" ,fontsize=12)
#         ax.set_xlabel('')

#         # plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.5, hspace=1)
    
#     # plt.tight_layout()
#     # fig.savefig(f"{ux.make_filepath(path_to_github,'visualization')}/{cavename}_plot2d.png", format='png')
#     return fig
#     # fig.savefig(f'C:/Users/celia/github/erc-karst-repositories/networks_datasets/png_figures/plot2d_{number}_{cavename}{section_name}.png', dpi=300, format='png',bbox_inches='tight')