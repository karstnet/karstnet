import os
import networkx as nx
import geopandas as gpd
from shapely.geometry import Point, LineString


# def networkx_to_shapefile(G,
#               crs,
#               outputdir='',
#               name=''):
#     '''
#     Transform graph data into an esri shapefile. Function written by Ana Tanaka
#     crs OxBelHa: 'EPSG:32616'
#     '''


#     #create shapefile with points
#     positions = nx.get_node_attributes(G, 'pos')
#     node_data = {'id': list(positions.keys()), 'geometry': [Point(pos) for pos in positions.values()]}
#     node_gdf = gpd.GeoDataFrame(node_data, crs=crs)
#     node_gdf.to_file(os.path.join(outputdir,f'{name}nodes.shp'))

#     #create associated shapefile with lines connecting the points.
#     edge_data = []
#     for u, v in G.edges():
#         edge_data.append({'geometry': LineString([positions[u], positions[v]]), 'source': u, 'target': v})

#     edge_gdf = gpd.GeoDataFrame(edge_data, crs=crs)   
#     edge_gdf.to_file(os.path.join(outputdir,f'{name}edges.shp'))


def networkx_to_shapefile(G,
                          crs,
                          name='cave',
                          outputdir='',
                          line=True,
                          point=True,
                          selected_nodes= None,
                          selected_edges = None,
                          pos_attr='pos',
                          ):
    """Export nodes and edges to two separate shapefiles for Gis visualization. 
    Requires the field 'pos' as an attribute of the graph. pos must be a list of three values x,y,z, or at least x,y.
    Modified from Transform graph data into an esri shapefile. Function written by Ana Tanaka

    Parameters
    ----------
    G : networkx.Graph
        networkx graph (that contains at least one node property 
        corresponding to the position of the nodes)
    crs : string
        Georeferencing code
        ex: crs OxBelHa: 'EPSG:32616'
    name : str, optional
        name of the file, by default 'cave'.
    outputdir : str, optional
        path to the folder where the file will be saved, by default ''.
    line : Bolean
        if True: create a line shapefile representing the edges
        if False: does not create a line shapefile
        By default None.
    point : Bolean
        if True: create a point shapefile representing the nodes
        if False: does not create a point shapefile
        By default None.
    selected_nodes: list
        list of nodes to export. By default None.
    selected_edges: list of list or list of tuple
        list of edges to export. By default None.
    pos_attr: string
        name of the dictionnary key containing the coordinates. By default 'pos'. 
        




    """

    #format node data for shapefile export, and save a list of edges
    if selected_nodes is None:
        #take all the node position and edges for the shapefile export
        positions = nx.get_node_attributes(G, pos_attr)
        list_edges = G.edges()
    else:
        # here there is an option to create a shapefile with just a selection of nodes
        positions = nx.get_node_attributes(nx.subgraph(G,selected_nodes), pos_attr)
        list_edges = nx.subgraph(G,selected_nodes).edges()
    
    node_data = {'id': list(positions.keys()), 'geometry': [Point(pos) for pos in positions.values()]}
    
    
    # save point shapefile
    if point == True:
        node_gdf = gpd.GeoDataFrame(node_data, crs=crs)
        node_gdf.to_file(os.path.join(outputdir,f'{name}_nodes.shp'))

    # save line shapefile
    #but first, check if there is any edges. selecting nodes can create a situation with no edges
    if list_edges and line == True:
        edge_data = []

        #choose if all a only part of the edges will be exported
        if selected_edges == None:
            for u, v in list_edges:
                # could add the flags here or any optional values
                edge_data.append({'geometry': LineString([positions[u], positions[v]]), 'source': u, 'target': v})
        else:
            for u, v in selected_edges:
                edge_data.append({'geometry': LineString([positions[u], positions[v]]), 'source': u, 'target': v})

        edge_gdf = gpd.GeoDataFrame(edge_data, crs=crs)   
        edge_gdf.to_file(os.path.join(outputdir,f'{name}_edges.shp'))


def networkx_from_shapefile(basename,
                        #   pos_attr='pos',
                          precision=1,
                        #   elevation = None
                          ):
    
    gdf = gpd.read_file(basename)
 
    #Create graph from shapefile
    G = nx.Graph()
    coord_to_label = {}

    for index, row in gdf.iterrows():
        line = row['geometry']
        if line is not None:
            coordinates = [(round(coord[0], precision), round(coord[1], precision)) for coord in line.coords] 

            for i, pos in enumerate(coordinates):
                if pos not in coord_to_label:
                    label_index = len(coord_to_label)
                    coord_to_label[pos] = label_index
                    G.add_node(label_index, pos=list(pos))

            edges = [(coord_to_label[coordinates[i]], coord_to_label[coordinates[i + 1]]) for i in range(len(coordinates) - 1)]
            G.add_edges_from(edges)

    return G
