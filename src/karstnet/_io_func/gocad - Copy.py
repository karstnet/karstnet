import networkx as nx
import karstnet as kn
import numpy as np

def networkx_to_gocad(G, 
                    data_type = 'lines', #or 'points'
                    properties = [], #
                    nodata_value = '-999999999',
                    name = 'graph_gocad_export',
                    node_id = [], #list of nodes id to export
                    pos_attr = 'pos'
                    ):
    """Export networkx graph data into a format readable by Gocad (.pl).
    This export requires x,y,z coordinate information attached as an attribute on each node on graph as a list of [x,y,z]
    It is possible to only give a few points or the entire dataset. It is also possible to add several properties.

    idea for improvement: 
    1. instead of having to attached x,y,z to the graph, give a dictionnary into the function??
    2. 


    Parameters
    ----------
    G : networkx graph
        Graph with coordinates position
    data_type : str, optional
        Choose data type between 'lines' and 'points', by default 'lines'
    properties : list of string(s)
        List containing the name of all the graph attributes to add to the Gocad output, by default []
        For now, the properties have to be a single value per node. 
    nodata_value : string
        by default '-999999999'
    name : str, optional
        path and name of the file to save. The extension is .pl , by default 'graph_gocad_export'
    node_id : list of node id, optional
        List of node id to use for the export, by default []
        If [] then the entire graph is selected
    pos_attr : str, optional
        name of the attribute containing the list of coordinates, by default 'pos'

    returns:
    --------
    saves a .pl file readable in Gocad
    """


    #the properties have to be a single value per node
    nodata_string = 'NO_DATA_VALUES ' + " ".join(str(item) for item in [nodata_value] * (len(properties)+3))
    
    if data_type=='points':
        header = 'GOCAD VSet 1.0'
        header_dataset = 'SUBVSET'
        
        #nodes
        if node_id == []:
            nodes = G.nodes
        else:
            nodes = node_id

        dataset = []

        for node in nodes:               
            #write coordinate string ex: '-6.53 -10.24 -28.11'
            # print(G.nodes('pos')[node])
            string_coordinates = " ".join(str(item) for item in G.nodes(pos_attr)[node])
            #write attribute string ex: 'att1 att2 atti'
            string_attribute = ''

            # string_id = " ".join(str(item) for item in nodes)
            properties_list = []
            for attribute in properties:
                if G.nodes(attribute)[node] is not None:
                    if type(G.nodes(attribute)[node])==list:
                        if len(G.nodes(attribute)[node])==2:
                            # for example csdim
                            string_attribute += str(G.nodes(attribute)[node][0]) + ' '
                            string_attribute += str(G.nodes(attribute)[node][1]) + ' '
                            properties_list.append(f'{attribute}_1')
                            properties_list.append(f'{attribute}_2')

                        if len(G.nodes(attribute)[node])!=2:
                            print(f'The chosen properties --{attribute}-- cannot be written to Gocad file')
                    else:
                        string_attribute += str(G.nodes(attribute)[node]) + ' '
                        properties_list.append(attribute)
                else:
                    # in the case there is no data
                    string_attribute += nodata_value
            
            #create a list of lines containing the notes attributes  
            dataset.append('PVRTX '+ str(node) +  ' ' + string_coordinates + ' ' + str(node) + ' ' + string_attribute)

        #write the lines
        lines = [   header, 
                'HEADER{',
                'name:' + name ,
                '}',
                'GOCAD_ORIGINAL_COORDINATE_SYSTEM',
                'ZPOSITIVE Elevation',
                'END_ORIGINAL_COORDINATE_SYSTEM',
                'PROPERTIES ID ' + " ".join(str(item) for item in properties_list),  
                nodata_string,
                header_dataset] + dataset + ['END']
    
    if data_type=='lines':
        header = 'GOCAD PLine 1'
        header_dataset = ''
        
        H = kn.utils.graph_to_branches(G)
        
        #write lines containing pline information
        #Gocad only read the pline right by single branch, with points in the right order
        dataset = []
        if nx.is_connected(H) == False:
            #iterate through the connectec components to find nodes where disconnection occured
            #search for the nodes that used to be degree >1 and are now degree 1.
            for i, subgraph_index in enumerate(nx.connected_components(H)):
                subgraph = nx.subgraph(H, subgraph_index)
                # print(subgraph)
                #segments
                seg = []
                nodes_start_end = [node for node,degree in dict(subgraph.degree()).items() if degree ==1]
                # print(nodes_start_end)

                #make sure its not a loop
                #in the case of a segment in the form of a loop, 
                # here will be no node of degree 1, so we replace by the first node in the list
                if nodes_start_end:
                    for path in nx.all_simple_edge_paths(subgraph, nodes_start_end[0],nodes_start_end[1]):
                        for edge in path:   
                            # print(edge) 
                            seg.append(" ".join(str(item) for item in edge))              
                            #nodes
                else:
                    #create a new node at the intersection
                    #take a random node in the loop:
                    intersection = list(subgraph.nodes())[0]
                    max_value = np.array(G.nodes()).max()
                    subgraph = subgraph.copy()
                    subgraph.add_node(max_value,   pos=subgraph.nodes(pos_attr)[intersection], 
                                            connected_component_number=subgraph.nodes('connected_component_number')[intersection],
                                            intersection=intersection)  
                    #create an edge that connects the new node
                    node = list(subgraph.neighbors(intersection))[0]
                    subgraph.add_edge(max_value,node) 
                    #remove the old ege
                    subgraph.remove_edge(intersection,node)  

                    seg = []
                    nodes_start_end = [node for node,degree in dict(subgraph.degree()).items() if degree ==1]

                    for path in nx.all_simple_edge_paths(subgraph, nodes_start_end[0],nodes_start_end[1]):
                        for edge in path:   
                            # print(edge) 
                            seg.append(" ".join(str(item) for item in edge))              
                            #nodes

                
                pvrtx = []
                for node in nx.shortest_path(subgraph, source=nodes_start_end[0], target=nodes_start_end[1]):               
                    #write coordinate string ex: '-6.53 -10.24 -28.11'
                    string_coordinates = " ".join(str(item) for item in subgraph.nodes(pos_attr)[node])
                    #write attribute string ex: 'att1 att2 atti'
                    string_attribute = ''

                    properties_list = []              
                    for attribute in properties:
                        if subgraph.nodes(attribute)[node] is not None:

                            if type(subgraph.nodes(attribute)[node])==list:
                                if len(subgraph.nodes(attribute)[node])==2:
                                    # for example csdim
                                    string_attribute += str(G.nodes(attribute)[node][0]) + ' '
                                    string_attribute += str(G.nodes(attribute)[node][1]) + ' '
                                    properties_list.append(f'{attribute}_1')
                                    properties_list.append(f'{attribute}_2')
                                if len(subgraph.nodes(attribute)[node])!=2:
                                    print(f'The chosen properties --{attribute}-- cannot be written to Gocad file')
                            else:
                                string_attribute += str(subgraph.nodes(attribute)[node]) + ' '
                                properties_list.append(attribute)
                        else:
                            # in the case there is no data
                            print('no properties for node: ', node )
                            string_attribute += nodata_value
                    
                    #create a list of lines  
                    pvrtx.append('PVRTX '+ str(node) +  ' ' + string_coordinates + ' ' + string_attribute)
                # create each branch text    
                dataset += ['ILINE'] + pvrtx + seg
           
     
        lines = [   header, 
                    'HEADER{',
                    'name:' + name ,
                    '}',
                    'GOCAD_ORIGINAL_COORDINATE_SYSTEM',
                    'ZPOSITIVE Elevation',
                    'END_ORIGINAL_COORDINATE_SYSTEM',
                    'PROPERTIES ' + " ".join(str(item) for item in properties_list),  
                    nodata_string,
                    header_dataset] + dataset + ['END']
               
    with open(name + '.pl', 'w') as f:
        f.write('\n'.join(lines))