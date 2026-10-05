import pandas as pd
import os
import networkx as nx
import karstnet as kn


def networkx_from_kndata_csv(inputpath,
                         node_attributes=['pos','csdim']):
    """Loads the cave network from the Github repository <https://github.com/ERC-Karst/KNdata-public> and returns a networkx graph. The cave graph is stored in a folder containing .csv files for the edges and the node attributes.

    Parameters
    ----------
    inputpath : string
        path to the folder containing the .csv files. In the folder, there should be one file for the edges and one file for each node attribute. The files should be named as follows:
            ID_Cavename_Subset_edges.csv: list of links between the nodes. [ID_from,ID_to]
            ID_Cavename_Subset_node_pos: list of the nodes ID with x,y,z coordinates. [ID,x,y,z]
            (Optional) ID_Cavename_Subset_node_csdim: list of the nodes ID with cross-sectional width dimension. [ID,CS_height,CS_width]


    node_attributes: list of string
        list of the node attribute names to attach to the graph (details about the node attributes can be found in the documentation of the KNdata-public repository).
        by default: ['pos','csdim'] 
            'pos': list of 3 floats, [easting,northing,elevation], 3D coordinates of the nodes in specific coordinate system
            'csdim': list of 2 floats, [Width,Height] cross-sectional dimensions of the nodes in m


    Returns
    -------
    networkx graph
    """

    G = nx.Graph()

    

    for file in os.listdir(inputpath):
        sep='' if inputpath.endswith('/') else '/'

        # load EDGES
        if file.endswith('edges.csv'):
            print('loading', file)
            df = pd.read_csv(f'{inputpath}{sep}{file}', delimiter=';')
            #create the graph with the edges
            G.add_edges_from(zip(df.from_id,df.to_id))

        # load node attributes
        else:
            attribute_name = file.split('.')[0].split('_')[-1]
            df = pd.read_csv(f'{inputpath}{sep}{file}', delimiter=';')

            if attribute_name in node_attributes and attribute_name == 'pos':
                print(f'loading {file}')
                #use dict(zip()) when unique entries per node
                dict_attribute = dict(zip(df.id,zip(df.x,df.y,df.z))) 
                nx.set_node_attributes(G,dict_attribute,attribute_name)             

            elif attribute_name in node_attributes and attribute_name == 'csdim':
                print(f'loading {file}')
                dict_attribute = dict(zip(df.id,zip(df.cswidth,df.csheight))) 
                nx.set_node_attributes(G,dict_attribute,attribute_name)          

    return G






# #------------------------------------------------------------------
# #function specific to sql import

def networkx_to_kndata_csv(G, outputpath, make_clean_folder_path = True, cavename='cavename'):
    """_summary_

    Args:
        G (_type_): _description_
        outputpath (_type_): _description_
        make_clean_folder_path (bool, optional): _description_. Defaults to True.
        cavename (str, optional): _description_. Defaults to 'cavename'.
    """
    #check and create new directory if necessary
    if make_clean_folder_path:
        filepath = kn.utils.make_filepath(outputpath,'clean_graph_csv')
    else:
        filepath = outputpath + '/'
    # cavename = G.graph['cavename']
       
    #save edges and flags
    # nx.write_edgelist(G, filepath + cavename + '_edges.csv', data=False)
    pd.DataFrame(G.edges()).to_csv(filepath + cavename + '_edges.csv', index=False, header = ['from_id','to_id'], sep=';', encoding = 'utf-8-sig')    
    print('saved edge list to ', filepath , cavename , '_edges.csv')
    
    #load and save all the existing NODE attribute names for this graph
    attribute_names = kn.utils.get_node_attribute_names(G)
    if len(attribute_names)==1:
        attribute_names = [attribute_names]
    for attribute_name in attribute_names:        
        df = kn.utils.attribute_dict_to_df(G,attribute_name, 'node')
        if attribute_name=='pos':
            df.to_csv(filepath + cavename + '_node_' + attribute_name + '.csv', index=False, header=['id','x','y','z'], sep=';', encoding = 'utf-8-sig')   
        elif attribute_name=='csdim':
            df.to_csv(filepath + cavename + '_node_' + attribute_name + '.csv', index=False, header=['id','cswidth','csheight'], sep=';', encoding = 'utf-8-sig')   



  
        print('saved edge attribute', attribute_name, ' to ', filepath, cavename, '_edge_', attribute_name, '.csv')



