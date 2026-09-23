import pandas as pd
import os
import networkx as nx
import karstnet as kn


def networkx_from_kndata_csv(inputpath,
                         node_attributes=['pos','csdim']):
    """Loads the cave graph from the Github repository erc-karst-repositories/networks_datasets

    Parameters
    ----------
    inputpath : string
        path to the folder containing the .csv files 
    node_attributes: list of string
        list of the node attribute names to attach to the graph
        by default: ['pos','csdim']

        

    Returns
    -------
    networkx graph
        clean graph of the cave
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



  
        # print('saved edge attribute', attribute_name, ' to ', filepath, cavename, '_edge_', attribute_name, '.csv')



# #----------------------------
# #not sure this one works



# def load_clean_graph_csv(inputpath,
#                          node_attributes=['pos','csdim','flags','fulladdress','idsql','splays','comments'], 
#                          edge_attributes=['flags']):
#     """Loads the cave graph from the Github repository erc-karst-repositories/networks_datasets

#     Parameters
#     ----------
#     inputpath : string
#         path to the folder containing the .csv files 
#     node_attributes: list of string
#         list of the node attribute names to attach to the graph
#         by default: ['pos','csdim','flags','fulladdress','idsql','splays','comments']
#     edge_attributes: list of string
#         list of the edge attribute names to attach to the graph
#         by default: ['flags']

        

#     Returns
#     -------
#     networkx graph
#         clean graph of the cave
#     """

#     G = nx.Graph()

    

#     for file in os.listdir(inputpath):
#         sep='' if inputpath.endswith('/') else '/'

#         # load EDGES
#         if file.endswith('edges.csv'):
#             print('loading', file)
#             df = pd.read_csv(f'{inputpath}{sep}{file}', delimiter=';')
#             #create the graph with the edges
#             G.add_edges_from(zip(df.from_id,df.to_id))

#         # Load Edges flags (if present)
#         elif file.endswith('edge_flags.csv') and 'flags' in edge_attributes:
#             print('loading',file)
#             df = pd.read_csv(f'{inputpath}{sep}{file}', delimiter=';')
#             #create the graph with the edges
#             dict_attribute = kn.utils.list2dict(list(zip(df.from_id,df.to_id)),df.flags)
#             nx.set_edge_attributes(G, dict_attribute,'flags')

#         # load node attributes
#         else:
#             attribute_name = file.split('.')[0].split('_')[-1]
#             df = pd.read_csv(f'{inputpath}{sep}{file}', delimiter=';')

#             if attribute_name in node_attributes and attribute_name == 'pos':
#                 print(f'loading {file}')
#                 #use dict(zip()) when unique entries per node
#                 dict_attribute = dict(zip(df.id,zip(df.x,df.y,df.z))) 
#                 nx.set_node_attributes(G,dict_attribute,attribute_name)             

#             elif attribute_name in node_attributes and attribute_name == 'csdim':
#                 print(f'loading {file}')
#                 dict_attribute = dict(zip(df.id,zip(df.cswidth,df.csheight))) 
#                 nx.set_node_attributes(G,dict_attribute,attribute_name) 

#             elif attribute_name in node_attributes and attribute_name == 'flags':
#                 print(f'loading {file}')
#                 #use kn.utils.list2dict() when mulitple entries per nodes
#                 dict_attribute = kn.utils.list2dict(df.id,df.flags)   
#                 nx.set_node_attributes(G,dict_attribute,attribute_name)

#             elif attribute_name in node_attributes and attribute_name == 'fulladdress':
#                 print(f'loading {file}')
#                 dict_attribute = kn.utils.list2dict(df.id,df.fulladdress)  
#                 nx.set_node_attributes(G,dict_attribute,attribute_name) 

#             elif attribute_name in node_attributes and attribute_name == 'idsql':
#                 print(f'loading {file}')
#                 dict_attribute = kn.utils.list2dict(df.id,df.idsql)   
#                 nx.set_node_attributes(G,dict_attribute,attribute_name)

#             elif attribute_name in node_attributes and attribute_name == 'splays':
#                 print(f'loading {file}')
#                 dict_attribute = kn.utils.list2dict(df.id,list(zip(df.x,df.y,df.z))) 
#                 nx.set_node_attributes(G,dict_attribute,attribute_name)  

#             elif attribute_name in node_attributes and attribute_name == 'comments':
#                 print(f'loading {file}')
#                 dict_attribute = kn.utils.list2dict(df.id,df.comment)   
#                 nx.set_node_attributes(G,dict_attribute,attribute_name)
        
              

#     return G






# # #------------------------------------------------------------------
# # #function specific to sql import

# def save_attribrutes_df_to_csv(G, outputpath, cavename='cavename'):
#     #check and create new directory if necessary
#     filepath = kn.utils.make_filepath(outputpath,'clean_graph_csv')
#     # cavename = G.graph['cavename']
       
#     #save edges and flags
#     # nx.write_edgelist(G, filepath + cavename + '_edges.csv', data=False)
#     pd.DataFrame(G.edges()).to_csv(filepath + cavename + '_edges.csv', index=False, header = ['from_id','to_id'], sep=';')    
#     print('saved edge list to ', filepath , cavename , '_edges.csv')
    
#     #load and save all the existing NODE attribute names for this graph
#     attribute_names = kn.utils.get_node_attribute_names(G)
#     if len(attribute_names)==1:
#         attribute_names = [attribute_names]
#     for attribute_name in attribute_names:        
#         df = kn.utils.attribute_dict_to_df(G,attribute_name, 'node')
#         if attribute_name=='pos' or attribute_name=='splays' or attribute_name=='splaylegs' or attribute_name=='splaylrud' :
#             df.to_csv(filepath + cavename + '_node_' + attribute_name + '.csv', index=False, header=['id','x','y','z'], sep=';')   
#         elif attribute_name=='csdim':
#             df.to_csv(filepath + cavename + '_node_' + attribute_name + '.csv', index=False, header=['id','cswidth','csheight'], sep=';')   
#         elif attribute_name=='flags':
#             df.to_csv(filepath + cavename + '_node_' + attribute_name + '.csv', index=False, header=['id','flag'], sep=';')   
#         elif attribute_name=='comments':
#             df.to_csv(filepath + cavename + '_node_' + attribute_name + '.csv', index=False, header=['id','comment'], sep=';')   
#         elif attribute_name=='idsql':
#             df.to_csv(filepath + cavename + '_node_' + attribute_name + '.csv', index=False, header=['id','idsql'], sep=';')    
#         elif attribute_name=='fulladdress':
#             df.to_csv(filepath + cavename + '_node_' + attribute_name + '.csv', index=False, header=['id','fulladdress'], sep=';')    
#         # print('saved node attribute', attribute_name, ' to ', filepath, cavename, '_node_' , attribute_name , '.csv')
    
#     #load and save all the existing EDGE attribute names for this graph
#     attribute_names = kn.utils.get_edge_attribute_names(G)
#     for attribute_name in attribute_names:
#         df = kn.utils.attribute_dict_to_df(G,attribute_name, 'edge')
#         if attribute_name=='flags':
#             df.to_csv(filepath + cavename + '_edge_' + attribute_name + '.csv', index=False, header=['from_id','to_id','flag'], sep=';')   
  
#         # print('saved edge attribute', attribute_name, ' to ', filepath, cavename, '_edge_', attribute_name, '.csv')



# # #----------------------------
# # #not sure this one works

