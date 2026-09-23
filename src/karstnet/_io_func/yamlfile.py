#%%
import networkx as nx
import yaml
import karstnet as kn



#!!! to do: transform to list all the tuples
def networkx_to_yaml(G,filepath):
    """Saves graph to a yaml file. This is a compact way to store the graph, 
    but still save it in a text readable format.

        structure of the dictionnary:
        {metadata:  {cavename: string
                     crs: string
                     original_data_rights: string
                     citation: string
                     ...}
        edges: list of list of edge
        
        edge_attributes:
            {flag:list of string
            comment:string
            ...}
        node_attributes:
            {'pos': list of float
            'fulladdress': list of string
            'idsql': list of int
            'flag': list of string
            'comment': list of string
            'splays': list of list of float
            'csdim': list of float
            ...}}}

    Parameters
    ----------
    G : networkx graph
        The networkx graph to save in yaml format
    filepath : string
        Path to save the file, with our without the .yaml extension
    """


    # initiate dictionnary
    dict_graph=dict()
    # extract metadata dictionnary attached to the graph:
    dict_graph['metadata'] = G.graph

    # exctract and attache edges (and make sure its a list of list, instead of a list of tuple)
    # this looks betters in the yaml file
    dict_graph['edges'] = list(map(list, list(G.edges())))

    # extract edges attributes

    list_name_edge_attributes = kn.utils.get_edge_attribute_names(G)
    if list_name_edge_attributes:
        #if there is any edges, initiate dictionnary and add to the dictionnary
        dict_graph['edge_attributes']={}
        for attribute_name in list_name_edge_attributes:
            attribute = nx.get_edge_attributes(G,attribute_name)
            # # create list of edges instead of tuples
            # for attribute_key in attribute:
            #     attribute[list(attribute_key)] = attribute.pop(attribute_key)

            #check if node attribute are tuples. if they are transform into lists    
            if type(attribute[list(attribute.keys())[0]]) == tuple:
                for attribute_key in attribute:
                    attribute[attribute_key]=list(attribute[attribute_key])   
            dict_graph['edge_attributes'][attribute_name] = attribute


    # extract all the node attributes
    list_name_node_attributes = list(kn.utils.get_node_attribute_names(G))
    if list_name_node_attributes:
        dict_graph['node_attributes']={}
        for attribute_name in list_name_node_attributes:
            attribute = nx.get_node_attributes(G,attribute_name)
            #check if node attribute are tuples. if they are transform into lists       
            if type(attribute[list(attribute.keys())[0]]) == tuple:
                for attribute_key in attribute:
                    attribute[attribute_key]=list(attribute[attribute_key])     
            dict_graph['node_attributes'][attribute_name] = attribute

    # # # Writing the data to a YAML file
    #check first if the extension exisits in path name. if not, just add it.
    if '.yaml' not in filepath:
        filepath = filepath + '.yaml'

    with open(filepath, 'w') as file:
        yaml.dump(dict_graph, file, default_flow_style=False, sort_keys=False)

    print(f"Graph saved to yaml here: {filepath}")



def networkx_from_yaml(filepath,add_edge_attr=True,add_node_attr=True):
    """_summary_

    Parameters
    ----------
    filepath : string
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

    Returns
    -------
    networkx graph
        Return a graph networkx with all the node and edge attributes selected 
        attached to the graph as attributes.

    """

    with open(filepath) as f:
        # yml = yaml.safe_load(f) 
        yml = yaml.load(f, Loader=yaml.Loader) 

    G = nx.Graph()

    # add metadata
    if 'metadata' in yml.keys():
        G.graph=yml['metadata']

    # add edges
    G.add_edges_from(yml['edges'])

    if add_edge_attr == False:
        pass

    # add edge attributes
    elif add_edge_attr == True:
        if 'edge_attributes' in yml.keys():
            for attribute_name in yml['edge_attributes']:
                nx.set_edge_attributes(G, yml['edge_attributes'][attribute_name],attribute_name)     
    elif add_edge_attr and type(add_edge_attr) == list:
        for attribute_name in add_edge_attr:
            if ('edge_attributes' in yml.keys()) & (attribute_name in yml['edge_attributes'].keys()): 
                nx.set_edge_attributes(G, yml['edge_attributes'][attribute_name],attribute_name)  
            else: 
                print(f'Warning load networkx from yaml: edge attribute {attribute_name} does not exit in graph, and was not added')
    # if only one attribute is entered
    elif add_edge_attr and type(add_edge_attr) == str:
        attribute_name = add_edge_attr
        if ('edge_attributes' in yml.keys()) & (attribute_name in yml['edge_attributes'].keys()): 
            nx.set_edge_attributes(G, yml['edge_attributes'][attribute_name],attribute_name)  
        else: 
            print(f'Warning load networkx from yaml: edge attribute {attribute_name} does not exit in graph, and was not added')

    # add node attributes
    if add_node_attr == False:
        pass

    elif add_node_attr == True:
        if 'node_attributes' in yml.keys():
            # print('adding node attributes')
            for attribute_name in yml['node_attributes']:
                nx.set_node_attributes(G,yml['node_attributes'][attribute_name],attribute_name) 

    elif add_node_attr and type(add_node_attr) == list:
        for attribute_name in add_node_attr:
            if ('node_attributes' in yml.keys()) & (attribute_name in yml['node_attributes'].keys()): 
                nx.set_node_attributes(G, yml['node_attributes'][attribute_name],attribute_name)  
            else: 
                print(f'Warning load networkx from yaml: node attribute {attribute_name} does not exit in graph, and was not added')
    # if only one attribute is entered
    elif add_node_attr and type(add_node_attr) == str:
        if ('node_attributes' in yml.keys()) & (attribute_name in yml['node_attributes'].keys()): 
            attribute_name = add_node_attr
        else: 
            print(f'Warning load networkx from yaml: node attribute {attribute_name} does not exit in graph, and was not added')

        nx.set_node_attributes(G, yml['node_attributes'][attribute_name],attribute_name)  
    return G 



# %%
