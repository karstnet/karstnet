# list of fictif surveys to represent the connected components (dates enable color visualisation)
# probably a loop

# raw survey data, flagged as duplicate and without any date, to appear gray and dashed.

#%%
# import numpy as np
# import pandas as pd
import networkx as nx
import datetime
# import os
import subprocess
# from os.path import normpath
import karstnet as kn

#!!! todo
# def write_graph_to_plt(graph: nx.Graph, basename: str, disco_keys=None):
#     """
#     Reads a networkx.Graph object and writes it to a plt file.
#     developped by Tanguy Racine, UNINE
#     """
#     coords = nx.get_node_attributes(graph, "coord")
#     plt_command = """{type} {y} {x} {z} S{station} P -9 -9 -9 -9"""
#     plt_file = """G18S\nNX D 1 1 1 C{name}.plt\n{points}"""
#     plt_lines = []

#     for line in graph.edges:
#         x1, y1, z1 = coords[line[0]]
#         x2, y2, z2 = coords[line[1]]

#         plt_lines.append(plt_command.format(
#             type="M", x=x1, y=y1, z = z1, station=line[0]))

#         plt_lines.append(plt_command.format(
#             type="D", x=x2, y=y2, z = z2, station=line[1]))
        
#     if disco_keys != None:
#         for key in disco_keys:
#             x, y, z = coords[key]

#         plt_lines.append(plt_command.format(
#             type="M", x=x, y=y, z=z, station=key))            

#     with open(f"{basename}.plt", "w+") as f:
#         f.write(
#             plt_file.format(
#                 name=basename,
#                 points="\n".join(plt_lines),
#             )
#         # )


def networkx_to_avend3d(G,H, outputpath, cavename):
   #nx2aven3d
    """G the original graph, H the clean graph

    Returns
    -------
    _type_
        _description_
    """
    pos3d = nx.get_node_attributes(G,'pos')

    # if manual_edges_file:
    #     #Extract data from manual edge file
    #     # after manual observations, you can add or remove edges
    #     # create a text file that contains 4 columns ['from', 'to', 'action','comments]
    #     ####################################################################

    #     #read file
    #     manual_edges =  pd.read_csv(manual_edges_file, delimiter=';')
    #     edges_to_add = manual_edges.loc[manual_edges['action']=='add'][['from', 'to']].apply(tuple, axis=1)
    #     edges_to_remove = manual_edges.loc[manual_edges['action']=='remove'][['from', 'to']].apply(tuple, axis=1)
    #     #extract comments and flags to add to the graph. this will help identifying nodes that are added
    #     #comments = dict(zip(edges_to_add,manual_edges.loc[manual_edges['action']=='add', 'comments']))
    #     flags = dict(zip(edges_to_add, ['manual_input']*len(edges_to_add)))
    #     reconnected_nodes_id = list(np.concatenate((np.array(manual_edges.loc[manual_edges['action']=='add']['from']),np.array(manual_edges.loc[manual_edges['action']=='add']['to'])), axis=None))

    #     #check if the edge to add is between already existing nodes. if not, remove the edge and give a warning
    #     for i,[node1,node2] in edges_to_add.items():
    #         if H.has_node(node1)==False or  H.has_node(node2)==False:
    #             print('one of the node of the edge to add is not in the graph without duplicates')
    #             print(node1, 'is in cleaned graph: ', H.has_node(node1))
    #             print(node2, 'is in cleaned graph: ', H.has_node(node2))        
    #             print('dropped the tuple', node1,node2)
    #             edges_to_add = edges_to_add.drop(i)

    #     #add something similar to check if the edges to remove exist??

    #     #Add and remove edges from H
    #     H.add_edges_from(edges_to_add)#, label='manually added')
    #     #nx.set_edge_attributes(H, comments, 'comments') 
    #     nx.set_edge_attributes(H, flags, 'flags')
    #     H.remove_edges_from(edges_to_remove)#, label='manually added')
    #     H.remove_nodes_from(list(nx.isolates(H))) 

    #     #RECONNECTED EDGES AND DISCONNECTED EDGES?
    #     strings_reconnected_edges = []
    #     for edge in edges_to_add:
    #         strings_reconnected_edges.append('RN_%s RN_%s'%(edge[0],edge[1]))

    #     #  RECONNECTED NODES
    #     ###################################
    #     strings_reconnected_nodes = []

    #     for id in reconnected_nodes_id:
    #         #create a string with the coordinates from the points
    #         position = " ".join([str(element) for element in pos3d[id]])
    #         #create each line of text -- fix DN_id x y z --
    #         strings_reconnected_nodes.append('fix RN_%s %s'%(id, position) )
            

    # else: 
    #     strings_reconnected_edges = ['']
    #     strings_reconnected_nodes = ['']   
    
    #DISCONNECTED NODES
    ######################
    # portion of the text containing the point for the disconnected nodes
    #comment how I find disconnected nodes. what is the logic
    disconnected_nodes_id = kn.clean.find_disconnected_node(G, H)
    if disconnected_nodes_id is not None:
        strings_disconnected_nodes = []

        for id_node in disconnected_nodes_id:
            #create a string with the coordinates from the points
            position = " ".join([str(element) for element in pos3d[id_node]])
            #create each line of text -- fix DN_id x y z --
            strings_disconnected_nodes.append(f'fix DN_{id_node} {position}' )
    else:
        strings_disconnected_nodes = ['']

    #SPLAYS
    ########################
    pos_splay = []
    edge_splay = []
    dict_splays = nx.get_node_attributes(G,'splays')
    for key,list_splays in dict_splays.items():
        for i,splay in enumerate(list_splays):
            pos_splay.append(f'fix splay-{key}-{i} {splay[0]} {splay[1]} {splay[2]}')
            edge_splay.append(f'{key} splay-{key}-{i}')
            



    #ENTIRE SURVEY
    ########################
    strings_all_nodes = []
    for id_node in G.nodes:
        position = " ".join([str(element) for element in pos3d[id_node]])
        strings_all_nodes.append(f'fix {id_node} {position}')

    strings_all_edges = []
    for edge in G.edges:
        strings_all_edges.append(f'{edge[0]} {edge[1]}')

    # text containing a list of survey, numbered from 1 to n, based on the connected components
    # the fake date enable coloring of the connected component in aven

    #CONNECTED COMPONENTS FROM CLEAN SURVEY
    #############################################

    text_clean_survey = []
    # dict_address = nx.get_node_attributes(G,'fulladdress')
    for i,cc in enumerate(nx.connected_components(H)):
        Hsub = H.subgraph(cc)

        strings_nodes = []
        for id_node in Hsub.nodes:
            position = " ".join([str(element) for element in pos3d[id_node]])
            strings_nodes.append(f'fix {id_node} {position}')

        strings_edges = []
        for edge in Hsub.edges:
            strings_edges.append(f'{edge[0]} {edge[1]}')

        text_clean_survey.extend([f'survey {i}',
                                'centreline',
                                f'date {(datetime.date(1000,1,1) + datetime.timedelta(days = i)).strftime("%Y.%m.%d")}'])
        text_clean_survey.extend(strings_nodes)
        text_clean_survey.extend(['',
                                'data nosurvey from to'])
        text_clean_survey.extend(strings_edges)
        text_clean_survey.extend(['endcentreline',
                                  'endsurvey',
                                  '',
                                  ''])
        
    # #attempt to display fulladdress, but not working
    # text_clean_survey = []
    # dict_address = nx.get_node_attributes(G,'fulladdress')
    # for i,cc in enumerate(nx.connected_components(H)):
    #     Hsub = H.subgraph(cc)

    #     strings_nodes = []
    #     for id_node in Hsub.nodes:
    #         position = " ".join([str(element) for element in pos3d[id_node]])
    #         strings_nodes.append(f'fix {dict_address[id_node][0]} {position}')

    #     strings_edges = []
    #     for edge in Hsub.edges:
    #         strings_edges.append(f'{dict_address[edge[0]][0]}{ dict_address[edge[1]][0]}')

    #     text_clean_survey.extend([f'survey {i}',
    #                             'centreline',
    #                             f'date {(datetime.date(1000,1,1) + datetime.timedelta(days = i)).strftime('%Y.%m.%d')}'])
    #     text_clean_survey.extend(strings_nodes)
    #     text_clean_survey.extend(['',
    #                             'data nosurvey from to'])
    #     text_clean_survey.extend(strings_edges)
    #     text_clean_survey.extend(['endcentreline',
    #                               'endsurvey',
    #                               '',
    #                               ''])



    main_text = []
    main_text.extend(['survey '+ cavename]                  
                    + ['']
                    + ['#this is a python compilation of a therion file to plot connected components and disconnected nodes']
                    + [ '', ''] 
                    #DISCONNECTED NODES
                    + ['centreline'] 
                    + strings_disconnected_nodes
                    + ['endcentreline']
                    + [ '', ''] 
                    # #RECONNECTED EDGES (GREY - NODATE)
                    # + ['survey reconnection', 'centreline']
                    # + strings_reconnected_nodes
                    # + ['' + 'data nosurvey from to']
                    # + ['flags splay']
                    # + strings_reconnected_edges
                    # + ['flags not splay']
                    # + ['endcentreline']
                    # + ['endsurvey']
                    # + [ '', '']
                    #CONNECTED COMPONENTS (COLORED BY DATES, SPLAY)
                    + text_clean_survey
                    #ALL DATA BACKGROUND (GREY - NODATE, DUPLICATE)
                    + ['survey all_data']
                    + ['centreline']
                    + strings_all_nodes
                    + pos_splay
                    + ['' + 'data nosurvey from to']
                    #background full survey before removing duplicates
                    + ['flags duplicate']
                    + strings_all_edges
                    + ['flags not duplicate']
                    # splays
                    + [ '']
                    + ['flags splay']
                    + edge_splay
                    + ['flags not splay']
                    + ['endcentreline']
                    + ['endsurvey']
                    + [ '', '']
                    + ['endsurvey'])

    # create .th file to run in therion to export .3d file    

    # sometimes the folder path ends with a slash and sometimes not. this way, it does not matter what the person inputs
    # even better, we should write a code that creates a folder ??

    # #check and create new directory if necessary
    # sep = '' if outputpath.endswith('/') else '/'
    # filepath = outputpath + sep + cavename 
    # isExist = os.path.exists(filepath)
    # if not isExist:
    #    os.makedirs(filepath)

    filepath = kn.utils.make_filepath(outputpath,'visualization/connectected_components_3d')


    # if outputpath[-1] == '/':
    #     separation = ''
    # else:
    #     separation = '/'

    with open(filepath + cavename + '_CC3d.th', 'w') as f:
        f.write('\n'.join(main_text))

    with open(filepath + cavename +'_CC3d.thconfig', 'w') as f:
        f.write('\n'.join(['source ' + cavename + '_CC3d.th', 
                        'export model -o ' + cavename + '_CC.3d']))



    thconfig_filename = '"' + filepath + cavename + '_CC3d.thconfig"'
    
    #chemin = "C:/Users/celia/OneDrive - unine.ch/in_progress/2_SystemMigovec/networkx_export/3d_visualization/system_migovec_CC3d.thconfig"
    #path_process = 'C:/Users/celia/github/erc-karst-repositories/networks_code_repo'
    # test_chemin = "therion '%s'"%path_process
    print(f"therion {thconfig_filename}")
    # print(normpath(test_chemin))
    log = subprocess.check_output(f"therion {thconfig_filename}", shell=True)

    #export file with original station stations name
    # for key in H.nodes:
    
    #     assert H.nodes('pos')[key] != 'NoneType', f'{key} has no position'

    #return the cleanned graph with flagged additional edges
    return H


#--------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------

def networkx_to_aven3dclean(G,filepath, cavename):
    #nx2clean3d
    print(f'running nx2clean3d')

    import subprocess
    #from os.path import normpath


    pos3d = nx.get_node_attributes(G,'pos')
    node_pos_all = []

    for id in G.nodes:
        #add station name and coordinates
        position = " ".join([str(element) for element in pos3d[id]])
        node_pos_all.append('fix %s %s'%(id, position))

    duplicate_edges = []
    surface_edges = []
    station_edges = []
    for edge in G.edges(data=True):
        
        if edge[2] == {}:
            station_edges.append('%s %s'%(edge[0],edge[1]))
        
        elif 'srf' in list(edge[2].values())[0]: #edge[2] == {'flags': 'srf'}:
            surface_edges.append('%s %s'%(edge[0],edge[1]))
            
        elif ('dpl' in list(edge[2].values())[0]) or ('rmv' in list(edge[2].values())[0]): #edge[2] == {'flags': 'dpl'}:
            duplicate_edges.append('%s %s'%(edge[0],edge[1]))
        else:
            #for all the other flags that are not listed
            station_edges.append('%s %s'%(edge[0],edge[1]))
        

    splay_edges = []
    for key,items in dict(G.nodes('splays')).items():
        #print(key)    
        if items is not None:
            i = 0
            #for item in items:   !!! toporobot export is list in list, but not here.. fix this         
            for pos in items:
                i=i+1
                # print(key)
                # print(pos)
                splay_name = 'splay_%s_%s'%(i,key)
                position = " ".join([str(element) for element in pos])
                node_pos_all.append('fix %s %s'%(splay_name, position))
            
                splay_edges.append('%s %s'%(key,splay_name))


        
    # strings_all_nodes = []
    # for id in G.nodes:
    #     position = " ".join([str(element) for element in pos3d[id]])
    #     strings_all_nodes.append('fix %s %s'%(id, position))

    # #ENTRANCES
    # entrances = []
    # for key,items in dict(G.nodes('flag')).items():   
    #     if 'ent' in items:
    #         entrances.append(f'station {key} "{key}" entrance')


    #write text file
    text=[]
    text.extend(['survey ' + cavename]+
                ['']+
                ['']+
                ['survey plot']+
                ['centreline']+
                ['']+
                node_pos_all+
                ['']+
                ['data nosurvey from to']+
                station_edges+
                ['']+
                ['flags splay']+
                splay_edges+
                ['flags not splay']+
                ['']+
                ['flags duplicate']+
                duplicate_edges+
                ['flags not duplicate']+
                ['']+
                ['flags surface']+
                surface_edges+
                ['flags not surface']+
                # ['']+
                # entrances+
                ['']+
                ['endcentreline']+
                ['endsurvey']+
                ['']+
                ['endsurvey'])


    #create output folder if it does not exist yet

    print(f'saving th files here: {filepath}{cavename}_clean.th')
    with open(f'{filepath}{cavename}_clean.th', 'w') as f:
        f.write('\n'.join(text))

    with open(f'{filepath}{cavename}_clean.thconfig', 'w') as f:
        f.write('\n'.join([f'source {cavename}_clean.th', 
                        f'export model -o  {cavename}_clean.3d',
                        f'export model -o  {cavename}_clean.lox',
                        # f'export model -o  output/{cavename}_clean.kml',
                        f'export model -o  {cavename}_shapefiles_clean -fmt esri '
                        ]))



    thconfig_filename = f'"{filepath}{cavename}_clean.thconfig"'

    #chemin = "C:/Users/celia/OneDrive - unine.ch/in_progress/2_SystemMigovec/networkx_export/3d_visualization/system_migovec_CC3d.thconfig"
    #path_process = 'C:/Users/celia/github/erc-karst-repositories/networks_code_repo'
    # test_chemin = "therion '%s'"%path_process
    print("therion %s"%thconfig_filename)
    # print(normpath(test_chemin))
    log = subprocess.check_output("therion %s"%thconfig_filename, shell=True)


#--------------------------------------------------------
# EXPORT TO .3D
##############
def nx2section3d(G,outputpath,list_sections=[]):

    import subprocess
    #from os.path import normpath

    #extract information about the different sections of the cave
    sections = nx.get_node_attributes(G,'section')

    cavename = cavename

    pos3d = nx.get_node_attributes(G,'pos')
    
    survey_text = []
    
    # unique_section_id = list(set(list(nx.get_node_attributes(G,'section').values())))

    

    for section_id in list_sections:
        print(f'section {section_id}')
        list_section_keys = [k for k, v in sections.items() if section_id in v]
        #print('length list_section_keys=',len(list_section_keys))
        H = G.subgraph(list_section_keys)
        node_pos_all = []
        
        for id in H.nodes:
            
            #add station name and coordinates
            position = " ".join([str(element) for element in pos3d[id]])
            node_pos_all.append('fix %s %s'%(id, position))

        duplicate_edges = []
        surface_edges = []
        station_edges = []
        for edge in H.edges(data=True):
            
            if edge[2] == {}:
                station_edges.append('%s %s'%(edge[0],edge[1]))
            
            elif 'srf' in list(edge[2].values())[0]: #edge[2] == {'flags': 'srf'}:
                surface_edges.append('%s %s'%(edge[0],edge[1]))
                
            elif 'dpl' in list(edge[2].values())[0]: #edge[2] == {'flags': 'dpl'}:
                duplicate_edges.append('%s %s'%(edge[0],edge[1]))
            

        splay_edges = []
        for key,items in dict(H.nodes('splays')).items():
            #print(key)    
            if items is not None:
                i = 0
                #for item in items:   !!! toporobot export is list in list, but not here.. fix this         
                for pos in items:
                    i=i+1
                    # print(key)
                    # print(pos)
                    splay_name = 'splay_%s_%s'%(i,key)
                    position = " ".join([str(element) for element in pos])
                    node_pos_all.append('fix %s %s'%(splay_name, position))
                
                    splay_edges.append('%s %s'%(key,splay_name))


            
        # strings_all_nodes = []
        # for id in H.nodes:
        #     position = " ".join([str(element) for element in pos3d[id]])
        #     strings_all_nodes.append('fix %s %s'%(id, position))

        # #ENTRANCES
        # entrances = []
        # for key,items in dict(H.nodes('flag')).items():   
        #     if 'ent' in items:
        #         entrances.append(f'station {key} "{key}" entrance')


        #write text file
        survey_text.extend(
                ['']+
                [f'survey {section_id}']+
                ['centreline']+
                ['']+
                node_pos_all+
                ['']+
                ['data nosurvey from to']+
                station_edges+
                ['']+
                ['flags splay']+
                splay_edges+
                ['flags not splay']+
                ['']+
                ['flags duplicate']+
                duplicate_edges+
                ['flags not duplicate']+
                ['']+
                ['flags surface']+
                surface_edges+
                ['flags not surface']+
                # ['']+
                # entrances+
                ['']+
                ['endcentreline']+
                ['endsurvey']
                )
        
    text = []
    text.extend([f'survey {cavename}']+
        ['']+
        survey_text+
        ['']+
        ['endsurvey'])

    

    filepath = kn.utils.make_filepath(outputpath,'visualization/sections')

    print('write .th file')
    with open(filepath + cavename + '_sections.th', 'w') as f:
        f.write('\n'.join(text))
    print('write .thconfig')
    with open(filepath + cavename + '_sections.thconfig', 'w') as f:
        f.write('\n'.join(['source ' + cavename + '_sections.th', 
                        'export model -o ' + cavename + '_sections.3d']))

    print('compile therion project')
    thconfig_filename = '"' + filepath + cavename + '_sections.thconfig"'

    #chemin = "C:/Users/celia/OneDrive - unine.ch/in_progress/2_SystemMigovec/networkx_export/3d_visualization/system_migovec_CC3d.thconfig"
    #path_process = 'C:/Users/celia/github/erc-karst-repositories/networks_code_repo'
    # test_chemin = "therion '%s'"%path_process
    print("therion %s"%thconfig_filename)
    # print(normpath(test_chemin))
    log = subprocess.check_output("therion %s"%thconfig_filename, shell=True)


    ##################
    ##################


# %%
