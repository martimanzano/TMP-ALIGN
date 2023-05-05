from networkx import DiGraph, ancestors, write_dot, read_dot, set_node_attributes
from networkx import *
import PathSummary
import subprocess
import settings

__author__ = 'marti'


def hypergraph_to_reaction_graph(rpath):
    reaction_graph = DiGraph()                          # reactions graph: simple digraph (without hyperedges)
    for node in rpath.dg.nodes():
        if 'n' in node:                             # iter through original hyperedge's nodes
            process_intermediate_node(node, rpath.dg, reaction_graph)
    rpath.reaction_graph = reaction_graph               # save it to the path object


def process_intermediate_node(int_node, dg, reaction_graph):
    for inc_reaction in dg.predecessors(int_node):      # hyperedeges incidents to THIS node
        if not reaction_graph.has_node(haTOn(inc_reaction)):
            reaction_graph.add_node(haTOn(inc_reaction))
        for out_reaction in dg.neighbors(int_node):     # hyperedges accesible from THIS node
            process_out_reaction(haTOn(out_reaction), haTOn(inc_reaction), reaction_graph)

def process_out_reaction(out_reaction, inc_reaction, reaction_graph):
    """
    Function that adds two reactions as nodes to the new digraph, then adds an edge to connect them
    """
    if not reaction_graph.has_node(out_reaction):
        reaction_graph.add_node(out_reaction)
    if not reaction_graph.has_edge(inc_reaction, out_reaction) and not isReversedEdge(inc_reaction, out_reaction):
        reaction_graph.add_edge(inc_reaction, out_reaction)

def isReversedEdge(inc_reaction, out_reaction):
    if inc_reaction == str(out_reaction + "rev") or str(inc_reaction+"rev") == str(out_reaction):
    #if inc_reaction == str(out_reaction + "rev"):
        return True
    return False

def haTOn(ha):
    return ha.replace("he","")

global_branch_index = 0

def reaction_graph_to_path_tree(rpath):
    """
    Function that converts an input reaction graph to a tree graph containing all the longest paths.
    First of all, the reaction graph is converted to a DAG to clean cycles.
    If there are more than one initial node, then the resulting graph is a forest, not a tree
    """
    longest_paths_list = []
    paths_tree = DiGraph()
    rpath.non_accessible_reactions = rpath.react_dict.copy()
    try:
        global global_branch_index
        global_branch_index = 0
        initial_nodes = []
        react_graph = rpath.reaction_graph
        #react_graph = test_four_complete2x()
        #graph_to_dot(react_graph, rpath.results_path, "test")
        initial_nodes = get_initial_nodes(react_graph)                       # Nodes without predecessors
        initial_nodes_attribute = ''

        for i_node in initial_nodes:
            initial_nodes_attribute += str(i_node) + '_' + str(global_branch_index) + ','
            visit_rg_node(react_graph, paths_tree, i_node, [], rpath.summary, longest_paths_list, rpath.non_accessible_reactions)
            #global_branch_index += 1
    except ValueError as x:
        print x.message
        rpath.summary.paths_limit_raised = True
    rpath.summary.pathsNumber = len(longest_paths_list)                  # or global_branch_list - 1
    rpath.longest_paths = longest_paths_list
    print "Total paths number: " + str(rpath.summary.pathsNumber)
    return paths_tree

def reaction_graph_to_longest_paths(rpath):
    longest_paths_list = []
    rpath.non_accessible_reactions = rpath.react_dict.copy()
    try:
        global global_branch_index
        global_branch_index = 0
        react_graph = rpath.reaction_graph
        initial_nodes = get_initial_nodes(react_graph)                        # Nodes without predecessors
        for i_node in initial_nodes:
            visit_rg_node_lite(react_graph, i_node, [], rpath.summary, longest_paths_list, rpath.non_accessible_reactions)
    except ValueError as x:
        rpath.summary.paths_limit_raised = True
        print x.message
    rpath.summary.pathsNumber = len(longest_paths_list)                  # or global_branch_list - 1
    print "Total paths number: " + str(rpath.summary.pathsNumber)
    return longest_paths_list

def get_initial_nodes(reaction_graph):
    # Initial node: node with input degree = 0. With the exception of initial # nodes that have a reversed reaction.
    # Reversed reaction's associated edges can be initial nodes if it's allowed in settings file.
    initial_nodes = []
    for node in reaction_graph.nodes():
        valid = True
        #if (type(node) is str and not "rev" in node) or (type(node) is int):
        if ("rev" in node and settings.reversed_reactions_canbe_initial_nodes) or (not "rev" in node):
                for incident in reaction_graph.predecessors(node):
                    if ((type(node) is str) and not (incident == str(node + "rev"))) or (type(node) is int):
                        valid = False
                        break
        else:
            valid = False

        if valid:
            initial_nodes.append(node)
    print "Initial nodes: " + str(initial_nodes)
    return initial_nodes

def get_node_branch(node, branch_index):
    return str(node) + "_" + str(branch_index)

def visit_rg_node(react_graph, paths_forest, node, partial_path, summary, longest_paths_list, non_accessible_reacts):
    global global_branch_index
    if global_branch_index - 1 > settings.MAX_PATHS:
        raise ValueError("Raised path's limit: " + str(settings.MAX_PATHS) + ", consider increasing MAX_PATHS variable or not allowing "
                                                                              "reversed reactions as initial nodes")
    local_branch_index = global_branch_index
    if not paths_forest.has_node(node):
        paths_forest.add_node(get_node_branch(node, local_branch_index))
    valid_successors = set(react_graph.neighbors(node)) - set(partial_path)
    for successor in valid_successors:
        paths_forest.add_node(get_node_branch(successor, global_branch_index))
        paths_forest.add_edge(get_node_branch(node, local_branch_index), get_node_branch(successor, global_branch_index))
        updated_path = partial_path + [node]
        visit_rg_node(react_graph, paths_forest, successor, updated_path, summary, longest_paths_list, non_accessible_reacts)
    if len(valid_successors) == 0:
        partial_path.append(node)
        global_branch_index += 1
        summary.pathsLong.append(len(partial_path))
        longest_paths_list.append(partial_path)
    non_accessible_reacts.pop(node, None)

def visit_rg_node_lite(react_graph, node, partial_path, summary, longest_paths_list, non_accessible_reacts):
    global global_branch_index
    if global_branch_index - 1 > settings.MAX_PATHS:
        raise ValueError("Raised path's limit: " + str(settings.MAX_PATHS) + ", consider increasing MAX_PATHS variable or not allowing "
                                                                              "reversed reactions as initial nodes")
    valid_successors = set(react_graph.neighbors(node)) - set(partial_path)
    for successor in valid_successors:
        updated_path = partial_path + [node]
        visit_rg_node_lite(react_graph, successor, updated_path, summary, longest_paths_list, non_accessible_reacts)
    if len(valid_successors) == 0:
        global_branch_index += 1
        partial_path.append(node)
        summary.pathsLong.append(len(partial_path))
        longest_paths_list.append(partial_path)
    non_accessible_reacts.pop(node, None)


def graph_to_dot(hgraph, result_path, namefile):
        save_file = result_path + namefile + '.dot'
        with open(save_file, 'w') as f:           # save graph to a .dot file
            try:
                #f.write(write(hgraph))
                write_dot(hgraph, f)
                f.close()
                #FNULL = open(os.devnull, 'w')    # then represent it as a image
                subprocess.Popen(['dot', save_file, '-Tpng', '-o', str(save_file+'.png')], stdout=None,
                                 stderr=None, stdin= None, close_fds= True)
                #os.remove(save_file)
            except IOError as e:
                print "I/O error when saving {2} png image ({0}): {1}".format(e.errno, e.strerror, namefile)
                f.close()

#tests
def test_digraph():
    testdg = DiGraph()
    testdg.add_nodes_from(['42','84','82','79','77'])
    testdg.add_node('84rev')
    testdg.add_node('82rev')
    testdg.add_edge('42','84')
    testdg.add_edge('84','84rev')
    testdg.add_edge('84rev','84')
    testdg.add_edge('84rev','82')
    testdg.add_edge('82','79')
    testdg.add_edge('82','82rev')
    testdg.add_edge('79','84rev')
    testdg.add_edge('79','82rev')
    testdg.add_edge('79','80')
    testdg.add_edge('79','77')
    testdg.add_edge('77','79')
    testdg.add_edge('82rev','79')
    testdg.add_edge('82rev','82')
    testdg.add_edge('82rev','84')
    return testdg

def test_complete():
    testdg = DiGraph()
    testdg.add_nodes_from(['i', '1', '2', '3', '4', '5'])
    testdg.add_edge('i','1')
    testdg.add_edge('1','2')
    testdg.add_edge('1','3')
    testdg.add_edge('1','4')
    testdg.add_edge('1','5')
    testdg.add_edge('2','1')
    testdg.add_edge('2','3')
    testdg.add_edge('2','4')
    testdg.add_edge('2','5')
    testdg.add_edge('3','1')
    testdg.add_edge('3','2')
    testdg.add_edge('3','4')
    testdg.add_edge('3','5')
    testdg.add_edge('4','1')
    testdg.add_edge('4','2')
    testdg.add_edge('4','3')
    testdg.add_edge('4','5')
    testdg.add_edge('5','1')
    testdg.add_edge('5','2')
    testdg.add_edge('5','3')
    testdg.add_edge('5','4')
    return testdg

def test_four_complete2x():
    testdg = DiGraph()
    testdg.add_nodes_from(['i', '1', '2', '3', '4', 'i2', '5', '6', '7', '8'])
    testdg.add_edge('i','1')
    testdg.add_edge('1','2')
    testdg.add_edge('1','3')
    testdg.add_edge('1','4')
    testdg.add_edge('2','1')
    testdg.add_edge('2','3')
    testdg.add_edge('2','4')
    testdg.add_edge('3','1')
    testdg.add_edge('3','2')
    testdg.add_edge('3','4')
    testdg.add_edge('4','1')
    testdg.add_edge('4','2')
    testdg.add_edge('4','3')

    testdg.add_edge('i2','5')
    testdg.add_edge('5','6')
    testdg.add_edge('5','7')
    testdg.add_edge('5','8')
    testdg.add_edge('6','5')
    testdg.add_edge('6','7')
    testdg.add_edge('6','8')
    testdg.add_edge('7','5')
    testdg.add_edge('7','6')
    testdg.add_edge('7','8')
    testdg.add_edge('8','5')
    testdg.add_edge('8','6')
    testdg.add_edge('8','7')

    testdg.add_edge('2', '7')
    return testdg
