import xml.etree.ElementTree as ET
import os, subprocess, sys
import networkx as nx
from graphTools import graph_to_dot
from PathSummary import PathSummary
import settings


class Rpath:
    def __init__(self, name, kgml_file, ecfile, directory):
        self.name = name
        self.file = kgml_file
        self.ec_dict = {}                                       # enzime dictionary id -> name
        self.d_react = []									    # lineal list to save the reactions
        self.d_comp = {}    									# hash list indexed by id's to save the compounds
        self.react_dict = {}                                    # hash list to reaction
        self.d_comp_lin = []									# lineal list of id's to iterate over compounds
        self.ind_react = -1								        # number of reactions and compounds to build the matrices
        self.ind_comp = -1
        self.ecfile = ecfile                                    # enzime's file
        self.results_path = directory
        self.hgfile = self.results_path+'hypergraph_'+name+'.dot'
        self._incidence_matrix = []
        self.invariants = []
        self.longest_paths = []
        self.reaction_graph = nx.DiGraph()
        self.dg = nx.DiGraph()
        self.paths_tree = nx.DiGraph()
        self.summary = PathSummary()
        self.non_accessible_reactions = {}
            
    def parse(self):
        ectree = ET.parse(self.ecfile)
        enz_root = ectree.getroot()
        for enzime in enz_root.iter('entry'):		                # save enzime's name to a hash dict
            self.ec_dict[enzime.get('id')] = enzime.get('name')

        tree = ET.parse(self.file)  # build an element self.tree object from the XML input self.fileetroot()
        root = tree.getroot()
        for root_child in root.iter():
            if root_child.tag == 'entry':
                if root_child.get('type') == 'compound':
                    self.ind_comp += 1
                    name = root_child.get('name')				        # for the compounds only name and id are saved
                    id = root_child.get('id')
                    self.d_comp[id] = (name, self.ind_comp)			    # hashing by id compound
                    self.d_comp_lin.append(id)					        # lineal structure [index] = id to easily iterate through the compounds
                    #self.dg.add_node(str("n"+id))
            elif root_child.tag == 'reaction':
                self.ind_react += 1
                sub_list = []                 
                id = root_child.get('id')						                    # id and name of the reaction
                name = root_child.get('name')
                reversible = (root_child.get('type') == 'reversible')
                hedge_id = str("he"+id)                                             # hedge ID = reaction ID
                if reversible:
                    revhedge_id = hedge_id + 'rev'                                  # reversed reaction
                self.dg.add_node(hedge_id, shape = 'box')
                for sub in root_child.iter('substrate'):		        # sublists of every reaction's substractes and products
                    sub_list.append(sub.get('id'))
                    if not self.dg.has_edge(str('n'+(sub.get('id'))), hedge_id):
                        self.dg.add_edge(str('n'+(sub.get('id'))), hedge_id )       # substractes -> hyperedge
                    if reversible and not self.dg.has_edge(revhedge_id, str('n'+(sub.get('id')))):
                        if not self.dg.has_node(revhedge_id):
                            self.dg.add_node(revhedge_id, shape = 'box')            # reverse reaction's hyperedge
                        self.dg.add_edge(revhedge_id, str('n'+(sub.get('id'))))    # rev. hyperedge -> substractes

                prd_list = []
                for prd in root_child.iter('product'):
                    prd_list.append(prd.get('id'))
                    if not self.dg.has_edge(hedge_id, str('n'+prd.get('id'))):
                        self.dg.add_edge(hedge_id, str('n'+prd.get('id')))          # hyperedge -> products
                    if reversible and not self.dg.has_edge(str('n'+prd.get('id')), revhedge_id):
                        self.dg.add_edge(str('n'+prd.get('id')), revhedge_id)       # products -> rev. hyperedge
                              
                details = (name, id, sub_list, prd_list, self.ec_dict.get(id, "ec:0.0.0.0"))	# all the elements are in a tuple
                self.d_react.append(details)							            # simple list indexed by an integer, not hashing
                self.react_dict[id] = details
                        
                if reversible:                  		# if the reaction is reversible, create and save the reverse reaction
                    self.ind_react += 1
                    details = (str(name+'#rev'), id, prd_list, sub_list, self.ec_dict.get(id, "ec:0.0.0.0")) # substracte and product lists are reversed
                    self.d_react.append(details)					    # save it as a different reaction
                    self.react_dict[id+'rev'] = details

        if len(self.d_react) > 0:
            graph_to_dot(self.dg, self.results_path, 'hypergraph')
            return True
        else:
            print "Detected pathway without reactions. Aborting process for this one..."
            return False

    def incidence_matrix(self):
        input_matrix = [[0 for x in range(self.ind_comp+1)] for y in range(self.ind_react+1)] 	        # files-columnes
        output_matrix = [[0 for x in range(self.ind_comp+1)] for y in range(self.ind_react+1)] 	        # files-columnes
        self._incidence_matrix = [[0 for x in range(self.ind_comp+1)] for y in range(self.ind_react+1)]  # files-columne

        for idx, cols in enumerate(self.d_react):
            for subs in cols[2]:
                index = self.d_comp[str(subs)][1]
                input_matrix[idx][index] = 1
                                                            # index: reaction -> compound (column-row)
            for prods in cols[3]:
                index = self.d_comp[str(prods)][1]
                output_matrix[idx][index] = 1
        
        for r in range(len(self._incidence_matrix)):		# incidence matrix = output matrix - input matrix
            for c in range(len(self._incidence_matrix[r])):
                self._incidence_matrix[r][c] = output_matrix[r][c]-input_matrix[r][c]


    def build_inc_matrix_file(self):
        outfile = open(str(self.results_path+self.name+'.txt'),"wb")	     # same name but .txt extension
        outfile.write(str(len(self._incidence_matrix[0]))+' '+str(len(self._incidence_matrix))+'\n')	# incidence matrix length
            
        for r in range(len(self._incidence_matrix[0])):   # ITERAR PER CADA PRIMER ELEMENT, CADA SEGON DE SA LLISTA, ETC
            outfile.write('\n')
            for c in range(len(self._incidence_matrix)):
                temp = str(self._incidence_matrix[c][r])
                while len(temp) < 2:				     # double spaces for 0 and 1, one space for -1 to keep it alineated
                    temp = " " + temp
                outfile.write(temp + " ")			
                
        outfile.write('\n\nCompounds:')				    # save compound's list
        for idx, cmp in enumerate(self.d_comp_lin):		# accessing the hash list through the lineal structure of cmp's ids
            outfile.write('\n'+str(idx+1)+' '+self.d_comp[str(cmp)][0]+'-'+cmp)	
            
        outfile.write('\n\nTransitions:')			    # finally save reaction's list
        for idx, reac in enumerate(self.d_react):
            if len(reac) > 4:
                outfile.write('\n'+str(idx+1)+' '+reac[0]+'-'+reac[1]+' '+reac[4])	# reaction's structure is lineal so we only need to iterate over it
            else:
                outfile.write('\n'+str(idx+1)+' '+reac[0]+'-'+reac[1])
        outfile.close()

    def invariants_info(self, library_path):
        if settings.align_tinvariants_reactions == True:
            tinvs = []
            return
        if sys.platform.startswith('linux'):
            #print "Linux system detected. Running 4ti2, please wait..."
            hilbert_loc = os.getcwd()+"/"+library_path + "hilbert_linux/hilbert"
        elif sys.platform.startswith('darwin'):
            #print "MAC OS X system detected. Running 4ti2, please wait..."
            hilbert_loc = os.getcwd()+"/"+library_path + "hilbert_macOS/hilbert"
        elif sys.platform.startswith('win32'):
            #print "Windows system detected. Running 4ti2, please wait..."
            hilbert_loc = os.getcwd()+"/"+library_path + "hilbert_win32/hilbert"
        else:
            print("4ti2 library not present for "+sys.platform+" aborting...")
            return 0
        inc_matrix_file_path = str(os.getcwd()+"/"+self.results_path+self.name+'.txt')
        FNULL = open(os.devnull, 'w')
        subprocess.call([hilbert_loc, inc_matrix_file_path], stdout=FNULL, stderr=FNULL)
        # -----------------------synrchronous----------------------------------#
        
        hfilename = str(self.results_path+self.name+'.txt.hil')
        hilfile = open(hfilename, "a+")				# pointer at 0 for reading, current EOF (end of file) for writing
        lines_list = hilfile.readlines()			# save all the self.file's lines
        rows, cols = (int(val) for val in lines_list[0].split())	            # rows and cols are in the first line
        data = [[int(val) for val in line.split()] for line in lines_list[2:]]  # T-INVs begins at row no. 4
            
        tinvs = []
        reversible = False
        for r_index, r in enumerate(data):                              # iterate through all the rows
            temp_inv = []                                               # variable to store indexes of every t-inv
            size_inv = 0                                                # and it's number of elements
            for react,c in enumerate(r):                                # iterate through columns
                if c == 1:                                              # '1' found, there is a t-inv
                    last_reaction = self.d_react[react][0]              # save it's name
                    temp_inv.append(react)                              # and it's index
                    if react+1 < cols:                                  # let's see what's next...
                        if data[r_index][react+1] == 1 and last_reaction+'#rev' == self.d_react[react+1][0]:  # is it a trival t-inv??
                            reversible = True                            # yes it is
                        else:
                            reversible = False                          # or not
                    if not reversible:
                        size_inv += 1                         # if it's not, then it's definetively a non trivial t-inv
            if size_inv > 1:
                tinvs.append((r_index, temp_inv))                      # save all the t-inv and it's row number only if it's not trivial
                        
        hilfile.write('Non-trivial T-invariants: \n')
                                                                # Write them to the inv file and to the self information
        for invariants in tinvs:                                       # iterate for all the t-inv
            hilfile.write('Row '+str(invariants[0]+1)+': ')             
            last_reaction = 'none'
            tinv_elements = []                                         # empty the elements info of every t-inv
            for idx, reactions in enumerate(invariants[1]):            # and iterate through every reaction of the t-inv
                #print self.d_react[reactions][0], last_reaction+'#rev'
                #print self.d_react[reactions]
                tinv_elements.append(self.d_react[reactions])          # save the reaction info
                if self.d_react[reactions][0] == last_reaction+'#rev':  
                    hilfile.write('(or #rev)\n')                         # t-inv formed by a trivial t-inv + a non triv. one
                else:
                    details = self.d_react[reactions][0] + ' ' + self.d_react[reactions][4]
                    if idx > 0:
                        details = ' \\ ' + details + '\n'
                    hilfile.write(details)					           # write the details of the t-inv
                last_reaction = self.d_react[reactions][0]
            self.invariants.append(tinv_elements)                      # and store the t-inv to the self var
        self.summary.noOfInvariants = len(self.invariants)
        hilfile.close()

        #print self.invariants
        # self._expand_all_cycles()                                    # get all the t-inv's combinations
        #print '----------------------------------------------------------------------'
        # print self.expanded_cicles


# --------deprecated section---------
    def invariants_to_sets(self):   # function that transforms our T-invariants to MPAlign ones and returns
        set_list = []               # a list of T-Invariants sets
        for invariant in self.invariants:
            invariant_set = set()
            for reaction in invariant:
                react_name = reaction[0]
                react_id = reaction[1]
                react_enz = reaction[4]
                reacts = react_name.split()
                MPAName = ""
                for idx, react in enumerate(reacts):
                    if "#rev" in react and idx == len(reacts) - 1:
                        react = react.replace("#rev", "")
                        react += "-" + react_id
                        react += "#rev"
                    elif idx == len(reacts) - 1:
                        react += "-" + react_id
                    MPAName += react + " "
                MPAName += react_enz.upper()
                invariant_set.add(MPAName)
            set_list.append(invariant_set)
        return set_list


def find_invariants_paths(invs_set, paths_set, coinc_file):  # search in which path is every T-invariant
    coinc_number = 0
    for inv_set in invs_set:                                 # for every T-invariant
        for path_set in paths_set:                           # compare it to every path until it's found
            if inv_set.issubset(path_set) or path_set.issubset(inv_set):  # is it inside the path?
                coinc_number += 1
                coinc_file.write("Coincidence found: T-Invariant - Path\n")
                coinc_file.writelines(x + ", " for x in list(inv_set))
                coinc_file.write("\n\n")
                coinc_file.writelines(x + ", " for x in list(path_set))
                coinc_file.write("\n-----------------------------------------------------------------------------\n")
                break                                        # stop searching for this T-Invariant
    coinc_file.write("Coincidences found: "+str(coinc_number) + " of " + str(len(invs_set)) + " total T-Invariants")

def _expand_all_cycles(self):
    for s_cycle in self.invariants:                                # for every cycle, generate all it's combinations
        expand_singular_cycle(s_cycle, self.expanded_cicles)       # and save it to the self variable


def expand_singular_cycle(s_cycle, expanded):
    for idx, elem in enumerate(s_cycle):                            # generate a new cycle from every cycle's element
        permutation = build_permutation(s_cycle, idx)               # specifically, from the element index in the cycle
        expanded.append(permutation)                                # and save it


def build_permutation(s_cycle, idx2):
    out = [s_cycle[idx2]]                                           # begin from the current element pointed by the index
    for it in range(len(s_cycle)-1):                                # and iterate cycle's lenght - 1 times
        idx2 = (idx2+1) % len(s_cycle)                              # through the cycle's elements
        out.append(s_cycle[idx2])                                   # append it to the "new" cycle and return it
    return out
# --------end deprecated section---------