from Libraries.mwmatching import maxWeightMatching
from Libraries.string_align import sequence_align
from Libraries.munkres import Munkres, print_matrix_to_file, make_cost_matrix
from scoreInfo import scoreInfo
import os
import sys
import pickle
from Util import getCompoundInfo
import settings

try:
    os.makedirs(settings.cache_folder)                       # Check and create cache dir/files
except OSError:
    if not os.path.isdir(settings.cache_folder):
        print 'Error creating cache directory'
        raise
if not os.path.isfile(settings.cache_folder + settings.comp_cache):
    print "Creating " + settings.comp_cache + " ..."
    f = open(settings.cache_folder + settings.comp_cache, "wb")
    pickle.dump({}, f)
    f.close()
if not os.path.isfile(settings.cache_folder + settings.react_cache):
    print "Creating " + settings.react_cache + " ..."
    f = open(settings.cache_folder + settings.react_cache, "wb")
    pickle.dump({}, f)
    f.close()


class Alignment:
    def __init__(self, path1, path2):
        self.path1 = path1
        self.path2 = path2
        # Compare
        self.gap_pen = 0
        self.compscore = {}                         # Compound score cache

        self._load_score()

        self._rscore = {}                           # Reaction score cache
        self._path_align = None
        self._reaction_align = None
        self._reaction_align_notinvariants = None
        self.tinv_align_score = None
        self.reaction_align_score = None
        self.reaction_align_score_notinvariants = None
        self.reactions_aligned_atTinvs = []


        self._max_score_compound = 0.99             # Maximum allowed score between different compounds

    def _load_score(self):
        f = open(settings.cache_folder + settings.comp_cache, "rb")
        self.compscore = pickle.load(f)
        f.close()

        f = open(settings.cache_folder + settings.react_cache, "rb")
        self._rscore = pickle.load(f)
        f.close()

    def _update(self, f="both"):
        if f == "both" or f == "compound":
            f = open(settings.cache_folder + settings.comp_cache, "wb")
            pickle.dump(self.compscore, f)
            f.close()

        if f == "both" or f == "reaction":
            f = open(settings.cache_folder + settings.react_cache, "wb")
            pickle.dump(self._rscore, f)
            f.close()

    def path_align(self): # Function originally from MP-Align (Adapted)
        """
        Devuelve el mejor alineamiento de caminos.
        reacion_set() --> {[i,j]: score}
            Devuelve el score al comparar
                route1.reaction_set()[i]
                    con
                route2.reaction_set()[j]
        """
        if not self._path_align is None:
            return self._path_align

        set1 = self.path1.longest_paths
        set2 = self.path2.longest_paths
        graph = []
        sscore = {}

        if len(set1) == 0 or len(set2) == 0:
            self._path_align = {}
            return {}

        len1 = len(set1)

        self._max_path = max(len(set1), len(set2))
        print "Alignment work to do: " + str(len(set1) * len(set2))

        for index1 in xrange(len(set1)):
            path1 = [self.path1.react_dict[x] for x in set1[index1]]
            for index2 in xrange(len(set2)):
                path2 = [self.path2.react_dict[x] for x in set2[index2]]
                #print "COMPARING PATH: \n"+str(path1)+" \nWITH\n "+str(path2)
                tmp = self._compare_path(path1, path2)
                #print "\nRESULT: \n"+str(tmp)+"\n---------------------------"
                graph.append([index1, index2+len1, tmp[2]])     # edge's cost. I.E [0, 4, 0.9]
                sscore[(index1, index2+len1)] = tmp             # 0: path of route1, 1: path of route2, 2: score
                                                # nodes have to be unique so absolute values are used for path2 nodes
        # Paths alignment
        match = maxWeightMatching(graph, maxcardinality=True)   # list where index = node from and content = node to
        self._path_align = {}
        for i in xrange(len1):
            if match[i] != -1:                                  # if the path has a matching other
                self._path_align[(i, match[i]-len1)] = sscore[(i, match[i])]

        return self._path_align

    def _compare_path(self, path1, path2):
        if path1 == path2:
            return [path1, path2, 1.0]
        # to compare paths we align the path seq. with reference to comp_reaction function
        return sequence_align(path1, path2, self._compare_reactions, self.gap_pen)


    def reaction_align(self, align_tinvariants): # Function originally from MP-Align (Adapted)
        if align_tinvariants and not self._reaction_align is None:
            return self._reaction_align
        if not align_tinvariants and not self.reaction_align_score_notinvariants is None:
            return self._reaction_align_notinvariants

        score = self.path_align()
        repeats_dict = {}
        graph = []
        reactions1 = self.path1.d_react
        reactions2 = self.path2.d_react
        indexesRoute1 = {}
        indexesRoute2 = {}
        # print len(reactions1), len(reactions2)
        # print len(self.path1.non_accessible_reactions), len(self.path2.non_accessible_reactions)
        #print self.path1.non_accessible_reactions[0]
        # print self.reactions_aligned_atTinvs
        # print self.path1.non_accessible_reactions
        for invReact in self.reactions_aligned_atTinvs:
            self.path1.non_accessible_reactions.pop(getReactionID(invReact), None)
            self.path2.non_accessible_reactions.pop(getReactionID(invReact), None)
        # print len(self.path1.non_accessible_reactions), len(self.path2.non_accessible_reactions)

        if not align_tinvariants:
            reactions1 = self.load_reactions_NOT_at_tinv_align(1)
            reactions2 = self.load_reactions_NOT_at_tinv_align(2)
        # print len(reactions1), len(reactions2)
        N1 = len(reactions1)
        # print score
        for match in score:
            # Iterate through already computed paths alignment
            path1 = score[match][0]
            path2 = score[match][1]

            for i in xrange(len(path1)):
                if align_tinvariants or (path1[i] in reactions1 and path2[i] in reactions2):
                    reaction1 = "-"
                    if path1[i] != "-":
                        reaction1 = path1[i][0] + "-" + path1[i][1] + " " + path1[i][4]
                        indexesRoute1[reaction1] = reactions1.index(path1[i])
                    reaction2 = "-"
                    if path2[i] != "-":
                        reaction2 = path2[i][0] + "-" + path2[i][1] + " " + path2[i][4]
                        indexesRoute2[reaction2] = reactions2.index(path2[i])
                    if not reaction1 in repeats_dict:
                        repeats_dict[reaction1] = {}
                    if not reaction2 in repeats_dict[reaction1]:
                        repeats_dict[reaction1][reaction2] = [0, 0.0, 0.0]

                    repeats_dict[reaction1][reaction2][0] += 1
                    repeats_dict[reaction1][reaction2][1] += score[match][2]
                    if repeats_dict[reaction1][reaction2][2] < score[match][2]:
                        repeats_dict[reaction1][reaction2][2] = score[match][2]
        # print repeats_dict
        for react1 in repeats_dict:
            for react2 in repeats_dict[react1]:
                if react1 != "-" and react2 != "-":
                    node_from = indexesRoute1[react1]   # reaction index inside the current path
                    node_to = N1 + indexesRoute2[react2]
                    graph.append([node_from, node_to, repeats_dict[react1][react2][0]]) # no. of times aligned

        match = maxWeightMatching(graph, maxcardinality=True)
        reaction_align = {}
        align_score = 0
        if len(match) > 0:
            for i in xrange(0, N1):
                if match[i] == -1:
                    reaction_align[reactionToString(reactions1[i])] = [None, 0]
                else:
                    jj = reactionToString(reactions1[i])
                    mm = reactionToString(reactions2[match[i]-N1])
                    if jj in repeats_dict and mm in repeats_dict[jj]: #and (align_tinvariants or reactions2[match[i]-N1] in reactions2):
                        score = 1.0 * repeats_dict[jj][mm][2]
                        reaction_align[reactionToString(reactions1[i])] = [mm, score]
                        align_score += score
        # Align reactions that are NOT in any path:
        reacts1_nonaccessible = [self.path1.react_dict[i] for i in self.path1.non_accessible_reactions.keys()]
        reacts2_nonaccessible = [self.path2.react_dict[i] for i in self.path2.non_accessible_reactions.keys()]
        if len(reacts1_nonaccessible) > 0 and len(reacts2_nonaccessible) > 0:
            sim_matrix_nonpaths = self._build_reaction_similarity_matrix(reacts1_nonaccessible, reacts2_nonaccessible)
            cost_matrix = make_cost_matrix(sim_matrix_nonpaths, lambda cost: 999.0 - cost)
            m = Munkres()
            #print sim_matrix_nonpaths
            indexes = m.compute(cost_matrix)                          # compute alignament
            for row, column in indexes:
                reaction_align[reactionToString(reacts1_nonaccessible[column])] = [reactionToString(reacts2_nonaccessible[row]), sim_matrix_nonpaths[row][column]]
                align_score += sim_matrix_nonpaths[row][column]

            #nonaccessible_max_len = max(len(reacts1_nonaccessible), len(reacts2_nonaccessible))

        if align_tinvariants:
            self.reaction_align_score = align_score/max(len(reactions1), len(reactions2))
            self._reaction_align = reaction_align
        else:
            self.reaction_align_score_notinvariants = align_score/max(len(reactions1), len(reactions2))
            self._reaction_align_notinvariants = reaction_align

        return reaction_align


    def load_reactions_NOT_at_tinv_align(self, num):
        if num == 1:
            reactions_NOTat_tinvariants_list = list(self.path1.d_react)
        else:
            reactions_NOTat_tinvariants_list = list(self.path2.d_react)

        for aligned_reaction in self.reactions_aligned_atTinvs:
            try:
                reactions_NOTat_tinvariants_list.remove(aligned_reaction)
            except ValueError:      # The reaction may be already deleted, so a ValueError exception is raised
                pass
        return reactions_NOTat_tinvariants_list


    def compare_all_tinvariants(self, fpath, score_threshold, scores):
        try:
            if len(self.path1.invariants) == 0 or len(self.path2.invariants) == 0:
                raise ValueError("Detected pathway with no non-trivial T-Invariants. Skipping T-invariant alignament...")
            f_namefile = fpath + 'T-invAlign-' + self.path1.name + '-' + self.path2.name   # file to save the alignament
            file_align = open(f_namefile, 'w')
            #file_align.write('///////T-INVARIANT PAIR ALIGNAMENT//////\n\n')
            #build the final alignament matrix structure
            tinv_alignament_matrix = [[0 for i in range(len(self.path2.invariants))] for j in range(len(self.path1.invariants))]
            for i1, inv1 in enumerate(self.path1.invariants):               # compare all with all the invariants
                self._check_invs_lenghts(inv1, scores)
                for i2, inv2 in enumerate(self.path2.invariants):
                    #file_align.write('T-Invariant 1:\n'),
                    #file_align.writelines(str(inv1))
                    #file_align.write('\n\nT-Invariant 2:\n'),
                    #file_align.writelines(str(inv2))
                    sim_matrix = self._build_reaction_similarity_matrix(inv1, inv2)     # compute the sim. matrix for inv1-2
                    tinv_alignament_matrix[i1][i2] = self._compute_alignament(sim_matrix, file_align, False) # align and save to the final sim. matrix
                    inv12_score = tinv_alignament_matrix[i1][i2]
                    #file_align.write('\n-------------------------------------------------------\n\n')

                    if inv12_score > score_threshold:
                        self._add_to_score_info(inv1, inv2, inv12_score, scores)        # save it to the info record
                    self._check_invs_lenghts(inv2, scores)

            file_align.write('///////T-INVARIANT FINAL ALIGNAMENT//////\n')
            total = self._compute_alignament(tinv_alignament_matrix, file_align, True)  # compute the final alignament and score
            self.tinv_align_score = total
            file_align.close()
            return True
        except ValueError:
            print "Error: " + str(sys.exc_info()[1])
            self.tinv_align_score = 0
            return False


    def _compute_alignament(self, alignament_matrix, file_align, final):
        cost_matrix = make_cost_matrix(alignament_matrix, lambda cost: 999.0 - cost)    # invert the simmilarity matrix, because Munkres needs a cost matrix
        m = Munkres()
        indexes = m.compute(cost_matrix)                                                # compute alignament
        total = 0
        if final:
            print_matrix_to_file(alignament_matrix, file_align, '\nScores matrix:\n')
            file_align.write("\n\n")
        for row, column in indexes:
            value = alignament_matrix[row][column]          # index the alignament matrix with the alignamed computed indexes by Munkres
            total += value                                  # increment the total score
            if final:
                file_align.write('T-Invariant 1:\n'),
                file_align.writelines([str(reactionToString(x) + " , ") for x in self.path1.invariants[row]])
                file_align.write('\n\nT-Invariant 2:\n'),
                file_align.writelines([str(reactionToString(x) + " , ") for x in self.path2.invariants[column]])
                self._write_tinv_reaction_alignment(self.path1.invariants[row], self.path2.invariants[column], file_align)
                file_align.write('\n-------------------------------------------------------\n\n')
        total /= max(len(alignament_matrix), len(alignament_matrix[0])) # and divide it by the max inv. lenght, so max score is 1
        if final:
            file_align.write('\nTotal: %f' % total)
        return total

    def _write_tinv_reaction_alignment(self, inv1, inv2, file_align):
        # this function is called during the t-inv final alignment writing and repeats the alignment to write also
        # how are the two t-invs invariants in terms of their reactions and save aligned reactions
        sim_matrix = self._build_reaction_similarity_matrix(inv1, inv2)
        cost_matrix = make_cost_matrix(sim_matrix, lambda cost: 999.0 - cost)
        m = Munkres()
        indexes = m.compute(cost_matrix)
        file_align.write('\n\nAligned this way:\n')
        for row, column in indexes:
            file_align.write(str(reactionToString(inv1[column]) + " with " + reactionToString(inv2[row]) + " - " + str(sim_matrix[row][column]) + "\n"))
            self.reactions_aligned_atTinvs.append(inv1[column])
            self.reactions_aligned_atTinvs.append(inv2[row])

    def _build_reaction_similarity_matrix(self, inv1, inv2):
        similarity_matrix = [[0 for i in range(len(inv1))] for j in range(len(inv2))]

        for i2, r2 in enumerate(inv2):
            for i1, r1 in enumerate(inv1):
                similarity_matrix[i2][i1] = self._compare_reactions(r2, r1)   # score is returned by comparing reactions
                #print similarity_matrix[i2][i1]
        return similarity_matrix

    def _add_to_score_info(self, inv1, inv2, inv12_score, scores):
        inv_from = tuple([i[0] + "," + i[4] for i in inv1])
        inv_to = tuple([i[0] + "," + i[4] for i in inv2])
        key = (inv_from, inv_to, inv12_score)
        path2_name = self.path2.name[0:3]
        if key in scores:                                       # existing score info with same invs. and score
            if not path2_name in scores[key].paths:
                scores[key].paths.append(path2_name)            # add the org. to the list
                #print 'AFEGIM UN ORG. A LA PUNTUACIO'
        elif not key in scores:                                 # non existing score for these invs
            # print 'NO TROBADA. INSERTAM... amb punt: ' + str(inv12_score)
            score = scoreInfo
            score.paths = [path2_name]
            scores[key] = score                                 # insert new score info

    def _check_invs_lenghts(self, inv, scores):
        sc = scoreInfo
        if len(inv) < sc.min_inv_lenght:
            sc.min_inv_lenght = len(inv)
        elif len(inv) > sc.max_inv_lenght:
            sc.max_inv_lenght = len(inv)

    def _compare_enzyme(self, enzyme1, enzyme2): # Function originally from MP-Align (Adapted and modified)
        # enzime score is not cached because it's faster to compute it every time
        tmp1 = enzyme1.split(":")[1].split(".")  # ec:1.1.1.1 --> [1, 1, 1, 1]
        tmp2 = enzyme2.split(":")[1].split(".")

        ret = 0.0
        if tmp1 == tmp2:
            ret = 1.0
        elif tmp1[0:2] == tmp2[0:2]:
            ret = 0.75
        elif tmp1[0:1] == tmp2[0:1]:
            ret = 0.5
        elif tmp1[0] == tmp2[0]:
            ret = 0.25

        return ret

    def _compare_compound(self, compound1, compound2):
        # cpd:C00100-12 --> C00100-12 --> C00100 --> 00100
        compound1 = compound1.split(":")[1].split("-")[0][1:]
        compound2 = compound2.split(":")[1].split("-")[0][1:]

        if compound1 == compound2:
            return 1.0      # return maximum and don't save it
        else:
            # Already computed?
            penalitzacio = self._max_score_compound
            if (compound1, compound2) in self.compscore:
                return penalitzacio*self.compscore[(compound1, compound2)]
            elif (compound2, compound1) in self.compscore:
                return penalitzacio*self.compscore[(compound2, compound1)]
            # if not, compute it
            ret = 0.0
            try:
                ret = simcompHelper(compound1, compound2)
            except:
                print "SIMPCOMP error. " + compound1 + " or " + compound2 + " may not be present at " + settings.comp_info_file
                print "Attempting recuperation by retrieval from KEGG. You must have internet access for this to work..."
                if not getCompoundInfo(str("C" + compound1), settings.comp_info_file) or not getCompoundInfo(str("C" + compound2), settings.comp_info_file):
                    sys.exit(1)
                ret = simcompHelper(compound1, compound2)
            finally:
                self.compscore[(compound1, compound2)] = ret  # save
                print str(ret) + ". Resuming alignment..."
                self._update("compound")
                return ret * penalitzacio


    def _compare_set_compound(self, set1, set2): # Function originally from MP-Align (Adapted)
        if max(len(set1), len(set2)) == 0:
            return 0.0
        bi_graph = []
        len1 = len(set1)
        for cpd1 in set1:
            index1 = set1.index(cpd1)
            for cpd2 in set2:
                bi_graph.append([index1, len1 + set2.index(cpd2), self._compare_compound(cpd1, cpd2)])

        match = maxWeightMatching(bi_graph)
        score = 0.0
        i = 0
        while (i < len(match)) and (i < min(len1, len(set2))):
            if match[i] != -1:
                score += self._compare_compound(set1[i], set2[match[i]-len1])
            i += 1

        return score/max(len(set1), len(set2))

    def _compare_reactions(self, reaction2, reaction1): # Function originally from MP-Align (Adapted)
        if reaction1[0] == reaction2[0]:
            return 1.0
        if (reaction1[0], reaction2[0]) in self._rscore:
            return self._rscore[(reaction1[0], reaction2[0])]

        dcmps1 = self.path1.d_comp
        dcmps2 = self.path2.d_comp
        # Entrada
        input1 = []
        for ids in reaction1[2]:
            input1.append(dcmps1[ids][0])
        #print input1
        input2 = []
        for ids in reaction2[2]:
            input2.append(dcmps2[ids][0])
        #print input2
        iscore = self._compare_set_compound(input1, input2)

        output1 = []
        for ids in reaction1[3]:
            output1.append(dcmps1[ids][0])
        #print output1
        output2 = []
        for ids in reaction2[3]:
            output2.append(dcmps2[ids][0])
        #print output2

        oscore = self._compare_set_compound(input1, input2)

        escore = self._compare_enzyme(reaction1[4], reaction2[4])

        reacts_score = escore * settings.simReact_enzimes_weight + iscore * settings.simReact_inputs_weight + oscore * settings.simReact_outputs_weight
        self._rscore[(reaction1[0], reaction2[0])] = reacts_score
        self._rscore[(reaction2[0], reaction1[0])] = reacts_score
        return reacts_score


def reactionToString(rTuple):
    return rTuple[0] + "-" + rTuple[1] + " " + rTuple[4]

def getReactionID(rTuple):
    if 'rev' in rTuple[0]:
        return rTuple[1] + "rev"
    else:
        return rTuple[1]

def simcompHelper(compound1, compound2):
    print "Retrieving cpd:C" + str(compound1) + " vs cpd:C"+ str(compound2)+" similarity from simcomp...",
    return float(os.popen("./simcomp " + settings.comp_info_file + " -e " + compound1 + \
                                 " -f " + compound2 + \
                                 " | awk {'print $6'}").read())
