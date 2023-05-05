#!/usr/bin/env python
from path import Rpath
from path import find_invariants_paths
from Alignmanent import Alignment, reactionToString
from scoreInfo import scoreInfo
from graphTools import *
from collections import OrderedDict
import subprocess
import urllib2
import os, shutil
import cPickle as pickle
import sys, time
from datetime import datetime as dt
from PathSummary import percentilePreSorted
import csv
import Util, settings


def _check_create_folder(folder):
    try:
        os.makedirs(folder)
    except OSError:
        if not os.path.isdir(folder):
            print 'error creating ' + folder + ' directory'
            raise

def _check_create_pickle_file(file):
    if not os.path.isfile(file):
        print "creating cache summary..."
        f = open(file, "wb")
        pickle.dump({}, f)
        f.close()


if len(sys.argv) >= 4:
    org_from = sys.argv[1]
    org_to = [sys.argv[2]]
    paths = [sys.argv[3]]
elif len(sys.argv) == 3:
    org_from = sys.argv[1]
    org_to = []
    paths = [sys.argv[2]]
else:
    org_from = settings.organism_from
    org_to = settings.organism_to
    paths = settings.paths

_check_create_folder(settings.results_path_part1)
_check_create_folder(settings.results_path)
_check_create_folder(settings.ec_path)
_check_create_folder(settings.alignment_path)
_check_create_folder(settings.alignment_info_path)
_check_create_pickle_file(settings.summary_file_cache)
tinv_repeats = {}

def main():
    print('////////T&MP-Align: Comparison, alignment and analysis of metabolic pathways////////')
    generalSummaryFile = open(settings.summary_file_cache, 'rb')
    summaryDict = pickle.load(generalSummaryFile)
    generalSummaryFile.close()
    for path in paths:
        path_from = org_from + path
        print "First pathway selected: " + path_from
        path_a = _process_path(org_from, path)
        # if hasattr(path_a, 'invariants'):
        #     add_to_repeats_dict(path_a.invariants, path)
        if path_a is not False:
            graphs_process(path_a, path, settings.compute_paths_tree)
            process_summary(path_a, summaryDict)
            scores = {}
            #align_success = False
            for t_org in org_to:
                try:
                    compute_align(path, path_from, t_org, path_a, summaryDict, scores)
                except:
                    timestr = time.strftime("%Y%m%d")
                    print "WARNING: an exception was caught and saved into a log file, skipping to the next alignment..."
                    with(open("Exceptions-" + timestr, "a+")) as exc_file:
                        exc_file.write("Exception while aligning " + path_from + " with " + t_org + path + " :")
                        exc_file.write(str(sys.exc_info()))
                        exc_file.write("\n\n")

            _process_tinv_alignment_info(scores, path)
    print "Processing summary..."
    generalSummaryFile = open(settings.summary_file_cache, 'wb')
    pickle.dump(summaryDict, generalSummaryFile)
    write_summary_to_file(summaryDict)
    # print_tinv_dict(2)
    print "Finished!"

def add_to_repeats_dict(tinvs, pathCode):
    for tinv in tinvs:
        tinv_tuple = tuple([x[0] for x in tinv])
        if tinv_tuple in tinv_repeats and not pathCode in tinv_repeats[tinv_tuple]:
            tinv_repeats[tinv_tuple].append(pathCode)
        else:
            tinv_repeats[tinv_tuple] = [pathCode]

def print_tinv_dict(threshold):
    with(open('tinvConservation', 'w')) as f:
        for key in tinv_repeats:
            if len(tinv_repeats[key]) >= threshold:
                f.write(str(key) + '\n' + str(tinv_repeats[key]))
                f.write('\n----------------\n')
        f.write('T-invariants analysed: ' + str(len(tinv_repeats)) + '\n')
        f.write('Paths analysed: ' + str(len(settings.paths)))

def compute_align(path, path_from, t_org, path_a, summaryDict, scores):
    path_to = t_org + path
    print "Second pathway selected: " + path_to
    path_b = _process_path(t_org, path)
    graphs_process(path_b, path, settings.compute_paths_tree)
    process_summary(path_b, summaryDict)
    if path_b is not False and not settings.skip_align:
        alignment = Alignment(path_a, path_b)
        if not settings.align_tinvariants_reactions:
            print "Aligning T-Invariants. Please wait..."
            align_success = alignment.compare_all_tinvariants(settings.alignment_path, settings.tinvScoreThreshold, scores)
        print "Aligning paths. Please wait..."
        path_align = alignment.path_align()
        process_paths_alignment(path_a.name, path_b.name, path_align)
        #reaction_align = alignment.reaction_align(align_tinvariants=True)
        #process_reactions_alignment(alignment, reaction_align, align_tinvariants=True)
        reaction_align = alignment.reaction_align(align_tinvariants=settings.align_tinvariants_reactions)
        process_reactions_alignment(alignment, reaction_align, settings.align_tinvariants_reactions)

        if not settings.align_tinvariants_reactions:
            print path_from + " vs " + path_to + " T-Invariant alignment score: %f" % alignment.tinv_align_score
            print path_from + " vs " + path_to + " reaction alignment score (Excluding T-Invariants reactions): " + str(alignment.reaction_align_score_notinvariants)
        else:
            print path_from + " vs " + path_to + " reaction alignment score (Like MP-Align): " + str(alignment.reaction_align_score)


def process_paths_alignment(path1, path2, path_align):
    path_align_filename = settings.alignment_path + "align_path-" + path1 + "-" + path2 + ".txt"
    out = open(path_align_filename, "w+")

    for b in path_align:
        r1 = path_align[b][0]
        r2 = path_align[b][1]
        out.write("------ " + str(path_align[b][2]) + " -----\n")
        for i in xrange(len(r1)):
            out.write(str(r1[i]) + "\t" + str(r2[i]) + "\n")
        out.write("\n\n")

    out.close()

def process_reactions_alignment(alignment, reaction_align, align_tinvariants):
    if align_tinvariants:
        filename_prefix = "align_reaction(MPAlign)-"
        reactions1 = alignment.path1.d_react
        score = str(alignment.reaction_align_score)
    else:
        filename_prefix = "align_reaction(T&MPAlign)-"
        reactions1 = alignment.load_reactions_NOT_at_tinv_align(1)
        score = str(alignment.reaction_align_score_notinvariants)
    reaction_align_filename = settings.alignment_path + filename_prefix + alignment.path1.name + "-" + \
                              alignment.path2.name + ".txt"
    out = open(reaction_align_filename, "w+")

    for b in reaction_align:
        if reaction_align[b][0] is not None:
            out.write(b + "\t---"  + str(reaction_align[b][1]) + "---\t" + reaction_align[b][0] + "\n")
        else:
            out.write(b + "\t---\t0\n")

    for c in reactions1:
        rname = reactionToString(c)
        if rname not in reaction_align:
            out.write(rname + "\t---\t0\n")
    out.write("\nReaction alignment final score: " + str(score))
    out.close()

def process_summary(path_x, summaryDict):
    orderedLongPaths = sorted(path_x.summary.pathsLong)
    if len(orderedLongPaths) > 0:
        path_x.summary.maxPathLong = orderedLongPaths[len(orderedLongPaths) - 1]
        path_x.summary.minPathLong = orderedLongPaths[0]
        path_x.summary.avrgPathLong = sum(orderedLongPaths) / float(len(orderedLongPaths))
        path_x.summary.percentilPath20 = percentilePreSorted(orderedLongPaths, 20)
        path_x.summary.percentilPath50 = percentilePreSorted(orderedLongPaths, 50)
        path_x.summary.percentilPath70 = percentilePreSorted(orderedLongPaths, 70)
        path_x.summary.percentilPath90 = percentilePreSorted(orderedLongPaths, 90)
    summaryDict[path_x.name] = path_x.summary


def write_summary_to_file(summaryDict):
    ordered_summary_dict = OrderedDict(sorted(summaryDict.items(), key=lambda t: t[1].pathsComputeTime))

    header = ["Org/Pathway", "Computing paths", "Path's No.", "T-Inv. No.", "Max Path Lenght", "Min Path Lenght",
                     "Avrg Path Lenght", "Path 20% perc.", "Path 50% perc.", "Path 70% perc.", "Path 90% perc."]

    with open(settings.summary_file, 'wb') as csvFile:
        csvWriter = csv.writer(csvFile, quotechar='|', quoting=csv.QUOTE_MINIMAL)
        csvWriter.writerow(header)

        for key in reversed(ordered_summary_dict):
            row = [key]
            row += ordered_summary_dict[key].to_list()
            csvWriter.writerow(row)


def graphs_process(path_x, pathway, compute_paths_tree = False):
    t1 = dt.now()
    hypergraph_to_reaction_graph(path_x)
    graph_to_dot(path_x.reaction_graph, path_x.results_path, 'reaction_graph')
    if compute_paths_tree:
        path_x.paths_tree = reaction_graph_to_path_tree(path_x)
        graph_to_dot(path_x.paths_tree, path_x.results_path, 'paths_tree')
    else:
        path_x.longest_paths = reaction_graph_to_longest_paths(path_x)

    t2 = dt.now()
    seconds = (t2 - t1).total_seconds()
    if path_x.summary.paths_limit_raised:
        seconds = sys.maxint
        path_x.summary.pathsNumber = sys.maxint
    path_x.summary.pathsComputeTime = seconds


def _process_tinv_alignment_info(scores, path):
    if len(scores) > 0:
        sc = scoreInfo
        namefile = 'Align-Summary_' + org_from + path
        f = open(settings.alignment_info_path + namefile, 'w')
        ordered_scores = OrderedDict(sorted(scores.items(), key=lambda t: t[0][2])) # sort scores
        f.write('///////T-INVARIANT ALIGNAMENT SUMMARY FOR ' + org_from + path + ' //////\n\n')
        f.write('Maximum T-invariant lenght: ' + str(sc.max_inv_lenght) +
                '  Minimum T-invariant lenght: ' + str(sc.min_inv_lenght) + '\n')
        f.write('Maximum alignament score: ' + str(ordered_scores.items()[len(ordered_scores) - 1][0][2]) +
                '  Minimum alignament score: ' + str(ordered_scores.items()[0][0][2]))
        f.write('\n\nSCORE       -       ORGS WITH THAT SCORE\n')
        for key in reversed(ordered_scores):
            f.write(str(round(key[2], 2)) + " - "),             # write score
            f.writelines(x + ", " for x in scores[key].paths)   # write organisms with that score and invariant
            f.write('\n')
            f.writelines(x + ", " for x in key[0])              # invariant from org_from
            if not key[2] == 1:
                f.write('\n')
                f.write('VS\n')                                 # if the score is < 1, then it's a different invariant
                f.writelines(x + ", " for x in key[1])          # so we append it
            f.write('\n\n')
        f.close()


def _process_path(org, pathway):
    path_folder = settings.results_path + org + pathway + '/'
    _check_create_folder(path_folder)

    kgml_file = _get_kgmlfile(path_folder, str(org+pathway))        # Download required files if they don't exist
    ec_file = _get_kgmlfile(settings.ec_path, str('ec' + pathway))
    rpath = Rpath(org+pathway, kgml_file, ec_file, path_folder)
    if not rpath.parse():                                           # Parse KGML pathway file
        return False
    rpath.incidence_matrix()
    rpath.build_inc_matrix_file()
    rpath.invariants_info(settings.libraries_path)                           # Find non-trivial T-Invariants
    if not settings.align_tinvariants_reactions:
        print "No. of T-invariants: " +str(len(rpath.invariants))
    else:
        print "Skipped T-invariants retrieval (align_tinvariants_reactions is set to True)"
    return rpath


def _get_kgmlfile(kgmldir, metab_path):
    try:
        path_kgmlfile = str(kgmldir+metab_path+".kgml")
        if not os.path.isfile(str(kgmldir+metab_path+".kgml")):
            print 'Downloading KGML file for '+metab_path+'...'
            url = urllib2.urlopen(str("http://www.kegg.jp/kegg-bin/download?entry="+metab_path+"&format=kgml"))
            kgmlfile = open(path_kgmlfile, 'w')
            content = url.read()
            if len(content) < 2:
                raise Exception
            kgmlfile.write(content)
            kgmlfile.close()
        return path_kgmlfile
    except:
        errorMessage = "Error downloading files. You may have entered an invalid pathway name or internet conection is lost"
        print errorMessage
        if os.path.exists(str(kgmldir)):
            shutil.rmtree(str(kgmldir))
        raise ValueError(errorMessage)

#----Deprecated section----
def _spider_script():
    valids = []
    organisms = ["afu", "ape", "hin", "lin", "mja", "mmu", "nme", "rno"]
    paths_file = open("pathways_search.txt", "r")
    lines_list = paths_file.readlines()
    pathways_d = [line.split() for line in lines_list]
    for indPath in pathways_d:
        ok = True
        if len(indPath) > 8:
            for org in organisms:
                path_folder = settings.results_path + org + indPath[8] + '/'
                if _get_kgmlfile(path_folder, str(org + indPath[8])) == False:
                    ok = False
                    break
            if ok:
                valids.append(indPath[8])
                print valids
    valids_file = open("valids.txt", "w")
    valids_file.writelines(valids)

    paths_file.close()


def _search_tinvariants_at_paths(a_name, path_a):
    # (DEPRECATED FUNCTION) Finds paths where T-Invariants are
    inv_paths_coinc_file = open(str("inv-paths-"+a_name), "w")

    print "MPAlign paths: "
    MPA_path_sets = _get_path_sets(a_name)
    _print_path_sets(MPA_path_sets)

    print "Normalized T-invariats to MPAlign format: "
    normalized_set_list = path_a.invariants_to_sets()
    _print_path_sets(normalized_set_list)

    find_invariants_paths(normalized_set_list, MPA_path_sets, inv_paths_coinc_file)
    inv_paths_coinc_file.close()


def _print_path_sets(path_sets):
    for path_set in path_sets:
        for reaction in path_set:
            print reaction + ",",
        print "\n"


def _get_path_sets(a_name):
    current_dir = os.getcwd()
    os.chdir("../MPAlign")
    pnml_name = str(a_name + ".pnml")
    subprocess.call(['python', '../MPAlign/scores.py', a_name[3:8], pnml_name, pnml_name, 'onlySets'])
    os.chdir = current_dir
    sets_file = open(str("../MPAlign/conjunts-" + a_name + ".p"), "r")
    paths = pickle.load(sets_file)
    sets = []
    for path in paths:
        sets.append(set(path))
    return sets


if __name__ == '__main__':
    main()
    #_spider_script()

