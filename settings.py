__author__ = 'marti'
import Util

#---General parameters---#
MAX_PATHS = 1000
tinvScoreThreshold = 0.8
compute_paths_tree = True
align_tinvariants_reactions = True # if set to true, this program will compute an alignment similar to the proper MP-Align
# If compute_paths_tree is False, paths'll be computed, but paths_tree won't be in the output
skip_align = False
reversed_reactions_canbe_initial_nodes = True
simReact_enzimes_weight = 0.4
simReact_inputs_weight = 0.3
simReact_outputs_weight = 0.3

#---Batch alignment (loaded if no arguments specified)---#
# To load organisms from CSV, insert Util.loadFromCSV("path/to/file.csv", "delimiterChar", columnNumber, rowsToSkip, maxNumber)
# Keep in mind that organism_from must be only one organism, so if you load from CSV, a list is returned and you have to index it
organism_from = 'hsa' #Util.loadFromCSV("Extras/Eukaryotes.csv", " ", 0, 1, 1)[0]
organism_to = ['oaa','xla','tgu','smm','bdi','mja','afu']
paths = ['00020']
# paths = ['00010', '00020', '00030', '00051', '00052', '00071', '00190',
#          '00780', '00790', '00860', '00900', '00910', '00920',
#          '00230', '00240', '00250', '00260', '00270', '00280', '00290',
#          '00300', '00330', '00350', '00360', '00400', '00450', '00480',
#          '00500', '00520', '00562', '00564', '00620', '00630', '00640',
#          '00650', '00670', '00730', '00740', '00750', '00760', '00770']

#---Directories---#
results_path_part1 = 'Results/'
results_path_part2 = 'Path_Info/'
results_path = results_path_part1 + results_path_part2
summary_file_cache = results_path_part1 + 'summaryCache.pkl'
summary_file = results_path_part1 + 'summary.csv'
ec_path = results_path_part1 + 'ECs/'
alignment_path = results_path_part1 + 'Align/'
alignment_info_path = alignment_path + 'T-inv-Alignament_SUMMARY/'
libraries_path = 'Libraries/'
cache_folder = results_path_part1 + "Scorecache/"
comp_cache = "compscore.pkl"
react_cache = "reactscore.pkl"
comp_info_file = "compound2"