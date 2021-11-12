from gd_ilp import run_optimizer
from gd_io import read_interval_file
import timeit
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("input_file", help="Default input is a file in RME format. The first line consist of the number "
                                       "of gene trees n. The next n lines represent gene trees. Then, the following "
                                       "line consist of the species tree.")
parser.add_argument('--intervallimit', dest="interval_Limit", type=int, default=0, nargs="?",
                    help="Set the maximal length of intervals. Default 0, and no intervals are shortened. If equals"
                         " x>0, then all intervals with length y>x are shortened to the size x. For an interval "
                         "<lca-mapping,b> nodes are removed from 'b' side. In other words we limit the distance between"
                         " the duplication mapping and lca-mapping not to exceed x. The value x = k+1, for the "
                         "definition of k-GME in the paper.  ")
parser.add_argument('--max_converted_spec', dest="max_converted_spec", type=int, default=-1, nargs="?",
                    help="Set the maximal number of speciations that can be converted into duplications.")
parser.add_argument('--time_limit', dest="time_limit", type=int, default=-1, nargs="?",
                    help="Set the maximal computation time in seconds.")


# Description
# def run_optimizer(stree,n,gf,debug=0,interval_limit=0,max_converted_spec=-1)
# stree - species tree (string)
# gtreesstr - list of gene trees (string)
# optional :
# debug=0 - if x>0 then shows additional info like ME score for all nodes 
# interval_limit=0  - if x>0 then shortens all intervals to the size==x
# max_converted_spec=-1 - if x>-1 then set the number of speciations to be converted into duplications

args = parser.parse_args()


print("Input file: " + args.input_file)
print("Interval Limit: " + str(args.interval_Limit))
print("Max number of converted speciations: " + str(args.max_converted_spec))

inputfile = args.input_file
# debug = 0  # to show ME score for particular nodes
interval_limit = args.interval_Limit
max_converted_spec = args.max_converted_spec
outfile = "results/out_" + inputfile[inputfile.rfind("/")+1:inputfile.find(".")] + "_" \
          + str(interval_limit) + "_" + str(max_converted_spec)  # if no extension like .txt, omit last letter

gtrees, st = read_interval_file(inputfile)
print("Output file: " + outfile + ".csv")

time_limit = args.time_limit
print("Time limit: ", time_limit)

print("\n\nRunning ILP model")
# NOTE: input trees are strings
gtreesstr = [str(x) for x in gtrees]
start = timeit.default_timer()
run_optimizer(str(st), gtreesstr, interval_limit, max_converted_spec, time_limit, outfile)
stop = timeit.default_timer()
print('Time: ', stop - start) 
