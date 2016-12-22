import pprocess
import wrapper
import time

def run(start,stop):
    # Parallel computation:
    indices = [i for i in range(start,stop)]

    nproc = len(indices)  # maximum number of simultaneous processes desired
    results = pprocess.Map(limit=nproc, reuse=1)
    parallel_function = results.manage(pprocess.MakeReusable(wrapper.run))
    tic = time.time()
    [parallel_function(ind) for ind in indices];  # Start computing things
    parallel_results = results[0:3]
    print "%f s for parallel computation." % (time.time() - tic)



if __name__ == "__main__":

    import sys, getopt
    try:
        args = sys.argv[1:]

        opt, arg = getopt.getopt(
            args, "ci",
            longopts=["start=","stop="])

    except getopt.GetoptError as err:
        print "No command line arguments"


    goodtogo1 = False
    goodtogo2 = False

    for o, a in opt:
        if o in ["--start"]:
            start = int(a)
            goodtogo1 = True
        if o in ["--stop"]:
            stop = int(a)
            goodtogo2 = True
    if goodtogo1:
        if goodtogo2:
            run(start,stop)
