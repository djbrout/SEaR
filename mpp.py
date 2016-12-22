import multiprocessing
import wrapper

def run(start,stop):
    # Parallel computation:

    jobs = []
    for i in range(start,stop):
        p = multiprocessing.Process(target=wrapper.run, args=(i,))
        jobs.append(p)
        p.start()

    for j in jobs:
        j.join()
        print '%s.exitcode = %s' % (j.name, j.exitcode)

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
