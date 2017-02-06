from sys import exit
import argparse
import logging
import os

from L1Analysis import L1Ana, L1Ntuple

def parse_options_and_init_log(loglevel=logging.INFO):
    """
    Adds often used options to the OptionParser...
    """
    parser = argparse.ArgumentParser(description="L1 Analysis Framework macro", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-f", "--fname", dest="fname", default="", type=str, help="A root file containing L1Ntuples.")
    parser.add_argument("-l", "--flist", dest="flist", default="", type=str, help="A txt file containing list of L1Ntuple files, one file per line.")
    parser.add_argument("-q" ,"--queue", dest="queue", default="1nh", help="Queue to submit the jobs to.")
    parser.add_argument("-j", "--njobs", dest="njobs", default=1, type=int, help="Number of jobs to submit")
    parser.add_argument("-n", "--nevents", dest="nevents", default=-1, type=int, help="Number of events to run, -1 for all [default: %default]")
    parser.add_argument("-w", "--workdir", dest="workdir", default='job', type=str, help="Work directory to create scripts and store output in [default: %default]")
    parser.add_argument("-o", "--outname", dest="outname", default="output", type=str, help="File name for output files. .root added automatically [default: %default]")
    parser.add_argument("-s", "--scriptname", dest="scriptname", default="ntuple.py", type=str, help="Script to run [default: %default]")
    parser.add_argument("-p", "--subparser", dest="subparser", default="ntuple", type=str, help="Subparser for script [default: %default]")
    parser.add_argument("--split_by_file", dest="split_by_file", action="store_true", help="File based splitting instead of event number based splitting")
    parser.add_argument("--json", dest="json", type=str, default=None, help="json file with good lumi sections")
    parser.add_argument("--runs", dest="runs", type=str, default=None, help="list of runs to analyse")
    parser.add_argument("--pos-side", dest="pos_side", default=False, action="store_true", help="Positive detector side only.")
    parser.add_argument("--neg-side", dest="neg_side", default=False, action="store_true", help="Negative detector side only.")
    parser.add_argument("--pos-charge", dest="pos_charge", default=False, action="store_true", help="Positive probe muon charge only.")
    parser.add_argument("--neg-charge", dest="neg_charge", default=False, action="store_true", help="Negative probe muon charge only.")
    parser.add_argument("--use-inv-mass-cut", dest="invmass", default=False, action="store_true", help="Use an invariant mass range for the tag and probe pair.")
    parser.add_argument("--use-extra-coord", dest="extraCoord", default=False, action="store_true", help="Use L1 extrapolated eta and phi coordinates.")
    parser.add_argument("--eta-restricted", dest="etarestricted", type=float, default=None, help="Upper eta value for isolation.")
    parser.add_argument("--emul", dest="emul", default=False, action="store_true", help="Make emulator histograms.")
    parser.add_argument("--legacy", dest="legacy", default=False, action="store_true", help="Use legacy muons translated to upgrade format.")
    parser.add_argument("--prefix", dest="prefix", type=str, default=None, help="Prefix for histogram names")
    parser.add_argument("--tftype", dest="tftype", type=str, default=None, help="Fill L1 muons from one TF.")
    parser.add_argument("--eta-bits", dest="etabits", type=int, default=None, help="Number of eta input bits for extrapolation LUT.")
    parser.add_argument("--iso-method", dest="isomethod", type=str, default=None, help="Isolation method. ['abs', 'rel', 'inner', 'outovertot', 'inner2x2', 'outovertot2x2']")
    parser.add_argument("--nvtx-min", dest="nvtxmin", type=int, default=None, help="Minimum number of vertices.")

    opts, unknown = parser.parse_known_args()
    if opts.fname == "" and opts.flist == "":
        from L1Analysis import L1Ana
        L1Ana.init_logging("L1Analysis", loglevel)
        L1Ana.log.fatal("No inputs specified!")
        parser.print_help()
        exit(0)
    else:
        from L1Analysis import L1Ana
        L1Ana.init_logging("L1Analysis", loglevel)

    return opts


def main():
    L1Ana.init_l1_analysis()
    print ""

    opts = parse_options_and_init_log()

    if opts.split_by_file:
        n_per_job = -1
        with open(opts.flist) as flistfile:
            lines = flistfile.readlines()
            nLines = len(lines)
            files_per_job = nLines / opts.njobs
    else:
        ntuple = L1Ntuple(opts.nevents)
        if opts.flist:
            ntuple.open_with_file_list(opts.flist)
        if opts.fname:
            ntuple.open_with_file(opts.fname)

        n_per_job = ntuple.nevents / opts.njobs

    # first make sure the directories are created
    if not os.path.exists(opts.workdir+"/out"):
        os.makedirs(opts.workdir+"/out")
    if not os.path.exists(opts.workdir+"/scripts"):
        os.makedirs(opts.workdir+"/scripts")
    if opts.flist and opts.split_by_file and not os.path.exists(opts.workdir+"/filelists"):
        os.makedirs(opts.workdir+"/filelists")

    start_up = "cd {cmssw_dir}/src\n".format(cmssw_dir=os.environ["CMSSW_BASE"])
    start_up += "eval `scram runtime -sh`\n"
    start_up += "cd {pwd}\n".format(pwd=os.getcwd())

    submission_string = ""
    hadd_string = "hadd " + os.path.abspath(opts.outname+"/out/{oname}_comb.root".format(oname=opts.outname))+" "
    # job_dir = os.path.abspath(opts.outname+"/scripts/")

    job_dir = os.path.abspath(opts.workdir+"/scripts/")
    out_dir = os.path.abspath(opts.workdir+"/out/")

    for i in range(opts.njobs):
        with open(opts.workdir+"/scripts/job_{i}.sh".format(i=i), "w") as job_script:
            job_script.write(start_up)
            outfile = opts.workdir+"/out/{name}_{n}.root".format(name=opts.outname, n=i)
            hadd_string += outfile + " "
            if opts.flist:
                if opts.split_by_file:
                    flistpath = opts.workdir+"/filelists/flist_{i}.txt".format(i=i)
                    with open(flistpath, "w") as flistfile_out:
                        for j in range(files_per_job):
                            flistfile_out.write(lines[i*files_per_job + j])
                        flistfile_out.close()
                    py_string = "python {script} -l {flist} {subparser} -o {out}"
                    py_string = py_string.format(script=opts.scriptname, flist=flistpath, subparser=opts.subparser, out=outfile)
                else:
                    py_string = "python {script} -l {flist} -n {n} -s {start} {subparser} -o {out}"
                    py_string = py_string.format(script=opts.scriptname, flist=opts.flist, n=n_per_job, start=i*n_per_job, subparser=opts.subparser, out=outfile)
            if opts.fname:
                py_string = "python {script} -f {fname} -n {n} -s {start} {subparser} -o {out}"
                py_string = py_string.format(script=opts.scriptname, fname=opts.fname, n=n_per_job, start=i*n_per_job, subparser=opts.subparser, out=outfile)
            if opts.json:
                py_string += " --json {json}".format(json=opts.json)
            if opts.runs:
                py_string += " --runs {runs}".format(runs=opts.runs)
            if opts.pos_side:
                py_string += " --pos-side"
            if opts.neg_side:
                py_string += " --neg-side"
            if opts.pos_charge:
                py_string += " --pos-charge"
            if opts.neg_charge:
                py_string += " --neg-charge"
            if opts.invmass:
                py_string += " --use-inv-mass-cut"
            if opts.extraCoord:
                py_string += " --use-extra-coord"
            if opts.etarestricted:
                py_string += " --eta-restricted {eta}".format(eta=opts.etarestricted)
            if opts.emul:
                py_string += " --emul"
            if opts.legacy:
                py_string += " --legacy"
            if opts.prefix:
                py_string += " --prefix {pref}".format(pref=opts.prefix)
            if opts.tftype:
                py_string += " --tftype {tf}".format(tf=opts.tftype)
            if opts.etabits:
                py_string += " --eta-bits {b}".format(b=opts.etabits)
            if opts.isomethod:
                py_string += " --iso-method {i}".format(i=opts.isomethod)
            if opts.nvtxmin:
                py_string += " --nvtx-min {i}".format(i=opts.nvtxmin)
            py_string += "\n"
            job_script.write(py_string)
            sub_string = "bsub -q {queue} -cwd {cwd} -J job_{i} {dir}/job_{i}.sh\n".format(queue=opts.queue, cwd=out_dir, dir=job_dir, i=i)
            submission_string += sub_string
        os.system('chmod 744 {dir}/job_{i}.sh'.format(dir=job_dir, i=i))

    with open(opts.workdir+"/submit.sh", "w") as submitfile:
        submitfile.write(submission_string)
    os.system('chmod 744 {dir}/submit.sh'.format(dir=opts.workdir))
    with open(opts.workdir+"/combine.sh", "w") as combfile:
        combfile.write(hadd_string)
    os.system('chmod 744 {dir}/combine.sh'.format(dir=opts.workdir))

    print "Will process", n_per_job, "events per job"
    print "execute", opts.workdir+"/submit.sh", "to submit"
    print "after completion run", opts.workdir+"/combine.sh", "to combine the ntuples."

if __name__ == "__main__":
    main()
