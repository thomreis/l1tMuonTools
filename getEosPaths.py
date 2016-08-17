import argparse
import subprocess
import commands

def parse_options_and_init_log():
    """
    Adds often used options to the OptionParser...
    """
    parser = argparse.ArgumentParser(description="EOS indexer macro", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-p", "--path", dest="path", default="", type=str, help="An EOS path where files should be searched.")
    parser.add_argument("-f", "--file-pattern", dest="filePattern", default="L1Ntuple", type=str, help="File name pattern to look for.")
    parser.add_argument("-o", "--out", dest="outfile", default="files.txt", type=str, help="Output file.")
    parser.add_argument("-a", "--append", dest="append", default=False, action="store_true", help="Append to file instead of overwrite.")

    opts, unknown = parser.parse_known_args()

    return opts


def main():
    opts = parse_options_and_init_log()

    eosPath = opts.path
    fpattern = opts.filePattern

    dateTime = long(0)
    latestGridSubTimeStr = ''

    #print 'eos find -d {path}/'.format(path=eosPath)
    proc = subprocess.Popen(['eos', 'find', '-d', '{path}/'.format(path=eosPath)], stdout=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    print stderr
    print stdout
    allPaths = stdout.split('\n')
    writeMode = 'w'
    if opts.append:
        writeMode = 'a'
    with open(opts.outfile, writeMode) as outfile:
        # first find the path with the latest grid submission time in case there are several
        for pathStr in allPaths:
            if len(pathStr) > 0:
                gridSubTimeStrPos = pathStr.find('/16')
                if gridSubTimeStrPos == -1:
                    continue
                gridSubTimeStr = pathStr[gridSubTimeStrPos+1:gridSubTimeStrPos+14]
                dateTimeStrs = gridSubTimeStr.split('_')
                dateTimeStr = dateTimeStrs[0]+dateTimeStrs[1]
                thisDateTime = long(dateTimeStr)
                if thisDateTime <= dateTime:
                    continue
                else:
                    dateTime = thisDateTime
                    latestGridSubTimeStr = gridSubTimeStr

        print 'Using files from the grid submission '+latestGridSubTimeStr

        # find the files from the latest grid submission                
        for pathStr in allPaths:
            if pathStr.find(latestGridSubTimeStr) == -1:
                continue
            if len(pathStr) > 0:
                spacePos = pathStr.find(' ')
                path = pathStr[0:spacePos]
                if path.rfind('/failed/') != -1:
                    continue
                proc = subprocess.Popen(['eos', 'ls', '-l', '{path}/'.format(path=path)], stdout=subprocess.PIPE)
                stdout, stderr = proc.communicate()
                lines = stdout.split('\n')
                for line in lines:
                    if len(line) > 0:
                        if line[0] != 'd':
                            spacePos = line.rfind(' ')
                            filestr = line[spacePos+1:]
                            if filestr.find(fpattern) != -1:
                                outfile.write('root://eoscms.cern.ch/'+path+filestr+'\n')

        outfile.close()

if __name__ == "__main__":
    main()
