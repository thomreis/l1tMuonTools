#!/usr/bin/env python
import argparse
import json
from operator import itemgetter
from itertools import groupby

"""
Generates the AND or OR of good lumi section json files
"""

def parse_options():
    desc = "Interface for JSON merger."
    parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--in1", dest="inPath2", type=str, help="Input file name 1.", default=None)
    parser.add_argument("--in2", dest="inPath1", type=str, help="Input file name 2.", default=None)
    parser.add_argument("--out", dest="outPath", type=str, help="Output file name.", default='merged.json')
    parser.add_argument("--mode", dest="mode", type=str, help="AND or OR.", default='AND')
    
    return parser.parse_args()

def and_ls_jsons(json1, json2):
    merged_json = {}
    keys2 = json2.keys()
    for key1 in json1.keys():
        # if both jsons contain this run number
        if key1 in keys2:
            # expand the ranges
            valuess1 = json1[key1]
            valuess2 = json2[key1]
            lss1 = []
            lss2 = []
            for values1 in valuess1:
                lss1 += range(values1[0], values1[1]+1)
            for values2 in valuess2:
                lss2 += range(values2[0], values2[1]+1)

            # make the intersection of the two lists. Works best with sets
            if len(lss1) > len(lss2):
                common_lss = sorted(list(set(lss2).intersection(lss1)))
            else:
                common_lss = sorted(list(set(lss1).intersection(lss2)))

            common_lss_grouped = []
            # from itertools documentation
            for k, g in groupby(enumerate(common_lss), lambda (i,x):i-x):
                lsrange = map(itemgetter(1), g)
                if len(lsrange) > 1:
                    common_lss_grouped.append([lsrange[0], lsrange[-1]])
                else:
                    common_lss_grouped.append([lsrange[0], lsrange[0]])

            # put everything in the new dictionary if at least one range is present
            if len(common_lss_grouped) > 0:
                merged_json[key1] = common_lss_grouped
    return merged_json

def main():
    options = parse_options()

    if options.inPath1 and options.inPath2:
        with open(options.inPath1) as json_file1:
            with open(options.inPath2) as json_file2:
                if options.mode == 'AND':
                    merged_dict = and_ls_jsons(json.load(json_file1), json.load(json_file2))
                elif options.mode == 'OR':
                    print 'OR is not yet implemented.'
                    exit()
                else:
                    print "Unknown mode '"+options.mode+"'. Choose 'AND' or 'OR'."
                    exit()

                with open(options.outPath, 'w') as outfile:
                    json.dump(merged_dict, outfile)

if __name__ == "__main__":
    main()

