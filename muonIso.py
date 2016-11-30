#!/usr/bin/env python
from ToolBox import parse_options_and_init_log
# have to do this first or ROOT masks the -h messages
opts, parser = parse_options_and_init_log()

from L1Analysis import L1Ana, L1Ntuple
from analysis_tools.plotting import HistManager, HistManager2d
from analysis_tools.selections import MuonSelections, Matcher
import exceptions
import json
import ROOT as root

def parse_options_upgradeMuonHistos(parser):
    """
    Adds often used options to the OptionParser...
    """
    parsers = parser.add_subparsers()
    sub_parser = parsers.add_parser("muonIso")
    sub_parser.add_argument("-o", "--outname", dest="outname", default="./muon_iso_histos.root", type=str, help="A root file name where to save the histograms.")
    sub_parser.add_argument("-j", "--json", dest="json", type=str, default=None, help="A json file with good lumi sections per run.")
    sub_parser.add_argument("-r", "--runs", dest="runs", type=str, default=None, help="A string of runs to check.")

    opts, unknown = parser.parse_known_args()
    return opts

def get_tftype(tf_muon_index):
    if tf_muon_index > 35 and tf_muon_index < 72:
        return 0 # BMTF
    elif tf_muon_index > 17 and tf_muon_index < 90:
        return 1 # OMTF
    else:
        return 2 # EMTF

def book_histograms():

    vars_bins = [['iet', 256, 0, 256], ['ieta', 83, -41, 42], ['iphi', 72, 0, 73], ['qual', 16, 0, 15], ['n', 2017, 0, 4034]]
    x_title_vars = {'iet':'iE_{T}', 'ieta':'i#eta', 'iphi':'i#phi', 'qual':'qual', 'n':'# towers'}
    x_title_units = {'iet':None, 'ieta':None, 'iphi':None, 'qual':None, 'n':None}

    x_vars_bins_2d = [['ieta', 83, -41, 42]]
    y_vars_bins_2d = [['iphi', 72, 0, 73]]

    x_title_vars_2d = {'ieta':'i#eta'}
    y_title_vars_2d = {'iphi':'i#phi'}

    x_title_units_2d = {'ieta':None}
    y_title_units_2d = {'iphi':None}

    varnames = []
    binnings = {}
    varnames2d = []
    binnings2d = {}

    for var_bin in vars_bins:
        varnames.append('l1_caloTower.{var}'.format(var=var_bin[0]))
        binnings['l1_caloTower.{var}'.format(var=var_bin[0])] = var_bin[1:]+[x_title_vars[var_bin[0]], x_title_units[var_bin[0]]]

    for var_bin_2d_x, var_bin_2d_y in zip(x_vars_bins_2d, y_vars_bins_2d):
        varnames2d.append('2d_caloTower.{xvar}_{yvar}'.format(xvar=var_bin_2d_x[0], yvar=var_bin_2d_y[0]))
        binnings2d['2d_caloTower.{xvar}_{yvar}'.format(xvar=var_bin_2d_x[0], yvar=var_bin_2d_y[0])] = [var_bin_2d_x[1:]+[x_title_vars_2d[var_bin_2d_x[0]], x_title_units_2d[var_bin_2d_x[0]]], var_bin_2d_y[1:]+[y_title_vars_2d[var_bin_2d_y[0]], y_title_units_2d[var_bin_2d_y[0]]]]

    return HistManager(list(set(varnames)), binnings), HistManager2d(list(set(varnames2d)), binnings2d)

def analyse(evt, hm, hm2d):
    
    l1MuColl = evt.upgrade
    l1CaloTwrColl = evt.caloTowers

    #bx_min = 0
    #bx_max = 0

    histoprefix = 'l1_caloTower'
    histoprefix2d = '2d_caloTower'

    nCaloTwr = l1CaloTwrColl.nTower
    hm.fill(histoprefix+'.n', nCaloTwr)

    for idx in range(nCaloTwr):
        hm.fill(histoprefix+'.iet', l1CaloTwrColl.iet[idx])
        hm.fill(histoprefix+'.ieta', l1CaloTwrColl.ieta[idx])
        hm.fill(histoprefix+'.iphi', l1CaloTwrColl.iphi[idx])

        hm2d.fill(histoprefix2d+'.ieta_iphi', l1CaloTwrColl.ieta[idx], l1CaloTwrColl.iphi[idx])


def save_histos(hm, hm2d, outfile):
    '''
    save all histograms in hm to outfile
    '''
    outfile.mkdir('all_runs')
    outfile.cd('all_runs')
    for varname in hm.get_varnames():
        hm.get(varname).Write()
    for varname in hm2d.get_varnames():
        hm2d.get(varname).Write()
        

def main():
    L1Ana.init_l1_analysis()
    opts = parse_options_upgradeMuonHistos(parser)
    print ""

    # book the histograms
    L1Ana.log.info("Booking combined run histograms.")
    hm, hm2d = book_histograms()

    ntuple = L1Ntuple(opts.nevents)

    if opts.flist:
        ntuple.open_with_file_list(opts.flist)
    if opts.fname:
        ntuple.open_with_file(opts.fname)

    # good runs from json file
    good_ls = None
    if opts.json:
        with open(opts.json) as json_file:    
            good_ls = json.load(json_file)

    # list of runs to run on
    runs_list = []
    if opts.runs:
        runs_list = [int(r) for r in opts.runs.split(",")]
        L1Ana.log.info("Processing only runs:")
        print runs_list

    start_evt = opts.start_event
    end_evt = opts.start_event+ntuple.nevents
    analysed_evt_ctr = 0
    try:
        for i in range(start_evt, end_evt):
            event = ntuple[i]
            if (i+1) % 1000 == 0:
                L1Ana.log.info("Processing event: {n}. Analysed events from selected runs/LS until now: {nAna}".format(n=i+1, nAna=analysed_evt_ctr))

            # if given, only process selected runs
            runnr = event.event.run
            if len(runs_list) > 0 and not runnr in runs_list:
                continue

            # apply json file if loaded
            if good_ls:
                analyze_this_ls = False
                ls = event.event.lumi
                if str(runnr) in good_ls:
                    for ls_list in good_ls[str(runnr)]:
                        if ls >= ls_list[0] and ls <= ls_list[1]:
                            analyze_this_ls = True
                            break
                if not analyze_this_ls:
                    continue

            # now do the analysis
            analyse(event, hm, hm2d)
            analysed_evt_ctr += 1
    except KeyboardInterrupt:
        L1Ana.log.info("Analysis interrupted after {n} events".format(n=i))

    L1Ana.log.info("Analysis of {nAna} events in selected runs/LS finished.".format(nAna=analysed_evt_ctr))

    # save histos to root file
    if saveHistos:
        output = root.TFile(opts.outname, 'recreate')
        output.cd()
        save_histos(hm, hm2d, output)
        output.Close()

if __name__ == "__main__":
    saveHistos = True
    main()

