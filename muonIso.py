#!/usr/bin/env python
from ToolBox import parse_options_and_init_log
# have to do this first or ROOT masks the -h messages
opts, parser = parse_options_and_init_log()

from L1Analysis import L1Ana, L1Ntuple
from analysis_tools.plotting import HistManager, HistManager2d
from analysis_tools.selections import MuonSelections, Matcher
from analysis_tools.isolation.caloTowerIso import CaloTowerIsolator
from math import floor, ceil
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

    vars_bins = [['iet', 256, 0, 256], ['ieta', 83, -41, 42], ['iphi', 73, 0, 73], ['iqual', 16, 0, 16], ['n', 2017, 0, 4034], ['total_cone_iet', 256, 0, 256], ['inner_cone_iet', 256, 0, 256], ['outer_cone_iet', 256, 0, 256], ['outer_over_total_cone_iet', 51, 0, 1.02]]
    x_title_vars = {'iet':'iE_{T}', 'ieta':'i#eta', 'iphi':'i#phi', 'iqual':'qual', 'n':'# towers', 'total_cone_iet':'iE_{T}^{total}', 'inner_cone_iet':'iE_{T}^{in}', 'outer_cone_iet':'iE_{T}^{out}', 'outer_over_total_cone_iet':'iE_{T}^{out}'}
    x_title_units = {'iet':None, 'ieta':None, 'iphi':None, 'iqual':None, 'n':None, 'total_cone_iet':None, 'inner_cone_iet':None, 'outer_cone_iet':None, 'outer_over_total_cone_iet':None}

    x_vars_bins_2d = [['ieta', 83, -41, 42], ['iet_ieta', 83, -41, 42], ['iet_ietarel', 165, -82, 83], ['iet_ietarel_red', 31, -15, 16]]
    y_vars_bins_2d = [['iphi', 73, 0, 73], ['iet_iphi', 73, 0, 73], ['iet_iphirel', 73, -36, 37], ['iet_iphirel_red', 31, -15, 16]]

    x_title_vars_2d = {'ieta':'i#eta', 'iet_ieta':'i#eta', 'iet_ietarel':'i#eta_{tower} - i#eta_{#mu}', 'iet_ietarel_red':'i#eta_{tower} - i#eta_{#mu}'}
    y_title_vars_2d = {'iphi':'i#phi', 'iet_iphi':'i#phi', 'iet_iphirel':'i#phi_{tower} - i#phi_{#mu}', 'iet_iphirel_red':'i#phi_{tower} - i#phi_{#mu}'}

    x_title_units_2d = {'ieta':None, 'iet_ieta':None, 'iet_ietarel':None, 'iet_ietarel_red':None}
    y_title_units_2d = {'iphi':None, 'iet_iphi':None, 'iet_iphirel':None, 'iet_iphirel_red':None}

    # the 2d histograms that should be tprofiles
    #profile_2d = {'2d_caloTower.iet_ieta_iet_iphi':True, '2d_caloTower.iet_ietarel_iet_iphirel':True, '2d_caloTower.iet_ietarel_red_iet_iphirel_red':True}

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

    bx_min = 0
    bx_max = 0

    # calo tower histograms
    histoprefix = 'l1_caloTower'
    histoprefix2d = '2d_caloTower'

    nCaloTwr = l1CaloTwrColl.nTower
    hm.fill(histoprefix+'.n', nCaloTwr)
    for i in range(nCaloTwr):
        hm.fill(histoprefix+'.iet', l1CaloTwrColl.iet[i])
        hm.fill(histoprefix+'.ieta', l1CaloTwrColl.ieta[i])
        hm.fill(histoprefix+'.iphi', l1CaloTwrColl.iphi[i])
        hm.fill(histoprefix+'.iqual', l1CaloTwrColl.iqual[i])

        hm2d.fill(histoprefix2d+'.ieta_iphi', l1CaloTwrColl.ieta[i], l1CaloTwrColl.iphi[i])
        hm2d.fill(histoprefix2d+'.iet_ieta_iet_iphi', l1CaloTwrColl.ieta[i], l1CaloTwrColl.iphi[i], l1CaloTwrColl.iet[i])

    # calo towers around L1 muons
    l1_muon_idcs = MuonSelections.select_ugmt_muons(l1MuColl, pt_min=0.5, bx_min=bx_min, bx_max=bx_max)

    for idx in l1_muon_idcs:
        muIEta = l1MuColl.muonIEta[idx]
        muIPhi = l1MuColl.muonIPhi[idx]

        muInCaloTowerIEta = CaloTowerIsolator.calc_muon_calo_tower_ieta(muIEta)
        muInCaloTowerIPhi = CaloTowerIsolator.calc_muon_calo_tower_iphi(muIPhi)

        # get towers relative to muon position and energy sums around
        relTwrs, iEtSums = CaloTowerIsolator.calc_calo_tower_sums(l1CaloTwrColl, muInCaloTowerIEta, muInCaloTowerIPhi, [(1, 1), (5, 5)])
        inner_cone_iet = iEtSums[1]
        total_cone_iet = iEtSums[2]

        for relTwr in relTwrs:
            hm2d.fill(histoprefix2d+'.iet_ietarel_iet_iphirel', relTwr[0], relTwr[1], relTwr[2])
            if abs(relTwr[0]) < 16 and abs(relTwr[1]) < 16:
                hm2d.fill(histoprefix2d+'.iet_ietarel_red_iet_iphirel_red', relTwr[0], relTwr[1], relTwr[2])

        hm.fill(histoprefix+'.total_cone_iet', total_cone_iet)
        hm.fill(histoprefix+'.inner_cone_iet', inner_cone_iet)
        hm.fill(histoprefix+'.outer_cone_iet', total_cone_iet - inner_cone_iet)
        if total_cone_iet > 0:
            hm.fill(histoprefix+'.outer_over_total_cone_iet', (total_cone_iet - inner_cone_iet) / float(total_cone_iet))

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
    muEtaScale = 0.010875
    saveHistos = True
    main()

