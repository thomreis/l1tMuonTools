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
    sub_parser = parsers.add_parser("muonDataEmulComp")
    sub_parser.add_argument("-o", "--outname", dest="outname", default="./muon_data_vs_emul_histos.root", type=str, help="A root file name where to save the histograms.")
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
    # define pt binning
    pt_bins = range(0, 60, 2)
    pt_bins += range(60, 80, 5)
    pt_bins += range(80, 100, 10)
    pt_bins += range(100, 200, 25)
    pt_bins += range(200, 300, 50)
    pt_bins += range(300, 500, 100)
    pt_bins += range(500, 1000, 250)
    pt_bins += range(1000, 1010, 10) # for overflow bin
    pt_bins.append(1010)

    vars_bins = [['pt', -1]+pt_bins, ['eta', 100, -2.5, 2.5], ['phi', 70, -3.5, 3.5], ['charge', 3, -1, 2], ['qual', 16, 0, 15], ['tfMuonIdx', 108, 0, 107], ['n', 9, 0, 9]]
    x_title_vars = {'pt':'p_{T}', 'eta':'#eta', 'phi':'#phi', 'charge':'charge', 'qual':'qual', 'tfMuonIdx':'TF muon index', 'n':'# muons'}
    x_title_units = {'pt':'GeV/c', 'eta':None, 'phi':None, 'charge':None, 'qual':None, 'tfMuonIdx':None, 'n':None}
    x_vars_bins_2d = [['pt', 150, 0, 300]+pt_bins, ['eta', 200, -2.5, 2.5], ['phi', 280, -3.5, 3.5], ['charge', 3, -1, 2], ['qual', 16, 0, 15], ['tfMuonIdx', 108, 0, 107], ['n', 9, 0, 9]]
    y_vars_bins_2d = [['pt', 150, 0, 300]+pt_bins, ['eta', 200, -2.5, 2.5], ['phi', 280, -3.5, 3.5], ['charge', 3, -1, 2], ['qual', 16, 0, 15], ['tfMuonIdx', 108, 0, 107], ['n', 9, 0, 9]]
    x_title_vars_2d = {'pt':'p_{T}^{data}', 'eta':'#eta_{data}', 'phi':'#phi_{data}', 'charge':'charge_{data}', 'qual':'qual_{data}', 'tfMuonIdx':'TF muon index_{data}', 'n':'# muons_{data}'}
    y_title_vars_2d = {'pt':'p_{T}^{emu}', 'eta':'#eta_{emu}', 'phi':'#phi_{emu}', 'charge':'charge_{emu}', 'qual':'qual_{emu}', 'tfMuonIdx':'TF muon index_{emu}', 'n':'# muons_{emul}'}
    x_title_units_2d = {'pt':'GeV/c', 'eta':None, 'phi':None, 'charge':None, 'qual':None, 'tfMuonIdx':None, 'n':None}
    y_title_units_2d = {'pt':'GeV/c', 'eta':None, 'phi':None, 'charge':None, 'qual':None, 'tfMuonIdx':None, 'n':None}

    varnames = []
    binnings = {}
    varnames2d = []
    binnings2d = {}

    for var_bin in vars_bins:
        varnames.append('data_muon.'+var_bin[0])
        varnames.append('emul_muon.'+var_bin[0])

        binnings['data_muon.'+var_bin[0]] = var_bin[1:]+[x_title_vars[var_bin[0]], x_title_units[var_bin[0]]]
        binnings['emul_muon.'+var_bin[0]] = var_bin[1:]+[x_title_vars[var_bin[0]], x_title_units[var_bin[0]]]

    varnames.append('data_emul_matched_muon.dr')
    varnames.append('data_emul_matched_muon.deta')
    varnames.append('data_emul_matched_muon.dphi')
    varnames.append('data_emul_othermatched_muons.dr')
    varnames.append('data_emul_othermatched_muons.deta')
    varnames.append('data_emul_othermatched_muons.dphi')
    binnings['data_emul_matched_muon.dr'] = [800, 0, 8, '#Delta R']
    binnings['data_emul_matched_muon.deta'] = [480, 0, 4.8, '#Delta #eta']
    binnings['data_emul_matched_muon.dphi'] = [630, 0, 6.3, '#Delta #phi']
    binnings['data_emul_othermatched_muons.dr'] = [800, 0, 8, '#Delta R']
    binnings['data_emul_othermatched_muons.deta'] = [480, 0, 4.8, '#Delta #eta']
    binnings['data_emul_othermatched_muons.dphi'] = [630, 0, 6.3, '#Delta #phi']

    for var_bin_2d in x_vars_bins_2d:
        varnames2d.append('2d_data_emul_muon.'+var_bin_2d[0])
        binnings2d['2d_data_emul_muon.'+var_bin_2d[0]] = [var_bin_2d[1:]+[x_title_vars_2d[var_bin_2d[0]], x_title_units_2d[var_bin_2d[0]]], var_bin_2d[1:]+[y_title_vars_2d[var_bin_2d[0]], y_title_units_2d[var_bin_2d[0]]]]
        if var_bin_2d[0] != 'n':
            varnames2d.append('2d_data_emul_matched_muon.'+var_bin_2d[0])
            binnings2d['2d_data_emul_matched_muon.'+var_bin_2d[0]] = [var_bin_2d[1:]+[x_title_vars_2d[var_bin_2d[0]], x_title_units_2d[var_bin_2d[0]]], var_bin_2d[1:]+[y_title_vars_2d[var_bin_2d[0]], y_title_units_2d[var_bin_2d[0]]]]

    return HistManager(list(set(varnames)), binnings), HistManager2d(list(set(varnames2d)), binnings2d)

def analyse(evt, hm, hm2d):
    dataColl = evt.upgrade
    emulColl = evt.upgradeEmu

    bx_min = 0
    bx_max = 0
    data_muon_idcs = MuonSelections.select_ugmt_muons(dataColl, pt_min=0.5, bx_min=bx_min, bx_max=bx_max)
    emul_muon_idcs = MuonSelections.select_ugmt_muons(emulColl, pt_min=0.5, bx_min=bx_min, bx_max=bx_max)

    hm.fill('data_muon.n', len(data_muon_idcs))
    hm.fill('emul_muon.n', len(emul_muon_idcs))
    hm2d.fill('2d_data_emul_muon.n', len(data_muon_idcs), len(emul_muon_idcs))

    for idx in data_muon_idcs:
        hm.fill('data_muon.pt', dataColl.muonEt[idx])
        hm.fill('data_muon.eta', dataColl.muonEta[idx])
        hm.fill('data_muon.phi', dataColl.muonPhi[idx])
        hm.fill('data_muon.charge', dataColl.muonChg[idx])
        hm.fill('data_muon.qual', dataColl.muonQual[idx])
        hm.fill('data_muon.tfMuonIdx', dataColl.muonTfMuonIdx[idx])

    for idx in emul_muon_idcs:
        hm.fill('emul_muon.pt', emulColl.muonEt[idx])
        hm.fill('emul_muon.eta', emulColl.muonEta[idx])
        hm.fill('emul_muon.phi', emulColl.muonPhi[idx])
        hm.fill('emul_muon.charge', emulColl.muonChg[idx])
        hm.fill('emul_muon.qual', emulColl.muonQual[idx])
        hm.fill('emul_muon.tfMuonIdx', emulColl.muonTfMuonIdx[idx])

    for dIdx, eIdx in zip(data_muon_idcs, emul_muon_idcs):
        hm2d.fill('2d_data_emul_muon.pt', dataColl.muonEt[dIdx], emulColl.muonEt[eIdx])
        hm2d.fill('2d_data_emul_muon.eta', dataColl.muonEta[dIdx], emulColl.muonEta[eIdx])
        hm2d.fill('2d_data_emul_muon.phi', dataColl.muonPhi[dIdx], emulColl.muonPhi[eIdx])
        hm2d.fill('2d_data_emul_muon.charge', dataColl.muonChg[dIdx], emulColl.muonChg[eIdx])
        hm2d.fill('2d_data_emul_muon.qual', dataColl.muonQual[dIdx], emulColl.muonQual[eIdx])
        hm2d.fill('2d_data_emul_muon.tfMuonIdx', dataColl.muonTfMuonIdx[dIdx], emulColl.muonTfMuonIdx[eIdx])

    # matching the muons
    matched_muons = Matcher.match_dr(dataColl.muonEta, dataColl.muonPhi, emulColl.muonEta, emulColl.muonPhi, cut=9999, idcs1=data_muon_idcs, idcs2=emul_muon_idcs)
    dataMuonsUsed = []
    emulMuonsUsed = []
    global matched_muon_ctr
    global uncancelled_muon_ctr
    for i in range(len(matched_muons)):
        if matched_muons[i][0] not in dataMuonsUsed and matched_muons[i][1] not in emulMuonsUsed:
            matched_muon_ctr += 1
            dataMuonsUsed.append(matched_muons[i][0])
            emulMuonsUsed.append(matched_muons[i][1])
            hm2d.fill('2d_data_emul_matched_muon.pt', dataColl.muonEt[dataMuonsUsed[-1]], emulColl.muonEt[emulMuonsUsed[-1]])
            hm2d.fill('2d_data_emul_matched_muon.eta', dataColl.muonEta[dataMuonsUsed[-1]], emulColl.muonEta[emulMuonsUsed[-1]])
            hm2d.fill('2d_data_emul_matched_muon.phi', dataColl.muonPhi[dataMuonsUsed[-1]], emulColl.muonPhi[emulMuonsUsed[-1]])
            hm2d.fill('2d_data_emul_matched_muon.charge', dataColl.muonChg[dataMuonsUsed[-1]], emulColl.muonChg[emulMuonsUsed[-1]])
            hm2d.fill('2d_data_emul_matched_muon.qual', dataColl.muonQual[dataMuonsUsed[-1]], emulColl.muonQual[emulMuonsUsed[-1]])
            hm2d.fill('2d_data_emul_matched_muon.tfMuonIdx', dataColl.muonTfMuonIdx[dataMuonsUsed[-1]], emulColl.muonTfMuonIdx[emulMuonsUsed[-1]])
            hm.fill('data_emul_matched_muon.dr', matched_muons[i][2])
            hm.fill('data_emul_matched_muon.deta', matched_muons[i][3])
            hm.fill('data_emul_matched_muon.dphi', matched_muons[i][4])
        else:
            hm.fill('data_emul_othermatched_muons.dr', matched_muons[i][2])
            hm.fill('data_emul_othermatched_muons.deta', matched_muons[i][3])
            hm.fill('data_emul_othermatched_muons.dphi', matched_muons[i][4])
            #if matched_muons[i][2] < 1.e-1:
            #    uncancelled_muon_ctr += 1
            #    dIdx = matched_muons[i][0]
            #    eIdx = matched_muons[i][1]
            #    print 'Found a further match with dR={dr}:'.format(dr=matched_muons[i][2])
            #    print '    data: pt={pt}, eta={eta}, phi={phi}, qual={qual}, muIdx={muIdx}'.format(pt=dataColl.muonEt[dIdx], eta=dataColl.muonEta[dIdx], phi=dataColl.muonPhi[dIdx], qual=dataColl.muonQual[dIdx], muIdx=dataColl.muonTfMuonIdx[dIdx])
            #    print '    emul: pt={pt}, eta={eta}, phi={phi}, qual={qual}, muIdx={muIdx}'.format(pt=emulColl.muonEt[eIdx], eta=emulColl.muonEta[eIdx], phi=emulColl.muonPhi[eIdx], qual=emulColl.muonQual[eIdx], muIdx=emulColl.muonTfMuonIdx[eIdx])
            #    print '    previous matches:'
            #    ctr = 0
            #    for dIdx, eIdx in zip(dataMuonsUsed, emulMuonsUsed):
            #        ctr += 1
            #        print '    {ctr}:'.format(ctr=ctr)
            #        print '    data: pt={pt}, eta={eta}, phi={phi}, qual={qual}, muIdx={muIdx}'.format(pt=dataColl.muonEt[dIdx], eta=dataColl.muonEta[dIdx], phi=dataColl.muonPhi[dIdx], qual=dataColl.muonQual[dIdx], muIdx=dataColl.muonTfMuonIdx[dIdx])
            #        print '    emul: pt={pt}, eta={eta}, phi={phi}, qual={qual}, muIdx={muIdx}'.format(pt=emulColl.muonEt[eIdx], eta=emulColl.muonEta[eIdx], phi=emulColl.muonPhi[eIdx], qual=emulColl.muonQual[eIdx], muIdx=emulColl.muonTfMuonIdx[eIdx])


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
    global matched_muon_ctr
    global uncancelled_muon_ctr 
    L1Ana.log.info("Found {mm} matched muons and {ucm} uncancelled muons with dR < 0.1.".format(mm=matched_muon_ctr, ucm=uncancelled_muon_ctr))

    # save histos to root file
    if saveHistos:
        output = root.TFile(opts.outname, 'recreate')
        output.cd()
        save_histos(hm, hm2d, output)
        output.Close()

if __name__ == "__main__":
    matched_muon_ctr = 0
    uncancelled_muon_ctr = 0
    saveHistos = True
    main()

