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
    sub_parser = parsers.add_parser("muonKinematics")
    sub_parser.add_argument("-o", "--outname", dest="outname", default="./muon_kinematics_histos.root", type=str, help="A root file name where to save the histograms.")
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
    pt_bins.append(300)

    vars_bins = [['pt', -1]+pt_bins, ['eta', 100, -2.5, 2.5], ['etaVtx', 100, -2.5, 2.5], ['phi', 70, -3.5, 3.5], ['phiVtx', 70, -3.5, 3.5], ['charge', 3, -1, 2], ['qual', 16, 0, 15], ['tfMuonIdx', 108, 0, 107],
                 ['mu1pt', -1]+pt_bins, ['mu1eta', 100, -2.5, 2.5], ['mu1etaVtx', 100, -2.5, 2.5], ['mu1phi', 70, -3.5, 3.5], ['mu1phiVtx', 70, -3.5, 3.5], ['mu1charge', 3, -1, 2], ['mu1qual', 16, 0, 15], ['mu1tfMuonIdx', 108, 0, 107],
                 ['mu2pt', -1]+pt_bins, ['mu2eta', 100, -2.5, 2.5], ['mu2etaVtx', 100, -2.5, 2.5], ['mu2phi', 70, -3.5, 3.5], ['mu2phiVtx', 70, -3.5, 3.5], ['mu2charge', 3, -1, 2], ['mu2qual', 16, 0, 15], ['mu2tfMuonIdx', 108, 0, 107],
                 ['mu3pt', -1]+pt_bins, ['mu3eta', 100, -2.5, 2.5], ['mu3etaVtx', 100, -2.5, 2.5], ['mu3phi', 70, -3.5, 3.5], ['mu3phiVtx', 70, -3.5, 3.5], ['mu3charge', 3, -1, 2], ['mu3qual', 16, 0, 15], ['mu3tfMuonIdx', 108, 0, 107],
                 ['n', 9, 0, 9], ['dpt', 300, 0, 300], ['dptoverpt', 50, 0, 1], ['dr', 600, 0, 6], ['deta', 480, 0, 4.8], ['dphi', 320, 0, 3.2]]

    x_title_vars = {'pt':'p_{T}', 'eta':'#eta', 'etaVtx':'#eta_{Vtx}', 'phi':'#phi', 'phiVtx':'#phi_{Vtx}', 'charge':'charge', 'qual':'qual', 'tfMuonIdx':'TF muon index',
                    'mu1pt':'p_{T}^{#mu1}', 'mu1eta':'#eta^{#mu1}', 'mu1etaVtx':'#eta_{Vtx}^{#mu1}', 'mu1phi':'#phi^{#mu1}', 'mu1phiVtx':'#phi_{Vtx}^{#mu1}', 'mu1charge':'charge^{#mu1}', 'mu1qual':'qual^{#mu1}', 'mu1tfMuonIdx':'TF muon index^{#mu1}',
                    'mu2pt':'p_{T}^{#mu2}', 'mu2eta':'#eta^{#mu2}', 'mu2etaVtx':'#eta_{Vtx}^{#mu2}', 'mu2phi':'#phi^{#mu2}', 'mu2phiVtx':'#phi_{Vtx}^{#mu2}', 'mu2charge':'charge^{#mu2}', 'mu2qual':'qual^{#mu2}', 'mu2tfMuonIdx':'TF muon index^{#mu2}',
                    'mu3pt':'p_{T}^{#mu3}', 'mu3eta':'#eta^{#mu3}', 'mu3etaVtx':'#eta_{Vtx}^{#mu3}', 'mu3phi':'#phi^{#mu3}', 'mu3phiVtx':'#phi_{Vtx}^{#mu3}', 'mu3charge':'charge^{#mu3}', 'mu3qual':'qual^{#mu3}', 'mu3tfMuonIdx':'TF muon index^{#mu3}',
                    'n':'# muons', 'dpt':'p_{T}^{#mu1} - p_{T}^{#mu2}', 'dptoverpt':'(p_{T}^{#mu1} - p_{T}^{#mu2}) / p_{T}^{#mu1}', 'dr':'#Delta R (#mu1, #mu2)', 'deta':'|#Delta#eta (#mu1, #mu2)|', 'dphi':'#Delta#phi (#mu1, #mu2)'}

    x_title_units = {'pt':'GeV/c', 'eta':None, 'etaVtx':None, 'phi':None, 'phiVtx':None, 'charge':None, 'qual':None, 'tfMuonIdx':None,
                     'mu1pt':'GeV/c', 'mu1eta':None, 'mu1etaVtx':None, 'mu1phi':None, 'mu1phiVtx':None, 'mu1charge':None, 'mu1qual':None, 'mu1tfMuonIdx':None,
                     'mu2pt':'GeV/c', 'mu2eta':None, 'mu2etaVtx':None, 'mu2phi':None, 'mu2phiVtx':None, 'mu2charge':None, 'mu2qual':None, 'mu2tfMuonIdx':None,
                     'mu3pt':'GeV/c', 'mu3eta':None, 'mu3etaVtx':None, 'mu3phi':None, 'mu3phiVtx':None, 'mu3charge':None, 'mu3qual':None, 'mu3tfMuonIdx':None,
                     'n':None, 'dpt':'GeV/c', 'dptoverpt':None, 'dr':None, 'deta':None, 'dphi':None}

    x_vars_bins_2d = [['pt', 150, 0, 300]+pt_bins, ['eta', 100, -2.5, 2.5], ['phi', 140, -3.5, 3.5], ['charge', 3, -1, 2], ['qual', 16, 0, 15], ['tfMuonIdx', 108, 0, 107]]
    y_vars_bins_2d = [['pt', 150, 0, 300]+pt_bins, ['eta', 100, -2.5, 2.5], ['phi', 140, -3.5, 3.5], ['charge', 3, -1, 2], ['qual', 16, 0, 15], ['tfMuonIdx', 108, 0, 107]]

    x_title_vars_2d = {'pt':'p_{T}^{#mu1}', 'eta':'#eta_{#mu1}', 'phi':'#phi_{#mu1}', 'charge':'charge_{#mu1}', 'qual':'qual_{#mu1}', 'tfMuonIdx':'TF muon index_{#mu1}'}
    y_title_vars_2d = {'pt':'p_{T}^{#mu2}', 'eta':'#eta_{#mu2}', 'phi':'#phi_{#mu2}', 'charge':'charge_{#mu2}', 'qual':'qual_{#mu2}', 'tfMuonIdx':'TF muon index_{#mu2}'}

    x_title_units_2d = {'pt':'GeV/c', 'eta':None, 'phi':None, 'charge':None, 'qual':None, 'tfMuonIdx':None}
    y_title_units_2d = {'pt':'GeV/c', 'eta':None, 'phi':None, 'charge':None, 'qual':None, 'tfMuonIdx':None}

    varnames = []
    binnings = {}
    varnames2d = []
    binnings2d = {}

    quals = [4, 8, 12]

    for qual in quals:
        for var_bin in vars_bins:
            varnames.append('l1_muon_q{q}.{var}'.format(q=qual, var=var_bin[0]))
            varnames.append('l1_muon_qmin{q}.{var}'.format(q=qual, var=var_bin[0]))
            binnings['l1_muon_q{q}.{var}'.format(q=qual, var=var_bin[0])] = var_bin[1:]+[x_title_vars[var_bin[0]], x_title_units[var_bin[0]]]
            binnings['l1_muon_qmin{q}.{var}'.format(q=qual, var=var_bin[0])] = var_bin[1:]+[x_title_vars[var_bin[0]], x_title_units[var_bin[0]]]

        for var_bin_2d in x_vars_bins_2d:
            varnames2d.append('2d_muon_qmin{q}.{var}'.format(q=qual, var=var_bin_2d[0]))
            binnings2d['2d_muon_qmin{q}.{var}'.format(q=qual, var=var_bin_2d[0])] = [var_bin_2d[1:]+[x_title_vars_2d[var_bin_2d[0]], x_title_units_2d[var_bin_2d[0]]], var_bin_2d[1:]+[y_title_vars_2d[var_bin_2d[0]], y_title_units_2d[var_bin_2d[0]]]]

    return HistManager(list(set(varnames)), binnings), HistManager2d(list(set(varnames2d)), binnings2d)

def analyse(evt, hm, hm2d):
    l1Coll = evt.upgrade

    bx_min = 0
    bx_max = 0

    quals = [4, 8, 12]

    for qual in quals:
        l1_muon_idcs = MuonSelections.select_ugmt_muons(l1Coll, pt_min=0.5, qual_min=qual, qual_max=qual, bx_min=bx_min, bx_max=bx_max)
        l1_muon_idcs_qmin = MuonSelections.select_ugmt_muons(l1Coll, pt_min=0.5, qual_min=qual, bx_min=bx_min, bx_max=bx_max)

        histoprefix = 'l1_muon_q{q}'.format(q=qual)
        histoprefixqmin = 'l1_muon_qmin{q}'.format(q=qual)

        nMu = len(l1_muon_idcs)
        nMuQmin = len(l1_muon_idcs_qmin)
        hm.fill(histoprefix+'.n', nMu)

        for idx in l1_muon_idcs:
            hm.fill(histoprefix+'.pt', l1Coll.muonEt[idx])
            hm.fill(histoprefix+'.eta', l1Coll.muonEta[idx])
            hm.fill(histoprefix+'.phi', l1Coll.muonPhi[idx])
            hm.fill(histoprefix+'.charge', l1Coll.muonChg[idx])
            hm.fill(histoprefix+'.qual', l1Coll.muonQual[idx])
            hm.fill(histoprefix+'.tfMuonIdx', l1Coll.muonTfMuonIdx[idx])

        for nMin in range(3):
            if nMu > nMin:
                hm.fill(histoprefix+'.mu{nMu}pt'.format(nMu=nMin+1), l1Coll.muonEt[l1_muon_idcs[nMin]])
                hm.fill(histoprefix+'.mu{nMu}eta'.format(nMu=nMin+1), l1Coll.muonEta[l1_muon_idcs[nMin]])
                hm.fill(histoprefix+'.mu{nMu}phi'.format(nMu=nMin+1), l1Coll.muonPhi[l1_muon_idcs[nMin]])
                hm.fill(histoprefix+'.mu{nMu}charge'.format(nMu=nMin+1), l1Coll.muonChg[l1_muon_idcs[nMin]])
                hm.fill(histoprefix+'.mu{nMu}qual'.format(nMu=nMin+1), l1Coll.muonQual[l1_muon_idcs[nMin]])
                hm.fill(histoprefix+'.mu{nMu}tfMuonIdx'.format(nMu=nMin+1), l1Coll.muonTfMuonIdx[l1_muon_idcs[nMin]])
            if nMuQmin > nMin:
                hm.fill(histoprefixqmin+'.mu{nMu}pt'.format(nMu=nMin+1), l1Coll.muonEt[l1_muon_idcs_qmin[nMin]])
                hm.fill(histoprefixqmin+'.mu{nMu}eta'.format(nMu=nMin+1), l1Coll.muonEta[l1_muon_idcs_qmin[nMin]])
                hm.fill(histoprefixqmin+'.mu{nMu}phi'.format(nMu=nMin+1), l1Coll.muonPhi[l1_muon_idcs_qmin[nMin]])
                hm.fill(histoprefixqmin+'.mu{nMu}charge'.format(nMu=nMin+1), l1Coll.muonChg[l1_muon_idcs_qmin[nMin]])
                hm.fill(histoprefixqmin+'.mu{nMu}qual'.format(nMu=nMin+1), l1Coll.muonQual[l1_muon_idcs_qmin[nMin]])
                hm.fill(histoprefixqmin+'.mu{nMu}tfMuonIdx'.format(nMu=nMin+1), l1Coll.muonTfMuonIdx[l1_muon_idcs_qmin[nMin]])
 
        if nMu > 1:
            hm.fill(histoprefix+'.dpt', l1Coll.muonEt[l1_muon_idcs[0]] - l1Coll.muonEt[l1_muon_idcs[1]])
            hm.fill(histoprefix+'.dptoverpt', (l1Coll.muonEt[l1_muon_idcs[0]] - l1Coll.muonEt[l1_muon_idcs[1]]) / l1Coll.muonEt[l1_muon_idcs[0]])
            hm.fill(histoprefix+'.dr', Matcher.delta_r(l1Coll.muonPhi[l1_muon_idcs[0]], l1Coll.muonEta[l1_muon_idcs[0]], l1Coll.muonPhi[l1_muon_idcs[1]], l1Coll.muonEta[l1_muon_idcs[1]]))
            hm.fill(histoprefix+'.deta', abs(l1Coll.muonEta[l1_muon_idcs[0]] - l1Coll.muonEta[l1_muon_idcs[1]]))
            hm.fill(histoprefix+'.dphi', Matcher.delta_phi(l1Coll.muonPhi[l1_muon_idcs[0]], l1Coll.muonPhi[l1_muon_idcs[1]]))

        histoprefix2d = '2d_muon_qmin{q}'.format(q=qual)
        if len(l1_muon_idcs_qmin) > 1:
            hm2d.fill(histoprefix2d+'.pt', l1Coll.muonEt[l1_muon_idcs_qmin[0]], l1Coll.muonEt[l1_muon_idcs_qmin[1]])
            hm2d.fill(histoprefix2d+'.eta', l1Coll.muonEta[l1_muon_idcs_qmin[0]], l1Coll.muonEta[l1_muon_idcs_qmin[1]])
            hm2d.fill(histoprefix2d+'.phi', l1Coll.muonPhi[l1_muon_idcs_qmin[0]], l1Coll.muonPhi[l1_muon_idcs_qmin[1]])
            hm2d.fill(histoprefix2d+'.charge', l1Coll.muonChg[l1_muon_idcs_qmin[0]], l1Coll.muonChg[l1_muon_idcs_qmin[1]])
            hm2d.fill(histoprefix2d+'.qual', l1Coll.muonQual[l1_muon_idcs_qmin[0]], l1Coll.muonQual[l1_muon_idcs_qmin[1]])
            hm2d.fill(histoprefix2d+'.tfMuonIdx', l1Coll.muonTfMuonIdx[l1_muon_idcs_qmin[0]], l1Coll.muonTfMuonIdx[l1_muon_idcs_qmin[1]])


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

