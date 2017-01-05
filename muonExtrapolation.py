#!/usr/bin/env python
from ToolBox import parse_options_and_init_log
# have to do this first or ROOT masks the -h messages
opts, parser = parse_options_and_init_log()

from L1Analysis import L1Ana, L1Ntuple
from analysis_tools.plotting import HistManager, HistManager2d
from analysis_tools.selections import MuonSelections, Matcher
import exceptions
import ROOT as root

def parse_options_upgradeMuonHistos(parser):
    """
    Adds often used options to the OptionParser...
    """
    parsers = parser.add_subparsers()
    sub_parser = parsers.add_parser("muonExtrapolation")
    sub_parser.add_argument("-o", "--outname", dest="outname", default="./muon_extrapolation_histos.root", type=str, help="A root file name where to save the histograms.")
    sub_parser.add_argument("--pos-charge", dest="pos_charge", default=False, action="store_true", help="Positive probe charge only.")
    sub_parser.add_argument("--neg-charge", dest="neg_charge", default=False, action="store_true", help="Negative probe charge only.")

    opts, unknown = parser.parse_known_args()
    return opts

def book_histograms(eta_ranges):
    # define pt binning
    pt_bins = range(0, 40, 1)
    pt_bins += range(40, 60, 2)
    pt_bins += range(60, 80, 5)
    pt_bins += range(80, 100, 10)
    pt_bins.append(100)

    # 1d histograms and tprofiles
    vars_bins = [['pt_dpt', -1]+pt_bins, ['pt_deta', -1]+pt_bins, ['pt_dphi', -1]+pt_bins]

    x_title_vars = {'pt_dpt':'p_{T}', 'pt_deta':'p_{T}', 'pt_dphi':'p_{T}'}

    x_title_units = {'pt_dpt':'GeV/c', 'pt_deta':'GeV/c', 'pt_dphi':'GeV/c'}

    profile_vars = ['pt_dpt', 'pt_deta', 'pt_dphi']

    # 2d histograms
    x_vars_bins_2d = [['pt', 150, 0, 300]]
    y_vars_bins_2d = [['dcharge', 5, -2, 3]]

    x_title_vars_2d = {'pt':'p_{T}'}
    y_title_vars_2d = {'dcharge':'charge_{L1} - charge_{GEN}'}

    x_title_units_2d = {'pt':'GeV/c'}
    y_title_units_2d = {'dcharge':None}

    varnames = []
    binnings = {}
    profiles = {}
    varnames2d = []
    binnings2d = {}

    for i, eta_range in enumerate(eta_ranges):
        eta_min = eta_range[0]
        eta_max = eta_range[1]
        histoprefix = 'l1_muon_absEtaMin{etaMin}_absEtaMax{etaMax}'.format(etaMin=eta_min, etaMax=eta_max)
        histoprefix2d = '2d_muon_absEtaMin{etaMin}_absEtaMax{etaMax}'.format(etaMin=eta_min, etaMax=eta_max)
        for var_bin in vars_bins:
            if i < 7 or var_bin[0] == 'pt_deta' or var_bin[0] == 'pt_dphi':
                varnames.append(histoprefix+'.{var}'.format(var=var_bin[0]))
                binnings[histoprefix+'.{var}'.format(var=var_bin[0])] = var_bin[1:]+[x_title_vars[var_bin[0]], x_title_units[var_bin[0]]]
                if var_bin[0] in profile_vars:
                    profiles[histoprefix+'.{var}'.format(var=var_bin[0])] = True
        if i < 7:
            for x_var_bin_2d, y_var_bin_2d in zip(x_vars_bins_2d, y_vars_bins_2d):
                varnames2d.append(histoprefix2d+'.{varx}_{vary}'.format(varx=x_var_bin_2d[0], vary=y_var_bin_2d[0]))
                binnings2d[histoprefix2d+'.{varx}_{vary}'.format(varx=x_var_bin_2d[0], vary=y_var_bin_2d[0])] = [x_var_bin_2d[1:]+[x_title_vars_2d[x_var_bin_2d[0]], x_title_units_2d[x_var_bin_2d[0]]], y_var_bin_2d[1:]+[y_title_vars_2d[y_var_bin_2d[0]], y_title_units_2d[y_var_bin_2d[0]]]]

    return HistManager(list(set(varnames)), binnings, profiles), HistManager2d(list(set(varnames2d)), binnings2d)

def analyse(evt, hm, hm2d, eta_ranges):
    genColl = evt.gen
    l1Coll = evt.upgrade

    bx_min = 0
    bx_max = 0

    # eta enlargement of the window for the GEN muons for matching
    gen_extra_eta_range = 0.0435

    gen_muon_idcs = MuonSelections.select_gen_muons(genColl, pt_min=0.5)
    l1_muon_idcs = MuonSelections.select_ugmt_muons(l1Coll, pt_min=0.5, bx_min=bx_min, pos_charge=pos_charge, neg_charge=neg_charge, bx_max=bx_max)

    for i, eta_range in enumerate(eta_ranges):
        eta_min = eta_range[0]
        eta_max = eta_range[1]
        # open the window around the L1 eta range for the GEN muons for matching
        gen_eta_min = eta_min - gen_extra_eta_range
        if gen_eta_min < 0.:
            gen_eta_min = 0.
        gen_eta_max = eta_max + gen_extra_eta_range

        eta_gen_muon_idcs = MuonSelections.select_gen_muons(genColl, abs_eta_min=gen_eta_min, abs_eta_max=gen_eta_max, idcs=gen_muon_idcs)
        eta_l1_muon_idcs = MuonSelections.select_ugmt_muons(l1Coll, abs_eta_min=eta_min, abs_eta_max=eta_max, idcs=l1_muon_idcs)

        matched_muons = Matcher.match_dr(genColl.partEta, genColl.partPhi, l1Coll.muonEta, l1Coll.muonPhi, cut=2., idcs1=eta_gen_muon_idcs, idcs2=eta_l1_muon_idcs)

        histoprefix = 'l1_muon_absEtaMin{etaMin}_absEtaMax{etaMax}'.format(etaMin=eta_min, etaMax=eta_max)
        if i < 7: # only for the first eta ranges
            histoprefix2d = '2d_muon_absEtaMin{etaMin}_absEtaMax{etaMax}'.format(etaMin=eta_min, etaMax=eta_max)
        
        genMuonsUsed = []
        for match in matched_muons:
            if match[0] not in genMuonsUsed:
                genMuonsUsed.append(match[0])
                hm.fill(histoprefix+'.pt_deta', l1Coll.muonEt[match[1]], abs(match[3]))
                hm.fill(histoprefix+'.pt_dphi', l1Coll.muonEt[match[1]], abs(match[4]))
                if i < 7: # fill only for the first eta ranges since histograms are not used for LUT generation
                    hm.fill(histoprefix+'.pt_dpt', l1Coll.muonEt[match[1]], abs(l1Coll.muonEt[match[1]] - genColl.partPt[match[0]]))
                    hm2d.fill(histoprefix2d+'.pt_dcharge', l1Coll.muonEt[match[1]], l1Coll.muonChg[match[1]] - genColl.partCh[match[0]])


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

    global pos_charge
    global neg_charge
    if opts.pos_charge and not opts.neg_charge:
        L1Ana.log.info("Only positive charge requested.")
        pos_charge = True
        neg_charge = False
    elif opts.neg_charge and not opts.pos_charge:
        L1Ana.log.info("Only negative charge requested.")
        pos_charge = False
        neg_charge = True
    elif opts.pos_charge and opts.neg_charge:
        L1Ana.log.warning("Only positive and only negative requested. Will include both charge options.")
        pos_charge = True
        neg_charge = True

#    eta_ranges = [[0, 2.4], [0, 0.83], [0.83, 1.24], [1.24, 2.4]]
    eta_ranges = [[0, 2.4], [0, 0.83], [0.83, 1.24], [1.24, 2.4], [1.2, 1.55], [1.55, 1.85], [1.85, 2.4]]

    # calculate eta ranges
    # The LUT uses a reduced eta coordinate with the two LSBs removed and the MSB masked.
    eta_scale = 0.010875
    eta_bits = 8
    red_eta_bits = 6
    red_eta_scale = 2**(eta_bits - red_eta_bits) * eta_scale
    for red_hw_eta in range(2**red_eta_bits):
        eta_ranges.append((red_hw_eta*red_eta_scale, (red_hw_eta+1)*red_eta_scale))

    # book the histograms
    L1Ana.log.info("Booking combined run histograms.")
    hm, hm2d = book_histograms(eta_ranges)

    ntuple = L1Ntuple(opts.nevents)

    if opts.flist:
        ntuple.open_with_file_list(opts.flist)
    if opts.fname:
        ntuple.open_with_file(opts.fname)

    start_evt = opts.start_event
    end_evt = opts.start_event+ntuple.nevents
    analysed_evt_ctr = 0
    try:
        for i in range(start_evt, end_evt):
            event = ntuple[i]
            if (i+1) % 1000 == 0:
                L1Ana.log.info("Processing event: {n}. Analysed events from selected runs/LS until now: {nAna}".format(n=i+1, nAna=analysed_evt_ctr))

            # now do the analysis
            analyse(event, hm, hm2d, eta_ranges)
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
    pos_charge = True
    neg_charge = True
    saveHistos = True
    main()

