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
    sub_parser.add_argument("--pos-eta", dest="pos_eta", default=False, action="store_true", help="Positive GEN eta only.")
    sub_parser.add_argument("--neg-eta", dest="neg_eta", default=False, action="store_true", help="Negative GEN eta only.")
    sub_parser.add_argument("--pos-charge", dest="pos_charge", default=False, action="store_true", help="Positive GEN charge only.")
    sub_parser.add_argument("--neg-charge", dest="neg_charge", default=False, action="store_true", help="Negative GEN charge only.")
    sub_parser.add_argument("--eta-bits", dest="etabits", type=int, default=6, help="Number of eta input bits.")
    sub_parser.add_argument("--tftype", dest="tf", type=str, default='boe', help="Use L1 muons from these TF. BMTF (b), OMTF (o), EMTF (e).")
    sub_parser.add_argument("--emul", dest="emul", default=False, action="store_true", help="Make emulator plots.")

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
    vars_bins = [['deta', 101, -1., 1.], ['dphi', 101, -1., 1.], ['pt_dpt', -1]+pt_bins, ['pt_deta', -1]+pt_bins, ['pt_dphi', -1]+pt_bins, ['pt_absdeta', -1]+pt_bins, ['pt_absdphi', -1]+pt_bins]

    x_title_vars = {'deta':'#eta_{L1} - #eta_{GEN}', 'dphi':'#phi_{L1} - #phi_{GEN}', 'pt_dpt':'p_{T}^{L1}', 'pt_deta':'p_{T}^{L1}', 'pt_dphi':'p_{T}^{L1}', 'pt_absdeta':'p_{T}^{L1}', 'pt_absdphi':'p_{T}^{L1}'}

    x_title_units = {'deta':None, 'dphi':None, 'pt_dpt':'GeV/c', 'pt_deta':'GeV/c', 'pt_dphi':'GeV/c', 'pt_absdeta':'GeV/c', 'pt_absdphi':'GeV/c'}

    profile_vars = ['pt_dpt', 'pt_deta', 'pt_dphi', 'pt_absdeta', 'pt_absdphi']

    # 2d histograms
    x_vars_bins_2d = [['pt', 150, 0, 300], ['pt', 150, 0, 300], ['pt', 150, 0, 300]]
    y_vars_bins_2d = [['dcharge', 5, -2, 3], ['deta', 100, -0.2, 0.2], ['dphi', 100, -1., 1.]]

    x_title_vars_2d = {'pt':'p_{T}^{L1}', 'pt':'p_{T}^{L1}', 'pt':'p_{T}^{L1}'}
    y_title_vars_2d = {'dcharge':'charge_{L1} - charge_{GEN}', 'deta':'#eta_{L1} - #eta_{GEN}', 'dphi':'#phi_{L1} - #phi_{GEN}'}

    x_title_units_2d = {'pt':'GeV/c', 'pt':'GeV/c', 'pt':'GeV/c'}
    y_title_units_2d = {'dcharge':None, 'deta':None, 'dphi':None}

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
        histoprefix_extrapol = 'l1_muon_extrapol_absEtaMin{etaMin}_absEtaMax{etaMax}'.format(etaMin=eta_min, etaMax=eta_max)
        histoprefix2d_extrapol = '2d_muon_extrapol_absEtaMin{etaMin}_absEtaMax{etaMax}'.format(etaMin=eta_min, etaMax=eta_max)
        for var_bin in vars_bins:
            if i < 7 or var_bin[0] == 'pt_deta' or var_bin[0] == 'pt_dphi' or var_bin[0] == 'pt_absdeta' or var_bin[0] == 'pt_absdphi':
                varnames.append(histoprefix+'.{var}'.format(var=var_bin[0]))
                varnames.append(histoprefix_extrapol+'.{var}'.format(var=var_bin[0]))
                binnings[histoprefix+'.{var}'.format(var=var_bin[0])] = var_bin[1:]+[x_title_vars[var_bin[0]], x_title_units[var_bin[0]]]
                binnings[histoprefix_extrapol+'.{var}'.format(var=var_bin[0])] = var_bin[1:]+[x_title_vars[var_bin[0]], x_title_units[var_bin[0]]]
                if var_bin[0] in profile_vars:
                    profiles[histoprefix+'.{var}'.format(var=var_bin[0])] = True
                    profiles[histoprefix_extrapol+'.{var}'.format(var=var_bin[0])] = True
        if i < 7:
            for x_var_bin_2d, y_var_bin_2d in zip(x_vars_bins_2d, y_vars_bins_2d):
                varnames2d.append(histoprefix2d+'.{varx}_{vary}'.format(varx=x_var_bin_2d[0], vary=y_var_bin_2d[0]))
                varnames2d.append(histoprefix2d_extrapol+'.{varx}_{vary}'.format(varx=x_var_bin_2d[0], vary=y_var_bin_2d[0]))
                binnings2d[histoprefix2d+'.{varx}_{vary}'.format(varx=x_var_bin_2d[0], vary=y_var_bin_2d[0])] = [x_var_bin_2d[1:]+[x_title_vars_2d[x_var_bin_2d[0]], x_title_units_2d[x_var_bin_2d[0]]], y_var_bin_2d[1:]+[y_title_vars_2d[y_var_bin_2d[0]], y_title_units_2d[y_var_bin_2d[0]]]]
                binnings2d[histoprefix2d_extrapol+'.{varx}_{vary}'.format(varx=x_var_bin_2d[0], vary=y_var_bin_2d[0])] = [x_var_bin_2d[1:]+[x_title_vars_2d[x_var_bin_2d[0]], x_title_units_2d[x_var_bin_2d[0]]], y_var_bin_2d[1:]+[y_title_vars_2d[y_var_bin_2d[0]], y_title_units_2d[y_var_bin_2d[0]]]]

    return HistManager(list(set(varnames)), binnings, profiles), HistManager2d(list(set(varnames2d)), binnings2d)

def analyse(evt, hm, hm2d, eta_ranges, emul=False, tf='boe'):
    genColl = evt.gen
    if emul:
        l1Coll = evt.upgradeEmu
    else:
        l1Coll = evt.upgrade

    tftype = []
    if tf.find('b') != -1:
        tftype.append(0)
    if tf.find('o') != -1:
        tftype.append(1)
    if tf.find('e') != -1:
        tftype.append(2)

    bx_min = 0
    bx_max = 0

    # eta enlargement of the window for the GEN muons for matching
    gen_extra_eta_range = 0.0435

    neg_eta_only = neg_eta and not pos_eta

    gen_muon_idcs = MuonSelections.select_gen_muons(genColl, pt_min=0.5, pos_eta=pos_eta, neg_eta=neg_eta, pos_charge=pos_charge, neg_charge=neg_charge)
    #l1_muon_idcs = MuonSelections.select_ugmt_muons(l1Coll, pt_min=0.5, bx_min=bx_min, bx_max=bx_max, tftype=tftype, qual_min=8)
    l1_muon_idcs = MuonSelections.select_ugmt_muons(l1Coll, pt_min=0.5, bx_min=bx_min, bx_max=bx_max, tftype=tftype)

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
        eta_extrapol_l1_muon_idcs = MuonSelections.select_ugmt_muons(l1Coll, abs_eta_min=eta_min, abs_eta_max=eta_max, idcs=l1_muon_idcs, useVtxExtraCoord=True)

        #if len(eta_l1_muon_idcs) != 2 or len(l1Coll.muonEta) != 2:
        #    continue

        #if l1Coll.muonEta[0] == l1Coll.muonEta[1]:
        #    continue

        matched_muons = Matcher.match_dr(l1Coll.muonEta, l1Coll.muonPhi, genColl.partEta, genColl.partPhi, cut=2., idcs1=eta_l1_muon_idcs, idcs2=eta_gen_muon_idcs)
        matched_muons_extrapol = Matcher.match_dr(l1Coll.muonEtaAtVtx, l1Coll.muonPhiAtVtx, genColl.partEta, genColl.partPhi, cut=2., idcs1=eta_extrapol_l1_muon_idcs, idcs2=eta_gen_muon_idcs)

        histoprefix = 'l1_muon_absEtaMin{etaMin}_absEtaMax{etaMax}'.format(etaMin=eta_min, etaMax=eta_max)
        histoprefix_extrapol = 'l1_muon_extrapol_absEtaMin{etaMin}_absEtaMax{etaMax}'.format(etaMin=eta_min, etaMax=eta_max)
        if i < 7: # only for the first eta ranges
            histoprefix2d = '2d_muon_absEtaMin{etaMin}_absEtaMax{etaMax}'.format(etaMin=eta_min, etaMax=eta_max)
            histoprefix2d_extrapol = '2d_muon_extrapol_absEtaMin{etaMin}_absEtaMax{etaMax}'.format(etaMin=eta_min, etaMax=eta_max)

        genMuonsUsed = []
        for match in matched_muons:
            if match[1] not in genMuonsUsed:
                genMuonsUsed.append(match[1])
                if neg_eta_only:
                    hm.fill(histoprefix+'.pt_deta', l1Coll.muonEt[match[0]], -1*match[3])
                else:
                    hm.fill(histoprefix+'.pt_deta', l1Coll.muonEt[match[0]], match[3])
                hm.fill(histoprefix+'.pt_dphi', l1Coll.muonEt[match[0]], match[4])
                hm.fill(histoprefix+'.pt_absdeta', l1Coll.muonEt[match[0]], abs(match[3]))
                hm.fill(histoprefix+'.pt_absdphi', l1Coll.muonEt[match[0]], abs(match[4]))
                #if l1Coll.muonEt[match[0]] < 8. and abs(match[4]) < 0.05:
                #    print matched_muons
                #    for match2 in matched_muons:
                #        print 'strange muon: pT L1: {ptl}, pT GEN: {ptg}, dr: {dr}, tfi: {tfi}, eta L1: {etal}, eta GEN: {etag}'.format(ptl=l1Coll.muonEt[match2[0]], ptg=genColl.partPt[match2[1]], dr=match2[2], tfi=l1Coll.muonTfMuonIdx[match2[0]], etal=l1Coll.muonEta[match2[0]], etag=genColl.partEta[match2[1]])
                if i < 7: # fill only for the first eta ranges since histograms are not used for LUT generation
                    #if match[3] < -0.7+1000:
                    #    print matched_muons
                    #    for match2 in matched_muons:
                    #        print 'strange muon: pT L1: {ptl}, pT GEN: {ptg}, dr: {dr}, tfi: {tfi}, eta L1: {etal}, eta GEN: {etag}, phi L1: {phil}, phi GEN: {phig}'.format(ptl=l1Coll.muonEt[match2[0]], ptg=genColl.partPt[match2[1]], dr=match2[2], tfi=l1Coll.muonTfMuonIdx[match2[0]], etal=l1Coll.muonEta[match2[0]], etag=genColl.partEta[match2[1]], phil=l1Coll.muonPhi[match2[0]], phig=genColl.partPhi[match2[1]])
                    hm.fill(histoprefix+'.deta', match[3])
                    hm.fill(histoprefix+'.dphi', match[4])
                    hm.fill(histoprefix+'.pt_dpt', l1Coll.muonEt[match[0]], abs(l1Coll.muonEt[match[0]] - genColl.partPt[match[1]]))
                    hm2d.fill(histoprefix2d+'.pt_dcharge', l1Coll.muonEt[match[0]], l1Coll.muonChg[match[0]] - genColl.partCh[match[1]])
                    hm2d.fill(histoprefix2d+'.pt_deta', l1Coll.muonEt[match[0]], match[3])
                    hm2d.fill(histoprefix2d+'.pt_dphi', l1Coll.muonEt[match[0]], match[4])

        genMuonsUsed = []
        for match in matched_muons_extrapol:
            if match[1] not in genMuonsUsed:
                genMuonsUsed.append(match[1])
                if neg_eta_only:
                    hm.fill(histoprefix_extrapol+'.pt_deta', l1Coll.muonEt[match[0]], -1*match[3])
                else:
                    hm.fill(histoprefix_extrapol+'.pt_deta', l1Coll.muonEt[match[0]], match[3])
                hm.fill(histoprefix_extrapol+'.pt_dphi', l1Coll.muonEt[match[0]], match[4])
                hm.fill(histoprefix_extrapol+'.pt_absdeta', l1Coll.muonEt[match[0]], abs(match[3]))
                hm.fill(histoprefix_extrapol+'.pt_absdphi', l1Coll.muonEt[match[0]], abs(match[4]))
                if i < 7: # fill only for the first eta ranges since histograms are not used for LUT generation
                    hm.fill(histoprefix_extrapol+'.deta', match[3])
                    hm.fill(histoprefix_extrapol+'.dphi', match[4])
                    hm.fill(histoprefix_extrapol+'.pt_dpt', l1Coll.muonEt[match[0]], abs(l1Coll.muonEt[match[0]] - genColl.partPt[match[1]]))
                    hm2d.fill(histoprefix2d_extrapol+'.pt_dcharge', l1Coll.muonEt[match[0]], l1Coll.muonChg[match[0]] - genColl.partCh[match[1]])
                    hm2d.fill(histoprefix2d_extrapol+'.pt_deta', l1Coll.muonEt[match[0]], match[3])
                    hm2d.fill(histoprefix2d_extrapol+'.pt_dphi', l1Coll.muonEt[match[0]], match[4])


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

    global pos_eta
    global neg_eta
    if opts.pos_eta and not opts.neg_eta:
        L1Ana.log.info("Only positive eta requested.")
        pos_eta = True
        neg_eta = False
    elif opts.neg_eta and not opts.pos_eta:
        L1Ana.log.info("Only negative eta requested.")
        pos_eta = False
        neg_eta = True
    elif opts.pos_eta and opts.neg_eta:
        L1Ana.log.warning("Only positive and only negative requested. Will include both eta options.")
        pos_eta = True
        neg_eta = True

    emul = opts.emul

#    eta_ranges = [[0, 2.4], [0, 0.83], [0.83, 1.24], [1.24, 2.4]]
    eta_ranges = [[0, 2.4], [0, 0.83], [0.83, 1.24], [1.24, 2.4], [1.2, 1.55], [1.55, 1.85], [1.85, 2.4]]
    #eta_ranges = [[0, 2.4]]

    # calculate eta ranges
    # The LUT uses a reduced eta coordinate with the two LSBs removed and the MSB masked.
    eta_scale = 0.010875
    eta_bits = 8
    red_eta_bits = opts.etabits
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
            analyse(event, hm, hm2d, eta_ranges, emul, opts.tf)
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
    pos_eta = True
    neg_eta = True
    pos_charge = True
    neg_charge = True
    saveHistos = True
    main()

