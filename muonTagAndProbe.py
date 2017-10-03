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
    sub_parser = parsers.add_parser("muonTagAndProbe")
    sub_parser.add_argument("-o", "--outname", dest="outname", default="./ugmt_tandp_eff_histos.root", type=str, help="A root file name where to save the histograms.")
    sub_parser.add_argument("-j", "--json", dest="json", type=str, default=None, help="A json file with good lumi sections per run.")
    sub_parser.add_argument("-r", "--runs", dest="runs", type=str, default=None, help="A string of runs to check.")
    sub_parser.add_argument("-p", "--pos-side", dest="pos_side", default=False, action="store_true", help="Positive detector side only.")
    sub_parser.add_argument("-n", "--neg-side", dest="neg_side", default=False, action="store_true", help="Negative detector side only.")
    sub_parser.add_argument("--pos-charge", dest="pos_charge", default=False, action="store_true", help="Positive probe charge only.")
    sub_parser.add_argument("--neg-charge", dest="neg_charge", default=False, action="store_true", help="Negative probe charge only.")
    sub_parser.add_argument("--use-inv-mass-cut", dest="invmasscut", default=False, action="store_true", help="Use an invariant mass range for the tag and probe pair.")
    sub_parser.add_argument("--use-l1-extra-coord", dest="l1extraCoord", default=False, action="store_true", help="Use L1 extrapolated eta and phi coordinates.")
    sub_parser.add_argument("--use-reco-extra-station", dest="recoExtraStation", type=int, default=0, help="Extrapolated reco muon coordinates. 0=Vertex, 1=1st muon station, 2=2nd muon station.")
    sub_parser.add_argument("--emul", dest="emul", default=False, action="store_true", help="Make emulator plots.")
    sub_parser.add_argument("--legacy", dest="legacy", default=False, action="store_true", help="Use legacy muons translated to upgrade format.")
    sub_parser.add_argument("--pa", dest="pa_run", default=False, action="store_true", help="Setup for pA run.")
    sub_parser.add_argument("--prefix", dest="prefix", type=str, default='', help="A prefix for the histogram names.")
    sub_parser.add_argument("--tftype", dest="tftype", type=str, default='', help="Fill L1 muons from one TF.")
    sub_parser.add_argument("--era", dest="era", type=str, default='2017pp', help="Era to select run numbers for plots.")
    sub_parser.add_argument("--pt-ranges", dest="ptranges", type=str, default='standard', help="A set of pT cuts to make plots for ['standard', 'extended'].")
    sub_parser.add_argument("--eta-ranges", dest="etaranges", type=str, default='standard', help="A set of eta ranges to make plots for ['minimal', 'standard', 'extended'].")

    opts, unknown = parser.parse_known_args()
    return opts

def get_tftype(tf_muon_index):
    if tf_muon_index > 35 and tf_muon_index < 72:
        return 0 # BMTF
    elif tf_muon_index > 17 and tf_muon_index < 90:
        return 1 # OMTF
    else:
        return 2 # EMTF

def book_histograms(eta_ranges, qual_ptmins_dict, res_probe_ptmins, match_deltas, emul=False, legacy=False):
    # define pt binning
    pt_bins = range(0, 30, 1)
    pt_bins = range(30, 50, 2)
    pt_bins += range(50, 80, 5)
    pt_bins += range(80, 100, 10)
    pt_bins += range(100, 200, 25)
    pt_bins += range(200, 300, 50)
    pt_bins += range(300, 500, 100)
    pt_bins += range(500, 1000, 250)
    pt_bins += range(1000, 1010, 10) # for overflow bin
    pt_bins.append(1010)

    p_bins = pt_bins

    vars_bins = [['pt', -1]+pt_bins, ['eta', 100, -2.5, 2.5], ['phi', 70, -3.5, 3.5], ['charge', 3, -1, 2], ['vtx', 60, 0, 60]]
    if era == '2016pp':
        vars_bins.append(['run', 13100, 271000, 284100])
    elif era == '2016pPb':
        vars_bins.append(['run', 2500, 284100, 286600])
    elif era == '2017pp':
        vars_bins.append(['run', 25000, 294645, 319645])
    else:
        vars_bins.append(['run', 10, 0, 10])
    x_title_vars = {'pt':'p_{T}', 'eta':'#eta', 'phi':'#phi', 'charge':'charge', 'vtx':'PU', 'run':'run number'}
    x_title_units = {'pt':'GeV/c', 'eta':None, 'phi':None, 'charge':None, 'vtx':None, 'run':None}
    probe_vars_bins = vars_bins + [['p', -1]+p_bins]
    probe_x_title_vars = {'p':'p'}
    probe_x_title_vars.update(x_title_vars)
    probe_x_title_units = {'p':'GeV/c'}
    probe_x_title_units.update(x_title_units)
    res_vars_bins = [['dpt', 100, -50, 50], ['dinvpt', 200, -2., 2.], ['deta', 100, -0.1, 0.1], ['dphi', 100, -0.2, 0.2], ['dcharge', 5, -2., 3]]
    res_x_title_vars = {'dpt':'p_{T}^{L1} - p_{T}^{reco}', 'dinvpt':'(p_{T}^{reco} - p_{T}^{L1}) / p_{T}^{L1}', 'deta':'#eta_{L1} - #eta_{reco}', 'dphi':'#phi_{L1} - #phi_{reco}', 'dcharge':'charge^{L1} - charge^{reco}'}
    res_x_title_units = {'dpt':'GeV', 'dinvpt':'GeV', 'deta':None, 'dphi':None, 'dcharge':None}
    #x_vars_bins_2d = [['pt', -1]+pt_bins, ['eta', 100, -2.5, 2.5], ['phi', 70, -3.5, 3.5]]
    #y_vars_bins_2d = [['pt', -1]+pt_bins, ['eta', 100, -2.5, 2.5], ['phi', 70, -3.5, 3.5]]
    #x_vars_bins_2d = [['pt', 150, 0, 300]+pt_bins, ['eta', 100, -2.5, 2.5], ['phi', 70, -3.5, 3.5]]
    #y_vars_bins_2d = [['pt', 150, 0, 300]+pt_bins, ['eta', 100, -2.5, 2.5], ['phi', 70, -3.5, 3.5]]
    x_vars_bins_2d = [['pt', 150, 0, 300]+pt_bins, ['eta', 200, -2.5, 2.5], ['phi', 280, -3.5, 3.5], ['charge', 3, -1, 2]]
    y_vars_bins_2d = [['pt', 150, 0, 300]+pt_bins, ['eta', 200, -2.5, 2.5], ['phi', 280, -3.5, 3.5], ['charge', 3, -1, 2]]
    x_title_vars_2d = {'pt':'p_{T}^{reco}', 'eta':'#eta_{reco}', 'phi':'#phi_{reco}', 'charge':'charge_{reco}'}
    y_title_vars_2d = {'pt':'p_{T}^{L1}', 'eta':'#eta_{L1}', 'phi':'#phi_{L1}', 'charge':'charge_{L1}'}
    x_title_units_2d = {'pt':'GeV/c', 'eta':None, 'phi':None, 'charge':None}
    y_title_units_2d = {'pt':'GeV/c', 'eta':None, 'phi':None, 'charge':None}

    namePrefix = prefix
    if tftype is 0:
        namePrefix += 'bmtf_only_'
    elif tftype is 1:
        namePrefix += 'omtf_only_'
    elif tftype is 2:
        namePrefix += 'emtf_only_'
    if emul:
        namePrefix += 'emu_'
    if legacy:
        namePrefix += 'legacy_'

    varnames = []
    binnings = {}
    varnames2d = []
    binnings2d = {}

    for var_bin in vars_bins:
        varnames.append(namePrefix+'tag_'+var_bin[0])
        binnings[namePrefix+'tag_'+var_bin[0]] = var_bin[1:]+[probe_x_title_vars[var_bin[0]]+'^{reco}', probe_x_title_units[var_bin[0]]]

    for eta_range in eta_ranges:
        eta_min = eta_range[0]
        eta_max = eta_range[1]
        eta_min_str = '_absEtaMin'+str(eta_min).replace('.', 'p')
        eta_max_str = '_absEtaMax'+str(eta_max).replace('.', 'p')
        eta_title = '{eMin} < |#eta| < {eMax}'.format(eMin=eta_min, eMax=eta_max)

        for q in range(16):
            if q in qual_ptmins_dict:
                ptmins_list = qual_ptmins_dict[q]
                qual_min_str = '_qualMin'+str(q)

                for ptmins in ptmins_list:
                    probe_pt_min = ptmins[0]
                    probe_ptmin_str = '_ptmin'+str(probe_pt_min).replace('.', 'p')
                    probe_cut_title = ' p_{T} > '+str(probe_pt_min)+' GeV/c'
                    for pt_min in ptmins[1]:
                        ptmin_str = '_ptmin'+str(pt_min).replace('.', 'p')
                        cut_title = ' p_{T} > '+str(pt_min)+' GeV/c'

                        varnames.append(namePrefix+'n_probes'+eta_min_str+eta_max_str+probe_ptmin_str)
                        binnings[namePrefix+'n_probes'+eta_min_str+eta_max_str+probe_ptmin_str] = [10, 0, 10, '# probes'+probe_cut_title]
                        for var_bin in vars_bins:
                            varnames.append(namePrefix+'l1_muon'+eta_min_str+eta_max_str+qual_min_str+ptmin_str+'_'+var_bin[0])
                            binnings[namePrefix+'l1_muon'+eta_min_str+eta_max_str+qual_min_str+ptmin_str+'_'+var_bin[0]] = var_bin[1:]+[x_title_vars[var_bin[0]]+'^{uGMT}', x_title_units[var_bin[0]]]
                        for var_bin in probe_vars_bins:
                            varnames.append(namePrefix+'probe'+eta_min_str+eta_max_str+probe_ptmin_str+'_'+var_bin[0])
                            binnings[namePrefix+'probe'+eta_min_str+eta_max_str+probe_ptmin_str+'_'+var_bin[0]] = var_bin[1:]+[probe_x_title_vars[var_bin[0]]+'^{reco}', probe_x_title_units[var_bin[0]]]

                        for delta_type in match_deltas:
                            match_delta = match_deltas[delta_type]
                            delta_str = '_'+delta_type+str(match_delta).replace('.', 'p')

                            varnames.append(namePrefix+'n_probe'+probe_ptmin_str+delta_str+'_matched_l1_muons'+eta_min_str+eta_max_str+qual_min_str+ptmin_str)
                            varnames.append(namePrefix+'n_l1_muons'+qual_min_str+ptmin_str+delta_str+'_matched_to_a_probe'+eta_min_str+eta_max_str+probe_ptmin_str)
                            for var_bin in vars_bins:
                                varnames.append(namePrefix+'best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+delta_str+'_matched_l1_muon'+qual_min_str+ptmin_str+'_'+var_bin[0])
                            for var_bin in probe_vars_bins:
                                varnames.append(namePrefix+'best_l1_muon'+qual_min_str+ptmin_str+delta_str+'_matched_probe'+eta_min_str+eta_max_str+probe_ptmin_str+'_'+var_bin[0])
                            for res_var_bin in res_vars_bins:
                                varnames.append(namePrefix+'res_best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+delta_str+'_matched_l1_muon'+qual_min_str+ptmin_str+'_'+res_var_bin[0])
                            for var_bin_2d in x_vars_bins_2d:
                                varnames2d.append(namePrefix+'2d_best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+delta_str+'_matched_l1_muon'+qual_min_str+ptmin_str+'_'+var_bin_2d[0])

                            varnames.append(namePrefix+'best_l1_muon'+qual_min_str+ptmin_str+delta_str+'_matched_probe'+eta_min_str+eta_max_str+probe_ptmin_str+'_'+delta_type)

                            # binnings
                            binnings[namePrefix+'n_probe'+probe_ptmin_str+delta_str+'_matched_l1_muons'+eta_min_str+eta_max_str+qual_min_str+ptmin_str] = [10, 0, 10, '# uGMT #mu'+cut_title+' #Delta R matched to probe'+probe_cut_title]
                            binnings[namePrefix+'n_l1_muons'+qual_min_str+ptmin_str+delta_str+'_matched_to_a_probe'+eta_min_str+eta_max_str+probe_ptmin_str+''] = [5, 0, 5, '# uGMT #mu'+cut_title+' #Delta R matched to one probe'+probe_cut_title]
                            for var_bin in vars_bins:
                                binnings[namePrefix+'best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+delta_str+'_matched_l1_muon'+qual_min_str+ptmin_str+'_'+var_bin[0]] = vars_bins[0][1:]+[x_title_vars[vars_bins[0][0]]+'^{L1}', x_title_units[vars_bins[0][0]]]
                            for var_bin in probe_vars_bins:
                                binnings[namePrefix+'probe'+eta_min_str+eta_max_str+probe_ptmin_str+'_'+var_bin[0]] = var_bin[1:]+[probe_x_title_vars[var_bin[0]]+'^{reco}', probe_x_title_units[var_bin[0]]]
                                binnings[namePrefix+'best_l1_muon'+qual_min_str+ptmin_str+delta_str+'_matched_probe'+eta_min_str+eta_max_str+probe_ptmin_str+'_'+var_bin[0]] = var_bin[1:]+[probe_x_title_vars[var_bin[0]]+'^{reco}', probe_x_title_units[var_bin[0]]]
                            for res_var_bin in res_vars_bins:
                                binnings[namePrefix+'res_best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+delta_str+'_matched_l1_muon'+qual_min_str+ptmin_str+'_'+res_var_bin[0]] = res_var_bin[1:]+[res_x_title_vars[res_var_bin[0]], res_x_title_units[res_var_bin[0]]]

                            for var_bin_2d in x_vars_bins_2d:
                                binnings2d[namePrefix+'2d_best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+delta_str+'_matched_l1_muon'+qual_min_str+ptmin_str+'_'+var_bin_2d[0]] = [var_bin_2d[1:]+[x_title_vars_2d[var_bin_2d[0]], x_title_units_2d[var_bin_2d[0]]], var_bin_2d[1:]+[y_title_vars_2d[var_bin_2d[0]], y_title_units_2d[var_bin_2d[0]]]]

                            binnings[namePrefix+'best_l1_muon'+qual_min_str+ptmin_str+delta_str+'_matched_probe'+eta_min_str+eta_max_str+probe_ptmin_str+'_'+delta_type] = [60, 0., 0.6, delta_type]

                # for resolution plots by probe pT range
                for i, probe_pt_min in enumerate(res_probe_ptmins):
                    probe_ptmin_str = '_ptmin'+str(probe_pt_min).replace('.', 'p')
                    if i < len(res_probe_ptmins)-1:
                        probe_ptmax_str = '_ptmax'+str(res_probe_ptmins[i+1]).replace('.', 'p')
                    else:
                        probe_ptmax_str = '_ptmax'

                    for delta_type in match_deltas:
                        match_delta = match_deltas[delta_type]
                        delta_str = '_'+delta_type+str(match_delta).replace('.', 'p')
                        for res_var_bin in res_vars_bins:
                            varnames.append(namePrefix+'res_best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+probe_ptmax_str+delta_str+'_matched_l1_muon'+qual_min_str+'_'+res_var_bin[0])
                            binnings[namePrefix+'res_best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+probe_ptmax_str+delta_str+'_matched_l1_muon'+qual_min_str+'_'+res_var_bin[0]] = res_var_bin[1:]+[res_x_title_vars[res_var_bin[0]], res_x_title_units[res_var_bin[0]]]

    return HistManager(list(set(varnames)), binnings), HistManager2d(list(set(varnames2d)), binnings2d)

def analyse(evt, hms, hms2d, eta_ranges, qual_ptmins_dict, res_probe_ptmins, match_deltas, emul=False, pp_run=True, legacy=False):
    recoColl = evt.recoMuon

    # at least 2 reco muons for tag and probe
    if recoColl.nMuons < 2:
        return

    # get tag muon indices
    if pp_run:
        tag_idcs = MuonSelections.select_tag_muons(recoColl, pt_min=30., abs_eta_max=2.4, pp_run=pp_run)
    else:
        tag_idcs = MuonSelections.select_tag_muons(recoColl, pt_min=18., abs_eta_max=2.4, pp_run=pp_run)
    if len(tag_idcs) < 1:
        return

    # get all probe muon indices
    all_probe_idcs = MuonSelections.select_probe_muons(recoColl, pt_min=0., pos_eta=pos_eta, neg_eta=neg_eta, pos_charge=pos_charge, neg_charge=neg_charge, extrapolated=recoExtraStation)
    #all_probe_idcs = MuonSelections.select_probe_muons(recoColl, pt_min=0., pt_max=15., pos_eta=pos_eta, neg_eta=neg_eta, pos_charge=pos_charge, neg_charge=neg_charge, extrapolated=recoExtraStation)
 
    # get all ugmt l1 muons
    bx_min = 0
    bx_max = 0
    # decide which ntuples to use
    namePrefix = prefix
    if tftype is 0:
        namePrefix += 'bmtf_only_'
    elif tftype is 1:
        namePrefix += 'omtf_only_'
    elif tftype is 2:
        namePrefix += 'emtf_only_'
    if emul:
        l1Coll = evt.upgradeEmu
        namePrefix += 'emu_'
    else:
        l1Coll = evt.upgrade

    # use translated legacy muons
    if legacy:
        l1Coll = evt.legacyGmtEmu
        namePrefix += 'legacy_'

    l1_muon_idcs = MuonSelections.select_ugmt_muons(l1Coll, pt_min=0.5, bx_min=bx_min, bx_max=bx_max, pos_eta=pos_eta, neg_eta=neg_eta, useVtxExtraCoord=useVtxExtraCoord)

    # vertex information
    nVtx = evt.recoVertex.nVtx
    # run number
    runnr = evt.event.run

    probeMomentumDict = {}

    # eta and phi variables to use depending on selected options
    zeros = [0.] * 50 # list with zeros for deta and dphi matching with the Matcher.match_dr function
    # match selected l1 muons to selected probes
    if useVtxExtraCoord:
        etas = l1Coll.muonEtaAtVtx
        phis = l1Coll.muonPhiAtVtx
    else:
        etas = l1Coll.muonEta
        phis = l1Coll.muonPhi
    # select reco muon coordinates at vertex or at 1st or 2nd muon station
    if recoExtraStation == 1:
        probeEtas = recoColl.etaSt1
        probePhis = recoColl.phiSt1
    elif recoExtraStation == 2:
        probeEtas = recoColl.etaSt2
        probePhis = recoColl.phiSt2
    else:
        probeEtas = recoColl.eta
        probePhis = recoColl.phi

    # loop over tags
    for tag_idx in tag_idcs:
        # fill tag kinematic plots
        for hm in hms:
            hm.fill(namePrefix+'tag_pt', recoColl.pt[tag_idx])
            hm.fill(namePrefix+'tag_eta', recoColl.eta[tag_idx])
            hm.fill(namePrefix+'tag_phi', recoColl.phi[tag_idx])
            hm.fill(namePrefix+'tag_charge', recoColl.charge[tag_idx])
        # remove the current tag from the list of probes
        probe_idcs = [idx for idx in all_probe_idcs if idx != tag_idx]
        # remove probes that are too close to the tag
        probe_idcs = [idx for idx in probe_idcs if Matcher.delta_r(recoColl.phi[idx], recoColl.eta[idx], recoColl.phi[tag_idx], recoColl.eta[tag_idx]) > 0.5]
        # select tag-probe pairs with invariant mass in selected window
        if useInvMassCut:
            tagLV = root.TLorentzVector()
            tagLV.SetPtEtaPhiM(recoColl.pt[tag_idx], recoColl.eta[tag_idx], recoColl.phi[tag_idx], 0.106)
        probeLV = root.TLorentzVector()
        invmass_probe_idcs = []
        for idx in probe_idcs:
            probeLV.SetPtEtaPhiM(recoColl.pt[idx], recoColl.eta[idx], recoColl.phi[idx], 0.106)
            probeMomentumDict[idx] = probeLV.P()
            if useInvMassCut:
                invMass = (tagLV + probeLV).M()
                if invMass < invMassMin or invMass > invMassMax:
                    continue
            invmass_probe_idcs.append(idx)

        # for all defined eta ranges
        for eta_range in eta_ranges:
            eta_min = eta_range[0]
            eta_max = eta_range[1]
            eta_min_str = '_absEtaMin'+str(eta_min).replace('.', 'p')
            eta_max_str = '_absEtaMax'+str(eta_max).replace('.', 'p')

            eta_probe_idcs = MuonSelections.select_reco_muons(recoColl, abs_eta_min=eta_min, abs_eta_max=eta_max, extrapolated=recoExtraStation, idcs=invmass_probe_idcs)
            eta_l1_muon_idcs = MuonSelections.select_ugmt_muons(l1Coll, abs_eta_min=eta_min, abs_eta_max=eta_max, idcs=l1_muon_idcs, useVtxExtraCoord=useVtxExtraCoord)

            # keep probe pt cuts in a list to not fill the histograms several times if two quality cuts use the same probe pt cut
            probe_pt_mins = []
            probe_ptmin_str_dict = {}
            eta_thr_probe_idcs_dict = {}
            # for all defined min quality ptmin_list combinations
            for q in range(16):
                if q in qual_ptmins_dict:
                    ptmins_list = qual_ptmins_dict[q]
                    qual_min_str = '_qualMin'+str(q)
                    eta_q_l1_muon_idcs = MuonSelections.select_ugmt_muons(l1Coll, qual_min=q, idcs=eta_l1_muon_idcs, useVtxExtraCoord=useVtxExtraCoord)

                    # for all defined min probe pt
                    for ptmins in ptmins_list:
                        if not ptmins[0] in probe_pt_mins: # fill only once for each pt min value
                            probe_pt_mins.append(ptmins[0])
                            probe_ptmin_str = '_ptmin'+str(ptmins[0]).replace('.', 'p')
                            probe_ptmin_str_dict[ptmins[0]] = probe_ptmin_str

                            eta_thr_probe_idcs = MuonSelections.select_reco_muons(recoColl, pt_min=probe_pt_mins[-1], idcs=eta_probe_idcs)
                            eta_thr_probe_idcs_dict[ptmins[0]] = eta_thr_probe_idcs
                            # fill the histograms with the probe kinematics
                            for hm in hms:
                                hm.fill(namePrefix+'n_probes'+eta_min_str+eta_max_str+probe_ptmin_str, len(eta_thr_probe_idcs))
                                for i in eta_thr_probe_idcs:
                                    hm.fill(namePrefix+'probe'+eta_min_str+eta_max_str+probe_ptmin_str+'_pt', recoColl.pt[i])
                                    hm.fill(namePrefix+'probe'+eta_min_str+eta_max_str+probe_ptmin_str+'_p', probeMomentumDict[i])
                                    hm.fill(namePrefix+'probe'+eta_min_str+eta_max_str+probe_ptmin_str+'_eta', recoColl.eta[i])
                                    hm.fill(namePrefix+'probe'+eta_min_str+eta_max_str+probe_ptmin_str+'_phi', recoColl.phi[i])
                                    hm.fill(namePrefix+'probe'+eta_min_str+eta_max_str+probe_ptmin_str+'_charge', recoColl.charge[i])
                                    hm.fill(namePrefix+'probe'+eta_min_str+eta_max_str+probe_ptmin_str+'_vtx', nVtx)
                                    hm.fill(namePrefix+'probe'+eta_min_str+eta_max_str+probe_ptmin_str+'_run', runnr)

                        probe_ptmin_str = probe_ptmin_str_dict[ptmins[0]]
                        eta_thr_probe_idcs = eta_thr_probe_idcs_dict[ptmins[0]]
                        # for all defined min l1 muon pt
                        for pt_min in ptmins[1]:
                            ptmin_str = '_ptmin'+str(pt_min).replace('.', 'p')

                            q_thr_l1_muon_idcs = MuonSelections.select_ugmt_muons(l1Coll, qual_min=q, pt_min=pt_min, idcs=l1_muon_idcs, useVtxExtraCoord=useVtxExtraCoord)
                            eta_q_thr_l1_muon_idcs = MuonSelections.select_ugmt_muons(l1Coll, pt_min=pt_min, idcs=eta_q_l1_muon_idcs, useVtxExtraCoord=useVtxExtraCoord)
                            # fill the histograms with the l1 muon kinematics
                            for hm in hms:
                                hm.fill(namePrefix+'n_probes'+eta_min_str+eta_max_str+probe_ptmin_str, len(eta_thr_probe_idcs))
                            for i in eta_q_thr_l1_muon_idcs:
                                if tftype == -1 or tftype == get_tftype(l1Coll.muonTfMuonIdx[i]):
                                    if useVtxExtraCoord:
                                        eta = l1Coll.muonEtaAtVtx[i]
                                        phi = l1Coll.muonPhiAtVtx[i]
                                    else:
                                        eta = l1Coll.muonEta[i]
                                        phi = l1Coll.muonPhi[i]

                                    for hm in hms:
                                        hm.fill(namePrefix+'l1_muon'+eta_min_str+eta_max_str+qual_min_str+ptmin_str+'_pt', l1Coll.muonEt[i])
                                        hm.fill(namePrefix+'l1_muon'+eta_min_str+eta_max_str+qual_min_str+ptmin_str+'_eta', eta)
                                        hm.fill(namePrefix+'l1_muon'+eta_min_str+eta_max_str+qual_min_str+ptmin_str+'_phi', phi)
                                        hm.fill(namePrefix+'l1_muon'+eta_min_str+eta_max_str+qual_min_str+ptmin_str+'_charge', l1Coll.muonChg[i])
                                        hm.fill(namePrefix+'l1_muon'+eta_min_str+eta_max_str+qual_min_str+ptmin_str+'_vtx', nVtx)
                                        hm.fill(namePrefix+'l1_muon'+eta_min_str+eta_max_str+qual_min_str+ptmin_str+'_run', runnr)

                            if len(q_thr_l1_muon_idcs) > 0:
                                # for all matching methodes (dr, deta, dphi)
                                for delta_type in match_deltas:
                                    match_delta = match_deltas[delta_type]
                                    delta_str = '_'+delta_type+str(match_delta).replace('.', 'p')
                                    if delta_type == 'dr':
                                        matched_l1_muons = Matcher.match_dr(etas, phis, probeEtas, probePhis, cut=match_delta, idcs1=q_thr_l1_muon_idcs, idcs2=eta_thr_probe_idcs) # match in delta R
                                    elif delta_type == 'deta':
                                        matched_l1_muons = Matcher.match_dr(etas, zeros, probeEtas, zeros, cut=match_delta, idcs1=q_thr_l1_muon_idcs, idcs2=eta_thr_probe_idcs) # match in delta eta
                                    elif delta_type == 'dphi':
                                        matched_l1_muons = Matcher.match_dr(zeros, phis, zeros, probePhis, cut=match_delta, idcs1=q_thr_l1_muon_idcs, idcs2=eta_thr_probe_idcs) # match in delta phi
                                    else:
                                        continue

                                    for hm in hms:
                                        hm.fill(namePrefix+'n_probe'+probe_ptmin_str+delta_str+'_matched_l1_muons'+eta_min_str+eta_max_str+qual_min_str+ptmin_str, len(matched_l1_muons))

                                    # how many l1 matches did we find for each probe muon
                                    for probe_idx in invmass_probe_idcs:
                                        l1_muon_cntr = 0
                                        histo_filled = False
                                        # fill matched histograms
                                        for i in range(len(matched_l1_muons)):
                                            if probe_idx == matched_l1_muons[i][1]:
                                                l1_muon_cntr += 1
                                                # fill muon values only for the first (and therefore best) match to this probe muon
                                                if not histo_filled:
                                                    if tftype == -1 or tftype == get_tftype(l1Coll.muonTfMuonIdx[matched_l1_muons[i][0]]):
                                                        eta = etas[matched_l1_muons[i][0]]
                                                        phi = phis[matched_l1_muons[i][0]]
                                                        for hm in hms:
                                                            hm.fill(namePrefix+'best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+delta_str+'_matched_l1_muon'+qual_min_str+ptmin_str+'_pt', l1Coll.muonEt[matched_l1_muons[i][0]])
                                                            hm.fill(namePrefix+'best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+delta_str+'_matched_l1_muon'+qual_min_str+ptmin_str+'_eta', eta)
                                                            hm.fill(namePrefix+'best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+delta_str+'_matched_l1_muon'+qual_min_str+ptmin_str+'_phi', phi)
                                                            hm.fill(namePrefix+'best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+delta_str+'_matched_l1_muon'+qual_min_str+ptmin_str+'_charge', l1Coll.muonChg[matched_l1_muons[i][0]])
                                                            hm.fill(namePrefix+'best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+delta_str+'_matched_l1_muon'+qual_min_str+ptmin_str+'_vtx', nVtx)
                                                            hm.fill(namePrefix+'best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+delta_str+'_matched_l1_muon'+qual_min_str+ptmin_str+'_run', runnr)
                                                            hm.fill(namePrefix+'best_l1_muon'+qual_min_str+ptmin_str+delta_str+'_matched_probe'+eta_min_str+eta_max_str+probe_ptmin_str+'_pt', recoColl.pt[probe_idx])
                                                            hm.fill(namePrefix+'best_l1_muon'+qual_min_str+ptmin_str+delta_str+'_matched_probe'+eta_min_str+eta_max_str+probe_ptmin_str+'_p', probeMomentumDict[probe_idx])
                                                            hm.fill(namePrefix+'best_l1_muon'+qual_min_str+ptmin_str+delta_str+'_matched_probe'+eta_min_str+eta_max_str+probe_ptmin_str+'_eta', recoColl.eta[probe_idx])
                                                            hm.fill(namePrefix+'best_l1_muon'+qual_min_str+ptmin_str+delta_str+'_matched_probe'+eta_min_str+eta_max_str+probe_ptmin_str+'_phi', recoColl.phi[probe_idx])
                                                            hm.fill(namePrefix+'best_l1_muon'+qual_min_str+ptmin_str+delta_str+'_matched_probe'+eta_min_str+eta_max_str+probe_ptmin_str+'_charge', recoColl.charge[probe_idx])
                                                            hm.fill(namePrefix+'best_l1_muon'+qual_min_str+ptmin_str+delta_str+'_matched_probe'+eta_min_str+eta_max_str+probe_ptmin_str+'_vtx', nVtx)
                                                            hm.fill(namePrefix+'best_l1_muon'+qual_min_str+ptmin_str+delta_str+'_matched_probe'+eta_min_str+eta_max_str+probe_ptmin_str+'_run', runnr)
                                                            hm.fill(namePrefix+'best_l1_muon'+qual_min_str+ptmin_str+delta_str+'_matched_probe'+eta_min_str+eta_max_str+probe_ptmin_str+'_'+delta_type, matched_l1_muons[i][2])
                                                            hm.fill(namePrefix+'res_best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+delta_str+'_matched_l1_muon'+qual_min_str+ptmin_str+'_dpt', l1Coll.muonEt[matched_l1_muons[i][0]] - recoColl.pt[probe_idx])
                                                            hm.fill(namePrefix+'res_best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+delta_str+'_matched_l1_muon'+qual_min_str+ptmin_str+'_dinvpt', (recoColl.pt[probe_idx] - l1Coll.muonEt[matched_l1_muons[i][0]]) / l1Coll.muonEt[matched_l1_muons[i][0]])
                                                            hm.fill(namePrefix+'res_best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+delta_str+'_matched_l1_muon'+qual_min_str+ptmin_str+'_deta', eta - probeEtas[probe_idx])
                                                            hm.fill(namePrefix+'res_best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+delta_str+'_matched_l1_muon'+qual_min_str+ptmin_str+'_dphi', phi - probePhis[probe_idx])
                                                            hm.fill(namePrefix+'res_best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+delta_str+'_matched_l1_muon'+qual_min_str+ptmin_str+'_dcharge', l1Coll.muonChg[matched_l1_muons[i][0]] - recoColl.charge[probe_idx])
                                                        for hm2d in hms2d:
                                                            hm2d.fill(namePrefix+'2d_best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+delta_str+'_matched_l1_muon'+qual_min_str+ptmin_str+'_pt', recoColl.pt[probe_idx], l1Coll.muonEt[matched_l1_muons[i][0]])
                                                            hm2d.fill(namePrefix+'2d_best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+delta_str+'_matched_l1_muon'+qual_min_str+ptmin_str+'_eta', probeEtas[probe_idx], eta)
                                                            hm2d.fill(namePrefix+'2d_best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+delta_str+'_matched_l1_muon'+qual_min_str+ptmin_str+'_phi', probePhis[probe_idx], phi)
                                                            hm2d.fill(namePrefix+'2d_best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+delta_str+'_matched_l1_muon'+qual_min_str+ptmin_str+'_charge', recoColl.charge[probe_idx], l1Coll.muonChg[matched_l1_muons[i][0]])
                                                    histo_filled = True
                                        for hm in hms:
                                            hm.fill(namePrefix+'n_l1_muons'+qual_min_str+ptmin_str+delta_str+'_matched_to_a_probe'+eta_min_str+eta_max_str+probe_ptmin_str+'', l1_muon_cntr)

                    # for resolution plots by pT range
                    for j, probe_pt_min in enumerate(res_probe_ptmins):
                        probe_ptmin_str = '_ptmin'+str(probe_pt_min).replace('.', 'p')
                        if j < len(res_probe_ptmins)-1:
                            probe_pt_max = res_probe_ptmins[j+1]
                            probe_ptmax_str = '_ptmax'+str(res_probe_ptmins[j+1]).replace('.', 'p')
                        else:
                            probe_pt_max = 1e99
                            probe_ptmax_str = '_ptmax'

                        eta_thr_probe_idcs = MuonSelections.select_reco_muons(recoColl, pt_min=probe_pt_min, pt_max=probe_pt_max, idcs=eta_probe_idcs)
                        eta_thr_probe_idcs_dict[ptmins[0]] = eta_thr_probe_idcs

                        ## fill the histograms with the probe kinematics
                        #for hm in hms:
                        #    hm.fill(namePrefix+'n_probes'+eta_min_str+eta_max_str+probe_ptmin_str+probe_ptmax_str, len(eta_thr_probe_idcs))
                        #    for i in eta_thr_probe_idcs:
                        #        hm.fill(namePrefix+'probe'+eta_min_str+eta_max_str+probe_ptmin_str+probe_ptmax_str+'_pt', recoColl.pt[i])
                        #        hm.fill(namePrefix+'probe'+eta_min_str+eta_max_str+probe_ptmin_str+probe_ptmax_str+'_p', probeMomentumDict[i])
                        #        hm.fill(namePrefix+'probe'+eta_min_str+eta_max_str+probe_ptmin_str+probe_ptmax_str+'_eta', recoColl.eta[i])
                        #        hm.fill(namePrefix+'probe'+eta_min_str+eta_max_str+probe_ptmin_str+probe_ptmax_str+'_phi', recoColl.phi[i])
                        #        hm.fill(namePrefix+'probe'+eta_min_str+eta_max_str+probe_ptmin_str+probe_ptmax_str+'_charge', recoColl.charge[i])
                        #        hm.fill(namePrefix+'probe'+eta_min_str+eta_max_str+probe_ptmin_str+probe_ptmax_str+'_vtx', nVtx)
                        #        hm.fill(namePrefix+'probe'+eta_min_str+eta_max_str+probe_ptmin_str+probe_ptmax_str+'_run', runnr)

                        q_l1_muon_idcs = MuonSelections.select_ugmt_muons(l1Coll, qual_min=q, idcs=l1_muon_idcs, useVtxExtraCoord=useVtxExtraCoord)

                        if len(q_l1_muon_idcs) > 0:
                            for delta_type in match_deltas:
                                match_delta = match_deltas[delta_type]
                                delta_str = '_'+delta_type+str(match_delta).replace('.', 'p')
                                if delta_type == 'dr':
                                    matched_l1_muons = Matcher.match_dr(etas, phis, probeEtas, probePhis, cut=match_delta, idcs1=q_l1_muon_idcs, idcs2=eta_thr_probe_idcs) # match in delta R
                                elif delta_type == 'deta':
                                    matched_l1_muons = Matcher.match_dr(etas, zeros, probeEtas, zeros, cut=match_delta, idcs1=q_l1_muon_idcs, idcs2=eta_thr_probe_idcs) # match in delta eta
                                elif delta_type == 'dphi':
                                    matched_l1_muons = Matcher.match_dr(zeros, phis, zeros, probePhis, cut=match_delta, idcs1=q_l1_muon_idcs, idcs2=eta_thr_probe_idcs) # match in delta phi
                                else:
                                    continue

                                #for hm in hms:
                                #    hm.fill(namePrefix+'n_probe'+probe_ptmin_str+probe_ptmax_str+delta_str+'_matched_l1_muons'+eta_min_str+eta_max_str+qual_min_str, len(matched_l1_muons))

                                # how many l1 matches did we find for each probe muon
                                for probe_idx in invmass_probe_idcs:
                                    l1_muon_cntr = 0
                                    histo_filled = False
                                    # fill matched histograms
                                    for i in range(len(matched_l1_muons)):
                                        if probe_idx == matched_l1_muons[i][1]:
                                            l1_muon_cntr += 1
                                            # fill muon values only for the first (and therefore best) match to this probe muon
                                            if not histo_filled:
                                                if tftype == -1 or tftype == get_tftype(l1Coll.muonTfMuonIdx[matched_l1_muons[i][0]]):
                                                    eta = etas[matched_l1_muons[i][0]]
                                                    phi = phis[matched_l1_muons[i][0]]
                                                    for hm in hms:
                                                        hm.fill(namePrefix+'res_best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+probe_ptmax_str+delta_str+'_matched_l1_muon'+qual_min_str+'_dpt', l1Coll.muonEt[matched_l1_muons[i][0]] - recoColl.pt[probe_idx])
                                                        hm.fill(namePrefix+'res_best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+probe_ptmax_str+delta_str+'_matched_l1_muon'+qual_min_str+'_dinvpt', (recoColl.pt[probe_idx] - l1Coll.muonEt[matched_l1_muons[i][0]]) / l1Coll.muonEt[matched_l1_muons[i][0]])
                                                        hm.fill(namePrefix+'res_best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+probe_ptmax_str+delta_str+'_matched_l1_muon'+qual_min_str+'_deta', eta - probeEtas[probe_idx])
                                                        hm.fill(namePrefix+'res_best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+probe_ptmax_str+delta_str+'_matched_l1_muon'+qual_min_str+'_dphi', phi - probePhis[probe_idx])
                                                        hm.fill(namePrefix+'res_best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+probe_ptmax_str+delta_str+'_matched_l1_muon'+qual_min_str+'_dcharge', l1Coll.muonChg[matched_l1_muons[i][0]] - recoColl.charge[probe_idx])

                                                histo_filled = True
                                    #for hm in hms:
                                    #    hm.fill(namePrefix+'n_l1_muons'+qual_min_str+delta_str+'_matched_to_a_probe'+eta_min_str+eta_max_str+probe_ptmin_str+probe_ptmax_str+'', l1_muon_cntr)


def save_histos(hm, hm2d, hm_runs, hm2d_runs, outfile):
    '''
    save all histograms in hm to outfile
    '''
    outfile.mkdir('all_runs')
    outfile.cd('all_runs')
    for varname in hm.get_varnames():
        hm.get(varname).Write()
    for varname in hm2d.get_varnames():
        hm2d.get(varname).Write()
    if perRunHistos:
        for runnr, hm_run in hm_runs.items():
            outfile.mkdir(str(runnr))
            outfile.cd('/'+str(runnr))
            for varname in hm_run.get_varnames():
                hm_run.get(varname).Write()
        for runnr, hm2d_run in hm2d_runs.items():
            if outfile.GetDirectory(str(runnr)) == 0:
                outfile.mkdir(str(runnr))
            outfile.cd('/'+str(runnr))
            for varname in hm2d_run.get_varnames():
                hm2d_run.get(varname).Write()
        

def main():
    L1Ana.init_l1_analysis()
    opts = parse_options_upgradeMuonHistos(parser)
    print ""

    global pos_eta
    global neg_eta
    if opts.pos_side and not opts.neg_side:
        L1Ana.log.info("Only positive eta side requested.")
        pos_eta = True
        neg_eta = False
    elif opts.neg_side and not opts.pos_side:
        L1Ana.log.info("Only negative eta side requested.")
        pos_eta = False
        neg_eta = True
    elif opts.pos_side and opts.neg_side:
        L1Ana.log.warning("Only positive side and only negative side requested. Will include both sides.")
        pos_eta = True
        neg_eta = True

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

    global useInvMassCut
    useInvMassCut = opts.invmasscut

    global useVtxExtraCoord
    useVtxExtraCoord = opts.l1extraCoord

    global recoExtraStation
    if opts.recoExtraStation == 1 or opts.recoExtraStation == 2:
        recoExtraStation = opts.recoExtraStation
    else:
        recoExtraStation = 0

    global prefix
    prefix = opts.prefix

    global tftype
    if opts.tftype == 'bmtf':
        tftype = 0
    elif opts.tftype == 'omtf':
        tftype = 1
    elif opts.tftype == 'emtf':
        tftype = 2

    global era
    era = opts.era

    emul = opts.emul
    legacy = opts.legacy
    pp_run = not opts.pa_run

    # combinations of probe_pt_min and the corresponding pt_min values for a quality
    # the first line defines which thresholds are going to be used for unmatched histograms
    if opts.ptranges == 'standard':
        ptmins_list_q12 = [[0.5, [0.5, 25]],
                           [33, [25]],
                          ]

        ptmins_list_q8 = [[0.5, [0.5, 7, 15]],
                          [10, [7]],
                          [20, [15]],
                          ]

        ptmins_list_q4 = [[0.5, [0.5, 3]],
                          [5, [3]],
                          ]
    elif opts.etaranges == 'extended':
        ptmins_list_q12 = [[0.5, [0.5, 3, 5, 7, 15, 25]],
                           [5, [3]],
                           [8, [5]],
                           [10, [7]],
                           [20, [15]],
                           [33, [25]],
                          ]

        ptmins_list_q8 = [[0.5, [0.5, 5, 7, 15, 25]],
                          [5, [3]],
                          [10, [7]],
                          [20, [15]],
                          [33, [25]],
                          ]

        ptmins_list_q4 = [[0.5, [0.5, 3, 5, 25]],
                          [5, [3]],
                          [8, [5]],
                          [33, [25]],
                          ]
    else:
        ptmins_list_q12 = [[0.5, [0.5]],
                          ]

        ptmins_list_q8 = [[0.5, [0.5]],
                          ]

        ptmins_list_q4 = [[0.5, [0.5]],
                          ]

    res_probe_ptmins = [0.5, 20, 30, 40, 50, 60, 100, 150]

    if opts.etaranges == 'standard':
        eta_ranges = [[0, 2.4], [0, 0.83], [0.83, 1.24], [1.24, 2.4]]
    elif opts.etaranges == 'extended':
        eta_ranges = [[0, 2.4], [0, 0.83], [0.83, 1.24], [1.24, 2.4], [1.24, 1.55], [1.55, 1.85], [1.85, 2.4]]
    else:
        eta_ranges = [[0, 2.4]]

    if legacy:
        qual_ptmins_dict = {5:ptmins_list_q12, 3:ptmins_list_q8, 2:ptmins_list_q4}
    else:
        qual_ptmins_dict = {12:ptmins_list_q12, 8:ptmins_list_q8, 4:ptmins_list_q4}
    match_deltas = {'dr':0.5, 'deta':0.5, 'dphi':0.5} # max deltas for matching

    # book the histograms
    L1Ana.log.info("Booking combined run histograms.")
    hm, hm2d = book_histograms(eta_ranges, qual_ptmins_dict, res_probe_ptmins, match_deltas, emul=emul, legacy=legacy)
    # histogram dicts per run
    hm_runs = {}
    hm2d_runs = {}

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

            # book histograms for this event's run number if not already done
            if perRunHistos and not runnr in hm_runs:
                L1Ana.log.info("Booking histograms for run {r}.".format(r=runnr))
                hm_runs[runnr], hm2d_runs[runnr] = book_histograms(eta_ranges, qual_ptmins_dict, res_probe_ptmins, match_deltas, emul=emul, legacy=legacy)

            # now do the analysis for all pt cut combinations
            if perRunHistos:
                analyse(event, [hm, hm_runs[runnr]], [hm2d, hm2d_runs[runnr]], eta_ranges, qual_ptmins_dict, res_probe_ptmins, match_deltas, emul=emul, pp_run=pp_run, legacy=legacy)
            else:
                analyse(event, [hm], [hm2d], eta_ranges, qual_ptmins_dict, res_probe_ptmins, match_deltas, emul=emul, pp_run=pp_run, legacy=legacy)
            analysed_evt_ctr += 1
    except KeyboardInterrupt:
        L1Ana.log.info("Analysis interrupted after {n} events".format(n=i))

    L1Ana.log.info("Analysis of {nAna} events in selected runs/LS finished.".format(nAna=analysed_evt_ctr))

    # save histos to root file
    if saveHistos:
        output = root.TFile(opts.outname, 'recreate')
        output.cd()
        save_histos(hm, hm2d, hm_runs, hm2d_runs, output)
        output.Close()

if __name__ == "__main__":
    pos_eta = True
    neg_eta = True
    pos_charge = True
    neg_charge = True
    useInvMassCut = False
    invMassMin = 71
    invMassMax = 111
    useVtxExtraCoord = False
    recoExtraStation = 0
    prefix = ''
    tftype = -1
    era = ''
    saveHistos = True
    best_only = False
    perRunHistos = False
    main()

