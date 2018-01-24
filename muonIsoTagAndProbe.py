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
    sub_parser = parsers.add_parser("muonIsoTagAndProbe")
    sub_parser.add_argument("-o", "--outname", dest="outname", default="./ugmt_tandp_eff_histos.root", type=str, help="A root file name where to save the histograms.")
    sub_parser.add_argument("-j", "--json", dest="json", type=str, default=None, help="A json file with good lumi sections per run.")
    sub_parser.add_argument("-r", "--runs", dest="runs", type=str, default=None, help="A string of runs to check.")
    sub_parser.add_argument("-p", "--pos-side", dest="pos_side", default=False, action="store_true", help="Positive detector side only.")
    sub_parser.add_argument("-n", "--neg-side", dest="neg_side", default=False, action="store_true", help="Negative detector side only.")
    sub_parser.add_argument("--pos-charge", dest="pos_charge", default=False, action="store_true", help="Positive probe charge only.")
    sub_parser.add_argument("--neg-charge", dest="neg_charge", default=False, action="store_true", help="Negative probe charge only.")
    sub_parser.add_argument("--use-inv-mass-cut", dest="invmasscut", default=False, action="store_true", help="Use an invariant mass range for the tag and probe pair.")
    sub_parser.add_argument("--use-extra-coord", dest="extraCoord", default=False, action="store_true", help="Use L1 extrapolated eta and phi coordinates.")
    sub_parser.add_argument("--eta-restricted", dest="etarestricted", type=float, default=3., help="Upper eta value for isolation.")
    sub_parser.add_argument("--emul", dest="emul", default=False, action="store_true", help="Make emulator plots.")
    sub_parser.add_argument("--prefix", dest="prefix", type=str, default='', help="A prefix for the histogram names.")
    sub_parser.add_argument("--tftype", dest="tftype", type=str, default='', help="Fill L1 muons from one TF.")
    sub_parser.add_argument("--iso-method", dest="isomethod", type=str, default='abs', help="Isolation method. ['abs', 'rel', 'inner', 'outovertot, 'inner2x2', 'outovertot2x2', 'mipptadjust']")
    sub_parser.add_argument("--nvtx-min", dest="nvtxmin", type=int, default=None, help="Minimum number of vertices.")
    sub_parser.add_argument("--nvtx-max", dest="nvtxmax", type=int, default=None, help="Maximum number of vertices.")

    opts, unknown = parser.parse_known_args()
    return opts

def get_tftype(tf_muon_index):
    if tf_muon_index > 35 and tf_muon_index < 72:
        return 0 # BMTF
    elif tf_muon_index > 17 and tf_muon_index < 90:
        return 1 # OMTF
    else:
        return 2 # EMTF

def book_histograms(eta_ranges, qual_ptmins_dict, match_deltas, iso_wps, emul=False):
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

    p_bins = pt_bins

    dr_str = '_dr'+str(match_deltas['dr'])

    vars_bins = [['pass', 1, 1, 2], ['pt', -1]+pt_bins, ['eta', 100, -2.5, 2.5], ['phi', 70, -3.5, 3.5], ['charge', 3, -1, 2], ['vtx', 100, 0, 100], ['run', 17000, 271725, 288725]]
    x_title_vars = {'pass':'pass', 'pt':'p_{T}', 'eta':'#eta', 'phi':'#phi', 'charge':'charge', 'vtx':'PU', 'run':'run number'}
    x_title_units = {'pass':None, 'pt':'GeV/c', 'eta':None, 'phi':None, 'charge':None, 'vtx':None, 'run':None}
    probe_vars_bins = vars_bins + [['p', -1]+p_bins]
    probe_x_title_vars = {'p':'p'}
    probe_x_title_vars.update(x_title_vars)
    probe_x_title_units = {'p':'GeV/c'}
    probe_x_title_units.update(x_title_units)
    res_vars_bins = [['dpt', 100, -50, 50], ['dinvpt', 200, -0.2, 0.2], ['deta', 100, -0.1, 0.1], ['dphi', 100, -0.2, 0.2]]
    res_x_title_vars = {'dpt':'p_{T}^{reco} - p_{T}^{L1}', 'dinvpt':'1/p_{T}^{reco} - 1/p_{T}^{L1}', 'deta':'#eta_{reco} - #eta_{L1}', 'dphi':'#phi_{reco} - #phi_{L1}'}
    res_x_title_units = {'dpt':'GeV', 'dinvpt':'GeV', 'deta':None, 'dphi':None}
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

    varnames = []
    binnings = {}
    varnames2d = []
    binnings2d = {}

    for iso_wp in iso_wps:
        iso_wp_str = '_isoMax{iso:.3f}'.format(iso=iso_wp)
        for eta_range in eta_ranges:
            eta_min = eta_range[0]
            eta_max = eta_range[1]
            eta_min_str = '_absEtaMin'+str(eta_min)
            eta_max_str = '_absEtaMax'+str(eta_max)
            eta_title = '{eMin} < |#eta| < {eMax}'.format(eMin=eta_min, eMax=eta_max)

            for q in range(16):
                if q in qual_ptmins_dict:
                    ptmins_list = qual_ptmins_dict[q]
                    qual_min_str = '_qualMin'+str(q)

                    for ptmins in ptmins_list:
                        probe_pt_min = ptmins[0]
                        probe_ptmin_str = '_ptmin'+str(probe_pt_min)
                        probe_cut_title = ' p_{T} > '+str(probe_pt_min)+' GeV/c'
                        for pt_min in ptmins[1]:
                            ptmin_str = '_ptmin'+str(pt_min)
                            cut_title = ' p_{T} > '+str(pt_min)+' GeV/c'

                            varnames.append(namePrefix+'n_probes'+eta_min_str+eta_max_str+probe_ptmin_str)
                            varnames.append(namePrefix+'n_probe'+probe_ptmin_str+dr_str+'_matched_l1_muons'+eta_min_str+eta_max_str+qual_min_str+ptmin_str+iso_wp_str)
                            varnames.append(namePrefix+'n_l1_muons'+qual_min_str+ptmin_str+iso_wp_str+dr_str+'_matched_to_a_probe'+eta_min_str+eta_max_str+probe_ptmin_str)
                            for var_bin in vars_bins:
                                varnames.append(namePrefix+'l1_muon'+eta_min_str+eta_max_str+qual_min_str+ptmin_str+iso_wp_str+'.'+var_bin[0])
                                #varnames.append(namePrefix+'best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+dr_str+'_matched_l1_muon'+qual_min_str+ptmin_str+iso_wp_str+'.'+var_bin[0])
                            for var_bin in probe_vars_bins:
                                varnames.append(namePrefix+'probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.'+var_bin[0])
                                varnames.append(namePrefix+'best_l1_muon'+qual_min_str+ptmin_str+iso_wp_str+dr_str+'_matched_probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.'+var_bin[0])
                            #for res_var_bin in res_vars_bins:
                            #    varnames.append(namePrefix+'res_best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+dr_str+'_matched_l1_muon'+qual_min_str+ptmin_str+iso_wp_str+'.'+res_var_bin[0])
                            #for var_bin_2d in x_vars_bins_2d:
                            #    varnames2d.append(namePrefix+'2d_best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+dr_str+'_matched_l1_muon'+qual_min_str+ptmin_str+iso_wp_str+'.'+var_bin_2d[0])

                            varnames.append(namePrefix+'best_l1_muon'+qual_min_str+ptmin_str+iso_wp_str+dr_str+'_matched_probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.dr')

                            # binnings
                            binnings[namePrefix+'n_probes'+eta_min_str+eta_max_str+probe_ptmin_str] = [10, 0, 10, '# probes'+probe_cut_title]
                            binnings[namePrefix+'n_probe'+probe_ptmin_str+dr_str+'_matched_l1_muons'+eta_min_str+eta_max_str+qual_min_str+ptmin_str+iso_wp_str] = [10, 0, 10, '# uGMT #mu'+cut_title+' #Delta R matched to probe'+probe_cut_title]
                            binnings[namePrefix+'n_l1_muons'+qual_min_str+ptmin_str+iso_wp_str+dr_str+'_matched_to_a_probe'+eta_min_str+eta_max_str+probe_ptmin_str+''] = [5, 0, 5, '# uGMT #mu'+cut_title+' #Delta R matched to one probe'+probe_cut_title]
                            for var_bin in vars_bins:
                                binnings[namePrefix+'l1_muon'+eta_min_str+eta_max_str+qual_min_str+ptmin_str+iso_wp_str+'.'+var_bin[0]] = var_bin[1:]+[x_title_vars[var_bin[0]]+'^{uGMT}', x_title_units[var_bin[0]]]
                                #binnings[namePrefix+'best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+dr_str+'_matched_l1_muon'+qual_min_str+ptmin_str+iso_wp_str+'.'+var_bin[0]] = vars_bins[0][1:]+[x_title_vars[vars_bins[0][0]]+'^{L1}', x_title_units[vars_bins[0][0]]]
                            for var_bin in probe_vars_bins:
                                binnings[namePrefix+'probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.'+var_bin[0]] = var_bin[1:]+[probe_x_title_vars[var_bin[0]]+'^{reco}', probe_x_title_units[var_bin[0]]]
                                binnings[namePrefix+'best_l1_muon'+qual_min_str+ptmin_str+iso_wp_str+dr_str+'_matched_probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.'+var_bin[0]] = var_bin[1:]+[probe_x_title_vars[var_bin[0]]+'^{reco}', probe_x_title_units[var_bin[0]]]
                            #for res_var_bin in res_vars_bins:
                            #    binnings[namePrefix+'res_best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+dr_str+'_matched_l1_muon'+qual_min_str+ptmin_str+iso_wp_str+'.'+res_var_bin[0]] = res_var_bin[1:]+[res_x_title_vars[res_var_bin[0]], res_x_title_units[res_var_bin[0]]]
                            #for var_bin_2d in x_vars_bins_2d:
                            #    binnings2d[namePrefix+'2d_best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+dr_str+'_matched_l1_muon'+qual_min_str+ptmin_str+iso_wp_str+'.'+var_bin_2d[0]] = [var_bin_2d[1:]+[x_title_vars_2d[var_bin_2d[0]], x_title_units_2d[var_bin_2d[0]]], var_bin_2d[1:]+[y_title_vars_2d[var_bin_2d[0]], y_title_units_2d[var_bin_2d[0]]]]

                            binnings[namePrefix+'best_l1_muon'+qual_min_str+ptmin_str+iso_wp_str+dr_str+'_matched_probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.dr'] = [60, 0., 0.6, '#Delta R']

    return HistManager(list(set(varnames)), binnings), HistManager2d(list(set(varnames2d)), binnings2d)

def fill_matched_muons(evt, hm, matched_muons, muon_type='', eta_strs = ['', ''], ptmin_strs=['', '']):
    eta_min_str = eta_strs[0]
    eta_max_str = eta_strs[1]
    probe_ptmin_str = ptmin_strs[0]
    ptmin_str = ptmin_strs[1]

    recoColl = evt.recoMuon
    quals = evt.upgrade.muonQual
    if muon_type == 'u': # uGMT
        muon_str = 'l1_muon'

    for matched_muon in matched_muons:
        # fill histograms with matched reco muon data
        hm.fill(muon_str+qual_min_str+ptmin_str+dr_str+'_matched_probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.pt', recoColl.pt[matched_muon[1]])
        hm.fill(muon_str+qual_min_str+ptmin_str+dr_str+'_matched_probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.eta', recoColl.eta[matched_muon[1]])
        hm.fill(muon_str+qual_min_str+ptmin_str+dr_str+'_matched_probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.phi', recoColl.phi[matched_muon[1]])
        hm.fill(muon_str+qual_min_str+ptmin_str+dr_str+'_matched_probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.dr', matched_muon[2])

def analyse(evt, hm, hm2d, hm_run, hm2d_run, eta_ranges, qual_ptmins_dict, match_deltas, iso_wps, emul=False):
    recoColl = evt.recoMuon

    # at least 2 reco muons for tag and probe
    if recoColl.nMuons < 2:
        return

    # get tag muon indices
    tag_idcs = MuonSelections.select_tag_muons(recoColl, pt_min=27., abs_eta_max=2.4)
    if len(tag_idcs) < 1:
        return

    # get all probe muon indices
    all_probe_idcs = MuonSelections.select_probe_muons(recoColl, pt_min=0., pos_eta=pos_eta, neg_eta=neg_eta, pos_charge=pos_charge, neg_charge=neg_charge)
    #all_probe_idcs = MuonSelections.select_probe_muons(recoColl, pt_min=0., pt_max=15., pos_eta=pos_eta, neg_eta=neg_eta, pos_charge=pos_charge, neg_charge=neg_charge)
 
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
        l1CaloTowerColl = evt.caloTowersEmu
        namePrefix += 'emu_'
    else:
        l1Coll = evt.upgrade
        l1CaloTowerColl = evt.caloTowers
    l1_muon_idcs = MuonSelections.select_ugmt_muons(l1Coll, pt_min=0.5, bx_min=bx_min, bx_max=bx_max, pos_eta=pos_eta, neg_eta=neg_eta, useVtxExtraCoord=useVtxExtraCoord)

    # vertex information
    nVtx = evt.recoVertex.nVtx
    # run number
    runnr = evt.event.run

    matchdr = match_deltas['dr']
    dr_str = '_dr'+str(matchdr)

    probeMomentumDict = {}

    # loop over tags
    for tag_idx in tag_idcs:
        # TODO fill tag kinematic plots
        # remove the current tag from the list of probes
        probe_idcs = [idx for idx in all_probe_idcs if idx != tag_idx]
        # remove probes that are too close to the tag
        probe_idcs = [idx for idx in probe_idcs if Matcher.delta_r(recoColl.phi[idx], recoColl.eta[idx], recoColl.phi[tag_idx] ,recoColl.eta[tag_idx]) > 0.5]
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

        probesFilled = False
        for iso_wp in iso_wps:
            iso_wp_str = '_isoMax{iso:.3f}'.format(iso=iso_wp)
            # remove non-isolated muons
            if iso_type == 2 or iso_type == 4:
                iso_min=iso_wp
                iso_max=999
            else:
                iso_min=0.
                iso_max=iso_wp
            if iso_type == 6:
                l1_iso_muon_idcs = l1_muon_idcs 
                l1_iso_muon_idcs_corr = MuonSelections.select_iso_ugmt_muons(l1Coll, l1CaloTowerColl, iso_min=iso_min, iso_max=iso_max, iso_eta_max=isoEtaMax, idcs=l1_muon_idcs, useVtxExtraCoord=useVtxExtraCoord, iso_type=iso_type)
            elif iso_type == 7:
                iso_min=iso_wp
                iso_max=1.
                l1_iso_muon_idcs = l1_muon_idcs 
                l1_iso_muon_idcs_corr = MuonSelections.select_iso_ugmt_muons(l1Coll, l1CaloTowerColl, iso_min=iso_min, iso_max=iso_max, iso_eta_max=isoEtaMax, idcs=l1_muon_idcs, useVtxExtraCoord=useVtxExtraCoord, iso_type=iso_type)
            else:
                l1_iso_muon_idcs = MuonSelections.select_iso_ugmt_muons(l1Coll, l1CaloTowerColl, iso_min=iso_min, iso_max=iso_max, iso_eta_max=isoEtaMax, idcs=l1_muon_idcs, useVtxExtraCoord=useVtxExtraCoord, iso_type=iso_type)

            # for all defined eta ranges
            for eta_range in eta_ranges:
                eta_min = eta_range[0]
                eta_max = eta_range[1]
                eta_min_str = '_absEtaMin'+str(eta_min)
                eta_max_str = '_absEtaMax'+str(eta_max)

                eta_probe_idcs = MuonSelections.select_reco_muons(recoColl, abs_eta_min=eta_min, abs_eta_max=eta_max, idcs=invmass_probe_idcs)
                eta_l1_iso_muon_idcs = MuonSelections.select_ugmt_muons(l1Coll, abs_eta_min=eta_min, abs_eta_max=eta_max, idcs=l1_iso_muon_idcs, useVtxExtraCoord=useVtxExtraCoord)

                # keep probe pt cuts in a list to not fill the histograms several times if two quality cuts use the same probe pt cut
                probe_pt_mins = []
                probe_ptmin_str_dict = {}
                eta_thr_probe_idcs_dict = {}
                # for all defined min quality ptmin_list combinations
                for q in range(16):
                    if q in qual_ptmins_dict:
                        ptmins_list = qual_ptmins_dict[q]
                        qual_min_str = '_qualMin'+str(q)
                        eta_q_l1_iso_muon_idcs = MuonSelections.select_ugmt_muons(l1Coll, qual_min=q, idcs=eta_l1_iso_muon_idcs, useVtxExtraCoord=useVtxExtraCoord)

                        # for all defined min probe pt
                        for ptmins in ptmins_list:
                            if not ptmins[0] in probe_pt_mins: # fill only once for each pt min value
                                probe_pt_mins.append(ptmins[0])
                                probe_ptmin_str = '_ptmin'+str(ptmins[0])
                                probe_ptmin_str_dict[ptmins[0]] = probe_ptmin_str

                                eta_thr_probe_idcs = MuonSelections.select_reco_muons(recoColl, pt_min=probe_pt_mins[-1], idcs=eta_probe_idcs)
                                eta_thr_probe_idcs_dict[ptmins[0]] = eta_thr_probe_idcs
                                # fill the histograms with the probe kinematics if not already done for previous iso threshold
                                if not probesFilled:
                                    hm.fill(namePrefix+'n_probes'+eta_min_str+eta_max_str+probe_ptmin_str, len(eta_thr_probe_idcs))
                                    if perRunHistos:
                                        hm_run.fill(namePrefix+'n_probes'+eta_min_str+eta_max_str+probe_ptmin_str, len(eta_thr_probe_idcs))
                                    for i in eta_thr_probe_idcs:
                                        hm.fill(namePrefix+'probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.pass', 1)
                                        hm.fill(namePrefix+'probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.pt', recoColl.pt[i])
                                        hm.fill(namePrefix+'probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.p', probeMomentumDict[i])
                                        hm.fill(namePrefix+'probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.eta', recoColl.eta[i])
                                        hm.fill(namePrefix+'probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.phi', recoColl.phi[i])
                                        hm.fill(namePrefix+'probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.charge', recoColl.charge[i])
                                        hm.fill(namePrefix+'probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.vtx', nVtx)
                                        hm.fill(namePrefix+'probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.run', runnr)
                                        if perRunHistos:
                                            hm_run.fill(namePrefix+'probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.pass', 1)
                                            hm_run.fill(namePrefix+'probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.pt', recoColl.pt[i])
                                            hm_run.fill(namePrefix+'probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.p', probeMomentumDict[i])
                                            hm_run.fill(namePrefix+'probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.eta', recoColl.eta[i])
                                            hm_run.fill(namePrefix+'probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.phi', recoColl.phi[i])
                                            hm_run.fill(namePrefix+'probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.charge', recoColl.charge[i])
                                            hm_run.fill(namePrefix+'probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.vtx', nVtx)
                                            hm_run.fill(namePrefix+'probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.run', runnr)

                            probe_ptmin_str = probe_ptmin_str_dict[ptmins[0]]
                            eta_thr_probe_idcs = eta_thr_probe_idcs_dict[ptmins[0]]
                            # for all defined min l1 muon pt
                            for pt_min in ptmins[1]:
                                ptmin_str = '_ptmin'+str(pt_min)

                                q_thr_l1_iso_muon_idcs = MuonSelections.select_ugmt_muons(l1Coll, qual_min=q, pt_min=pt_min, idcs=l1_iso_muon_idcs, corr_idcs=l1_iso_muon_idcs_corr, pt_corr=ptCorrVal, useVtxExtraCoord=useVtxExtraCoord)
                                eta_q_thr_l1_iso_muon_idcs = MuonSelections.select_ugmt_muons(l1Coll, pt_min=pt_min, idcs=eta_q_l1_iso_muon_idcs, corr_idcs=l1_iso_muon_idcs_corr, pt_corr=ptCorrVal, useVtxExtraCoord=useVtxExtraCoord)
                                # fill the histograms with the l1 muon kinematics
                                for i in eta_q_thr_l1_iso_muon_idcs:
                                    if tftype == -1 or tftype == get_tftype(l1Coll.muonTfMuonIdx[i]):
                                        if useVtxExtraCoord:
                                            eta = l1Coll.muonEtaAtVtx[i]
                                            phi = l1Coll.muonPhiAtVtx[i]
                                        else:
                                            eta = l1Coll.muonEta[i]
                                            phi = l1Coll.muonPhi[i]

                                        ptCorr = 1.
                                        if (iso_type == 6 or iso_type == 7) and i in l1_iso_muon_idcs_corr:
                                            ptCorr = ptCorrVal

                                        hm.fill(namePrefix+'l1_muon'+eta_min_str+eta_max_str+qual_min_str+ptmin_str+iso_wp_str+'.pass', 1)
                                        hm.fill(namePrefix+'l1_muon'+eta_min_str+eta_max_str+qual_min_str+ptmin_str+iso_wp_str+'.pt', l1Coll.muonEt[i] * ptCorr)
                                        hm.fill(namePrefix+'l1_muon'+eta_min_str+eta_max_str+qual_min_str+ptmin_str+iso_wp_str+'.eta', eta)
                                        hm.fill(namePrefix+'l1_muon'+eta_min_str+eta_max_str+qual_min_str+ptmin_str+iso_wp_str+'.phi', phi)
                                        hm.fill(namePrefix+'l1_muon'+eta_min_str+eta_max_str+qual_min_str+ptmin_str+iso_wp_str+'.charge', l1Coll.muonChg[i])
                                        hm.fill(namePrefix+'l1_muon'+eta_min_str+eta_max_str+qual_min_str+ptmin_str+iso_wp_str+'.vtx', nVtx)
                                        hm.fill(namePrefix+'l1_muon'+eta_min_str+eta_max_str+qual_min_str+ptmin_str+iso_wp_str+'.run', runnr)
                                        if perRunHistos:
                                            hm_run.fill(namePrefix+'l1_muon'+eta_min_str+eta_max_str+qual_min_str+ptmin_str+iso_wp_str+'.pass', 1)
                                            hm_run.fill(namePrefix+'l1_muon'+eta_min_str+eta_max_str+qual_min_str+ptmin_str+iso_wp_str+'.pt', l1Coll.muonEt[i] * ptCorr)
                                            hm_run.fill(namePrefix+'l1_muon'+eta_min_str+eta_max_str+qual_min_str+ptmin_str+iso_wp_str+'.eta', eta)
                                            hm_run.fill(namePrefix+'l1_muon'+eta_min_str+eta_max_str+qual_min_str+ptmin_str+iso_wp_str+'.phi', phi)
                                            hm_run.fill(namePrefix+'l1_muon'+eta_min_str+eta_max_str+qual_min_str+ptmin_str+iso_wp_str+'.charge', l1Coll.muonChg[i])
                                            hm_run.fill(namePrefix+'l1_muon'+eta_min_str+eta_max_str+qual_min_str+ptmin_str+iso_wp_str+'.vtx', nVtx)
                                            hm_run.fill(namePrefix+'l1_muon'+eta_min_str+eta_max_str+qual_min_str+ptmin_str+iso_wp_str+'.run', runnr)


                                if len(q_thr_l1_iso_muon_idcs) > 0:
                                    zeros = [0.] * 50 # list with zeros for deta and dphi matching with the Matcher.match_dr function
                                    # match selected l1 muons to selected probes
                                    if useVtxExtraCoord:
                                        etas = l1Coll.muonEtaAtVtx
                                        phis = l1Coll.muonPhiAtVtx
                                    else:
                                        etas = l1Coll.muonEta
                                        phis = l1Coll.muonPhi
                                    matched_l1_muons = Matcher.match_dr(etas, phis, recoColl.eta, recoColl.phi, cut=matchdr, idcs1=q_thr_l1_iso_muon_idcs, idcs2=eta_thr_probe_idcs) # match in delta R
                                    hm.fill(namePrefix+'n_probe'+probe_ptmin_str+dr_str+'_matched_l1_muons'+eta_min_str+eta_max_str+qual_min_str+ptmin_str+iso_wp_str, len(matched_l1_muons))
                                    if perRunHistos:
                                        hm_run.fill(namePrefix+'n_probe'+probe_ptmin_str+dr_str+'_matched_l1_muons'+eta_min_str+eta_max_str+qual_min_str+ptmin_str+iso_wp_str, len(matched_l1_muons))

                                    # how many l1 matches did we find for each probe muon
                                    for probe_idx in invmass_probe_idcs:
                                        l1_muon_cntr = 0
                                        histo_filled = False
                                        # fill dR matched histograms
                                        for i in range(len(matched_l1_muons)):
                                            if probe_idx == matched_l1_muons[i][1]:
                                                l1_muon_cntr += 1
                                                # fill muon values only for the first (and therefore best) match to this probe muon
                                                if not histo_filled:
                                                    if tftype == -1 or tftype == get_tftype(l1Coll.muonTfMuonIdx[matched_l1_muons[i][0]]):
                                                        if useVtxExtraCoord:
                                                            eta = l1Coll.muonEtaAtVtx[matched_l1_muons[i][0]]
                                                            phi = l1Coll.muonPhiAtVtx[matched_l1_muons[i][0]]
                                                        else:
                                                            eta = l1Coll.muonEta[matched_l1_muons[i][0]]
                                                            phi = l1Coll.muonPhi[matched_l1_muons[i][0]]
                                                        #hm.fill(namePrefix+'best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+dr_str+'_matched_l1_muon'+qual_min_str+ptmin_str+iso_wp_str+'.pass', 1)
                                                        #hm.fill(namePrefix+'best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+dr_str+'_matched_l1_muon'+qual_min_str+ptmin_str+iso_wp_str+'.pt', l1Coll.muonEt[matched_l1_muons[i][0]] * ptCorr)
                                                        #hm.fill(namePrefix+'best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+dr_str+'_matched_l1_muon'+qual_min_str+ptmin_str+iso_wp_str+'.eta', eta)
                                                        #hm.fill(namePrefix+'best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+dr_str+'_matched_l1_muon'+qual_min_str+ptmin_str+iso_wp_str+'.phi', phi)
                                                        #hm.fill(namePrefix+'best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+dr_str+'_matched_l1_muon'+qual_min_str+ptmin_str+iso_wp_str+'.charge', l1Coll.muonChg[matched_l1_muons[i][0]])
                                                        #hm.fill(namePrefix+'best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+dr_str+'_matched_l1_muon'+qual_min_str+ptmin_str+iso_wp_str+'.vtx', nVtx)
                                                        #hm.fill(namePrefix+'best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+dr_str+'_matched_l1_muon'+qual_min_str+ptmin_str+iso_wp_str+'.run', runnr)
                                                        hm.fill(namePrefix+'best_l1_muon'+qual_min_str+ptmin_str+iso_wp_str+dr_str+'_matched_probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.pass', 1)
                                                        hm.fill(namePrefix+'best_l1_muon'+qual_min_str+ptmin_str+iso_wp_str+dr_str+'_matched_probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.pt', recoColl.pt[probe_idx])
                                                        hm.fill(namePrefix+'best_l1_muon'+qual_min_str+ptmin_str+iso_wp_str+dr_str+'_matched_probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.p', probeMomentumDict[probe_idx])
                                                        hm.fill(namePrefix+'best_l1_muon'+qual_min_str+ptmin_str+iso_wp_str+dr_str+'_matched_probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.eta', recoColl.eta[probe_idx])
                                                        hm.fill(namePrefix+'best_l1_muon'+qual_min_str+ptmin_str+iso_wp_str+dr_str+'_matched_probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.phi', recoColl.phi[probe_idx])
                                                        hm.fill(namePrefix+'best_l1_muon'+qual_min_str+ptmin_str+iso_wp_str+dr_str+'_matched_probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.charge', recoColl.charge[probe_idx])
                                                        hm.fill(namePrefix+'best_l1_muon'+qual_min_str+ptmin_str+iso_wp_str+dr_str+'_matched_probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.vtx', nVtx)
                                                        hm.fill(namePrefix+'best_l1_muon'+qual_min_str+ptmin_str+iso_wp_str+dr_str+'_matched_probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.run', runnr)
                                                        hm.fill(namePrefix+'best_l1_muon'+qual_min_str+ptmin_str+iso_wp_str+dr_str+'_matched_probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.dr', matched_l1_muons[i][2])
                                                        #hm.fill(namePrefix+'res_best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+dr_str+'_matched_l1_muon'+qual_min_str+ptmin_str+iso_wp_str+'.dpt', recoColl.pt[probe_idx] - (l1Coll.muonEt[matched_l1_muons[i][0]] * ptCorr))
                                                        #hm.fill(namePrefix+'res_best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+dr_str+'_matched_l1_muon'+qual_min_str+ptmin_str+iso_wp_str+'.dinvpt', 1./recoColl.pt[probe_idx] - 1./l1Coll.muonEt[matched_l1_muons[i][0]])
                                                        #hm.fill(namePrefix+'res_best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+dr_str+'_matched_l1_muon'+qual_min_str+ptmin_str+iso_wp_str+'.deta', recoColl.eta[probe_idx] - eta)
                                                        #hm.fill(namePrefix+'res_best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+dr_str+'_matched_l1_muon'+qual_min_str+ptmin_str+iso_wp_str+'.dphi', recoColl.phi[probe_idx] - phi)
                                                        #hm2d.fill(namePrefix+'2d_best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+dr_str+'_matched_l1_muon'+qual_min_str+ptmin_str+iso_wp_str+'.pt', recoColl.pt[probe_idx], l1Coll.muonEt[matched_l1_muons[i][0]] * ptCorr)
                                                        #hm2d.fill(namePrefix+'2d_best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+dr_str+'_matched_l1_muon'+qual_min_str+ptmin_str+iso_wp_str+'.eta', recoColl.eta[probe_idx], eta)
                                                        #hm2d.fill(namePrefix+'2d_best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+dr_str+'_matched_l1_muon'+qual_min_str+ptmin_str+iso_wp_str+'.phi', recoColl.phi[probe_idx], phi)
                                                        #hm2d.fill(namePrefix+'2d_best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+dr_str+'_matched_l1_muon'+qual_min_str+ptmin_str+iso_wp_str+'.charge', recoColl.charge[probe_idx], l1Coll.muonChg[matched_l1_muons[i][0]])

                                                        if perRunHistos:
                                                            #hm_run.fill(namePrefix+'best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+dr_str+'_matched_l1_muon'+qual_min_str+ptmin_str+iso_wp_str+'.pass', 1)
                                                            #hm_run.fill(namePrefix+'best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+dr_str+'_matched_l1_muon'+qual_min_str+ptmin_str+iso_wp_str+'.pt', l1Coll.muonEt[matched_l1_muons[i][0]] * ptCorr)
                                                            #hm_run.fill(namePrefix+'best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+dr_str+'_matched_l1_muon'+qual_min_str+ptmin_str+iso_wp_str+'.eta', eta)
                                                            #hm_run.fill(namePrefix+'best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+dr_str+'_matched_l1_muon'+qual_min_str+ptmin_str+iso_wp_str+'.phi', phi)
                                                            #hm_run.fill(namePrefix+'best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+dr_str+'_matched_l1_muon'+qual_min_str+ptmin_str+iso_wp_str+'.charge', l1Coll.muonChg[matched_l1_muons[i][0]])
                                                            #hm_run.fill(namePrefix+'best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+dr_str+'_matched_l1_muon'+qual_min_str+ptmin_str+iso_wp_str+'.vtx', nVtx)
                                                            #hm_run.fill(namePrefix+'best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+dr_str+'_matched_l1_muon'+qual_min_str+ptmin_str+iso_wp_str+'.run', runnr)
                                                            hm_run.fill(namePrefix+'best_l1_muon'+qual_min_str+ptmin_str+iso_wp_str+dr_str+'_matched_probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.pass', 1)
                                                            hm_run.fill(namePrefix+'best_l1_muon'+qual_min_str+ptmin_str+iso_wp_str+dr_str+'_matched_probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.pt', recoColl.pt[probe_idx])
                                                            hm_run.fill(namePrefix+'best_l1_muon'+qual_min_str+ptmin_str+iso_wp_str+dr_str+'_matched_probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.p', probeMomentumDict[probe_idx])
                                                            hm_run.fill(namePrefix+'best_l1_muon'+qual_min_str+ptmin_str+iso_wp_str+dr_str+'_matched_probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.eta', recoColl.eta[probe_idx])
                                                            hm_run.fill(namePrefix+'best_l1_muon'+qual_min_str+ptmin_str+iso_wp_str+dr_str+'_matched_probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.phi', recoColl.phi[probe_idx])
                                                            hm_run.fill(namePrefix+'best_l1_muon'+qual_min_str+ptmin_str+iso_wp_str+dr_str+'_matched_probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.charge', recoColl.charge[probe_idx])
                                                            hm_run.fill(namePrefix+'best_l1_muon'+qual_min_str+ptmin_str+iso_wp_str+dr_str+'_matched_probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.vtx', nVtx)
                                                            hm_run.fill(namePrefix+'best_l1_muon'+qual_min_str+ptmin_str+iso_wp_str+dr_str+'_matched_probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.run', runnr)
                                                            hm_run.fill(namePrefix+'best_l1_muon'+qual_min_str+ptmin_str+iso_wp_str+dr_str+'_matched_probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.dr', matched_l1_muons[i][2])
                                                            #hm_run.fill(namePrefix+'res_best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+dr_str+'_matched_l1_muon'+qual_min_str+ptmin_str+iso_wp_str+'.dpt', recoColl.pt[probe_idx] - (l1Coll.muonEt[matched_l1_muons[i][0]] * ptCorr))
                                                            #hm_run.fill(namePrefix+'res_best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+dr_str+'_matched_l1_muon'+qual_min_str+ptmin_str+iso_wp_str+'.dinvpt', 1./recoColl.pt[probe_idx] - 1./l1Coll.muonEt[matched_l1_muons[i][0]])
                                                            #hm_run.fill(namePrefix+'res_best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+dr_str+'_matched_l1_muon'+qual_min_str+ptmin_str+iso_wp_str+'.deta', recoColl.eta[probe_idx] - eta)
                                                            #hm_run.fill(namePrefix+'res_best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+dr_str+'_matched_l1_muon'+qual_min_str+ptmin_str+iso_wp_str+'.dphi', recoColl.phi[probe_idx] - phi)
                                                            #hm2d_run.fill(namePrefix+'2d_best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+dr_str+'_matched_l1_muon'+qual_min_str+ptmin_str+iso_wp_str+'.pt', recoColl.pt[probe_idx], l1Coll.muonEt[matched_l1_muons[i][0]] + ptCorr)
                                                            #hm2d_run.fill(namePrefix+'2d_best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+dr_str+'_matched_l1_muon'+qual_min_str+ptmin_str+iso_wp_str+'.eta', recoColl.eta[probe_idx], eta)
                                                            #hm2d_run.fill(namePrefix+'2d_best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+dr_str+'_matched_l1_muon'+qual_min_str+ptmin_str+iso_wp_str+'.phi', recoColl.phi[probe_idx], phi)
                                                            #hm2d_run.fill(namePrefix+'2d_best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+dr_str+'_matched_l1_muon'+qual_min_str+ptmin_str+iso_wp_str+'.charge', recoColl.charge[probe_idx], l1Coll.muonChg[matched_l1_muons[i][0]])
                                                    histo_filled = True
                                        hm.fill(namePrefix+'n_l1_muons'+qual_min_str+ptmin_str+iso_wp_str+dr_str+'_matched_to_a_probe'+eta_min_str+eta_max_str+probe_ptmin_str+'', l1_muon_cntr)
                                        if perRunHistos:
                                            hm_run.fill(namePrefix+'n_l1_muons'+qual_min_str+ptmin_str+iso_wp_str+dr_str+'_matched_to_a_probe'+eta_min_str+eta_max_str+probe_ptmin_str+'', l1_muon_cntr)

                                # fill all matched l1 muons
                                #fill_matched_muons(evt, hm, matched_l1_muons, 'u', eta_strs=[eta_min_str, eta_max_str], ptmin_strs=[probe_ptmin_str, ptmin_str])
            # probes need to be filled only once for the iso thresholds
            probesFilled = True

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
    useVtxExtraCoord = opts.extraCoord

    global iso_type
    if opts.isomethod == 'rel':
        iso_type = 1
    elif opts.isomethod == 'inner':
        iso_type = 2
    elif opts.isomethod == 'outovertot':
        iso_type = 3
    elif opts.isomethod == 'inner2x2':
        iso_type = 4
    elif opts.isomethod == 'outovertot2x2':
        iso_type = 5
    elif opts.isomethod == 'mipptadjust':
        iso_type = 6
    elif opts.isomethod == 'mipptadjust2':
        iso_type = 7
    else:
        iso_type = 0

    global isoEtaMax
    isoEtaMax = opts.etarestricted

    global prefix
    prefix = opts.prefix

    global tftype
    if opts.tftype == 'bmtf':
        tftype = 0
    elif opts.tftype == 'omtf':
        tftype = 1
    elif opts.tftype == 'emtf':
        tftype = 2

    emul = opts.emul

    # combinations of probe_pt_min and the corresponding pt_min values for a quality
    # the first line defines which thresholds are going to be used for unmatched histograms
    ptmins_list_q12 = [[0.5, [18, 20, 22, 24]],
                       [26, [18]],
                       [28, [20]],
                       [30, [22]],
                       [32, [24]],
                      ]

    ptmins_list_q8 = [[0.5, [3, 5, 8, 10, 12]],
                      [5, [3]],
                      [8, [5]],
                      [12, [8]],
                      [14, [10]],
                      [16, [12]],
                      ]

    ptmins_list_q4 = [[0.5, [3, 5, 8, 10, 12]],
                      [5, [3]],
                      [8, [5]],
                      [12, [8]],
                      [14, [10]],
                      [16, [12]],
                      ]

    eta_ranges = [[0, 2.4]]
#    eta_ranges = [[0, 2.4], [0, 0.83], [0.83, 1.24], [1.24, 2.4]]
#    eta_ranges = [[0, 2.4], [0, 0.83], [0.83, 1.24], [1.24, 2.4], [1.2, 1.55], [1.55, 1.85], [1.85, 2.4]]
    qual_ptmins_dict = {12:ptmins_list_q12, 8:ptmins_list_q8, 4:ptmins_list_q4}
    match_deltas = {'dr':0.5, 'deta':0.5, 'dphi':0.5} # max deltas for matching

    if iso_type == 0: # absolute isolation
        iso_wps = [0, 1, 3, 5, 7, 9, 11, 15, 20, 25, 28, 30, 31]
    elif iso_type == 1: # relative isolation
        iso_wps = [0., 1/2., 1/3., 2/3., 3/4., 4/5., 5/6., 6/7., 1., 1.5, 3., 5., 10., 31., 62.]
    elif iso_type == 3: # outer cone over total cone
        iso_wps = [0., 1/1., 1/2., 1/3., 2/3., 3/4., 4/5., 5/6., 6/7., 7/8., 8/9., 9/10., 19/20., 30/31., 99/100.]
    elif iso_type == 5: # outer cone over total cone 2x2
        iso_wps = [0., 1/1., 1/2., 1/3., 2/3., 3/4., 4/5., 5/6., 6/7., 7/8., 8/9., 9/10., 19/20., 30/31., 99/100.]
    elif iso_type == 7: # MIP pt adjust 2
        iso_wps = [2., 1., 0.5, 0.3, 0.1]
    else: # inner cone
        #iso_wps = [0., 1., 2., 3., 4., 5., 6., 7., 8., 9.]
        iso_wps = [0., 1., 2., 4., 6., 8., 10., 15.]

    # book the histograms
    L1Ana.log.info("Booking combined run histograms.")
    hm, hm2d = book_histograms(eta_ranges, qual_ptmins_dict, match_deltas, iso_wps, emul=emul)
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

            # process only if event has minimum or less than a maximum number of vertices
            if opts.nvtxmin:
                if event.recoVertex.nVtx < opts.nvtxmin:
                    continue
            if opts.nvtxmax:
                if event.recoVertex.nVtx > opts.nvtxmax:
                    continue

            # book histograms for this event's run number if not already done
            if perRunHistos and not runnr in hm_runs:
                L1Ana.log.info("Booking histograms for run {r}.".format(r=runnr))
                hm_runs[runnr], hm2d_runs[runnr] = book_histograms(eta_ranges, qual_ptmins_dict, match_deltas, iso_wps, emul=emul)

            # now do the analysis for all pt cut combinations
            if perRunHistos:
                analyse(event, hm, hm2d, hm_runs[runnr], hm2d_runs[runnr], eta_ranges, qual_ptmins_dict, match_deltas, iso_wps, emul=emul)
            else:
                analyse(event, hm, hm2d, None, None, eta_ranges, qual_ptmins_dict, match_deltas, iso_wps, emul=emul)
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
    iso_type = 0
    isoEtaMax = 3.
    prefix = ''
    tftype = -1
    ptCorrVal = 0.9
    saveHistos = True
    best_only = False
    perRunHistos = False
    main()

