#!/usr/bin/env python
from ToolBox import parse_options_and_init_log
# have to do this first or ROOT masks the -h messages
opts, parser = parse_options_and_init_log()

from L1Analysis import L1Ana, L1Ntuple
from analysis_tools.plotting import HistManager
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

    opts, unknown = parser.parse_known_args()
    return opts

def get_tftype(tf_muon_index):
    if tf_muon_index > 35 and tf_muon_index < 72:
        return 0 # BMTF
    elif tf_muon_index > 17 and tf_muon_index < 90:
        return 1 # OMTF
    else:
        return 2 # EMTF

def book_histograms(eta_ranges, qual_ptmins_dict, match_deltas, emul=False):
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
    dphi_str = '_dphi'+str(match_deltas['dphi'])
    deta_str = '_deta'+str(match_deltas['deta'])

    vars_bins = [['pt', -1]+pt_bins, ['eta', 100, -2.5, 2.5], ['phi', 70, -3.5, 3.5], ['charge', 3, -1, 2], ['vtx', 60, 0, 60], ['run', 17000, 271725, 288725]]
    x_title_vars = {'pt':'p_{T}', 'eta':'#eta', 'phi':'#phi', 'charge':'charge', 'vtx':'PU', 'run':'run number'}
    x_title_units = {'pt':'GeV/c', 'eta':None, 'phi':None, 'charge':None, 'vtx':None, 'run':None}
    probe_vars_bins = vars_bins + [['p', -1]+p_bins]
    probe_x_title_vars = {'p':'p'}
    probe_x_title_vars.update(x_title_vars)
    probe_x_title_units = {'p':'GeV/c'}
    probe_x_title_units.update(x_title_units)
    res_vars_bins = [['dpt', 100, -50, 50], ['dinvpt', 200, -0.2, 0.2], ['deta', 100, -0.1, 0.1], ['dphi', 100, -0.2, 0.2]]
    res_x_title_vars = {'dpt':'p_{T}^{reco} - p_{T}^{L1}', 'dinvpt':'1/p_{T}^{reco} - 1/p_{T}^{L1}', 'deta':'#eta_{reco} - #eta_{L1}', 'dphi':'#phi_{reco} - #phi_{L1}'}
    res_x_title_units = {'dpt':'GeV', 'dinvpt':'GeV', 'deta':None, 'dphi':None}

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
                        varnames.append(namePrefix+'n_probe'+probe_ptmin_str+dr_str+'_matched_l1_muons'+eta_min_str+eta_max_str+qual_min_str+ptmin_str)
                        varnames.append(namePrefix+'n_l1_muons'+qual_min_str+ptmin_str+dr_str+'_matched_to_a_probe'+eta_min_str+eta_max_str+probe_ptmin_str)
                        for var_bin in vars_bins:
                            varnames.append(namePrefix+'l1_muon'+eta_min_str+eta_max_str+qual_min_str+ptmin_str+'.'+var_bin[0])
                            varnames.append(namePrefix+'best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+dr_str+'_matched_l1_muon'+qual_min_str+ptmin_str+'.'+var_bin[0])
                        for var_bin in probe_vars_bins:
                            varnames.append(namePrefix+'probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.'+var_bin[0])
                            varnames.append(namePrefix+'best_l1_muon'+qual_min_str+ptmin_str+dr_str+'_matched_probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.'+var_bin[0])
                        for res_var_bin in res_vars_bins:
                            varnames.append(namePrefix+'res_best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+dr_str+'_matched_l1_muon'+qual_min_str+ptmin_str+'.'+res_var_bin[0])

                        varnames.append(namePrefix+'best_l1_muon'+qual_min_str+ptmin_str+dr_str+'_matched_probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.dr')

                        # binnings
                        binnings[namePrefix+'n_probes'+eta_min_str+eta_max_str+probe_ptmin_str] = [10, 0, 10, '# probes'+probe_cut_title]
                        binnings[namePrefix+'n_probe'+probe_ptmin_str+dr_str+'_matched_l1_muons'+eta_min_str+eta_max_str+qual_min_str+ptmin_str] = [10, 0, 10, '# uGMT #mu'+cut_title+' #Delta R matched to probe'+probe_cut_title]
                        binnings[namePrefix+'n_l1_muons'+qual_min_str+ptmin_str+dr_str+'_matched_to_a_probe'+eta_min_str+eta_max_str+probe_ptmin_str+''] = [5, 0, 5, '# uGMT #mu'+cut_title+' #Delta R matched to one probe'+probe_cut_title]
                        for var_bin in vars_bins:
                            binnings[namePrefix+'l1_muon'+eta_min_str+eta_max_str+qual_min_str+ptmin_str+'.'+var_bin[0]] = var_bin[1:]+[x_title_vars[var_bin[0]]+'^{uGMT}', x_title_units[var_bin[0]]]
                            binnings[namePrefix+'best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+dr_str+'_matched_l1_muon'+qual_min_str+ptmin_str+'.'+var_bin[0]] = vars_bins[0][1:]+[x_title_vars[vars_bins[0][0]]+'^{L1}', x_title_units[vars_bins[0][0]]]
                        for var_bin in probe_vars_bins:
                            binnings[namePrefix+'probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.'+var_bin[0]] = var_bin[1:]+[probe_x_title_vars[var_bin[0]]+'^{reco}', probe_x_title_units[var_bin[0]]]
                            binnings[namePrefix+'best_l1_muon'+qual_min_str+ptmin_str+dr_str+'_matched_probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.'+var_bin[0]] = var_bin[1:]+[probe_x_title_vars[var_bin[0]]+'^{reco}', probe_x_title_units[var_bin[0]]]
                        for res_var_bin in res_vars_bins:
                            binnings[namePrefix+'res_best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+dr_str+'_matched_l1_muon'+qual_min_str+ptmin_str+'.'+res_var_bin[0]] = res_var_bin[1:]+[res_x_title_vars[res_var_bin[0]], res_x_title_units[res_var_bin[0]]]

                        binnings[namePrefix+'best_l1_muon'+qual_min_str+ptmin_str+dr_str+'_matched_probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.dr'] = [60, 0., 0.6, '#Delta R']

    return HistManager(list(set(varnames)), binnings)

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

def analyse(evt, hm, eta_ranges, qual_ptmins_dict, match_deltas, emul=False):
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
    l1_muon_idcs = MuonSelections.select_ugmt_muons(l1Coll, pt_min=0.5, bx_min=bx_min, bx_max=bx_max, pos_eta=pos_eta, neg_eta=neg_eta)

    # vertex information
    nVtx = evt.recoVertex.nVtx
    # run number
    runnr = evt.event.run

    matchdr = match_deltas['dr']
    matchdphi = match_deltas['dphi']
    matchdeta = match_deltas['deta']
    dr_str = '_dr'+str(matchdr)
    dphi_str = '_dphi'+str(matchdphi)
    deta_str = '_deta'+str(matchdeta)

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

        # for all defined eta ranges
        for eta_range in eta_ranges:
            eta_min = eta_range[0]
            eta_max = eta_range[1]
            eta_min_str = '_absEtaMin'+str(eta_min)
            eta_max_str = '_absEtaMax'+str(eta_max)

            eta_probe_idcs = MuonSelections.select_reco_muons(recoColl, abs_eta_min=eta_min, abs_eta_max=eta_max, idcs=invmass_probe_idcs)
            eta_l1_muon_idcs = MuonSelections.select_ugmt_muons(l1Coll, abs_eta_min=eta_min, abs_eta_max=eta_max, idcs=l1_muon_idcs)

            # keep probe pt cuts in a list to not fill the histograms several times if two quality cuts use the same probe pt cut
            probe_pt_mins = []
            probe_ptmin_str_dict = {}
            eta_thr_probe_idcs_dict = {}
            # for all defined min quality ptmin_list combinations
            for q in range(16):
                if q in qual_ptmins_dict:
                    ptmins_list = qual_ptmins_dict[q]
                    qual_min_str = '_qualMin'+str(q)
                    eta_q_l1_muon_idcs = MuonSelections.select_ugmt_muons(l1Coll, qual_min=q, idcs=eta_l1_muon_idcs)

                    # for all defined min probe pt
                    for ptmins in ptmins_list:
                        if not ptmins[0] in probe_pt_mins: # fill only once for each pt min value
                            probe_pt_mins.append(ptmins[0])
                            probe_ptmin_str = '_ptmin'+str(ptmins[0])
                            probe_ptmin_str_dict[ptmins[0]] = probe_ptmin_str

                            eta_thr_probe_idcs = MuonSelections.select_reco_muons(recoColl, pt_min=probe_pt_mins[-1], idcs=eta_probe_idcs)
                            eta_thr_probe_idcs_dict[ptmins[0]] = eta_thr_probe_idcs
                            # fill the histograms with the probe kinematics
                            hm.fill(namePrefix+'n_probes'+eta_min_str+eta_max_str+probe_ptmin_str, len(eta_thr_probe_idcs))
                            for i in eta_thr_probe_idcs:
                                hm.fill(namePrefix+'probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.pt', recoColl.pt[i])
                                hm.fill(namePrefix+'probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.p', probeMomentumDict[i])
                                hm.fill(namePrefix+'probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.eta', recoColl.eta[i])
                                hm.fill(namePrefix+'probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.phi', recoColl.phi[i])
                                hm.fill(namePrefix+'probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.charge', recoColl.charge[i])
                                hm.fill(namePrefix+'probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.vtx', nVtx)
                                hm.fill(namePrefix+'probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.run', runnr)

                        probe_ptmin_str = probe_ptmin_str_dict[ptmins[0]]
                        eta_thr_probe_idcs = eta_thr_probe_idcs_dict[ptmins[0]]
                        # for all defined min l1 muon pt
                        for pt_min in ptmins[1]:
                            ptmin_str = '_ptmin'+str(pt_min)

                            q_thr_l1_muon_idcs = MuonSelections.select_ugmt_muons(l1Coll, qual_min=q, pt_min=pt_min, idcs=l1_muon_idcs)
                            eta_q_thr_l1_muon_idcs = MuonSelections.select_ugmt_muons(l1Coll, pt_min=pt_min, idcs=eta_q_l1_muon_idcs)
                            # fill the histograms with the l1 muon kinematics
                            hm.fill(namePrefix+'n_probes'+eta_min_str+eta_max_str+probe_ptmin_str, len(eta_thr_probe_idcs))
                            for i in eta_q_thr_l1_muon_idcs:
                                if tftype == -1 or tftype == get_tftype(l1Coll.muonTfMuonIdx[i]):
                                    hm.fill(namePrefix+'l1_muon'+eta_min_str+eta_max_str+qual_min_str+ptmin_str+'.pt', l1Coll.muonEt[i])
                                    hm.fill(namePrefix+'l1_muon'+eta_min_str+eta_max_str+qual_min_str+ptmin_str+'.eta', l1Coll.muonEta[i])
                                    hm.fill(namePrefix+'l1_muon'+eta_min_str+eta_max_str+qual_min_str+ptmin_str+'.phi', l1Coll.muonPhi[i])
                                    hm.fill(namePrefix+'l1_muon'+eta_min_str+eta_max_str+qual_min_str+ptmin_str+'.charge', l1Coll.muonChg[i])
                                    hm.fill(namePrefix+'l1_muon'+eta_min_str+eta_max_str+qual_min_str+ptmin_str+'.vtx', nVtx)
                                    hm.fill(namePrefix+'l1_muon'+eta_min_str+eta_max_str+qual_min_str+ptmin_str+'.run', runnr)


                            if len(q_thr_l1_muon_idcs) > 0:
                                zeros = [0.] * 50 # list with zeros for deta and dphi matching with the Matcher.match_dr function
                                # match selected l1 muons to selected probes
                                matched_l1_muons = Matcher.match_dr(l1Coll.muonEta, l1Coll.muonPhi, recoColl.eta, recoColl.phi, cut=matchdr, idcs1=q_thr_l1_muon_idcs, idcs2=eta_thr_probe_idcs) # match in delta R
                                hm.fill(namePrefix+'n_probe'+probe_ptmin_str+dr_str+'_matched_l1_muons'+eta_min_str+eta_max_str+qual_min_str+ptmin_str, len(matched_l1_muons))

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
                                                    hm.fill(namePrefix+'best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+dr_str+'_matched_l1_muon'+qual_min_str+ptmin_str+'.pt', l1Coll.muonEt[matched_l1_muons[i][0]])
                                                    hm.fill(namePrefix+'best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+dr_str+'_matched_l1_muon'+qual_min_str+ptmin_str+'.eta', l1Coll.muonEta[matched_l1_muons[i][0]])
                                                    hm.fill(namePrefix+'best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+dr_str+'_matched_l1_muon'+qual_min_str+ptmin_str+'.phi', l1Coll.muonPhi[matched_l1_muons[i][0]])
                                                    hm.fill(namePrefix+'best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+dr_str+'_matched_l1_muon'+qual_min_str+ptmin_str+'.charge', l1Coll.muonChg[matched_l1_muons[i][0]])
                                                    hm.fill(namePrefix+'best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+dr_str+'_matched_l1_muon'+qual_min_str+ptmin_str+'.vtx', nVtx)
                                                    hm.fill(namePrefix+'best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+dr_str+'_matched_l1_muon'+qual_min_str+ptmin_str+'.run', runnr)
                                                    hm.fill(namePrefix+'best_l1_muon'+qual_min_str+ptmin_str+dr_str+'_matched_probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.pt', recoColl.pt[probe_idx])
                                                    hm.fill(namePrefix+'best_l1_muon'+qual_min_str+ptmin_str+dr_str+'_matched_probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.p', probeMomentumDict[probe_idx])
                                                    hm.fill(namePrefix+'best_l1_muon'+qual_min_str+ptmin_str+dr_str+'_matched_probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.eta', recoColl.eta[probe_idx])
                                                    hm.fill(namePrefix+'best_l1_muon'+qual_min_str+ptmin_str+dr_str+'_matched_probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.phi', recoColl.phi[probe_idx])
                                                    hm.fill(namePrefix+'best_l1_muon'+qual_min_str+ptmin_str+dr_str+'_matched_probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.charge', recoColl.charge[probe_idx])
                                                    hm.fill(namePrefix+'best_l1_muon'+qual_min_str+ptmin_str+dr_str+'_matched_probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.vtx', nVtx)
                                                    hm.fill(namePrefix+'best_l1_muon'+qual_min_str+ptmin_str+dr_str+'_matched_probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.run', runnr)
                                                    hm.fill(namePrefix+'best_l1_muon'+qual_min_str+ptmin_str+dr_str+'_matched_probe'+eta_min_str+eta_max_str+probe_ptmin_str+'.dr', matched_l1_muons[i][2])
                                                    hm.fill(namePrefix+'res_best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+dr_str+'_matched_l1_muon'+qual_min_str+ptmin_str+'.dpt', recoColl.pt[probe_idx] - l1Coll.muonEt[matched_l1_muons[i][0]])
                                                    hm.fill(namePrefix+'res_best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+dr_str+'_matched_l1_muon'+qual_min_str+ptmin_str+'.dinvpt', 1./recoColl.pt[probe_idx] - 1./l1Coll.muonEt[matched_l1_muons[i][0]])
                                                    hm.fill(namePrefix+'res_best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+dr_str+'_matched_l1_muon'+qual_min_str+ptmin_str+'.deta', recoColl.eta[probe_idx] - l1Coll.muonEta[matched_l1_muons[i][0]])
                                                    hm.fill(namePrefix+'res_best_probe'+eta_min_str+eta_max_str+probe_ptmin_str+dr_str+'_matched_l1_muon'+qual_min_str+ptmin_str+'.dphi', recoColl.phi[probe_idx] - l1Coll.muonPhi[matched_l1_muons[i][0]])

                                                histo_filled = True
                                    hm.fill(namePrefix+'n_l1_muons'+qual_min_str+ptmin_str+dr_str+'_matched_to_a_probe'+eta_min_str+eta_max_str+probe_ptmin_str+'', l1_muon_cntr)



                                # fill all matched l1 muons
                                #fill_matched_muons(evt, hm, matched_l1_muons, 'u', eta_strs=[eta_min_str, eta_max_str], ptmin_strs=[probe_ptmin_str, ptmin_str])

def save_histos(hm, outfile):
    '''
    save all histograms in hm to outfile
    '''
    outfile.mkdir('all_runs')
    outfile.cd('all_runs')
    for varname in hm.get_varnames():
        hm.get(varname).Write()
        

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

    emul = opts.emul

    # combinations of probe_pt_min and the corresponding pt_min values for a quality
    # the first line defines which thresholds are going to be used for unmatched histograms
    ptmins_list_q12 = [[0.5, [0.5, 22]],
                       [30, [22]],
                       [100, [22]],
                      ]

    ptmins_list_q8 = [[0.5, [0.5, 5, 12, 22]],
                      [8, [5]],
                      [16, [12]],
                      [30, [22]],
                      ]

    ptmins_list_q4 = [[0.5, [0.5, 5, 12, 22]],
                      [8, [5]],
                      [16, [12]],
                      [30, [22]],
                      ]

#    eta_ranges = [[0, 2.4]]
    eta_ranges = [[0, 2.4], [0, 0.83], [0.83, 1.24], [1.24, 2.4]]
#    eta_ranges = [[0, 2.4], [0, 0.83], [0.83, 1.24], [1.24, 2.4], [1.2, 1.55], [1.55, 1.85], [1.85, 2.4]]
#    qual_ptmins_dict = {12:ptmins_list_q12, 8:ptmins_list_q8, 4:ptmins_list_q4}
    qual_ptmins_dict = {12:ptmins_list_q12}
    match_deltas = {'dr':0.5, 'deta':0.5, 'dphi':0.5} # max deltas for matching

    # book the histograms
    L1Ana.log.info("Booking combined run histograms.")
    hm = book_histograms(eta_ranges, qual_ptmins_dict, match_deltas, emul=emul)

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

            # now do the analysis for all pt cut combinations
            analyse(event, hm, eta_ranges, qual_ptmins_dict, match_deltas, emul=emul)
            analysed_evt_ctr += 1
    except KeyboardInterrupt:
        L1Ana.log.info("Analysis interrupted after {n} events".format(n=i))

    L1Ana.log.info("Analysis of {nAna} events in selected runs/LS finished.".format(nAna=analysed_evt_ctr))

    # save histos to root file
    if saveHistos:
        output = root.TFile(opts.outname, 'recreate')
        output.cd()
        save_histos(hm, output)
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
    saveHistos = True
    best_only = False
    main()

