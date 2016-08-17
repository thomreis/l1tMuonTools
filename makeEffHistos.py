#!/usr/bin/env python
from ToolBox import parse_options_and_init_log
# have to do this first or ROOT masks the -h messages
opts, parser = parse_options_and_init_log()

from L1Analysis import L1Ana, L1Ntuple
from analysis_tools.plotting import HistManager
from analysis_tools.selections import MuonSelections, Matcher
import ROOT as root

def parse_options_upgradeMuonHistos(parser):
    """
    Adds often used options to the OptionParser...
    """
    parsers = parser.add_subparsers()
    sub_parser = parsers.add_parser("makeEffHistos")
    sub_parser.add_argument("-o", "--outname", dest="outname", default="./ugmt_eff_histos.root", type=str, help="A root file name where to save the histograms.")

    opts, unknown = parser.parse_known_args()
    return opts

def book_histograms(eta_ranges, ptmins_list):
    # define pt binning
    pt_bins = range(0, 60, 2)
    pt_bins += range(60, 80, 5)
    pt_bins += range(80, 100, 10)
    pt_bins.append(100)

    vars_bins = [['pt', -1]+pt_bins, ['eta', 100, -2.5, 2.5], ['phi', 70, -3.5, 3.5]]
    x_title_vars = {'pt':'p_{T}', 'eta':'#eta', 'phi':'#phi'}
    x_title_units = {'pt':'GeV/c', 'eta':None, 'phi':None}

    varnames = []
    binnings = {}

    for eta_range in eta_ranges:
        eta_min = eta_range[0]
        eta_max = eta_range[1]
        eta_min_str = '_absEtaMin'+str(eta_min)
        eta_max_str = '_absEtaMax'+str(eta_max)
        eta_title = '{eMin} < |#eta| < {eMax}'.format(eMin=eta_min, eMax=eta_max)

        for ptmins in ptmins_list:
            reco_pt_min = ptmins[0]
            reco_ptmin_str = '_ptmin'+str(reco_pt_min)
            reco_cut_title = ' p_{T} > '+str(reco_pt_min)+' GeV/c'
            for pt_min in ptmins[1]:
                ptmin_str = '_ptmin'+str(pt_min)
                cut_title = ' p_{T} > '+str(pt_min)+' GeV/c'

                varnames.append('n_reco_muons'+eta_min_str+eta_max_str+reco_ptmin_str)
                varnames.append('n_ugmt_muons'+eta_min_str+eta_max_str+ptmin_str)
                varnames.append('n_reco'+reco_ptmin_str+'_matched_ugmt_muons'+eta_min_str+eta_max_str+ptmin_str)
                varnames.append('n_ugmt'+ptmin_str+'_matched_to_a_reco_muon'+eta_min_str+eta_max_str+reco_ptmin_str)
                for var_bin in vars_bins:
                    varnames.append('reco_muon'+eta_min_str+eta_max_str+reco_ptmin_str+'.'+var_bin[0])
                    varnames.append('ugmt_muon'+eta_min_str+eta_max_str+ptmin_str+'.'+var_bin[0])
                    varnames.append('best_ugmt'+ptmin_str+'_matched_reco_muon'+eta_min_str+eta_max_str+reco_ptmin_str+'.'+var_bin[0])
                    varnames.append('ugmt'+ptmin_str+'_matched_reco_muon'+eta_min_str+eta_max_str+reco_ptmin_str+'.'+var_bin[0])

                    for qual in range(16):
                        qual_str = '_q'+str(qual)

                        varnames.append('ugmt_muon'+eta_min_str+eta_max_str+ptmin_str+qual_str+'.'+var_bin[0])
                        varnames.append('best_ugmt'+ptmin_str+qual_str+'_matched_reco_muon'+eta_min_str+eta_max_str+reco_ptmin_str+'.'+var_bin[0])
                        varnames.append('ugmt'+ptmin_str+qual_str+'_matched_reco_muon'+eta_min_str+eta_max_str+reco_ptmin_str+'.'+var_bin[0])

                varnames.append('best_reco'+eta_min_str+eta_max_str+reco_ptmin_str+'_matched_ugmt_muon'+ptmin_str+'.pt')

                varnames.append('best_ugmt'+ptmin_str+'_matched_reco_muon'+eta_min_str+eta_max_str+reco_ptmin_str+'.dr')
                varnames.append('ugmt'+ptmin_str+'_matched_reco_muon'+eta_min_str+eta_max_str+reco_ptmin_str+'.dr')

                for qual in range(16):
                    qual_str = '_q'+str(qual)

                    varnames.append('best_ugmt'+ptmin_str+qual_str+'_matched_reco_muon'+eta_min_str+eta_max_str+reco_ptmin_str+'.dr')
                    varnames.append('ugmt'+ptmin_str+qual_str+'_matched_reco_muon'+eta_min_str+eta_max_str+reco_ptmin_str+'.dr')

                binnings['n_reco_muons'+eta_min_str+eta_max_str+reco_ptmin_str] = [10, 0, 10, '# reco #mu'+reco_cut_title]
                binnings['n_ugmt_muons'+eta_min_str+eta_max_str+ptmin_str] = [10, 0, 10, '# uGMT #mu'+cut_title]
                binnings['n_reco'+reco_ptmin_str+'_matched_ugmt_muons'+eta_min_str+eta_max_str+ptmin_str] = [10, 0, 10, '# uGMT #mu'+cut_title+' matched to reco #mu'+reco_cut_title]
                binnings['n_ugmt'+ptmin_str+'_matched_to_a_reco_muon'+eta_min_str+eta_max_str+reco_ptmin_str+''] = [5, 0, 5, '# uGMT #mu'+cut_title+' matched to one reco #mu'+reco_cut_title]
                for var_bin in vars_bins:
                    binnings['reco_muon'+eta_min_str+eta_max_str+reco_ptmin_str+'.'+var_bin[0]] = var_bin[1:]+[x_title_vars[var_bin[0]]+'^{reco}', x_title_units[var_bin[0]]]
                    binnings['ugmt_muon'+eta_min_str+eta_max_str+ptmin_str+'.'+var_bin[0]] = var_bin[1:]+[x_title_vars[var_bin[0]]+'^{uGMT}', x_title_units[var_bin[0]]]
                    binnings['best_ugmt'+ptmin_str+'_matched_reco_muon'+eta_min_str+eta_max_str+reco_ptmin_str+'.'+var_bin[0]] = var_bin[1:]+[x_title_vars[var_bin[0]]+'^{reco}', x_title_units[var_bin[0]]]
                    binnings['ugmt'+ptmin_str+'_matched_reco_muon'+eta_min_str+eta_max_str+reco_ptmin_str+'.'+var_bin[0]] = var_bin[1:]+[x_title_vars[var_bin[0]]+'^{reco}', x_title_units[var_bin[0]]]

                    for qual in range(16):
                        qual_str = '_q'+str(qual)
                        qualTitle = 'q'+str(qual)

                        binnings['ugmt_muon'+eta_min_str+eta_max_str+ptmin_str+qual_str+'.'+var_bin[0]] = var_bin[1:]+[x_title_vars[var_bin[0]]+'^{uGMT}', x_title_units[var_bin[0]]]
                        binnings['best_ugmt'+ptmin_str+qual_str+'_matched_reco_muon'+eta_min_str+eta_max_str+reco_ptmin_str+'.'+var_bin[0]] = var_bin[1:]+['('+qualTitle+') '+x_title_vars[var_bin[0]]+'^{reco}', x_title_units[var_bin[0]]]
                        binnings['ugmt'+ptmin_str+qual_str+'_matched_reco_muon'+eta_min_str+eta_max_str+reco_ptmin_str+'.'+var_bin[0]] = var_bin[1:]+['('+qualTitle+') '+x_title_vars[var_bin[0]]+'^{reco}', x_title_units[var_bin[0]]]

                binnings['best_reco'+eta_min_str+eta_max_str+reco_ptmin_str+'_matched_ugmt_muon'+ptmin_str+'.pt'] = vars_bins[0][1:]+[x_title_vars[vars_bins[0][0]]+'^{L1}', x_title_units[vars_bins[0][0]]]

                binnings['best_ugmt'+ptmin_str+'_matched_reco_muon'+eta_min_str+eta_max_str+reco_ptmin_str+'.dr'] = [60, 0., 0.6, '#Delta R']
                binnings['ugmt'+ptmin_str+'_matched_reco_muon'+eta_min_str+eta_max_str+reco_ptmin_str+'.dr'] = [60, 0., 0.6, '#Delta R']

                for qual in range(16):
                    qual_str = '_q'+str(qual)
                    qualTitle = 'q'+str(qual)

                    binnings['best_ugmt'+ptmin_str+qual_str+'_matched_reco_muon'+eta_min_str+eta_max_str+reco_ptmin_str+'.dr'] = [60, 0., 0.6, '('+qualTitle+') #Delta R']
                    binnings['ugmt'+ptmin_str+qual_str+'_matched_reco_muon'+eta_min_str+eta_max_str+reco_ptmin_str+'.dr'] = [60, 0., 0.6, '('+qualTitle+') #Delta R']

    return HistManager(list(set(varnames)), binnings)

def fill_matched_muons(evt, hm, matched_muons, muon_type='', eta_strs = ['', ''], ptmin_strs=['', '']):
    eta_min_str = eta_strs[0]
    eta_max_str = eta_strs[1]
    reco_ptmin_str = ptmin_strs[0]
    ptmin_str = ptmin_strs[1]

    recoColl = evt.recoMuon
    quals = evt.upgrade.muonQual
    if muon_type == 'u': # uGMT
        muon_str = 'ugmt'

    for matched_muon in matched_muons:
        # fill histograms with matched reco muon data
        hm.fill(muon_str+ptmin_str+'_matched_reco_muon'+eta_min_str+eta_max_str+reco_ptmin_str+'.pt', recoColl.pt[matched_muon[1]])
        hm.fill(muon_str+ptmin_str+'_matched_reco_muon'+eta_min_str+eta_max_str+reco_ptmin_str+'.eta', recoColl.eta[matched_muon[1]])
        hm.fill(muon_str+ptmin_str+'_matched_reco_muon'+eta_min_str+eta_max_str+reco_ptmin_str+'.phi', recoColl.phi[matched_muon[1]])
        hm.fill(muon_str+ptmin_str+'_matched_reco_muon'+eta_min_str+eta_max_str+reco_ptmin_str+'.dr', matched_muon[2])
        qual_str = '_q'+str(int(quals[matched_muon[0]]))
        hm.fill(muon_str+ptmin_str+qual_str+'_matched_reco_muon'+eta_min_str+eta_max_str+reco_ptmin_str+'.pt', recoColl.pt[matched_muon[1]])
        hm.fill(muon_str+ptmin_str+qual_str+'_matched_reco_muon'+eta_min_str+eta_max_str+reco_ptmin_str+'.eta', recoColl.eta[matched_muon[1]])
        hm.fill(muon_str+ptmin_str+qual_str+'_matched_reco_muon'+eta_min_str+eta_max_str+reco_ptmin_str+'.phi', recoColl.phi[matched_muon[1]])
        hm.fill(muon_str+ptmin_str+qual_str+'_matched_reco_muon'+eta_min_str+eta_max_str+reco_ptmin_str+'.dr', matched_muon[2])

def analyse(evt, hm, eta_ranges, ptmins_list, qualities):
    qual_min = qualities[0]

    bx_min = 0
    bx_max = 0

    # decide which ntuples to use
    emul = False
    if emul:
        ugmtColl = evt.upgradeEmu
    else:
        ugmtColl = evt.upgrade

    recoColl = evt.recoMuon
    vertexColl = evt.recoVertex

    reco_muon_idcs = MuonSelections.select_reco_muons(recoColl, pt_min=0.5, only_pos_eta=only_pos_eta)
    ugmt_muon_idcs = MuonSelections.select_ugmt_muons(ugmtColl, pt_min=0.5, qual_min=qual_min, bx_min=bx_min, bx_max=bx_max, only_pos_eta=only_pos_eta)

    thr_ugmt_muon_idcs_cache = {}

    for eta_range in eta_ranges:
        eta_min = eta_range[0]
        eta_max = eta_range[1]
        eta_min_str = '_absEtaMin'+str(eta_min)
        eta_max_str = '_absEtaMax'+str(eta_max)

        eta_reco_muon_idcs = MuonSelections.select_reco_muons(recoColl, abs_eta_min=eta_min, abs_eta_max=eta_max, idcs=reco_muon_idcs)
        eta_ugmt_muon_idcs = MuonSelections.select_ugmt_muons(ugmtColl, abs_eta_min=eta_min, abs_eta_max=eta_max, idcs=ugmt_muon_idcs)

        for ptmins in ptmins_list:
            reco_pt_min = ptmins[0]
            reco_ptmin_str = '_ptmin'+str(reco_pt_min)

            eta_thr_reco_muon_idcs = MuonSelections.select_reco_muons(recoColl, pt_min=reco_pt_min, idcs=eta_reco_muon_idcs)
            # fill the histograms with the recoerated muon values
            hm.fill('n_reco_muons'+eta_min_str+eta_max_str+reco_ptmin_str, len(eta_thr_reco_muon_idcs))
            for i in eta_thr_reco_muon_idcs:
                hm.fill('reco_muon'+eta_min_str+eta_max_str+reco_ptmin_str+'.pt', recoColl.pt[i])
                hm.fill('reco_muon'+eta_min_str+eta_max_str+reco_ptmin_str+'.eta', recoColl.eta[i])
                hm.fill('reco_muon'+eta_min_str+eta_max_str+reco_ptmin_str+'.phi', recoColl.phi[i])

            for pt_min in ptmins[1]:
                ptmin_str = '_ptmin'+str(pt_min)

                thr_ugmt_muon_idcs = MuonSelections.select_ugmt_muons(ugmtColl, pt_min=pt_min, idcs=ugmt_muon_idcs)
                thr_ugmt_muon_idcs_cache[pt_min] = thr_ugmt_muon_idcs

                eta_thr_ugmt_muon_idcs = MuonSelections.select_ugmt_muons(ugmtColl, pt_min=pt_min, idcs=eta_ugmt_muon_idcs)

                if reco_pt_min == ptmins_list[0][0]:
                    # fill the histograms with the ugmt muon values
                    hm.fill('n_ugmt_muons'+eta_min_str+eta_max_str+ptmin_str, len(eta_thr_ugmt_muon_idcs))
                    for i in eta_thr_ugmt_muon_idcs:
                        hm.fill('ugmt_muon'+eta_min_str+eta_max_str+ptmin_str+'.pt', ugmtColl.muonEt[i])
                        hm.fill('ugmt_muon'+eta_min_str+eta_max_str+ptmin_str+'.eta', ugmtColl.muonEta[i])
                        hm.fill('ugmt_muon'+eta_min_str+eta_max_str+ptmin_str+'.phi', ugmtColl.muonPhi[i])
                        qual_str = '_q'+str(int(ugmtColl.muonQual[i]))
                        hm.fill('ugmt_muon'+eta_min_str+eta_max_str+ptmin_str+qual_str+'.pt', ugmtColl.muonEt[i])
                        hm.fill('ugmt_muon'+eta_min_str+eta_max_str+ptmin_str+qual_str+'.eta', ugmtColl.muonEta[i])
                        hm.fill('ugmt_muon'+eta_min_str+eta_max_str+ptmin_str+qual_str+'.phi', ugmtColl.muonPhi[i])

                if len(eta_thr_reco_muon_idcs) > 0:
                    ##########################################################################
                    # match selected ugmt muons to selected recoerated muons
                    matched_ugmt_muons = Matcher.match_dr(ugmtColl.muonEta, ugmtColl.muonPhi, recoColl.eta, recoColl.phi, idcs1=thr_ugmt_muon_idcs, idcs2=eta_thr_reco_muon_idcs)
                    hm.fill('n_reco'+reco_ptmin_str+'_matched_ugmt_muons'+eta_min_str+eta_max_str+ptmin_str, len(matched_ugmt_muons))

                    # how many ugmt matches did we find for each recoerated muon
                    for reco_muon_idx in eta_thr_reco_muon_idcs:
                        ugmt_muon_cntr = 0
                        histo_filled = False
                        for i in range(len(matched_ugmt_muons)):
                            if reco_muon_idx == matched_ugmt_muons[i][1]:
                                ugmt_muon_cntr += 1
                                # fill muon values only for the first (and therefore best) match to this reco muon
                                if not histo_filled:
                                    hm.fill('best_reco'+eta_min_str+eta_max_str+reco_ptmin_str+'_matched_ugmt_muon'+ptmin_str+'.pt', ugmtColl.muonEt[matched_ugmt_muons[i][0]])
                                    hm.fill('best_ugmt'+ptmin_str+'_matched_reco_muon'+eta_min_str+eta_max_str+reco_ptmin_str+'.pt', recoColl.pt[reco_muon_idx])
                                    hm.fill('best_ugmt'+ptmin_str+'_matched_reco_muon'+eta_min_str+eta_max_str+reco_ptmin_str+'.eta', recoColl.eta[reco_muon_idx])
                                    hm.fill('best_ugmt'+ptmin_str+'_matched_reco_muon'+eta_min_str+eta_max_str+reco_ptmin_str+'.phi', recoColl.phi[reco_muon_idx])
                                    hm.fill('best_ugmt'+ptmin_str+'_matched_reco_muon'+eta_min_str+eta_max_str+reco_ptmin_str+'.dr', matched_ugmt_muons[i][2])
                                    qual_str = '_q'+str(int(ugmtColl.muonQual[matched_ugmt_muons[i][0]]))
                                    hm.fill('best_ugmt'+ptmin_str+qual_str+'_matched_reco_muon'+eta_min_str+eta_max_str+reco_ptmin_str+'.pt', recoColl.pt[reco_muon_idx])
                                    hm.fill('best_ugmt'+ptmin_str+qual_str+'_matched_reco_muon'+eta_min_str+eta_max_str+reco_ptmin_str+'.eta', recoColl.eta[reco_muon_idx])
                                    hm.fill('best_ugmt'+ptmin_str+qual_str+'_matched_reco_muon'+eta_min_str+eta_max_str+reco_ptmin_str+'.phi', recoColl.phi[reco_muon_idx])
                                    hm.fill('best_ugmt'+ptmin_str+qual_str+'_matched_reco_muon'+eta_min_str+eta_max_str+reco_ptmin_str+'.dr', matched_ugmt_muons[i][2])
                                    histo_filled = True
                        hm.fill('n_ugmt'+ptmin_str+'_matched_to_a_reco_muon'+eta_min_str+eta_max_str+reco_ptmin_str+'', ugmt_muon_cntr)

                    # fill all matched ugmt muons
                    fill_matched_muons(evt, hm, matched_ugmt_muons, 'u', eta_strs=[eta_min_str, eta_max_str], ptmin_strs=[reco_ptmin_str, ptmin_str])

def save_histos(hm, outfile):
    '''
    save all histograms in hm to outfile
    '''
    outfile.cd()
    for varname in hm.get_varnames():
        hm.get(varname).Write()

def main():
    L1Ana.init_l1_analysis()
    opts = parse_options_upgradeMuonHistos(parser)
    print ""

    # combinations of reco_pt_min and the corresponding pt_min values
    # the first line defines which thresholds are going to be used for unmatched histograms
    ptmins_list = [[0.5, [0.5, 12, 16, 20, 24, 30]],
#                  [12, [12]],
#                  [16, [12, 16]],
#                  [20, [12, 16, 20]],
                  [24, [16, 20, 24]],
#                  [30, [20, 24, 30]],
                 ]

    eta_ranges = [[0, 2.5], [0, 0.83], [0.83, 1.24], [1.24, 2.5]]
    qualities = [12, 8] # [uGMT, GMT]

    # book the histograms
    hm = book_histograms(eta_ranges, ptmins_list)

    ntuple = L1Ntuple(opts.nevents)

    if opts.flist:
        ntuple.open_with_file_list(opts.flist)
    if opts.fname:
        ntuple.open_with_file(opts.fname)

    start_evt = opts.start_event
    end_evt = opts.start_event+ntuple.nevents
    for i in range(start_evt, end_evt):
        event = ntuple[i]
        if (i+1) % 1000 == 0:
            L1Ana.log.info("Processing event: {n}".format(n=i+1))
        # now do the analysis for all pt cut combinations
        analyse(event, hm, eta_ranges, ptmins_list, qualities)

    # save histos to root file
    if saveHistos:
        output = root.TFile(opts.outname, 'recreate')
        output.cd()
        save_histos(hm, output)
        output.Close()

if __name__ == "__main__":
    only_pos_eta = False
    saveHistos = True
    best_only = False
    main()

