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

ptScale = 0.5
etaScale = 0.010875
phiScale = 0.010908

def parse_options_upgradeRateHistos(parser):
    """
    Adds often used options to the OptionParser...
    """
    parsers = parser.add_subparsers()
    sub_parser = parsers.add_parser("makeRateHistos")
    sub_parser.add_argument("-o", "--outname", dest="outname", default="./ugmt_rate_histos.root", type=str, help="A root file name where to save the histograms.")
    sub_parser.add_argument("-j", "--json", dest="json", type=str, default=None, help="A json file with good lumi sections per run.")
    sub_parser.add_argument("-e", "--emul", dest="emul", action='store_true', help="Use emulated collections instead of unpacked ones.")
    sub_parser.add_argument("--use-extra-coord", dest="extraCoord", default=False, action="store_true", help="Use L1 extrapolated eta and phi coordinates.")
    sub_parser.add_argument("--eta-restricted", dest="etarestricted", type=float, default=3., help="Upper eta value for isolation.")
    sub_parser.add_argument("--iso-method", dest="isomethod", type=str, default='abs', help="Isolation method. ['abs', 'rel', 'inner', 'outovertot', 'inner2x2', 'outovertot2x2']")

    opts, unknown = parser.parse_known_args()
    return opts

def book_histograms(eta_ranges, thresholds, qualities, iso_wps):
    varnames = []
    binnings = {}

    varnames.append('n_evts_analysed')
    binnings['n_evts_analysed'] = [1, 0, 1, '']

    # define pt binning for GMT and OMTF
    pt_bins = [-1] # -1 indicates a variable binning
    pt_bins += range(8)
    pt_bins += range(8, 20, 2)
    pt_bins += range(20, 50, 5)
    pt_bins += range(50, 110, 10)

    # define eta binning
    eta_bins = [-1] # -1 indicates a variable binning
    eta_bins += [0.01*i for i in range(-250, -170, 5)]
    eta_bins += [0.1*i for i in range(-17, 0)]
    eta_bins += [0.1*i for i in range(18)]
    eta_bins += [0.01*i for i in range(175, 255, 5)]

    for iso_wp in iso_wps:
        iso_wp_str = '_isoMax{iso:.3f}'.format(iso=iso_wp)
        for eta_range in eta_ranges:
            eta_min = eta_range[0]
            eta_max = eta_range[1]
            eta_min_str = '_absEtaMin'+str(eta_min)
            eta_max_str = '_absEtaMax'+str(eta_max)
            eta_title = '{eMin} < |#eta| < {eMax}'.format(eMin=eta_min, eMax=eta_max)

            for threshold in thresholds:
                thr_str = '_ptmin'+str(threshold)
                thr_title = 'p_{T} > '+str(threshold)+' GeV/c'

                #varnames.append('gmt_muon'+eta_min_str+eta_max_str+thr_str+iso_wp_str+'_qual')
                #varnames.append('ugmt_muon'+eta_min_str+eta_max_str+thr_str+iso_wp_str+'_qual')
                #varnames.append('bmtf_ugmt_muon'+eta_min_str+eta_max_str+thr_str+iso_wp_str+'_qual')
                #varnames.append('omtf_ugmt_muon'+eta_min_str+eta_max_str+thr_str+iso_wp_str+'_qual')
                #varnames.append('emtf_ugmt_muon'+eta_min_str+eta_max_str+thr_str+iso_wp_str+'_qual')
                #binnings['gmt_muon'+eta_min_str+eta_max_str+thr_str+iso_wp_str+'_qual'] = [16, 0, 16, '('+eta_title+', '+thr_title+') quality']
                #binnings['ugmt_muon'+eta_min_str+eta_max_str+thr_str+iso_wp_str+'_qual'] = [16, 0, 16, '('+eta_title+', '+thr_title+') quality']
                #binnings['bmtf_ugmt_muon'+eta_min_str+eta_max_str+thr_str+iso_wp_str+'_qual'] = [16, 0, 16, '('+eta_title+', '+thr_title+') quality']
                #binnings['omtf_ugmt_muon'+eta_min_str+eta_max_str+thr_str+iso_wp_str+'_qual'] = [16, 0, 16, '('+eta_title+', '+thr_title+') quality']
                #binnings['emtf_ugmt_muon'+eta_min_str+eta_max_str+thr_str+iso_wp_str+'_qual'] = [16, 0, 16, '('+eta_title+', '+thr_title+') quality']

                for qMin in qualities:
                    qMin_str = '_qmin'+str(qMin)
                    qTitle = 'q>'+str(qMin)

                    #varnames.append('n_gmt_muons'+eta_min_str+eta_max_str+thr_str+iso_wp_str+qMin_str)
                    varnames.append('n_ugmt_muons'+eta_min_str+eta_max_str+thr_str+iso_wp_str+qMin_str)
                    #varnames.append('n_bmtf_ugmt_muons'+eta_min_str+eta_max_str+thr_str+iso_wp_str+qMin_str)
                    #varnames.append('n_omtf_ugmt_muons'+eta_min_str+eta_max_str+thr_str+iso_wp_str+qMin_str)
                    #varnames.append('n_emtf_ugmt_muons'+eta_min_str+eta_max_str+thr_str+iso_wp_str+qMin_str)
                    #binnings['n_gmt_muons'+eta_min_str+eta_max_str+thr_str+iso_wp_str+qMin_str] = [10, 0, 10, '# GMT #mu '+eta_title+', '+thr_title+', '+qTitle]
                    binnings['n_ugmt_muons'+eta_min_str+eta_max_str+thr_str+iso_wp_str+qMin_str] = [10, 0, 10, '# uGMT #mu '+eta_title+', '+thr_title+', '+qTitle]
                    #binnings['n_bmtf_ugmt_muons'+eta_min_str+eta_max_str+thr_str+iso_wp_str+qMin_str] = [10, 0, 10, '# BMTF uGMT #mu '+eta_title+', '+thr_title+', '+qTitle]
                    #binnings['n_omtf_ugmt_muons'+eta_min_str+eta_max_str+thr_str+iso_wp_str+qMin_str] = [10, 0, 10, '# OMTF uGMT #mu '+eta_title+', '+thr_title+', '+qTitle]
                    #binnings['n_emtf_ugmt_muons'+eta_min_str+eta_max_str+thr_str+iso_wp_str+qMin_str] = [10, 0, 10, '# EMTF uGMT #mu '+eta_title+', '+thr_title+', '+qTitle]

                    #varnames.append('gmt_muon'+eta_min_str+eta_max_str+thr_str+iso_wp_str+qMin_str+'_phi')
                    #varnames.append('ugmt_muon'+eta_min_str+eta_max_str+thr_str+iso_wp_str+qMin_str+'_phi')
                    #varnames.append('bmtf_ugmt_muon'+eta_min_str+eta_max_str+thr_str+iso_wp_str+qMin_str+'_phi')
                    #varnames.append('omtf_ugmt_muon'+eta_min_str+eta_max_str+thr_str+iso_wp_str+qMin_str+'_phi')
                    #varnames.append('emtf_ugmt_muon'+eta_min_str+eta_max_str+thr_str+iso_wp_str+qMin_str+'_phi')
                    #binnings['gmt_muon'+eta_min_str+eta_max_str+thr_str+iso_wp_str+qMin_str+'_phi'] = [35, -3.5, 3.5, 'GMT #mu ('+eta_title+', '+thr_title+', '+qTitle+') #phi']
                    #binnings['ugmt_muon'+eta_min_str+eta_max_str+thr_str+iso_wp_str+qMin_str+'_phi'] = [35, -3.5, 3.5, 'uGMT #mu ('+eta_title+', '+thr_title+', '+qTitle+') #phi']
                    #binnings['bmtf_ugmt_muon'+eta_min_str+eta_max_str+thr_str+iso_wp_str+qMin_str+'_phi'] = [35, -3.5, 3.5, 'BMTF uGMT #mu ('+eta_title+', '+thr_title+', '+qTitle+') #phi']
                    #binnings['omtf_ugmt_muon'+eta_min_str+eta_max_str+thr_str+iso_wp_str+qMin_str+'_phi'] = [35, -3.5, 3.5, 'OMTF uGMT #mu ('+eta_title+', '+thr_title+', '+qTitle+') #phi']
                    #binnings['emtf_ugmt_muon'+eta_min_str+eta_max_str+thr_str+iso_wp_str+qMin_str+'_phi'] = [35, -3.5, 3.5, 'EMTF uGMT #mu ('+eta_title+', '+thr_title+', '+qTitle+') #phi']

                #for qual in range(16):
                #    qual_str = '_q'+str(qual)
                #    qualTitle = 'q'+str(qual)

                #    #varnames.append('gmt_muon'+eta_min_str+eta_max_str+thr_str+iso_wp_str+qual_str+'_phi')
                #    varnames.append('ugmt_muon'+eta_min_str+eta_max_str+thr_str+iso_wp_str+qual_str+'_phi')
                #    #varnames.append('bmtf_ugmt_muon'+eta_min_str+eta_max_str+thr_str+iso_wp_str+qual_str+'_phi')
                #    #varnames.append('omtf_ugmt_muon'+eta_min_str+eta_max_str+thr_str+iso_wp_str+qual_str+'_phi')
                #    #varnames.append('emtf_ugmt_muon'+eta_min_str+eta_max_str+thr_str+iso_wp_str+qual_str+'_phi')
                #    #binnings['gmt_muon'+eta_min_str+eta_max_str+thr_str+iso_wp_str+qual_str+'_phi'] = [35, -3.5, 3.5, 'GMT #mu ('+eta_title+', '+thr_title+', '+qualTitle+') #phi']
                #    binnings['ugmt_muon'+eta_min_str+eta_max_str+thr_str+iso_wp_str+qual_str+'_phi'] = [35, -3.5, 3.5, 'uGMT #mu ('+eta_title+', '+thr_title+', '+qualTitle+') #phi']
                #    #binnings['bmtf_ugmt_muon'+eta_min_str+eta_max_str+thr_str+iso_wp_str+qual_str+'_phi'] = [35, -3.5, 3.5, 'BMTF uGMT #mu ('+eta_title+', '+thr_title+', '+qualTitle+') #phi']
                #    #binnings['omtf_ugmt_muon'+eta_min_str+eta_max_str+thr_str+iso_wp_str+qual_str+'_phi'] = [35, -3.5, 3.5, 'OMTF uGMT #mu ('+eta_title+', '+thr_title+', '+qualTitle+') #phi']
                #    #binnings['emtf_ugmt_muon'+eta_min_str+eta_max_str+thr_str+iso_wp_str+qual_str+'_phi'] = [35, -3.5, 3.5, 'EMTF uGMT #mu ('+eta_title+', '+thr_title+', '+qualTitle+') #phi']

            for qMin in qualities:
                qMin_str = '_qmin'+str(qMin)
                qTitle = 'q>'+str(qMin)

                #varnames.append('gmt_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_pt')
                varnames.append('ugmt_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_pt')
                #varnames.append('bmtf_ugmt_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_pt')
                #varnames.append('omtf_ugmt_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_pt')
                #varnames.append('emtf_ugmt_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_pt')
                #binnings['gmt_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_pt'] = [100, 0, 100, '('+eta_title+', '+qTitle+') p_{T}', 'GeV/c']
                binnings['ugmt_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_pt'] = [100, 0, 100, '('+eta_title+', '+qTitle+') p_{T}', 'GeV/c']
                #binnings['bmtf_ugmt_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_pt'] = [100, 0, 100, '('+eta_title+', '+qTitle+') p_{T}', 'GeV/c']
                #binnings['omtf_ugmt_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_pt'] = [100, 0, 100, '('+eta_title+', '+qTitle+') p_{T}', 'GeV/c']
                #binnings['emtf_ugmt_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_pt'] = [100, 0, 100, '('+eta_title+', '+qTitle+') p_{T}', 'GeV/c']

                #varnames.append('gmt_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_varBin_pt')
                #varnames.append('ugmt_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_varBin_pt')
                #varnames.append('bmtf_ugmt_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_varBin_pt')
                #varnames.append('omtf_ugmt_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_varBin_pt')
                #varnames.append('emtf_ugmt_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_varBin_pt')
                #binnings['gmt_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_varBin_pt'] = pt_bins+['('+eta_title+', '+qTitle+') p_{T}', 'GeV/c']
                #binnings['ugmt_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_varBin_pt'] = pt_bins+['('+eta_title+', '+qTitle+') p_{T}', 'GeV/c']
                #binnings['bmtf_ugmt_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_varBin_pt'] = pt_bins+['('+eta_title+', '+qTitle+') p_{T}', 'GeV/c']
                #binnings['omtf_ugmt_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_varBin_pt'] = pt_bins+['('+eta_title+', '+qTitle+') p_{T}', 'GeV/c']
                #binnings['emtf_ugmt_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_varBin_pt'] = pt_bins+['('+eta_title+', '+qTitle+') p_{T}', 'GeV/c']

                #varnames.append('gmt_highest_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_pt')
                varnames.append('ugmt_highest_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_pt')
                #varnames.append('bmtf_ugmt_highest_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_pt')
                #varnames.append('omtf_ugmt_highest_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_pt')
                #varnames.append('emtf_ugmt_highest_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_pt')
                #varnames.append('bmtf_highest_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_pt')
                #varnames.append('omtf_highest_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_pt')
                #varnames.append('emtf_highest_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_pt')
                #binnings['gmt_highest_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_pt'] = [100, 0, 100, '('+eta_title+', '+qTitle+') p_{T}', 'GeV/c']
                binnings['ugmt_highest_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_pt'] = [100, 0, 100, '('+eta_title+', '+qTitle+') p_{T}', 'GeV/c']
                #binnings['bmtf_ugmt_highest_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_pt'] = [100, 0, 100, '('+eta_title+', '+qTitle+') p_{T}', 'GeV/c']
                #binnings['omtf_ugmt_highest_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_pt'] = [100, 0, 100, '('+eta_title+', '+qTitle+') p_{T}', 'GeV/c']
                #binnings['emtf_ugmt_highest_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_pt'] = [100, 0, 100, '('+eta_title+', '+qTitle+') p_{T}', 'GeV/c']
                #binnings['bmtf_highest_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_pt'] = [100, 0, 100, '('+eta_title+', '+qTitle+') p_{T}', 'GeV/c']
                #binnings['omtf_highest_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_pt'] = [100, 0, 100, '('+eta_title+', '+qTitle+') p_{T}', 'GeV/c']
                #binnings['emtf_highest_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_pt'] = [100, 0, 100, '('+eta_title+', '+qTitle+') p_{T}', 'GeV/c']

                #varnames.append('gmt_highest_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_varBin_pt')
                #varnames.append('ugmt_highest_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_varBin_pt')
                #varnames.append('bmtf_ugmt_highest_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_varBin_pt')
                #varnames.append('omtf_ugmt_highest_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_varBin_pt')
                #varnames.append('emtf_ugmt_highest_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_varBin_pt')
                #varnames.append('bmtf_highest_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_varBin_pt')
                #varnames.append('omtf_highest_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_varBin_pt')
                #varnames.append('emtf_highest_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_varBin_pt')
                #binnings['gmt_highest_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_varBin_pt'] = pt_bins+['('+eta_title+', '+qTitle+') p_{T}', 'GeV/c']
                #binnings['ugmt_highest_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_varBin_pt'] = pt_bins+['('+eta_title+', '+qTitle+') p_{T}', 'GeV/c']
                #binnings['bmtf_ugmt_highest_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_varBin_pt'] = pt_bins+['('+eta_title+', '+qTitle+') p_{T}', 'GeV/c']
                #binnings['omtf_ugmt_highest_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_varBin_pt'] = pt_bins+['('+eta_title+', '+qTitle+') p_{T}', 'GeV/c']
                #binnings['emtf_ugmt_highest_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_varBin_pt'] = pt_bins+['('+eta_title+', '+qTitle+') p_{T}', 'GeV/c']
                #binnings['bmtf_highest_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_varBin_pt'] = pt_bins+['('+eta_title+', '+qTitle+') p_{T}', 'GeV/c']
                #binnings['omtf_highest_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_varBin_pt'] = pt_bins+['('+eta_title+', '+qTitle+') p_{T}', 'GeV/c']
                #binnings['emtf_highest_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_varBin_pt'] = pt_bins+['('+eta_title+', '+qTitle+') p_{T}', 'GeV/c']

            #for qual in range(16):
            #    qual_str = '_q'+str(qual)
            #    qualTitle = 'q'+str(qual)

            #    #varnames.append('gmt_muon'+eta_min_str+eta_max_str+qual_str+iso_wp_str+'_pt')
            #    varnames.append('ugmt_muon'+eta_min_str+eta_max_str+qual_str+iso_wp_str+'_pt')
            #    #varnames.append('bmtf_ugmt_muon'+eta_min_str+eta_max_str+qual_str+iso_wp_str+'_pt')
            #    #varnames.append('omtf_ugmt_muon'+eta_min_str+eta_max_str+qual_str+iso_wp_str+'_pt')
            #    #varnames.append('emtf_ugmt_muon'+eta_min_str+eta_max_str+qual_str+iso_wp_str+'_pt')
            #    #binnings['gmt_muon'+eta_min_str+eta_max_str+qual_str+iso_wp_str+'_pt'] = [100, 0, 100, '('+eta_title+', '+qualTitle+') p_{T}', 'GeV/c']
            #    binnings['ugmt_muon'+eta_min_str+eta_max_str+qual_str+iso_wp_str+'_pt'] = [100, 0, 100, '('+eta_title+', '+qualTitle+') p_{T}', 'GeV/c']
            #    #binnings['bmtf_ugmt_muon'+eta_min_str+eta_max_str+qual_str+iso_wp_str+'_pt'] = [100, 0, 100, '('+eta_title+', '+qualTitle+') p_{T}', 'GeV/c']
            #    #binnings['omtf_ugmt_muon'+eta_min_str+eta_max_str+qual_str+iso_wp_str+'_pt'] = [100, 0, 100, '('+eta_title+', '+qualTitle+') p_{T}', 'GeV/c']
            #    #binnings['emtf_ugmt_muon'+eta_min_str+eta_max_str+qual_str+iso_wp_str+'_pt'] = [100, 0, 100, '('+eta_title+', '+qualTitle+') p_{T}', 'GeV/c']

            #    #varnames.append('gmt_muon'+eta_min_str+eta_max_str+qual_str+iso_wp_str+'_varBin_pt')
            #    #varnames.append('ugmt_muon'+eta_min_str+eta_max_str+qual_str+iso_wp_str+'_varBin_pt')
            #    #varnames.append('bmtf_ugmt_muon'+eta_min_str+eta_max_str+qual_str+iso_wp_str+'_varBin_pt')
            #    #varnames.append('omtf_ugmt_muon'+eta_min_str+eta_max_str+qual_str+iso_wp_str+'_varBin_pt')
            #    #varnames.append('emtf_ugmt_muon'+eta_min_str+eta_max_str+qual_str+iso_wp_str+'_varBin_pt')
            #    #binnings['gmt_muon'+eta_min_str+eta_max_str+qual_str+iso_wp_str+'_varBin_pt'] = pt_bins+['('+eta_title+', '+qualTitle+') p_{T}', 'GeV/c']
            #    #binnings['ugmt_muon'+eta_min_str+eta_max_str+qual_str+iso_wp_str+'_varBin_pt'] = pt_bins+['('+eta_title+', '+qualTitle+') p_{T}', 'GeV/c']
            #    #binnings['bmtf_ugmt_muon'+eta_min_str+eta_max_str+qual_str+iso_wp_str+'_varBin_pt'] = pt_bins+['('+eta_title+', '+qualTitle+') p_{T}', 'GeV/c']
            #    #binnings['omtf_ugmt_muon'+eta_min_str+eta_max_str+qual_str+iso_wp_str+'_varBin_pt'] = pt_bins+['('+eta_title+', '+qualTitle+') p_{T}', 'GeV/c']
            #    #binnings['emtf_ugmt_muon'+eta_min_str+eta_max_str+qual_str+iso_wp_str+'_varBin_pt'] = pt_bins+['('+eta_title+', '+qualTitle+') p_{T}', 'GeV/c']

        ##########################################################################
        #for threshold in thresholds:
        #    thr_str = '_ptmin'+str(threshold)
        #    thr_title = 'p_{T} > '+str(threshold)+' GeV/c'

        #    for qMin in qualities:
        #        qMin_str = '_qmin'+str(qMin)
        #        qTitle = 'q>'+str(qMin)

        #        #varnames.append('gmt_muon'+thr_str+iso_wp_str+qMin_str+'_eta')
        #        varnames.append('ugmt_muon'+thr_str+iso_wp_str+qMin_str+'_eta')
        #        #varnames.append('bmtf_ugmt_muon'+thr_str+iso_wp_str+qMin_str+'_eta')
        #        #varnames.append('omtf_ugmt_muon'+thr_str+iso_wp_str+qMin_str+'_eta')
        #        #varnames.append('emtf_ugmt_muon'+thr_str+iso_wp_str+qMin_str+'_eta')
        #        #binnings['gmt_muon'+thr_str+iso_wp_str+qMin_str+'_eta'] = eta_bins+['GMT #mu ('+thr_title+', '+qTitle+') #eta']
        #        binnings['ugmt_muon'+thr_str+iso_wp_str+qMin_str+'_eta'] = eta_bins+['uGMT #mu ('+thr_title+', '+qTitle+') #eta']
        #        #binnings['bmtf_ugmt_muon'+thr_str+iso_wp_str+qMin_str+'_eta'] = eta_bins+['BMTF uGMT #mu ('+thr_title+', '+qTitle+') #eta']
        #        #binnings['omtf_ugmt_muon'+thr_str+iso_wp_str+qMin_str+'_eta'] = eta_bins+['OMTF uGMT #mu ('+thr_title+', '+qTitle+') #eta']
        #        #binnings['emtf_ugmt_muon'+thr_str+iso_wp_str+qMin_str+'_eta'] = eta_bins+['EMTF uGMT #mu ('+thr_title+', '+qTitle+') #eta']

        #    #for qual in range(16):
        #    #    qual_str = '_q'+str(qual)
        #    #    qualTitle = 'q'+str(qual)

        #    #    #varnames.append('gmt_muon'+thr_str+iso_wp_str+qual_str+'_eta')
        #    #    varnames.append('ugmt_muon'+thr_str+iso_wp_str+qual_str+'_eta')
        #    #    #varnames.append('bmtf_ugmt_muon'+thr_str+iso_wp_str+qual_str+'_eta')
        #    #    #varnames.append('omtf_ugmt_muon'+thr_str+iso_wp_str+qual_str+'_eta')
        #    #    #varnames.append('emtf_ugmt_muon'+thr_str+iso_wp_str+qual_str+'_eta')
        #    #    #binnings['gmt_muon'+thr_str+iso_wp_str+qual_str+'_eta'] = eta_bins+['GMT #mu ('+thr_title+', '+qualTitle+') #eta']
        #    #    binnings['ugmt_muon'+thr_str+iso_wp_str+qual_str+'_eta'] = eta_bins+['uGMT #mu ('+thr_title+', '+qualTitle+') #eta']
        #    #    #binnings['bmtf_ugmt_muon'+thr_str+iso_wp_str+qual_str+'_eta'] = eta_bins+['BMTF uGMT #mu ('+thr_title+', '+qualTitle+') #eta']
        #    #    #binnings['omtf_ugmt_muon'+thr_str+iso_wp_str+qual_str+'_eta'] = eta_bins+['OMTF uGMT #mu ('+thr_title+', '+qualTitle+') #eta']
        #    #    #binnings['emtf_ugmt_muon'+thr_str+iso_wp_str+qual_str+'_eta'] = eta_bins+['EMTF uGMT #mu ('+thr_title+', '+qualTitle+') #eta']

    return HistManager(list(set(varnames)), binnings)

def get_highest_pt(candColl, idcs, gmt=False, tf=False):
    ptList = []
    for i in idcs:
        if not gmt:
            if not tf:
                ptList.append(candColl.muonEt[i])
            else:
                ptList.append(candColl.tfMuonHwPt[i] * ptScale)
        else:
            ptList.append(candColl.Pt[i])
    ptList.sort()
    return ptList.pop()


def get_highest_pt_idx(candColl, idcs, gmt=False, tf=False):
    highestPt = 0.
    highestPtIdx = -1
    for i in idcs:
        if not gmt:
            if not tf:
                pt = candColl.muonEt[i]
            else:
                pt = candColl.tfMuonHwPt[i] * ptScale
        else:
            pt = candColl.Pt[i]
        if pt > highestPt:
            highestPtIdx = i
            highestPt = pt
    return highestPtIdx


def analyse(evt, hm, eta_ranges, thresholds, qualities, iso_wps, emul):
    # USER HOOK
    # do what you want to do with the ntuples here
    # example:
    #print "GMT:", evt.gmt.N
    #print "UGMT:", evt.upgrade.n

    #if evt.upgrade.nMuons < 2:
    #    return

    bx_min = 0
    bx_max = 0

    # decide which ntuples to use
    if emul:
        ugmtColl = evt.upgradeEmu
        l1CaloTowerColl = evt.caloTowersEmu
    else:
        ugmtColl = evt.upgrade
        l1CaloTowerColl = evt.caloTowers

    #gmt_muon_idcs = MuonSelections.select_gmt_muons(evt.gmt, pt_min=0.5, bx_min=bx_min, bx_max=bx_max, pos_eta=pos_eta, neg_eta=neg_eta)
    gmt_muon_idcs = [] # don't fill GMT histograms
    ugmt_muon_idcs = MuonSelections.select_ugmt_muons(ugmtColl, pt_min=0.5, bx_min=bx_min, bx_max=bx_max, pos_eta=pos_eta, neg_eta=neg_eta, useVtxExtraCoord=useVtxExtraCoord)

    for iso_wp in iso_wps:
        iso_wp_str = '_isoMax{iso:.3f}'.format(iso=iso_wp)
        # remove non-isolated muons
        if iso_type == 2:
            iso_min=iso_wp
            iso_max=999
        else:
            iso_min=0.
            iso_max=iso_wp
        ugmt_iso_muon_idcs = MuonSelections.select_iso_ugmt_muons(ugmtColl, l1CaloTowerColl, iso_min=iso_min, iso_max=iso_max, iso_eta_max=isoEtaMax, idcs=ugmt_muon_idcs, useVtxExtraCoord=useVtxExtraCoord, iso_type=iso_type)

        for eta_range in eta_ranges:
            eta_min = eta_range[0]
            eta_max = eta_range[1]
            eta_min_str = '_absEtaMin'+str(eta_min)
            eta_max_str = '_absEtaMax'+str(eta_max)

            eta_gmt_muon_idcs = MuonSelections.select_gmt_muons(evt.gmt, abs_eta_min=eta_min, abs_eta_max=eta_max, idcs=gmt_muon_idcs)
            eta_ugmt_iso_muon_idcs = MuonSelections.select_ugmt_muons(ugmtColl, abs_eta_min=eta_min, abs_eta_max=eta_max, idcs=ugmt_iso_muon_idcs, useVtxExtraCoord=useVtxExtraCoord)

            for threshold in thresholds:
                thr_str = '_ptmin'+str(threshold)

                eta_thr_gmt_muon_idcs = MuonSelections.select_gmt_muons(evt.gmt, pt_min=threshold, idcs=eta_gmt_muon_idcs)
                eta_thr_ugmt_iso_muon_idcs = MuonSelections.select_ugmt_muons(ugmtColl, pt_min=threshold, idcs=eta_ugmt_iso_muon_idcs, useVtxExtraCoord=useVtxExtraCoord)

                #for i in eta_thr_gmt_muon_idcs:
                #    hm.fill('gmt_muon'+eta_min_str+eta_max_str+thr_str+iso_wp_str+'_qual', evt.gmt.Qual[i])
                #for i in eta_thr_ugmt_iso_muon_idcs:
                #    hm.fill('ugmt_muon'+eta_min_str+eta_max_str+thr_str+iso_wp_str+'_qual', ugmtColl.muonQual[i])
                #    #tftype = MuonSelections.getTfTypeFromTfMuonIdx(ugmtColl.muonTfMuonIdx[i])
                #    #if tftype is 0:
                #    #    hm.fill('bmtf_ugmt_muon'+eta_min_str+eta_max_str+thr_str+iso_wp_str+'_qual', ugmtColl.muonQual[i])
                #    #elif tftype is 1:
                #    #    hm.fill('omtf_ugmt_muon'+eta_min_str+eta_max_str+thr_str+iso_wp_str+'_qual', ugmtColl.muonQual[i])
                #    #elif tftype is 2:
                #    #    hm.fill('emtf_ugmt_muon'+eta_min_str+eta_max_str+thr_str+iso_wp_str+'_qual', ugmtColl.muonQual[i])

                for qMin in qualities:
                    qMin_str = '_qmin'+str(qMin)

                    eta_thr_q_gmt_muon_idcs = MuonSelections.select_gmt_muons(evt.gmt, qual_min=qMin, idcs=eta_thr_gmt_muon_idcs)
                    eta_thr_q_ugmt_iso_muon_idcs = MuonSelections.select_ugmt_muons(ugmtColl, qual_min=qMin, idcs=eta_thr_ugmt_iso_muon_idcs, useVtxExtraCoord=useVtxExtraCoord)

                    bmtfUgmtCtr = 0
                    omtfUgmtCtr = 0
                    emtfUgmtCtr = 0

                    #for i in eta_thr_q_gmt_muon_idcs:
                    #    hm.fill('gmt_muon'+eta_min_str+eta_max_str+thr_str+iso_wp_str+qMin_str+'_phi', Matcher.norm_phi(evt.gmt.Phi[i]))
                    #for i in eta_thr_q_ugmt_iso_muon_idcs:
                    #    if useVtxExtraCoord:
                    #        phi = ugmtColl.muonPhiAtVtx[i]
                    #    else:
                    #        phi = ugmtColl.muonPhi[i]

                    #    hm.fill('ugmt_muon'+eta_min_str+eta_max_str+thr_str+iso_wp_str+qMin_str+'_phi', phi)
                    #    #tftype = MuonSelections.getTfTypeFromTfMuonIdx(ugmtColl.muonTfMuonIdx[i])
                    #    #if tftype is 0:
                    #    #    hm.fill('bmtf_ugmt_muon'+eta_min_str+eta_max_str+thr_str+iso_wp_str+qMin_str+'_phi', phi)
                    #    #    bmtfUgmtCtr += 1
                    #    #elif tftype is 1:
                    #    #    hm.fill('omtf_ugmt_muon'+eta_min_str+eta_max_str+thr_str+iso_wp_str+qMin_str+'_phi', phi)
                    #    #    omtfUgmtCtr += 1
                    #    #elif tftype is 2:
                    #    #    hm.fill('emtf_ugmt_muon'+eta_min_str+eta_max_str+thr_str+iso_wp_str+qMin_str+'_phi', phi)
                    #    #    emtfUgmtCtr += 1

                    #hm.fill('n_gmt_muons'+eta_min_str+eta_max_str+thr_str+iso_wp_str+qMin_str, len(eta_thr_q_gmt_muon_idcs))
                    hm.fill('n_ugmt_muons'+eta_min_str+eta_max_str+thr_str+iso_wp_str+qMin_str, len(eta_thr_q_ugmt_iso_muon_idcs))
                    #hm.fill('n_bmtf_ugmt_muons'+eta_min_str+eta_max_str+thr_str+iso_wp_str+qMin_str, bmtfUgmtCtr)
                    #hm.fill('n_omtf_ugmt_muons'+eta_min_str+eta_max_str+thr_str+iso_wp_str+qMin_str, omtfUgmtCtr)
                    #hm.fill('n_emtf_ugmt_muons'+eta_min_str+eta_max_str+thr_str+iso_wp_str+qMin_str, emtfUgmtCtr)

                #for qual in range(16):
                #    qual_str = '_q'+str(qual)

                #    #for i in eta_thr_gmt_muon_idcs:
                #    #    if evt.gmt.Qual[i] == qual:
                #    #        hm.fill('gmt_muon'+eta_min_str+eta_max_str+thr_str+iso_wp_str+qual_str+'_phi', Matcher.norm_phi(evt.gmt.Phi[i]))
                #    for i in eta_thr_ugmt_iso_muon_idcs:
                #        if ugmtColl.muonQual[i] == qual:
                #            if useVtxExtraCoord:
                #                phi = ugmtColl.muonPhiAtVtx[i]
                #            else:
                #                phi = ugmtColl.muonPhi[i]

                #            hm.fill('ugmt_muon'+eta_min_str+eta_max_str+thr_str+iso_wp_str+qual_str+'_phi', phi)
                #            #tftype = MuonSelections.getTfTypeFromTfMuonIdx(ugmtColl.muonTfMuonIdx[i])
                #            #if tftype is 0:
                #            #    hm.fill('bmtf_ugmt_muon'+eta_min_str+eta_max_str+thr_str+iso_wp_str+qual_str+'_phi', phi)
                #            #elif tftype is 1:
                #            #    hm.fill('omtf_ugmt_muon'+eta_min_str+eta_max_str+thr_str+iso_wp_str+qual_str+'_phi', phi)
                #            #elif tftype is 2:
                #            #    hm.fill('emtf_ugmt_muon'+eta_min_str+eta_max_str+thr_str+iso_wp_str+qual_str+'_phi', phi)

            for qMin in qualities:
                qMin_str = '_qmin'+str(qMin)

                eta_q_gmt_muon_idcs = MuonSelections.select_gmt_muons(evt.gmt, qual_min=qMin, idcs=eta_gmt_muon_idcs)
                eta_q_ugmt_iso_muon_idcs = MuonSelections.select_ugmt_muons(ugmtColl, qual_min=qMin, idcs=eta_ugmt_iso_muon_idcs, useVtxExtraCoord=useVtxExtraCoord)

                #for i in eta_q_gmt_muon_idcs:
                #    hm.fill('gmt_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_pt', evt.gmt.Pt[i])
                #    hm.fill('gmt_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_varBin_pt', evt.gmt.Pt[i])
                for i in eta_q_ugmt_iso_muon_idcs:
                    hm.fill('ugmt_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_pt', ugmtColl.muonEt[i])
                    #hm.fill('ugmt_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_varBin_pt', ugmtColl.muonEt[i])
                    #tftype = MuonSelections.getTfTypeFromTfMuonIdx(ugmtColl.muonTfMuonIdx[i])
                    #if tftype is 0:
                    #    hm.fill('bmtf_ugmt_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_pt', ugmtColl.muonEt[i])
                    #    hm.fill('bmtf_ugmt_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_varBin_pt', ugmtColl.muonEt[i])
                    #elif tftype is 1:
                    #    hm.fill('omtf_ugmt_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_pt', ugmtColl.muonEt[i])
                    #    hm.fill('omtf_ugmt_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_varBin_pt', ugmtColl.muonEt[i])
                    #elif tftype is 2:
                    #    hm.fill('emtf_ugmt_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_pt', ugmtColl.muonEt[i])
                    #    hm.fill('emtf_ugmt_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_varBin_pt', ugmtColl.muonEt[i])

                #If len(eta_q_gmt_muon_idcs):
                #    highestPt = get_highest_pt(evt.gmt, eta_q_gmt_muon_idcs, gmt=True)
                #    hm.fill('gmt_highest_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_pt', highestPt)
                #    hm.fill('gmt_highest_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_varBin_pt', highestPt)
                if len(eta_q_ugmt_iso_muon_idcs):
                    highestPtIdx = get_highest_pt_idx(ugmtColl, eta_q_ugmt_iso_muon_idcs)
                    highestPt = ugmtColl.muonEt[highestPtIdx]
                    hm.fill('ugmt_highest_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_pt', highestPt)
                    #hm.fill('ugmt_highest_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_varBin_pt', highestPt)
                    #tftype = MuonSelections.getTfTypeFromTfMuonIdx(ugmtColl.muonTfMuonIdx[highestPtIdx])
                    #if tftype is 0:
                    #    hm.fill('bmtf_ugmt_highest_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_pt', highestPt)
                    #    hm.fill('bmtf_ugmt_highest_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_varBin_pt', highestPt)
                    #elif tftype is 1:
                    #    hm.fill('omtf_ugmt_highest_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_pt', highestPt)
                    #    hm.fill('omtf_ugmt_highest_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_varBin_pt', highestPt)
                    #elif tftype is 2:
                    #    hm.fill('emtf_ugmt_highest_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_pt', highestPt)
                    #    hm.fill('emtf_ugmt_highest_muon'+eta_min_str+eta_max_str+qMin_str+iso_wp_str+'_varBin_pt', highestPt)

            #for qual in range(16):
            #    qual_str = '_q'+str(qual)

            #    #for i in eta_gmt_muon_idcs:
            #    #    if evt.gmt.Qual[i] == qual:
            #    #        hm.fill('gmt_muon'+eta_min_str+eta_max_str+qual_str+iso_wp_str+'_pt', evt.gmt.Pt[i])
            #    #        hm.fill('gmt_muon'+eta_min_str+eta_max_str+qual_str+iso_wp_str+'_varBin_pt', evt.gmt.Pt[i])
            #    for i in eta_ugmt_iso_muon_idcs:
            #        if ugmtColl.muonQual[i] == qual:
            #            hm.fill('ugmt_muon'+eta_min_str+eta_max_str+qual_str+iso_wp_str+'_pt', ugmtColl.muonEt[i])
            #            #hm.fill('ugmt_muon'+eta_min_str+eta_max_str+qual_str+iso_wp_str+'_varBin_pt', ugmtColl.muonEt[i])
            #            #tftype = MuonSelections.getTfTypeFromTfMuonIdx(ugmtColl.muonTfMuonIdx[i])
            #            #if tftype is 0:
            #            #    hm.fill('bmtf_ugmt_muon'+eta_min_str+eta_max_str+qual_str+iso_wp_str+'_pt', ugmtColl.muonEt[i])
            #            #    hm.fill('bmtf_ugmt_muon'+eta_min_str+eta_max_str+qual_str+iso_wp_str+'_varBin_pt', ugmtColl.muonEt[i])
            #            #elif tftype is 1:
            #            #    hm.fill('omtf_ugmt_muon'+eta_min_str+eta_max_str+qual_str+iso_wp_str+'_pt', ugmtColl.muonEt[i])
            #            #    hm.fill('omtf_ugmt_muon'+eta_min_str+eta_max_str+qual_str+iso_wp_str+'_varBin_pt', ugmtColl.muonEt[i])
            #            #elif tftype is 2:
            #            #    hm.fill('emtf_ugmt_muon'+eta_min_str+eta_max_str+qual_str+iso_wp_str+'_pt', ugmtColl.muonEt[i])
            #            #    hm.fill('emtf_ugmt_muon'+eta_min_str+eta_max_str+qual_str+iso_wp_str+'_varBin_pt', ugmtColl.muonEt[i])

        #for threshold in thresholds:
        #    thr_str = '_ptmin'+str(threshold)

        #    thr_gmt_muon_idcs = MuonSelections.select_gmt_muons(evt.gmt, pt_min=threshold, idcs=gmt_muon_idcs)
        #    thr_ugmt_iso_muon_idcs = MuonSelections.select_ugmt_muons(ugmtColl, pt_min=threshold, idcs=ugmt_iso_muon_idcs, useVtxExtraCoord=useVtxExtraCoord)

        #    for qMin in qualities:
        #        qMin_str = '_qmin'+str(qMin)

        #        thr_q_gmt_muon_idcs = MuonSelections.select_gmt_muons(evt.gmt, pt_min=threshold, qual_min=qMin, idcs=gmt_muon_idcs)
        #        thr_q_ugmt_iso_muon_idcs = MuonSelections.select_ugmt_muons(ugmtColl, pt_min=threshold, qual_min=qMin, idcs=ugmt_iso_muon_idcs, useVtxExtraCoord=useVtxExtraCoord)
        #    
        #        #for i in thr_q_gmt_muon_idcs:
        #        #    hm.fill('gmt_muon'+thr_str+iso_wp_str+qMin_str+'_eta', evt.gmt.Eta[i])
        #        for i in thr_q_ugmt_iso_muon_idcs:
        #            if useVtxExtraCoord:
        #                eta = ugmtColl.muonEtaAtVtx[i]
        #            else:
        #                eta = ugmtColl.muonEta[i]

        #            hm.fill('ugmt_muon'+thr_str+iso_wp_str+qMin_str+'_eta', eta)
        #            #tftype = MuonSelections.getTfTypeFromTfMuonIdx(ugmtColl.muonTfMuonIdx[i])
        #            #if tftype is 0:
        #            #    hm.fill('bmtf_ugmt_muon'+thr_str+iso_wp_str+qMin_str+'_eta', eta)
        #            #elif tftype is 1:
        #            #    hm.fill('omtf_ugmt_muon'+thr_str+iso_wp_str+qMin_str+'_eta', eta)
        #            #elif tftype is 2:
        #            #    hm.fill('emtf_ugmt_muon'+thr_str+iso_wp_str+qMin_str+'_eta', eta)

        #    #for qual in range(16):
        #    #    qual_str = '_q'+str(qual)

        #    #    #for i in thr_gmt_muon_idcs:
        #    #    #    if evt.gmt.Qual[i] == qual:
        #    #    #        hm.fill('gmt_muon'+thr_str+iso_wp_str+qual_str+'_eta', evt.gmt.Eta[i])
        #    #    for i in thr_ugmt_iso_muon_idcs:
        #    #        if ugmtColl.muonQual[i] == qual:
        #    #            if useVtxExtraCoord:
        #    #                eta = ugmtColl.muonEtaAtVtx[i]
        #    #            else:
        #    #                eta = ugmtColl.muonEta[i]

        #    #            hm.fill('ugmt_muon'+thr_str+iso_wp_str+qual_str+'_eta', eta)
        #    #            #tftype = MuonSelections.getTfTypeFromTfMuonIdx(ugmtColl.muonTfMuonIdx[i])
        #    #            #if tftype is 0:
        #    #            #    hm.fill('bmtf_ugmt_muon'+thr_str+iso_wp_str+qual_str+'_eta', eta)
        #    #            #elif tftype is 1:
        #    #            #    hm.fill('omtf_ugmt_muon'+thr_str+iso_wp_str+qual_str+'_eta', eta)
        #    #            #elif tftype is 2:
        #    #            #    hm.fill('emtf_ugmt_muon'+thr_str+iso_wp_str+qual_str+'_eta', eta)

def save_histos(hm, outfile):
    '''
    save all histograms in hm to outfile
    '''
    outfile.cd()
    for varname in hm.get_varnames():
        hm.get(varname).Write()

def main():
    L1Ana.init_l1_analysis()
    opts = parse_options_upgradeRateHistos(parser)
    emulated = opts.emul
    print ""

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
    else:
        iso_type = 0

    global isoEtaMax
    isoEtaMax = opts.etarestricted

    #eta_ranges = [[0, 2.5], [0, 2.1], [0, 0.83], [0.83, 1.24], [1.24, 2.5], [1.24, 2.1]]
    #thresholds = [1, 5, 10, 12, 16, 20, 24, 30]
    #qualities = range(16)
    #eta_ranges = [[0, 2.5], [0, 2.1], [0, 0.83], [0.83, 1.24], [1.24, 2.5]]
    #thresholds = [0, 3, 5, 7, 12, 18, 22]
    #qualities = [0, 4, 8, 12]
    if iso_type == 0: # absolute isolation
        iso_wps = [0, 1, 3, 5, 7, 9, 11, 15, 20, 25, 28, 30, 31]
    elif iso_type == 1: # relative isolation
        iso_wps = [0., 1/2., 1/3., 2/3., 3/4., 4/5., 5/6., 6/7., 1., 1.5, 3., 5., 10., 31., 62.]
    elif iso_type == 3 or iso_type == 5: # outer cone over total cone
        iso_wps = [0., 1/1., 1/2., 1/3., 2/3., 3/4., 4/5., 5/6., 6/7., 7/8., 8/9., 9/10., 19/20., 30/31.]
    else: # inner cone
        iso_wps = [0., 1., 2., 3., 4., 5., 6., 7., 8., 9.]

    eta_ranges = [[0, 2.5]]
    thresholds = [0]
    qualities = [4, 8, 12]

    # book the histograms
    hm = book_histograms(eta_ranges, thresholds, qualities, iso_wps)

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

    start_evt = opts.start_event
    end_evt = opts.start_event+ntuple.nevents
    analysed_evt_ctr = 0
    try:
        for i in range(start_evt, end_evt):
            event = ntuple[i]
            if (i+1) % 1000 == 0:
                L1Ana.log.info("Processing event: {n}. Analysed events from selected runs/LS until now: {nAna}".format(n=i+1, nAna=analysed_evt_ctr))

            runnr = event.event.run
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

            analyse(event, hm, eta_ranges, thresholds, qualities, iso_wps, emulated)
            hm.fill('n_evts_analysed', 0.5)
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
    useVtxExtraCoord = False
    iso_type = 0
    isoEtaMax = 3.
    saveHistos = True
    main()

