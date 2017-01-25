#!/usr/bin/env python
from ToolBox import parse_options_and_init_log
# have to do this first or ROOT masks the -h messages
opts, parser = parse_options_and_init_log()

from L1Analysis import L1Ana, L1Ntuple
from analysis_tools.plotting import HistManager
from analysis_tools.selections import MuonSelections, Matcher
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

    opts, unknown = parser.parse_known_args()
    return opts

def book_histograms(eta_ranges, thresholds, qualities):
    varnames = []
    binnings = {}

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

    for eta_range in eta_ranges:
        eta_min = eta_range[0]
        eta_max = eta_range[1]
        eta_min_str = '_absEtaMin'+str(eta_min)
        eta_max_str = '_absEtaMax'+str(eta_max)
        eta_title = '{eMin} < |#eta| < {eMax}'.format(eMin=eta_min, eMax=eta_max)

        for threshold in thresholds:
            thr_str = '_ptmin'+str(threshold)
            thr_title = 'p_{T} > '+str(threshold)+' GeV/c'

            varnames.append('gmt_muon'+eta_min_str+eta_max_str+thr_str+'_qual')
            varnames.append('ugmt_muon'+eta_min_str+eta_max_str+thr_str+'_qual')
            varnames.append('bmtf_ugmt_muon'+eta_min_str+eta_max_str+thr_str+'_qual')
            varnames.append('omtf_ugmt_muon'+eta_min_str+eta_max_str+thr_str+'_qual')
            varnames.append('emtf_ugmt_muon'+eta_min_str+eta_max_str+thr_str+'_qual')
            varnames.append('bmtf_muon'+eta_min_str+eta_max_str+thr_str+'_qual')
            varnames.append('omtf_muon'+eta_min_str+eta_max_str+thr_str+'_qual')
            varnames.append('emtf_muon'+eta_min_str+eta_max_str+thr_str+'_qual')
            binnings['gmt_muon'+eta_min_str+eta_max_str+thr_str+'_qual'] = [16, 0, 16, '('+eta_title+', '+thr_title+') quality']
            binnings['ugmt_muon'+eta_min_str+eta_max_str+thr_str+'_qual'] = [16, 0, 16, '('+eta_title+', '+thr_title+') quality']
            binnings['bmtf_ugmt_muon'+eta_min_str+eta_max_str+thr_str+'_qual'] = [16, 0, 16, '('+eta_title+', '+thr_title+') quality']
            binnings['omtf_ugmt_muon'+eta_min_str+eta_max_str+thr_str+'_qual'] = [16, 0, 16, '('+eta_title+', '+thr_title+') quality']
            binnings['emtf_ugmt_muon'+eta_min_str+eta_max_str+thr_str+'_qual'] = [16, 0, 16, '('+eta_title+', '+thr_title+') quality']
            binnings['bmtf_muon'+eta_min_str+eta_max_str+thr_str+'_qual'] = [16, 0, 16, '('+eta_title+', '+thr_title+') quality']
            binnings['omtf_muon'+eta_min_str+eta_max_str+thr_str+'_qual'] = [16, 0, 16, '('+eta_title+', '+thr_title+') quality']
            binnings['emtf_muon'+eta_min_str+eta_max_str+thr_str+'_qual'] = [16, 0, 16, '('+eta_title+', '+thr_title+') quality']

            for qMin in qualities:
                qMin_str = '_qmin'+str(qMin)
                qTitle = 'q>'+str(qMin)

                varnames.append('n_gmt_muons'+eta_min_str+eta_max_str+thr_str+qMin_str)
                varnames.append('n_ugmt_muons'+eta_min_str+eta_max_str+thr_str+qMin_str)
                varnames.append('n_bmtf_ugmt_muons'+eta_min_str+eta_max_str+thr_str+qMin_str)
                varnames.append('n_omtf_ugmt_muons'+eta_min_str+eta_max_str+thr_str+qMin_str)
                varnames.append('n_emtf_ugmt_muons'+eta_min_str+eta_max_str+thr_str+qMin_str)
                varnames.append('n_bmtf_muons'+eta_min_str+eta_max_str+thr_str+qMin_str)
                varnames.append('n_omtf_muons'+eta_min_str+eta_max_str+thr_str+qMin_str)
                varnames.append('n_emtf_muons'+eta_min_str+eta_max_str+thr_str+qMin_str)
                binnings['n_gmt_muons'+eta_min_str+eta_max_str+thr_str+qMin_str] = [10, 0, 10, '# GMT #mu '+eta_title+', '+thr_title+', '+qTitle]
                binnings['n_ugmt_muons'+eta_min_str+eta_max_str+thr_str+qMin_str] = [10, 0, 10, '# uGMT #mu '+eta_title+', '+thr_title+', '+qTitle]
                binnings['n_bmtf_ugmt_muons'+eta_min_str+eta_max_str+thr_str+qMin_str] = [10, 0, 10, '# BMTF uGMT #mu '+eta_title+', '+thr_title+', '+qTitle]
                binnings['n_omtf_ugmt_muons'+eta_min_str+eta_max_str+thr_str+qMin_str] = [10, 0, 10, '# OMTF uGMT #mu '+eta_title+', '+thr_title+', '+qTitle]
                binnings['n_emtf_ugmt_muons'+eta_min_str+eta_max_str+thr_str+qMin_str] = [10, 0, 10, '# EMTF uGMT #mu '+eta_title+', '+thr_title+', '+qTitle]
                binnings['n_bmtf_muons'+eta_min_str+eta_max_str+thr_str+qMin_str] = [10, 0, 10, '# BMTF #mu '+eta_title+', '+thr_title+', '+qTitle]
                binnings['n_omtf_muons'+eta_min_str+eta_max_str+thr_str+qMin_str] = [10, 0, 10, '# OMTF #mu '+eta_title+', '+thr_title+', '+qTitle]
                binnings['n_emtf_muons'+eta_min_str+eta_max_str+thr_str+qMin_str] = [10, 0, 10, '# EMTF #mu '+eta_title+', '+thr_title+', '+qTitle]

                varnames.append('gmt_muon'+eta_min_str+eta_max_str+thr_str+qMin_str+'_phi')
                varnames.append('ugmt_muon'+eta_min_str+eta_max_str+thr_str+qMin_str+'_phi')
                varnames.append('bmtf_ugmt_muon'+eta_min_str+eta_max_str+thr_str+qMin_str+'_phi')
                varnames.append('omtf_ugmt_muon'+eta_min_str+eta_max_str+thr_str+qMin_str+'_phi')
                varnames.append('emtf_ugmt_muon'+eta_min_str+eta_max_str+thr_str+qMin_str+'_phi')
                varnames.append('bmtf_muon'+eta_min_str+eta_max_str+thr_str+qMin_str+'_phi')
                varnames.append('omtf_muon'+eta_min_str+eta_max_str+thr_str+qMin_str+'_phi')
                varnames.append('emtf_muon'+eta_min_str+eta_max_str+thr_str+qMin_str+'_phi')
                binnings['gmt_muon'+eta_min_str+eta_max_str+thr_str+qMin_str+'_phi'] = [35, -3.5, 3.5, 'GMT #mu ('+eta_title+', '+thr_title+', '+qTitle+') #phi']
                binnings['ugmt_muon'+eta_min_str+eta_max_str+thr_str+qMin_str+'_phi'] = [35, -3.5, 3.5, 'uGMT #mu ('+eta_title+', '+thr_title+', '+qTitle+') #phi']
                binnings['bmtf_ugmt_muon'+eta_min_str+eta_max_str+thr_str+qMin_str+'_phi'] = [35, -3.5, 3.5, 'BMTF uGMT #mu ('+eta_title+', '+thr_title+', '+qTitle+') #phi']
                binnings['omtf_ugmt_muon'+eta_min_str+eta_max_str+thr_str+qMin_str+'_phi'] = [35, -3.5, 3.5, 'OMTF uGMT #mu ('+eta_title+', '+thr_title+', '+qTitle+') #phi']
                binnings['emtf_ugmt_muon'+eta_min_str+eta_max_str+thr_str+qMin_str+'_phi'] = [35, -3.5, 3.5, 'EMTF uGMT #mu ('+eta_title+', '+thr_title+', '+qTitle+') #phi']
                binnings['bmtf_muon'+eta_min_str+eta_max_str+thr_str+qMin_str+'_phi'] = [35, -3.5, 3.5, 'BMTF #mu ('+eta_title+', '+thr_title+', '+qTitle+') #phi']
                binnings['omtf_muon'+eta_min_str+eta_max_str+thr_str+qMin_str+'_phi'] = [35, -3.5, 3.5, 'OMTF #mu ('+eta_title+', '+thr_title+', '+qTitle+') #phi']
                binnings['emtf_muon'+eta_min_str+eta_max_str+thr_str+qMin_str+'_phi'] = [35, -3.5, 3.5, 'EMTF #mu ('+eta_title+', '+thr_title+', '+qTitle+') #phi']

            for qual in range(16):
                qual_str = '_q'+str(qual)
                qualTitle = 'q'+str(qual)

                varnames.append('gmt_muon'+eta_min_str+eta_max_str+thr_str+qual_str+'_phi')
                varnames.append('ugmt_muon'+eta_min_str+eta_max_str+thr_str+qual_str+'_phi')
                varnames.append('bmtf_ugmt_muon'+eta_min_str+eta_max_str+thr_str+qual_str+'_phi')
                varnames.append('omtf_ugmt_muon'+eta_min_str+eta_max_str+thr_str+qual_str+'_phi')
                varnames.append('emtf_ugmt_muon'+eta_min_str+eta_max_str+thr_str+qual_str+'_phi')
                varnames.append('bmtf_muon'+eta_min_str+eta_max_str+thr_str+qual_str+'_phi')
                varnames.append('omtf_muon'+eta_min_str+eta_max_str+thr_str+qual_str+'_phi')
                varnames.append('emtf_muon'+eta_min_str+eta_max_str+thr_str+qual_str+'_phi')
                binnings['gmt_muon'+eta_min_str+eta_max_str+thr_str+qual_str+'_phi'] = [35, -3.5, 3.5, 'GMT #mu ('+eta_title+', '+thr_title+', '+qualTitle+') #phi']
                binnings['ugmt_muon'+eta_min_str+eta_max_str+thr_str+qual_str+'_phi'] = [35, -3.5, 3.5, 'uGMT #mu ('+eta_title+', '+thr_title+', '+qualTitle+') #phi']
                binnings['bmtf_ugmt_muon'+eta_min_str+eta_max_str+thr_str+qual_str+'_phi'] = [35, -3.5, 3.5, 'BMTF uGMT #mu ('+eta_title+', '+thr_title+', '+qualTitle+') #phi']
                binnings['omtf_ugmt_muon'+eta_min_str+eta_max_str+thr_str+qual_str+'_phi'] = [35, -3.5, 3.5, 'OMTF uGMT #mu ('+eta_title+', '+thr_title+', '+qualTitle+') #phi']
                binnings['emtf_ugmt_muon'+eta_min_str+eta_max_str+thr_str+qual_str+'_phi'] = [35, -3.5, 3.5, 'EMTF uGMT #mu ('+eta_title+', '+thr_title+', '+qualTitle+') #phi']
                binnings['bmtf_muon'+eta_min_str+eta_max_str+thr_str+qual_str+'_phi'] = [35, -3.5, 3.5, 'BMTF #mu ('+eta_title+', '+thr_title+', '+qualTitle+') #phi']
                binnings['omtf_muon'+eta_min_str+eta_max_str+thr_str+qual_str+'_phi'] = [35, -3.5, 3.5, 'OMTF #mu ('+eta_title+', '+thr_title+', '+qualTitle+') #phi']
                binnings['emtf_muon'+eta_min_str+eta_max_str+thr_str+qual_str+'_phi'] = [35, -3.5, 3.5, 'EMTF #mu ('+eta_title+', '+thr_title+', '+qualTitle+') #phi']

        for qMin in qualities:
            qMin_str = '_qmin'+str(qMin)
            qTitle = 'q>'+str(qMin)

            varnames.append('gmt_muon'+eta_min_str+eta_max_str+qMin_str+'_pt')
            varnames.append('ugmt_muon'+eta_min_str+eta_max_str+qMin_str+'_pt')
            varnames.append('bmtf_ugmt_muon'+eta_min_str+eta_max_str+qMin_str+'_pt')
            varnames.append('omtf_ugmt_muon'+eta_min_str+eta_max_str+qMin_str+'_pt')
            varnames.append('emtf_ugmt_muon'+eta_min_str+eta_max_str+qMin_str+'_pt')
            varnames.append('bmtf_muon'+eta_min_str+eta_max_str+qMin_str+'_pt')
            varnames.append('omtf_muon'+eta_min_str+eta_max_str+qMin_str+'_pt')
            varnames.append('emtf_muon'+eta_min_str+eta_max_str+qMin_str+'_pt')
            binnings['gmt_muon'+eta_min_str+eta_max_str+qMin_str+'_pt'] = [100, 0, 100, '('+eta_title+', '+qTitle+') p_{T}', 'GeV/c']
            binnings['ugmt_muon'+eta_min_str+eta_max_str+qMin_str+'_pt'] = [100, 0, 100, '('+eta_title+', '+qTitle+') p_{T}', 'GeV/c']
            binnings['bmtf_ugmt_muon'+eta_min_str+eta_max_str+qMin_str+'_pt'] = [100, 0, 100, '('+eta_title+', '+qTitle+') p_{T}', 'GeV/c']
            binnings['omtf_ugmt_muon'+eta_min_str+eta_max_str+qMin_str+'_pt'] = [100, 0, 100, '('+eta_title+', '+qTitle+') p_{T}', 'GeV/c']
            binnings['emtf_ugmt_muon'+eta_min_str+eta_max_str+qMin_str+'_pt'] = [100, 0, 100, '('+eta_title+', '+qTitle+') p_{T}', 'GeV/c']
            binnings['bmtf_muon'+eta_min_str+eta_max_str+qMin_str+'_pt'] = [100, 0, 100, '('+eta_title+', '+qTitle+') p_{T}', 'GeV/c']
            binnings['omtf_muon'+eta_min_str+eta_max_str+qMin_str+'_pt'] = [100, 0, 100, '('+eta_title+', '+qTitle+') p_{T}', 'GeV/c']
            binnings['emtf_muon'+eta_min_str+eta_max_str+qMin_str+'_pt'] = [100, 0, 100, '('+eta_title+', '+qTitle+') p_{T}', 'GeV/c']

            varnames.append('gmt_muon'+eta_min_str+eta_max_str+qMin_str+'_varBin_pt')
            varnames.append('ugmt_muon'+eta_min_str+eta_max_str+qMin_str+'_varBin_pt')
            varnames.append('bmtf_ugmt_muon'+eta_min_str+eta_max_str+qMin_str+'_varBin_pt')
            varnames.append('omtf_ugmt_muon'+eta_min_str+eta_max_str+qMin_str+'_varBin_pt')
            varnames.append('emtf_ugmt_muon'+eta_min_str+eta_max_str+qMin_str+'_varBin_pt')
            varnames.append('bmtf_muon'+eta_min_str+eta_max_str+qMin_str+'_varBin_pt')
            varnames.append('omtf_muon'+eta_min_str+eta_max_str+qMin_str+'_varBin_pt')
            varnames.append('emtf_muon'+eta_min_str+eta_max_str+qMin_str+'_varBin_pt')
            binnings['gmt_muon'+eta_min_str+eta_max_str+qMin_str+'_varBin_pt'] = pt_bins+['('+eta_title+', '+qTitle+') p_{T}', 'GeV/c']
            binnings['ugmt_muon'+eta_min_str+eta_max_str+qMin_str+'_varBin_pt'] = pt_bins+['('+eta_title+', '+qTitle+') p_{T}', 'GeV/c']
            binnings['bmtf_ugmt_muon'+eta_min_str+eta_max_str+qMin_str+'_varBin_pt'] = pt_bins+['('+eta_title+', '+qTitle+') p_{T}', 'GeV/c']
            binnings['omtf_ugmt_muon'+eta_min_str+eta_max_str+qMin_str+'_varBin_pt'] = pt_bins+['('+eta_title+', '+qTitle+') p_{T}', 'GeV/c']
            binnings['emtf_ugmt_muon'+eta_min_str+eta_max_str+qMin_str+'_varBin_pt'] = pt_bins+['('+eta_title+', '+qTitle+') p_{T}', 'GeV/c']
            binnings['bmtf_muon'+eta_min_str+eta_max_str+qMin_str+'_varBin_pt'] = pt_bins+['('+eta_title+', '+qTitle+') p_{T}', 'GeV/c']
            binnings['omtf_muon'+eta_min_str+eta_max_str+qMin_str+'_varBin_pt'] = pt_bins+['('+eta_title+', '+qTitle+') p_{T}', 'GeV/c']
            binnings['emtf_muon'+eta_min_str+eta_max_str+qMin_str+'_varBin_pt'] = pt_bins+['('+eta_title+', '+qTitle+') p_{T}', 'GeV/c']

            varnames.append('gmt_highest_muon'+eta_min_str+eta_max_str+qMin_str+'_pt')
            varnames.append('ugmt_highest_muon'+eta_min_str+eta_max_str+qMin_str+'_pt')
            varnames.append('bmtf_ugmt_highest_muon'+eta_min_str+eta_max_str+qMin_str+'_pt')
            varnames.append('omtf_ugmt_highest_muon'+eta_min_str+eta_max_str+qMin_str+'_pt')
            varnames.append('emtf_ugmt_highest_muon'+eta_min_str+eta_max_str+qMin_str+'_pt')
            varnames.append('bmtf_highest_muon'+eta_min_str+eta_max_str+qMin_str+'_pt')
            varnames.append('omtf_highest_muon'+eta_min_str+eta_max_str+qMin_str+'_pt')
            varnames.append('emtf_highest_muon'+eta_min_str+eta_max_str+qMin_str+'_pt')
            binnings['gmt_highest_muon'+eta_min_str+eta_max_str+qMin_str+'_pt'] = [100, 0, 100, '('+eta_title+', '+qTitle+') p_{T}', 'GeV/c']
            binnings['ugmt_highest_muon'+eta_min_str+eta_max_str+qMin_str+'_pt'] = [100, 0, 100, '('+eta_title+', '+qTitle+') p_{T}', 'GeV/c']
            binnings['bmtf_ugmt_highest_muon'+eta_min_str+eta_max_str+qMin_str+'_pt'] = [100, 0, 100, '('+eta_title+', '+qTitle+') p_{T}', 'GeV/c']
            binnings['omtf_ugmt_highest_muon'+eta_min_str+eta_max_str+qMin_str+'_pt'] = [100, 0, 100, '('+eta_title+', '+qTitle+') p_{T}', 'GeV/c']
            binnings['emtf_ugmt_highest_muon'+eta_min_str+eta_max_str+qMin_str+'_pt'] = [100, 0, 100, '('+eta_title+', '+qTitle+') p_{T}', 'GeV/c']
            binnings['bmtf_highest_muon'+eta_min_str+eta_max_str+qMin_str+'_pt'] = [100, 0, 100, '('+eta_title+', '+qTitle+') p_{T}', 'GeV/c']
            binnings['omtf_highest_muon'+eta_min_str+eta_max_str+qMin_str+'_pt'] = [100, 0, 100, '('+eta_title+', '+qTitle+') p_{T}', 'GeV/c']
            binnings['emtf_highest_muon'+eta_min_str+eta_max_str+qMin_str+'_pt'] = [100, 0, 100, '('+eta_title+', '+qTitle+') p_{T}', 'GeV/c']

            varnames.append('gmt_highest_muon'+eta_min_str+eta_max_str+qMin_str+'_varBin_pt')
            varnames.append('ugmt_highest_muon'+eta_min_str+eta_max_str+qMin_str+'_varBin_pt')
            varnames.append('bmtf_ugmt_highest_muon'+eta_min_str+eta_max_str+qMin_str+'_varBin_pt')
            varnames.append('omtf_ugmt_highest_muon'+eta_min_str+eta_max_str+qMin_str+'_varBin_pt')
            varnames.append('emtf_ugmt_highest_muon'+eta_min_str+eta_max_str+qMin_str+'_varBin_pt')
            varnames.append('bmtf_highest_muon'+eta_min_str+eta_max_str+qMin_str+'_varBin_pt')
            varnames.append('omtf_highest_muon'+eta_min_str+eta_max_str+qMin_str+'_varBin_pt')
            varnames.append('emtf_highest_muon'+eta_min_str+eta_max_str+qMin_str+'_varBin_pt')
            binnings['gmt_highest_muon'+eta_min_str+eta_max_str+qMin_str+'_varBin_pt'] = pt_bins+['('+eta_title+', '+qTitle+') p_{T}', 'GeV/c']
            binnings['ugmt_highest_muon'+eta_min_str+eta_max_str+qMin_str+'_varBin_pt'] = pt_bins+['('+eta_title+', '+qTitle+') p_{T}', 'GeV/c']
            binnings['bmtf_ugmt_highest_muon'+eta_min_str+eta_max_str+qMin_str+'_varBin_pt'] = pt_bins+['('+eta_title+', '+qTitle+') p_{T}', 'GeV/c']
            binnings['omtf_ugmt_highest_muon'+eta_min_str+eta_max_str+qMin_str+'_varBin_pt'] = pt_bins+['('+eta_title+', '+qTitle+') p_{T}', 'GeV/c']
            binnings['emtf_ugmt_highest_muon'+eta_min_str+eta_max_str+qMin_str+'_varBin_pt'] = pt_bins+['('+eta_title+', '+qTitle+') p_{T}', 'GeV/c']
            binnings['bmtf_highest_muon'+eta_min_str+eta_max_str+qMin_str+'_varBin_pt'] = pt_bins+['('+eta_title+', '+qTitle+') p_{T}', 'GeV/c']
            binnings['omtf_highest_muon'+eta_min_str+eta_max_str+qMin_str+'_varBin_pt'] = pt_bins+['('+eta_title+', '+qTitle+') p_{T}', 'GeV/c']
            binnings['emtf_highest_muon'+eta_min_str+eta_max_str+qMin_str+'_varBin_pt'] = pt_bins+['('+eta_title+', '+qTitle+') p_{T}', 'GeV/c']

        for qual in range(16):
            qual_str = '_q'+str(qual)
            qualTitle = 'q'+str(qual)

            varnames.append('gmt_muon'+eta_min_str+eta_max_str+qual_str+'_pt')
            varnames.append('ugmt_muon'+eta_min_str+eta_max_str+qual_str+'_pt')
            varnames.append('bmtf_ugmt_muon'+eta_min_str+eta_max_str+qual_str+'_pt')
            varnames.append('omtf_ugmt_muon'+eta_min_str+eta_max_str+qual_str+'_pt')
            varnames.append('emtf_ugmt_muon'+eta_min_str+eta_max_str+qual_str+'_pt')
            varnames.append('bmtf_muon'+eta_min_str+eta_max_str+qual_str+'_pt')
            varnames.append('omtf_muon'+eta_min_str+eta_max_str+qual_str+'_pt')
            varnames.append('emtf_muon'+eta_min_str+eta_max_str+qual_str+'_pt')
            binnings['gmt_muon'+eta_min_str+eta_max_str+qual_str+'_pt'] = [100, 0, 100, '('+eta_title+', '+qualTitle+') p_{T}', 'GeV/c']
            binnings['ugmt_muon'+eta_min_str+eta_max_str+qual_str+'_pt'] = [100, 0, 100, '('+eta_title+', '+qualTitle+') p_{T}', 'GeV/c']
            binnings['bmtf_ugmt_muon'+eta_min_str+eta_max_str+qual_str+'_pt'] = [100, 0, 100, '('+eta_title+', '+qualTitle+') p_{T}', 'GeV/c']
            binnings['omtf_ugmt_muon'+eta_min_str+eta_max_str+qual_str+'_pt'] = [100, 0, 100, '('+eta_title+', '+qualTitle+') p_{T}', 'GeV/c']
            binnings['emtf_ugmt_muon'+eta_min_str+eta_max_str+qual_str+'_pt'] = [100, 0, 100, '('+eta_title+', '+qualTitle+') p_{T}', 'GeV/c']
            binnings['bmtf_muon'+eta_min_str+eta_max_str+qual_str+'_pt'] = [100, 0, 100, '('+eta_title+', '+qualTitle+') p_{T}', 'GeV/c']
            binnings['omtf_muon'+eta_min_str+eta_max_str+qual_str+'_pt'] = [100, 0, 100, '('+eta_title+', '+qualTitle+') p_{T}', 'GeV/c']
            binnings['emtf_muon'+eta_min_str+eta_max_str+qual_str+'_pt'] = [100, 0, 100, '('+eta_title+', '+qualTitle+') p_{T}', 'GeV/c']

            varnames.append('gmt_muon'+eta_min_str+eta_max_str+qual_str+'_varBin_pt')
            varnames.append('ugmt_muon'+eta_min_str+eta_max_str+qual_str+'_varBin_pt')
            varnames.append('bmtf_ugmt_muon'+eta_min_str+eta_max_str+qual_str+'_varBin_pt')
            varnames.append('omtf_ugmt_muon'+eta_min_str+eta_max_str+qual_str+'_varBin_pt')
            varnames.append('emtf_ugmt_muon'+eta_min_str+eta_max_str+qual_str+'_varBin_pt')
            varnames.append('bmtf_muon'+eta_min_str+eta_max_str+qual_str+'_varBin_pt')
            varnames.append('omtf_muon'+eta_min_str+eta_max_str+qual_str+'_varBin_pt')
            varnames.append('emtf_muon'+eta_min_str+eta_max_str+qual_str+'_varBin_pt')
            binnings['gmt_muon'+eta_min_str+eta_max_str+qual_str+'_varBin_pt'] = pt_bins+['('+eta_title+', '+qualTitle+') p_{T}', 'GeV/c']
            binnings['ugmt_muon'+eta_min_str+eta_max_str+qual_str+'_varBin_pt'] = pt_bins+['('+eta_title+', '+qualTitle+') p_{T}', 'GeV/c']
            binnings['bmtf_ugmt_muon'+eta_min_str+eta_max_str+qual_str+'_varBin_pt'] = pt_bins+['('+eta_title+', '+qualTitle+') p_{T}', 'GeV/c']
            binnings['omtf_ugmt_muon'+eta_min_str+eta_max_str+qual_str+'_varBin_pt'] = pt_bins+['('+eta_title+', '+qualTitle+') p_{T}', 'GeV/c']
            binnings['emtf_ugmt_muon'+eta_min_str+eta_max_str+qual_str+'_varBin_pt'] = pt_bins+['('+eta_title+', '+qualTitle+') p_{T}', 'GeV/c']
            binnings['bmtf_muon'+eta_min_str+eta_max_str+qual_str+'_varBin_pt'] = pt_bins+['('+eta_title+', '+qualTitle+') p_{T}', 'GeV/c']
            binnings['omtf_muon'+eta_min_str+eta_max_str+qual_str+'_varBin_pt'] = pt_bins+['('+eta_title+', '+qualTitle+') p_{T}', 'GeV/c']
            binnings['emtf_muon'+eta_min_str+eta_max_str+qual_str+'_varBin_pt'] = pt_bins+['('+eta_title+', '+qualTitle+') p_{T}', 'GeV/c']

    ##########################################################################
    for threshold in thresholds:
        thr_str = '_ptmin'+str(threshold)
        thr_title = 'p_{T} > '+str(threshold)+' GeV/c'

        for qMin in qualities:
            qMin_str = '_qmin'+str(qMin)
            qTitle = 'q>'+str(qMin)

            varnames.append('gmt_muon'+thr_str+qMin_str+'_eta')
            varnames.append('ugmt_muon'+thr_str+qMin_str+'_eta')
            varnames.append('bmtf_ugmt_muon'+thr_str+qMin_str+'_eta')
            varnames.append('omtf_ugmt_muon'+thr_str+qMin_str+'_eta')
            varnames.append('emtf_ugmt_muon'+thr_str+qMin_str+'_eta')
            varnames.append('bmtf_muon'+thr_str+qMin_str+'_eta')
            varnames.append('omtf_muon'+thr_str+qMin_str+'_eta')
            varnames.append('emtf_muon'+thr_str+qMin_str+'_eta')
            binnings['gmt_muon'+thr_str+qMin_str+'_eta'] = eta_bins+['GMT #mu ('+thr_title+', '+qTitle+') #eta']
            binnings['ugmt_muon'+thr_str+qMin_str+'_eta'] = eta_bins+['uGMT #mu ('+thr_title+', '+qTitle+') #eta']
            binnings['bmtf_ugmt_muon'+thr_str+qMin_str+'_eta'] = eta_bins+['BMTF uGMT #mu ('+thr_title+', '+qTitle+') #eta']
            binnings['omtf_ugmt_muon'+thr_str+qMin_str+'_eta'] = eta_bins+['OMTF uGMT #mu ('+thr_title+', '+qTitle+') #eta']
            binnings['emtf_ugmt_muon'+thr_str+qMin_str+'_eta'] = eta_bins+['EMTF uGMT #mu ('+thr_title+', '+qTitle+') #eta']
            binnings['bmtf_muon'+thr_str+qMin_str+'_eta'] = eta_bins+['BMTF #mu ('+thr_title+', '+qTitle+') #eta']
            binnings['omtf_muon'+thr_str+qMin_str+'_eta'] = eta_bins+['OMTF #mu ('+thr_title+', '+qTitle+') #eta']
            binnings['emtf_muon'+thr_str+qMin_str+'_eta'] = eta_bins+['EMTF #mu ('+thr_title+', '+qTitle+') #eta']

        for qual in range(16):
            qual_str = '_q'+str(qual)
            qualTitle = 'q'+str(qual)

            varnames.append('gmt_muon'+thr_str+qual_str+'_eta')
            varnames.append('ugmt_muon'+thr_str+qual_str+'_eta')
            varnames.append('bmtf_ugmt_muon'+thr_str+qual_str+'_eta')
            varnames.append('omtf_ugmt_muon'+thr_str+qual_str+'_eta')
            varnames.append('emtf_ugmt_muon'+thr_str+qual_str+'_eta')
            varnames.append('bmtf_muon'+thr_str+qual_str+'_eta')
            varnames.append('omtf_muon'+thr_str+qual_str+'_eta')
            varnames.append('emtf_muon'+thr_str+qual_str+'_eta')
            binnings['gmt_muon'+thr_str+qual_str+'_eta'] = eta_bins+['GMT #mu ('+thr_title+', '+qualTitle+') #eta']
            binnings['ugmt_muon'+thr_str+qual_str+'_eta'] = eta_bins+['uGMT #mu ('+thr_title+', '+qualTitle+') #eta']
            binnings['bmtf_ugmt_muon'+thr_str+qual_str+'_eta'] = eta_bins+['BMTF uGMT #mu ('+thr_title+', '+qualTitle+') #eta']
            binnings['omtf_ugmt_muon'+thr_str+qual_str+'_eta'] = eta_bins+['OMTF uGMT #mu ('+thr_title+', '+qualTitle+') #eta']
            binnings['emtf_ugmt_muon'+thr_str+qual_str+'_eta'] = eta_bins+['EMTF uGMT #mu ('+thr_title+', '+qualTitle+') #eta']
            binnings['bmtf_muon'+thr_str+qual_str+'_eta'] = eta_bins+['BMTF #mu ('+thr_title+', '+qualTitle+') #eta']
            binnings['omtf_muon'+thr_str+qual_str+'_eta'] = eta_bins+['OMTF #mu ('+thr_title+', '+qualTitle+') #eta']
            binnings['emtf_muon'+thr_str+qual_str+'_eta'] = eta_bins+['EMTF #mu ('+thr_title+', '+qualTitle+') #eta']

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


def analyse(evt, hm, eta_ranges, thresholds, qualities):
    # USER HOOK
    # do what you want to do with the ntuples here
    # example:
    #print "GMT:", evt.gmt.N
    #print "UGMT:", evt.upgrade.n

    bx_min = 0
    bx_max = 0

    #gmt_muon_idcs = MuonSelections.select_gmt_muons(evt.gmt, pt_min=0.5, bx_min=bx_min, bx_max=bx_max, pos_eta=pos_eta, neg_eta=neg_eta)
    gmt_muon_idcs = [] # don't fill GMT histograms
    ugmt_muon_idcs = MuonSelections.select_ugmt_muons(evt.upgrade, pt_min=0.5, bx_min=bx_min, bx_max=bx_max, pos_eta=pos_eta, neg_eta=neg_eta)
    #bmtf_muon_idcs = MuonSelections.select_tf_muons(evt.upgradeBmtf, pt_min=0.5, tftype=0, pos_eta=pos_eta, neg_eta=neg_eta)
    #omtf_muon_idcs = MuonSelections.select_tf_muons(evt.upgradeOmtf, pt_min=0.5, tftype=1, pos_eta=pos_eta, neg_eta=neg_eta)
    #emtf_muon_idcs = MuonSelections.select_tf_muons(evt.upgradeEmtf, pt_min=0.5, tftype=2, pos_eta=pos_eta, neg_eta=neg_eta)
    bmtf_muon_idcs = MuonSelections.select_tf_muons(evt.upgradeBmtf, pt_min=0.5, pos_eta=pos_eta, neg_eta=neg_eta)
    omtf_muon_idcs = MuonSelections.select_tf_muons(evt.upgradeOmtf, pt_min=0.5, pos_eta=pos_eta, neg_eta=neg_eta)
    emtf_muon_idcs = MuonSelections.select_tf_muons(evt.upgradeEmtf, pt_min=0.5, pos_eta=pos_eta, neg_eta=neg_eta)

    for eta_range in eta_ranges:
        eta_min = eta_range[0]
        eta_max = eta_range[1]
        eta_min_str = '_absEtaMin'+str(eta_min)
        eta_max_str = '_absEtaMax'+str(eta_max)

        eta_gmt_muon_idcs = MuonSelections.select_gmt_muons(evt.gmt, abs_eta_min=eta_min, abs_eta_max=eta_max, idcs=gmt_muon_idcs)
        eta_ugmt_muon_idcs = MuonSelections.select_ugmt_muons(evt.upgrade, abs_eta_min=eta_min, abs_eta_max=eta_max, idcs=ugmt_muon_idcs)
        eta_bmtf_muon_idcs = MuonSelections.select_tf_muons(evt.upgradeBmtf, abs_eta_min=eta_min, abs_eta_max=eta_max, idcs=bmtf_muon_idcs)
        eta_omtf_muon_idcs = MuonSelections.select_tf_muons(evt.upgradeOmtf, abs_eta_min=eta_min, abs_eta_max=eta_max, idcs=omtf_muon_idcs)
        eta_emtf_muon_idcs = MuonSelections.select_tf_muons(evt.upgradeEmtf, abs_eta_min=eta_min, abs_eta_max=eta_max, idcs=emtf_muon_idcs)

        for threshold in thresholds:
            thr_str = '_ptmin'+str(threshold)

            eta_thr_gmt_muon_idcs = MuonSelections.select_gmt_muons(evt.gmt, pt_min=threshold, idcs=eta_gmt_muon_idcs)
            eta_thr_ugmt_muon_idcs = MuonSelections.select_ugmt_muons(evt.upgrade, pt_min=threshold, idcs=eta_ugmt_muon_idcs)
            eta_thr_bmtf_muon_idcs = MuonSelections.select_tf_muons(evt.upgradeBmtf, pt_min=threshold, idcs=eta_bmtf_muon_idcs)
            eta_thr_omtf_muon_idcs = MuonSelections.select_tf_muons(evt.upgradeOmtf, pt_min=threshold, idcs=eta_omtf_muon_idcs)
            eta_thr_emtf_muon_idcs = MuonSelections.select_tf_muons(evt.upgradeEmtf, pt_min=threshold, idcs=eta_emtf_muon_idcs)

            for i in eta_thr_gmt_muon_idcs:
                hm.fill('gmt_muon'+eta_min_str+eta_max_str+thr_str+'_qual', evt.gmt.Qual[i])
            for i in eta_thr_ugmt_muon_idcs:
                hm.fill('ugmt_muon'+eta_min_str+eta_max_str+thr_str+'_qual', evt.upgrade.muonQual[i])
                tftype = MuonSelections.getTfTypeFromTfMuonIdx(evt.upgrade.muonTfMuonIdx[i])
                if tftype is 0:
                    hm.fill('bmtf_ugmt_muon'+eta_min_str+eta_max_str+thr_str+'_qual', evt.upgrade.muonQual[i])
                elif tftype is 1:
                    hm.fill('omtf_ugmt_muon'+eta_min_str+eta_max_str+thr_str+'_qual', evt.upgrade.muonQual[i])
                elif tftype is 2:
                    hm.fill('emtf_ugmt_muon'+eta_min_str+eta_max_str+thr_str+'_qual', evt.upgrade.muonQual[i])
            for i in eta_thr_bmtf_muon_idcs:
                hm.fill('bmtf_muon'+eta_min_str+eta_max_str+thr_str+'_qual', evt.upgradeBmtf.tfMuonHwQual[i])
            for i in eta_thr_omtf_muon_idcs:
                hm.fill('omtf_muon'+eta_min_str+eta_max_str+thr_str+'_qual', evt.upgradeOmtf.tfMuonHwQual[i])
            for i in eta_thr_emtf_muon_idcs:
                hm.fill('emtf_muon'+eta_min_str+eta_max_str+thr_str+'_qual', evt.upgradeEmtf.tfMuonHwQual[i])

            for qMin in qualities:
                qMin_str = '_qmin'+str(qMin)

                eta_thr_q_gmt_muon_idcs = MuonSelections.select_gmt_muons(evt.gmt, qual_min=qMin, idcs=eta_thr_gmt_muon_idcs)
                eta_thr_q_ugmt_muon_idcs = MuonSelections.select_ugmt_muons(evt.upgrade, qual_min=qMin, idcs=eta_thr_ugmt_muon_idcs)
                eta_thr_q_bmtf_muon_idcs = MuonSelections.select_tf_muons(evt.upgradeBmtf, qual_min=qMin, idcs=eta_thr_bmtf_muon_idcs)
                eta_thr_q_omtf_muon_idcs = MuonSelections.select_tf_muons(evt.upgradeOmtf, qual_min=qMin, idcs=eta_thr_omtf_muon_idcs)
                eta_thr_q_emtf_muon_idcs = MuonSelections.select_tf_muons(evt.upgradeEmtf, qual_min=qMin, idcs=eta_thr_emtf_muon_idcs)

                bmtfUgmtCtr = 0
                omtfUgmtCtr = 0
                emtfUgmtCtr = 0

                for i in eta_thr_q_gmt_muon_idcs:
                    hm.fill('gmt_muon'+eta_min_str+eta_max_str+thr_str+qMin_str+'_phi', Matcher.norm_phi(evt.gmt.Phi[i]))
                for i in eta_thr_q_ugmt_muon_idcs:
                    hm.fill('ugmt_muon'+eta_min_str+eta_max_str+thr_str+qMin_str+'_phi', evt.upgrade.muonPhi[i])
                    tftype = MuonSelections.getTfTypeFromTfMuonIdx(evt.upgrade.muonTfMuonIdx[i])
                    if tftype is 0:
                        hm.fill('bmtf_ugmt_muon'+eta_min_str+eta_max_str+thr_str+qMin_str+'_phi', evt.upgrade.muonPhi[i])
                        bmtfUgmtCtr += 1
                    elif tftype is 1:
                        hm.fill('omtf_ugmt_muon'+eta_min_str+eta_max_str+thr_str+qMin_str+'_phi', evt.upgrade.muonPhi[i])
                        omtfUgmtCtr += 1
                    elif tftype is 2:
                        hm.fill('emtf_ugmt_muon'+eta_min_str+eta_max_str+thr_str+qMin_str+'_phi', evt.upgrade.muonPhi[i])
                        emtfUgmtCtr += 1
                for i in eta_thr_q_bmtf_muon_idcs:
                    hm.fill('bmtf_muon'+eta_min_str+eta_max_str+thr_str+qMin_str+'_phi', Matcher.norm_phi(evt.upgradeBmtf.tfMuonHwPhi[i]*phiScale))
                for i in eta_thr_q_omtf_muon_idcs:
                    hm.fill('omtf_muon'+eta_min_str+eta_max_str+thr_str+qMin_str+'_phi', Matcher.norm_phi(evt.upgradeOmtf.tfMuonHwPhi[i]*phiScale))
                for i in eta_thr_q_emtf_muon_idcs:
                    hm.fill('emtf_muon'+eta_min_str+eta_max_str+thr_str+qMin_str+'_phi', Matcher.norm_phi(evt.upgradeEmtf.tfMuonHwPhi[i]*phiScale))

                hm.fill('n_gmt_muons'+eta_min_str+eta_max_str+thr_str+qMin_str, len(eta_thr_q_gmt_muon_idcs))
                hm.fill('n_ugmt_muons'+eta_min_str+eta_max_str+thr_str+qMin_str, len(eta_thr_q_ugmt_muon_idcs))
                hm.fill('n_bmtf_ugmt_muons'+eta_min_str+eta_max_str+thr_str+qMin_str, bmtfUgmtCtr)
                hm.fill('n_omtf_ugmt_muons'+eta_min_str+eta_max_str+thr_str+qMin_str, omtfUgmtCtr)
                hm.fill('n_emtf_ugmt_muons'+eta_min_str+eta_max_str+thr_str+qMin_str, emtfUgmtCtr)
                hm.fill('n_bmtf_muons'+eta_min_str+eta_max_str+thr_str+qMin_str, len(eta_thr_q_bmtf_muon_idcs))
                hm.fill('n_omtf_muons'+eta_min_str+eta_max_str+thr_str+qMin_str, len(eta_thr_q_omtf_muon_idcs))
                hm.fill('n_emtf_muons'+eta_min_str+eta_max_str+thr_str+qMin_str, len(eta_thr_q_emtf_muon_idcs))

            for qual in range(16):
                qual_str = '_q'+str(qual)

                for i in eta_thr_gmt_muon_idcs:
                    if evt.gmt.Qual[i] == qual:
                        hm.fill('gmt_muon'+eta_min_str+eta_max_str+thr_str+qual_str+'_phi', Matcher.norm_phi(evt.gmt.Phi[i]))
                for i in eta_thr_ugmt_muon_idcs:
                    if evt.upgrade.muonQual[i] == qual:
                        hm.fill('ugmt_muon'+eta_min_str+eta_max_str+thr_str+qual_str+'_phi', evt.upgrade.muonPhi[i])
                        tftype = MuonSelections.getTfTypeFromTfMuonIdx(evt.upgrade.muonTfMuonIdx[i])
                        if tftype is 0:
                            hm.fill('bmtf_ugmt_muon'+eta_min_str+eta_max_str+thr_str+qual_str+'_phi', evt.upgrade.muonPhi[i])
                        elif tftype is 1:
                            hm.fill('omtf_ugmt_muon'+eta_min_str+eta_max_str+thr_str+qual_str+'_phi', evt.upgrade.muonPhi[i])
                        elif tftype is 2:
                            hm.fill('emtf_ugmt_muon'+eta_min_str+eta_max_str+thr_str+qual_str+'_phi', evt.upgrade.muonPhi[i])
                for i in eta_thr_bmtf_muon_idcs:
                    if evt.upgradeBmtf.tfMuonHwQual[i] == qual:
                        hm.fill('bmtf_muon'+eta_min_str+eta_max_str+thr_str+qual_str+'_phi', Matcher.norm_phi(evt.upgradeBmtf.tfMuonHwPhi[i]*phiScale))
                for i in eta_thr_omtf_muon_idcs:
                    if evt.upgradeOmtf.tfMuonHwQual[i] == qual:
                        hm.fill('omtf_muon'+eta_min_str+eta_max_str+thr_str+qual_str+'_phi', Matcher.norm_phi(evt.upgradeOmtf.tfMuonHwPhi[i]*phiScale))
                for i in eta_thr_emtf_muon_idcs:
                    if evt.upgradeEmtf.tfMuonHwQual[i] == qual:
                        hm.fill('emtf_muon'+eta_min_str+eta_max_str+thr_str+qual_str+'_phi', Matcher.norm_phi(evt.upgradeEmtf.tfMuonHwPhi[i]*phiScale))

        for qMin in qualities:
            qMin_str = '_qmin'+str(qMin)

            eta_q_gmt_muon_idcs = MuonSelections.select_gmt_muons(evt.gmt, qual_min=qMin, idcs=eta_gmt_muon_idcs)
            eta_q_ugmt_muon_idcs = MuonSelections.select_ugmt_muons(evt.upgrade, qual_min=qMin, idcs=eta_ugmt_muon_idcs)
            eta_q_bmtf_muon_idcs = MuonSelections.select_tf_muons(evt.upgradeBmtf, qual_min=qMin, idcs=eta_bmtf_muon_idcs)
            eta_q_omtf_muon_idcs = MuonSelections.select_tf_muons(evt.upgradeOmtf, qual_min=qMin, idcs=eta_omtf_muon_idcs)
            eta_q_emtf_muon_idcs = MuonSelections.select_tf_muons(evt.upgradeEmtf, qual_min=qMin, idcs=eta_emtf_muon_idcs)

            for i in eta_q_gmt_muon_idcs:
                hm.fill('gmt_muon'+eta_min_str+eta_max_str+qMin_str+'_pt', evt.gmt.Pt[i])
                hm.fill('gmt_muon'+eta_min_str+eta_max_str+qMin_str+'_varBin_pt', evt.gmt.Pt[i])
            for i in eta_q_ugmt_muon_idcs:
                hm.fill('ugmt_muon'+eta_min_str+eta_max_str+qMin_str+'_pt', evt.upgrade.muonEt[i])
                hm.fill('ugmt_muon'+eta_min_str+eta_max_str+qMin_str+'_varBin_pt', evt.upgrade.muonEt[i])
                tftype = MuonSelections.getTfTypeFromTfMuonIdx(evt.upgrade.muonTfMuonIdx[i])
                if tftype is 0:
                    hm.fill('bmtf_ugmt_muon'+eta_min_str+eta_max_str+qMin_str+'_pt', evt.upgrade.muonEt[i])
                    hm.fill('bmtf_ugmt_muon'+eta_min_str+eta_max_str+qMin_str+'_varBin_pt', evt.upgrade.muonEt[i])
                elif tftype is 1:
                    hm.fill('omtf_ugmt_muon'+eta_min_str+eta_max_str+qMin_str+'_pt', evt.upgrade.muonEt[i])
                    hm.fill('omtf_ugmt_muon'+eta_min_str+eta_max_str+qMin_str+'_varBin_pt', evt.upgrade.muonEt[i])
                elif tftype is 2:
                    hm.fill('emtf_ugmt_muon'+eta_min_str+eta_max_str+qMin_str+'_pt', evt.upgrade.muonEt[i])
                    hm.fill('emtf_ugmt_muon'+eta_min_str+eta_max_str+qMin_str+'_varBin_pt', evt.upgrade.muonEt[i])
            for i in eta_q_bmtf_muon_idcs:
                hm.fill('bmtf_muon'+eta_min_str+eta_max_str+qMin_str+'_pt', evt.upgradeBmtf.tfMuonHwPt[i]*ptScale)
                hm.fill('bmtf_muon'+eta_min_str+eta_max_str+qMin_str+'_varBin_pt', evt.upgradeBmtf.tfMuonHwPt[i]*ptScale)
            for i in eta_q_omtf_muon_idcs:
                hm.fill('omtf_muon'+eta_min_str+eta_max_str+qMin_str+'_pt', evt.upgradeOmtf.tfMuonHwPt[i]*ptScale)
                hm.fill('omtf_muon'+eta_min_str+eta_max_str+qMin_str+'_varBin_pt', evt.upgradeOmtf.tfMuonHwPt[i]*ptScale)
            for i in eta_q_emtf_muon_idcs:
                hm.fill('emtf_muon'+eta_min_str+eta_max_str+qMin_str+'_pt', evt.upgradeEmtf.tfMuonHwPt[i]*ptScale)
                hm.fill('emtf_muon'+eta_min_str+eta_max_str+qMin_str+'_varBin_pt', evt.upgradeEmtf.tfMuonHwPt[i]*ptScale)

            if len(eta_q_gmt_muon_idcs):
                highestPt = get_highest_pt(evt.gmt, eta_q_gmt_muon_idcs, gmt=True)
                hm.fill('gmt_highest_muon'+eta_min_str+eta_max_str+qMin_str+'_pt', highestPt)
                hm.fill('gmt_highest_muon'+eta_min_str+eta_max_str+qMin_str+'_varBin_pt', highestPt)
            if len(eta_q_ugmt_muon_idcs):
                highestPtIdx = get_highest_pt_idx(evt.upgrade, eta_q_ugmt_muon_idcs)
                highestPt = evt.upgrade.muonEt[highestPtIdx]
                hm.fill('ugmt_highest_muon'+eta_min_str+eta_max_str+qMin_str+'_pt', highestPt)
                hm.fill('ugmt_highest_muon'+eta_min_str+eta_max_str+qMin_str+'_varBin_pt', highestPt)
                tftype = MuonSelections.getTfTypeFromTfMuonIdx(evt.upgrade.muonTfMuonIdx[highestPtIdx])
                if tftype is 0:
                    hm.fill('bmtf_ugmt_highest_muon'+eta_min_str+eta_max_str+qMin_str+'_pt', highestPt)
                    hm.fill('bmtf_ugmt_highest_muon'+eta_min_str+eta_max_str+qMin_str+'_varBin_pt', highestPt)
                elif tftype is 1:
                    hm.fill('omtf_ugmt_highest_muon'+eta_min_str+eta_max_str+qMin_str+'_pt', highestPt)
                    hm.fill('omtf_ugmt_highest_muon'+eta_min_str+eta_max_str+qMin_str+'_varBin_pt', highestPt)
                elif tftype is 2:
                    hm.fill('emtf_ugmt_highest_muon'+eta_min_str+eta_max_str+qMin_str+'_pt', highestPt)
                    hm.fill('emtf_ugmt_highest_muon'+eta_min_str+eta_max_str+qMin_str+'_varBin_pt', highestPt)
            if len(eta_q_bmtf_muon_idcs):
                highestPt = get_highest_pt(evt.upgradeBmtf, eta_q_bmtf_muon_idcs, tf=True)
                hm.fill('bmtf_highest_muon'+eta_min_str+eta_max_str+qMin_str+'_pt', highestPt)
                hm.fill('bmtf_highest_muon'+eta_min_str+eta_max_str+qMin_str+'_varBin_pt', highestPt)
            if len(eta_q_omtf_muon_idcs):
                highestPt = get_highest_pt(evt.upgradeOmtf, eta_q_omtf_muon_idcs, tf=True)
                hm.fill('omtf_highest_muon'+eta_min_str+eta_max_str+qMin_str+'_pt', highestPt)
                hm.fill('omtf_highest_muon'+eta_min_str+eta_max_str+qMin_str+'_varBin_pt', highestPt)
            if len(eta_q_emtf_muon_idcs):
                highestPt = get_highest_pt(evt.upgradeEmtf, eta_q_emtf_muon_idcs, tf=True)
                hm.fill('emtf_highest_muon'+eta_min_str+eta_max_str+qMin_str+'_pt', highestPt)
                hm.fill('emtf_highest_muon'+eta_min_str+eta_max_str+qMin_str+'_varBin_pt', highestPt)

        for qual in range(16):
            qual_str = '_q'+str(qual)

            for i in eta_gmt_muon_idcs:
                if evt.gmt.Qual[i] == qual:
                    hm.fill('gmt_muon'+eta_min_str+eta_max_str+qual_str+'_pt', evt.gmt.Pt[i])
                    hm.fill('gmt_muon'+eta_min_str+eta_max_str+qual_str+'_varBin_pt', evt.gmt.Pt[i])
            for i in eta_ugmt_muon_idcs:
                if evt.upgrade.muonQual[i] == qual:
                    hm.fill('ugmt_muon'+eta_min_str+eta_max_str+qual_str+'_pt', evt.upgrade.muonEt[i])
                    hm.fill('ugmt_muon'+eta_min_str+eta_max_str+qual_str+'_varBin_pt', evt.upgrade.muonEt[i])
                    tftype = MuonSelections.getTfTypeFromTfMuonIdx(evt.upgrade.muonTfMuonIdx[i])
                    if tftype is 0:
                        hm.fill('bmtf_ugmt_muon'+eta_min_str+eta_max_str+qual_str+'_pt', evt.upgrade.muonEt[i])
                        hm.fill('bmtf_ugmt_muon'+eta_min_str+eta_max_str+qual_str+'_varBin_pt', evt.upgrade.muonEt[i])
                    elif tftype is 1:
                        hm.fill('omtf_ugmt_muon'+eta_min_str+eta_max_str+qual_str+'_pt', evt.upgrade.muonEt[i])
                        hm.fill('omtf_ugmt_muon'+eta_min_str+eta_max_str+qual_str+'_varBin_pt', evt.upgrade.muonEt[i])
                    elif tftype is 2:
                        hm.fill('emtf_ugmt_muon'+eta_min_str+eta_max_str+qual_str+'_pt', evt.upgrade.muonEt[i])
                        hm.fill('emtf_ugmt_muon'+eta_min_str+eta_max_str+qual_str+'_varBin_pt', evt.upgrade.muonEt[i])
            for i in eta_bmtf_muon_idcs:
                if evt.upgradeBmtf.tfMuonHwQual[i] == qual:
                    hm.fill('bmtf_muon'+eta_min_str+eta_max_str+qual_str+'_pt', evt.upgradeBmtf.tfMuonHwPt[i]*ptScale)
                    hm.fill('bmtf_muon'+eta_min_str+eta_max_str+qual_str+'_varBin_pt', evt.upgradeBmtf.tfMuonHwPt[i]*ptScale)
            for i in eta_omtf_muon_idcs:
                if evt.upgradeOmtf.tfMuonHwQual[i] == qual:
                    hm.fill('omtf_muon'+eta_min_str+eta_max_str+qual_str+'_pt', evt.upgradeOmtf.tfMuonHwPt[i]*ptScale)
                    hm.fill('omtf_muon'+eta_min_str+eta_max_str+qual_str+'_varBin_pt', evt.upgradeOmtf.tfMuonHwPt[i]*ptScale)
            for i in eta_emtf_muon_idcs:
                if evt.upgradeEmtf.tfMuonHwQual[i] == qual:
                    hm.fill('emtf_muon'+eta_min_str+eta_max_str+qual_str+'_pt', evt.upgradeEmtf.tfMuonHwPt[i]*ptScale)
                    hm.fill('emtf_muon'+eta_min_str+eta_max_str+qual_str+'_varBin_pt', evt.upgradeEmtf.tfMuonHwPt[i]*ptScale)

    for threshold in thresholds:
        thr_str = '_ptmin'+str(threshold)

        thr_gmt_muon_idcs = MuonSelections.select_gmt_muons(evt.gmt, pt_min=threshold, idcs=gmt_muon_idcs)
        thr_ugmt_muon_idcs = MuonSelections.select_ugmt_muons(evt.upgrade, pt_min=threshold, idcs=ugmt_muon_idcs)
        thr_bmtf_muon_idcs = MuonSelections.select_tf_muons(evt.upgradeBmtf, pt_min=threshold, idcs=bmtf_muon_idcs)
        thr_omtf_muon_idcs = MuonSelections.select_tf_muons(evt.upgradeOmtf, pt_min=threshold, idcs=omtf_muon_idcs)
        thr_emtf_muon_idcs = MuonSelections.select_tf_muons(evt.upgradeEmtf, pt_min=threshold, idcs=emtf_muon_idcs)

        for qMin in qualities:
            qMin_str = '_qmin'+str(qMin)

            thr_q_gmt_muon_idcs = MuonSelections.select_gmt_muons(evt.gmt, pt_min=threshold, qual_min=qMin, idcs=gmt_muon_idcs)
            thr_q_ugmt_muon_idcs = MuonSelections.select_ugmt_muons(evt.upgrade, pt_min=threshold, qual_min=qMin, idcs=ugmt_muon_idcs)
            thr_q_bmtf_muon_idcs = MuonSelections.select_tf_muons(evt.upgradeBmtf, pt_min=threshold, qual_min=qMin, idcs=bmtf_muon_idcs)
            thr_q_omtf_muon_idcs = MuonSelections.select_tf_muons(evt.upgradeOmtf, pt_min=threshold, qual_min=qMin, idcs=omtf_muon_idcs)
            thr_q_emtf_muon_idcs = MuonSelections.select_tf_muons(evt.upgradeEmtf, pt_min=threshold, qual_min=qMin, idcs=emtf_muon_idcs)
        
            for i in thr_q_gmt_muon_idcs:
                hm.fill('gmt_muon'+thr_str+qMin_str+'_eta', evt.gmt.Eta[i])
            for i in thr_q_ugmt_muon_idcs:
                hm.fill('ugmt_muon'+thr_str+qMin_str+'_eta', evt.upgrade.muonEta[i])
                tftype = MuonSelections.getTfTypeFromTfMuonIdx(evt.upgrade.muonTfMuonIdx[i])
                if tftype is 0:
                    hm.fill('bmtf_ugmt_muon'+thr_str+qMin_str+'_eta', evt.upgrade.muonEta[i])
                elif tftype is 1:
                    hm.fill('omtf_ugmt_muon'+thr_str+qMin_str+'_eta', evt.upgrade.muonEta[i])
                elif tftype is 2:
                    hm.fill('emtf_ugmt_muon'+thr_str+qMin_str+'_eta', evt.upgrade.muonEta[i])
            for i in thr_q_bmtf_muon_idcs:
                hm.fill('bmtf_muon'+thr_str+qMin_str+'_eta', evt.upgradeBmtf.tfMuonHwEta[i]*etaScale)
            for i in thr_q_omtf_muon_idcs:
                hm.fill('omtf_muon'+thr_str+qMin_str+'_eta', evt.upgradeOmtf.tfMuonHwEta[i]*etaScale)
            for i in thr_q_emtf_muon_idcs:
                hm.fill('emtf_muon'+thr_str+qMin_str+'_eta', evt.upgradeEmtf.tfMuonHwEta[i]*etaScale)

        for qual in range(16):
            qual_str = '_q'+str(qual)

            for i in thr_gmt_muon_idcs:
                if evt.gmt.Qual[i] == qual:
                    hm.fill('gmt_muon'+thr_str+qual_str+'_eta', evt.gmt.Eta[i])
            for i in thr_ugmt_muon_idcs:
                if evt.upgrade.muonQual[i] == qual:
                    hm.fill('ugmt_muon'+thr_str+qual_str+'_eta', evt.upgrade.muonEta[i])
                    tftype = MuonSelections.getTfTypeFromTfMuonIdx(evt.upgrade.muonTfMuonIdx[i])
                    if tftype is 0:
                        hm.fill('bmtf_ugmt_muon'+thr_str+qual_str+'_eta', evt.upgrade.muonEta[i])
                    elif tftype is 1:
                        hm.fill('omtf_ugmt_muon'+thr_str+qual_str+'_eta', evt.upgrade.muonEta[i])
                    elif tftype is 2:
                        hm.fill('emtf_ugmt_muon'+thr_str+qual_str+'_eta', evt.upgrade.muonEta[i])
            for i in thr_bmtf_muon_idcs:
                if evt.upgradeBmtf.tfMuonHwQual[i] == qual:
                    hm.fill('bmtf_muon'+thr_str+qual_str+'_eta', evt.upgradeBmtf.tfMuonHwEta[i]*etaScale)
            for i in thr_omtf_muon_idcs:
                if evt.upgradeOmtf.tfMuonHwQual[i] == qual:
                    hm.fill('omtf_muon'+thr_str+qual_str+'_eta', evt.upgradeOmtf.tfMuonHwEta[i]*etaScale)
            for i in thr_emtf_muon_idcs:
                if evt.upgradeEmtf.tfMuonHwQual[i] == qual:
                    hm.fill('emtf_muon'+thr_str+qual_str+'_eta', evt.upgradeEmtf.tfMuonHwEta[i]*etaScale)

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
    print ""

    #eta_ranges = [[0, 2.5], [0, 2.1], [0, 0.83], [0.83, 1.24], [1.24, 2.5], [1.24, 2.1]]
    #thresholds = [1, 5, 10, 12, 16, 20, 24, 30]
    #qualities = range(16)
    eta_ranges = [[0, 2.5], [0, 2.1], [0, 0.83], [0.83, 1.24], [1.24, 2.5]]
    thresholds = [0, 18]
    qualities = [0, 4, 8, 12]
    # book the histograms
    hm = book_histograms(eta_ranges, thresholds, qualities)

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
        analyse(event, hm, eta_ranges, thresholds, qualities)

    # save histos to root file
    if saveHistos:
        output = root.TFile(opts.outname, 'recreate')
        output.cd()
        save_histos(hm, output)
        output.Close()

if __name__ == "__main__":
    pos_eta = True
    neg_eta = True
    saveHistos = True
    main()

