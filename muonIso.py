#!/usr/bin/env python
from ToolBox import parse_options_and_init_log
# have to do this first or ROOT masks the -h messages
opts, parser = parse_options_and_init_log()

from L1Analysis import L1Ana, L1Ntuple
from analysis_tools.plotting import HistManager, HistManager2d
from analysis_tools.selections import MuonSelections, Matcher
from analysis_tools.isolation.caloTowerIso import CaloTowerIsolator
from math import floor, ceil
import exceptions
import json
import ROOT as root

def parse_options_upgradeMuonHistos(parser):
    """
    Adds often used options to the OptionParser...
    """
    parsers = parser.add_subparsers()
    sub_parser = parsers.add_parser("muonIso")
    sub_parser.add_argument("-o", "--outname", dest="outname", default="./muon_iso_histos.root", type=str, help="A root file name where to save the histograms.")
    sub_parser.add_argument("-j", "--json", dest="json", type=str, default=None, help="A json file with good lumi sections per run.")
    sub_parser.add_argument("-r", "--runs", dest="runs", type=str, default=None, help="A string of runs to check.")
    sub_parser.add_argument("--use-extra-coord", dest="extraCoord", default=False, action="store_true", help="Use L1 extrapolated eta and phi coordinates.")
    sub_parser.add_argument("--use-reco-muon", dest="use_reco", default=False, action="store_true", help="Use RECO muons.")
    sub_parser.add_argument("--emul", dest="emul", default=False, action="store_true", help="Make emulator plots.")
    sub_parser.add_argument("--pos-charge", dest="pos_charge", default=False, action="store_true", help="Positive L1 charge only.")
    sub_parser.add_argument("--neg-charge", dest="neg_charge", default=False, action="store_true", help="Negative L1 charge only.")
    sub_parser.add_argument("--tftype", dest="tftype", type=str, default='', help="Fill L1 muons from one TF.")
    sub_parser.add_argument("--match-to-gen", dest="matchtogen", default=False, action="store_true", help="Match to GEN muons.")

    opts, unknown = parser.parse_known_args()
    return opts

def get_tftype(tf_muon_index):
    if tf_muon_index > 35 and tf_muon_index < 72:
        return 0 # BMTF
    elif tf_muon_index > 17 and tf_muon_index < 90:
        return 1 # OMTF
    else:
        return 2 # EMTF

def book_histograms(l1PtMins, emul=False):
    namePrefix = ''
    if emul:
        namePrefix += 'emu_'

    # define pt binning
    pt_bins = range(0, 60, 1)
    pt_bins += range(60, 80, 2)
    pt_bins += range(80, 100, 5)
    pt_bins += range(100, 200, 10)
    pt_bins += range(200, 300, 50)
    pt_bins.append(300)

    vars_bins = [['iet', 256, 0, 256], ['ieta', 83, -41, 42], ['iphi', 73, 0, 73], ['iqual', 16, 0, 16], ['n', 2017, 0, 4034]]
    x_title_vars = {'iet':'iE_{T}', 'ieta':'i#eta', 'iphi':'i#phi', 'iqual':'qual', 'n':'# towers'}
    x_title_units = {'iet':None, 'ieta':None, 'iphi':None, 'iqual':None, 'n':None}

    x_vars_bins_2d = [['ieta', 83, -41, 42], ['iet_ieta', 83, -41, 42]]
    y_vars_bins_2d = [['iphi', 73, 0, 73], ['iet_iphi', 73, 0, 73]]

    x_title_vars_2d = {'ieta':'i#eta', 'iet_ieta':'i#eta'}
    y_title_vars_2d = {'iphi':'i#phi', 'iet_iphi':'i#phi'}

    x_title_units_2d = {'ieta':None, 'iet_ieta':None}
    y_title_units_2d = {'iphi':None, 'iet_iphi':None}

    vars_bins_mu = [['n_mu', 15, 0, 15],
                    ['mu_pt', -1]+pt_bins,
                    ['area_1x1_iet', 256, 0, 256],
                    ['area_1x3_iet', 256, 0, 256],
                    ['area_1x5_iet', 256, 0, 256],
                    ['area_1xm2to0_iet', 256, 0, 256],
                    ['area_1x0top2_iet', 256, 0, 256],
                    ['area_3x3_iet', 256, 0, 256],
                    ['area_3xm3to0_iet', 256, 0, 256],
                    ['area_3x0top3_iet', 256, 0, 256],
                    ['area_3xm7to0_iet', 256, 0, 256],
                    ['area_3x0top7_iet', 256, 0, 256],
                    ['area_11x11_iet', 256, 0, 256],
                    ['area_11x11-1x1_iet', 256, 0, 256],
                    ['area_11x11-1x3_iet', 256, 0, 256],
                    ['area_11x11-1x5_iet', 256, 0, 256],
                    ['area_11x11-3x3_iet', 256, 0, 256],
                    ['area_1xm2to0_iet_minus_area_1x0top2_iet', 511, -255, 256],
                    ['area_1xm2to0_iet_minus_area_1x0top2_iet_mu_chg_pos', 511, -255, 256],
                    ['area_1xm2to0_iet_minus_area_1x0top2_iet_mu_chg_neg', 511, -255, 256],
                    ['area_3xm3to0_iet_over_area_3x0top3_iet', 100, 0, 10.],
                    ['area_3xm3to0_iet_over_area_3x0top3_iet_mu_chg_pos', 100, 0, 10.],
                    ['area_3xm3to0_iet_over_area_3x0top3_iet_mu_chg_neg', 100, 0, 10.],
                    ['area_3xm3to0_iet_minus_area_3x0top3_iet', 511, -255, 256],
                    ['area_3xm3to0_iet_minus_area_3x0top3_iet_mu_chg_pos', 511, -255, 256],
                    ['area_3xm3to0_iet_minus_area_3x0top3_iet_mu_chg_neg', 511, -255, 256],
                    ['area_3xm3to0_iet_minus_area_3x0top3_iet_over_area_3xm3to0_iet_plus_area_3x0top3_iet', 101, -1., 1.02],
                    ['area_3xm3to0_iet_minus_area_3x0top3_iet_over_area_3xm3to0_iet_plus_area_3x0top3_iet_times_mu_chg', 101, -1., 1.02],
                    ['area_3xm3to0_iet_minus_area_3x0top3_iet_over_area_3xm3to0_iet_plus_area_3x0top3_iet_mu_chg_pos', 101, -1., 1.02],
                    ['area_3xm3to0_iet_minus_area_3x0top3_iet_over_area_3xm3to0_iet_plus_area_3x0top3_iet_mu_chg_neg', 101, -1., 1.02],
                    ['area_3xm7to0_iet_minus_area_3x0top7_iet', 511, -255, 256],
                    ['area_3xm7to0_iet_minus_area_3x0top7_iet_mu_chg_pos', 511, -255, 256],
                    ['area_3xm7to0_iet_minus_area_3x0top7_iet_mu_chg_neg', 511, -255, 256],
                    ['area_3xm7to0_iet_minus_area_3x0top7_iet_over_area_3xm7to0_iet_plus_area_3x0top7_iet', 101, -1., 1.02],
                    ['area_3xm7to0_iet_minus_area_3x0top7_iet_over_area_3xm7to0_iet_plus_area_3x0top7_iet_times_mu_chg', 101, -1., 1.02],
                    ['area_3xm7to0_iet_minus_area_3x0top7_iet_over_area_3xm7to0_iet_plus_area_3x0top7_iet_mu_chg_pos', 101, -1., 1.02],
                    ['area_3xm7to0_iet_minus_area_3x0top7_iet_over_area_3xm7to0_iet_plus_area_3x0top7_iet_mu_chg_neg', 101, -1., 1.02],
                    ['area_11x11-1x1_over_area_11x11_iet', 51, 0, 1.02],
                    ['area_11x11-1x3_over_area_11x11_iet', 51, 0, 1.02],
                    ['area_11x11-1x5_over_area_11x11_iet', 51, 0, 1.02],
                    ['area_11x11-3x3_over_area_11x11_iet', 51, 0, 1.02],
                    ['twobytwo_area_1x1_iet', 256, 0, 256],
                    ['twobytwo_area_1xm1to0_iet', 256, 0, 256],
                    ['twobytwo_area_1x0top1_iet', 256, 0, 256],
                    ['twobytwo_area_5x5_iet', 256, 0, 256],
                    ['twobytwo_area_5x5-1x1_iet', 256, 0, 256],
                    ['twobytwo_area_1xm1to0_iet_minus_twobytwo_area_1x0top1_iet', 511, -255, 256],
                    ['twobytwo_area_1xm1to0_iet_minus_twobytwo_area_1x0top1_iet_mu_chg_pos', 511, -255, 256],
                    ['twobytwo_area_1xm1to0_iet_minus_twobytwo_area_1x0top1_iet_mu_chg_neg', 511, -255, 256],
                    ['twobytwo_area_5x5-1x1_over_twobytwo_area_5x5_iet', 51, 0, 1.02],
                    ['twobytwo_area_5x5_iet_over_mu_ipt', 100, 0, 1.],
                    ['twobytwo_area_5x5-1x1_iet_over_mu_ipt', 100, 0, 3.],
                    ['twobytwo_area_5x5_5bit_iet_over_mu_ipt', 100, 0, 3.],
                    ['twobytwo_area_5x5-1x1_5bit_iet_over_mu_ipt', 100, 0, 3.]]
    x_title_vars_mu = {'n_mu':'# L1 muons',
                       'mu_pt':'p_{T}',
                       'area_1x1_iet':'iE_{T}^{1x1}',
                       'area_1x3_iet':'iE_{T}^{1x3}',
                       'area_1x5_iet':'iE_{T}^{1x5}',
                       'area_1xm2to0_iet':'iE_{T}^{1x-2to0}',
                       'area_1x0top2_iet':'iE_{T}^{1x0to+2}',
                       'area_3x3_iet':'iE_{T}^{3x3}',
                       'area_3xm3to0_iet':'iE_{T}^{3x-3to0}',
                       'area_3x0top3_iet':'iE_{T}^{3x0to+3}',
                       'area_3xm7to0_iet':'iE_{T}^{3x-7to0}',
                       'area_3x0top7_iet':'iE_{T}^{3x0to+7}',
                       'area_11x11_iet':'iE_{T}^{11x11}',
                       'area_11x11-1x1_iet':'iE_{T}^{11x11-1x1}',
                       'area_11x11-1x3_iet':'iE_{T}^{11x11-1x3}',
                       'area_11x11-1x5_iet':'iE_{T}^{11x11-1x5}',
                       'area_11x11-3x3_iet':'iE_{T}^{11x11-3x3}',
                       'area_1xm2to0_iet_minus_area_1x0top2_iet':'iE_{T}^{1x-2to0} - iE_{T}^{1x0to+2}',
                       'area_1xm2to0_iet_minus_area_1x0top2_iet_mu_chg_pos':'positive #mu charge iE_{T}^{1x-2to0} - iE_{T}^{1x0to+2}',
                       'area_1xm2to0_iet_minus_area_1x0top2_iet_mu_chg_neg':'negative #mu charge iE_{T}^{1x-2to0} - iE_{T}^{1x0to+2}',
                       'area_3xm3to0_iet_over_area_3x0top3_iet':'iE_{T}^{3x-3to0} / iE_{T}^{3x0to+3}',
                       'area_3xm3to0_iet_over_area_3x0top3_iet_mu_chg_pos':'positive #mu charge iE_{T}^{3x-3to0} / iE_{T}^{3x0to+3}',
                       'area_3xm3to0_iet_over_area_3x0top3_iet_mu_chg_neg':'negative #mu charge iE_{T}^{3x-3to0} / iE_{T}^{3x0to+3}',
                       'area_3xm3to0_iet_minus_area_3x0top3_iet':'iE_{T}^{3x-3to0} - iE_{T}^{3x0to+3}',
                       'area_3xm3to0_iet_minus_area_3x0top3_iet_mu_chg_pos':'positive #mu charge iE_{T}^{3x-3to0} - iE_{T}^{3x0to+3}',
                       'area_3xm3to0_iet_minus_area_3x0top3_iet_mu_chg_neg':'negative #mu charge iE_{T}^{3x-3to0} - iE_{T}^{3x0to+3}',
                       'area_3xm3to0_iet_minus_area_3x0top3_iet_over_area_3xm3to0_iet_plus_area_3x0top3_iet':'(iE_{T}^{3x-3to0} - iE_{T}^{3x0to+3}) / (iE_{T}^{3x-3to0} + iE_{T}^{3x0to+3})',
                       'area_3xm3to0_iet_minus_area_3x0top3_iet_over_area_3xm3to0_iet_plus_area_3x0top3_iet_times_mu_chg':'#mu charge * (iE_{T}^{3x-3to0} - iE_{T}^{3x0to+3}) / (iE_{T}^{3x-3to0} + iE_{T}^{3x0to+3})',
                       'area_3xm3to0_iet_minus_area_3x0top3_iet_over_area_3xm3to0_iet_plus_area_3x0top3_iet_mu_chg_pos':'positive #mu charge (iE_{T}^{3x-3to0} - iE_{T}^{3x0to+3}) / (charge iE_{T}^{3x-3to0} - iE_{T}^{3x0to+3})',
                       'area_3xm3to0_iet_minus_area_3x0top3_iet_over_area_3xm3to0_iet_plus_area_3x0top3_iet_mu_chg_neg':'negative #mu charge (iE_{T}^{3x-3to0} - iE_{T}^{3x0to+3}) / (charge iE_{T}^{3x-3to0} - iE_{T}^{3x0to+3})',
                       'area_3xm7to0_iet_minus_area_3x0top7_iet':'iE_{T}^{3x-3to0} - iE_{T}^{3x0to+3}',
                       'area_3xm7to0_iet_minus_area_3x0top7_iet_mu_chg_pos':'positive #mu charge iE_{T}^{3x-7to0} - iE_{T}^{3x0to+7}',
                       'area_3xm7to0_iet_minus_area_3x0top7_iet_mu_chg_neg':'negative #mu charge iE_{T}^{3x-7to0} - iE_{T}^{3x0to+7}',
                       'area_3xm7to0_iet_minus_area_3x0top7_iet_over_area_3xm7to0_iet_plus_area_3x0top7_iet':'(iE_{T}^{3x-7to0} - iE_{T}^{3x0to+7}) / (iE_{T}^{3x-7to0} + iE_{T}^{3x0to+7})',
                       'area_3xm7to0_iet_minus_area_3x0top7_iet_over_area_3xm7to0_iet_plus_area_3x0top7_iet_times_mu_chg':'#mu charge * (iE_{T}^{3x-7to0} - iE_{T}^{3x0to+7}) / (iE_{T}^{3x-7to0} + iE_{T}^{3x0to+7})',
                       'area_3xm7to0_iet_minus_area_3x0top7_iet_over_area_3xm7to0_iet_plus_area_3x0top7_iet_mu_chg_pos':'positive #mu charge (iE_{T}^{3x-7to0} - iE_{T}^{3x0to+7}) / (charge iE_{T}^{3x-7to0} - iE_{T}^{3x0to+7})',
                       'area_3xm7to0_iet_minus_area_3x0top7_iet_over_area_3xm7to0_iet_plus_area_3x0top7_iet_mu_chg_neg':'negative #mu charge (iE_{T}^{3x-7to0} - iE_{T}^{3x0to+7}) / (charge iE_{T}^{3x-7to0} - iE_{T}^{3x0to+7})',
                       'area_11x11-1x1_over_area_11x11_iet':'iE_{T}^{11x11-1x1}/iE_{T}^{11x11}',
                       'area_11x11-1x3_over_area_11x11_iet':'iE_{T}^{11x11-1x3}/iE_{T}^{11x11}',
                       'area_11x11-1x5_over_area_11x11_iet':'iE_{T}^{11x11-1x5}/iE_{T}^{11x11}',
                       'area_11x11-3x3_over_area_11x11_iet':'iE_{T}^{11x11-3x3}/iE_{T}^{11x11}',
                       'twobytwo_area_1x1_iet':'2x2 area iE_{T}^{1x1}',
                       'twobytwo_area_1xm1to0_iet':'2x2 area iE_{T}^{1x-2to0}',
                       'twobytwo_area_1x0top1_iet':'2x2 area iE_{T}^{1x0to+2}',
                       'twobytwo_area_5x5_iet':'2x2 area iE_{T}^{5x5}',
                       'twobytwo_area_5x5-1x1_iet':'2x2 area iE_{T}^{5x5-1x1}',
                       'twobytwo_area_1xm1to0_iet_minus_twobytwo_area_1x0top1_iet':'2x2 area iE_{T}^{1x-2to0} - iE_{T}^{1x0to+2}',
                       'twobytwo_area_1xm1to0_iet_minus_twobytwo_area_1x0top1_iet_mu_chg_pos':'positive #mu charge 2x2 area iE_{T}^{1x-2to0} - iE_{T}^{1x0to+2}',
                       'twobytwo_area_1xm1to0_iet_minus_twobytwo_area_1x0top1_iet_mu_chg_neg':'negative #mu charge 2x2 area iE_{T}^{1x-2to0} - iE_{T}^{1x0to+2}',
                       'twobytwo_area_5x5-1x1_over_twobytwo_area_5x5_iet':'2x2 area iE_{T}^{5x5-1x1}/iE_{T}^{5x5}',
                       'twobytwo_area_5x5_iet_over_mu_ipt':'2x2 area iE_{T}^{5x5} / p_{T}^{#mu}',
                       'twobytwo_area_5x5-1x1_iet_over_mu_ipt':'2x2 area iE_{T}^{5x5-1x1} / p_{T}^{#mu}',
                       'twobytwo_area_5x5_5bit_iet_over_mu_ipt':'5 bit (iE_{T}^{max}=31) 2x2 area iE_{T}^{5x5} / p_{T}^{#mu}',
                       'twobytwo_area_5x5-1x1_5bit_iet_over_mu_ipt':'5 bit (iE_{T}^{max}=31) 2x2 area iE_{T}^{5x5-1x1} / p_{T}^{#mu}'}
    x_title_units_mu = {'n_mu':None,
                        'mu_pt':'GeV',
                        'area_1x1_iet':None,
                        'area_1x3_iet':None,
                        'area_1x5_iet':None,
                        'area_1xm2to0_iet':None,
                        'area_1x0top2_iet':None,
                        'area_3xm3to0_iet':None,
                        'area_3x0top3_iet':None,
                        'area_3xm7to0_iet':None,
                        'area_3x0top7_iet':None,
                        'area_3x3_iet':None,
                        'area_11x11_iet':None,
                        'area_11x11-1x1_iet':None,
                        'area_11x11-1x3_iet':None,
                        'area_11x11-1x5_iet':None,
                        'area_11x11-3x3_iet':None,
                        'area_1xm2to0_iet_minus_area_1x0top2_iet':None,
                        'area_1xm2to0_iet_minus_area_1x0top2_iet_mu_chg_pos':None,
                        'area_1xm2to0_iet_minus_area_1x0top2_iet_mu_chg_neg':None,
                        'area_3xm3to0_iet_over_area_3x0top3_iet':None,
                        'area_3xm3to0_iet_over_area_3x0top3_iet_mu_chg_pos':None,
                        'area_3xm3to0_iet_over_area_3x0top3_iet_mu_chg_neg':None,
                        'area_3xm3to0_iet_minus_area_3x0top3_iet':None,
                        'area_3xm3to0_iet_minus_area_3x0top3_iet_mu_chg_pos':None,
                        'area_3xm3to0_iet_minus_area_3x0top3_iet_mu_chg_neg':None,
                        'area_3xm3to0_iet_minus_area_3x0top3_iet_over_area_3xm3to0_iet_plus_area_3x0top3_iet':None,
                        'area_3xm3to0_iet_minus_area_3x0top3_iet_over_area_3xm3to0_iet_plus_area_3x0top3_iet_times_mu_chg':None,
                        'area_3xm3to0_iet_minus_area_3x0top3_iet_over_area_3xm3to0_iet_plus_area_3x0top3_iet_mu_chg_pos':None,
                        'area_3xm3to0_iet_minus_area_3x0top3_iet_over_area_3xm3to0_iet_plus_area_3x0top3_iet_mu_chg_neg':None,
                        'area_3xm7to0_iet_minus_area_3x0top7_iet':None,
                        'area_3xm7to0_iet_minus_area_3x0top7_iet_mu_chg_pos':None,
                        'area_3xm7to0_iet_minus_area_3x0top7_iet_mu_chg_neg':None,
                        'area_3xm7to0_iet_minus_area_3x0top7_iet_over_area_3xm7to0_iet_plus_area_3x0top7_iet':None,
                        'area_3xm7to0_iet_minus_area_3x0top7_iet_over_area_3xm7to0_iet_plus_area_3x0top7_iet_times_mu_chg':None,
                        'area_3xm7to0_iet_minus_area_3x0top7_iet_over_area_3xm7to0_iet_plus_area_3x0top7_iet_mu_chg_pos':None,
                        'area_3xm7to0_iet_minus_area_3x0top7_iet_over_area_3xm7to0_iet_plus_area_3x0top7_iet_mu_chg_neg':None,
                        'area_11x11-1x1_over_area_11x11_iet':None,
                        'area_11x11-1x3_over_area_11x11_iet':None,
                        'area_11x11-1x5_over_area_11x11_iet':None,
                        'area_11x11-3x3_over_area_11x11_iet':None,
                        'twobytwo_area_1x1_iet':None,
                        'twobytwo_area_1xm1to0_iet':None,
                        'twobytwo_area_1x0top1_iet':None,
                        'twobytwo_area_5x5_iet':None,
                        'twobytwo_area_5x5-1x1_iet':None,
                        'twobytwo_area_1xm1to0_iet_minus_twobytwo_area_1x0top1_iet':None,
                        'twobytwo_area_1xm1to0_iet_minus_twobytwo_area_1x0top1_iet_mu_chg_pos':None,
                        'twobytwo_area_1xm1to0_iet_minus_twobytwo_area_1x0top1_iet_mu_chg_neg':None,
                        'twobytwo_area_5x5-1x1_over_twobytwo_area_5x5_iet':None,
                        'twobytwo_area_5x5_iet_over_mu_ipt':None,
                        'twobytwo_area_5x5-1x1_iet_over_mu_ipt':None,
                        'twobytwo_area_5x5_5bit_iet_over_mu_ipt':None,
                        'twobytwo_area_5x5-1x1_5bit_iet_over_mu_ipt':None}

    x_vars_bins_2d_mu = [['iet_ietarel', 165, -82, 83],
                         ['iet_ietarel_red', 31, -15, 16],
                         ['area_11x11_iet', 256, 0, 256],
                         ['area_11x11_iet', 256, 0, 256],
                         ['area_11x11_iet', 256, 0, 256],
                         ['area_11x11_iet', 256, 0, 256],
                         ['area_11x11-1x1_iet', 256, 0, 256],
                         ['area_11x11-1x3_iet', 256, 0, 256],
                         ['area_11x11-1x5_iet', 256, 0, 256],
                         ['area_11x11-3x3_iet', 256, 0, 256],
                         ['mu_pt', -1]+pt_bins,
                         ['mu_pt', -1]+pt_bins,
                         ['mu_pt', -1]+pt_bins,
                         ['mu_pt', -1]+pt_bins,
                         ['mu_pt', -1]+pt_bins,
                         ['mu_pt', -1]+pt_bins,
                         ['mu_pt', -1]+pt_bins,
                         ['mu_pt', -1]+pt_bins,
                         ['mu_pt', -1]+pt_bins,
                         ['mu_pt', -1]+pt_bins]
    y_vars_bins_2d_mu = [['iet_iphirel', 73, -36, 37],
                         ['iet_iphirel_red', 31, -15, 16],
                         ['area_1x1_iet', 256, 0, 256],
                         ['area_1x3_iet', 256, 0, 256],
                         ['area_1x5_iet', 256, 0, 256],
                         ['area_3x3_iet', 256, 0, 256],
                         ['area_1x1_iet', 256, 0, 256],
                         ['area_1x3_iet', 256, 0, 256],
                         ['area_1x5_iet', 256, 0, 256],
                         ['area_3x3_iet', 256, 0, 256],
                         ['area_1x1_iet', 256, 0, 256],
                         ['area_1x3_iet', 256, 0, 256],
                         ['area_1x5_iet', 256, 0, 256],
                         ['area_3x3_iet', 256, 0, 256],
                         ['area_1xm2to0_iet_minus_area_1x0top2_iet', 511, -255, 256],
                         ['area_1xm2to0_iet_minus_area_1x0top2_iet_mu_chg_pos', 511, -255, 256],
                         ['area_1xm2to0_iet_minus_area_1x0top2_iet_mu_chg_neg', 511, -255, 256],
                         ['twobytwo_area_1x1_iet', 256, 0, 256],
                         ['twobytwo_area_5x5_iet', 256, 0, 256],
                         ['twobytwo_area_5x5-1x1_iet', 256, 0, 256]]

    x_title_vars_2d_mu = {'iet_ietarel':'i#eta_{tower} - i#eta_{#mu}',
                          'iet_ietarel_red':'i#eta_{tower} - i#eta_{#mu}',
                          'area_11x11_iet':'iE_{T}^{11x11}',
                          'area_11x11_iet':'iE_{T}^{11x11}',
                          'area_11x11_iet':'iE_{T}^{11x11}',
                          'area_11x11_iet':'iE_{T}^{11x11}',
                          'area_11x11-1x1_iet':'iE_{T}^{11x11-1x1}',
                          'area_11x11-1x3_iet':'iE_{T}^{11x11-1x3}',
                          'area_11x11-1x5_iet':'iE_{T}^{11x11-1x5}',
                          'area_11x11-3x3_iet':'iE_{T}^{11x11-3x3}',
                          'mu_pt':'p_{T}',
                          'mu_pt':'p_{T}',
                          'mu_pt':'p_{T}',
                          'mu_pt':'p_{T}',
                          'mu_pt':'p_{T}',
                          'mu_pt':'p_{T}',
                          'mu_pt':'p_{T}',
                          'mu_pt':'p_{T}',
                          'mu_pt':'p_{T}',
                          'mu_pt':'p_{T}'}
    y_title_vars_2d_mu = {'iet_iphirel':'i#phi_{tower} - i#phi_{#mu}',
                          'iet_iphirel_red':'i#phi_{tower} - i#phi_{#mu}',
                          'area_1x1_iet':'iE_{T}^{1x1}',
                          'area_1x3_iet':'iE_{T}^{1x3}',
                          'area_1x5_iet':'iE_{T}^{1x5}',
                          'area_3x3_iet':'iE_{T}^{3x3}',
                          'area_1x1_iet':'iE_{T}^{1x1}',
                          'area_1x3_iet':'iE_{T}^{1x3}',
                          'area_1x5_iet':'iE_{T}^{1x5}',
                          'area_3x3_iet':'iE_{T}^{3x3}',
                          'area_1x1_iet':'iE_{T}^{1x1}',
                          'area_1x3_iet':'iE_{T}^{1x3}',
                          'area_1x5_iet':'iE_{T}^{1x5}',
                          'area_3x3_iet':'iE_{T}^{3x3}',
                          'area_1xm2to0_iet_minus_area_1x0top2_iet':'iE_{T}^{1x-2to0} - iE_{T}^{1x0to+2}',
                          'area_1xm2to0_iet_minus_area_1x0top2_iet_mu_chg_pos':'positive #mu charge iE_{T}^{1x-2to0} - iE_{T}^{1x0to+2}',
                          'area_1xm2to0_iet_minus_area_1x0top2_iet_mu_chg_neg':'negative #mu charge iE_{T}^{1x-2to0} - iE_{T}^{1x0to+2}',
                          'twobytwo_area_1x1_iet':'2x2 area iE_{T}^{1x1}',
                          'twobytwo_area_5x5_iet':'2x2 area iE_{T}^{5x5}',
                          'twobytwo_area_5x5-1x1_iet':'2x2 area iE_{T}^{5x5-1x1}'}

    x_title_units_2d_mu = {'iet_ietarel':None,
                           'iet_ietarel_red':None,
                           'area_11x11_iet':None,
                           'area_11x11_iet':None,
                           'area_11x11_iet':None,
                           'area_11x11_iet':None,
                           'area_11x11-1x1_iet':None,
                           'area_11x11-1x3_iet':None,
                           'area_11x11-1x5_iet':None,
                           'area_11x11-3x3_iet':None,
                           'mu_pt':'GeV',
                           'mu_pt':'GeV',
                           'mu_pt':'GeV',
                           'mu_pt':'GeV',
                           'mu_pt':'GeV',
                           'mu_pt':'GeV',
                           'mu_pt':'GeV',
                           'mu_pt':'GeV',
                           'mu_pt':'GeV',
                           'mu_pt':'GeV'}
    y_title_units_2d_mu = {'iet_iphirel':None,
                           'iet_iphirel_red':None,
                           'area_1x1_iet':None,
                           'area_1x3_iet':None,
                           'area_1x5_iet':None,
                           'area_3x3_iet':None,
                           'area_1x1_iet':None,
                           'area_1x3_iet':None,
                           'area_1x5_iet':None,
                           'area_3x3_iet':None,
                           'area_1x1_iet':None,
                           'area_1x3_iet':None,
                           'area_1x5_iet':None,
                           'area_3x3_iet':None,
                           'area_1xm2to0_iet_minus_area_1x0top2_iet':None,
                           'area_1xm2to0_iet_minus_area_1x0top2_iet_mu_chg_pos':None,
                           'area_1xm2to0_iet_minus_area_1x0top2_iet_mu_chg_neg':None,
                           'twobytwo_area_1x1_iet':None,
                           'twobytwo_area_5x5_iet':None,
                           'twobytwo_area_5x5-1x1_iet':None}

    # the 2d histograms that should be tprofiles
    #profile_2d = {'2d_caloTower.iet_ieta_iet_iphi':True, '2d_caloTower.iet_ietarel_iet_iphirel':True, '2d_caloTower.iet_ietarel_red_iet_iphirel_red':True}

    varnames = []
    binnings = {}
    varnames2d = []
    binnings2d = {}

    for var_bin in vars_bins:
        varnames.append(namePrefix+'l1_caloTower.{var}'.format(var=var_bin[0]))
        binnings[namePrefix+'l1_caloTower.{var}'.format(var=var_bin[0])] = var_bin[1:]+[x_title_vars[var_bin[0]], x_title_units[var_bin[0]]]

    for var_bin_2d_x, var_bin_2d_y in zip(x_vars_bins_2d, y_vars_bins_2d):
        varnames2d.append(namePrefix+'2d_caloTower.{xvar}_{yvar}'.format(xvar=var_bin_2d_x[0], yvar=var_bin_2d_y[0]))
        binnings2d[namePrefix+'2d_caloTower.{xvar}_{yvar}'.format(xvar=var_bin_2d_x[0], yvar=var_bin_2d_y[0])] = [var_bin_2d_x[1:]+[x_title_vars_2d[var_bin_2d_x[0]], x_title_units_2d[var_bin_2d_x[0]]], var_bin_2d_y[1:]+[y_title_vars_2d[var_bin_2d_y[0]], y_title_units_2d[var_bin_2d_y[0]]]]

    for l1PtMin in l1PtMins:
        ptMinStr = '_l1ptmin{ptmin}'.format(ptmin=l1PtMin)

        for var_bin_mu in vars_bins_mu:
            varnames.append(namePrefix+'l1_caloTower'+ptMinStr+'.{var}'.format(var=var_bin_mu[0]))
            binnings[namePrefix+'l1_caloTower'+ptMinStr+'.{var}'.format(var=var_bin_mu[0])] = var_bin_mu[1:]+[x_title_vars_mu[var_bin_mu[0]], x_title_units_mu[var_bin_mu[0]]]

        for var_bin_2d_mu_x, var_bin_2d_mu_y in zip(x_vars_bins_2d_mu, y_vars_bins_2d_mu):
            varnames2d.append(namePrefix+'2d_caloTower'+ptMinStr+'.{xvar}_{yvar}'.format(xvar=var_bin_2d_mu_x[0], yvar=var_bin_2d_mu_y[0]))
            binnings2d[namePrefix+'2d_caloTower'+ptMinStr+'.{xvar}_{yvar}'.format(xvar=var_bin_2d_mu_x[0], yvar=var_bin_2d_mu_y[0])] = [var_bin_2d_mu_x[1:]+[x_title_vars_2d_mu[var_bin_2d_mu_x[0]], x_title_units_2d_mu[var_bin_2d_mu_x[0]]], var_bin_2d_mu_y[1:]+[y_title_vars_2d_mu[var_bin_2d_mu_y[0]], y_title_units_2d_mu[var_bin_2d_mu_y[0]]]]

    return HistManager(list(set(varnames)), binnings), HistManager2d(list(set(varnames2d)), binnings2d)

def analyse(evt, hm, hm2d, l1PtMins, emul=False):
    if emul:
        l1MuColl = evt.upgradeEmu
        l1CaloTwrColl = evt.caloTowersEmu
        emuPrefix = 'emu_'
    else:
        l1MuColl = evt.upgrade
        l1CaloTwrColl = evt.caloTowers
        emuPrefix = ''

    if useRecoMuon:
        recoMuonColl = evt.recoMuon

    bx_min = 0
    bx_max = 0

    # calo tower histograms
    histoprefix = emuPrefix+'l1_caloTower'
    histoprefix2d = emuPrefix+'2d_caloTower'

    nCaloTwr = l1CaloTwrColl.nTower
    hm.fill(histoprefix+'.n', nCaloTwr)
    for i in range(nCaloTwr):
        hm.fill(histoprefix+'.iet', l1CaloTwrColl.iet[i])
        hm.fill(histoprefix+'.ieta', l1CaloTwrColl.ieta[i])
        hm.fill(histoprefix+'.iphi', l1CaloTwrColl.iphi[i])
        hm.fill(histoprefix+'.iqual', l1CaloTwrColl.iqual[i])

        hm2d.fill(histoprefix2d+'.ieta_iphi', l1CaloTwrColl.ieta[i], l1CaloTwrColl.iphi[i])
        hm2d.fill(histoprefix2d+'.iet_ieta_iet_iphi', l1CaloTwrColl.ieta[i], l1CaloTwrColl.iphi[i], l1CaloTwrColl.iet[i])

    if useRecoMuon:
        reco_muon_idcs = MuonSelections.select_reco_muons(recoMuonColl, pos_charge=pos_charge, neg_charge=neg_charge)
        l1_muon_idcs = MuonSelections.select_ugmt_muons(l1MuColl, pt_min=0, bx_min=bx_min, bx_max=bx_max, pos_charge=pos_charge, neg_charge=neg_charge, tftype=tftype, useVtxExtraCoord=useVtxExtraCoord)
        l1matched_reco_muon_idcs = []
        for l1_idx in l1_muon_idcs:
            matched_idcs = [idx for idx in reco_muon_idcs if Matcher.delta_r(recoMuonColl.phiSt2[idx], recoMuonColl.etaSt2[idx], l1MuColl.muonPhi[l1_idx], l1MuColl.muonEta[l1_idx]) < 0.5]
            l1matched_reco_muon_idcs += matched_idcs
        l1_muon_idcs = list(set(l1matched_reco_muon_idcs))
    else:
        l1_muon_idcs = MuonSelections.select_ugmt_muons(l1MuColl, pt_min=0, bx_min=bx_min, bx_max=bx_max, pos_charge=pos_charge, neg_charge=neg_charge, tftype=tftype, useVtxExtraCoord=useVtxExtraCoord)

    for idx in l1_muon_idcs:
        if useRecoMuon:
            muIEta = int(91.954 * recoMuonColl.eta[idx])
            muIPhi = int(91.673 * recoMuonColl.phi[idx])
            if muIPhi < 0:
                muIPhi += 576
            elif muIPhi >= 576:
                muIPhi -= 576
        else:
            muIEta = l1MuColl.muonIEta[idx]
            muIPhi = l1MuColl.muonIPhi[idx]
            if useVtxExtraCoord:
                muIEta += l1MuColl.muonIDEta[idx]
                muIPhi += l1MuColl.muonIDPhi[idx]
                if muIPhi < 0:
                    muIPhi += 576
                elif muIPhi >= 576:
                    muIPhi -= 576

        muInCaloTowerIEta = CaloTowerIsolator.calc_muon_calo_tower_ieta(muIEta)
        muInCaloTowerIPhi = CaloTowerIsolator.calc_muon_calo_tower_iphi(muIPhi)

        # get towers relative to muon position and energy sums around
        relTwrs, iEtSums, iEmSums, iHadSums = CaloTowerIsolator.calc_calo_tower_sums(l1CaloTwrColl, muInCaloTowerIEta, muInCaloTowerIPhi, [(0, 0, 0, 0), (0, 0, -1, 1), (0, 0, -2, 2), (0, 0, -2, 0), (0, 0, 0, 2), (-1, 1, -1, 1), (-1, 1, -3, 0), (-1, 1, 0, 3), (-1, 1, -7, 0), (-1, 1, 0, 7), (-5, 5, -5, 5)])
        area_1x1_iet = iEtSums[0]
        area_1x3_iet = iEtSums[1]
        area_1x5_iet = iEtSums[2]
        area_1xm2to0_iet = iEtSums[3]
        area_1x0top2_iet = iEtSums[4]
        area_3x3_iet = iEtSums[5]
        area_3xm3to0_iet = iEtSums[6]
        area_3x0top3_iet = iEtSums[7]
        area_3xm7to0_iet = iEtSums[8]
        area_3x0top7_iet = iEtSums[9]
        area_11x11_iet = iEtSums[10]
        # energy sums of 2x2 towers
        iEtSums2x2, iEmSums2x2, iHadSums2x2 = CaloTowerIsolator.calc_calo_tower_2x2_sums(l1CaloTwrColl, muInCaloTowerIEta, muInCaloTowerIPhi, [(0, 0, 0, 0), (0, 0, -1, 0), (0, 0, 0, 1), (-2, 2, -2, 2)])
        twobytwo_area_1x1_iet = iEtSums2x2[0]
        twobytwo_area_1xm1to0_iet = iEtSums2x2[1]
        twobytwo_area_1x0top1_iet = iEtSums2x2[2]
        twobytwo_area_5x5_iet = iEtSums2x2[3]
        twobytwo_area_5x5_minus_1x1_iet = twobytwo_area_5x5_iet - twobytwo_area_1x1_iet

        twobytwo_area_5x5_5bit_iet = twobytwo_area_5x5_iet
        if twobytwo_area_5x5_5bit_iet > 31:
            twobytwo_area_5x5_5bit_iet = 31
        twobytwo_area_5x5_minus_1x1_5bit_iet = twobytwo_area_5x5_minus_1x1_iet
        if twobytwo_area_5x5_minus_1x1_5bit_iet > 31:
            twobytwo_area_5x5_minus_1x1_5bit_iet = 31
        #print '{i} {t} {i2} {t2}'.format(i=area_3x3_iet, t=area_11x11_iet, i2=twobytwo_area_1x1_iet, t2=twobytwo_area_5x5_iet)

        for l1PtMin in l1PtMins:
            # calo towers around L1 muons
            if useRecoMuon:
                ptmin_l1_muon_idcs = MuonSelections.select_reco_muons(recoMuonColl, pt_min=l1PtMin, idcs=l1_muon_idcs)
            else:
                ptmin_l1_muon_idcs = MuonSelections.select_ugmt_muons(l1MuColl, pt_min=l1PtMin, idcs=l1_muon_idcs)

            if idx not in ptmin_l1_muon_idcs:
                continue

            ptMinStr = '_l1ptmin{ptmin}'.format(ptmin=l1PtMin)

            if useRecoMuon:
                muEt = recoMuonColl.pt[idx]
                muIEt = 1 + 2 * int(muEt)
                muChg = recoMuonColl.charge[idx]
            else:
                muEt = l1MuColl.muonEt[idx]
                muIEt = l1MuColl.muonIEt[idx]
                muChg = l1MuColl.muonChg[idx]

            hm.fill(histoprefix+ptMinStr+'.n_mu', len(ptmin_l1_muon_idcs))
            hm.fill(histoprefix+ptMinStr+'.mu_pt', muEt)

            for relTwr in relTwrs:
                hm2d.fill(histoprefix2d+ptMinStr+'.iet_ietarel_iet_iphirel', relTwr[0], relTwr[1], relTwr[2])
                if abs(relTwr[0]) < 16 and abs(relTwr[1]) < 16:
                    hm2d.fill(histoprefix2d+ptMinStr+'.iet_ietarel_red_iet_iphirel_red', relTwr[0], relTwr[1], relTwr[2])

            hm.fill(histoprefix+ptMinStr+'.area_1x1_iet', area_1x1_iet)
            hm.fill(histoprefix+ptMinStr+'.area_1x3_iet', area_1x3_iet)
            hm.fill(histoprefix+ptMinStr+'.area_1x5_iet', area_1x5_iet)
            hm.fill(histoprefix+ptMinStr+'.area_1xm2to0_iet', area_1xm2to0_iet)
            hm.fill(histoprefix+ptMinStr+'.area_1x0top2_iet', area_1x0top2_iet)
            hm.fill(histoprefix+ptMinStr+'.area_3x3_iet', area_3x3_iet)
            hm.fill(histoprefix+ptMinStr+'.area_3xm3to0_iet', area_3xm3to0_iet)
            hm.fill(histoprefix+ptMinStr+'.area_3x0top3_iet', area_3x0top3_iet)
            hm.fill(histoprefix+ptMinStr+'.area_3xm7to0_iet', area_3xm7to0_iet)
            hm.fill(histoprefix+ptMinStr+'.area_3x0top7_iet', area_3x0top7_iet)
            hm.fill(histoprefix+ptMinStr+'.area_11x11_iet', area_11x11_iet)
            hm.fill(histoprefix+ptMinStr+'.area_11x11-1x1_iet', area_11x11_iet - area_1x1_iet)
            hm.fill(histoprefix+ptMinStr+'.area_11x11-1x3_iet', area_11x11_iet - area_1x3_iet)
            hm.fill(histoprefix+ptMinStr+'.area_11x11-1x5_iet', area_11x11_iet - area_1x5_iet)
            hm.fill(histoprefix+ptMinStr+'.area_11x11-3x3_iet', area_11x11_iet - area_3x3_iet)
            if area_11x11_iet > 0:
                hm.fill(histoprefix+ptMinStr+'.area_11x11-1x1_over_area_11x11_iet', (area_11x11_iet - area_1x1_iet) / float(area_11x11_iet))
                hm.fill(histoprefix+ptMinStr+'.area_11x11-1x3_over_area_11x11_iet', (area_11x11_iet - area_1x3_iet) / float(area_11x11_iet))
                hm.fill(histoprefix+ptMinStr+'.area_11x11-1x5_over_area_11x11_iet', (area_11x11_iet - area_1x5_iet) / float(area_11x11_iet))
                hm.fill(histoprefix+ptMinStr+'.area_11x11-3x3_over_area_11x11_iet', (area_11x11_iet - area_3x3_iet) / float(area_11x11_iet))

            hm.fill(histoprefix+ptMinStr+'.twobytwo_area_1x1_iet', twobytwo_area_1x1_iet)
            hm.fill(histoprefix+ptMinStr+'.twobytwo_area_1xm1to0_iet', twobytwo_area_1xm1to0_iet)
            hm.fill(histoprefix+ptMinStr+'.twobytwo_area_1x0top1_iet', twobytwo_area_1x0top1_iet)
            hm.fill(histoprefix+ptMinStr+'.twobytwo_area_5x5_iet', twobytwo_area_5x5_iet)
            hm.fill(histoprefix+ptMinStr+'.twobytwo_area_5x5-1x1_iet', twobytwo_area_5x5_minus_1x1_iet)
            if twobytwo_area_5x5_iet > 0:
                hm.fill(histoprefix+ptMinStr+'.twobytwo_area_5x5-1x1_over_twobytwo_area_5x5_iet', twobytwo_area_5x5_minus_1x1_iet / float(twobytwo_area_5x5_iet))
            if muIEt > 0:
                hm.fill(histoprefix+ptMinStr+'.twobytwo_area_5x5_iet_over_mu_ipt', twobytwo_area_5x5_iet / float(muIEt))
                hm.fill(histoprefix+ptMinStr+'.twobytwo_area_5x5_5bit_iet_over_mu_ipt', twobytwo_area_5x5_5bit_iet / float(muIEt))
                hm.fill(histoprefix+ptMinStr+'.twobytwo_area_5x5-1x1_iet_over_mu_ipt', twobytwo_area_5x5_minus_1x1_iet / float(muIEt))
                hm.fill(histoprefix+ptMinStr+'.twobytwo_area_5x5-1x1_5bit_iet_over_mu_ipt', twobytwo_area_5x5_minus_1x1_5bit_iet / float(muIEt))

            hm.fill(histoprefix+ptMinStr+'.area_1xm2to0_iet_minus_area_1x0top2_iet', area_1xm2to0_iet - area_1x0top2_iet)
            hm.fill(histoprefix+ptMinStr+'.area_3xm3to0_iet_minus_area_3x0top3_iet', area_3xm3to0_iet - area_3x0top3_iet)
            hm.fill(histoprefix+ptMinStr+'.area_3xm7to0_iet_minus_area_3x0top7_iet', area_3xm7to0_iet - area_3x0top7_iet)
            if area_3x0top3_iet != 0:
                hm.fill(histoprefix+ptMinStr+'.area_3xm3to0_iet_over_area_3x0top3_iet', area_3xm3to0_iet / float(area_3x0top3_iet))
            if area_3xm3to0_iet > 0 or area_3x0top3_iet > 0:
                hm.fill(histoprefix+ptMinStr+'.area_3xm3to0_iet_minus_area_3x0top3_iet_over_area_3xm3to0_iet_plus_area_3x0top3_iet', (area_3xm3to0_iet - area_3x0top3_iet) / float(area_3xm3to0_iet + area_3x0top3_iet))
                hm.fill(histoprefix+ptMinStr+'.area_3xm3to0_iet_minus_area_3x0top3_iet_over_area_3xm3to0_iet_plus_area_3x0top3_iet_times_mu_chg', muChg * (area_3xm3to0_iet - area_3x0top3_iet) / float(area_3xm3to0_iet + area_3x0top3_iet))
            if area_3xm7to0_iet > 0 or area_3x0top7_iet > 0:
                hm.fill(histoprefix+ptMinStr+'.area_3xm7to0_iet_minus_area_3x0top7_iet_over_area_3xm7to0_iet_plus_area_3x0top7_iet', (area_3xm7to0_iet - area_3x0top7_iet) / float(area_3xm7to0_iet + area_3x0top7_iet))
                hm.fill(histoprefix+ptMinStr+'.area_3xm7to0_iet_minus_area_3x0top7_iet_over_area_3xm7to0_iet_plus_area_3x0top7_iet_times_mu_chg', muChg * (area_3xm7to0_iet - area_3x0top7_iet) / float(area_3xm7to0_iet + area_3x0top7_iet))
            hm.fill(histoprefix+ptMinStr+'.twobytwo_area_1xm1to0_iet_minus_twobytwo_area_1x0top1_iet', twobytwo_area_1xm1to0_iet - twobytwo_area_1x0top1_iet)
            if muChg > 0:
                hm.fill(histoprefix+ptMinStr+'.area_1xm2to0_iet_minus_area_1x0top2_iet_mu_chg_pos', area_1xm2to0_iet - area_1x0top2_iet)
                hm.fill(histoprefix+ptMinStr+'.area_3xm3to0_iet_minus_area_3x0top3_iet_mu_chg_pos', area_3xm3to0_iet - area_3x0top3_iet)
                hm.fill(histoprefix+ptMinStr+'.area_3xm7to0_iet_minus_area_3x0top7_iet_mu_chg_pos', area_3xm7to0_iet - area_3x0top7_iet)
                if area_3x0top3_iet != 0:
                    hm.fill(histoprefix+ptMinStr+'.area_3xm3to0_iet_over_area_3x0top3_iet_mu_chg_pos', area_3xm3to0_iet / float(area_3x0top3_iet))
                if area_3xm3to0_iet > 0 or area_3x0top3_iet > 0:
                    hm.fill(histoprefix+ptMinStr+'.area_3xm3to0_iet_minus_area_3x0top3_iet_over_area_3xm3to0_iet_plus_area_3x0top3_iet_mu_chg_pos', (area_3xm3to0_iet - area_3x0top3_iet) / float(area_3xm3to0_iet + area_3x0top3_iet))
                if area_3xm7to0_iet > 0 or area_3x0top7_iet > 0:
                    hm.fill(histoprefix+ptMinStr+'.area_3xm7to0_iet_minus_area_3x0top7_iet_over_area_3xm7to0_iet_plus_area_3x0top7_iet_mu_chg_pos', (area_3xm7to0_iet - area_3x0top7_iet) / float(area_3xm7to0_iet + area_3x0top7_iet))
                hm.fill(histoprefix+ptMinStr+'.twobytwo_area_1xm1to0_iet_minus_twobytwo_area_1x0top1_iet_mu_chg_pos', twobytwo_area_1xm1to0_iet - twobytwo_area_1x0top1_iet)
            else:
                hm.fill(histoprefix+ptMinStr+'.area_1xm2to0_iet_minus_area_1x0top2_iet_mu_chg_neg', area_1xm2to0_iet - area_1x0top2_iet)
                hm.fill(histoprefix+ptMinStr+'.area_3xm3to0_iet_minus_area_3x0top3_iet_mu_chg_neg', area_3xm3to0_iet - area_3x0top3_iet)
                hm.fill(histoprefix+ptMinStr+'.area_3xm7to0_iet_minus_area_3x0top7_iet_mu_chg_neg', area_3xm7to0_iet - area_3x0top7_iet)
                if area_3x0top3_iet != 0:
                    hm.fill(histoprefix+ptMinStr+'.area_3xm3to0_iet_over_area_3x0top3_iet_mu_chg_neg', area_3xm3to0_iet / float(area_3x0top3_iet))
                if area_3xm3to0_iet > 0 or area_3x0top3_iet > 0:
                    hm.fill(histoprefix+ptMinStr+'.area_3xm3to0_iet_minus_area_3x0top3_iet_over_area_3xm3to0_iet_plus_area_3x0top3_iet_mu_chg_neg', (area_3xm3to0_iet - area_3x0top3_iet) / float(area_3xm3to0_iet + area_3x0top3_iet))
                if area_3xm7to0_iet > 0 or area_3x0top7_iet > 0:
                    hm.fill(histoprefix+ptMinStr+'.area_3xm7to0_iet_minus_area_3x0top7_iet_over_area_3xm7to0_iet_plus_area_3x0top7_iet_mu_chg_neg', (area_3xm7to0_iet - area_3x0top7_iet) / float(area_3xm7to0_iet + area_3x0top7_iet))
                hm.fill(histoprefix+ptMinStr+'.twobytwo_area_1xm1to0_iet_minus_twobytwo_area_1x0top1_iet_mu_chg_neg', twobytwo_area_1xm1to0_iet - twobytwo_area_1x0top1_iet)

            hm2d.fill(histoprefix2d+ptMinStr+'.area_11x11_iet_area_1x1_iet', area_11x11_iet, area_1x1_iet)
            hm2d.fill(histoprefix2d+ptMinStr+'.area_11x11_iet_area_1x3_iet', area_11x11_iet, area_1x3_iet)
            hm2d.fill(histoprefix2d+ptMinStr+'.area_11x11_iet_area_1x5_iet', area_11x11_iet, area_1x5_iet)
            hm2d.fill(histoprefix2d+ptMinStr+'.area_11x11_iet_area_3x3_iet', area_11x11_iet, area_3x3_iet)
            hm2d.fill(histoprefix2d+ptMinStr+'.area_11x11-1x1_iet_area_1x1_iet', area_11x11_iet - area_1x1_iet, area_1x1_iet)
            hm2d.fill(histoprefix2d+ptMinStr+'.area_11x11-1x3_iet_area_1x3_iet', area_11x11_iet - area_1x3_iet, area_1x3_iet)
            hm2d.fill(histoprefix2d+ptMinStr+'.area_11x11-1x5_iet_area_1x5_iet', area_11x11_iet - area_1x5_iet, area_1x5_iet)
            hm2d.fill(histoprefix2d+ptMinStr+'.area_11x11-3x3_iet_area_3x3_iet', area_11x11_iet - area_3x3_iet, area_3x3_iet)
            if l1PtMin == 0:
                hm2d.fill(histoprefix2d+ptMinStr+'.mu_pt_area_1x1_iet', muEt, area_1x1_iet)
                hm2d.fill(histoprefix2d+ptMinStr+'.mu_pt_area_1x3_iet', muEt, area_1x3_iet)
                hm2d.fill(histoprefix2d+ptMinStr+'.mu_pt_area_1x5_iet', muEt, area_1x5_iet)
                hm2d.fill(histoprefix2d+ptMinStr+'.mu_pt_area_3x3_iet', muEt, area_3x3_iet)
                hm2d.fill(histoprefix2d+ptMinStr+'.mu_pt_area_1xm2to0_iet_minus_area_1x0top2_iet', muEt, area_1xm2to0_iet - area_1x0top2_iet)
                if muChg > 0:
                    hm2d.fill(histoprefix2d+ptMinStr+'.mu_pt_area_1xm2to0_iet_minus_area_1x0top2_iet_mu_chg_pos', muEt, area_1xm2to0_iet - area_1x0top2_iet)
                else:
                    hm2d.fill(histoprefix2d+ptMinStr+'.mu_pt_area_1xm2to0_iet_minus_area_1x0top2_iet_mu_chg_neg', muEt, area_1xm2to0_iet - area_1x0top2_iet)
                hm2d.fill(histoprefix2d+ptMinStr+'.mu_pt_twobytwo_area_1x1_iet', muEt, twobytwo_area_1x1_iet)
                hm2d.fill(histoprefix2d+ptMinStr+'.mu_pt_twobytwo_area_5x5_iet', muEt, twobytwo_area_5x5_iet)
                hm2d.fill(histoprefix2d+ptMinStr+'.mu_pt_twobytwo_area_5x5-1x1_iet', muEt, twobytwo_area_5x5_iet - twobytwo_area_1x1_iet)

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

    emul = opts.emul

    global useVtxExtraCoord
    useVtxExtraCoord = opts.extraCoord

    global useRecoMuon
    useRecoMuon = opts.use_reco

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

    global tftype
    if opts.tftype == 'bmtf':
        tftype = 0
    elif opts.tftype == 'omtf':
        tftype = 1
    elif opts.tftype == 'emtf':
        tftype = 2

    l1PtMins = [0, 5, 12, 20]
    #l1PtMins = [0, 5, 12, 20, 25]

    # book the histograms
    L1Ana.log.info("Booking combined run histograms.")
    hm, hm2d = book_histograms(l1PtMins, emul=emul)

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
            analyse(event, hm, hm2d, l1PtMins, emul=emul)
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
    muEtaScale = 0.010875
    useVtxExtraCoord = False
    useRecoMuon = False
    pos_charge = True
    neg_charge = True
    tftype = None
    saveHistos = True
    main()

