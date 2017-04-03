#!/usr/bin/env python
from ToolBox import parse_options_and_init_log
# have to do this first or ROOT masks the -h messages
opts, parser = parse_options_and_init_log()

from L1Analysis import L1Ana
from analysis_tools.plotting import HistManager, HistManager2d
from analysis_tools.plottools import *
import ROOT as root
import os

def parse_options_plotRates(parser):
    """
    Adds often used options to the OptionParser...
    """
    parsers = parser.add_subparsers()
    sub_parser = parsers.add_parser("plotExtrapolation")
    sub_parser.add_argument("-i", "--interactive", dest="interactive", action='store_false', help="Draw plots on screen.")
    sub_parser.add_argument("--delta", dest="delta", default=False, action='store_true', help="Plot Gen delta TProfile plots.")
    sub_parser.add_argument("--2d", dest="twod", default=False, action='store_true', help="Plot 2D L1 muon1 vs. muon2 plots.")
    sub_parser.add_argument("--extrapolated", dest="extrapolated", default=False, action='store_true', help="Plots from extrapolated coordinates.")
    sub_parser.add_argument("--public", dest="public", default=False, action='store_true', help="Plot style for publication.")

    opts, unknown = parser.parse_known_args()
    return opts

def plot_tprof_muons(hm, hName, xTitle='', yTitle='# muons', eta_ranges=[], threshold=False, stacked=False, normToBinWidth=False, normalise=False, xMax=None, reg='', scaleFactor=1., data=False, rebin=1):
    hDefs = []
    if reg == '':
        for i, eta_range in enumerate(eta_ranges):
            eta_min = eta_range[0]
            eta_max = eta_range[1]
            ugmt_dict = {'num':hName.replace('XXXX', '{etaMin}'.format(etaMin=eta_min)).replace('YYYY', '{etaMax}'.format(etaMax=eta_max)), 'den':None, 'drawopt':'', 'stacked':False, 'err':True}
            # set legend text
            style = hist_style('etarange_{i}'.format(i=i), marker=True)
            style['legtext'] = '{etaMin} #leq |#eta| < {etaMax}'.format(etaMin=eta_min, etaMax=eta_max)
            ugmt_dict.update(style)
            hDefs.append(ugmt_dict)

    logy = False

    if normalise:
        yTitle='A.U.'

    xBase = 0.01
    yBase = 0.78
    notes = extract_notes_from_name(hName, xBase, yBase, qualTxt=True)

    prefix = 'mucomp_'

    return plot_hists(hm, hDefs, xTitle, yTitle, threshold, normToBinWidth, normalise, xMax, prefix, notes, scaleFactor, False, logy, data=data, rebin=rebin)

def main():
    opts = parse_options_plotRates(parser)
    batchRun = opts.interactive
    if batchRun:
        root.gROOT.SetBatch(True)

    extrapol_str = ''
    if opts.extrapolated:
        extrapol_str = '_extrapol'

    isData = False 

#    eta_ranges = [[0, 0.83], [0.83, 1.24], [1.24, 2.4]]
#    eta_ranges = [[1.2, 1.55], [1.55, 1.85], [1.85, 2.4]]
    eta_ranges = [[0, 2.4], [0, 0.83], [0.83, 1.24], [1.24, 2.4], [1.2, 1.55], [1.55, 1.85], [1.85, 2.4]]

    hm = HistManager(filename=opts.fname, subdir='all_runs')

    L1Ana.init_l1_analysis()
    print ""

    # holds the canvases, histograms, etc.
    objects = []

    ##########################################################################
    #rebinPt = 1
    rebinPt = 2
    # L1 muonu kinematic variables
    if opts.delta:
        objects.append(plot_tprof_muons(hm, 'l1_muon'+extrapol_str+'_absEtaMinXXXX_absEtaMaxYYYY.pt_dpt', xTitle='p_{T}^{L1} (GeV/c)', yTitle='<|p_{T}^{L1} - p_{T}^{GEN}|>', eta_ranges=eta_ranges, data=isData, rebin=rebinPt))
        objects.append(plot_tprof_muons(hm, 'l1_muon'+extrapol_str+'_absEtaMinXXXX_absEtaMaxYYYY.pt_deta', xTitle='p_{T}^{L1} (GeV/c)', yTitle='<#eta^{L1} - #eta^{GEN}>', eta_ranges=eta_ranges, data=isData, xMax=80., rebin=rebinPt))
        objects.append(plot_tprof_muons(hm, 'l1_muon'+extrapol_str+'_absEtaMinXXXX_absEtaMaxYYYY.pt_dphi', xTitle='p_{T}^{L1} (GeV/c)', yTitle='<#phi^{L1} - #phi^{GEN}>', eta_ranges=eta_ranges, data=isData, xMax=80., rebin=rebinPt))
        objects.append(plot_tprof_muons(hm, 'l1_muon'+extrapol_str+'_absEtaMinXXXX_absEtaMaxYYYY.pt_absdeta', xTitle='p_{T}^{L1} (GeV/c)', yTitle='<|#eta^{L1} - #eta^{GEN}|>', eta_ranges=eta_ranges, data=isData, xMax=80., rebin=rebinPt))
        objects.append(plot_tprof_muons(hm, 'l1_muon'+extrapol_str+'_absEtaMinXXXX_absEtaMaxYYYY.pt_absdphi', xTitle='p_{T}^{L1} (GeV/c)', yTitle='<|#phi^{L1} - #phi^{GEN}|>', eta_ranges=eta_ranges, data=isData, xMax=80., rebin=rebinPt))

    # 2d reco vs. L1 plots
    if opts.twod:
        hm2d = HistManager2d(filename=opts.fname, subdir='all_runs')
        for eta_range in eta_ranges:
            eta_min = eta_range[0]
            eta_max = eta_range[1]
            histoprefix2d = '2d_muon'+extrapol_str+'_absEtaMin{etaMin}_absEtaMax{etaMax}'.format(etaMin=eta_min, etaMax=eta_max)

            objects.append(plot_2dhist(hm2d, histoprefix2d+'.pt_dcharge', drawDiag=False, data=isData))
            objects.append(plot_2dhist(hm2d, histoprefix2d+'.pt_deta', drawDiag=False, data=isData, xMax=80.))
            objects.append(plot_2dhist(hm2d, histoprefix2d+'.pt_dphi', drawDiag=False, data=isData, xMax=80.))

    ##########################################################################
    # save plots to root file
    if savePlots:
        plotdir = 'plots_'+opts.fname.replace('.root','').partition('/')[0]
        if opts.public:
            plotdir += '_public'
        if not os.path.exists(plotdir):
            os.makedirs(plotdir)
        fname_extra_str = ''
        if opts.extrapolated:
            fname_extra_str = '_extrapolated'
        output = root.TFile('./'+plotdir+'/l1_muon_gencomp'+fname_extra_str+'_plots.root', 'recreate')
        output.cd()
        for obj in objects:
            c = obj[0]

            c.Write(c.GetName())
            if opts.public:
                c.Print('./'+plotdir+'/'+c.GetName()+'.pdf', '.pdf')
            c.Print('./'+plotdir+'/'+c.GetName()+'.png', '.png')
        print 'get the plots with: scp -r lxplus:{pwd}/{plotdir}/ .'.format(pwd=os.getcwd(), plotdir=plotdir)
        output.Close()

    # wait
    if not batchRun:
        raw_input("Press ENTER to quit.")

if __name__ == "__main__":
    savePlots = True
    batchRun = True
    font = 42
    fontSize = 0.04
    main()

