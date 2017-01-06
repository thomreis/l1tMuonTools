#!/usr/bin/env python
from ToolBox import parse_options_and_init_log
# have to do this first or ROOT masks the -h messages
opts, parser = parse_options_and_init_log()

from L1Analysis import L1Ana, L1Ntuple
from analysis_tools.plotting import HistManager, HistManager2d
from analysis_tools.plottools import *
import ROOT as root
import re
import os

def parse_options_plotRates(parser):
    """
    Adds often used options to the OptionParser...
    """
    parsers = parser.add_subparsers()
    sub_parser = parsers.add_parser("plotKinematics")
    sub_parser.add_argument("-i", "--interactive", dest="interactive", action='store_false', help="Draw plots on screen.")
    sub_parser.add_argument("--mukin", dest="mukin", default=False, action='store_true', help="Plot kinematics plots for the first three muons.")
    sub_parser.add_argument("--qstack", dest="qstack", default=False, action='store_true', help="Plot a quality stack plot for muon kinematics.")
    sub_parser.add_argument("--2d", dest="twod", default=False, action='store_true', help="Plot 2D L1 muon1 vs. muon2 plots.")
    sub_parser.add_argument("--public", dest="public", default=False, action='store_true', help="Plot style for publication.")

    opts, unknown = parser.parse_known_args()
    return opts


def plot_hists_standard(hm, hName, den=None, xTitle='', yTitle='# muons', threshold=False, stacked=False, normToBinWidth=False, xMax=None, reg='', scaleFactor=1., data=False):
    ugmt_dict = {'num':hName, 'den':den, 'drawopt':'', 'stacked':False, 'err':True}
    ugmt_dict.update(hist_style('data', lw=2))
    hDefs = []
    if reg == '':
        hDefs.append(ugmt_dict)

    if (hName.find('.pt') != -1 or hName.find('.n') != -1) and not den:
        logy = True
    else:
        logy = False

    xBase = 0.17
    yBase = 0.56
    notes = extract_notes_from_name(hName, xBase, yBase)
    if den:
        den_eta_number_strs = re.findall(r'[\d\.\d]+', den[den.find('EtaMin')+6:den.find('EtaMax')+12])
        if len(den_eta_number_strs) > 1 and eta_number_strs != den_eta_number_strs:
            den_note_str = den_eta_number_strs[0]+' < |#eta^{GMT}| < '+den_eta_number_strs[1]
            notes.append([xBase, yBase-0.05, den_note_str, True])

    prefix = 'c_'

    return plot_hists(hm, hDefs, xTitle, yTitle, threshold, normToBinWidth, False, xMax, prefix, notes, scaleFactor, False, logy, data)

def plot_hists_qstack(hm, hName, den=None, xTitle='', yTitle='# muons', threshold=False, stacked=False, normToBinWidth=False, xMax=None, reg='', scaleFactor=1., data=False):
    hDefs = []
    if stacked:
        stylePref = 'data_q'
        lw=1
    else:
        stylePref = 'data_qmin'
        lw=2
    if reg == '':
        for q in reversed(range(4,16,4)):
            ugmt_q_dict = {'num':hName.replace('XX', '{q}'.format(q=q)), 'den':den, 'drawopt':'hist', 'stacked':stacked, 'err':False}
            ugmt_q_dict.update(hist_style(stylePref+str(q), filled=stacked, lw=lw))
            hDefs.append(ugmt_q_dict)

    if (hName.find('.pt') != -1 or hName.find('.n') != -1) and not den:
        logy = True
    else:
        logy = False

    xBase = 0.07
    yBase = 0.56
    notes = extract_notes_from_name(hName, xBase, yBase, qualTxt=False)
    if den:
        den_eta_number_strs = re.findall(r'[\d\.\d]+', den[den.find('EtaMin')+6:den.find('EtaMax')+12])
        if len(den_eta_number_strs) > 1 and eta_number_strs != den_eta_number_strs:
            den_note_str = den_eta_number_strs[0]+' < |#eta^{GMT}| < '+den_eta_number_strs[1]
            notes.append([xBase, yBase-0.05, den_note_str, True])

    prefix = 'qualcomp_'

    return plot_hists(hm, hDefs, xTitle, yTitle, threshold, normToBinWidth, False, xMax, prefix, notes, scaleFactor, False, logy, data=data)

def plot_hists_muons(hm, hName, xTitle='', yTitle='# muons', threshold=False, stacked=False, normToBinWidth=False, normalise=False, xMax=None, reg='', scaleFactor=1., data=False):
    hDefs = []
    if reg == '':
        for i in reversed(range(1, 4)):
            ugmt_q_dict = {'num':hName.replace('muX', 'mu{i}'.format(i=i)), 'den':None, 'drawopt':'', 'stacked':False, 'err':True}
            ugmt_q_dict.update(hist_style('mu_{i}'.format(i=i), marker=True))
            hDefs.append(ugmt_q_dict)

    if hName[-2:] == 'pt':
        logy = True
    else:
        logy = False

    if normalise:
        yTitle='A.U.'

    xBase = 0.01
    yBase = 0.78
    notes = extract_notes_from_name(hName, xBase, yBase, qualTxt=True)

    prefix = 'mucomp_'

    return plot_hists(hm, hDefs, xTitle, yTitle, threshold, normToBinWidth, normalise, xMax, prefix, notes, scaleFactor, False, logy, data=data)

def main():
    opts = parse_options_plotRates(parser)
    batchRun = opts.interactive
    if batchRun:
        root.gROOT.SetBatch(True)


    quals = [4, 8, 12]
    isData = True

    hm = HistManager(filename=opts.fname, subdir='all_runs')

    L1Ana.init_l1_analysis()
    print ""

    # holds the canvases, histograms, etc.
    objects = []

    ##########################################################################
    # L1 muonu kinematic variables
    if opts.mukin:
        objects.append(plot_hists_muons(hm, 'l1_muon_qmin12.muXpt', xTitle='p_{T} (GeV/c)', yTitle='# muons/(GeV/c)', normToBinWidth=True, normalise=False, data=isData))
        objects.append(plot_hists_muons(hm, 'l1_muon_qmin12.muXeta', xTitle='#eta', normalise=True, data=isData))
        objects.append(plot_hists_muons(hm, 'l1_muon_qmin12.muXphi', xTitle='#phi', normalise=True, data=isData))
        #objects.append(plot_hists_muons(hm, 'l1_muon_qmin12.muXqual', xTitle='quality', normalise=True, data=isData))
        objects.append(plot_hists_muons(hm, 'l1_muon_qmin12.muXcharge', xTitle='charge', normalise=True, data=isData))
        objects.append(plot_hists_muons(hm, 'l1_muon_qmin12.muXtfMuonIdx', xTitle='TF muon index', normalise=True, data=isData))

        objects.append(plot_hists_muons(hm, 'l1_muon_qmin8.muXpt', xTitle='p_{T} (GeV/c)', yTitle='# muons/(GeV/c)', normToBinWidth=True, normalise=False, data=isData))
        objects.append(plot_hists_muons(hm, 'l1_muon_qmin8.muXeta', xTitle='#eta', normalise=True, data=isData))
        objects.append(plot_hists_muons(hm, 'l1_muon_qmin8.muXphi', xTitle='#phi', normalise=True, data=isData))
        objects.append(plot_hists_muons(hm, 'l1_muon_qmin8.muXqual', xTitle='quality', normalise=True, data=isData))
        objects.append(plot_hists_muons(hm, 'l1_muon_qmin8.muXcharge', xTitle='charge', normalise=True, data=isData))
        objects.append(plot_hists_muons(hm, 'l1_muon_qmin8.muXtfMuonIdx', xTitle='TF muon index', normalise=True, data=isData))

        objects.append(plot_hists_muons(hm, 'l1_muon_qmin4.muXpt', xTitle='p_{T} (GeV/c)', yTitle='# muons/(GeV/c)', normToBinWidth=True, normalise=False, data=isData))
        objects.append(plot_hists_muons(hm, 'l1_muon_qmin4.muXeta', xTitle='#eta', normalise=True, data=isData))
        objects.append(plot_hists_muons(hm, 'l1_muon_qmin4.muXphi', xTitle='#phi', normalise=True, data=isData))
        objects.append(plot_hists_muons(hm, 'l1_muon_qmin4.muXqual', xTitle='quality', normalise=True, data=isData))
        objects.append(plot_hists_muons(hm, 'l1_muon_qmin4.muXcharge', xTitle='charge', normalise=True, data=isData))
        objects.append(plot_hists_muons(hm, 'l1_muon_qmin4.muXtfMuonIdx', xTitle='TF muon index', normalise=True, data=isData))

    # quality stack
    if opts.qstack:
        objects.append(plot_hists_qstack(hm, 'l1_muon_qminXX.n', xTitle='#mu multiplicity', data=isData))
        objects.append(plot_hists_qstack(hm, 'l1_muon_qXX.pt', xTitle='p_{T} (GeV/c)', yTitle='# muons/(GeV/c)', normToBinWidth=True, stacked=True, data=isData))
        objects.append(plot_hists_qstack(hm, 'l1_muon_qXX.eta', stacked=True, data=isData))
        objects.append(plot_hists_qstack(hm, 'l1_muon_qXX.phi', stacked=True, data=isData))
        objects.append(plot_hists_qstack(hm, 'l1_muon_qXX.qual', stacked=True, data=isData))
        objects.append(plot_hists_qstack(hm, 'l1_muon_qXX.charge', stacked=True, data=isData))
        objects.append(plot_hists_qstack(hm, 'l1_muon_qXX.tfMuonIdx', stacked=True, data=isData))
        objects.append(plot_hists_qstack(hm, 'l1_muon_qXX.dpt', stacked=True, data=isData))
        objects.append(plot_hists_qstack(hm, 'l1_muon_qXX.dpt', stacked=True, data=isData, xMax=25))
        objects.append(plot_hists_qstack(hm, 'l1_muon_qXX.dptoverpt', stacked=True, data=isData))
        objects.append(plot_hists_qstack(hm, 'l1_muon_qXX.dr', stacked=True, data=isData))
        objects.append(plot_hists_qstack(hm, 'l1_muon_qXX.dr', stacked=True, data=isData, xMax=0.5))
        objects.append(plot_hists_qstack(hm, 'l1_muon_qXX.deta', stacked=True, data=isData))
        objects.append(plot_hists_qstack(hm, 'l1_muon_qXX.deta', stacked=True, data=isData, xMax=0.5))
        objects.append(plot_hists_qstack(hm, 'l1_muon_qXX.dphi', stacked=True, data=isData))
        objects.append(plot_hists_qstack(hm, 'l1_muon_qXX.dphi', stacked=True, data=isData, xMax=0.5))

    # 2d reco vs. L1 plots
    if opts.twod:
        hm2d = HistManager2d(filename=opts.fname, subdir='all_runs')

        for qual in quals:
            histoprefix2d = '2d_muon_qmin{q}'.format(q=qual)
            objects.append(plot_2dhist(hm2d, histoprefix2d+'.pt', drawDiag=False, data=isData))
            objects.append(plot_2dhist(hm2d, histoprefix2d+'.eta', drawDiag=False, data=isData))
            objects.append(plot_2dhist(hm2d, histoprefix2d+'.phi', drawDiag=False, data=isData))
            objects.append(plot_2dhist(hm2d, histoprefix2d+'.qual', drawDiag=False, data=isData))
            objects.append(plot_2dhist(hm2d, histoprefix2d+'.charge', drawDiag=False, data=isData))
            objects.append(plot_2dhist(hm2d, histoprefix2d+'.tfMuonIdx', drawDiag=False, data=isData))

    ##########################################################################
    # save plots to root file
    if savePlots:
        plotdir = 'plots_'+opts.fname.replace('.root','').partition('/')[0]
        if opts.public:
            plotdir += '_public'
        if not os.path.exists(plotdir):
            os.makedirs(plotdir)
        output = root.TFile('./'+plotdir+'/l1_muon_kinematic_plots.root', 'recreate')
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

