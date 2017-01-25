#!/usr/bin/env python
from L1Analysis import L1Ana
from analysis_tools.plotting import HistManager
from analysis_tools.plottools import *
import ROOT as root
import os
import argparse

def parse_options():
    """
    Command line option parser
    """
    parser = argparse.ArgumentParser(description="Isolation ROC curve plot macro", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--feff", dest="feff", default="", type=str, help="Root file containing efficiency histograms.")
    parser.add_argument("--frate", dest="frate", default="", type=str, help="Root file containing rate histograms.")
    parser.add_argument("-i", "--interactive", dest="interactive", action='store_false', help="Draw plots on screen.")
    parser.add_argument("-b", "--bunches", dest="bunches", default=0, type=int, help="Number of colliding bunches")
    parser.add_argument("--pu", dest="pu", default=20, type=int, help="Average PU. default=20")
    parser.add_argument("--xsect", dest="xsect", default=80, type=float, help="Total cross section in mb. default=80 mb")
    parser.add_argument("--instlumi", dest="instlumi", default=1.2e34, type=float, help="Instantaneous luminosity. default=1.2e-34 cm-2s-1")
    parser.add_argument("--scale", dest="scale", default=1., type=float, help="Additional scale factor for rate calculation")
    parser.add_argument("-l", "--legacy", dest="legacy", action='store_true', help="Draw plots relative to legacy.")
    parser.add_argument("--public", dest="public", default=False, action='store_true', help="Plot style for publication.")
    parser.add_argument("--iso-method", dest="isomethod", type=str, default='abs', help="Isolation method. ['abs', 'rel', 'inner', 'outovertot', 'inner2x2', 'outovertot2x2']")

    opts, unknown = parser.parse_known_args()
    if opts.feff == "" or opts.frate == "":
        parser.print_help()
        exit(0)
    return opts

def get_eff(hm, hNames):
    varNames = hm.get_varnames()
    if hNames['num'] not in varNames:
        print 'Error: Numerator histogram '+hNames['num']+' not found.'
        return (0., 0., 0.)
    if hNames['den'] not in varNames:
        print 'Error: Denominator histogram '+hNames['den']+' not found.'
        return (0., 0., 0.)

    effHist = hm.get_efficiency(hNames['num'], hNames['den']).Clone()

    eff = effHist.GetEfficiency(1)
    effErrLow = effHist.GetEfficiencyErrorLow(1)
    effErrUp = effHist.GetEfficiencyErrorUp(1)
    return (eff, effErrLow, effErrUp)

def get_relative_eff(hm, hNames, reference):
    eff = get_eff(hm, hNames)
    relEff = 0.
    if reference > 0.:
        relEff = eff[0] / reference
    return relEff

def get_rate(hm, hName, threshold=0., scaleFactor=1.):
    if hName not in hm.get_varnames():
        print 'Error: Rate histogram '+hName+' not found.'
        return (0., 0.)

    h = hm.get_threshold_hist(hName).Clone()
    if scaleFactor != 1.:
        h.Scale(scaleFactor)
    binNr = h.FindBin(threshold)
    binCont = h.GetBinContent(binNr)
    binErr = h.GetBinError(binNr)
    return (binCont, binErr)

def get_relative_rate(hm, hName, reference, threshold=0., scaleFactor=1.):
    rate = get_rate(hm, hName, threshold, scaleFactor)
    relRate = 0.
    if reference > 0.:
        relRate = rate[0] / reference
    return relRate

def print_rates(hm, hName, scaleFactor=1.):
    print '===== Rates ====='
    print hName
    print ''
    print '\nThreshold              uGMT'
    for threshold in [0, 3, 5, 7, 10, 12, 14, 16, 18, 20, 22, 25, 30, 40, 50, 60]:
        rate = get_rate(hm, hName, threshold, scaleFactor)
        print '{threshold:>3} GeV:   {rate:>8.3f} +/- {err:>5.3f} kHz'.format(threshold=threshold, rate=rate[0], err=rate[1])

    print '================='

def print_efficiencies(hm, hNames):
    print '===== Efficiencies ====='
    print hNames['num']
    print hNames['den']
    print ''
    print '\nThreshold              uGMT'
    eff = get_eff(hm, hNames)
    print '{threshold:>3} GeV:   {eff:>8.3f} +{errUp:>5.3f} -{errLow:>5.3f}'.format(threshold=22, eff=eff[0], errUp=eff[2], errLow=eff[1])

    print '================='

def plot_roc_curves(thresholdAndPointss, qual):
    canvas_name = 'roc_curves_qualMin{q}'.format(q=qual)
    c = root.TCanvas(canvas_name, canvas_name, 100, 100, 600, 600)
    c.cd()
    set_root_style()

    colors = [root.kRed+1, root.kGreen+1, root.kCyan+1, root.kBlue+1, root.kMagenta+1]

    # setup legend according to how many graphs are in the plot
    legYmin = 0.5-0.04*len(thresholdAndPointss)
    legXmin = 0.68
    legXmax = 0.9
    legend = root.TLegend(legXmin, legYmin, legXmax, 0.5)
    legend.SetTextFont(font)
    legend.SetTextSize(fontSize)
    legend.SetBorderSize(0)
    legend.SetFillColor(19)
    legend.SetFillStyle(0)
    #legend.SetNColumns(2)
    legEntries = []

    # two points to mark the edges
    roc_curves = [root.TGraph(2)]
    roc_curves[-1].SetPoint(0, 0., 0.)
    roc_curves[-1].SetPoint(1, 1., 1.)
    roc_curves[-1].SetMarkerStyle(root.kDot)
    roc_curves[-1].Draw('AP')
 
    # actual ROC curves
    for i, thresholdAndPoints in enumerate(thresholdAndPointss):
        thr = thresholdAndPoints[0]
        roc_curves.append(root.TGraph(len(thresholdAndPoints[1])))
        roc_curves[-1].SetName('roc_graph_ptmin{thr}'.format(thr=thr))
        roc_curves[-1].SetMarkerStyle(root.kFullCircle)
        roc_curves[-1].SetMarkerColor(colors[i])
        roc_curves[-1].SetLineColor(colors[i])
        for j, p in enumerate(thresholdAndPoints[1]):
            roc_curves[-1].SetPoint(j, p[0], p[1])
        roc_curves[-1].Draw('LP')
        legend.AddEntry(roc_curves[-1], '{thr} GeV'.format(thr=thr), 'lp')

    legend.Draw('same')

    xAxis = roc_curves[0].GetXaxis()
    yAxis = roc_curves[0].GetYaxis()
    xAxis.SetTitle('relative rate')
    xAxis.SetTitleFont(font)
    xAxis.SetLabelFont(font)
    xAxis.SetLabelSize(fontSize)
    xAxis.SetRangeUser(0., 1.)
    yAxis.SetTitle('relative efficiency')
    yAxis.SetTitleOffset(1.5)
    yAxis.SetTitleFont(font)
    yAxis.SetLabelFont(font)
    yAxis.SetLabelSize(fontSize)
    yAxis.SetRangeUser(0., 1.)

    # diagonal line
    line = root.TLine(0., 0., 1., 1.)
    line.SetLineStyle(root.kDashed)
    line.SetLineColor(root.kMagenta)
    line.Draw('same')

    notes = [[0.17, 0.86, 'quality #geq {q}'.format(q=qual), False]]
    notes.append([0.53, 0.93, 'CMS internal, 13 TeV', False])
    text = add_text(notes)

    c.Modified()
    c.Update()

    return c, roc_curves, legend, line, text

def main():
    opts = parse_options()
    plotLegacy = opts.legacy
    batchRun = opts.interactive
    if batchRun:
        root.gROOT.SetBatch(True)

    hmEff = HistManager(filename=opts.feff, subdir='all_runs')
    hmRate = HistManager(filename=opts.frate)

    nEvtsAna = hmRate.get('n_evts_analysed').GetBinContent(1)
    print '{n} events have been analysed.'.format(n=nEvtsAna)

    # calculate the scale factor for rate in Hz
    orbitFreq = 11245.6
    nCollBunches = opts.bunches
    nZeroBiasEvents = nEvtsAna
    crossSect = opts.xsect
    instLumi = opts.instlumi
    pu = opts.pu
    thisIsData=True
    # determine that this is MC if there is no number of colliding bunches given (defaults to 0 then)
    if nCollBunches == 0:
        print "No number of colliding bunches given. Assuming this is MC"
        print "Using {instLumi} cm-2s-1 as instantaneous luminosity, {crossSect} mb as cross section, and {pu} as average number of pileup to determine number of colliding bunches.".format(instLumi=instLumi, crossSect=crossSect, pu=pu)
        nCollBunches = round(instLumi * crossSect*1e-27 / (pu * orbitFreq))
        thisIsData=False
    else:
        print "Assuming this is data"
    convFactorToHz = orbitFreq * nCollBunches / nZeroBiasEvents
    print 'Conversion factor to rate in Hz with {orbitFreq} Hz orbit frequency, {nCollBunches} colliding bunches and {nZeroBiasEvents} analyzed zero bias events: {convFactorToHz}'.format(orbitFreq=orbitFreq, nCollBunches=nCollBunches, nZeroBiasEvents=nZeroBiasEvents, convFactorToHz=convFactorToHz)
    if opts.scale != 1.:
        convFactorToHz *= opts.scale
        print 'Conversion factor after applying additinoal scale factor of {sf}: {convFactorToHz}'.format(sf=opts.scale, convFactorToHz=convFactorToHz)

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

    #L1Ana.init_l1_analysis()
    print ""

    # holds the canvases, histograms, etc.
    objects = []

    ##########################################################################
    plotSetups = [{'q':12, 'thr':[18, 20, 22, 24], 'pthr':[26, 28, 30, 32]},
                 {'q':8, 'thr':[3, 5, 8, 10, 12], 'pthr':[5, 8, 12, 14, 16]},
                 {'q':4, 'thr':[3, 5, 8, 10, 12], 'pthr':[5, 8, 12, 14, 16]}]
    for plotSetup in plotSetups:
        qualMin = plotSetup['q']
        thresholdAndPointss = []
        for threshold, probeThreshold in zip(plotSetup['thr'], plotSetup['pthr']):
            rate_str = 'ugmt_highest_muon_absEtaMin0_absEtaMax2.5_qmin{qmin}'.format(qmin=qualMin)
            eff_l1_str = 'emu_best_l1_muon_qualMin{qmin}_ptmin{thr}_isoMax{iso:.3f}_dr0.5'
            eff_probe_str = 'probe_absEtaMin0_absEtaMax2.4_ptmin{pthr}.pass'.format(pthr=probeThreshold)

            if iso_type == 0:
                ref_iso = 31.
                iso_wps = [0, 1, 3, 5, 7, 9, 11, 15, 20, 25, 31]
            elif iso_type == 1:
                ref_iso = 62.
                iso_wps = [0., 1/2., 1/3., 2/3., 3/4., 4/5., 5/6., 6/7., 1., 1.5, 3., 5., 10., 31., 62.]
            elif iso_type == 3:
                ref_iso = 1.
                iso_wps = [0., 1/2., 1/3., 2/3., 3/4., 4/5., 5/6., 6/7., 7/8., 8/9., 9/10., 19/20., 29/30., 99/100.]
            elif iso_type == 5:
                ref_iso = 1.
                iso_wps = [0., 1/2., 1/3., 2/3., 3/4., 4/5., 5/6., 6/7., 7/8., 8/9., 9/10., 19/20., 30/31.]
            else:
                ref_iso = 0.
                iso_wps = [1., 2., 3., 4., 5., 6., 7., 8., 9.]
            iso_wps.sort()
            hNamesRef = {'num':eff_l1_str.format(qmin=qualMin, thr=threshold, iso=ref_iso)+'_matched_'+eff_probe_str,
                         'den':'emu_'+eff_probe_str}
            referenceRate = get_rate(hmRate, rate_str+'_isoMax{iso:.3f}_pt'.format(iso=ref_iso), threshold=threshold, scaleFactor=convFactorToHz / 1000.)
            referenceEff = get_eff(hmEff, hNamesRef)
            points = []
            print '==================='
            print 'L1 threshold: {thr} GeV, probe threshold: {pthr} GeV, quality >= {qmin}'.format(thr=threshold, pthr=probeThreshold, qmin=qualMin)
            print '-------------------'
            print ' iso   rel. rate   rel. eff'
            print '-------------------'
            for iso_wp in iso_wps:
                iso_wp_str = '_isoMax{iso:.3f}'.format(iso=iso_wp)
                #print_rates(hmRate, rate_str+iso_wp_str+'_pt', scaleFactor=convFactorToHz / 1000.)
                hNames = {'num':eff_l1_str.format(qmin=qualMin, thr=threshold, iso=iso_wp)+'_matched_'+eff_probe_str, 'den':'emu_'+eff_probe_str}
                #print_efficiencies(hmEff, hNames)
                relRate = get_relative_rate(hmRate, rate_str+iso_wp_str+'_pt', referenceRate[0], threshold=threshold, scaleFactor=convFactorToHz / 1000.)
                relEff = get_relative_eff(hmEff, hNames, referenceEff[0])
                print '{iso:.3f}    {rate:.4f}     {eff:.4f}'.format(iso=iso_wp, rate=relRate, eff=relEff)
                points.append((relRate, relEff))
            thresholdAndPointss.append((threshold, points))
            print '==================='

        objects.append(plot_roc_curves(thresholdAndPointss, qualMin))

    ##########################################################################
    # save plots to root file
    if savePlots:
        effFileString = opts.feff.replace('.root','')
        effFileString = effFileString.replace('ugmt_', '')
        effFileString = effFileString.replace('histos_', '')
        rateFileString = opts.frate.replace('.root','')
        rateFileString = rateFileString.replace('ugmt_iso_', '')
        rateFileString = rateFileString.replace('histos_', '')
        plotdir = 'plots_ROC_'+effFileString.partition('/')[0]+'_'+rateFileString.partition('/')[0]
        if opts.public:
            plotdir += '_public'
        if not os.path.exists(plotdir):
            os.makedirs(plotdir)
        output = root.TFile('./'+plotdir+'/ugmt_iso_roc_plots.root', 'recreate')
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
    plotLegacy = False
    font = 42
    fontSize = 0.04
    main()

