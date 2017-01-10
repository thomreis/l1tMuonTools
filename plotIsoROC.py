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

def plot_roc_curve(points):
    c = root.TCanvas('roc_curve', 'roc_curve', 100, 100, 600, 600)
    c.cd()
    set_root_style()

    roc_curve = root.TGraph(len(points))
    roc_curve.SetName('roc_graph')
    roc_curve.SetMarkerStyle(root.kFullCircle)
    for i, p in enumerate(points):
        roc_curve.SetPoint(i, p[0], p[1])

    roc_curve.Draw('ALP')

    xAxis = roc_curve.GetXaxis()
    yAxis = roc_curve.GetYaxis()
    xAxis.SetTitle('relative rate')
    xAxis.SetTitleFont(font)
    xAxis.SetLabelFont(font)
    xAxis.SetLabelSize(fontSize)
    yAxis.SetTitle('relative efficiency')
    yAxis.SetTitleOffset(1.5)
    yAxis.SetTitleFont(font)
    yAxis.SetLabelFont(font)
    yAxis.SetLabelSize(fontSize)

    line = root.TLine(0., 0., 1., 1.)
    line.SetLineStyle(root.kDashed)
    line.SetLineColor(root.kMagenta)
    line.Draw('same')

    c.Modified()
    c.Update()

    return c, roc_curve, line

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

    #L1Ana.init_l1_analysis()
    print ""

    # holds the canvases, histograms, etc.
    objects = []

    ##########################################################################
    ## uGMT kinematic variables
    #objects.append(plot_hists_standard(hm, 'muon_absEtaMin0_absEtaMax2.5_qmin12_varBin_pt', xTitle='p_{T} (GeV/c)', yTitle='# muons/(GeV/c)', stacked=True, normToBinWidth=True, data=thisIsData))
    #iso_wps = [0., 1/1., 1/2., 1/3., 2/3.]

    threshold = 12
    probeThreshold = 16
    rate_str = 'ugmt_highest_muon_absEtaMin0_absEtaMax2.5_qmin12'
    eff_l1_str = 'emu_best_l1_muon_qualMin12_ptmin{thr}_isoMax{iso:.3f}_dr0.5'
    eff_probe_str = 'probe_absEtaMin0_absEtaMax2.4_ptmin{pthr}.pass'.format(pthr=probeThreshold)
    hNamesRef = {'num':eff_l1_str.format(thr=threshold, iso=1.)+'_matched_'+eff_probe_str,
                 'den':'emu_'+eff_probe_str}

    referenceRate = get_rate(hmRate, rate_str+'_isoMax1.000_pt', threshold=threshold, scaleFactor=convFactorToHz / 1000.)
    referenceEff = get_eff(hmEff, hNamesRef)

    iso_wps = [0., 1/1., 1/2., 1/3., 2/3., 3/4., 4/5., 5/6., 6/7., 7/8., 8/9., 9/10., 19/20., 29/30., 99/100., 499/500., 999/1000.]
    iso_wps.sort()
    points = []
    for iso_wp in iso_wps:
        iso_wp_str = '_isoMax{iso:.3f}'.format(iso=iso_wp)
        #print_rates(hmRate, rate_str+iso_wp_str+'_pt', scaleFactor=convFactorToHz / 1000.)
        hNames = {'num':eff_l1_str.format(thr=threshold, iso=iso_wp)+'_matched_'+eff_probe_str, 'den':'emu_'+eff_probe_str}
        #print_efficiencies(hmEff, hNames)
        relRate = get_relative_rate(hmRate, rate_str+iso_wp_str+'_pt', referenceRate[0], threshold=threshold, scaleFactor=convFactorToHz / 1000.)
        relEff = get_relative_eff(hmEff, hNames, referenceEff[0])
        print 'iso {iso:.3f}: rate={rate}, eff={eff}'.format(iso=iso_wp, rate=relRate, eff=relEff)
        points.append((relRate, relEff))

    objects.append(plot_roc_curve(points))

    ##########################################################################
    # save plots to root file
    if savePlots:
        plotdir = 'plots_ROC_'+opts.feff.replace('.root','').partition('/')[0]
        if opts.public:
            plotdir += '_public'
        if not os.path.exists(plotdir):
            os.makedirs(plotdir)
        output = root.TFile('./'+plotdir+'/ugmt_rate_plots.root', 'recreate')
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

