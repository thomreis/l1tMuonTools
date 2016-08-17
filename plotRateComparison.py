#!/usr/bin/env python
from analysis_tools.plotting import HistManager
from analysis_tools.selections import MuonSelections, Matcher
import argparse
import ROOT as root
import re

def parse_options():
    """
    Adds often used options to the OptionParser...
    """
    parser = argparse.ArgumentParser(description="L1 Analysis Framework macro", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--interactive", dest="interactive", action='store_false', help="Draw plots on screen.")
    parser.add_argument("--fOne", dest="f1", default="", type=str, help="A root file containing reference histograms.")
    parser.add_argument("--fTwo", dest="f2", default="", type=str, help="A root file containing histograms to compare.")
    parser.add_argument("--nOne", dest="nevents1", default=1, type=int, help="Total nmumber of events for file 1")
    parser.add_argument("--bOne", dest="bunches1", default=0, type=int, help="Number of colliding bunches for file 1")
    parser.add_argument("--nTwo", dest="nevents2", default=1, type=int, help="Total nmumber of events for file 2")
    parser.add_argument("--bTwo", dest="bunches2", default=0, type=int, help="Number of colliding bunches for file 2")
    parser.add_argument("--pu", dest="pu", default=20, type=int, help="Average PU. default=20")
    parser.add_argument("--xsect", dest="xsect", default=80, type=float, help="Total cross section in mb. default=80 mb")
    parser.add_argument("--instlumi", dest="instlumi", default=1.2e34, type=float, help="Instantaneous luminosity. default=1.2e-34 cm-2s-1")
    parser.add_argument("--scale", dest="scale", default=1., type=float, help="Additional scale factor for rate calculate")
    parser.add_argument("-l", "--legacy", dest="legacy", action='store_true', help="Draw plots relative to legacy.")

    opts, unknown = parser.parse_known_args()
    return opts

def set_root_style():
    root.gStyle.SetTitleFont(font)
    root.gStyle.SetStatFont(font)
    root.gStyle.SetTextFont(font)
    root.gStyle.SetLabelFont(font)
    root.gStyle.SetLegendFont(font)
    root.gStyle.SetMarkerStyle(20)
    root.gStyle.SetOptStat(0)
    root.gStyle.SetOptFit(0)
    root.gStyle.SetOptTitle(0)
    root.gPad.SetTopMargin(0.08)
    root.gPad.SetLeftMargin(0.14)
    root.gPad.SetRightMargin(0.06)
    root.gPad.SetTickx(1)
    root.gPad.SetTicky(1)

def plot_hists_comp(hDefs, hName=None, xTitle=None, yTitle='# muons', threshold=False, normToBinWidth=False, canvasPrefix='', notes=None, data=False):
    scaled = False
    ratioPadRelSize = 0.3
    histoPadFontSizeFactor = 1 / (1. - ratioPadRelSize)
    ratioPadFontSizeFactor = 1 / ratioPadRelSize
    textFontSize = 0.04
    axisLabelSize = 0.035
    axisTitleSize = 0.04
    leftMargin = 0.12
    rightMargin = 0.06

    if hName != None:
        for hDef in hDefs:
            hDef['name'] = hDef['namePrefix'] + hName

    name = canvasPrefix+hDefs[0]['name']
    if normToBinWidth and not threshold:
        name = 'normToBinWidth_'+name

    # setup legend according to how many histograms are in the plot
    legYmin = 0.9-0.08*len(hDefs)
    legXmin = 0.48
    legXmax = 0.8
    canvWidth = 600
    if legYmin < 0.6:
        legXmin = 0.48
        #legXmax = 1.
        #canvWidth = 730
    legend = root.TLegend(legXmin, legYmin, legXmax, 0.9)
    legend.SetTextFont(font)
    legend.SetTextSize(textFontSize * histoPadFontSizeFactor)
    legend.SetBorderSize(0)
    legend.SetFillColor(19)
    legend.SetFillStyle(0)
    #legend.SetNColumns(2)
    legEntries = []

    hs = []
    hrs = []
    hStack = root.THStack()
    maxBinValue = 0
    # get all the histograms and set their plot style
    for i, hDef in enumerate(hDefs):
        hm = hDef['hm']
        if threshold:
            h = hm.get_threshold_hist(hDef['name']).Clone()
        else:
            h = hm.get(hDef['name']).Clone()

        if normToBinWidth and not threshold:
            for bin in range(1, h.GetNbinsX()+1):
               h.SetBinContent(bin, h.GetBinContent(bin) / h.GetBinWidth(bin))
               h.SetBinError(bin, h.GetBinError(bin) / h.GetBinWidth(bin))
        elif normToBinWidth:
            print 'Ignoring normToBinWidth flag for threshold'

        if hDef['sf'] != 1.:
            h.Scale(hDef['sf'])
            scaled = True
        h.SetLineColor(hDef['lc'])
        h.SetLineStyle(hDef['ls'])
        h.SetLineWidth(2)
        h.GetXaxis().SetLabelSize(axisLabelSize * histoPadFontSizeFactor)
        h.GetYaxis().SetLabelSize(axisLabelSize * histoPadFontSizeFactor)
        h.GetXaxis().SetTitleSize(axisTitleSize * histoPadFontSizeFactor)
        h.GetYaxis().SetTitleSize(axisTitleSize * histoPadFontSizeFactor)
        h.GetXaxis().SetTitleOffset(0.8 * histoPadFontSizeFactor)
        h.GetYaxis().SetTitleOffset(1.)
        legStyle = 'l'
        if hDef['fc']:
            h.SetFillColor(hDef['fc'])
            h.SetLineWidth(1)
            legStyle = 'f'
            hStack.Add(h)
        legEntries.append(legend.AddEntry(h, hDef['legtext'], legStyle))
        if h.GetBinContent(h.GetMaximumBin()) > maxBinValue:
            maxBinValue = h.GetBinContent(h.GetMaximumBin())
        hs.append(h)

        # ratio histograms
        if i == 0:
            hDen = h.Clone()
        hRatio = h.Clone('hRatio_h{nom}_over_h0'.format(nom=i))
        hRatio.Divide(hRatio, hDen, 1, 1, "b")
        hRatio.GetXaxis().SetLabelSize(axisLabelSize * ratioPadFontSizeFactor)
        hRatio.GetYaxis().SetLabelSize(axisLabelSize * ratioPadFontSizeFactor)
        hRatio.GetXaxis().SetTitleSize(axisTitleSize * ratioPadFontSizeFactor)
        hRatio.GetYaxis().SetTitleSize(axisTitleSize * ratioPadFontSizeFactor)
        hRatio.GetXaxis().SetTitleOffset(1.)
        hRatio.GetYaxis().SetTitleOffset(1.)
        hRatio.GetYaxis().SetNdivisions(205)
        hrs.append(hRatio)

    # replace histograms to be stacked with stack histograms
    if hStack.GetNhists() > 0:
        styles = hist_styles(True)
        canvas_name = 'c_rates_stacked_'+name
        stackHistos = hStack.GetStack()
        j = len(stackHistos)-1
        hStackSum = stackHistos[j].Clone()
        for i, hDef in enumerate(hDefs):
            if hDef['fc']:
                hs[i] = stackHistos[j].Clone()
                j -= 1
        hStackRatio = hStackSum.Clone('hStackRatio_hStack_over_h0')
        hStackRatio.Divide(hStackRatio, hDen, 1, 1, "b")
        hStackRatio.SetLineColor(styles['emul_ratio']['lc'])
        hStackRatio.SetLineStyle(styles['emul_ratio']['ls'])
        hStackRatio.SetLineWidth(2)
        hStackRatio.GetXaxis().SetLabelSize(axisLabelSize * ratioPadFontSizeFactor)
        hStackRatio.GetYaxis().SetLabelSize(axisLabelSize * ratioPadFontSizeFactor)
        hStackRatio.GetXaxis().SetTitleSize(axisTitleSize * ratioPadFontSizeFactor)
        hStackRatio.GetYaxis().SetTitleSize(axisTitleSize * ratioPadFontSizeFactor)
        hStackRatio.GetXaxis().SetTitleOffset(1.)
        hStackRatio.GetYaxis().SetTitleOffset(1.)
        hStackRatio.GetYaxis().SetNdivisions(205)
    else:
        canvas_name = 'c_rates_'+name

    if scaled:
        canvas_name += '_scaled'

    # create canvas and draw on it
    canvas_title = canvas_name
    c = root.TCanvas(canvas_name, canvas_title, 100, 100, canvWidth, 600)
    p1 = root.TPad(canvas_name + '_p1', canvas_title + '_p1', 0., ratioPadRelSize, 1., 1.)
    p2 = root.TPad(canvas_name + '_p2', canvas_title + '_p2', 0., 0., 1., ratioPadRelSize)
    p1.Draw()
    p2.Draw()

    # Pad 1: histograms
    p1.cd()
    p1.SetBottomMargin(0.05)
    p1.SetLeftMargin(leftMargin)
    p1.SetRightMargin(rightMargin)
    p1.SetTicks()

    if name[-2:] == 'pt':
        p1.SetLogy(True)

    set_root_style()

    #if legYmin < 0.6:
    #    root.gPad.SetRightMargin(0.2)

    if xTitle:
        hs[0].GetXaxis().SetTitle(xTitle)
    hs[0].GetYaxis().SetTitleOffset(1.6 / histoPadFontSizeFactor)
    hs[0].GetYaxis().SetTitle(yTitle)
    if not p1.GetLogy():
        yMax = 1.1*maxBinValue
        #if maxBinValue <= 1.:
        #    yMax = 1.3
        hs[0].GetYaxis().SetRangeUser(0., yMax)
    # draw
    hs[0].SetLineWidth(2)
    legEntries[0].SetObject(hs[0])
    legEntries[0].SetOption(legEntries[0].GetOption()+'le')
    hs[0].Draw('pe0')
    for h in hs[1:]:
        h.Draw('histsame')
    hs[0].Draw('pe0same')
    hs[0].Draw('sameaxis')

    # draw vertical lines to mark TF boundaries
    lines = []
    if name[-3:] == 'eta':
        lines.append(root.TLine(-0.83, 0., -0.83, yMax))
        lines[-1].SetLineStyle(root.kDotted)
        lines[-1].Draw('same')
        lines.append(root.TLine(-1.24, 0., -1.24, yMax))
        lines[-1].SetLineStyle(root.kDotted)
        lines[-1].Draw('same')
        lines.append(root.TLine(0.83, 0., 0.83, yMax))
        lines[-1].SetLineStyle(root.kDotted)
        lines[-1].Draw('same')
        lines.append(root.TLine(1.24, 0., 1.24, yMax))
        lines[-1].SetLineStyle(root.kDotted)
        lines[-1].Draw('same')

    legend.Draw('same')

    tex = root.TLatex()
    tex.SetNDC()
    tex.SetTextFont(font)
    tex.SetTextSize(0.04 * histoPadFontSizeFactor)
    #tex.DrawLatex(0.484, 0.93, 'Simulation, 13 TeV')
    if canvWidth > 600:
        if data:
            #tex.DrawLatex(0.48, 0.93, 'CMS preliminary, 13 TeV')
            tex.DrawLatex(0.48, 0.93, 'CMS internal, 13 TeV')
        else:
            tex.DrawLatex(0.484, 0.93, 'CMS Simulation, 13 TeV')
    else:
        if data:
            #tex.DrawLatex(0.551, 0.93, 'CMS preliminary, 13 TeV')
            tex.DrawLatex(0.551, 0.93, 'CMS internal, 13 TeV')
        else:
            tex.DrawLatex(0.555, 0.93, 'CMS Simulation, 13 TeV')
    if notes:
        tex.SetTextSize(0.035)
        for note in notes:
            tex.DrawLatex(note[0], note[1], note[2])

    # Pad 2: ratio
    p2.cd()
    p2.SetTopMargin(0.1)
    p2.SetBottomMargin(0.3)
    p2.SetLeftMargin(leftMargin)
    p2.SetRightMargin(rightMargin)
    p2.SetTicks()
    p2.SetGridy(True)
    set_root_style()

    yMin = 0.5
    yMax = 1.5
    if hStack.GetNhists() > 0:
        if xTitle:
            hStackRatio.GetXaxis().SetTitle(xTitle)
        hStackRatio.GetYaxis().SetTitleOffset(1.6 / ratioPadFontSizeFactor)
        hStackRatio.GetYaxis().SetTitle('ratio')
        hStackRatio.GetYaxis().SetRangeUser(yMin, yMax)
        # draw
        hStackRatio.SetLineWidth(2)
        hStackRatio.Draw()
        hStackRatio.Draw('sameaxis')
    else:
        # first histo ratio to plot (0 would be hDen/hDen = 1)
        firstHist = 1

        if xTitle:
            hrs[firstHist].GetXaxis().SetTitle(xTitle)
        hrs[firstHist].GetYaxis().SetTitleOffset(1.6 / ratioPadFontSizeFactor)
        hrs[firstHist].GetYaxis().SetTitle('ratio')
        hrs[firstHist].GetYaxis().SetRangeUser(yMin, yMax)
        # draw
        hrs[firstHist].SetLineWidth(2)
        hrs[firstHist].Draw()
        for h in hrs[2:]:
            h.Draw('same')
        hrs[firstHist].Draw('same')
        hrs[firstHist].Draw('sameaxis')

    c.Modified()
    c.Update()

    return [c, p1, p2, hs, hrs, hStackRatio, legend, lines, tex]

def print_rates(hm, hName, scaleFactor=1.):
    hNames = ['gmt_'+hName.replace('qmin12', 'qmin8'), 'ugmt_'+hName, 'bmtf_ugmt_'+hName, 'omtf_ugmt_'+hName, 'emtf_ugmt_'+hName]

    print '===== Rates ====='
    print hName
    print ''
    histos = []
    print 'System        16 GeV        20 GeV        25 GeV'
    for name in hNames:
        histos.append(hm.get_threshold_hist(name).Clone())
        if scaleFactor != 1.:
            histos[-1].Scale(scaleFactor)

        bin16 = histos[-1].FindBin(16)
        bin20 = histos[-1].FindBin(20)
        bin25 = histos[-1].FindBin(25)
        print '{name} rate: {sixteengev:>7.2f} kHz   {twentygev:>7.2f} kHz   {twentyfivegev:>7.2f} kHz'.format(name=name.split('_')[0], sixteengev=histos[-1].GetBinContent(bin16), twentygev=histos[-1].GetBinContent(bin20), twentyfivegev=histos[-1].GetBinContent(bin25))

    print '\nThreshold           GMT                   uGMT             ratio'
    for threshold in [0, 3, 5, 7, 10, 12, 14, 16, 18, 20, 22, 25, 30, 40, 50, 60]:
        gmtBinNr = histos[0].FindBin(threshold)
        ugmtBinNr = histos[1].FindBin(threshold)
        gmtCont = histos[0].GetBinContent(gmtBinNr)
        ugmtCont = histos[1].GetBinContent(ugmtBinNr)
        gmtErr = histos[0].GetBinError(gmtBinNr)
        ugmtErr = histos[1].GetBinError(ugmtBinNr)
        ratio = -1.
        if gmtCont != 0:
            ratio = ugmtCont/gmtCont
        print '{threshold:>3} GeV:   {gmt:>8.3f} +/- {gmterr:>5.3f} kHz   {ugmt:>8.3f} +/- {ugmterr:>5.3f} kHz   {ratio:>8.3f}'.format(threshold=threshold, gmt=gmtCont, gmterr=gmtErr, ugmt=ugmtCont, ugmterr=ugmtErr, ratio=ratio)

    print '================='

def hist_styles(stacked=False):
    styles = {}
    styles['gmt'] = {'lc':root.kCyan, 'ls':root.kSolid, 'fc':None, 'legtext':'GMT'}
    styles['ugmt'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':None, 'legtext':'uGMT'}
    styles['unpacked'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':None, 'legtext':'unpacked'}
    styles['unpacked_bmtf'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':None, 'legtext':'unpacked BMTF'}
    styles['unpacked_omtf'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':None, 'legtext':'unpacked OMTF'}
    styles['unpacked_emtf'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':None, 'legtext':'unpacked EMTF'}
    styles['emul'] = {'lc':root.kMagenta-4, 'ls':root.kSolid, 'fc':root.kMagenta-7, 'legtext':'emulated'}
    styles['emul_ratio'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':None, 'legtext':'emulated / unpacked'}
    styles['emul_bmtf'] = {'lc':root.kBlue-4, 'ls':root.kSolid, 'fc':root.kBlue-7, 'legtext':'emulated BMTF'}
    styles['emul_omtf'] = {'lc':root.kGreen-4, 'ls':root.kSolid, 'fc':root.kGreen-7, 'legtext':'emulated OMTF'}
    styles['emul_emtf'] = {'lc':root.kRed-4, 'ls':root.kSolid, 'fc':root.kRed-7, 'legtext':'emulated EMTF'}
    if stacked:
        styles['bmtf_ugmt'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':root.kBlue-4, 'legtext':'BMTF uGMT'}
        styles['omtf_ugmt'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':root.kGreen-4, 'legtext':'OMTF uGMT'}
        styles['emtf_ugmt'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':root.kRed-4, 'legtext':'EMTF uGMT'}
        styles['bmtf'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':root.kBlue-4, 'legtext':'BMTF'}
        styles['omtf'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':root.kGreen-4, 'legtext':'OMTF'}
        styles['emtf'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':root.kRed-4, 'legtext':'EMTF'}
    else:
        styles['bmtf_ugmt'] = {'lc':root.kBlue-4, 'ls':root.kSolid, 'fc':None, 'legtext':'BMTF uGMT'}
        styles['omtf_ugmt'] = {'lc':root.kGreen-4, 'ls':root.kSolid, 'fc':None, 'legtext':'OMTF uGMT'}
        styles['emtf_ugmt'] = {'lc':root.kRed-4, 'ls':root.kSolid, 'fc':None, 'legtext':'EMTF uGMT'}
        styles['bmtf_ugmt_q'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':None, 'legtext':'BMTF uGMT'}
        styles['omtf_ugmt_q'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':None, 'legtext':'OMTF uGMT'}
        styles['emtf_ugmt_q'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':None, 'legtext':'EMTF uGMT'}
        styles['bmtf'] = {'lc':root.kBlue-4, 'ls':root.kSolid, 'fc':None, 'legtext':'BMTF'}
        styles['omtf'] = {'lc':root.kGreen-4, 'ls':root.kSolid, 'fc':None, 'legtext':'OMTF'}
        styles['emtf'] = {'lc':root.kRed-4, 'ls':root.kSolid, 'fc':None, 'legtext':'EMTF'}
        styles['bmtf_q'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':None, 'legtext':'BMTF'}
        styles['omtf_q'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':None, 'legtext':'OMTF'}
        styles['emtf_q'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':None, 'legtext':'EMTF'}

    styles['ugmt_q0'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':root.kRed+3, 'legtext':'uGMT q0'}
    styles['ugmt_q1'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':root.kRed, 'legtext':'uGMT q1'}
    styles['ugmt_q2'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':root.kOrange+8, 'legtext':'uGMT q2'}
    styles['ugmt_q3'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':root.kOrange, 'legtext':'uGMT q3'}
    styles['ugmt_q4'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':root.kYellow, 'legtext':'uGMT q4'}
    styles['ugmt_q5'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':root.kGreen, 'legtext':'uGMT q5'}
    styles['ugmt_q6'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':root.kGreen+3, 'legtext':'uGMT q6'}
    styles['ugmt_q7'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':root.kCyan, 'legtext':'uGMT q7'}
    styles['ugmt_q8'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':root.kCyan+3, 'legtext':'uGMT q8'}
    styles['ugmt_q9'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':root.kAzure+7, 'legtext':'uGMT q9'}
    styles['ugmt_q10'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':root.kBlue, 'legtext':'uGMT q10'}
    styles['ugmt_q11'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':root.kBlue+3, 'legtext':'uGMT q11'}
    styles['ugmt_q12'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':root.kViolet, 'legtext':'uGMT q12'}
    styles['ugmt_q13'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':root.kMagenta, 'legtext':'uGMT q13'}
    styles['ugmt_q14'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':root.kMagenta+3, 'legtext':'uGMT q14'}
    styles['ugmt_q15'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':root.kViolet+3, 'legtext':'uGMT q15'}

    styles['tf_q0'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':root.kRed+3, 'legtext':'TF q0'}
    styles['tf_q1'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':root.kRed, 'legtext':'TF q1'}
    styles['tf_q2'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':root.kOrange+8, 'legtext':'TF q2'}
    styles['tf_q3'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':root.kOrange, 'legtext':'TF q3'}
    styles['tf_q4'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':root.kYellow, 'legtext':'TF q4'}
    styles['tf_q5'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':root.kGreen, 'legtext':'TF q5'}
    styles['tf_q6'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':root.kGreen+3, 'legtext':'TF q6'}
    styles['tf_q7'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':root.kCyan, 'legtext':'TF q7'}
    styles['tf_q8'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':root.kCyan+3, 'legtext':'TF q8'}
    styles['tf_q9'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':root.kAzure+7, 'legtext':'TF q9'}
    styles['tf_q10'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':root.kBlue, 'legtext':'TF q10'}
    styles['tf_q11'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':root.kBlue+3, 'legtext':'TF q11'}
    styles['tf_q12'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':root.kViolet, 'legtext':'TF q12'}
    styles['tf_q13'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':root.kMagenta, 'legtext':'TF q13'}
    styles['tf_q14'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':root.kMagenta+3, 'legtext':'TF q14'}
    styles['tf_q15'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':root.kViolet+3, 'legtext':'TF q15'}

    return styles

def plot_hists_standard(hm, hName, den=None, xTitle='', yTitle='# muons', threshold=False, stacked=False, normToBinWidth=False, tfMuonOrig='ugmt', reg='', scaleFactor=1., data=False):
    styles = hist_styles(stacked)

    if tfMuonOrig == 'ugmt':
        ugmt_str = '_ugmt'
        prefix = ''
    elif tfMuonOrig == 'tf':
        ugmt_str = ''
        prefix = 'tf_'

    ugmt_dict = {'num':'ugmt_'+hName, 'den':den}
    bmtf_dict = {'num':'bmtf'+ugmt_str+'_'+hName, 'den':den}
    omtf_dict = {'num':'omtf'+ugmt_str+'_'+hName, 'den':den}
    emtf_dict = {'num':'emtf'+ugmt_str+'_'+hName, 'den':den}
    ugmt_dict.update(styles['ugmt'])
    bmtf_dict.update(styles['bmtf'+ugmt_str])
    omtf_dict.update(styles['omtf'+ugmt_str])
    emtf_dict.update(styles['emtf'+ugmt_str])
    hDefs = []
    if reg == '':
        hDefs.append(ugmt_dict)
        hDefs.append(bmtf_dict)
        hDefs.append(omtf_dict)
        hDefs.append(emtf_dict)
    elif reg == 'b':
        hDefs.append(bmtf_dict)
        prefix += 'bmtf_'
    elif reg == 'o':
        hDefs.append(omtf_dict)
        prefix += 'omtf_'
    elif reg == 'e':
        hDefs.append(emtf_dict)
        prefix += 'emtf_'
    if plotLegacy:
        if den:
            gmt_dict = {'num':den, 'den':den}
            gmt_dict.update(styles['gmt'])
        else:
            gmt_dict = {'num':'gmt_'+hName.replace('qmin12', 'qmin8'), 'den':den}
            gmt_dict.update(styles['gmt'])
        hDefs.append(gmt_dict)

    # extract eta range from histogram name
    eta_number_strs = re.findall(r'[\d\.\d]+', hName[hName.find('EtaMin')+6:hName.find('EtaMax')+12])
    if len(eta_number_strs) > 1:
        note_str = eta_number_strs[0]+' < |#eta| < '+eta_number_strs[1]
        notes = [[0.17, 0.86, note_str]]
        if den:
            den_eta_number_strs = re.findall(r'[\d\.\d]+', den[den.find('EtaMin')+6:den.find('EtaMax')+12])
            if len(den_eta_number_strs) > 1 and eta_number_strs != den_eta_number_strs:
                den_note_str = den_eta_number_strs[0]+' < |#eta^{GMT}| < '+den_eta_number_strs[1]
                notes.append([0.17, 0.81, den_note_str])
    else:
        notes = None

    return plot_hists(hm, hDefs, xTitle, yTitle, threshold, normToBinWidth, prefix, notes, scaleFactor, data)

def main():
    plotLegacy = opts.legacy
    batchRun = opts.interactive
    if batchRun:
        root.gROOT.SetBatch(True)

    # calculate the scale factor for rate in Hz
    orbitFreq = 11245.6
    nCollBunches1 = opts.bunches1
    nZeroBiasEvents1 = opts.nevents1
    nCollBunches2 = opts.bunches2
    nZeroBiasEvents2 = opts.nevents2
    crossSect = opts.xsect
    instLumi = opts.instlumi
    pu = opts.pu
    thisIsData=True
    # determine that this is MC if there is no number of colliding bunches given (defaults to 0 then)
    if nCollBunches1 == 0:
        print "No number of colliding bunches given. Assuming this is MC"
        print "Using {instLumi} cm-2s-1 as instantaneous luminosity, {crossSect} mb as cross section, and {pu} as average number of pileup to determine number of colliding bunches.".format(instLumi=instLumi, crossSect=crossSect, pu=pu)
        nCollBunches1 = round(instLumi * crossSect*1e-27 / (pu * orbitFreq))
        nCollBunches2 = nCollBunches1
        thisIsData=False
    else:
        print "Assuming this is data"
    convFactorToHz1 = orbitFreq * nCollBunches1 / nZeroBiasEvents1
    convFactorToHz2 = orbitFreq * nCollBunches2 / nZeroBiasEvents2
    print 'Conversion factor to rate in Hz for file 1 with {orbitFreq} Hz orbit frequency, {nCollBunches1} colliding bunches and {nZeroBiasEvents1} analyzed zero bias events: {convFactorToHz1}'.format(orbitFreq=orbitFreq, nCollBunches1=nCollBunches1, nZeroBiasEvents1=nZeroBiasEvents1, convFactorToHz1=convFactorToHz1)
    print 'Conversion factor to rate in Hz for file 2 with {orbitFreq} Hz orbit frequency, {nCollBunches2} colliding bunches and {nZeroBiasEvents2} analyzed zero bias events: {convFactorToHz2}'.format(orbitFreq=orbitFreq, nCollBunches2=nCollBunches2, nZeroBiasEvents2=nZeroBiasEvents2, convFactorToHz2=convFactorToHz2)
    if opts.scale != 1.:
        convFactorToHz1 *= opts.scale
        print 'Conversion factor for file 1 after applying additinoal scale factor of {sf}: {convFactorToHz1}'.format(sf=opts.scale, convFactorToHz1=convFactorToHz1)
        print 'Conversion factor for file 2 after applying additinoal scale factor of {sf}: {convFactorToHz2}'.format(sf=opts.scale, convFactorToHz2=convFactorToHz2)

    hm1 = HistManager(filename=opts.f1)
    hm2 = HistManager(filename=opts.f2)

    convFactorTokHz1 = convFactorToHz1 / 1000.
    convFactorTokHz2 = convFactorToHz2 / 1000.

    styles = hist_styles(False)
    hName = 'dummyName'
    dict1      = {'hm':hm1, 'sf':convFactorTokHz1, 'name':hName, 'namePrefix':''}
    bmtf_dict1 = {'hm':hm1, 'sf':convFactorTokHz1, 'name':hName, 'namePrefix':'bmtf_'}
    omtf_dict1 = {'hm':hm1, 'sf':convFactorTokHz1, 'name':hName, 'namePrefix':'omtf_'}
    emtf_dict1 = {'hm':hm1, 'sf':convFactorTokHz1, 'name':hName, 'namePrefix':'emtf_'}
    dict2      = {'hm':hm2, 'sf':convFactorTokHz2, 'name':hName, 'namePrefix':''}
    bmtf_dict2 = {'hm':hm2, 'sf':convFactorTokHz2, 'name':hName, 'namePrefix':'bmtf_'}
    omtf_dict2 = {'hm':hm2, 'sf':convFactorTokHz2, 'name':hName, 'namePrefix':'omtf_'}
    emtf_dict2 = {'hm':hm2, 'sf':convFactorTokHz2, 'name':hName, 'namePrefix':'emtf_'}
    dict1.update(styles['unpacked'])
    bmtf_dict1.update(styles['unpacked_bmtf'])
    omtf_dict1.update(styles['unpacked_omtf'])
    emtf_dict1.update(styles['unpacked_emtf'])
    dict2.update(styles['emul'])
    bmtf_dict2.update(styles['emul_bmtf'])
    omtf_dict2.update(styles['emul_omtf'])
    emtf_dict2.update(styles['emul_emtf'])

    hDefs = [dict1, dict2]
    hDefs_bmtf = [bmtf_dict1, bmtf_dict2]
    hDefs_omtf = [omtf_dict1, omtf_dict2]
    hDefs_emtf = [emtf_dict1, emtf_dict2]
    hDefs_withTfs = [dict1, bmtf_dict2, omtf_dict2, emtf_dict2]

    # holds the canvases, histograms, etc.
    objects = []

    ##########################################################################
    # uGMT kinematic variables

    # uGMT rates for regions
    objects.append(plot_hists_comp(hDefs, hName='ugmt_highest_muon_absEtaMin0_absEtaMax2.5_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='kHz', threshold=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_bmtf, hName='ugmt_highest_muon_absEtaMin0_absEtaMax2.5_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='kHz', threshold=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_omtf, hName='ugmt_highest_muon_absEtaMin0_absEtaMax2.5_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='kHz', threshold=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_emtf, hName='ugmt_highest_muon_absEtaMin0_absEtaMax2.5_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='kHz', threshold=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs, hName='ugmt_highest_muon_absEtaMin0_absEtaMax2.5_qmin8_pt', xTitle='p_{T} (GeV/c)', yTitle='kHz', threshold=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_bmtf, hName='ugmt_highest_muon_absEtaMin0_absEtaMax2.5_qmin8_pt', xTitle='p_{T} (GeV/c)', yTitle='kHz', threshold=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_omtf, hName='ugmt_highest_muon_absEtaMin0_absEtaMax2.5_qmin8_pt', xTitle='p_{T} (GeV/c)', yTitle='kHz', threshold=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_emtf, hName='ugmt_highest_muon_absEtaMin0_absEtaMax2.5_qmin8_pt', xTitle='p_{T} (GeV/c)', yTitle='kHz', threshold=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs, hName='ugmt_highest_muon_absEtaMin0_absEtaMax2.5_qmin4_pt', xTitle='p_{T} (GeV/c)', yTitle='kHz', threshold=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_bmtf, hName='ugmt_highest_muon_absEtaMin0_absEtaMax2.5_qmin4_pt', xTitle='p_{T} (GeV/c)', yTitle='kHz', threshold=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_omtf, hName='ugmt_highest_muon_absEtaMin0_absEtaMax2.5_qmin4_pt', xTitle='p_{T} (GeV/c)', yTitle='kHz', threshold=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_emtf, hName='ugmt_highest_muon_absEtaMin0_absEtaMax2.5_qmin4_pt', xTitle='p_{T} (GeV/c)', yTitle='kHz', threshold=True, data=thisIsData))

    objects.append(plot_hists_comp(hDefs_withTfs, hName='ugmt_muon_ptmin0_qmin12_eta', xTitle='#eta', yTitle='', normToBinWidth=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_bmtf, hName='ugmt_muon_ptmin0_qmin12_eta', xTitle='#eta', yTitle='', normToBinWidth=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_omtf, hName='ugmt_muon_ptmin0_qmin12_eta', xTitle='#eta', yTitle='', normToBinWidth=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_emtf, hName='ugmt_muon_ptmin0_qmin12_eta', xTitle='#eta', yTitle='', normToBinWidth=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_withTfs, hName='ugmt_muon_ptmin4_qmin12_eta', xTitle='#eta', yTitle='', normToBinWidth=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_bmtf, hName='ugmt_muon_ptmin4_qmin12_eta', xTitle='#eta', yTitle='', normToBinWidth=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_omtf, hName='ugmt_muon_ptmin4_qmin12_eta', xTitle='#eta', yTitle='', normToBinWidth=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_emtf, hName='ugmt_muon_ptmin4_qmin12_eta', xTitle='#eta', yTitle='', normToBinWidth=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_withTfs, hName='ugmt_muon_ptmin8_qmin12_eta', xTitle='#eta', yTitle='', normToBinWidth=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_bmtf, hName='ugmt_muon_ptmin8_qmin12_eta', xTitle='#eta', yTitle='', normToBinWidth=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_omtf, hName='ugmt_muon_ptmin8_qmin12_eta', xTitle='#eta', yTitle='', normToBinWidth=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_emtf, hName='ugmt_muon_ptmin8_qmin12_eta', xTitle='#eta', yTitle='', normToBinWidth=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_withTfs, hName='ugmt_muon_ptmin11_qmin12_eta', xTitle='#eta', yTitle='', normToBinWidth=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_bmtf, hName='ugmt_muon_ptmin11_qmin12_eta', xTitle='#eta', yTitle='', normToBinWidth=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_omtf, hName='ugmt_muon_ptmin11_qmin12_eta', xTitle='#eta', yTitle='', normToBinWidth=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_emtf, hName='ugmt_muon_ptmin11_qmin12_eta', xTitle='#eta', yTitle='', normToBinWidth=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_withTfs, hName='ugmt_muon_ptmin18_qmin12_eta', xTitle='#eta', yTitle='', normToBinWidth=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_bmtf, hName='ugmt_muon_ptmin18_qmin12_eta', xTitle='#eta', yTitle='', normToBinWidth=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_omtf, hName='ugmt_muon_ptmin18_qmin12_eta', xTitle='#eta', yTitle='', normToBinWidth=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_emtf, hName='ugmt_muon_ptmin18_qmin12_eta', xTitle='#eta', yTitle='', normToBinWidth=True, data=thisIsData))

    objects.append(plot_hists_comp(hDefs_withTfs, hName='ugmt_muon_ptmin0_qmin8_eta', xTitle='#eta', yTitle='', normToBinWidth=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_bmtf, hName='ugmt_muon_ptmin0_qmin8_eta', xTitle='#eta', yTitle='', normToBinWidth=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_omtf, hName='ugmt_muon_ptmin0_qmin8_eta', xTitle='#eta', yTitle='', normToBinWidth=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_emtf, hName='ugmt_muon_ptmin0_qmin8_eta', xTitle='#eta', yTitle='', normToBinWidth=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_withTfs, hName='ugmt_muon_ptmin4_qmin8_eta', xTitle='#eta', yTitle='', normToBinWidth=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_bmtf, hName='ugmt_muon_ptmin4_qmin8_eta', xTitle='#eta', yTitle='', normToBinWidth=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_omtf, hName='ugmt_muon_ptmin4_qmin8_eta', xTitle='#eta', yTitle='', normToBinWidth=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_emtf, hName='ugmt_muon_ptmin4_qmin8_eta', xTitle='#eta', yTitle='', normToBinWidth=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_withTfs, hName='ugmt_muon_ptmin8_qmin8_eta', xTitle='#eta', yTitle='', normToBinWidth=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_bmtf, hName='ugmt_muon_ptmin8_qmin8_eta', xTitle='#eta', yTitle='', normToBinWidth=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_omtf, hName='ugmt_muon_ptmin8_qmin8_eta', xTitle='#eta', yTitle='', normToBinWidth=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_emtf, hName='ugmt_muon_ptmin8_qmin8_eta', xTitle='#eta', yTitle='', normToBinWidth=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_withTfs, hName='ugmt_muon_ptmin11_qmin8_eta', xTitle='#eta', yTitle='', normToBinWidth=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_bmtf, hName='ugmt_muon_ptmin11_qmin8_eta', xTitle='#eta', yTitle='', normToBinWidth=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_omtf, hName='ugmt_muon_ptmin11_qmin8_eta', xTitle='#eta', yTitle='', normToBinWidth=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_emtf, hName='ugmt_muon_ptmin11_qmin8_eta', xTitle='#eta', yTitle='', normToBinWidth=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_withTfs, hName='ugmt_muon_ptmin18_qmin8_eta', xTitle='#eta', yTitle='', normToBinWidth=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_bmtf, hName='ugmt_muon_ptmin18_qmin8_eta', xTitle='#eta', yTitle='', normToBinWidth=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_omtf, hName='ugmt_muon_ptmin18_qmin8_eta', xTitle='#eta', yTitle='', normToBinWidth=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_emtf, hName='ugmt_muon_ptmin18_qmin8_eta', xTitle='#eta', yTitle='', normToBinWidth=True, data=thisIsData))

    objects.append(plot_hists_comp(hDefs_withTfs, hName='ugmt_muon_ptmin0_qmin4_eta', xTitle='#eta', yTitle='', normToBinWidth=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_bmtf, hName='ugmt_muon_ptmin0_qmin4_eta', xTitle='#eta', yTitle='', normToBinWidth=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_omtf, hName='ugmt_muon_ptmin0_qmin4_eta', xTitle='#eta', yTitle='', normToBinWidth=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_emtf, hName='ugmt_muon_ptmin0_qmin4_eta', xTitle='#eta', yTitle='', normToBinWidth=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_withTfs, hName='ugmt_muon_ptmin4_qmin4_eta', xTitle='#eta', yTitle='', normToBinWidth=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_bmtf, hName='ugmt_muon_ptmin4_qmin4_eta', xTitle='#eta', yTitle='', normToBinWidth=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_omtf, hName='ugmt_muon_ptmin4_qmin4_eta', xTitle='#eta', yTitle='', normToBinWidth=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_emtf, hName='ugmt_muon_ptmin4_qmin4_eta', xTitle='#eta', yTitle='', normToBinWidth=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_withTfs, hName='ugmt_muon_ptmin8_qmin4_eta', xTitle='#eta', yTitle='', normToBinWidth=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_bmtf, hName='ugmt_muon_ptmin8_qmin4_eta', xTitle='#eta', yTitle='', normToBinWidth=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_omtf, hName='ugmt_muon_ptmin8_qmin4_eta', xTitle='#eta', yTitle='', normToBinWidth=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_emtf, hName='ugmt_muon_ptmin8_qmin4_eta', xTitle='#eta', yTitle='', normToBinWidth=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_withTfs, hName='ugmt_muon_ptmin11_qmin4_eta', xTitle='#eta', yTitle='', normToBinWidth=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_bmtf, hName='ugmt_muon_ptmin11_qmin4_eta', xTitle='#eta', yTitle='', normToBinWidth=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_omtf, hName='ugmt_muon_ptmin11_qmin4_eta', xTitle='#eta', yTitle='', normToBinWidth=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_emtf, hName='ugmt_muon_ptmin11_qmin4_eta', xTitle='#eta', yTitle='', normToBinWidth=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_withTfs, hName='ugmt_muon_ptmin18_qmin4_eta', xTitle='#eta', yTitle='', normToBinWidth=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_bmtf, hName='ugmt_muon_ptmin18_qmin4_eta', xTitle='#eta', yTitle='', normToBinWidth=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_omtf, hName='ugmt_muon_ptmin18_qmin4_eta', xTitle='#eta', yTitle='', normToBinWidth=True, data=thisIsData))
    objects.append(plot_hists_comp(hDefs_emtf, hName='ugmt_muon_ptmin18_qmin4_eta', xTitle='#eta', yTitle='', normToBinWidth=True, data=thisIsData))

    ##########################################################################
    # save plots to root file
    if savePlots:
        output = root.TFile('./ugmt_rate_comp_plots.root', 'recreate')
        output.cd()
        for obj in objects:
            c = obj[0]
            c.Write(c.GetName())
            c.Print('./plots/'+c.GetName()+'.pdf', '.pdf')
            c.Print('./plots/'+c.GetName()+'.png', '.png')
        output.Close()

    # wait
    if not batchRun:
        raw_input("Press ENTER to quit.")

if __name__ == "__main__":
    opts = parse_options()
    savePlots = True
    batchRun = True
    plotLegacy = False
    font = 42
    main()

