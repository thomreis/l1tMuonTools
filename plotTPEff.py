#!/usr/bin/env python
from ToolBox import parse_options_and_init_log
# have to do this first or ROOT masks the -h messages
opts, parser = parse_options_and_init_log()

from L1Analysis import L1Ana, L1Ntuple
from analysis_tools.plotting import HistManager, HistManager2d
from analysis_tools.selections import MuonSelections, Matcher
import ROOT as root
import re
import os

def parse_options_plotEff(parser):
    """
    Adds often used options to the OptionParser...
    """
    parsers = parser.add_subparsers()
    sub_parser = parsers.add_parser("plotTPEff")
    sub_parser.add_argument("-i", "--interactive", dest="interactive", action='store_false', help="Draw plots on screen.")
    sub_parser.add_argument("-r", "--run", dest="runnr", default="all_runs", help="Run number for plots.")
    sub_parser.add_argument("-e", "--emul", dest="emul", default=False, action='store_true', help="Use emulator plots.")
    sub_parser.add_argument("--legacy", dest="legacy", default=False, action='store_true', help="Use legacy plots.")
    # histograms to produce
    sub_parser.add_argument("--eff", dest="eff", default=False, action='store_true', help="Plot efficiencies.")
    sub_parser.add_argument("--eff-tf", dest="efftf", default=False, action='store_true', help="Plot efficiencies for TF regions.")
    sub_parser.add_argument("--qualcomp", dest="qualcomp", default=False, action='store_true', help="Plot efficiencies for different qualities.")
    sub_parser.add_argument("--delta", dest="delta", default=False, action='store_true', help="Plot L1 - RECO difference plots.")
    sub_parser.add_argument("--fit-delta", dest="fitdelta", default=False, action='store_true', help="Fit L1 - RECO difference plots.")
    sub_parser.add_argument("--fit-dphi", dest="fitdphi", default=False, action='store_true', help="Fit phi_L1 - phi_RECO difference plot.")
    sub_parser.add_argument("--2d", dest="twod", default=False, action='store_true', help="Plot 2D L1 vs. RECO plots.")
    sub_parser.add_argument("--data-emul", dest="dataemul", default=False, action='store_true', help="Plot data and emulator efficiencies.")
    sub_parser.add_argument("--upgrade-legacy", dest="upgradelegacy", default=False, action='store_true', help="Plot upgrade and legacy efficiencies.")
    sub_parser.add_argument("--by-charge", dest="bycharge", default=False, action='store_true', help="Plot efficiencies by charge.")
    sub_parser.add_argument("--by-tf", dest="bytf", default=False, action='store_true', help="Plot efficiencies by track finder.")
    sub_parser.add_argument("--public", dest="public", default=False, action='store_true', help="Plot style for publication.")
    sub_parser.add_argument("--control", dest="control", default=False, action='store_true', help="Plot control histograms.")
    # options to compare histograms from two files/runs
    sub_parser.add_argument("--fname2", dest="fname2", default=None, help="Second file to take reference histograms from.")
    sub_parser.add_argument("--run2", dest="runnr2", default=None, help="Run number for reference plots. If none is given the same as for the test plot is taken.")
    sub_parser.add_argument("--leg-txt1", dest="legtxt1", default='hm1', help="Legend text for test histograms.")
    sub_parser.add_argument("--leg-txt2", dest="legtxt2", default='hm2', help="Legend text for reference histograms.")
    # plot options
    sub_parser.add_argument("--mc", dest="mc", default=False, action='store_true', help="This is MC.")
    sub_parser.add_argument("--preliminary", dest="prelim", default=False, action='store_true', help="This is a preliminary result.")
    sub_parser.add_argument("--year", dest="year", default=None, help="Year of data taken.")
    sub_parser.add_argument("--lumi", dest="lumi", type=str, default=None, help="Integrated luminosity.")

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
    root.gStyle.SetOptFit(1)
    root.gStyle.SetOptTitle(0)
    root.gStyle.SetNumberContours(99)
    root.gPad.SetTopMargin(0.08)
    root.gPad.SetLeftMargin(0.14)
    root.gPad.SetRightMargin(0.10)
    root.gPad.SetTickx(1)
    root.gPad.SetTicky(1)
    root.gPad.SetGridx(True)
    root.gPad.SetGridy(True)

def plot_2dhist(hm2d, hName, drawDiag=True, clOpts=None):
    canvas_name = 'c_'+hName

    # create canvas and draw on it
    canvas_title = canvas_name
    c = root.TCanvas(canvas_name, canvas_title, 100, 100, 600, 600)
    c.cd()
    set_root_style()
    root.gPad.SetRightMargin(0.14)

    if hName not in hm2d.get_varnames():
        return [c]
    h = hm2d.get(hName).Clone()

    h.GetXaxis().SetTitleFont(font)
    h.GetXaxis().SetLabelFont(font)
    h.GetXaxis().SetLabelSize(fontSize)
    h.GetXaxis().SetNoExponent()
    h.GetYaxis().SetTitleOffset(1.5)
    h.GetYaxis().SetTitleFont(font)
    h.GetYaxis().SetLabelFont(font)
    h.GetYaxis().SetLabelSize(fontSize)
    h.GetZaxis().SetTitleFont(font)
    h.GetZaxis().SetLabelFont(font)
    h.GetZaxis().SetLabelSize(fontSize)
    h.Draw('colz')

    xBase = 0.56
    yBase = 0.13
    notes = extract_notes_from_name(hName, xBase, yBase)

    tex = add_text(notes=notes)
    tex = add_cms_header_text(tex, clOpts, rPos=0.86)

    lines = draw_tf_eta_regions(hName=hName, xMin=h.GetXaxis().GetXmin(), yMin=h.GetYaxis().GetXmin(), xMax=h.GetXaxis().GetXmax(), yMax=h.GetYaxis().GetXmax(), twoD=True, drawDiag=drawDiag)

    c.Modified()
    c.Update()
    
    return [c, h, tex, lines]


def plot_hists(hDefs, xTitle=None, yTitle='# muons', normToBinWidth=False, normToNEntries=0, prefix='', notes=None, autoZoomX=False, xMin=None, xMax=None, logY=False, addOverflow=False, rebin=1, clOpts=None):
    name = prefix
    if hDefs[0]['den']:
        name += hDefs[0]['num']+'_over_'+hDefs[0]['den']
        #legYmin = 0.3-0.04*len(hDefs)
        legYmin = 0.34
    else:
        name = hDefs[0]['num']
        legYmin = 0.81
    if addOverflow:
        name += '_withOverflow'
    if normToBinWidth:
        name += '_normToBinWidth'

    # setup legend according to how many histograms are in the plot
    legend = setup_legend(xMin=0.6, yMin=legYmin, xMax=0.83, yMax=legYmin+0.05*len(hDefs))
    legEntries = []

    if (xMin or xMax or autoZoomX) and addOverflow:
        print 'No overflow bin if xMin/xMax or autoZoomX are set.'
        addOverflow = False

    intL = 0
    intR = 0
    hs = []
    hStack = root.THStack()
    # get all the histograms and set their plot style
    for hDef in hDefs:
        hm = hDef['hm']
        if hDef['den']:
            if hDef['num'] not in hm.get_varnames() or hDef['den'] not in hm.get_varnames():
                print 'Error: ' + hDef['num'] + ' or ' + hDef['den'] + ' not found.'
                continue
            h = hm.get_ratio(hDef['num'], hDef['den'], addoverflow=addOverflow, rebin=rebin).Clone()
        else:
            if hDef['num'] not in hm.get_varnames():
                print 'Error: ' + hDef['num'] + ' not found.'
                continue
            h = hm.get(hDef['num'], addoverflow=addOverflow, rebin=rebin).Clone()
            # find out where the weight of the distributio is to place the legend and notes
            nBinsX = h.GetNbinsX()
            intL += h.Integral(0, nBinsX/3)
            intR += h.Integral(nBinsX/3+1, nBinsX)
            if normToBinWidth:
                for bin in range(1, h.GetNbinsX()+1):
                   h.SetBinContent(bin, h.GetBinContent(bin) / h.GetBinWidth(bin))
                   h.SetBinError(bin, h.GetBinError(bin) / h.GetBinWidth(bin))
            if normToNEntries > 0.:
                integral = h.Integral()
                if integral > 0.:
                    scalefactor = normToNEntries / integral
                else:
                    scalefactor = 1.
                h.Scale(scalefactor)

        h.SetLineColor(hDef['lc'])
        h.SetLineStyle(hDef['ls'])
        h.SetLineWidth(2)
        h.SetMarkerColor(hDef['mc'])
        h.SetMarkerStyle(hDef['ms'])
        h.SetMarkerSize(0.75)
        legStyle = 'lep'
        if hDef['fc']:
            h.SetFillColor(hDef['fc'])
            h.SetLineWidth(1)
            legStyle = 'f'
            # if a fill colour is defined stack this histogram with others
            hStack.Add(h)
        legEntries.append(legend.AddEntry(h, hDef['legtext'], legStyle))
        hs.append(h)

    # replace histograms to be stacked with stack histograms
    if hStack.GetNhists() > 0:
        canvas_name = 'c_'+name+'_stacked'
        stackHistos = hStack.GetStack()
        j = len(stackHistos)-1
        for i, hDef in enumerate(hDefs):
            if hDef['fc']:
                hs[i] = stackHistos[j].Clone()
                j -= 1
    else:
        canvas_name = 'c_'+name

    # create canvas and draw on it
    canvas_title = canvas_name
    canvWidth = 600
    c = root.TCanvas(canvas_name, canvas_title, 100, 100, canvWidth, 600)
    c.cd()
    if (name[-3:] == '.pt' and not hDefs[0]['den']) or logY:
        c.SetLogy(True)

    set_root_style()

    if len(hs) == 0:
        return [c]

    # axis
    drawOpts = ''
    xAxis = hs[0].GetXaxis()
    yAxis = hs[0].GetYaxis()
    maxBinValue = hs[0].GetBinContent(hs[0].GetMaximumBin())
    if xTitle:
        xAxis.SetTitle(xTitle)
    xAxis.SetTitleFont(font)
    xAxis.SetLabelFont(font)
    xAxis.SetLabelSize(fontSize)
    xAxis.SetNoExponent()
    yAxis.SetTitleOffset(1.5)
    yAxis.SetTitleFont(font)
    yAxis.SetTitle(yTitle)
    yAxis.SetLabelFont(font)
    yAxis.SetLabelSize(fontSize)
    # customise axis ranges
    #if name[-7:] == '_dinvpt':
    #    xAxis.SetRangeUser(-0.05, 0.05)
    # find the lowest and highest filled bin
    if autoZoomX:
        minBin = xAxis.GetNbins()
        maxBin = 1
        for hist in hs:
            for b in range(1, hist.GetNbinsX()):
                if hist.GetBinContent(b) == 0:
                    continue
                if b < minBin:
                    minBin = b
                if b > maxBin:
                    maxBin = b
        xAxis.SetRange(minBin, maxBin)
        xAxis.SetLabelSize(0.03)
        xAxis.SetNdivisions(508)
    if xMin and not xMax:
        xAxis.SetRangeUser(xMin, xAxis.GetBinUpEdge(xAxis.GetNbins()))
    elif xMax and not xMin:
        xAxis.SetRangeUser(xAxis.GetBinLowEdge(1), xMax)
    elif xMin and xMax:
        xAxis.SetRangeUser(xMin, xMax)
    yMax = yAxis.GetXmax()
    # set y axis
    if not c.GetLogy():
        #yMax = 1.05*maxBinValue
        if maxBinValue < 1.05:
            yMax = 1.05
            yAxis.SetRangeUser(0., yMax)

    # draw
    hs[0].SetLineWidth(2)
    legEntries[0].SetObject(hs[0])
    legEntries[0].SetOption(legEntries[0].GetOption()+'le')
    hs[0].Draw(drawOpts)
    for h in hs[1:]:
        h.Draw(drawOpts+'same')
    hs[0].Draw(drawOpts+'same')
    hs[0].Draw(drawOpts+'sameaxis')

    # draw vertical lines to mark TF boundaries
    lines = draw_tf_eta_regions(hName=name, yMax=yMax)

    legend.Draw('same')

    tex = add_text(notes=notes, placeRight=(not hDefs[0]['den'] and intL > intR), addOverflow=addOverflow)
    tex = add_cms_header_text(tex, clOpts)

    c.Modified()
    c.Update()

    return [c, hs, legend, lines, tex]


def plot_effs(hDefs, xTitle=None, yTitle='# muons', prefix='', notes=None, autoZoomX=False, xMax=None, addOverflow=False, rebin=1, clOpts=None):
    name = prefix
    name += hDefs[0]['num']+'_over_'+hDefs[0]['den']
    if xMax:
        name += '_xMax{max}'.format(max=xMax)
    if addOverflow:
        name += '_withOverflow'

    draw_legend = False
    # setup legend according to how many histograms are in the plot
    legYmin = 0.34
    legend = setup_legend(xMin=0.5, yMin=legYmin, xMax=0.73, yMax=legYmin+0.05*len(hDefs))
    legEntries = []

    if autoZoomX and addOverflow:
        print 'No overflow bin if autoZoomX is set.'
        addOverflow = False

    drawOpts = 'PZ0'
    effs = []
    effGraphs = []
    # get all the histograms and set their plot style
    for hDef in hDefs:
        hm = hDef['hm']
        if hDef['num'] not in hm.get_varnames() or hDef['den'] not in hm.get_varnames():
            print 'Error: ' + hDef['num'] + ' or ' + hDef['den'] + ' not found.'
            continue
        if addOverflow and xMax != None:
            eff = hm.get_efficiency_int(hDef['num'], hDef['den'], integrateToFromRight=xMax, rebin=rebin, removeProbeBinsNEntriesBelow=0).Clone()
        else:
            eff = hm.get_efficiency(hDef['num'], hDef['den'], addoverflow=addOverflow, rebin=rebin, removeProbeBinsNEntriesBelow=0).Clone()

        eff.SetLineColor(hDef['lc'])
        eff.SetLineStyle(hDef['ls'])
        eff.SetLineWidth(2)
        eff.SetMarkerColor(hDef['mc'])
        eff.SetMarkerStyle(hDef['ms'])
        eff.SetMarkerSize(0.75)
        legStyle = 'lep'
        if hDef['fc']:
            eff.SetFillColor(hDef['fc'])
            eff.SetLineWidth(1)
            legStyle = 'f'
        effs.append(eff)
        effGraphs.append(eff.CreateGraph('A'+drawOpts))
        if hDef['legtext']:
            draw_legend = True
            legEntries.append(legend.AddEntry(effGraphs[-1], hDef['legtext'], legStyle))

    # create canvas and draw on it
    canvas_name = 'c_'+name
    canvas_title = canvas_name
    canvWidth = 600
    c = root.TCanvas(canvas_name, canvas_title, 100, 100, canvWidth, 600)
    c.cd()

    set_root_style()

    # return if no efficiency graph was created
    if len(effGraphs) == 0:
        return [c]

    # axis
    totHist0 = effs[0].GetTotalHistogram()
    xRangeLo = totHist0.GetXaxis().GetBinLowEdge(1)
    xRangeHi = totHist0.GetXaxis().GetBinUpEdge(totHist0.GetNbinsX())
    effGraphs[0].GetXaxis().SetLimits(xRangeLo, xRangeHi)
    axisHisto = effGraphs[0].GetHistogram()
    xAxis = axisHisto.GetXaxis()
    yAxis = axisHisto.GetYaxis()
    if xTitle:
        xAxis.SetTitle(xTitle)
    xAxis.SetTitleFont(font)
    xAxis.SetLabelFont(font)
    xAxis.SetLabelSize(fontSize)
    xAxis.SetNoExponent()
    xAxis.SetRangeUser(xRangeLo, xRangeHi)
    yAxis.SetTitleOffset(1.5)
    yAxis.SetTitleFont(font)
    yAxis.SetTitle(yTitle)
    yAxis.SetLabelFont(font)
    yAxis.SetLabelSize(fontSize)
    # find the lowest and highest filled bin
    if autoZoomX:
        minBin = xAxis.GetNbins()
        maxBin = 1
        for eg in effGraphs:
            xBuff = eg.GetX()
            xBuff.SetSize(eg.GetN())
            xList = list(xBuff)
            if len(xList) > 0:
                minBin = xAxis.FindBin(xList[0])
                maxBin = xAxis.FindBin(xList[-1])
        xAxis.SetRange(minBin, maxBin)
        xAxis.SetLabelSize(0.03)
        xAxis.SetNdivisions(508)
    if xMax:
        xAxis.SetRangeUser(xAxis.GetBinLowEdge(1), xMax)
    # set y axis
    yMax = 1.05
    yAxis.SetRangeUser(0., yMax)

    # draw
    #effGraphs[0].SetLineWidth(2)
    #legEntries[0].SetObject(effGraphs[0])
    effGraphs[0].Draw('A'+drawOpts)
    revEffGraphs = effGraphs
    revEffGraphs.reverse()
    for eg in revEffGraphs:
        eg.Draw(drawOpts+'same')
    axisHisto.Draw('sameaxis')

    # draw vertical lines to mark TF boundaries
    lines = None
    if not clOpts.public:
        lines = draw_tf_eta_regions(hName=name, yMax=yMax)

    if draw_legend:
        legend.Draw('same')

    tex = add_text(notes=notes, addOverflow=addOverflow)
    tex = add_cms_header_text(tex, clOpts)

    c.Modified()
    c.Update()

    return [c, effGraphs, axisHisto, legend, lines, tex]


def setup_legend(xMin=0.5, yMin=0.34, xMax=0.73, yMax=0.39):
    legend = root.TLegend(xMin, yMin, xMax, yMax)
    legend.SetTextFont(font)
    legend.SetTextSize(0.035)
    legend.SetBorderSize(0)
    legend.SetFillColor(19)
    legend.SetFillStyle(0)
    #legend.SetNColumns(2)
    return legend


# add text to plots 
def add_text(notes=None, placeRight=False, addOverflow=False):
    tex = root.TLatex()
    tex.SetNDC()
    tex.SetTextFont(font)
    tex.SetTextSize(fontSize)
    if notes:
        xOffset = 0
        yOffset = 0
        if placeRight:
            xOffset = 0.45
            yOffset = -0.06
        tex.SetTextSize(0.035)
        for note in notes:
            if note[3]: # defines if changing the note position is allowed
                tex.DrawLatex(note[0]+xOffset, note[1]+yOffset, note[2])
            else:
                tex.DrawLatex(note[0], note[1], note[2])
    #if addOverflow:
    #    tex.SetTextAngle(90)
    #    tex.DrawLatex(0.87, 0.20, 'highest bin includes overflow')
    #    tex.SetTextAngle(0)
    return tex

# add the CMS header text to the plot
def add_cms_header_text(tex=None, clOpts=None, lPos=0.14, rPos=0.9):
    if not tex:
        tex = root.TLatex()
        tex.SetNDC()

    energy='13 TeV'

    cmsTextFont = 61
    cmsTextSize = 0.05
    extraTextFont = 52
    extraOverCmsTextSize = 0.8

    cmsText = 'CMS'
    extraText = ''
    lumiText = ''
    if clOpts:
        if clOpts.prelim:
            extraText = "preliminary"
        elif not clOpts.public:
            extraText   = "internal"

        lumiText = ''
        if clOpts.year:
            lumiText += '{y}'.format(y=clOpts.year)
        if clOpts.mc:
            lumiText += ' simulation'
        else:
            lumiText += ' data'
        if clOpts.lumi:
            lumiText += '  {l}'.format(l=clOpts.lumi)
        lumiText += ' ('+energy+')'
 
    tex.SetTextFont(cmsTextFont)
    tex.SetTextSize(cmsTextSize)
    tex.DrawLatex(lPos, 0.93, cmsText)
    tex.SetTextFont(extraTextFont)
    tex.SetTextSize(cmsTextSize*extraOverCmsTextSize)
    tex.DrawLatex(lPos+0.12, 0.93, extraText)
    tex.SetTextFont(font)
    tex.SetTextAlign(31)
    tex.DrawLatex(rPos, 0.93, lumiText)
    return tex

# draw vertical lines to mark TF boundaries
def draw_tf_eta_regions(hName='', xMin=0., yMin=0., xMax=1., yMax=1., twoD=False, drawDiag=False):
    lines = []
    if hName[-4:] == '_eta':
        lines.append(root.TLine(0.83, yMin, 0.83, yMax))
        lines[-1].SetLineStyle(root.kDashed)
        lines[-1].Draw('same')
        lines.append(root.TLine(-0.83, yMin, -0.83, yMax))
        lines[-1].SetLineStyle(root.kDashed)
        lines[-1].Draw('same')
        lines.append(root.TLine(1.24, yMin, 1.24, yMax))
        lines[-1].SetLineStyle(root.kDashed)
        lines[-1].Draw('same')
        lines.append(root.TLine(-1.24, yMin, -1.24, yMax))
        lines[-1].SetLineStyle(root.kDashed)
        lines[-1].Draw('same')

        if twoD:
            lines.append(root.TLine(xMin, 0.83, xMax, 0.83))
            lines[-1].SetLineStyle(root.kDashed)
            lines[-1].Draw('same')
            lines.append(root.TLine(xMin, -0.83, xMax, -0.83))
            lines[-1].SetLineStyle(root.kDashed)
            lines[-1].Draw('same')
            lines.append(root.TLine(xMin, 1.24, xMax, 1.24))
            lines[-1].SetLineStyle(root.kDashed)
            lines[-1].Draw('same')
            lines.append(root.TLine(xMin, -1.24, xMax, -1.24))
            lines[-1].SetLineStyle(root.kDashed)
            lines[-1].Draw('same')

    if drawDiag:
        lines.append(root.TLine(xMin, yMin, xMax, yMax))
        lines[-1].SetLineStyle(root.kSolid)
        lines[-1].SetLineColor(root.kMagenta)
        lines[-1].Draw('same')
    return lines
    

def hist_styles(stacked=False):
    styles = {}
    styles['gmt'] = {'lc':root.kCyan, 'ls':root.kSolid, 'fc':None, 'mc':root.kCyan, 'ms':root.kCircle, 'legtext':'GMT'}
    styles['ugmt'] = {'lc':root.kRed, 'ls':root.kSolid, 'fc':None, 'mc':root.kRed, 'ms':root.kDot, 'legtext':'uGMT'}
    styles['data'] = {'lc':root.kRed, 'ls':root.kSolid, 'fc':None, 'mc':root.kRed, 'ms':root.kFullTriangleUp, 'legtext':'data'}
    styles['data_q12'] = {'lc':root.kRed, 'ls':root.kSolid, 'fc':None, 'mc':root.kRed, 'ms':root.kOpenTriangleUp, 'legtext':'data Q #geq 12'}
    styles['data_q8'] = {'lc':root.kBlue, 'ls':root.kSolid, 'fc':None, 'mc':root.kBlue, 'ms':root.kOpenCircle, 'legtext':'data Q #geq 8'}
    styles['data_q4'] = {'lc':root.kGreen, 'ls':root.kSolid, 'fc':None, 'mc':root.kGreen, 'ms':root.kFullSquare, 'legtext':'data Q #geq 4'}
    styles['data_alleta'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':None, 'mc':root.kBlack, 'ms':root.kFullCircle, 'legtext':'0 < |#eta^{reco}| < 2.4'}
    styles['data_bmtfeta'] = {'lc':root.kBlue, 'ls':root.kSolid, 'fc':None, 'mc':root.kBlue, 'ms':root.kOpenCircle, 'legtext':'0 < |#eta^{reco}| < 0.83'}
    styles['data_omtfeta'] = {'lc':root.kGreen+1, 'ls':root.kSolid, 'fc':None, 'mc':root.kGreen+1, 'ms':root.kOpenSquare, 'legtext':'0.83 < |#eta^{reco}| < 1.24'}
    styles['data_emtfeta'] = {'lc':root.kRed, 'ls':root.kSolid, 'fc':None, 'mc':root.kRed, 'ms':root.kOpenTriangleUp, 'legtext':'1.24 < |#eta^{reco}| < 2.4'}
    styles['emul'] = {'lc':root.kBlue, 'ls':root.kSolid, 'fc':None, 'mc':root.kBlue, 'ms':root.kFullCircle, 'legtext':'emul'}
    styles['emul_q12'] = {'lc':root.kRed, 'ls':root.kSolid, 'fc':None, 'mc':root.kRed, 'ms':root.kOpenTriangleUp, 'legtext':'emul Q #geq 12'}
    styles['emul_q8'] = {'lc':root.kBlue, 'ls':root.kSolid, 'fc':None, 'mc':root.kBlue, 'ms':root.kOpenCircle, 'legtext':'emul Q #geq 8'}
    styles['emul_q4'] = {'lc':root.kGreen, 'ls':root.kSolid, 'fc':None, 'mc':root.kGreen, 'ms':root.kFullSquare, 'legtext':'emul Q #geq 4'}
    styles['upgrade'] = {'lc':root.kRed, 'ls':root.kSolid, 'fc':None, 'mc':root.kRed, 'ms':root.kFullTriangleUp, 'legtext':'upgrade'}
    styles['upgrade_q12'] = {'lc':root.kRed, 'ls':root.kSolid, 'fc':None, 'mc':root.kRed, 'ms':root.kOpenTriangleUp, 'legtext':'upgrade Q #geq 12'}
    styles['upgrade_q8'] = {'lc':root.kBlue, 'ls':root.kSolid, 'fc':None, 'mc':root.kBlue, 'ms':root.kOpenCircle, 'legtext':'upgrade Q #geq 8'}
    styles['upgrade_q4'] = {'lc':root.kGreen, 'ls':root.kSolid, 'fc':None, 'mc':root.kGreen, 'ms':root.kFullSquare, 'legtext':'upgrade Q #geq 4'}
    styles['legacy'] = {'lc':root.kBlue, 'ls':root.kSolid, 'fc':None, 'mc':root.kBlue, 'ms':root.kFullCircle, 'legtext':'legacy'}
    styles['legacy_q5'] = {'lc':root.kRed, 'ls':root.kSolid, 'fc':None, 'mc':root.kRed, 'ms':root.kOpenTriangleUp, 'legtext':'legacy Q #geq 5'}
    styles['legacy_q3'] = {'lc':root.kBlue, 'ls':root.kSolid, 'fc':None, 'mc':root.kBlue, 'ms':root.kOpenCircle, 'legtext':'legacy Q #geq 3'}
    styles['legacy_q2'] = {'lc':root.kGreen, 'ls':root.kSolid, 'fc':None, 'mc':root.kGreen, 'ms':root.kFullSquare, 'legtext':'legacy Q #geq 2'}
    styles['data_pub'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':None, 'mc':root.kBlack, 'ms':root.kFullCircle, 'legtext':''}
    styles['data_0'] = {'lc':root.kRed, 'ls':root.kSolid, 'fc':None, 'mc':root.kRed, 'ms':root.kOpenTriangleUp, 'legtext':'L1 p_{T} #geq XX GeV'}
    styles['data_1'] = {'lc':root.kBlue, 'ls':root.kSolid, 'fc':None, 'mc':root.kBlue, 'ms':root.kOpenCircle, 'legtext':'L1 p_{T} #geq XX GeV'}
    styles['data_2'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':None, 'mc':root.kBlack, 'ms':root.kOpenSquare, 'legtext':'L1 p_{T} #geq XX GeV'}
    styles['emul_pub'] = {'lc':root.kBlue, 'ls':root.kSolid, 'fc':None, 'mc':root.kBlue, 'ms':root.kFullCircle, 'legtext':''}
    styles['old_sel'] = {'lc':root.kRed, 'ls':root.kSolid, 'fc':None, 'mc':root.kRed, 'ms':root.kFullTriangleUp, 'legtext':'old selection'}
    styles['new_sel'] = {'lc':root.kBlue, 'ls':root.kSolid, 'fc':None, 'mc':root.kBlue, 'ms':root.kFullCircle, 'legtext':'new selection'}
    styles['pos_charge'] = {'lc':root.kRed, 'ls':root.kSolid, 'fc':None, 'mc':root.kRed, 'ms':root.kFullTriangleUp, 'legtext':'pos. charge'}
    styles['neg_charge'] = {'lc':root.kBlue, 'ls':root.kSolid, 'fc':None, 'mc':root.kBlue, 'ms':root.kFullCircle, 'legtext':'neg. charge'}
    styles['data_all'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':None, 'mc':root.kBlack, 'ms':root.kOpenSquare, 'legtext':'data'}
    styles['data_bmtf'] = {'lc':root.kBlue, 'ls':root.kSolid, 'fc':None, 'mc':root.kBlue, 'ms':root.kFullCircle, 'legtext':'data BMTF'}
    styles['data_omtf'] = {'lc':root.kGreen, 'ls':root.kSolid, 'fc':None, 'mc':root.kGreen, 'ms':root.kFullTriangleUp, 'legtext':'data OMTF'}
    styles['data_emtf'] = {'lc':root.kRed, 'ls':root.kSolid, 'fc':None, 'mc':root.kRed, 'ms':root.kFullTriangleDown, 'legtext':'data EMTF'}
    styles['emul_all'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':None, 'mc':root.kBlack, 'ms':root.kFullSquare, 'legtext':'emul'}
    styles['emul_bmtf'] = {'lc':root.kBlue, 'ls':root.kSolid, 'fc':None, 'mc':root.kBlue, 'ms':root.kOpenCircle, 'legtext':'emul BMTF'}
    styles['emul_omtf'] = {'lc':root.kGreen, 'ls':root.kSolid, 'fc':None, 'mc':root.kGreen, 'ms':root.kOpenTriangleUp, 'legtext':'emul OMTF'}
    styles['emul_emtf'] = {'lc':root.kRed, 'ls':root.kSolid, 'fc':None, 'mc':root.kRed, 'ms':root.kOpenTriangleDown, 'legtext':'emul EMTF'}
    styles['hm1'] = {'lc':root.kRed, 'ls':root.kSolid, 'fc':None, 'mc':root.kRed, 'ms':root.kFullTriangleUp, 'legtext':'hm1'}
    styles['hm2'] = {'lc':root.kBlue, 'ls':root.kSolid, 'fc':None, 'mc':root.kBlue, 'ms':root.kFullCircle, 'legtext':'hm2'}
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
        styles['bmtf_ugmt_reg'] = {'lc':root.kBlue-4, 'ls':root.kSolid, 'fc':None, 'legtext':'BMTF uGMT (0 < |#eta^{reco}| < 0.83)'}
        styles['omtf_ugmt_reg'] = {'lc':root.kGreen-4, 'ls':root.kSolid, 'fc':None, 'legtext':'OMTF uGMT (0.83 < |#eta^{reco}| < 1.24)'}
        styles['emtf_ugmt_reg'] = {'lc':root.kRed-4, 'ls':root.kSolid, 'fc':None, 'legtext':'EMTF uGMT (1.24 < |#eta^{reco}| < 2.5)'}
        styles['bmtf_reg'] = {'lc':root.kBlue-4, 'ls':root.kSolid, 'fc':None, 'legtext':'BMTF (0 < |#eta^{reco}| < 0.83)'}
        styles['omtf_reg'] = {'lc':root.kGreen-4, 'ls':root.kSolid, 'fc':None, 'legtext':'OMTF (0.83 < |#eta^{reco}| < 1.24)'}
        styles['emtf_reg'] = {'lc':root.kRed-4, 'ls':root.kSolid, 'fc':None, 'legtext':'EMTF (1.24 < |#eta^{reco}| < 2.5)'}
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

def extract_notes_from_name(name, xBase, yBase, etaTxt=True, qualTxt=True, ptTxt=True, public=False):
    notes = []
    # extract eta range
    if etaTxt and name.rfind('_eta') == -1:
        eta_number_strs = re.findall(r'[\d\p\d]+', name[name.find('EtaMin')+6:name.find('EtaMax')+12])
        if len(eta_number_strs) > 1:
            if public:
                if len(eta_number_strs[0]) > 3:
                    eta_number_strs[0] = eta_number_strs[0][0:3]
                if len(eta_number_strs[1]) > 3:
                    eta_number_strs[1] = eta_number_strs[1][0:3]
                note_str = eta_number_strs[0].replace('p', '.')+' < |#eta| < '+eta_number_strs[1].replace('p', '.')
            else:
                note_str = eta_number_strs[0].replace('p', '.')+' < |#eta^{reco}| < '+eta_number_strs[1].replace('p', '.')
            notes.append([xBase, yBase+0.15, note_str, True])
    # extract quality
    if qualTxt and name.find('qualMin') != -1:
        qualPos = name.find('qualMin')
        qual_strs = re.findall(r'[\d\p\d]+', name[qualPos+7:qualPos+9])
        if len(qual_strs) > 0:
            if public:
                if qual_strs[0] == '12':
                    qual_note_str = 'Tight L1 quality'
                elif qual_strs[0] == '8':
                    qual_note_str = 'Loose L1 quality'
                elif qual_strs[0] == '4':
                    qual_note_str = 'Very loose L1 quality'
                else:
                    qual_note_str = ''
            else:
                qual_note_str = 'Quality #geq '+qual_strs[0].replace('p', '.')
            notes.append([xBase, yBase+0.1, qual_note_str, True])
    # extract pt range
    if ptTxt:
        overPos = name.find('_over_')
        if overPos != -1:
            ptminPos = name.find('ptmin')
            l1_ptmin_strs = re.findall(r'\d+\p?\d*', name[ptminPos+5:ptminPos+9])
            if len(l1_ptmin_strs) > 0:
                if l1_ptmin_strs[0][-1] == '.':
                    l1_ptmin_strs[0] = l1_ptmin_strs[0][0:-1]
                if public:
                    notes.append([xBase, yBase+0.05, 'L1 p_{T} #geq '+l1_ptmin_strs[0].replace('p', '.')+' GeV', True])
                else:
                    notes.append([xBase, yBase+0.05, 'p_{T}^{L1} #geq '+l1_ptmin_strs[0].replace('p', '.')+' GeV', True])
            ptminPos = overPos + name[overPos:].rfind('ptmin')
            probe_ptmin_strs = re.findall(r'\d+\p?\d*', name[ptminPos+5:ptminPos+9])
            if len(probe_ptmin_strs) > 0:
                if probe_ptmin_strs[0][-1] == 'p':
                    probe_ptmin_strs[0] = probe_ptmin_strs[0][0:-1]
                if public:
                    notes.append([xBase, yBase, 'Offline p_{T} > '+probe_ptmin_strs[0].replace('p', '.')+' GeV', True])
                else:
                    notes.append([xBase, yBase, 'p_{T}^{reco} > '+probe_ptmin_strs[0].replace('p', '.')+' GeV', True])
        elif name.find('_matched_') != -1:
            substrings = name.split('_matched_')
            for substr in substrings:
                if substr.find('probe') != -1 and not substr[-3:] == '_pt':
                    ptminPos = substr.find('ptmin')
                    ptmin_strs = re.findall(r'\d+\p?\d*', substr[ptminPos+5:ptminPos+9])
                    if len(ptmin_strs) > 0:
                        if ptmin_strs[0][-1] == 'p':
                            ptmin_strs[0] = ptmin_strs[0][0:-1]
                        if public:
                            if ptmin_strs[0][0] != '0':
                                notes.append([xBase, yBase, 'Offline p_{T} > '+ptmin_strs[0].replace('p', '.')+' GeV', True])
                        else:
                            notes.append([xBase, yBase, 'p_{T}^{reco} > '+ptmin_strs[0].replace('p', '.')+' GeV', True])
                    ptmaxPos = substr.find('ptmax')
                    ptmax_strs = re.findall(r'\d+\p?\d*', substr[ptmaxPos+5:ptmaxPos+8])
                    if len(ptmax_strs) > 0:
                        if ptmax_strs[0][-1] == 'p':
                            ptmax_strs[0] = ptmax_strs[0][0:-1]
                        if public:
                            if ptmax_strs[0][0] != '0':
                                notes.append([xBase, yBase-0.05, 'Offline p_{T} < '+ptmax_strs[0].replace('p', '.')+' GeV', True])
                        else:
                            notes.append([xBase, yBase-0.05, 'p_{T}^{reco} < '+ptmax_strs[0].replace('p', '.')+' GeV', True])
                if substr.find('l1_muon') != -1:
                    ptminPos = substr.find('ptmin')
                    ptmin_strs = re.findall(r'\d+\p?\d*', substr[ptminPos+5:ptminPos+9])
                    if len(ptmin_strs) > 0:
                        if ptmin_strs[0][-1] == 'p':
                            ptmin_strs[0] = ptmin_strs[0][0:-1]
                        if public:
                            notes.append([xBase, yBase+0.05, 'L1 p_{T} #geq '+ptmin_strs[0].replace('p', '.')+' GeV', True])
                        else:
                            notes.append([xBase, yBase+0.05, 'p_{T}^{L1} #geq '+ptmin_strs[0].replace('p', '.')+' GeV', True])
    return notes

# plot histogram
def plot_hists_standard(hm, hName, den=None, hNamePrefix='', xTitle='', yTitle='', stacked=False, normToBinWidth=False, normToNEntries=0, autoZoomX=False, xMin=None, xMax=None, logY=False, addOverflow=False, rebin=1, clOpts=None):
    styles = hist_styles(stacked)

    h_dict = {'hm':hm, 'num':hNamePrefix+hName, 'den':den}
    h_dict.update(styles['data'])
    hDefs = [h_dict]

    xBase = 0.17
    yBase = 0.71
    if den:
        xBase = 0.5
        yBase = 0.13
    notes = extract_notes_from_name(hName, xBase, yBase)

    return plot_hists(hDefs, xTitle, yTitle, normToBinWidth, normToNEntries, '', notes, autoZoomX, xMin, xMax, logY, addOverflow, rebin, clOpts)

# plot data and emulator histogram
def plot_hists_data_emul(hm, hName, den=None, hNamePrefix='', xTitle='', yTitle='', normToBinWidth=False, normToNEntries=0, autoZoomX=False, xMin=None, xMax=None, logY=False, addOverflow=False, rebin=1, clOpts=None):
    styles = hist_styles(False)

    data_dict = {'hm':hm, 'num':hNamePrefix+hName, 'den':den}
    data_dict.update(styles['data'])
    emul_dict = {'hm':hm, 'num':'emu_'+hNamePrefix+hName, 'den':den}
    emul_dict.update(styles['emul'])
    hDefs = [data_dict]
    hDefs.append(emul_dict)

    xBase = 0.17
    yBase = 0.71
    if den:
        xBase = 0.5
        yBase = 0.13
    notes = extract_notes_from_name(hName, xBase, yBase)

    return plot_hists(hDefs, xTitle, yTitle, normToBinWidth, normToNEntries, 'dataEmulHisto_', notes, autoZoomX, xMin, xMax, logY, addOverflow, rebin, clOpts)

# plot efficiency
def plot_eff_standard(hm, hName, den, hNamePrefix='', xTitle='', yTitle='', autoZoomX=False, xMax=None, addOverflow=False, rebin=1, clOpts=None):
    styles = hist_styles(False)

    styleKey = 'data'
    denPrefix = ''
    if clOpts:
        if clOpts.emul:
            hNamePrefix = 'emu_'+hNamePrefix
            denPrefix = 'emu_'
            styleKey = 'emul'
        if clOpts.legacy:
            hName = hName.replace('qualMin12', 'qualMin5')
            hName = hName.replace('qualMin8', 'qualMin3')
            hName = hName.replace('qualMin4', 'qualMin2')
            hNamePrefix = 'legacy_'+hNamePrefix
            denPrefix = 'legacy_'
            styleKey = 'legacy'

    h_dict = {'hm':hm, 'num':hNamePrefix+hName, 'den':denPrefix+den}
    h_dict.update(styles[styleKey])
    hDefs = [h_dict]

    xBase = 0.5
    yBase = 0.13
    notes = extract_notes_from_name(hName, xBase, yBase)

    return plot_effs(hDefs, xTitle, yTitle, '', notes, autoZoomX, xMax, addOverflow, rebin, clOpts)

# plot public style efficiency
def plot_eff_public(hm, hName, den, hNamePrefix='', xTitle='', yTitle='', autoZoomX=False, xMax=None, addOverflow=False, rebin=1, clOpts=None):
    styles = hist_styles(False)

    styleKey = 'data_pub'
    denPrefix = ''
    if clOpts:
        if clOpts.emul:
            hNamePrefix = 'emu_'+hNamePrefix
            denPrefix = 'emu_'
            styleKey = 'emul_pub'

    h_dict = {'hm':hm, 'num':hNamePrefix+hName, 'den':denPrefix+den}
    h_dict.update(styles[styleKey])
    hDefs = [h_dict]

    xBase = 0.6
    yBase = 0.14
    notes = extract_notes_from_name(hName, xBase, yBase, etaTxt=True, qualTxt=True, ptTxt=True, public=True)

    return plot_effs(hDefs, xTitle, yTitle, '', notes, autoZoomX, xMax, addOverflow, rebin, clOpts)

# plot public style efficiencies for thresholds 18 GeV, 22 GeV, and 25 GeV
def plot_eff_public_pt(hm, hName, den, hNamePrefix='', ptmins=[], xTitle='', yTitle='', autoZoomX=False, xMax=None, addOverflow=False, rebin=1, clOpts=None):
    styles = hist_styles(False)

    styleKey = 'data'
    denPrefix = ''
    if clOpts:
        if clOpts.emul:
            hNamePrefix = 'emu_'+hNamePrefix
            denPrefix = 'emu_'
            styleKey = 'emul'

    hDefs = []
    for i, ptmin in enumerate(ptmins):
        ptmin_dict = {'hm':hm, 'num':hNamePrefix+hName.replace('ptminXX_', 'ptmin'+str(ptmin).replace('.', 'p')+'_'), 'den':denPrefix+den}
        ptmin_dict.update(styles[styleKey+'_'+str(i)])
        ptmin_dict['legtext'] = ptmin_dict['legtext'].replace('XX', str(ptmin))
        hDefs.append(ptmin_dict)

    xBase = 0.558
    yBase = 0.07
    notes = extract_notes_from_name(hName, xBase, yBase, etaTxt=True, qualTxt=True, ptTxt=False, public=True)

    return plot_effs(hDefs, xTitle, yTitle, 'ptRange_', notes, autoZoomX, xMax, addOverflow, rebin, clOpts)

# plot data and emulator efficiencies
def plot_eff_data_emul(hm, hName, den, hNamePrefix='', xTitle='', yTitle='', autoZoomX=False, xMax=None, addOverflow=False, rebin=1, clOpts=None):
    styles = hist_styles(False)

    data_dict = {'hm':hm, 'num':hNamePrefix+hName, 'den':den}
    data_dict.update(styles['data'])
    emul_dict = {'hm':hm, 'num':'emu_'+hNamePrefix+hName, 'den':'emu_'+den}
    emul_dict.update(styles['emul'])
    #data_dict.update(styles['old_sel'])
    #emul_dict = {'num':'newSel_'+hNamePrefix+hName, 'den':'newSel_'+den}
    #emul_dict.update(styles['new_sel'])
    hDefs = [data_dict]
    hDefs.append(emul_dict)

    xBase = 0.5
    yBase = 0.13
    notes = extract_notes_from_name(hName, xBase, yBase)

    return plot_effs(hDefs, xTitle, yTitle, 'dataEmul_', notes, autoZoomX, xMax, addOverflow, rebin, clOpts)
    #return plot_effs(hDefs, xTitle, yTitle, 'selComp_', notes, autoZoomX, xMax, addOverflow, rebin, clOpts)

# plot data and emulator efficiencies
def plot_eff_upgrade_legacy(hm, hName, den, hNamePrefix='', xTitle='', yTitle='', autoZoomX=False, xMax=None, addOverflow=False, rebin=1, clOpts=None):
    styles = hist_styles(False)

    upgrade_dict = {'hm':hm, 'num':hNamePrefix+hName, 'den':den}
    upgrade_dict.update(styles['upgrade'])
    hNameLegacy = hName
    hNameLegacy = hNameLegacy.replace('qualMin12', 'qualMin5')
    hNameLegacy = hNameLegacy.replace('qualMin8', 'qualMin3')
    hNameLegacy = hNameLegacy.replace('qualMin4', 'qualMin2')
    legacy_dict = {'hm':hm, 'num':'legacy_'+hNamePrefix+hNameLegacy, 'den':'legacy_'+den}
    legacy_dict.update(styles['legacy'])
    #data_dict.update(styles['old_sel'])
    #emul_dict = {'num':'newSel_'+hNamePrefix+hName, 'den':'newSel_'+den}
    #emul_dict.update(styles['new_sel'])
    hDefs = [upgrade_dict]
    hDefs.append(legacy_dict)

    xBase = 0.5
    yBase = 0.13
    notes = extract_notes_from_name(hName, xBase, yBase)

    return plot_effs(hDefs, xTitle, yTitle, 'upgradeLegacy_', notes, autoZoomX, xMax, addOverflow, rebin, clOpts)
    #return plot_effs(hDefs, xTitle, yTitle, 'selComp_', notes, autoZoomX, xMax, addOverflow, rebin, clOpts)

# plot efficiencies for qualities 12, 8, and 4
def plot_eff_qual(hm, hName, den, hNamePrefix='', xTitle='', yTitle='', autoZoomX=False, xMax=None, addOverflow=False, rebin=1, clOpts=None):
    styles = hist_styles(False)

    styleKey = 'data'
    denPrefix = ''
    if clOpts:
        if clOpts.emul:
            hNamePrefix = 'emu_'+hNamePrefix
            denPrefix = 'emu_'
            styleKey = 'emul'

    q12_dict = {'hm':hm, 'num':hNamePrefix+hName.replace('qualMinXX_', 'qualMin12_'), 'den':denPrefix+den}
    q12_dict.update(styles[styleKey+'_q12'])
    q8_dict = {'hm':hm, 'num':hNamePrefix+hName.replace('qualMinXX_', 'qualMin8_'), 'den':denPrefix+den}
    q8_dict.update(styles[styleKey+'_q8'])
    q4_dict = {'hm':hm, 'num':hNamePrefix+hName.replace('qualMinXX_', 'qualMin4_'), 'den':denPrefix+den}
    q4_dict.update(styles[styleKey+'_q4'])
    hDefs = [q12_dict]
    hDefs.append(q8_dict)
    hDefs.append(q4_dict)

    xBase = 0.5
    yBase = 0.13
    notes = extract_notes_from_name(hName, xBase, yBase)

    return plot_effs(hDefs, xTitle, yTitle, 'qual_', notes, autoZoomX, xMax, addOverflow, rebin, clOpts)

# plot efficiencies for eta ranges on one plot
def plot_eff_eta_ranges(hm, hName, den, hNamePrefix='', xTitle='', yTitle='', autoZoomX=False, xMax=None, addOverflow=False, rebin=1, clOpts=None):
    styles = hist_styles(False)

    styleKey = 'data'
    denPrefix = ''
    if clOpts:
        if clOpts.emul:
            hNamePrefix = 'emu_'+hNamePrefix
            denPrefix = 'emu_'
            styleKey = 'emul'

    all_dict = {'hm':hm, 'num':hNamePrefix+hName.replace('absEtaMinXX', 'absEtaMin0').replace('absEtaMaxYY', 'absEtaMax2p4'), 'den':denPrefix+den.replace('absEtaMinXX', 'absEtaMin0').replace('absEtaMaxYY', 'absEtaMax2p4')}
    all_dict.update(styles[styleKey+'_alleta'])
    bmtfreg_dict = {'hm':hm, 'num':hNamePrefix+hName.replace('absEtaMinXX', 'absEtaMin0').replace('absEtaMaxYY', 'absEtaMax0p83'), 'den':denPrefix+den.replace('absEtaMinXX', 'absEtaMin0').replace('absEtaMaxYY', 'absEtaMax0p83')}
    bmtfreg_dict.update(styles[styleKey+'_bmtfeta'])
    omtfreg_dict = {'hm':hm, 'num':hNamePrefix+hName.replace('absEtaMinXX', 'absEtaMin0p83').replace('absEtaMaxYY', 'absEtaMax1p24'), 'den':denPrefix+den.replace('absEtaMinXX', 'absEtaMin0p83').replace('absEtaMaxYY', 'absEtaMax1p24')}
    omtfreg_dict.update(styles[styleKey+'_omtfeta'])
    emtfreg_dict = {'hm':hm, 'num':hNamePrefix+hName.replace('absEtaMinXX', 'absEtaMin1p24').replace('absEtaMaxYY', 'absEtaMax2p4'), 'den':denPrefix+den.replace('absEtaMinXX', 'absEtaMin1p24').replace('absEtaMaxYY', 'absEtaMax2p4')}
    emtfreg_dict.update(styles[styleKey+'_emtfeta'])
    hDefs = []
    hDefs.append(bmtfreg_dict)
    hDefs.append(omtfreg_dict)
    hDefs.append(emtfreg_dict)
    hDefs.append(all_dict)

    xBase = 0.5
    yBase = 0.13
    notes = extract_notes_from_name(hName, xBase, yBase)

    return plot_effs(hDefs, xTitle, yTitle, 'efftf_', notes, autoZoomX, xMax, addOverflow, rebin, clOpts)

# plot efficiency contribution from each TF and overall efficiency
def plot_eff_tf(hm, hName, den, hNamePrefix='', xTitle='', yTitle='', autoZoomX=False, xMax=None, addOverflow=False, rebin=1, clOpts=None):
    styles = hist_styles(False)

    styleKey = 'data'
    denPrefix = ''
    if clOpts:
        if clOpts.emul:
            hNamePrefix = 'emu_'+hNamePrefix
            denPrefix = 'emu_'
            styleKey = 'emul'

    all_dict = {'hm':hm, 'num':hNamePrefix+hName, 'den':denPrefix+den}
    all_dict.update(styles[styleKey+'_all'])
    bmtf_dict = {'hm':hm, 'num':'bmtf_only_'+hNamePrefix+hName, 'den':denPrefix+'bmtf_only_'+den}
    bmtf_dict.update(styles[styleKey+'_bmtf'])
    omtf_dict = {'hm':hm, 'num':'omtf_only_'+hNamePrefix+hName, 'den':denPrefix+'omtf_only_'+den}
    omtf_dict.update(styles[styleKey+'_omtf'])
    emtf_dict = {'hm':hm, 'num':'emtf_only_'+hNamePrefix+hName, 'den':denPrefix+'emtf_only_'+den}
    emtf_dict.update(styles[styleKey+'_emtf'])
    hDefs = [all_dict]
    hDefs.append(bmtf_dict)
    hDefs.append(omtf_dict)
    hDefs.append(emtf_dict)

    xBase = 0.5
    yBase = 0.13
    notes = extract_notes_from_name(hName, xBase, yBase)

    return plot_effs(hDefs, xTitle, yTitle, 'tf_', notes, autoZoomX, xMax, addOverflow, rebin, clOpts)

# plot efficiencies for positive and negative charge
def plot_eff_charge(hm, hName, den, hNamePrefix='', xTitle='', yTitle='', autoZoomX=False, xMax=None, addOverflow=False, rebin=1, clOpts=None):
    styles = hist_styles(False)

    pos_dict = {'hm':hm, 'num':'posCharge'+hNamePrefix+hName, 'den':'posCharge'+den}
    pos_dict.update(styles['pos_charge'])
    neg_dict = {'hm':hm, 'num':'negCharge'+hNamePrefix+hName, 'den':'negCharge'+den}
    neg_dict.update(styles['neg_charge'])
    hDefs = [pos_dict]
    hDefs.append(neg_dict)

    xBase = 0.5
    yBase = 0.13
    notes = extract_notes_from_name(hName, xBase, yBase)

    return plot_effs(hDefs, xTitle, yTitle, 'chargeComp_', notes, autoZoomX, xMax, addOverflow, rebin, clOpts)

# plot efficiencies with the same name from two different histogram managers (files, runs)
def plot_eff_hm_comp(hm1, hm2, hName, den, hNamePrefix='', xTitle='', yTitle='', autoZoomX=False, xMax=None, addOverflow=False, rebin=1, clOpts=None):
    styles = hist_styles(False)

    denPrefix = ''
    if clOpts:
        if clOpts.emul:
            hNamePrefix = 'emu_'+hNamePrefix
            denPrefix = 'emu_'
        legTxt1 = clOpts.legtxt1
        legTxt2 = clOpts.legtxt2
    else:
        legTxt1 = ''
        legTxt2 = ''

    hm1_dict = {'hm':hm1, 'num':hNamePrefix+hName, 'den':denPrefix+den}
    hm1_dict.update(styles['hm1'])
    hm1_dict['legtext'] = legTxt1
    hm2_dict = {'hm':hm2, 'num':hNamePrefix+hName, 'den':denPrefix+den}
    hm2_dict.update(styles['hm2'])
    hm2_dict['legtext'] = legTxt2
    hDefs = [hm1_dict]
    hDefs.append(hm2_dict)

    xBase = 0.5
    yBase = 0.13
    notes = extract_notes_from_name(hName, xBase, yBase)

    return plot_effs(hDefs, xTitle, yTitle, 'comp_', notes, autoZoomX, xMax, addOverflow, rebin, clOpts)

def fitHisto(h, fitfunc):
    f = root.TF1('f', fitfunc)
    if h.GetEntries() > 10: # need some entries at least
        f.SetLineWidth(1)
        f.SetLineColor(root.kBlue)
        peak = h.GetBinCenter(h.GetMaximumBin())
        f.SetParameter(1, peak)
        f.SetParameter(2, h.GetRMS())
        f.SetParameter(3, -1.)
        f.SetParameter(4, 5000.)
        factorFit1 = 1.
        factorFit2 = 4
        # raw and fine fit
        h.Fit(f, '', '', peak-factorFit1*h.GetRMS(), peak+factorFit1*h.GetRMS())
        h.Fit(f, '', '', f.GetParameter(1)-factorFit2*f.GetParameter(2), f.GetParameter(1)+factorFit2*f.GetParameter(2))

        # move the stats box
        statsBoxCoord = [0.63, 0.69, 0.94, 0.80]
        root.gPad.Update()
        statsBox = h.GetListOfFunctions().FindObject('stats')
        statsBox.SetX1NDC(statsBoxCoord[0])
        statsBox.SetY1NDC(statsBoxCoord[1])
        statsBox.SetX2NDC(statsBoxCoord[2])
        statsBox.SetY2NDC(statsBoxCoord[3])
        statsBox.SetLineColor(root.kBlue)
        statsBox.SetTextColor(root.kBlue)
    return f

def plot_res(graph, name, xTitle='', yTitle='', clOpts=None):
    # create canvas and draw on it
    canvas_name = 'c_'+name
    canvas_title = canvas_name
    canvWidth = 600
    c = root.TCanvas(canvas_name, canvas_title, 100, 100, canvWidth, 600)
    c.cd()
    set_root_style()

    graph.SetMarkerStyle(root.kFullTriangleUp)
    graph.SetMarkerColor(root.kRed)
    graph.SetLineColor(root.kRed)

    # axis
    xAxis = graph.GetXaxis()
    yAxis = graph.GetYaxis()
    xAxis.SetTitle(xTitle)
    xAxis.SetTitleFont(font)
    xAxis.SetLabelFont(font)
    xAxis.SetLabelSize(fontSize)
    xAxis.SetNoExponent()
    yAxis.SetTitleOffset(1.5)
    yAxis.SetTitleFont(font)
    yAxis.SetTitle(yTitle)
    yAxis.SetLabelFont(font)
    yAxis.SetLabelSize(fontSize)
    yAxis.SetNoExponent()

    graph.Draw('ap')

    xBase = 0.5
    yBase = 0.13
    notes = extract_notes_from_name(canvas_name, xBase, yBase)
    tex = add_text(notes=notes)
    tex = add_cms_header_text(tex, clOpts)

    c.Modified()
    c.Update()

    return [c, tex]

def main():
    opts = parse_options_plotEff(parser)
    batchRun = opts.interactive
    if batchRun:
        root.gROOT.SetBatch(True)

    prefix = ''
    emul = opts.emul
    if emul:
        prefix = 'emu_'
    legacy = opts.legacy
    if legacy:
        prefix = 'legacy_'

    hm = HistManager(filename=opts.fname, subdir=opts.runnr)

    objects = []

    # histogram namebits
    reco_0to2p5 = 'probe_absEtaMin0_absEtaMax2p5_ptmin'
    reco_0to2p4 = 'probe_absEtaMin0_absEtaMax2p4_ptmin'
    reco_0to2p1 = 'probe_absEtaMin0_absEtaMax2p1_ptmin'
    reco_0to0p83 = 'probe_absEtaMin0_absEtaMax0p83_ptmin'
    reco_0p83to1p24 = 'probe_absEtaMin0p83_absEtaMax1p24_ptmin'
    reco_1p24to2p5 = 'probe_absEtaMin1p24_absEtaMax2p5_ptmin'
    reco_1p24to2p4 = 'probe_absEtaMin1p24_absEtaMax2p4_ptmin'
    reco_0to0p7 = 'probe_absEtaMin0_absEtaMax0p7_ptmin'
    reco_0p8to0p85 = 'probe_absEtaMin0p8_absEtaMax0p85_ptmin'
    reco_1p24to1p7 = 'probe_absEtaMin1p24_absEtaMax1p7_ptmin'
    reco_1p2to1p55 = 'probe_absEtaMin1p2_absEtaMax1p55_ptmin'
    reco_1p24to1p55 = 'probe_absEtaMin1p24_absEtaMax1p55_ptmin'
    reco_1p55to1p85 = 'probe_absEtaMin1p55_absEtaMax1p85_ptmin'
    reco_1p55to2p4 = 'probe_absEtaMin1p55_absEtaMax2p4_ptmin'
    reco_1p7to1p85 = 'probe_absEtaMin1p7_absEtaMax1p85_ptmin'
    reco_1p85to2p5 = 'probe_absEtaMin1p85_absEtaMax2p5_ptmin'
    reco_1p85to2p4 = 'probe_absEtaMin1p85_absEtaMax2p4_ptmin'
    #etaRanges = [reco_0to2p5, reco_0to0p83, reco_0p83to1p24, reco_1p24to2p5, reco_0to0p7, reco_0p8to0p85, reco_1p24to1p7, reco_1p7to1p85, reco_1p85to2p5]
    #etaRanges = [reco_0to2p5, reco_0to0p83, reco_0p83to1p24, reco_1p24to2p5, reco_1p24to1p7, reco_1p7to1p85, reco_1p85to2p5]
    #etaRanges = [reco_0to2p5, reco_0to0p83, reco_0p83to1p24, reco_1p24to2p5, reco_1p2to1p55, reco_1p55to1p85, reco_1p85to2p5]
    #etaRanges = [reco_0to2p4, reco_0to0p83, reco_0p83to1p24, reco_1p24to2p4, reco_1p2to1p55, reco_1p55to1p85, reco_1p85to2p4]
    #etaRanges = [reco_1p24to2p4, reco_1p24to1p55, reco_1p55to1p85, reco_1p55to2p4, reco_1p85to2p4]
    etaRanges = [reco_0to2p4, reco_0to0p83, reco_0p83to1p24, reco_1p24to2p4]
    #tfEtaRanges = [reco_0to2p5, reco_0to0p83, reco_0p83to1p24, reco_1p24to2p5]
    tfEtaRanges = [reco_0to2p4, reco_0to0p83, reco_0p83to1p24, reco_1p24to2p4]
    #tfEtaRanges = etaRanges

    #drStr = 'dr0p5_'
    #detaStr = 'deta0p5_'
    #dphiStr = 'dphi0p5_'
    ## l1 at 2nd station
    #drStr = 'dr0p3_'
    #detaStr = 'deta0p15_'
    #dphiStr = 'dphi0p25_'
    # l1 at vertex
    drStr = 'dr0p2_'
    detaStr = 'deta0p15_'
    dphiStr = 'dphi0p15_'
    ## reco at 2nd station
    #drStr = 'dr0p1_'
    #detaStr = 'deta0p1_'
    #dphiStr = 'dphi0p025_'

    res_probe_ptmins = [0.5, 20, 30, 40, 50, 60, 100, 150]

    yTitle_eff = 'L1 muon efficiency'
    yTitle_nMatch = '# best matched probes'
    xMax=100
    rebinPt = 1
    rebinEta = 1
    rebinPhi = 1
    #rebinPt = 2
    #rebinEta = 2
    #rebinPhi = 2

    if opts.eff:
        for etaRange in etaRanges:
            ## pt plots
            # quality 12
            objects.append(plot_eff_standard(hm, 'l1_muon_qualMin12_ptmin25_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, xMax=xMax, addOverflow=True, rebin=rebinPt, clOpts=opts))
            objects.append(plot_eff_standard(hm, 'l1_muon_qualMin12_ptmin25_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, addOverflow=True, clOpts=opts))

            # quality 8
            objects.append(plot_eff_standard(hm, 'l1_muon_qualMin8_ptmin7_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, xMax=xMax, addOverflow=True, rebin=rebinPt, clOpts=opts))
            objects.append(plot_eff_standard(hm, 'l1_muon_qualMin8_ptmin7_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, addOverflow=True, clOpts=opts))

            objects.append(plot_eff_standard(hm, 'l1_muon_qualMin8_ptmin15_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, xMax=xMax, addOverflow=True, rebin=rebinPt, clOpts=opts))
            objects.append(plot_eff_standard(hm, 'l1_muon_qualMin8_ptmin15_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, addOverflow=True, clOpts=opts))

            # quality 4
            objects.append(plot_eff_standard(hm, 'l1_muon_qualMin4_ptmin3_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, xMax=xMax, addOverflow=True, rebin=rebinPt, clOpts=opts))
            objects.append(plot_eff_standard(hm, 'l1_muon_qualMin4_ptmin3_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, addOverflow=True, clOpts=opts))

            ## p plots
            objects.append(plot_eff_standard(hm, 'l1_muon_qualMin12_ptmin25_'+drStr+'matched_'+etaRange+'0p5_p', etaRange+'0p5_p', 'best_', xTitle='p_{reco} (GeV/c)', yTitle=yTitle_eff, addOverflow=True, clOpts=opts))
            objects.append(plot_eff_standard(hm, 'l1_muon_qualMin8_ptmin15_'+drStr+'matched_'+etaRange+'0p5_p', etaRange+'0p5_p', 'best_', xTitle='p_{reco} (GeV/c)', yTitle=yTitle_eff, addOverflow=True, clOpts=opts))
            objects.append(plot_eff_standard(hm, 'l1_muon_qualMin4_ptmin3_'+drStr+'matched_'+etaRange+'0p5_p', etaRange+'0p5_p', 'best_', xTitle='p_{reco} (GeV/c)', yTitle=yTitle_eff, addOverflow=True, clOpts=opts))

            # phi plots
            objects.append(plot_eff_standard(hm, 'l1_muon_qualMin12_ptmin25_'+drStr+'matched_'+etaRange+'33_phi', etaRange+'33_phi', 'best_', xTitle='#phi^{reco}', yTitle=yTitle_eff, rebin=rebinPhi, clOpts=opts))

            objects.append(plot_eff_standard(hm, 'l1_muon_qualMin8_ptmin7_'+drStr+'matched_'+etaRange+'10_phi', etaRange+'10_phi', 'best_', xTitle='#phi^{reco}', yTitle=yTitle_eff, rebin=rebinPhi, clOpts=opts))
            objects.append(plot_eff_standard(hm, 'l1_muon_qualMin8_ptmin15_'+drStr+'matched_'+etaRange+'20_phi', etaRange+'20_phi', 'best_', xTitle='#phi^{reco}', yTitle=yTitle_eff, rebin=rebinPhi, clOpts=opts))

            objects.append(plot_eff_standard(hm, 'l1_muon_qualMin4_ptmin3_'+drStr+'matched_'+etaRange+'5_phi', etaRange+'5_phi', 'best_', xTitle='#phi^{reco}', yTitle=yTitle_eff, rebin=rebinPhi, clOpts=opts))

        ## eta plots
        etaRange = reco_0to2p4
        objects.append(plot_eff_standard(hm, 'l1_muon_qualMin12_ptmin25_'+drStr+'matched_'+etaRange+'33_eta', etaRange+'33_eta', 'best_', xTitle='#eta^{reco}', yTitle=yTitle_eff, rebin=rebinEta, clOpts=opts))

        objects.append(plot_eff_standard(hm, 'l1_muon_qualMin8_ptmin7_'+drStr+'matched_'+etaRange+'10_eta', etaRange+'10_eta', 'best_', xTitle='#eta^{reco}', yTitle=yTitle_eff, rebin=rebinEta, clOpts=opts))
        objects.append(plot_eff_standard(hm, 'l1_muon_qualMin8_ptmin15_'+drStr+'matched_'+etaRange+'20_eta', etaRange+'20_eta', 'best_', xTitle='#eta^{reco}', yTitle=yTitle_eff, rebin=rebinEta, clOpts=opts))

        objects.append(plot_eff_standard(hm, 'l1_muon_qualMin4_ptmin3_'+drStr+'matched_'+etaRange+'5_eta', etaRange+'5_eta', 'best_', xTitle='#eta^{reco}', yTitle=yTitle_eff, rebin=rebinEta, clOpts=opts))

        for etaRange in tfEtaRanges:
            # nVtx plots
            objects.append(plot_eff_standard(hm, 'l1_muon_qualMin12_ptmin25_'+drStr+'matched_'+etaRange+'33_vtx', etaRange+'33_vtx', 'best_', xTitle='PU', yTitle=yTitle_eff, clOpts=opts))

            if opts.runnr == 'all_runs':
                # run plots
                objects.append(plot_eff_standard(hm, 'l1_muon_qualMin12_ptmin25_'+drStr+'matched_'+etaRange+'33_run', etaRange+'33_run', 'best_', xTitle='run number', yTitle=yTitle_eff, autoZoomX=True, clOpts=opts))

            # charge plots
            objects.append(plot_eff_standard(hm, 'l1_muon_qualMin12_ptmin25_'+drStr+'matched_'+etaRange+'33_charge', etaRange+'33_charge', 'best_', xTitle='charge_{reco}', yTitle=yTitle_eff, clOpts=opts))

    # efficiency plots with the curves of all TF eta regions plus the full eta range
    if opts.efftf:
        # pt plots
        objects.append(plot_eff_eta_ranges(hm, 'l1_muon_qualMin12_ptmin25_'+drStr+'matched_probe_absEtaMinXX_absEtaMaxYY_ptmin0p5_pt', 'probe_absEtaMinXX_absEtaMaxYY_ptmin0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, xMax=xMax, addOverflow=True, rebin=rebinPt, clOpts=opts))
        objects.append(plot_eff_eta_ranges(hm, 'l1_muon_qualMin12_ptmin25_'+drStr+'matched_probe_absEtaMinXX_absEtaMaxYY_ptmin0p5_pt', 'probe_absEtaMinXX_absEtaMaxYY_ptmin0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, addOverflow=True, clOpts=opts))
        objects.append(plot_eff_eta_ranges(hm, 'l1_muon_qualMin8_ptmin7_'+drStr+'matched_probe_absEtaMinXX_absEtaMaxYY_ptmin0p5_pt', 'probe_absEtaMinXX_absEtaMaxYY_ptmin0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, xMax=xMax, addOverflow=True, rebin=rebinPt, clOpts=opts))
        objects.append(plot_eff_eta_ranges(hm, 'l1_muon_qualMin8_ptmin7_'+drStr+'matched_probe_absEtaMinXX_absEtaMaxYY_ptmin0p5_pt', 'probe_absEtaMinXX_absEtaMaxYY_ptmin0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, addOverflow=True, clOpts=opts))
        objects.append(plot_eff_eta_ranges(hm, 'l1_muon_qualMin8_ptmin15_'+drStr+'matched_probe_absEtaMinXX_absEtaMaxYY_ptmin0p5_pt', 'probe_absEtaMinXX_absEtaMaxYY_ptmin0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, xMax=xMax, addOverflow=True, rebin=rebinPt, clOpts=opts))
        objects.append(plot_eff_eta_ranges(hm, 'l1_muon_qualMin8_ptmin15_'+drStr+'matched_probe_absEtaMinXX_absEtaMaxYY_ptmin0p5_pt', 'probe_absEtaMinXX_absEtaMaxYY_ptmin0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, addOverflow=True, clOpts=opts))
        # phi plots
        objects.append(plot_eff_eta_ranges(hm, 'l1_muon_qualMin12_ptmin25_'+drStr+'matched_probe_absEtaMinXX_absEtaMaxYY_ptmin33_phi', 'probe_absEtaMinXX_absEtaMaxYY_ptmin33_phi', 'best_', xTitle='#phi^{reco}', yTitle=yTitle_eff, rebin=rebinPhi, clOpts=opts))
        objects.append(plot_eff_eta_ranges(hm, 'l1_muon_qualMin8_ptmin7_'+drStr+'matched_probe_absEtaMinXX_absEtaMaxYY_ptmin10_phi', 'probe_absEtaMinXX_absEtaMaxYY_ptmin10_phi', 'best_', xTitle='#phi^{reco}', yTitle=yTitle_eff, rebin=rebinPhi, clOpts=opts))
        objects.append(plot_eff_eta_ranges(hm, 'l1_muon_qualMin8_ptmin15_'+drStr+'matched_probe_absEtaMinXX_absEtaMaxYY_ptmin20_phi', 'probe_absEtaMinXX_absEtaMaxYY_ptmin20_phi', 'best_', xTitle='#phi^{reco}', yTitle=yTitle_eff, rebin=rebinPhi, clOpts=opts))
        # run plots
        objects.append(plot_eff_eta_ranges(hm, 'l1_muon_qualMin12_ptmin25_'+drStr+'matched_probe_absEtaMinXX_absEtaMaxYY_ptmin33_run', 'probe_absEtaMinXX_absEtaMaxYY_ptmin33_run', 'best_', xTitle='run number', yTitle=yTitle_eff, autoZoomX=True, clOpts=opts))
        objects.append(plot_eff_eta_ranges(hm, 'l1_muon_qualMin8_ptmin7_'+drStr+'matched_probe_absEtaMinXX_absEtaMaxYY_ptmin10_run', 'probe_absEtaMinXX_absEtaMaxYY_ptmin10_run', 'best_', xTitle='run number', yTitle=yTitle_eff, autoZoomX=True, clOpts=opts))
        objects.append(plot_eff_eta_ranges(hm, 'l1_muon_qualMin8_ptmin15_'+drStr+'matched_probe_absEtaMinXX_absEtaMaxYY_ptmin20_run', 'probe_absEtaMinXX_absEtaMaxYY_ptmin20_run', 'best_', xTitle='run number', yTitle=yTitle_eff, autoZoomX=True, clOpts=opts))
        # vtx plots
        objects.append(plot_eff_eta_ranges(hm, 'l1_muon_qualMin12_ptmin25_'+drStr+'matched_probe_absEtaMinXX_absEtaMaxYY_ptmin33_vtx', 'probe_absEtaMinXX_absEtaMaxYY_ptmin33_vtx', 'best_', xTitle='PU', yTitle=yTitle_eff, clOpts=opts))
        objects.append(plot_eff_eta_ranges(hm, 'l1_muon_qualMin8_ptmin7_'+drStr+'matched_probe_absEtaMinXX_absEtaMaxYY_ptmin10_vtx', 'probe_absEtaMinXX_absEtaMaxYY_ptmin10_vtx', 'best_', xTitle='PU', yTitle=yTitle_eff, clOpts=opts))
        objects.append(plot_eff_eta_ranges(hm, 'l1_muon_qualMin8_ptmin15_'+drStr+'matched_probe_absEtaMinXX_absEtaMaxYY_ptmin20_vtx', 'probe_absEtaMinXX_absEtaMaxYY_ptmin20_vtx', 'best_', xTitle='PU', yTitle=yTitle_eff, clOpts=opts))

    if opts.qualcomp:
        for etaRange in etaRanges:
            ## pt plots
            objects.append(plot_eff_qual(hm, 'l1_muon_qualMinXX_ptmin25_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, addOverflow=True, clOpts=opts))

        ## eta plots
        etaRange = reco_0to2p4
        objects.append(plot_eff_qual(hm, 'l1_muon_qualMinXX_ptmin25_'+drStr+'matched_'+etaRange+'33_eta', etaRange+'33_eta', 'best_', xTitle='#eta^{reco}', yTitle=yTitle_eff, rebin=rebinEta, clOpts=opts))

    if opts.delta:
        etaRange = reco_0to2p4
        # delta R plots
        objects.append(plot_hists_data_emul(hm, 'l1_muon_qualMin12_ptmin25_'+drStr+'matched_'+etaRange+'0p5_dr', hNamePrefix='best_', normToNEntries=100, xTitle='#DeltaR', yTitle=yTitle_nMatch, clOpts=opts))

        # delta eta plots
        objects.append(plot_hists_data_emul(hm, 'l1_muon_qualMin12_ptmin25_'+detaStr+'matched_'+etaRange+'0p5_deta', hNamePrefix='best_', normToNEntries=100, xTitle='#Delta#eta', yTitle=yTitle_nMatch, clOpts=opts))

        # delta phi plots
        objects.append(plot_hists_data_emul(hm, 'l1_muon_qualMin12_ptmin25_'+dphiStr+'matched_'+etaRange+'0p5_dphi', hNamePrefix='best_', normToNEntries=100, xTitle='#Delta#phi', yTitle=yTitle_nMatch, clOpts=opts))

        resGraphs = []
        # L1 - reco plots
        for etaRange in tfEtaRanges:
            plotNames = [etaRange+'0p5_'+drStr+'matched_l1_muon_qualMin12_ptmin0p5']

            for plotName in plotNames:
                # inverse pt resolution
                objects.append(plot_hists_data_emul(hm, plotName+'_dinvpt', hNamePrefix=prefix+'res_best_', normToNEntries=100, xTitle='(p_{T}^{RECO} - p_{T}^{L1}) / p_{T}^{L1}', clOpts=opts))
                if opts.fitdelta and len(objects[-1]) > 1:
                    lastHs = objects[-1][1]
                    fitHisto(lastHs[0], 'gaus')
                # pt resolution
                objects.append(plot_hists_data_emul(hm, plotName+'_dpt', hNamePrefix=prefix+'res_best_', normToNEntries=100, xTitle='p_{T}^{L1} - p_{T}^{RECO}', clOpts=opts))
                if opts.fitdelta and len(objects[-1]) > 1:
                    lastHs = objects[-1][1]
                    fitHisto(lastHs[0], 'crystalball')
                # eta resolution
                objects.append(plot_hists_data_emul(hm, plotName+'_deta', hNamePrefix=prefix+'res_best_', normToNEntries=100, xTitle='#eta_{L1} - #eta_{RECO}', clOpts=opts))
                if opts.fitdelta and len(objects[-1]) > 1:
                    lastHs = objects[-1][1]
                    fitHisto(lastHs[0], 'gaus')
                # phi resolution
                objects.append(plot_hists_data_emul(hm, plotName+'_dphi', hNamePrefix=prefix+'res_best_', normToNEntries=100, xTitle='#phi_{L1} - #phi_{RECO}', clOpts=opts))
                if opts.fitdphi and len(objects[-1]) > 1:
                    lastHs = objects[-1][1]
                    fitHisto(lastHs[0], 'gaus')
                # charge matching
                objects.append(plot_hists_data_emul(hm, plotName+'_dcharge', hNamePrefix=prefix+'res_best_', normToNEntries=100, xTitle='charge^{L1} - charge^{RECO}', clOpts=opts))

            if opts.fitdelta:
                resGraphs.append(root.TGraphErrors(len(res_probe_ptmins)))
                resGraphs.append(root.TGraphErrors(len(res_probe_ptmins)))
                resGraphs.append(root.TGraphErrors(len(res_probe_ptmins)))
                resGraphs.append(root.TGraphErrors(len(res_probe_ptmins)))
                if opts.fitdphi:
                    resGraphs.append(root.TGraphErrors(len(res_probe_ptmins)))
                    resGraphs.append(root.TGraphErrors(len(res_probe_ptmins)))
                # for resolution plots by pT range
                for j, probe_pt_min in enumerate(res_probe_ptmins):
                    probe_ptmin_str = str(probe_pt_min).replace('.', 'p')
                    if j < len(res_probe_ptmins)-1:
                        probe_pt_max = res_probe_ptmins[j+1]
                        probe_ptmax_str = '_ptmax'+str(probe_pt_max).replace('.', 'p')
                    else:
                        probe_pt_max = 2*res_probe_ptmins[-1]
                        probe_ptmax_str = '_ptmax'

                    plotName = etaRange+probe_ptmin_str+probe_ptmax_str+'_'+drStr+'matched_l1_muon_qualMin12'

                    # inverse pt resolution
                    objects.append(plot_hists_data_emul(hm, plotName+'_dinvpt', hNamePrefix=prefix+'res_best_', normToNEntries=100, xTitle='(p_{T}^{RECO} - p_{T}^{L1}) / p_{T}^{L1}', clOpts=opts))
                    if len(objects[-1]) > 1:
                        lastHs = objects[-1][1]
                        #fit = fitHisto(lastHs[0], 'crystalball')
                        fit = fitHisto(lastHs[0], 'gaus')
                        resGraphs[-1].SetPoint(j, (probe_pt_max+probe_pt_min)/2., fit.GetParameter(1))
                        resGraphs[-1].SetPointError(j, (probe_pt_max-probe_pt_min)/2., fit.GetParError(1))
                        resGraphs[-2].SetPoint(j, (probe_pt_max+probe_pt_min)/2., fit.GetParameter(2))
                        resGraphs[-2].SetPointError(j, (probe_pt_max-probe_pt_min)/2., fit.GetParError(2))
                    else:
                        resGraphs[-1].RemovePoint(j)
                        resGraphs[-2].RemovePoint(j)
                    # eta resolution
                    objects.append(plot_hists_data_emul(hm, plotName+'_deta', hNamePrefix=prefix+'res_best_', normToNEntries=100, xTitle='#eta_{L1} - #eta_{RECO}', clOpts=opts))
                    if len(objects[-1]) > 1:
                        lastHs = objects[-1][1]
                        fit = fitHisto(lastHs[0], 'gaus')
                        resGraphs[-3].SetPoint(j, (probe_pt_max+probe_pt_min)/2., fit.GetParameter(1))
                        resGraphs[-3].SetPointError(j, (probe_pt_max-probe_pt_min)/2., fit.GetParError(1))
                        resGraphs[-4].SetPoint(j, (probe_pt_max+probe_pt_min)/2., fit.GetParameter(2))
                        resGraphs[-4].SetPointError(j, (probe_pt_max-probe_pt_min)/2., fit.GetParError(2))
                    else:
                        resGraphs[-3].RemovePoint(j)
                        resGraphs[-4].RemovePoint(j)
                    if opts.fitdphi:
                        # phi resolution
                        objects.append(plot_hists_data_emul(hm, plotName+'_dphi', hNamePrefix=prefix+'res_best_', normToNEntries=100, xTitle='#phi_{L1} - #phi_{RECO}', clOpts=opts))
                        if len(objects[-1]) > 1:
                            lastHs = objects[-1][1]
                            fit = fitHisto(lastHs[0], 'gaus')
                            resGraphs[-5].SetPoint(j, (probe_pt_max+probe_pt_min)/2., fit.GetParameter(1))
                            resGraphs[-5].SetPointError(j, (probe_pt_max-probe_pt_min)/2., fit.GetParError(1))
                            resGraphs[-6].SetPoint(j, (probe_pt_max+probe_pt_min)/2., fit.GetParameter(2))
                            resGraphs[-6].SetPointError(j, (probe_pt_max-probe_pt_min)/2., fit.GetParError(2))
                        else:
                            resGraphs[-5].RemovePoint(j)
                            resGraphs[-6].RemovePoint(j)
                objects.append(plot_res(resGraphs[-1], 'res_'+etaRange+'_dinvpt_mean', xTitle='p_{T}^{RECO}', yTitle='<(p_{T}^{RECO} - p_{T}^{L1}) / p_{T}^{L1}>', clOpts=opts))
                objects.append(plot_res(resGraphs[-2], 'res_'+etaRange+'_dinvpt_sigma', xTitle='p_{T}^{RECO}', yTitle='#sigma((p_{T}^{RECO} - p_{T}^{L1}) / p_{T}^{L1})', clOpts=opts))
                objects.append(plot_res(resGraphs[-3], 'res_'+etaRange+'_deta_mean', xTitle='p_{T}^{RECO}', yTitle='<#eta_{L1} - #eta_{RECO}>', clOpts=opts))
                objects.append(plot_res(resGraphs[-4], 'res_'+etaRange+'_deta_sigma', xTitle='p_{T}^{RECO}', yTitle='#sigma(#eta_{L1} - #eta_{RECO})', clOpts=opts))
                if opts.fitdphi:
                    objects.append(plot_res(resGraphs[-5], 'res_'+etaRange+'_dphi_mean', xTitle='p_{T}^{RECO}', yTitle='<#phi_{L1} - #phi_{RECO}>', clOpts=opts))
                    objects.append(plot_res(resGraphs[-6], 'res_'+etaRange+'_dphi_sigma', xTitle='p_{T}^{RECO}', yTitle='#sigma(#phi_{L1} - #phi_{RECO})', clOpts=opts))

    # 2d reco vs. L1 plots
    if opts.twod:
        hm2d = HistManager2d(filename=opts.fname, subdir=opts.runnr)
        etaRanges2d = [reco_0to0p83, reco_0p83to1p24, reco_1p24to2p4]
        etaRange = reco_0to2p4
        # quality 8
        objects.append(plot_2dhist(hm2d, prefix+'2d_best_'+etaRange+'0p5_'+drStr+'matched_l1_muon_qualMin8_ptmin0p5_pt', clOpts=opts))
        objects.append(plot_2dhist(hm2d, prefix+'2d_best_'+etaRange+'10_'+dphiStr+'matched_l1_muon_qualMin8_ptmin7_eta', clOpts=opts))
        objects.append(plot_2dhist(hm2d, prefix+'2d_best_'+etaRange+'10_'+detaStr+'matched_l1_muon_qualMin8_ptmin7_phi', clOpts=opts))

        for etaRange in etaRanges2d:
            objects.append(plot_2dhist(hm2d, prefix+'2d_best_'+etaRange+'0p5_'+drStr+'matched_l1_muon_qualMin8_ptmin0p5_pt', clOpts=opts))
            objects.append(plot_2dhist(hm2d, prefix+'2d_best_'+etaRange+'10_'+detaStr+'matched_l1_muon_qualMin8_ptmin7_phi', clOpts=opts))

        # quality 12
        etaRange = reco_0to2p4
        objects.append(plot_2dhist(hm2d, prefix+'2d_best_'+etaRange+'0p5_'+drStr+'matched_l1_muon_qualMin12_ptmin0p5_pt', clOpts=opts))
        objects.append(plot_2dhist(hm2d, prefix+'2d_best_'+etaRange+'33_'+dphiStr+'matched_l1_muon_qualMin12_ptmin25_eta', clOpts=opts))
        objects.append(plot_2dhist(hm2d, prefix+'2d_best_'+etaRange+'33_'+detaStr+'matched_l1_muon_qualMin12_ptmin25_phi', clOpts=opts))

        for etaRange in etaRanges2d:
            objects.append(plot_2dhist(hm2d, prefix+'2d_best_'+etaRange+'0p5_'+drStr+'matched_l1_muon_qualMin12_ptmin0p5_pt', clOpts=opts))
            objects.append(plot_2dhist(hm2d, prefix+'2d_best_'+etaRange+'33_'+detaStr+'matched_l1_muon_qualMin12_ptmin25_phi', clOpts=opts))

    if opts.dataemul:
        for etaRange in etaRanges:
            ## pt plots
            # quality 12
            objects.append(plot_eff_data_emul(hm, 'l1_muon_qualMin12_ptmin0p5_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, xMax=xMax, addOverflow=True, rebin=rebinPt, clOpts=opts))
            objects.append(plot_eff_data_emul(hm, 'l1_muon_qualMin12_ptmin0p5_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, addOverflow=True, clOpts=opts))

            objects.append(plot_eff_data_emul(hm, 'l1_muon_qualMin12_ptmin25_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, xMax=xMax, addOverflow=True, rebin=rebinPt, clOpts=opts))
            objects.append(plot_eff_data_emul(hm, 'l1_muon_qualMin12_ptmin25_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, addOverflow=True, clOpts=opts))

            # quality 8
            objects.append(plot_eff_data_emul(hm, 'l1_muon_qualMin8_ptmin0p5_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, xMax=xMax, addOverflow=True, rebin=rebinPt, clOpts=opts))
            objects.append(plot_eff_data_emul(hm, 'l1_muon_qualMin8_ptmin0p5_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, addOverflow=True, clOpts=opts))

            objects.append(plot_eff_data_emul(hm, 'l1_muon_qualMin8_ptmin7_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, xMax=xMax, addOverflow=True, rebin=rebinPt, clOpts=opts))
            objects.append(plot_eff_data_emul(hm, 'l1_muon_qualMin8_ptmin7_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, addOverflow=True, clOpts=opts))

            objects.append(plot_eff_data_emul(hm, 'l1_muon_qualMin8_ptmin15_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, xMax=xMax, addOverflow=True, rebin=rebinPt, clOpts=opts))
            objects.append(plot_eff_data_emul(hm, 'l1_muon_qualMin8_ptmin15_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, addOverflow=True, clOpts=opts))

            # quality 4
            objects.append(plot_eff_data_emul(hm, 'l1_muon_qualMin4_ptmin0p5_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, xMax=xMax, addOverflow=True, rebin=rebinPt, clOpts=opts))
            objects.append(plot_eff_data_emul(hm, 'l1_muon_qualMin4_ptmin0p5_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, addOverflow=True, clOpts=opts))

            objects.append(plot_eff_data_emul(hm, 'l1_muon_qualMin4_ptmin3_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, xMax=xMax, addOverflow=True, rebin=rebinPt, clOpts=opts))
            objects.append(plot_eff_data_emul(hm, 'l1_muon_qualMin4_ptmin3_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, addOverflow=True, clOpts=opts))

            ## p plots
            objects.append(plot_eff_data_emul(hm, 'l1_muon_qualMin12_ptmin25_'+drStr+'matched_'+etaRange+'0p5_p', etaRange+'0p5_p', 'best_', xTitle='p_{reco} (GeV/c)', yTitle=yTitle_eff, addOverflow=True, clOpts=opts))

            # phi plots
            objects.append(plot_eff_data_emul(hm, 'l1_muon_qualMin12_ptmin25_'+drStr+'matched_'+etaRange+'33_phi', etaRange+'33_phi', 'best_', xTitle='#phi^{reco}', yTitle=yTitle_eff, rebin=rebinPhi, clOpts=opts))

            objects.append(plot_eff_data_emul(hm, 'l1_muon_qualMin8_ptmin7_'+drStr+'matched_'+etaRange+'10_phi', etaRange+'10_phi', 'best_', xTitle='#phi^{reco}', yTitle=yTitle_eff, rebin=rebinPhi, clOpts=opts))
            objects.append(plot_eff_data_emul(hm, 'l1_muon_qualMin8_ptmin15_'+drStr+'matched_'+etaRange+'20_phi', etaRange+'20_phi', 'best_', xTitle='#phi^{reco}', yTitle=yTitle_eff, rebin=rebinPhi, clOpts=opts))

        ## eta plots
        etaRange = reco_0to2p4
        objects.append(plot_eff_data_emul(hm, 'l1_muon_qualMin12_ptmin25_'+drStr+'matched_'+etaRange+'33_eta', etaRange+'33_eta', 'best_', xTitle='#eta^{reco}', yTitle=yTitle_eff, rebin=rebinEta, clOpts=opts))

        objects.append(plot_eff_data_emul(hm, 'l1_muon_qualMin8_ptmin7_'+drStr+'matched_'+etaRange+'10_eta', etaRange+'10_eta', 'best_', xTitle='#eta^{reco}', yTitle=yTitle_eff, rebin=rebinEta, clOpts=opts))
        objects.append(plot_eff_data_emul(hm, 'l1_muon_qualMin8_ptmin15_'+drStr+'matched_'+etaRange+'20_eta', etaRange+'20_eta', 'best_', xTitle='#eta^{reco}', yTitle=yTitle_eff, rebin=rebinEta, clOpts=opts))

        objects.append(plot_eff_data_emul(hm, 'l1_muon_qualMin4_ptmin3_'+drStr+'matched_'+etaRange+'5_eta', etaRange+'5_eta', 'best_', xTitle='#eta^{reco}', yTitle=yTitle_eff, rebin=rebinEta, clOpts=opts))

        for etaRange in tfEtaRanges:
            # nVtx plots
            objects.append(plot_eff_data_emul(hm, 'l1_muon_qualMin12_ptmin25_'+drStr+'matched_'+etaRange+'33_vtx', etaRange+'33_vtx', 'best_', xTitle='PU', yTitle=yTitle_eff, clOpts=opts))

            objects.append(plot_eff_data_emul(hm, 'l1_muon_qualMin8_ptmin7_'+drStr+'matched_'+etaRange+'10_vtx', etaRange+'10_vtx', 'best_', xTitle='PU', yTitle=yTitle_eff, clOpts=opts))
            objects.append(plot_eff_data_emul(hm, 'l1_muon_qualMin8_ptmin15_'+drStr+'matched_'+etaRange+'20_vtx', etaRange+'20_vtx', 'best_', xTitle='PU', yTitle=yTitle_eff, clOpts=opts))

            objects.append(plot_eff_data_emul(hm, 'l1_muon_qualMin4_ptmin3_'+drStr+'matched_'+etaRange+'5_vtx', etaRange+'5_vtx', 'best_', xTitle='PU', yTitle=yTitle_eff, clOpts=opts))

            if opts.runnr == 'all_runs':
                # run plots
                objects.append(plot_eff_data_emul(hm, 'l1_muon_qualMin12_ptmin25_'+drStr+'matched_'+etaRange+'33_run', etaRange+'33_run', 'best_', xTitle='run number', yTitle=yTitle_eff, autoZoomX=True, clOpts=opts))

                objects.append(plot_eff_data_emul(hm, 'l1_muon_qualMin8_ptmin7_'+drStr+'matched_'+etaRange+'10_run', etaRange+'10_run', 'best_', xTitle='run number', yTitle=yTitle_eff, autoZoomX=True, clOpts=opts))
                objects.append(plot_eff_data_emul(hm, 'l1_muon_qualMin8_ptmin15_'+drStr+'matched_'+etaRange+'20_run', etaRange+'20_run', 'best_', xTitle='run number', yTitle=yTitle_eff, autoZoomX=True, clOpts=opts))

                objects.append(plot_eff_data_emul(hm, 'l1_muon_qualMin4_ptmin3_'+drStr+'matched_'+etaRange+'5_run', etaRange+'5_run', 'best_', xTitle='run number', yTitle=yTitle_eff, autoZoomX=True, clOpts=opts))

            # charge plots
            objects.append(plot_eff_data_emul(hm, 'l1_muon_qualMin12_ptmin25_'+drStr+'matched_'+etaRange+'33_charge', etaRange+'33_charge', 'best_', xTitle='charge_{reco}', yTitle=yTitle_eff, clOpts=opts))

            objects.append(plot_eff_data_emul(hm, 'l1_muon_qualMin8_ptmin7_'+drStr+'matched_'+etaRange+'10_charge', etaRange+'10_charge', 'best_', xTitle='charge_{reco}', yTitle=yTitle_eff, clOpts=opts))
            objects.append(plot_eff_data_emul(hm, 'l1_muon_qualMin8_ptmin15_'+drStr+'matched_'+etaRange+'20_charge', etaRange+'20_charge', 'best_', xTitle='charge_{reco}', yTitle=yTitle_eff, clOpts=opts))

            objects.append(plot_eff_data_emul(hm, 'l1_muon_qualMin4_ptmin3_'+drStr+'matched_'+etaRange+'5_charge', etaRange+'5_charge', 'best_', xTitle='charge_{reco}', yTitle=yTitle_eff, clOpts=opts))

    if opts.upgradelegacy:
        for etaRange in etaRanges:
            ## pt plots
            # quality 12
            objects.append(plot_eff_upgrade_legacy(hm, 'l1_muon_qualMin12_ptmin0p5_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, xMax=xMax, addOverflow=True, rebin=rebinPt, clOpts=opts))
            objects.append(plot_eff_upgrade_legacy(hm, 'l1_muon_qualMin12_ptmin0p5_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, addOverflow=True, clOpts=opts))

            objects.append(plot_eff_upgrade_legacy(hm, 'l1_muon_qualMin12_ptmin25_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, xMax=xMax, addOverflow=True, rebin=rebinPt, clOpts=opts))
            objects.append(plot_eff_upgrade_legacy(hm, 'l1_muon_qualMin12_ptmin25_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, addOverflow=True, clOpts=opts))

            # quality 8
            objects.append(plot_eff_upgrade_legacy(hm, 'l1_muon_qualMin8_ptmin0p5_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, xMax=xMax, addOverflow=True, rebin=rebinPt, clOpts=opts))
            objects.append(plot_eff_upgrade_legacy(hm, 'l1_muon_qualMin8_ptmin0p5_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, addOverflow=True, clOpts=opts))

            objects.append(plot_eff_upgrade_legacy(hm, 'l1_muon_qualMin8_ptmin7_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, xMax=xMax, addOverflow=True, rebin=rebinPt, clOpts=opts))
            objects.append(plot_eff_upgrade_legacy(hm, 'l1_muon_qualMin8_ptmin7_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, addOverflow=True, clOpts=opts))

            objects.append(plot_eff_upgrade_legacy(hm, 'l1_muon_qualMin8_ptmin15_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, xMax=xMax, addOverflow=True, rebin=rebinPt, clOpts=opts))
            objects.append(plot_eff_upgrade_legacy(hm, 'l1_muon_qualMin8_ptmin15_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, addOverflow=True, clOpts=opts))

            # quality 4
            objects.append(plot_eff_upgrade_legacy(hm, 'l1_muon_qualMin4_ptmin0p5_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, xMax=xMax, addOverflow=True, rebin=rebinPt, clOpts=opts))
            objects.append(plot_eff_upgrade_legacy(hm, 'l1_muon_qualMin4_ptmin0p5_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, addOverflow=True, clOpts=opts))

            objects.append(plot_eff_upgrade_legacy(hm, 'l1_muon_qualMin4_ptmin3_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, xMax=xMax, addOverflow=True, rebin=rebinPt, clOpts=opts))
            objects.append(plot_eff_upgrade_legacy(hm, 'l1_muon_qualMin4_ptmin3_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, addOverflow=True, clOpts=opts))

            ## p plots
            objects.append(plot_eff_upgrade_legacy(hm, 'l1_muon_qualMin12_ptmin25_'+drStr+'matched_'+etaRange+'0p5_p', etaRange+'0p5_p', 'best_', xTitle='p_{reco} (GeV/c)', yTitle=yTitle_eff, addOverflow=True, clOpts=opts))

            # phi plots
            objects.append(plot_eff_upgrade_legacy(hm, 'l1_muon_qualMin12_ptmin25_'+drStr+'matched_'+etaRange+'33_phi', etaRange+'33_phi', 'best_', xTitle='#phi^{reco}', yTitle=yTitle_eff, rebin=rebinPhi, clOpts=opts))

            objects.append(plot_eff_upgrade_legacy(hm, 'l1_muon_qualMin8_ptmin7_'+drStr+'matched_'+etaRange+'10_phi', etaRange+'10_phi', 'best_', xTitle='#phi^{reco}', yTitle=yTitle_eff, rebin=rebinPhi, clOpts=opts))
            objects.append(plot_eff_upgrade_legacy(hm, 'l1_muon_qualMin8_ptmin15_'+drStr+'matched_'+etaRange+'20_phi', etaRange+'20_phi', 'best_', xTitle='#phi^{reco}', yTitle=yTitle_eff, rebin=rebinPhi, clOpts=opts))

        ## eta plots
        etaRange = reco_0to2p4
        objects.append(plot_eff_upgrade_legacy(hm, 'l1_muon_qualMin12_ptmin25_'+drStr+'matched_'+etaRange+'33_eta', etaRange+'33_eta', 'best_', xTitle='#eta^{reco}', yTitle=yTitle_eff, rebin=rebinEta, clOpts=opts))

        objects.append(plot_eff_upgrade_legacy(hm, 'l1_muon_qualMin8_ptmin7_'+drStr+'matched_'+etaRange+'10_eta', etaRange+'10_eta', 'best_', xTitle='#eta^{reco}', yTitle=yTitle_eff, rebin=rebinEta, clOpts=opts))
        objects.append(plot_eff_upgrade_legacy(hm, 'l1_muon_qualMin8_ptmin15_'+drStr+'matched_'+etaRange+'20_eta', etaRange+'20_eta', 'best_', xTitle='#eta^{reco}', yTitle=yTitle_eff, rebin=rebinEta, clOpts=opts))

        objects.append(plot_eff_upgrade_legacy(hm, 'l1_muon_qualMin4_ptmin3_'+drStr+'matched_'+etaRange+'5_eta', etaRange+'5_eta', 'best_', xTitle='#eta^{reco}', yTitle=yTitle_eff, rebin=rebinEta, clOpts=opts))

        for etaRange in tfEtaRanges:
            # nVtx plots
            objects.append(plot_eff_upgrade_legacy(hm, 'l1_muon_qualMin12_ptmin25_'+drStr+'matched_'+etaRange+'33_vtx', etaRange+'33_vtx', 'best_', xTitle='PU', yTitle=yTitle_eff, clOpts=opts))

            objects.append(plot_eff_upgrade_legacy(hm, 'l1_muon_qualMin8_ptmin7_'+drStr+'matched_'+etaRange+'10_vtx', etaRange+'10_vtx', 'best_', xTitle='PU', yTitle=yTitle_eff, clOpts=opts))
            objects.append(plot_eff_upgrade_legacy(hm, 'l1_muon_qualMin8_ptmin15_'+drStr+'matched_'+etaRange+'20_vtx', etaRange+'20_vtx', 'best_', xTitle='PU', yTitle=yTitle_eff, clOpts=opts))

            objects.append(plot_eff_upgrade_legacy(hm, 'l1_muon_qualMin4_ptmin3_'+drStr+'matched_'+etaRange+'5_vtx', etaRange+'5_vtx', 'best_', xTitle='PU', yTitle=yTitle_eff, clOpts=opts))

            if opts.runnr == 'all_runs':
                # run plots
                objects.append(plot_eff_upgrade_legacy(hm, 'l1_muon_qualMin12_ptmin25_'+drStr+'matched_'+etaRange+'33_run', etaRange+'33_run', 'best_', xTitle='run number', yTitle=yTitle_eff, autoZoomX=True, clOpts=opts))

                objects.append(plot_eff_upgrade_legacy(hm, 'l1_muon_qualMin8_ptmin7_'+drStr+'matched_'+etaRange+'10_run', etaRange+'10_run', 'best_', xTitle='run number', yTitle=yTitle_eff, autoZoomX=True, clOpts=opts))
                objects.append(plot_eff_upgrade_legacy(hm, 'l1_muon_qualMin8_ptmin15_'+drStr+'matched_'+etaRange+'20_run', etaRange+'20_run', 'best_', xTitle='run number', yTitle=yTitle_eff, autoZoomX=True, clOpts=opts))

                objects.append(plot_eff_upgrade_legacy(hm, 'l1_muon_qualMin4_ptmin3_'+drStr+'matched_'+etaRange+'5_run', etaRange+'5_run', 'best_', xTitle='run number', yTitle=yTitle_eff, autoZoomX=True, clOpts=opts))

            # charge plots
            objects.append(plot_eff_upgrade_legacy(hm, 'l1_muon_qualMin12_ptmin25_'+drStr+'matched_'+etaRange+'33_charge', etaRange+'33_charge', 'best_', xTitle='charge_{reco}', yTitle=yTitle_eff, clOpts=opts))

            objects.append(plot_eff_upgrade_legacy(hm, 'l1_muon_qualMin8_ptmin7_'+drStr+'matched_'+etaRange+'10_charge', etaRange+'10_charge', 'best_', xTitle='charge_{reco}', yTitle=yTitle_eff, clOpts=opts))
            objects.append(plot_eff_upgrade_legacy(hm, 'l1_muon_qualMin8_ptmin15_'+drStr+'matched_'+etaRange+'20_charge', etaRange+'20_charge', 'best_', xTitle='charge_{reco}', yTitle=yTitle_eff, clOpts=opts))

            objects.append(plot_eff_upgrade_legacy(hm, 'l1_muon_qualMin4_ptmin3_'+drStr+'matched_'+etaRange+'5_charge', etaRange+'5_charge', 'best_', xTitle='charge_{reco}', yTitle=yTitle_eff, clOpts=opts))

    if opts.bytf:
        for etaRange in etaRanges:
            # quality 12
            objects.append(plot_eff_tf(hm, 'l1_muon_qualMin12_ptmin25_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, xMax=xMax, addOverflow=True, rebin=rebinPt, clOpts=opts))
            objects.append(plot_eff_tf(hm, 'l1_muon_qualMin12_ptmin25_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, addOverflow=True, clOpts=opts))
            # quality 8
            objects.append(plot_eff_tf(hm, 'l1_muon_qualMin8_ptmin7_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, xMax=xMax, addOverflow=True, rebin=rebinPt, clOpts=opts))
            objects.append(plot_eff_tf(hm, 'l1_muon_qualMin8_ptmin7_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, addOverflow=True, clOpts=opts))

            objects.append(plot_eff_tf(hm, 'l1_muon_qualMin8_ptmin15_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, xMax=xMax, addOverflow=True, rebin=rebinPt, clOpts=opts))
            objects.append(plot_eff_tf(hm, 'l1_muon_qualMin8_ptmin15_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, addOverflow=True, clOpts=opts))

            # quality 4
            objects.append(plot_eff_tf(hm, 'l1_muon_qualMin4_ptmin3_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, xMax=xMax, addOverflow=True, rebin=rebinPt, clOpts=opts))
            objects.append(plot_eff_tf(hm, 'l1_muon_qualMin4_ptmin3_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, addOverflow=True, clOpts=opts))

            ## p plots
            objects.append(plot_eff_tf(hm, 'l1_muon_qualMin12_ptmin25_'+drStr+'matched_'+etaRange+'0p5_p', etaRange+'0p5_p', 'best_', xTitle='p_{reco} (GeV/c)', yTitle=yTitle_eff, addOverflow=True, clOpts=opts))

            # phi plots
            objects.append(plot_eff_tf(hm, 'l1_muon_qualMin12_ptmin25_'+drStr+'matched_'+etaRange+'33_phi', etaRange+'33_phi', 'best_', xTitle='#phi^{reco}', yTitle=yTitle_eff, rebin=rebinPhi, clOpts=opts))
            objects.append(plot_eff_tf(hm, 'l1_muon_qualMin8_ptmin7_'+drStr+'matched_'+etaRange+'10_phi', etaRange+'10_phi', 'best_', xTitle='#phi^{reco}', yTitle=yTitle_eff, rebin=rebinPhi, clOpts=opts))
            objects.append(plot_eff_tf(hm, 'l1_muon_qualMin8_ptmin15_'+drStr+'matched_'+etaRange+'20_phi', etaRange+'20_phi', 'best_', xTitle='#phi^{reco}', yTitle=yTitle_eff, rebin=rebinPhi, clOpts=opts))

        ## eta plots
        etaRange = reco_0to2p4
        objects.append(plot_eff_tf(hm, 'l1_muon_qualMin12_ptmin25_'+drStr+'matched_'+etaRange+'33_eta', etaRange+'33_eta', 'best_', xTitle='#eta^{reco}', yTitle=yTitle_eff, rebin=rebinEta, clOpts=opts))

        objects.append(plot_eff_tf(hm, 'l1_muon_qualMin8_ptmin7_'+drStr+'matched_'+etaRange+'10_eta', etaRange+'10_eta', 'best_', xTitle='#eta^{reco}', yTitle=yTitle_eff, rebin=rebinEta, clOpts=opts))
        objects.append(plot_eff_tf(hm, 'l1_muon_qualMin8_ptmin15_'+drStr+'matched_'+etaRange+'20_eta', etaRange+'20_eta', 'best_', xTitle='#eta^{reco}', yTitle=yTitle_eff, rebin=rebinEta, clOpts=opts))

        objects.append(plot_eff_tf(hm, 'l1_muon_qualMin4_ptmin3_'+drStr+'matched_'+etaRange+'5_eta', etaRange+'5_eta', 'best_', xTitle='#eta^{reco}', yTitle=yTitle_eff, rebin=rebinEta, clOpts=opts))

        for etaRange in tfEtaRanges:
            # charge plots
            objects.append(plot_eff_tf(hm, 'l1_muon_qualMin12_ptmin25_'+drStr+'matched_'+etaRange+'33_charge', etaRange+'33_charge', 'best_', xTitle='charge_{reco}', yTitle=yTitle_eff, clOpts=opts))

    if opts.bycharge:
        for etaRange in etaRanges:
            # quality 12
            objects.append(plot_eff_charge(hm, 'l1_muon_qualMin12_ptmin25_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, xMax=xMax, addOverflow=True, rebin=rebinPt, clOpts=opts))
            objects.append(plot_eff_charge(hm, 'l1_muon_qualMin12_ptmin25_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, addOverflow=True, clOpts=opts))

            # quality 8
            objects.append(plot_eff_charge(hm, 'l1_muon_qualMin8_ptmin7_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, xMax=xMax, addOverflow=True, rebin=rebinPt, clOpts=opts))
            objects.append(plot_eff_charge(hm, 'l1_muon_qualMin8_ptmin7_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, addOverflow=True, clOpts=opts))

            objects.append(plot_eff_charge(hm, 'l1_muon_qualMin8_ptmin15_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, xMax=xMax, addOverflow=True, rebin=rebinPt, clOpts=opts))
            objects.append(plot_eff_charge(hm, 'l1_muon_qualMin8_ptmin15_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, addOverflow=True, clOpts=opts))

            # quality 4
            objects.append(plot_eff_charge(hm, 'l1_muon_qualMin4_ptmin3_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, xMax=xMax, addOverflow=True, rebin=rebinPt, clOpts=opts))
            objects.append(plot_eff_charge(hm, 'l1_muon_qualMin4_ptmin3_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, addOverflow=True, clOpts=opts))

            ## p plots
            objects.append(plot_eff_charge(hm, 'l1_muon_qualMin12_ptmin25_'+drStr+'matched_'+etaRange+'0p5_p', etaRange+'0p5_p', 'best_', xTitle='p_{reco} (GeV/c)', yTitle=yTitle_eff, addOverflow=True, clOpts=opts))

            # phi plots
            objects.append(plot_eff_charge(hm, 'l1_muon_qualMin12_ptmin25_'+drStr+'matched_'+etaRange+'33_phi', etaRange+'33_phi', 'best_', xTitle='#phi^{reco}', yTitle=yTitle_eff, rebin=rebinPhi, clOpts=opts))

            objects.append(plot_eff_charge(hm, 'l1_muon_qualMin8_ptmin7_'+drStr+'matched_'+etaRange+'10_phi', etaRange+'10_phi', 'best_', xTitle='#phi^{reco}', yTitle=yTitle_eff, rebin=rebinPhi, clOpts=opts))
            objects.append(plot_eff_charge(hm, 'l1_muon_qualMin8_ptmin15_'+drStr+'matched_'+etaRange+'20_phi', etaRange+'20_phi', 'best_', xTitle='#phi^{reco}', yTitle=yTitle_eff, rebin=rebinPhi, clOpts=opts))

        ## eta plots
        etaRange = reco_0to2p4
        objects.append(plot_eff_charge(hm, 'l1_muon_qualMin12_ptmin25_'+drStr+'matched_'+etaRange+'33_eta', etaRange+'33_eta', 'best_', xTitle='#eta^{reco}', yTitle=yTitle_eff, rebin=rebinEta, clOpts=opts))

        objects.append(plot_eff_charge(hm, 'l1_muon_qualMin8_ptmin7_'+drStr+'matched_'+etaRange+'10_eta', etaRange+'10_eta', 'best_', xTitle='#eta^{reco}', yTitle=yTitle_eff, rebin=rebinEta, clOpts=opts))
        objects.append(plot_eff_charge(hm, 'l1_muon_qualMin8_ptmin15_'+drStr+'matched_'+etaRange+'20_eta', etaRange+'20_eta', 'best_', xTitle='#eta^{reco}', yTitle=yTitle_eff, rebin=rebinEta, clOpts=opts))

        objects.append(plot_eff_charge(hm, 'l1_muon_qualMin4_ptmin3_'+drStr+'matched_'+etaRange+'5_eta', etaRange+'5_eta', 'best_', xTitle='#eta^{reco}', yTitle=yTitle_eff, rebin=rebinEta, clOpts=opts))

        for etaRange in tfEtaRanges:
            # charge plots
            objects.append(plot_eff_charge(hm, 'l1_muon_qualMin12_ptmin25_'+drStr+'matched_'+etaRange+'33_charge', etaRange+'33_charge', 'best_', xTitle='charge_{reco}', yTitle=yTitle_eff, clOpts=opts))

    if opts.public:
        yTitle_eff = 'L1 muon efficiency'
        for etaRange in etaRanges:
            # pt plots
            objects.append(plot_eff_public(hm, 'l1_muon_qualMin12_ptmin25_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T} (GeV)', yTitle=yTitle_eff, xMax=500, addOverflow=True, clOpts=opts))
            objects.append(plot_eff_public(hm, 'l1_muon_qualMin12_ptmin25_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T} (GeV)', yTitle=yTitle_eff, xMax=50, addOverflow=True, clOpts=opts))

            objects.append(plot_eff_public_pt(hm, 'l1_muon_qualMin8_ptminXX_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', ptmins=[7, 15], xTitle='p_{T} (GeV)', yTitle=yTitle_eff, xMax=500, addOverflow=True, clOpts=opts))
            objects.append(plot_eff_public_pt(hm, 'l1_muon_qualMin8_ptminXX_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', ptmins=[7, 15], xTitle='p_{T} (GeV)', yTitle=yTitle_eff, xMax=50, addOverflow=True, clOpts=opts))

            # phi plots
            objects.append(plot_eff_public(hm, 'l1_muon_qualMin12_ptmin25_'+drStr+'matched_'+etaRange+'33_phi', etaRange+'33_phi', 'best_', xTitle='#phi', yTitle=yTitle_eff, rebin=7, clOpts=opts))

            # run number
            objects.append(plot_eff_public(hm, 'l1_muon_qualMin12_ptmin25_'+drStr+'matched_'+etaRange+'33_run', etaRange+'33_run', 'best_', xTitle='run', yTitle=yTitle_eff, autoZoomX=True, clOpts=opts))

        etaRange = reco_0to2p4
        # eta plots
        objects.append(plot_eff_public(hm, 'l1_muon_qualMin12_ptmin25_'+drStr+'matched_'+etaRange+'33_eta', etaRange+'33_eta', 'best_', xTitle='#eta', yTitle=yTitle_eff, rebin=4, clOpts=opts))
        objects.append(plot_eff_public(hm, 'l1_muon_qualMin8_ptmin7_'+drStr+'matched_'+etaRange+'10_eta', etaRange+'10_eta', 'best_', xTitle='#eta', yTitle=yTitle_eff, rebin=4, clOpts=opts))
        objects.append(plot_eff_public(hm, 'l1_muon_qualMin8_ptmin15_'+drStr+'matched_'+etaRange+'20_eta', etaRange+'20_eta', 'best_', xTitle='#eta', yTitle=yTitle_eff, rebin=4, clOpts=opts))

    if opts.fname2:
        if opts.runnr2:
            runnr2 = opts.runnr2
        else:
            runnr2 = 'all_runs'
        hm2 = HistManager(filename=opts.fname2, subdir=runnr2)
        legTxt1 = opts.legtxt1
        legTxt2 = opts.legtxt2
        
        for etaRange in etaRanges:
            ## pt plots
            # quality 12
            objects.append(plot_eff_hm_comp(hm, hm2, 'l1_muon_qualMin12_ptmin25_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, xMax=xMax, addOverflow=True, rebin=rebinPt, clOpts=opts))
            objects.append(plot_eff_hm_comp(hm, hm2, 'l1_muon_qualMin12_ptmin25_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, addOverflow=True, clOpts=opts))

            # quality 8
            objects.append(plot_eff_hm_comp(hm, hm2, 'l1_muon_qualMin8_ptmin7_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, xMax=xMax, addOverflow=True, rebin=rebinPt, clOpts=opts))
            objects.append(plot_eff_hm_comp(hm, hm2, 'l1_muon_qualMin8_ptmin7_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, addOverflow=True, clOpts=opts))

            objects.append(plot_eff_hm_comp(hm, hm2, 'l1_muon_qualMin8_ptmin15_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, xMax=xMax, addOverflow=True, rebin=rebinPt, clOpts=opts))
            objects.append(plot_eff_hm_comp(hm, hm2, 'l1_muon_qualMin8_ptmin15_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, addOverflow=True, clOpts=opts))

            # quality 4
            objects.append(plot_eff_hm_comp(hm, hm2, 'l1_muon_qualMin4_ptmin3_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, xMax=xMax, addOverflow=True, rebin=rebinPt, clOpts=opts))
            objects.append(plot_eff_hm_comp(hm, hm2, 'l1_muon_qualMin4_ptmin3_'+drStr+'matched_'+etaRange+'0p5_pt', etaRange+'0p5_pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_eff, addOverflow=True, clOpts=opts))

            # phi plots
            objects.append(plot_eff_hm_comp(hm, hm2, 'l1_muon_qualMin12_ptmin25_'+drStr+'matched_'+etaRange+'33_phi', etaRange+'33_phi', 'best_', xTitle='#phi^{reco}', yTitle=yTitle_eff, rebin=rebinPhi, clOpts=opts))

            objects.append(plot_eff_hm_comp(hm, hm2, 'l1_muon_qualMin8_ptmin7_'+drStr+'matched_'+etaRange+'10_phi', etaRange+'10_phi', 'best_', xTitle='#phi^{reco}', yTitle=yTitle_eff, rebin=rebinPhi, clOpts=opts))
            objects.append(plot_eff_hm_comp(hm, hm2, 'l1_muon_qualMin8_ptmin15_'+drStr+'matched_'+etaRange+'20_phi', etaRange+'20_phi', 'best_', xTitle='#phi^{reco}', yTitle=yTitle_eff, rebin=rebinPhi, clOpts=opts))

            objects.append(plot_eff_hm_comp(hm, hm2, 'l1_muon_qualMin4_ptmin3_'+drStr+'matched_'+etaRange+'5_phi', etaRange+'5_phi', 'best_', xTitle='#phi^{reco}', yTitle=yTitle_eff, rebin=rebinPhi, clOpts=opts))

        ## eta plots
        etaRange = reco_0to2p4
        objects.append(plot_eff_hm_comp(hm, hm2, 'l1_muon_qualMin12_ptmin25_'+drStr+'matched_'+etaRange+'33_eta', etaRange+'33_eta', 'best_', xTitle='#eta^{reco}', yTitle=yTitle_eff, rebin=rebinEta, clOpts=opts))

        objects.append(plot_eff_hm_comp(hm, hm2, 'l1_muon_qualMin8_ptmin7_'+drStr+'matched_'+etaRange+'10_eta', etaRange+'10_eta', 'best_', xTitle='#eta^{reco}', yTitle=yTitle_eff, rebin=rebinEta, clOpts=opts))
        objects.append(plot_eff_hm_comp(hm, hm2, 'l1_muon_qualMin8_ptmin15_'+drStr+'matched_'+etaRange+'20_eta', etaRange+'20_eta', 'best_', xTitle='#eta^{reco}', yTitle=yTitle_eff, rebin=rebinEta, clOpts=opts))

        objects.append(plot_eff_hm_comp(hm, hm2, 'l1_muon_qualMin4_ptmin3_'+drStr+'matched_'+etaRange+'5_eta', etaRange+'5_eta', 'best_', xTitle='#eta^{reco}', yTitle=yTitle_eff, rebin=rebinEta, clOpts=opts))

        for etaRange in tfEtaRanges:
            # nVtx plots
            objects.append(plot_eff_hm_comp(hm, hm2, 'l1_muon_qualMin12_ptmin25_'+drStr+'matched_'+etaRange+'33_vtx', etaRange+'33_vtx', 'best_', xTitle='PU', yTitle=yTitle_eff, clOpts=opts))

            #if opts.runnr == 'all_runs':
            #    # run plots
            #    objects.append(plot_eff_hm_comp(hm, hm2, 'l1_muon_qualMin12_ptmin25_'+drStr+'matched_'+etaRange+'33_run', etaRange+'33_run', 'best_', xTitle='run number', yTitle=yTitle_eff, autoZoomX=True, clOpts=opts))

            # charge plots
            objects.append(plot_eff_hm_comp(hm, hm2, 'l1_muon_qualMin12_ptmin25_'+drStr+'matched_'+etaRange+'33_charge', etaRange+'33_charge', 'best_', xTitle='charge_{reco}', yTitle=yTitle_eff, clOpts=opts))

    if opts.control:
        eta_range_0to2p4 = '_absEtaMin0_absEtaMax2p4'
        eta_range_0to0p83 = '_absEtaMin0_absEtaMax0p83'
        eta_range_0p83to1p24 = '_absEtaMin0p83_absEtaMax1p24'
        eta_range_1p24to2p4 = '_absEtaMin1p24_absEtaMax2p4'
        yTitle_tags = '# tag muons'
        yTitle_probes = '# probe muons'
        yTitle_l1 = '# L1 muons'
        rebin_eta = 2
        rebin_phi = 2
        eta_ranges = [eta_range_0to2p4, eta_range_0to0p83, eta_range_0p83to1p24, eta_range_1p24to2p4]
        qualMins = [4, 8, 12]
        for eta_range in eta_ranges:
            objects.append(plot_hists_standard(hm, 'probe'+eta_range+'_ptmin0p5_pt', xTitle='p_{T}^{reco} (GeV)', yTitle=yTitle_probes+' / GeV', normToBinWidth=True, logY=True, addOverflow=True, clOpts=opts))
            objects.append(plot_hists_standard(hm, 'probe'+eta_range+'_ptmin0p5_p', xTitle='p^{reco} (GeV)', yTitle=yTitle_probes+' / GeV', normToBinWidth=True, logY=True, addOverflow=True, clOpts=opts))
            objects.append(plot_hists_standard(hm, 'probe'+eta_range+'_ptmin0p5_phi', xTitle='#phi^{reco}', yTitle=yTitle_probes, rebin=rebin_phi, clOpts=opts))
            objects.append(plot_hists_standard(hm, 'probe'+eta_range+'_ptmin0p5_charge', xTitle='#mu charge', yTitle=yTitle_probes, clOpts=opts))

            for qualMin in qualMins:
                objects.append(plot_hists_standard(hm, 'l1_muon'+eta_range+'_qualMin{q}_ptmin0p5_pt'.format(q=qualMin), xTitle='p_{T}^{L1} (GeV)', yTitle=yTitle_l1+' / GeV', normToBinWidth=True, xMax=300, addOverflow=True, logY=True, clOpts=opts))
                objects.append(plot_hists_standard(hm, 'l1_muon'+eta_range+'_qualMin{q}_ptmin0p5_phi'.format(q=qualMin), xTitle='#phi^{L1}', yTitle=yTitle_l1, rebin=rebin_phi, clOpts=opts))
                objects.append(plot_hists_standard(hm, 'l1_muon'+eta_range+'_qualMin{q}_ptmin0p5_charge'.format(q=qualMin), xTitle='L1 #mu charge', yTitle=yTitle_l1, clOpts=opts))

        objects.append(plot_hists_standard(hm, 'probe'+eta_range_0to2p4+'_ptmin0p5_eta', xTitle='#eta^{reco}', yTitle=yTitle_probes, rebin=rebin_eta, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'probe'+eta_range_0to2p4+'_ptmin0p5_vtx', xTitle='# vertex', yTitle=yTitle_probes, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'n_probes'+eta_range_0to2p4+'_ptmin0p5', xTitle='probe muon multiplicity', yTitle='# events', logY=True, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'l1_muon'+eta_range_0to2p4+'_qualMin4_ptmin0p5_eta', xTitle='#eta^{L1}', yTitle=yTitle_l1, rebin=rebin_eta, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'l1_muon'+eta_range_0to2p4+'_qualMin8_ptmin0p5_eta', xTitle='#eta^{L1}', yTitle=yTitle_l1, rebin=rebin_eta, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'l1_muon'+eta_range_0to2p4+'_qualMin12_ptmin0p5_eta', xTitle='#eta^{L1}', yTitle=yTitle_l1, rebin=rebin_eta, clOpts=opts))

        # tag plots
        objects.append(plot_hists_standard(hm, 'tag_pt', xTitle='tag p_{T}^{reco} (GeV)', yTitle=yTitle_tags+' / GeV', normToBinWidth=True, logY=True, addOverflow=True, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'tag_phi', xTitle='tag #phi^{reco}', yTitle=yTitle_tags, rebin=rebin_phi, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'tag_eta', xTitle='tag #eta^{reco}', yTitle=yTitle_tags, rebin=rebin_eta, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'tag_charge', xTitle='tag #mu charge', yTitle=yTitle_tags, clOpts=opts))

    # save canvases to root file
    if savePlots:
        plotdir = 'plots_'+opts.fname.replace('.root','').partition('/')[0] + '_' + opts.runnr
        #plotdir += '_oldSelVsNewSel'
        if opts.bycharge:
            plotdir += '_chargeComp'
        if opts.bytf:
            plotdir += '_tf'
        if opts.public:
            plotdir += '_public'
        if opts.fname2:
            plotdir += '_'+legTxt1+'_vs_'+legTxt2
        if not os.path.exists(plotdir):
            os.makedirs(plotdir)
        output = root.TFile('./'+plotdir+'/ugmt_tandp_eff_plots.root', 'recreate')
        output.cd()
        for obj in objects:
            c = obj[0]
            plotName = c.GetName()[2:]
            c.Write(plotName)
            if opts.public:
                c.Print('./'+plotdir+'/'+plotName+'.pdf', '.pdf')
            c.Print('./'+plotdir+'/'+plotName+'.png', '.png')
        print 'get the plots with: scp -r lxplus:{pwd}/{plotdir}/ .'.format(pwd=os.getcwd(), plotdir=plotdir)
        output.Close()

    # wait
    if not batchRun:
        raw_input("Press ENTER to quit.")

if __name__ == "__main__":
    savePlots = True
    batchRun = True
    best_only = False
    font = 42
    fontSize = 0.04
    main()

