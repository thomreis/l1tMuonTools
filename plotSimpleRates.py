#!/usr/bin/env python
from ToolBox import parse_options_and_init_log
# have to do this first or ROOT masks the -h messages
opts, parser = parse_options_and_init_log()

from L1Analysis import L1Ana, L1Ntuple
from analysis_tools.plotting import HistManager
from analysis_tools.selections import MuonSelections, Matcher
import ROOT as root
import re
import os

def parse_options_plotRates(parser):
    """
    Adds often used options to the OptionParser...
    """
    parsers = parser.add_subparsers()
    sub_parser = parsers.add_parser("plotRates")
    sub_parser.add_argument("-i", "--interactive", dest="interactive", action='store_false', help="Draw plots on screen.")
    sub_parser.add_argument("-b", "--bunches", dest="bunches", default=0, type=int, help="Number of colliding bunches")
    sub_parser.add_argument("--pu", dest="pu", default=20, type=int, help="Average PU. default=20")
    sub_parser.add_argument("--xsect", dest="xsect", default=80, type=float, help="Total cross section in mb. default=80 mb")
    sub_parser.add_argument("--instlumi", dest="instlumi", default=1.2e34, type=float, help="Instantaneous luminosity. default=1.2e-34 cm-2s-1")
    sub_parser.add_argument("--scale", dest="scale", default=1., type=float, help="Additional scale factor for rate calculation")
    sub_parser.add_argument("--scale-to-bunches", dest="scaleto", default=0, type=int, help="Scale to this many bunches")
    sub_parser.add_argument("-l", "--legacy", dest="legacy", action='store_true', help="Draw plots relative to legacy.")
    sub_parser.add_argument("--qstack", dest="qstack", default=False, action='store_true', help="Plot stacked rate by quality.")
    sub_parser.add_argument("--public", dest="public", default=False, action='store_true', help="Plot style for publication.")
    sub_parser.add_argument("--print-rates", dest="printrates", default=False, action='store_true', help="Print out table of rates.")
    # options to compare histograms from two files/runs
    sub_parser.add_argument("--fname2", dest="fname2", default=None, help="Second file to take reference histograms from.")
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
    root.gStyle.SetOptFit(0)
    root.gStyle.SetOptTitle(0)
    root.gPad.SetTopMargin(0.08)
    root.gPad.SetLeftMargin(0.14)
    root.gPad.SetRightMargin(0.06)
    root.gPad.SetTickx(1)
    root.gPad.SetTicky(1)
    root.gPad.SetGridx(1)
    root.gPad.SetGridy(1)

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
    if addOverflow:
        tex.SetTextAngle(90)
        tex.DrawLatex(0.87, 0.20, 'highest bin includes overflow')
        tex.SetTextAngle(0)
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
 
def extract_notes_from_name(name, xBase, yBase, etaTxt=True, qualTxt=True, ptTxt=True, public=False):
    notes = []
    # extract eta range
    if etaTxt and name.rfind('_eta') == -1:
        eta_number_strs = re.findall(r'[\d\.\d]+', name[name.find('EtaMin')+6:name.find('EtaMax')+12])
        if len(eta_number_strs) > 1:
            if public:
                if len(eta_number_strs[0]) > 3:
                    eta_number_strs[0] = eta_number_strs[0][0:3]
                if len(eta_number_strs[1]) > 3:
                    eta_number_strs[1] = eta_number_strs[1][0:3]
                note_str = eta_number_strs[0]+' < |#eta| < '+eta_number_strs[1]
            else:
                note_str = eta_number_strs[0]+' < |#eta| < '+eta_number_strs[1]
            notes.append([xBase, yBase+0.15, note_str, True])
    # extract quality
    if qualTxt and name.find('qmin') != -1:
        qualPos = name.find('qmin')
        qual_strs = re.findall(r'[\d\.\d]+', name[qualPos+4:qualPos+6])
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
                qual_note_str = 'Quality #geq '+qual_strs[0]
            notes.append([xBase, yBase+0.1, qual_note_str, True])
    # extract pt range
    if ptTxt:
        ptminPos = name.find('ptmin')
        l1_ptmin_strs = re.findall(r'\d+\.?\d*', name[ptminPos+5:ptminPos+9])
        if len(l1_ptmin_strs) > 0:
            if l1_ptmin_strs[0][-1] == '.':
                l1_ptmin_strs[0] = l1_ptmin_strs[0][0:-1]
            if public:
                notes.append([xBase, yBase+0.05, 'L1 p_{T} #geq '+l1_ptmin_strs[0]+' GeV', True])
            else:
                notes.append([xBase, yBase+0.05, 'p_{T}^{L1} #geq '+l1_ptmin_strs[0]+' GeV', True])
    return notes

def add_main_pad(canvas, hDefs):

    return canvas

def plot_hists2(hDefs, axisTitles=['', ''], threshold=False, normToBinWidth=False, canvasPrefix='', notes=None, scaleFactor=1., data=False, errorBars=True, clOpts=None):
    den = hDefs[0]['den']
    if den:
        name = canvasPrefix+hDefs[0]['num']+'_over_'+den
    else:
        name = canvasPrefix+hDefs[0]['num']
    if normToBinWidth and not threshold and not den:
        name = 'normToBinWidth_'+name

    # setup legend according to how many histograms are in the plot
    legYmin = 0.9-0.04*len(hDefs)
    legXmin = 0.68
    legXmax = 0.9
    canvWidth = 600
    if legYmin < 0.6:
        legXmin = 0.8
        legXmax = 1.
        canvWidth = 730
    legend = setup_legend(xMin=legXmin, yMin=legYmin, xMax=legXmax, yMax=0.9)
    legEntries = []

    hs = []
    hStack = root.THStack()
    # get all the histograms and set their plot style
    for hDef in hDefs:
        hm = hDef['hm']
        hName = hDef['num']
        # denominator
        if hDef['denhm']:
            denHm = hDef['denhm']
        else:
            denHm = hm
        hDenName = hDef['den']

        # get the histogram as a threshold histogram or normal histogram, with or without ratio
        if threshold:
            h = hm.get_threshold_hist(hName).Clone()
            if hDenName:
                hDen = denHm.get_threshold_hist(hDenName)
                h.Divide(h, hDen, 1, 1, "b")
        else:
            h = hm.get(hName).Clone()
            if hDenName:
                hDen = denHm.get(hDenName)
                h.Divide(h, hDen, 1, 1, "b")

        if normToBinWidth and not threshold and not hDenName:
            for bin in range(1, h.GetNbinsX()+1):
               h.SetBinContent(bin, h.GetBinContent(bin) / h.GetBinWidth(bin))
               h.SetBinError(bin, h.GetBinError(bin) / h.GetBinWidth(bin))
        elif normToBinWidth:
            print 'Ignoring normToBinWidth flag for threshold or ratio plots'

        if scaleFactor != 1.:
            h.Scale(scaleFactor)
        h.SetLineColor(hDef['lc'])
        h.SetLineStyle(hDef['ls'])
        if errorBars:
            h.SetLineWidth(2)
            legStyle = 'l'
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
        canvas_name = 'c_rates_stacked_'+name
        stackHistos = hStack.GetStack()
        j = len(stackHistos)-1
        for i, hDef in enumerate(hDefs):
            if hDef['fc']:
                hs[i] = stackHistos[j].Clone()
                j -= 1
    else:
        canvas_name = 'c_rates_'+name

    if scaleFactor != 1.:
        canvas_name += '_scaled'

    # create canvas and draw on it
    canvas_title = canvas_name
    c = root.TCanvas(canvas_name, canvas_title, 100, 100, canvWidth, 600)
    c.cd()
    if name[-2:] == 'pt' and not den:
        c.SetLogy(True)

    set_root_style()
    if legYmin < 0.6:
        root.gPad.SetRightMargin(0.2)

    xAxis = hs[0].GetXaxis()
    yAxis = hs[0].GetYaxis()
    maxBinValue = hs[0].GetBinContent(hs[0].GetMaximumBin())
    xAxis.SetTitle(axisTitles[0])
    xAxis.SetTitleFont(font)
    xAxis.SetLabelFont(font)
    xAxis.SetLabelSize(fontSize)
    xAxis.SetNoExponent()
    yAxis.SetTitleOffset(1.5)
    yAxis.SetTitle(axisTitles[1])
    yAxis.SetTitleFont(font)
    yAxis.SetLabelFont(font)
    yAxis.SetLabelSize(fontSize)
    yMax = yAxis.GetXmax()
    if not c.GetLogy():
        yMax = 1.2*maxBinValue
        if maxBinValue <= 1. and scaleFactor == 1.:
            yMax = 1.3
        yAxis.SetRangeUser(0., yMax)
    # draw
    if errorBars:
        hs[0].SetLineWidth(2)
    legEntries[0].SetObject(hs[0])
    if errorBars:
        legEntries[0].SetOption(legEntries[0].GetOption()+'le')
    hs[0].Draw('hist')
    for h in hs[1:]:
        h.Draw('histsame')
    if errorBars:
        hs[0].Draw('same')
    hs[0].Draw('sameaxis')

    # draw vertical lines to mark TF boundaries
    lines = None
    if not clOpts.public:
        lines = draw_tf_eta_regions(hName=name, yMax=yMax)

    legend.Draw('same')

    tex = add_text(notes, True, False)
    tex = add_cms_header_text(tex, clOpts, rPos=0.94)

    c.Modified()
    c.Update()

    return [c, hs, legend, lines, tex]

def plot_hists(hDefs, xTitle=None, yTitle='# muons', threshold=False, normToBinWidth=False, canvasPrefix='', notes=None, scaleFactor=1., data=False, errorBars=True, clOpts=None):
    den = hDefs[0]['den']
    if den:
        name = canvasPrefix+hDefs[0]['num']+'_over_'+den
    else:
        name = canvasPrefix+hDefs[0]['num']
    if normToBinWidth and not threshold and not den:
        name = 'normToBinWidth_'+name

    # setup legend according to how many histograms are in the plot
    legYmin = 0.9-0.04*len(hDefs)
    legXmin = 0.68
    legXmax = 0.9
    canvWidth = 600
    if legYmin < 0.6:
        legXmin = 0.8
        legXmax = 1.
        canvWidth = 730
    legend = setup_legend(xMin=legXmin, yMin=legYmin, xMax=legXmax, yMax=0.9)
    legEntries = []

    hs = []
    hStack = root.THStack()
    # get all the histograms and set their plot style
    for hDef in hDefs:
        hm = hDef['hm']
        hName = hDef['num']
        #if hName[:4] == 'gmt_':
        #    hName = hName.replace('qmin12', 'qmin5')
        #    hName = hName.replace('qmin8', 'qmin3')
        #    hName = hName.replace('qmin4', 'qmin2')

        if threshold:
            h = hm.get_threshold_hist(hName).Clone()
            if den:
                hDen = hm.get_threshold_hist(den)
                h.Divide(h, hDen, 1, 1, "b")
        else:
            if den:
                h = hm.get_ratio(hName, den).Clone()
            else:
                h = hm.get(hName).Clone()

        if normToBinWidth and not threshold and not den:
            for bin in range(1, h.GetNbinsX()+1):
               h.SetBinContent(bin, h.GetBinContent(bin) / h.GetBinWidth(bin))
               h.SetBinError(bin, h.GetBinError(bin) / h.GetBinWidth(bin))
        elif normToBinWidth:
            print 'Ignoring normToBinWidth flag for threshold or ratio plots'

        if scaleFactor != 1.:
            h.Scale(scaleFactor)
        h.SetLineColor(hDef['lc'])
        h.SetLineStyle(hDef['ls'])
        if errorBars:
            h.SetLineWidth(2)
            legStyle = 'l'
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
        canvas_name = 'c_rates_stacked_'+name
        stackHistos = hStack.GetStack()
        j = len(stackHistos)-1
        for i, hDef in enumerate(hDefs):
            if hDef['fc']:
                hs[i] = stackHistos[j].Clone()
                j -= 1
    else:
        canvas_name = 'c_rates_'+name

    if scaleFactor != 1.:
        canvas_name += '_scaled'

    # create canvas and draw on it
    canvas_title = canvas_name
    c = root.TCanvas(canvas_name, canvas_title, 100, 100, canvWidth, 600)
    c.cd()
    if name[-2:] == 'pt' and not den:
        c.SetLogy(True)

    set_root_style()
    if legYmin < 0.6:
        root.gPad.SetRightMargin(0.2)

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
    maxBinValue = hs[0].GetBinContent(hs[0].GetMaximumBin())
    yMax = yAxis.GetXmax()
    if not c.GetLogy():
        yMax = 1.2*maxBinValue
        if maxBinValue <= 1. and scaleFactor == 1.:
            yMax = 1.3
        yAxis.SetRangeUser(0., yMax)
    # draw
    if errorBars:
        hs[0].SetLineWidth(2)
    legEntries[0].SetObject(hs[0])
    if errorBars:
        legEntries[0].SetOption(legEntries[0].GetOption()+'le')
    hs[0].Draw('hist')
    for h in hs[1:]:
        h.Draw('histsame')
    if errorBars:
        hs[0].Draw('same')
    hs[0].Draw('sameaxis')

    # draw vertical lines to mark TF boundaries
    lines = None
    if not clOpts.public:
        lines = draw_tf_eta_regions(hName=name, yMax=yMax)

    legend.Draw('same')

    tex = add_text(notes, True, False)
    tex = add_cms_header_text(tex, clOpts, rPos=0.94)

    c.Modified()
    c.Update()

    return [c, hs, legend, lines, tex]

def print_rates(hm, hName, scaleFactor=1., includeLegacy=False):
    hNames = ['ugmt_'+hName, 'bmtf_ugmt_'+hName, 'omtf_ugmt_'+hName, 'emtf_ugmt_'+hName]
    titles = ['uGMT', 'BMTF', 'OMTF', 'EMTF']
    if includeLegacy:
        hNames.insert(0, 'gmt_'+hName.replace('qmin12', 'qmin5'))
        titles.insert(0, 'GMT')

    print '===== Rates ====='
    print hName
    print ''
    histos = {}
    for name in hNames:
        histos[name] = hm.get_threshold_hist(name).Clone()
        if scaleFactor != 1.:
            histos[name].Scale(scaleFactor)

    titleStr = 'Threshold'
    for i, title in enumerate(titles):
        titleStr += '          '
        if i > 0:
            titleStr += '                 '
        titleStr += title
    print titleStr
    for threshold in [0, 3, 5, 7, 10, 12, 14, 16, 18, 20, 22, 25, 30, 40, 50, 60, 70, 80, 90, 100]:
        printStr = '{threshold:>3} GeV:'.format(threshold=threshold)
        conts = []
        errs = []
        for name in hNames:
            binNr = histos[name].FindBin(threshold)
            conts.append(histos[name].GetBinContent(binNr))
            errs.append(histos[name].GetBinError(binNr))
            printStr += ' | {cont:>8.3f} +/- {err:>7.3f} kHz'.format(cont=conts[-1], err=errs[-1])
            if len(conts) > 1:
                ratio = -1.
                if conts[0] != 0:
                    ratio = conts[-1]/conts[0]
                printStr += ' {ratio:>6.3f}'.format(ratio=ratio)
        print printStr
    print '================='

def print_rates_comp(hm, hm2, hName, scaleFactor=1., legTxts=['', '']):
    hNames = ['ugmt_'+hName, 'emtf_ugmt_'+hName]
    titles = ['uGMT', 'EMTF']

    print '===== Rates ====='
    print hName
    print ''
    histos = {}
    histos2 = {}
    for name in hNames:
        histos[name] = hm.get_threshold_hist(name).Clone()
        histos2[name] = hm2.get_threshold_hist(name).Clone()
        if scaleFactor != 1.:
            histos[name].Scale(scaleFactor)
            histos2[name].Scale(scaleFactor)

    titleStr = 'Threshold'
    for i, title in enumerate(titles):
        titleStr += '      '
        if i > 0:
            titleStr += '                 '
        titleStr += '{title} {ltxt1}                    {title} {ltxt2}'.format(title=title, ltxt1=legTxts[0], ltxt2=legTxts[1])
    print titleStr
    for threshold in [0, 3, 5, 7, 10, 12, 14, 16, 18, 20, 22, 25, 30, 40, 50, 60, 70, 80, 90, 100]:
        printStr = '{threshold:>3} GeV:'.format(threshold=threshold)
        conts = []
        errs = []
        conts2 = []
        errs2 = []
        for name in hNames:
            binNr = histos[name].FindBin(threshold)
            conts.append(histos[name].GetBinContent(binNr))
            errs.append(histos[name].GetBinError(binNr))
            binNr2 = histos2[name].FindBin(threshold)
            conts2.append(histos2[name].GetBinContent(binNr2))
            errs2.append(histos2[name].GetBinError(binNr2))
            printStr += ' | {cont:>8.3f} +/- {err:>7.3f} kHz'.format(cont=conts[-1], err=errs[-1])
            printStr += ' | {cont:>8.3f} +/- {err:>7.3f} kHz'.format(cont=conts2[-1], err=errs2[-1])
            ratio = -1.
            if conts[-1] != 0:
                ratio = conts2[-1]/conts[-1]
            printStr += ' {ratio:>6.3f}'.format(ratio=ratio)
        print printStr
    print '================='

def hist_styles(stacked=False):
    styles = {}
    styles['gmt'] = {'lc':root.kCyan, 'ls':root.kSolid, 'fc':None, 'legtext':'GMT'}
    styles['ugmt'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':None, 'legtext':'uGMT'}
    styles['hm1'] = {'lc':root.kBlue-4, 'ls':root.kSolid, 'fc':None, 'legtext':'hm1'}
    styles['hm2'] = {'lc':root.kRed-4, 'ls':root.kSolid, 'fc':None, 'legtext':'hm2'}
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

def plot_hists_standard(hm, hName, den=None, xTitle='', yTitle='# muons', threshold=False, stacked=False, normToBinWidth=False, tfMuonOrig='ugmt', reg='', scaleFactor=1., data=False, clOpts=None):
    styles = hist_styles(stacked)

    if tfMuonOrig == 'ugmt':
        ugmt_str = '_ugmt'
        prefix = ''
    elif tfMuonOrig == 'tf':
        ugmt_str = ''
        prefix = 'tf_'

    ugmt_dict = {'hm':hm, 'num':'ugmt_'+hName, 'den':den}
    bmtf_dict = {'hm':hm, 'num':'bmtf'+ugmt_str+'_'+hName, 'den':den}
    omtf_dict = {'hm':hm, 'num':'omtf'+ugmt_str+'_'+hName, 'den':den}
    emtf_dict = {'hm':hm, 'num':'emtf'+ugmt_str+'_'+hName, 'den':den}
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
    if clOpts.legacy:
        if den:
            gmt_dict = {'num':den, 'den':den}
            gmt_dict.update(styles['gmt'])
        else:
            hName_gmt = hName#.replace('qmin12', 'qmin5')
            #hName_gmt = hName_gmt.replace('qmin8', 'qmin3')
            #hName_gmt = hName_gmt.replace('qmin4', 'qmin2')
            gmt_dict = {'num':'gmt_'+hName_gmt, 'den':den}
            gmt_dict.update(styles['gmt'])
        hDefs.append(gmt_dict)

    xBase = 0.17
    yBase = 0.56
    notes = extract_notes_from_name(hName, xBase, yBase, public=clOpts.public)
    if den:
        den_eta_number_strs = re.findall(r'[\d\.\d]+', den[den.find('EtaMin')+6:den.find('EtaMax')+12])
        #if len(den_eta_number_strs) > 1 and eta_number_strs != den_eta_number_strs:
        if len(den_eta_number_strs) > 1:
            den_note_str = den_eta_number_strs[0]+' < |#eta^{GMT}| < '+den_eta_number_strs[1]
            notes.append([xBase, yBase-0.05, den_note_str, True])

    return plot_hists(hDefs, xTitle, yTitle, threshold, normToBinWidth, prefix, notes, scaleFactor, data, clOpts=clOpts)

def plot_hists_qstack(hm, hName, den=None, xTitle='', yTitle='# muons', threshold=False, stacked=False, normToBinWidth=False, tfMuonOrig='ugmt', reg='', scaleFactor=1., data=False, clOpts=None):
    styles = hist_styles(False)

    if tfMuonOrig == 'ugmt':
        ugmt_str = '_ugmt'
        style_str = 'ugmt'
        prefix = 'q_'
    elif tfMuonOrig == 'tf':
        ugmt_str = ''
        style_str = 'tf'
        prefix = 'q_tf_'

    hDefs = []
    if reg == '':
        #ugmt_dict = {'hm':hm, 'num':'ugmt_'+hName, 'den':den}
        #ugmt_dict.update(styles['ugmt'])
        #hDefs.append(ugmt_dict)
        for q in reversed(range(16)):
            ugmt_q_dict = {'hm':hm, 'num':'ugmt_'+hName.replace('qmin12', 'q{q}'.format(q=q)), 'den':den}
            ugmt_q_dict.update(styles['ugmt_q{q}'.format(q=q)])
            hDefs.append(ugmt_q_dict)
    elif reg == 'b':
        #bmtf_dict = {'hm':hm, 'num':'bmtf'+ugmt_str+'_'+hName, 'den':den}
        #bmtf_dict.update(styles['bmtf'+ugmt_str+'_q'])
        #hDefs.append(bmtf_dict)
        for q in reversed(range(16)):
            bmtf_q_dict = {'hm':hm, 'num':'bmtf'+ugmt_str+'_'+hName.replace('qmin12', 'q{q}'.format(q=q)), 'den':den}
            bmtf_q_dict.update(styles[style_str+'_q{q}'.format(q=q)])
            hDefs.append(bmtf_q_dict)
        prefix += 'bmtf_'
    elif reg == 'o':
        #omtf_dict = {'hm':hm, 'num':'omtf'+ugmt_str+'_'+hName, 'den':den}
        #omtf_dict.update(styles['omtf'+ugmt_str+'_q'])
        #hDefs.append(omtf_dict)
        for q in reversed(range(16)):
            omtf_q_dict = {'hm':hm, 'num':'omtf'+ugmt_str+'_'+hName.replace('qmin12', 'q{q}'.format(q=q)), 'den':den}
            omtf_q_dict.update(styles[style_str+'_q{q}'.format(q=q)])
            hDefs.append(omtf_q_dict)
        prefix += 'omtf_'
    elif reg == 'e':
        #emtf_dict = {'hm':hm, 'num':'emtf'+ugmt_str+'_'+hName, 'den':den}
        #emtf_dict.update(styles['emtf'+ugmt_str+'_q'])
        #hDefs.append(emtf_dict)
        for q in reversed(range(16)):
            emtf_q_dict = {'hm':hm, 'num':'emtf'+ugmt_str+'_'+hName.replace('qmin12', 'q{q}'.format(q=q)), 'den':den}
            emtf_q_dict.update(styles[style_str+'_q{q}'.format(q=q)])
            hDefs.append(emtf_q_dict)
        prefix += 'emtf_'
    if clOpts.legacy:
        if den:
            gmt_dict = {'hm':hm, 'num':den, 'den':den}
            gmt_dict.update(styles['gmt'])
        else:
            gmt_dict = {'hm':hm, 'num':'gmt_'+hName, 'den':den}
            gmt_dict.update(styles['gmt'])
        hDefs.append(gmt_dict)

    xBase = 0.07
    yBase = 0.56
    notes = extract_notes_from_name(hName, xBase, yBase, qualTxt=False, public=clOpts.public)
    if den:
        den_eta_number_strs = re.findall(r'[\d\.\d]+', den[den.find('EtaMin')+6:den.find('EtaMax')+12])
        #if len(den_eta_number_strs) > 1 and eta_number_strs != den_eta_number_strs:
        if len(den_eta_number_strs) > 1:
            den_note_str = den_eta_number_strs[0]+' < |#eta^{GMT}| < '+den_eta_number_strs[1]
            notes.append([xBase, yBase-0.05, den_note_str, True])

    return plot_hists(hDefs, xTitle, yTitle, threshold, normToBinWidth, prefix, notes, scaleFactor, data=data, errorBars=False, clOpts=clOpts)

def plot_hists_comp(hm1, hm2, hName, xTitle='', yTitle='# muons', threshold=False, stacked=False, normToBinWidth=False, reg='', scaleFactor=1., data=False, legTxts=['hm1', 'hm2'], clOpts=None):
    styles = hist_styles(stacked)

    #hm1_dict = {'hm':hm1, 'num':hName, 'denhm':hm2, 'den':hName}
    hm1_dict = {'hm':hm1, 'num':hName, 'den':''}
    hm2_dict = {'hm':hm2, 'num':hName, 'den':''}
    hm1_dict.update(styles['hm1'])
    hm2_dict.update(styles['hm2'])
    hm1_dict['legtext'] = legTxts[0]
    hm2_dict['legtext'] = legTxts[1]
    hDefs = [hm1_dict, hm2_dict]

    xBase = 0.17
    yBase = 0.56
    notes = extract_notes_from_name(hName, xBase, yBase, public=clOpts.public)

    return plot_hists(hDefs, xTitle, yTitle, threshold, normToBinWidth, 'comp_', notes, scaleFactor, data, clOpts=clOpts)

def plot_hists_comp_ratio(hm1, hm2, hName, xTitle='', yTitle='# muons', threshold=False, stacked=False, normToBinWidth=False, reg='', scaleFactor=1., data=False, legTxts=['hm1', 'hm2'], clOpts=None):
    styles = hist_styles(stacked)

    hm_dict = {'hm':hm2, 'num':hName, 'denhm':hm1, 'den':hName}
    hm_dict.update(styles['hm2'])
    hm_dict['legtext'] = legTxts[1]+'/'+legTxts[0]
    hDefs = [hm_dict]

    xBase = 0.17
    yBase = 0.56
    notes = extract_notes_from_name(hName, xBase, yBase, public=clOpts.public)

    return plot_hists2(hDefs, [xTitle, yTitle], threshold, normToBinWidth, 'compratio_', notes, scaleFactor, data, clOpts=clOpts)

def plot_hists_comp2(hm1, hm2, hName, xTitle='', yTitle='# muons', threshold=False, stacked=False, normToBinWidth=False, reg='', scaleFactor=1., data=False, legTxts=['hm1', 'hm2'], clOpts=None):
    styles = hist_styles(stacked)

    hm1_dict = {'hm':hm1, 'num':hName, 'denhm':hm2, 'den':hName}
    hm2_dict = {'hm':hm2, 'num':hName, 'den':''}
    hm1_dict.update(styles['hm1'])
    hm2_dict.update(styles['hm2'])
    hm1_dict['legtext'] = legTxts[0]
    hm2_dict['legtext'] = legTxts[1]
    hDefs = [hm1_dict, hm2_dict]

    xBase = 0.17
    yBase = 0.56
    notes = extract_notes_from_name(hName, xBase, yBase, public=clOpts.public)

    return plot_hists(hDefs, xTitle, yTitle, threshold, normToBinWidth, 'comp_', notes, scaleFactor, data, clOpts=clOpts)

def main():
    opts = parse_options_plotRates(parser)
    batchRun = opts.interactive
    if batchRun:
        root.gROOT.SetBatch(True)

    hm = HistManager(filename=opts.fname)
    if opts.fname2:
        hm2 = HistManager(filename=opts.fname2)
        legTxt1 = opts.legtxt1
        legTxt2 = opts.legtxt2
 
    nEvtsAna = hm.get('n_evts_analysed').GetBinContent(1)
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
    if opts.scaleto > 0:
        convFactorToHz *= opts.scaleto / float(nCollBunches)
        print 'Conversion factor after scaling to {to} bunches: {convFactorToHz}'.format(to=opts.scaleto, convFactorToHz=convFactorToHz)
    if opts.scale != 1.:
        convFactorToHz *= opts.scale
        print 'Conversion factor after applying additinoal scale factor of {sf}: {convFactorToHz}'.format(sf=opts.scale, convFactorToHz=convFactorToHz)

    L1Ana.init_l1_analysis()
    print ""

    # holds the canvases, histograms, etc.
    objects = []

    ##########################################################################
    if opts.printrates:
        print 'Rates'
        print_rates(hm, 'highest_muon_absEtaMin0_absEtaMax2.5_qmin12_pt', scaleFactor=convFactorToHz / 1000.)
        print_rates(hm, 'highest_muon_absEtaMin0_absEtaMax2.1_qmin12_pt', scaleFactor=convFactorToHz / 1000.)
        #print_rates(hm, 'highest_muon_absEtaMin0_absEtaMax2.5_qmin4_pt', scaleFactor=convFactorToHz / 1000.)
        #print_rates(hm, 'highest_muon_absEtaMin0_absEtaMax2.1_qmin4_pt', scaleFactor=convFactorToHz / 1000.)

        print 'Rates scaled to 1 colliding bunch'
        print_rates(hm, 'highest_muon_absEtaMin0_absEtaMax2.5_qmin12_pt', scaleFactor=convFactorToHz / 1000. / nCollBunches)
        print_rates(hm, 'highest_muon_absEtaMin0_absEtaMax2.1_qmin12_pt', scaleFactor=convFactorToHz / 1000. / nCollBunches)
        #print_rates(hm, 'highest_muon_absEtaMin0_absEtaMax2.5_qmin4_pt', scaleFactor=convFactorToHz / 1000. / nCollBunches)
        #print_rates(hm, 'highest_muon_absEtaMin0_absEtaMax2.1_qmin4_pt', scaleFactor=convFactorToHz / 1000. / nCollBunches)

        if opts.fname2:
            print 'Rates'
            print_rates_comp(hm, hm2, 'highest_muon_absEtaMin0_absEtaMax2.5_qmin12_pt', scaleFactor=convFactorToHz / 1000., legTxts=[legTxt1, legTxt2])
            print_rates_comp(hm, hm2, 'highest_muon_absEtaMin0_absEtaMax2.1_qmin12_pt', scaleFactor=convFactorToHz / 1000., legTxts=[legTxt1, legTxt2])

    if opts.public:
        objects.append(plot_hists_standard(hm, 'highest_muon_absEtaMin0_absEtaMax2.5_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='Integrated # events', threshold=True, stacked=True, data=thisIsData, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'highest_muon_absEtaMin0_absEtaMax2.5_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='Integrated # events', threshold=True, data=thisIsData, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'highest_muon_absEtaMin0_absEtaMax2.5_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='kHz', threshold=True, scaleFactor=convFactorToHz / 1000., data=thisIsData, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'highest_muon_absEtaMin0_absEtaMax0.83_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='kHz', threshold=True, scaleFactor=convFactorToHz / 1000., data=thisIsData, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'highest_muon_absEtaMin0_absEtaMax0.83_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='Integrated # events', threshold=True, stacked=True, data=thisIsData, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'highest_muon_absEtaMin0.83_absEtaMax1.24_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='kHz', threshold=True, scaleFactor=convFactorToHz / 1000., data=thisIsData, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'highest_muon_absEtaMin0.83_absEtaMax1.24_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='Integrated # events', threshold=True, stacked=True, data=thisIsData, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'highest_muon_absEtaMin1.24_absEtaMax2.5_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='kHz', threshold=True, scaleFactor=convFactorToHz / 1000., data=thisIsData, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'highest_muon_absEtaMin1.24_absEtaMax2.5_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='Integrated # events', threshold=True, stacked=True, data=thisIsData, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'highest_muon_absEtaMin0_absEtaMax2.5_qmin8_pt', xTitle='p_{T} (GeV/c)', yTitle='kHz', threshold=True, scaleFactor=convFactorToHz / 1000., data=thisIsData, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'highest_muon_absEtaMin0_absEtaMax2.5_qmin4_pt', xTitle='p_{T} (GeV/c)', yTitle='kHz', threshold=True, scaleFactor=convFactorToHz / 1000., data=thisIsData, clOpts=opts))

        objects.append(plot_hists_standard(hm, 'muon_ptmin0_qmin12_eta', xTitle='#eta', yTitle='', stacked=True, normToBinWidth=True, data=thisIsData, scaleFactor=1./nZeroBiasEvents, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'muon_ptmin7_qmin12_eta', xTitle='#eta', yTitle='', stacked=True, normToBinWidth=True, data=thisIsData, scaleFactor=1./nZeroBiasEvents, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'muon_ptmin12_qmin12_eta', xTitle='#eta', yTitle='', stacked=True, normToBinWidth=True, data=thisIsData, scaleFactor=1./nZeroBiasEvents, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'muon_ptmin22_qmin12_eta', xTitle='#eta', yTitle='', stacked=True, normToBinWidth=True, data=thisIsData, scaleFactor=1./nZeroBiasEvents, clOpts=opts))

        if opts.legacy:
            # relative uGMT rates for regions single muon quality
            objects.append(plot_hists_standard(hm, 'highest_muon_absEtaMin0_absEtaMax2.5_qmin12_pt', 'gmt_highest_muon_absEtaMin0_absEtaMax2.5_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='Integrated # events / Integrated # GMT events', threshold=True, stacked=True, data=thisIsData, clOpts=opts))
            objects.append(plot_hists_standard(hm, 'highest_muon_absEtaMin0_absEtaMax0.83_qmin12_pt', 'gmt_highest_muon_absEtaMin0_absEtaMax0.83_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='Integrated # events / Integrated # GMT events', threshold=True, stacked=True, data=thisIsData, clOpts=opts))
            objects.append(plot_hists_standard(hm, 'highest_muon_absEtaMin0.83_absEtaMax1.24_qmin12_pt', 'gmt_highest_muon_absEtaMin0.83_absEtaMax1.24_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='Integrated # events / Integrated # GMT events', threshold=True, stacked=True, data=thisIsData, clOpts=opts))
            objects.append(plot_hists_standard(hm, 'highest_muon_absEtaMin1.24_absEtaMax2.5_qmin12_pt', 'gmt_highest_muon_absEtaMin1.24_absEtaMax2.5_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='Integrated # events / Integrated # GMT events', threshold=True, stacked=True, data=thisIsData, clOpts=opts))

            # relative uGMT eta distributions
            objects.append(plot_hists_standard(hm, 'muon_ptmin0_qmin12_eta', 'gmt_muon_ptmin0_qmin12_eta', xTitle='#eta', yTitle='# muons / # GMT muons', stacked=True, data=thisIsData, clOpts=opts))
            objects.append(plot_hists_standard(hm, 'muon_ptmin7_qmin12_eta', 'gmt_muon_ptmin7_qmin12_eta', xTitle='#eta', yTitle='# muons / # GMT muons', stacked=True, data=thisIsData, clOpts=opts))
            objects.append(plot_hists_standard(hm, 'muon_ptmin12_qmin12_eta', 'gmt_muon_ptmin12_qmin12_eta', xTitle='#eta', yTitle='# muons / # GMT muons', stacked=True, data=thisIsData, clOpts=opts))
            objects.append(plot_hists_standard(hm, 'muon_ptmin22_qmin12_eta', 'gmt_muon_ptmin22_qmin12_eta', xTitle='#eta', yTitle='# muons / # GMT muons', stacked=True, data=thisIsData, clOpts=opts))

    else:
        # uGMT kinematic variables
        objects.append(plot_hists_standard(hm, 'muon_absEtaMin0_absEtaMax2.5_qmin12_varBin_pt', xTitle='p_{T} (GeV/c)', yTitle='# muons/(GeV/c)', stacked=True, normToBinWidth=True, data=thisIsData, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'muon_absEtaMin0_absEtaMax0.83_qmin12_varBin_pt', xTitle='p_{T} (GeV/c)', yTitle='# muons/(GeV/c)', stacked=True, normToBinWidth=True, data=thisIsData, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'muon_absEtaMin0.83_absEtaMax1.24_qmin12_varBin_pt', xTitle='p_{T} (GeV/c)', yTitle='# muons/(GeV/c)', stacked=True, normToBinWidth=True, data=thisIsData, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'muon_absEtaMin1.24_absEtaMax2.5_qmin12_varBin_pt', xTitle='p_{T} (GeV/c)', yTitle='# muons/(GeV/c)', stacked=True, normToBinWidth=True, data=thisIsData, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'muon_ptmin0_qmin12_eta', xTitle='#eta', yTitle='# muons', stacked=True, data=thisIsData, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'muon_ptmin22_qmin12_eta', xTitle='#eta', yTitle='# muons', stacked=True, data=thisIsData, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'muon_ptmin0_qmin12_eta', xTitle='#eta', yTitle='', stacked=True, normToBinWidth=True, data=thisIsData, scaleFactor=1./nZeroBiasEvents, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'muon_ptmin3_qmin12_eta', xTitle='#eta', yTitle='', stacked=True, normToBinWidth=True, data=thisIsData, scaleFactor=1./nZeroBiasEvents, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'muon_ptmin5_qmin12_eta', xTitle='#eta', yTitle='', stacked=True, normToBinWidth=True, data=thisIsData, scaleFactor=1./nZeroBiasEvents, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'muon_ptmin7_qmin12_eta', xTitle='#eta', yTitle='', stacked=True, normToBinWidth=True, data=thisIsData, scaleFactor=1./nZeroBiasEvents, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'muon_ptmin12_qmin12_eta', xTitle='#eta', yTitle='', stacked=True, normToBinWidth=True, data=thisIsData, scaleFactor=1./nZeroBiasEvents, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'muon_ptmin18_qmin12_eta', xTitle='#eta', yTitle='', stacked=True, normToBinWidth=True, data=thisIsData, scaleFactor=1./nZeroBiasEvents, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'muon_ptmin22_qmin12_eta', xTitle='#eta', yTitle='', stacked=True, normToBinWidth=True, data=thisIsData, scaleFactor=1./nZeroBiasEvents, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'muon_ptmin0_qmin8_eta', xTitle='#eta', yTitle='', stacked=True, normToBinWidth=True, data=thisIsData, scaleFactor=1./nZeroBiasEvents, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'muon_ptmin3_qmin8_eta', xTitle='#eta', yTitle='', stacked=True, normToBinWidth=True, data=thisIsData, scaleFactor=1./nZeroBiasEvents, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'muon_ptmin5_qmin8_eta', xTitle='#eta', yTitle='', stacked=True, normToBinWidth=True, data=thisIsData, scaleFactor=1./nZeroBiasEvents, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'muon_ptmin7_qmin8_eta', xTitle='#eta', yTitle='', stacked=True, normToBinWidth=True, data=thisIsData, scaleFactor=1./nZeroBiasEvents, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'muon_ptmin12_qmin8_eta', xTitle='#eta', yTitle='', stacked=True, normToBinWidth=True, data=thisIsData, scaleFactor=1./nZeroBiasEvents, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'muon_ptmin18_qmin8_eta', xTitle='#eta', yTitle='', stacked=True, normToBinWidth=True, data=thisIsData, scaleFactor=1./nZeroBiasEvents, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'muon_ptmin22_qmin8_eta', xTitle='#eta', yTitle='', stacked=True, normToBinWidth=True, data=thisIsData, scaleFactor=1./nZeroBiasEvents, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'muon_ptmin0_qmin4_eta', xTitle='#eta', yTitle='', stacked=True, normToBinWidth=True, data=thisIsData, scaleFactor=1./nZeroBiasEvents, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'muon_ptmin3_qmin4_eta', xTitle='#eta', yTitle='', stacked=True, normToBinWidth=True, data=thisIsData, scaleFactor=1./nZeroBiasEvents, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'muon_ptmin5_qmin4_eta', xTitle='#eta', yTitle='', stacked=True, normToBinWidth=True, data=thisIsData, scaleFactor=1./nZeroBiasEvents, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'muon_ptmin7_qmin4_eta', xTitle='#eta', yTitle='', stacked=True, normToBinWidth=True, data=thisIsData, scaleFactor=1./nZeroBiasEvents, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'muon_ptmin12_qmin4_eta', xTitle='#eta', yTitle='', stacked=True, normToBinWidth=True, data=thisIsData, scaleFactor=1./nZeroBiasEvents, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'muon_ptmin18_qmin4_eta', xTitle='#eta', yTitle='', stacked=True, normToBinWidth=True, data=thisIsData, scaleFactor=1./nZeroBiasEvents, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'muon_ptmin22_qmin4_eta', xTitle='#eta', yTitle='', stacked=True, normToBinWidth=True, data=thisIsData, scaleFactor=1./nZeroBiasEvents, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'muon_absEtaMin0_absEtaMax2.5_ptmin0_qmin12_phi', xTitle='#phi', yTitle='# muons', stacked=True, data=thisIsData, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'muon_absEtaMin0_absEtaMax2.5_ptmin3_qmin12_phi', xTitle='#phi', yTitle='# muons', stacked=True, data=thisIsData, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'muon_absEtaMin0_absEtaMax2.5_ptmin5_qmin12_phi', xTitle='#phi', yTitle='# muons', stacked=True, data=thisIsData, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'muon_absEtaMin0_absEtaMax2.5_ptmin7_qmin12_phi', xTitle='#phi', yTitle='# muons', stacked=True, data=thisIsData, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'muon_absEtaMin0_absEtaMax2.5_ptmin12_qmin12_phi', xTitle='#phi', yTitle='# muons', stacked=True, data=thisIsData, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'muon_absEtaMin0_absEtaMax2.5_ptmin18_qmin12_phi', xTitle='#phi', yTitle='# muons', stacked=True, data=thisIsData, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'muon_absEtaMin0_absEtaMax2.5_ptmin22_qmin12_phi', xTitle='#phi', yTitle='# muons', stacked=True, data=thisIsData, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'muon_absEtaMin0_absEtaMax2.5_ptmin0_qual', xTitle='#mu quality', yTitle='# muons', stacked=True, data=thisIsData, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'muon_absEtaMin0_absEtaMax2.5_ptmin3_qual', xTitle='#mu quality', yTitle='# muons', stacked=True, data=thisIsData, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'muon_absEtaMin0_absEtaMax2.5_ptmin5_qual', xTitle='#mu quality', yTitle='# muons', stacked=True, data=thisIsData, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'muon_absEtaMin0_absEtaMax2.5_ptmin7_qual', xTitle='#mu quality', yTitle='# muons', stacked=True, data=thisIsData, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'muon_absEtaMin0_absEtaMax2.5_ptmin12_qual', xTitle='#mu quality', yTitle='# muons', stacked=True, data=thisIsData, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'muon_absEtaMin0_absEtaMax2.5_ptmin18_qual', xTitle='#mu quality', yTitle='# muons', stacked=True, data=thisIsData, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'muon_absEtaMin0_absEtaMax2.5_ptmin22_qual', xTitle='#mu quality', yTitle='# muons', stacked=True, data=thisIsData, clOpts=opts))

        # uGMT rates for regions
        objects.append(plot_hists_standard(hm, 'highest_muon_absEtaMin0_absEtaMax2.5_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='Integrated # events', threshold=True, stacked=True, data=thisIsData, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'highest_muon_absEtaMin0_absEtaMax2.5_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='Integrated # events', threshold=True, data=thisIsData, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'highest_muon_absEtaMin0_absEtaMax2.5_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='kHz', threshold=True, scaleFactor=convFactorToHz / 1000., data=thisIsData, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'highest_muon_absEtaMin0_absEtaMax2.1_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='kHz', threshold=True, scaleFactor=convFactorToHz / 1000., data=thisIsData, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'highest_muon_absEtaMin0_absEtaMax0.83_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='kHz', threshold=True, scaleFactor=convFactorToHz / 1000., data=thisIsData, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'highest_muon_absEtaMin0_absEtaMax0.83_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='Integrated # events', threshold=True, stacked=True, data=thisIsData, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'highest_muon_absEtaMin0.83_absEtaMax1.24_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='kHz', threshold=True, scaleFactor=convFactorToHz / 1000., data=thisIsData, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'highest_muon_absEtaMin0.83_absEtaMax1.24_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='Integrated # events', threshold=True, stacked=True, data=thisIsData, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'highest_muon_absEtaMin1.24_absEtaMax2.5_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='kHz', threshold=True, scaleFactor=convFactorToHz / 1000., data=thisIsData, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'highest_muon_absEtaMin1.24_absEtaMax2.5_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='Integrated # events', threshold=True, stacked=True, data=thisIsData, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'highest_muon_absEtaMin0_absEtaMax2.5_qmin8_pt', xTitle='p_{T} (GeV/c)', yTitle='kHz', threshold=True, scaleFactor=convFactorToHz / 1000., data=thisIsData, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'highest_muon_absEtaMin0_absEtaMax2.1_qmin8_pt', xTitle='p_{T} (GeV/c)', yTitle='kHz', threshold=True, scaleFactor=convFactorToHz / 1000., data=thisIsData, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'highest_muon_absEtaMin0_absEtaMax2.5_qmin4_pt', xTitle='p_{T} (GeV/c)', yTitle='kHz', threshold=True, scaleFactor=convFactorToHz / 1000., data=thisIsData, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'highest_muon_absEtaMin0_absEtaMax2.1_qmin4_pt', xTitle='p_{T} (GeV/c)', yTitle='kHz', threshold=True, scaleFactor=convFactorToHz / 1000., data=thisIsData, clOpts=opts))

        if opts.fname2:
            objects.append(plot_hists_comp(hm, hm2, 'ugmt_highest_muon_absEtaMin0_absEtaMax2.5_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='kHz', threshold=True, scaleFactor=convFactorToHz / 1000., data=thisIsData, legTxts=[legTxt1, legTxt2], clOpts=opts))
            objects.append(plot_hists_comp(hm, hm2, 'bmtf_ugmt_highest_muon_absEtaMin0_absEtaMax2.5_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='kHz', threshold=True, scaleFactor=convFactorToHz / 1000., data=thisIsData, legTxts=['BMTF '+legTxt1, 'BMTF '+legTxt2], clOpts=opts))
            objects.append(plot_hists_comp(hm, hm2, 'omtf_ugmt_highest_muon_absEtaMin0_absEtaMax2.5_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='kHz', threshold=True, scaleFactor=convFactorToHz / 1000., data=thisIsData, legTxts=['OMTF '+legTxt1, 'OMTF '+legTxt2], clOpts=opts))
            objects.append(plot_hists_comp(hm, hm2, 'emtf_ugmt_highest_muon_absEtaMin0_absEtaMax2.5_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='kHz', threshold=True, scaleFactor=convFactorToHz / 1000., data=thisIsData, legTxts=['EMTF '+legTxt1, 'EMTF '+legTxt2], clOpts=opts))
            objects.append(plot_hists_comp(hm, hm2, 'ugmt_highest_muon_absEtaMin0_absEtaMax2.5_qmin8_pt', xTitle='p_{T} (GeV/c)', yTitle='kHz', threshold=True, scaleFactor=convFactorToHz / 1000., data=thisIsData, legTxts=[legTxt1, legTxt2], clOpts=opts))
            objects.append(plot_hists_comp(hm, hm2, 'bmtf_ugmt_highest_muon_absEtaMin0_absEtaMax2.5_qmin8_pt', xTitle='p_{T} (GeV/c)', yTitle='kHz', threshold=True, scaleFactor=convFactorToHz / 1000., data=thisIsData, legTxts=['BMTF '+legTxt1, 'BMTF '+legTxt2], clOpts=opts))
            objects.append(plot_hists_comp(hm, hm2, 'omtf_ugmt_highest_muon_absEtaMin0_absEtaMax2.5_qmin8_pt', xTitle='p_{T} (GeV/c)', yTitle='kHz', threshold=True, scaleFactor=convFactorToHz / 1000., data=thisIsData, legTxts=['OMTF '+legTxt1, 'OMTF '+legTxt2], clOpts=opts))
            objects.append(plot_hists_comp(hm, hm2, 'emtf_ugmt_highest_muon_absEtaMin0_absEtaMax2.5_qmin8_pt', xTitle='p_{T} (GeV/c)', yTitle='kHz', threshold=True, scaleFactor=convFactorToHz / 1000., data=thisIsData, legTxts=['EMTF '+legTxt1, 'EMTF '+legTxt2], clOpts=opts))
            objects.append(plot_hists_comp(hm, hm2, 'ugmt_highest_muon_absEtaMin0_absEtaMax2.5_qmin4_pt', xTitle='p_{T} (GeV/c)', yTitle='kHz', threshold=True, scaleFactor=convFactorToHz / 1000., data=thisIsData, legTxts=[legTxt1, legTxt2], clOpts=opts))
            objects.append(plot_hists_comp(hm, hm2, 'bmtf_ugmt_highest_muon_absEtaMin0_absEtaMax2.5_qmin4_pt', xTitle='p_{T} (GeV/c)', yTitle='kHz', threshold=True, scaleFactor=convFactorToHz / 1000., data=thisIsData, legTxts=['BMTF '+legTxt1, 'BMTF '+legTxt2], clOpts=opts))
            objects.append(plot_hists_comp(hm, hm2, 'omtf_ugmt_highest_muon_absEtaMin0_absEtaMax2.5_qmin4_pt', xTitle='p_{T} (GeV/c)', yTitle='kHz', threshold=True, scaleFactor=convFactorToHz / 1000., data=thisIsData, legTxts=['OMTF '+legTxt1, 'OMTF '+legTxt2], clOpts=opts))
            objects.append(plot_hists_comp(hm, hm2, 'emtf_ugmt_highest_muon_absEtaMin0_absEtaMax2.5_qmin4_pt', xTitle='p_{T} (GeV/c)', yTitle='kHz', threshold=True, scaleFactor=convFactorToHz / 1000., data=thisIsData, legTxts=['EMTF '+legTxt1, 'EMTF '+legTxt2], clOpts=opts))

            objects.append(plot_hists_comp(hm, hm2, 'ugmt_muon_ptmin22_qmin12_eta', xTitle='#eta', yTitle='', scaleFactor=1./nZeroBiasEvents, data=thisIsData, legTxts=[legTxt1, legTxt2], clOpts=opts))
            objects.append(plot_hists_comp(hm, hm2, 'ugmt_muon_ptmin12_qmin8_eta', xTitle='#eta', yTitle='', scaleFactor=1./nZeroBiasEvents, data=thisIsData, legTxts=[legTxt1, legTxt2], clOpts=opts))

            objects.append(plot_hists_comp_ratio(hm, hm2, 'ugmt_highest_muon_absEtaMin0_absEtaMax2.5_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle=legTxt2+'/'+legTxt1, threshold=True, data=thisIsData, legTxts=[legTxt1, legTxt2], clOpts=opts))
            objects.append(plot_hists_comp_ratio(hm, hm2, 'bmtf_ugmt_highest_muon_absEtaMin0_absEtaMax2.5_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='BMTF '+legTxt2+'/BMTF_'+legTxt1, threshold=True, data=thisIsData, legTxts=['BMTF '+legTxt1, 'BMTF '+legTxt2], clOpts=opts))
            objects.append(plot_hists_comp_ratio(hm, hm2, 'omtf_ugmt_highest_muon_absEtaMin0_absEtaMax2.5_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='OMTF '+legTxt2+'/OMTF_'+legTxt1, threshold=True, data=thisIsData, legTxts=['OMTF '+legTxt1, 'OMTF '+legTxt2], clOpts=opts))
            objects.append(plot_hists_comp_ratio(hm, hm2, 'emtf_ugmt_highest_muon_absEtaMin0_absEtaMax2.5_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='EMTF '+legTxt2+'/EMTF_'+legTxt1, threshold=True, data=thisIsData, legTxts=['EMTF '+legTxt1, 'EMTF '+legTxt2], clOpts=opts))
            objects.append(plot_hists_comp_ratio(hm, hm2, 'ugmt_highest_muon_absEtaMin0_absEtaMax2.5_qmin8_pt', xTitle='p_{T} (GeV/c)', yTitle=legTxt2+'/'+legTxt1, threshold=True, data=thisIsData, legTxts=[legTxt1, legTxt2], clOpts=opts))
            objects.append(plot_hists_comp_ratio(hm, hm2, 'bmtf_ugmt_highest_muon_absEtaMin0_absEtaMax2.5_qmin8_pt', xTitle='p_{T} (GeV/c)', yTitle='BMTF '+legTxt2+'/BMTF_'+legTxt1, threshold=True, data=thisIsData, legTxts=['BMTF '+legTxt1, 'BMTF '+legTxt2], clOpts=opts))
            objects.append(plot_hists_comp_ratio(hm, hm2, 'omtf_ugmt_highest_muon_absEtaMin0_absEtaMax2.5_qmin8_pt', xTitle='p_{T} (GeV/c)', yTitle='OMTF '+legTxt2+'/OMTF_'+legTxt1, threshold=True, data=thisIsData, legTxts=['OMTF '+legTxt1, 'OMTF '+legTxt2], clOpts=opts))
            objects.append(plot_hists_comp_ratio(hm, hm2, 'emtf_ugmt_highest_muon_absEtaMin0_absEtaMax2.5_qmin8_pt', xTitle='p_{T} (GeV/c)', yTitle='EMTF '+legTxt2+'/EMTF_'+legTxt1, threshold=True, data=thisIsData, legTxts=['EMTF '+legTxt1, 'EMTF '+legTxt2], clOpts=opts))
            objects.append(plot_hists_comp_ratio(hm, hm2, 'ugmt_highest_muon_absEtaMin0_absEtaMax2.5_qmin4_pt', xTitle='p_{T} (GeV/c)', yTitle=legTxt2+'/'+legTxt1, threshold=True, data=thisIsData, legTxts=[legTxt1, legTxt2], clOpts=opts))
            objects.append(plot_hists_comp_ratio(hm, hm2, 'bmtf_ugmt_highest_muon_absEtaMin0_absEtaMax2.5_qmin4_pt', xTitle='p_{T} (GeV/c)', yTitle='BMTF '+legTxt2+'/BMTF_'+legTxt1, threshold=True, data=thisIsData, legTxts=['BMTF '+legTxt1, 'BMTF '+legTxt2], clOpts=opts))
            objects.append(plot_hists_comp_ratio(hm, hm2, 'omtf_ugmt_highest_muon_absEtaMin0_absEtaMax2.5_qmin4_pt', xTitle='p_{T} (GeV/c)', yTitle='OMTF '+legTxt2+'/OMTF_'+legTxt1, threshold=True, data=thisIsData, legTxts=['OMTF '+legTxt1, 'OMTF '+legTxt2], clOpts=opts))
            objects.append(plot_hists_comp_ratio(hm, hm2, 'emtf_ugmt_highest_muon_absEtaMin0_absEtaMax2.5_qmin4_pt', xTitle='p_{T} (GeV/c)', yTitle='EMTF '+legTxt2+'/EMTF_'+legTxt1, threshold=True, data=thisIsData, legTxts=['EMTF '+legTxt1, 'EMTF '+legTxt2], clOpts=opts))

            objects.append(plot_hists_comp_ratio(hm, hm2, 'ugmt_muon_ptmin22_qmin12_eta', xTitle='#eta', yTitle=legTxt2+'/'+legTxt1, data=thisIsData, legTxts=[legTxt1, legTxt2], clOpts=opts))
            objects.append(plot_hists_comp_ratio(hm, hm2, 'ugmt_muon_ptmin12_qmin8_eta', xTitle='#eta', yTitle=legTxt2+'/'+legTxt1, data=thisIsData, legTxts=[legTxt1, legTxt2], clOpts=opts))

        if opts.qstack:
            # q stack uGMT rates for regions
            objects.append(plot_hists_qstack(hm, 'muon_absEtaMin0_absEtaMax2.5_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='Integrated # events', threshold=True, stacked=True, data=thisIsData, clOpts=opts))
            objects.append(plot_hists_qstack(hm, 'muon_absEtaMin0_absEtaMax2.5_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='kHz', threshold=True, stacked=True, scaleFactor=convFactorToHz / 1000., data=thisIsData, clOpts=opts))
            objects.append(plot_hists_qstack(hm, 'muon_absEtaMin0_absEtaMax0.83_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='Integrated # events', threshold=True, stacked=True, data=thisIsData, clOpts=opts))
            objects.append(plot_hists_qstack(hm, 'muon_absEtaMin0.83_absEtaMax1.24_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='Integrated # events', threshold=True, stacked=True, data=thisIsData, clOpts=opts))
            objects.append(plot_hists_qstack(hm, 'muon_absEtaMin1.24_absEtaMax2.5_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='Integrated # events', threshold=True, stacked=True, data=thisIsData, clOpts=opts))
            objects.append(plot_hists_qstack(hm, 'muon_ptmin0_qmin12_eta', xTitle='#eta', yTitle='', stacked=True, normToBinWidth=True, data=thisIsData, clOpts=opts))
            objects.append(plot_hists_qstack(hm, 'muon_ptmin3_qmin12_eta', xTitle='#eta', yTitle='', stacked=True, normToBinWidth=True, data=thisIsData, clOpts=opts))
            objects.append(plot_hists_qstack(hm, 'muon_ptmin5_qmin12_eta', xTitle='#eta', yTitle='', stacked=True, normToBinWidth=True, data=thisIsData, clOpts=opts))
            objects.append(plot_hists_qstack(hm, 'muon_ptmin7_qmin12_eta', xTitle='#eta', yTitle='', stacked=True, normToBinWidth=True, data=thisIsData, clOpts=opts))
            objects.append(plot_hists_qstack(hm, 'muon_ptmin12_qmin12_eta', xTitle='#eta', yTitle='', stacked=True, normToBinWidth=True, data=thisIsData, clOpts=opts))
            objects.append(plot_hists_qstack(hm, 'muon_ptmin18_qmin12_eta', xTitle='#eta', yTitle='', stacked=True, normToBinWidth=True, data=thisIsData, clOpts=opts))
            objects.append(plot_hists_qstack(hm, 'muon_ptmin22_qmin12_eta', xTitle='#eta', yTitle='', stacked=True, normToBinWidth=True, data=thisIsData, clOpts=opts))

        if opts.legacy:
            # relative uGMT rates for regions single muon quality
            objects.append(plot_hists_standard(hm, 'highest_muon_absEtaMin0_absEtaMax2.5_qmin12_pt', 'gmt_highest_muon_absEtaMin0_absEtaMax2.5_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='Integrated # events / Integrated # GMT events', threshold=True, stacked=True, data=thisIsData, clOpts=opts))
            objects.append(plot_hists_standard(hm, 'highest_muon_absEtaMin0_absEtaMax0.83_qmin12_pt', 'gmt_highest_muon_absEtaMin0_absEtaMax0.83_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='Integrated # events / Integrated # GMT events', threshold=True, stacked=True, data=thisIsData, clOpts=opts))
            objects.append(plot_hists_standard(hm, 'highest_muon_absEtaMin0.83_absEtaMax1.24_qmin12_pt', 'gmt_highest_muon_absEtaMin0.83_absEtaMax1.24_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='Integrated # events / Integrated # GMT events', threshold=True, stacked=True, data=thisIsData, clOpts=opts))
            objects.append(plot_hists_standard(hm, 'highest_muon_absEtaMin1.24_absEtaMax2.5_qmin12_pt', 'gmt_highest_muon_absEtaMin1.24_absEtaMax2.5_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='Integrated # events / Integrated # GMT events', threshold=True, stacked=True, data=thisIsData, clOpts=opts))
            objects.append(plot_hists_standard(hm, 'highest_muon_absEtaMin0_absEtaMax0.83_qmin12_pt', 'gmt_highest_muon_absEtaMin0_absEtaMax0.83_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='Integrated # events / Integrated # GMT events', threshold=True, stacked=True, reg='b', data=thisIsData, clOpts=opts))
            objects.append(plot_hists_standard(hm, 'highest_muon_absEtaMin0.83_absEtaMax1.24_qmin12_pt', 'gmt_highest_muon_absEtaMin0.83_absEtaMax1.24_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='Integrated # events / Integrated # GMT events', threshold=True, stacked=True, reg='o', data=thisIsData, clOpts=opts))
            objects.append(plot_hists_standard(hm, 'highest_muon_absEtaMin1.24_absEtaMax2.5_qmin12_pt', 'gmt_highest_muon_absEtaMin1.24_absEtaMax2.5_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='Integrated # events / Integrated # GMT events', threshold=True, stacked=True, reg='e', data=thisIsData, clOpts=opts))

            # relative uGMT rates for regions double muon quality
            objects.append(plot_hists_standard(hm, 'highest_muon_absEtaMin0_absEtaMax2.5_qmin8_pt', 'gmt_highest_muon_absEtaMin0_absEtaMax2.5_qmin8_pt', xTitle='p_{T} (GeV/c)', yTitle='Integrated # events / Integrated # GMT events', threshold=True, stacked=True, data=thisIsData, clOpts=opts))
            objects.append(plot_hists_standard(hm, 'highest_muon_absEtaMin0_absEtaMax0.83_qmin8_pt', 'gmt_highest_muon_absEtaMin0_absEtaMax0.83_qmin8_pt', xTitle='p_{T} (GeV/c)', yTitle='Integrated # events / Integrated # GMT events', threshold=True, stacked=True, data=thisIsData, clOpts=opts))
            objects.append(plot_hists_standard(hm, 'highest_muon_absEtaMin0.83_absEtaMax1.24_qmin8_pt', 'gmt_highest_muon_absEtaMin0.83_absEtaMax1.24_qmin8_pt', xTitle='p_{T} (GeV/c)', yTitle='Integrated # events / Integrated # GMT events', threshold=True, stacked=True, data=thisIsData, clOpts=opts))
            objects.append(plot_hists_standard(hm, 'highest_muon_absEtaMin1.24_absEtaMax2.5_qmin8_pt', 'gmt_highest_muon_absEtaMin1.24_absEtaMax2.5_qmin8_pt', xTitle='p_{T} (GeV/c)', yTitle='Integrated # events / Integrated # GMT events', threshold=True, stacked=True, data=thisIsData, clOpts=opts))
            objects.append(plot_hists_standard(hm, 'highest_muon_absEtaMin0_absEtaMax0.83_qmin8_pt', 'gmt_highest_muon_absEtaMin0_absEtaMax0.83_qmin8_pt', xTitle='p_{T} (GeV/c)', yTitle='Integrated # events / Integrated # GMT events', threshold=True, stacked=True, reg='b', data=thisIsData, clOpts=opts))
            objects.append(plot_hists_standard(hm, 'highest_muon_absEtaMin0.83_absEtaMax1.24_qmin8_pt', 'gmt_highest_muon_absEtaMin0.83_absEtaMax1.24_qmin8_pt', xTitle='p_{T} (GeV/c)', yTitle='Integrated # events / Integrated # GMT events', threshold=True, stacked=True, reg='o', data=thisIsData, clOpts=opts))
            objects.append(plot_hists_standard(hm, 'highest_muon_absEtaMin1.24_absEtaMax2.5_qmin8_pt', 'gmt_highest_muon_absEtaMin1.24_absEtaMax2.5_qmin8_pt', xTitle='p_{T} (GeV/c)', yTitle='Integrated # events / Integrated # GMT events', threshold=True, stacked=True, reg='e', data=thisIsData, clOpts=opts))

            # relative uGMT rates for regions muon open quality
            objects.append(plot_hists_standard(hm, 'highest_muon_absEtaMin0_absEtaMax2.5_qmin4_pt', 'gmt_highest_muon_absEtaMin0_absEtaMax2.5_qmin4_pt', xTitle='p_{T} (GeV/c)', yTitle='Integrated # events / Integrated # GMT events', threshold=True, stacked=True, data=thisIsData, clOpts=opts))
            objects.append(plot_hists_standard(hm, 'highest_muon_absEtaMin0_absEtaMax0.83_qmin4_pt', 'gmt_highest_muon_absEtaMin0_absEtaMax0.83_qmin4_pt', xTitle='p_{T} (GeV/c)', yTitle='Integrated # events / Integrated # GMT events', threshold=True, stacked=True, data=thisIsData, clOpts=opts))
            objects.append(plot_hists_standard(hm, 'highest_muon_absEtaMin0.83_absEtaMax1.24_qmin4_pt', 'gmt_highest_muon_absEtaMin0.83_absEtaMax1.24_qmin4_pt', xTitle='p_{T} (GeV/c)', yTitle='Integrated # events / Integrated # GMT events', threshold=True, stacked=True, data=thisIsData, clOpts=opts))
            objects.append(plot_hists_standard(hm, 'highest_muon_absEtaMin1.24_absEtaMax2.5_qmin4_pt', 'gmt_highest_muon_absEtaMin1.24_absEtaMax2.5_qmin4_pt', xTitle='p_{T} (GeV/c)', yTitle='Integrated # events / Integrated # GMT events', threshold=True, stacked=True, data=thisIsData, clOpts=opts))
            objects.append(plot_hists_standard(hm, 'highest_muon_absEtaMin0_absEtaMax0.83_qmin4_pt', 'gmt_highest_muon_absEtaMin0_absEtaMax0.83_qmin4_pt', xTitle='p_{T} (GeV/c)', yTitle='Integrated # events / Integrated # GMT events', threshold=True, stacked=True, reg='b', data=thisIsData, clOpts=opts))
            objects.append(plot_hists_standard(hm, 'highest_muon_absEtaMin0.83_absEtaMax1.24_qmin4_pt', 'gmt_highest_muon_absEtaMin0.83_absEtaMax1.24_qmin4_pt', xTitle='p_{T} (GeV/c)', yTitle='Integrated # events / Integrated # GMT events', threshold=True, stacked=True, reg='o', data=thisIsData, clOpts=opts))
            objects.append(plot_hists_standard(hm, 'highest_muon_absEtaMin1.24_absEtaMax2.5_qmin4_pt', 'gmt_highest_muon_absEtaMin1.24_absEtaMax2.5_qmin4_pt', xTitle='p_{T} (GeV/c)', yTitle='Integrated # events / Integrated # GMT events', threshold=True, stacked=True, reg='e', data=thisIsData, clOpts=opts))

            if opts.qstack:
                # q stack relative uGMT rates for regions
                objects.append(plot_hists_qstack(hm, 'muon_absEtaMin0_absEtaMax2.5_qmin12_pt', 'gmt_muon_absEtaMin0_absEtaMax2.5_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='Integrated # events / Integrated # GMT events', threshold=True, stacked=True, data=thisIsData, clOpts=opts))
                objects.append(plot_hists_qstack(hm, 'muon_absEtaMin0_absEtaMax0.83_qmin12_pt', 'gmt_muon_absEtaMin0_absEtaMax0.83_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='Integrated # events / Integrated # GMT events', threshold=True, stacked=True, data=thisIsData, clOpts=opts))
                objects.append(plot_hists_qstack(hm, 'muon_absEtaMin0.83_absEtaMax1.24_qmin12_pt', 'gmt_muon_absEtaMin0.83_absEtaMax1.24_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='Integrated # events / Integrated # GMT events', threshold=True, stacked=True, data=thisIsData, clOpts=opts))
                objects.append(plot_hists_qstack(hm, 'muon_absEtaMin1.24_absEtaMax2.5_qmin12_pt', 'gmt_muon_absEtaMin1.24_absEtaMax2.5_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='Integrated # events / Integrated # GMT events', threshold=True, stacked=True, data=thisIsData, clOpts=opts))
                objects.append(plot_hists_qstack(hm, 'muon_absEtaMin0_absEtaMax0.83_qmin12_pt', 'gmt_muon_absEtaMin0_absEtaMax0.83_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='Integrated # events / Integrated # GMT events', threshold=True, stacked=True, reg='b', data=thisIsData, clOpts=opts))
                objects.append(plot_hists_qstack(hm, 'muon_absEtaMin0.83_absEtaMax1.24_qmin12_pt', 'gmt_muon_absEtaMin0.83_absEtaMax1.24_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='Integrated # events / Integrated # GMT events', threshold=True, stacked=True, reg='o', data=thisIsData, clOpts=opts))
                objects.append(plot_hists_qstack(hm, 'muon_absEtaMin1.24_absEtaMax2.5_qmin12_pt', 'gmt_muon_absEtaMin1.24_absEtaMax2.5_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='Integrated # events / Integrated # GMT events', threshold=True, stacked=True, reg='e', data=thisIsData, clOpts=opts))

        # uGMT TF rates for 0<|eta|<2.5
        objects.append(plot_hists_standard(hm, 'highest_muon_absEtaMin0_absEtaMax0.83_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='Integrated # events', threshold=True, stacked=True, reg='b', data=thisIsData, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'highest_muon_absEtaMin0.83_absEtaMax1.24_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='Integrated # events', threshold=True, stacked=True, reg='o', data=thisIsData, clOpts=opts))
        objects.append(plot_hists_standard(hm, 'highest_muon_absEtaMin1.24_absEtaMax2.5_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='Integrated # events', threshold=True, stacked=True, reg='e', data=thisIsData, clOpts=opts))

        if opts.qstack:
            # q stack uGMT TF rates for 0<|eta|<2.5
            objects.append(plot_hists_qstack(hm, 'muon_absEtaMin0_absEtaMax0.83_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='Integrated # events', threshold=True, stacked=True, reg='b', data=thisIsData, clOpts=opts))
            objects.append(plot_hists_qstack(hm, 'muon_absEtaMin0.83_absEtaMax1.24_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='Integrated # events', threshold=True, stacked=True, reg='o', data=thisIsData, clOpts=opts))
            objects.append(plot_hists_qstack(hm, 'muon_absEtaMin1.24_absEtaMax2.5_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='Integrated # events', threshold=True, stacked=True, reg='e', data=thisIsData, clOpts=opts))

        if opts.legacy:
            # relative uGMT TF rates for 0<|eta|<2.5
            objects.append(plot_hists_standard(hm, 'highest_muon_absEtaMin0_absEtaMax2.5_qmin12_pt', 'gmt_highest_muon_absEtaMin0_absEtaMax0.83_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='Integrated # events / Integrated # GMT events', threshold=True, stacked=True, reg='b', data=thisIsData, clOpts=opts))
            objects.append(plot_hists_standard(hm, 'highest_muon_absEtaMin0_absEtaMax2.5_qmin12_pt', 'gmt_highest_muon_absEtaMin0.83_absEtaMax1.24_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='Integrated # events / Integrated # GMT events', threshold=True, stacked=True, reg='o', data=thisIsData, clOpts=opts))
            objects.append(plot_hists_standard(hm, 'highest_muon_absEtaMin0_absEtaMax2.5_qmin12_pt', 'gmt_highest_muon_absEtaMin1.24_absEtaMax2.5_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='Integrated # events / Integrated # GMT events', threshold=True, stacked=True, reg='e', data=thisIsData, clOpts=opts))

            if opts.qstack:
                # q stack relative uGMT TF rates for 0<|eta|<2.5
                objects.append(plot_hists_qstack(hm, 'muon_absEtaMin0_absEtaMax2.5_qmin12_pt', 'gmt_muon_absEtaMin0_absEtaMax0.83_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='Integrated # events / Integrated # GMT events', threshold=True, stacked=True, reg='b', data=thisIsData, clOpts=opts))
                objects.append(plot_hists_qstack(hm, 'muon_absEtaMin0_absEtaMax2.5_qmin12_pt', 'gmt_muon_absEtaMin0.83_absEtaMax1.24_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='Integrated # events / Integrated # GMT events', threshold=True, stacked=True, reg='o', data=thisIsData, clOpts=opts))
                objects.append(plot_hists_qstack(hm, 'muon_absEtaMin0_absEtaMax2.5_qmin12_pt', 'gmt_muon_absEtaMin1.24_absEtaMax2.5_qmin12_pt', xTitle='p_{T} (GeV/c)', yTitle='Integrated # events / Integrated # GMT events', threshold=True, stacked=True, reg='e', data=thisIsData, clOpts=opts))

            # relative uGMT eta distributions
            objects.append(plot_hists_standard(hm, 'muon_ptmin0_qmin12_eta', 'gmt_muon_ptmin0_qmin12_eta', xTitle='#eta', yTitle='# muons / # GMT muons', stacked=True, data=thisIsData, clOpts=opts))
            objects.append(plot_hists_standard(hm, 'muon_ptmin0_qmin12_eta', 'gmt_muon_ptmin0_qmin12_eta', xTitle='#eta', yTitle='# muons / # GMT muons', stacked=False, data=thisIsData, clOpts=opts))
            objects.append(plot_hists_standard(hm, 'muon_ptmin12_qmin12_eta', 'gmt_muon_ptmin12_qmin12_eta', xTitle='#eta', yTitle='# muons / # GMT muons', stacked=True, data=thisIsData, clOpts=opts))
            objects.append(plot_hists_standard(hm, 'muon_ptmin12_qmin12_eta', 'gmt_muon_ptmin12_qmin12_eta', xTitle='#eta', yTitle='# muons / # GMT muons', stacked=False, data=thisIsData, clOpts=opts))
            objects.append(plot_hists_standard(hm, 'muon_ptmin22_qmin12_eta', 'gmt_muon_ptmin22_qmin12_eta', xTitle='#eta', yTitle='# muons / # GMT muons', stacked=True, data=thisIsData, clOpts=opts))
            objects.append(plot_hists_standard(hm, 'muon_ptmin22_qmin12_eta', 'gmt_muon_ptmin22_qmin12_eta', xTitle='#eta', yTitle='# muons / # GMT muons', stacked=False, data=thisIsData, clOpts=opts))

            if opts.qstack:
                objects.append(plot_hists_qstack  (hm, 'muon_ptmin0_qmin12_eta', 'gmt_muon_ptmin0_qmin12_eta', xTitle='#eta', yTitle='# muons / # GMT muons', stacked=True, data=thisIsData, clOpts=opts))
                objects.append(plot_hists_qstack  (hm, 'muon_ptmin12_qmin12_eta', 'gmt_muon_ptmin12_qmin12_eta', xTitle='#eta', yTitle='# muons / # GMT muons', stacked=True, data=thisIsData, clOpts=opts))
                objects.append(plot_hists_qstack  (hm, 'muon_ptmin22_qmin12_eta', 'gmt_muon_ptmin22_qmin12_eta', xTitle='#eta', yTitle='# muons / # GMT muons', stacked=True, data=thisIsData, clOpts=opts))

    ##########################################################################
    # save plots to root file
    if savePlots:
        plotdir = 'plots_'+opts.fname.replace('.root','').partition('/')[0]
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
    font = 42
    fontSize = 0.04
    main()

