#!/usr/bin/env python

import ROOT as root
import re

font = 42
fontSize = 0.04

# set the basic plotting style
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

# extract values from histogram name
def extract_notes_from_name(name, xBase, yBase, etaTxt=True, qualTxt=True, ptTxt=True):
    notes = []
    # extract eta range
    if etaTxt and name.rfind('.eta') == -1:
        eta_number_strs = re.findall(r'[\d\.\d]+', name[name.find('EtaMin')+6:name.find('EtaMax')+12])
        if len(eta_number_strs) > 1:
            # remove trailing dots
            for i in range(2):
                if eta_number_strs[i][-1] == '.':
                    eta_number_strs[i] = eta_number_strs[i][:-1]
            note_str = eta_number_strs[0]+' < |#eta| < '+eta_number_strs[1]
            notes.append([xBase, yBase+0.15, note_str, True])
    # extract quality
    if qualTxt and name.find('qmin') != -1:
        qualPos = name.find('qmin')
        qual_strs = re.findall(r'[\d]+', name[qualPos+4:qualPos+6])
        if len(qual_strs) > 0:
            qual_note_str = 'Quality #geq '+qual_strs[0]
            notes.append([xBase, yBase+0.1, qual_note_str, True])
    # extract pt range
    if ptTxt:
        ptminPos = name.find('ptmin')
        l1_ptmin_strs = re.findall(r'\d+\.?\d*', name[ptminPos+5:ptminPos+9])
        if len(l1_ptmin_strs) > 0:
            if l1_ptmin_strs[0][-1] == '.':
                l1_ptmin_strs[0] = l1_ptmin_strs[0][0:-1]
            notes.append([xBase, yBase+0.05, 'p_{T}^{L1} #geq '+l1_ptmin_strs[0]+' GeV', True])
    return notes

def plot_2dhist(hm2d, hName, drawDiag=True, data=False, xMax=None):
    canvas_name = hName
    if xMax:
        canvas_name += '_xMax{max}'.format(max=xMax)

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
    if xMax:
        h.GetXaxis().SetRangeUser(h.GetXaxis().GetBinLowEdge(1), xMax)
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

    if data:
        notes.append([0.53, 0.93, 'CMS internal, 13 TeV', False])
    else:
        notes.append([0.53, 0.93, 'CMS Simulation, 13 TeV', False])
    tex = add_text(notes=notes)

    lines = draw_tf_eta_regions(hName=hName, xMin=h.GetXaxis().GetXmin(), yMin=h.GetYaxis().GetXmin(), xMax=h.GetXaxis().GetXmax(), yMax=h.GetYaxis().GetXmax(), twoD=True, drawDiag=drawDiag)

    c.Modified()
    c.Update()
    
    return [c, h, tex, lines]


def plot_hists(hDefs, xTitle=None, yTitle='# muons', threshold=False, normToBinWidth=False, normalise=False, xMax=None, canvasPrefix='', notes=None, scaleFactor=1., logx=False, logy=False, data=False, rebin=1):
    den = hDefs[0]['den']
    if den:
        name = canvasPrefix+hDefs[0]['num']+'_over_'+den
    else:
        name = canvasPrefix+hDefs[0]['num']
    if normToBinWidth and not threshold and not den:
        name += '_normToBinWidth'
    if xMax:
        name += '_xMax{max}'.format(max=xMax)

    # setup legend according to how many histograms are in the plot
    legYmin = 0.9-0.04*len(hDefs)
    legXmin = 0.68
    legXmax = 0.9
    canvWidth = 600
    if legYmin < 0.6:
        legXmin = 0.8
        legXmax = 1.
        canvWidth = 730
    legend = root.TLegend(legXmin, legYmin, legXmax, 0.9)
    legend.SetTextFont(font)
    legend.SetTextSize(0.03)
    legend.SetBorderSize(0)
    legend.SetFillColor(19)
    legend.SetFillStyle(0)
    #legend.SetNColumns(2)
    legEntries = []

    hs = []
    hStack = root.THStack()
    # get all the histograms and set their plot style
    for hDef in hDefs:
        hm = hDef['hm']
        if 'den' in hDef:
            if 'denhm' in hDef:
                denhm = hDef['denhm']
            else:
                denhm = hm
        if threshold:
            h = hm.get_threshold_hist(hDef['num'], rebin=rebin).Clone()
            if den:
                hDen = denhm.get_threshold_hist(den, rebin=rebin)
                h.Divide(h, hDen, 1, 1, "b")
        else:
            if den:
                h = denhm.get_ratio(hDef['num'], den, rebin=rebin).Clone()
            else:
                h = hm.get(hDef['num'], rebin=rebin).Clone()

        if normalise:
            h.Scale(1/h.Integral())

        if normToBinWidth and not threshold and not den:
            for bin in range(1, h.GetNbinsX()+1):
               h.SetBinContent(bin, h.GetBinContent(bin) / h.GetBinWidth(bin))
               h.SetBinError(bin, h.GetBinError(bin) / h.GetBinWidth(bin))
        elif normToBinWidth:
            print 'Ignoring normToBinWidth flag for threshold or ratio plots'

        if scaleFactor != 1.:
            h.Scale(scaleFactor)

        # set plot style
        legStyle = ''
        if hDef['lc']:
            h.SetLineColor(hDef['lc'])
            legStyle += 'l'
        if hDef['ls']:
            h.SetLineStyle(hDef['ls'])
        if hDef['lw']:
            h.SetLineWidth(hDef['lw'])
        if hDef['fc']:
            h.SetFillColor(hDef['fc'])
            legStyle = 'f'
        if hDef['mc']:
            h.SetMarkerColor(hDef['mc'])
            h.SetMarkerSize(0.75)
            legStyle += 'p'
        if hDef['ms']:
            h.SetMarkerStyle(hDef['ms'])
        if hDef['stacked']:
            hStack.Add(h)
        if hDef['err']:
            legStyle += 'e'
        legEntries.append(legend.AddEntry(h, hDef['legtext'], legStyle))
        hs.append(h)

    # replace histograms to be stacked with stack histograms
    if hStack.GetNhists() > 0:
        canvas_name = name
        stackHistos = hStack.GetStack()
        j = len(stackHistos)-1
        for i, hDef in enumerate(hDefs):
            if hDef['stacked']:
                hs[i] = stackHistos[j].Clone()
                j -= 1
    else:
        canvas_name = name

    # find max and min bins > 0
    minBinValueNonZero = 1.
    maxBinValue = hs[0].GetBinContent(hs[0].GetMaximumBin())
    minBinValue = hs[0].GetBinContent(hs[0].GetMinimumBin())
    for h in hs[1:]:
        for b in range(1,h.GetNbinsX()+1):
            if h.GetBinContent(b) < minBinValueNonZero and h.GetBinContent(b) > 0.:
                minBinValueNonZero = h.GetBinContent(b)
        if h.GetBinContent(h.GetMaximumBin()) > maxBinValue:
            maxBinValue = h.GetBinContent(h.GetMaximumBin())

    if scaleFactor != 1.:
        canvas_name += '_scaled'

    # create canvas and draw on it
    canvas_title = canvas_name
    c = root.TCanvas(canvas_name, canvas_title, 100, 100, canvWidth, 600)
    c.cd()
    c.SetLogx(logx)
    c.SetLogy(logy)

    set_root_style()
    if legYmin < 0.6:
        root.gPad.SetRightMargin(0.2)

    # axis settings
    xAxis = hs[0].GetXaxis()
    yAxis = hs[0].GetYaxis()
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
    yMin = 0.
    if minBinValue < 0.:
        yMin = 1.2*minBinValue
    yMax = 1.2*maxBinValue
    if c.GetLogy():
        yMin = 0.5 * minBinValueNonZero
        yMax = 2. * maxBinValue
    yAxis.SetRangeUser(yMin, yMax)
    if xMax:
        xAxis.SetRangeUser(xAxis.GetBinLowEdge(1), xMax)

    # draw
    hs[0].Draw(hDefs[0]['drawopt'])
    for i, h in enumerate(hs[1:]):
        h.Draw(hDefs[i+1]['drawopt']+'same')
    hs[0].Draw('sameaxis')

    # draw vertical lines to mark TF boundaries
    lines = draw_lut_steps(hName=hDefs[0]['num'], xMin=h.GetXaxis().GetXmin(), yMin=h.GetYaxis().GetXmin(), xMax=h.GetXaxis().GetXmax(), yMax=h.GetYaxis().GetXmax())
    lines.append(draw_tf_eta_regions(hName=hDefs[0]['num'], xMin=h.GetXaxis().GetXmin(), yMin=h.GetYaxis().GetXmin(), xMax=h.GetXaxis().GetXmax(), yMax=h.GetYaxis().GetXmax()))

    legend.Draw('same')

    if canvWidth > 600:
        if data:
            #notes.append([0.48, 0.93, 'CMS preliminary, 13 TeV', False])
            notes.append([0.48, 0.93, 'CMS internal, 13 TeV', False])
        else:
            notes.append([0.484, 0.93, 'CMS Simulation, 13 TeV', False])
    else:
        if data:
            #notes.append([0.551, 0.93, 'CMS preliminary, 13 TeV', False])
            notes.append([0.551, 0.93, 'CMS internal, 13 TeV', False])
        else:
            notes.append([0.555, 0.93, 'CMS Simulation, 13 TeV', False])
    tex = add_text(notes, True, False)

    c.Modified()
    c.Update()

    return [c, hs, legend, lines, tex]

# draw vertical lines to mark TF boundaries
def draw_tf_eta_regions(hName='', xMin=0., yMin=0., xMax=1., yMax=1., twoD=False, drawDiag=False):
    lines = []
    if hName[-4:] == '.eta':
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

# draw steps to indicate LUT resolution for extrapolation
def draw_lut_steps(hName='', xMin=0., yMin=0., xMax=1., yMax=1.):
    lines = []
    step = 0.
    if hName[-5:] == '_deta' or hName[-8:] == '_absdeta':
       step = 0.010875
    if hName[-5:] == '_dphi' or hName[-8:] == '_absdphi':
       step = 0.0435
    if hName[-5:] == '_deta' or hName[-8:] == '_absdeta' or hName[-5:] == '_dphi' or hName[-8:] == '_absdphi':
        for i in range(1, 16):
            yCoord = i * step
            lines.append(root.TLine(0., yCoord, 63., yCoord))
            lines[-1].SetLineStyle(root.kDashed)
            lines[-1].Draw('same')
        lines.append(root.TLine(63., yMin, 63., yMax))
        lines[-1].SetLineStyle(root.kDashed)
        lines[-1].Draw('same')

    return lines
    

def define_styles():
    styles = {}
    styles['gmt'] = {'lc':root.kCyan, 'ls':root.kSolid, 'fc':root.kCyan, 'mc':root.kCyan, 'ms':root.kOpenCircle, 'legtext':'GMT'}
    styles['ugmt'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':root.kBlack, 'mc':root.kBlack, 'ms':root.kOpenCircle, 'legtext':'uGMT'}
    styles['data'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':root.kBlack, 'mc':root.kBlack, 'ms':root.kFullCircle, 'legtext':'data'}
    styles['etarange_0'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':root.kBlack, 'mc':root.kBlack, 'ms':root.kOpenTriangleDown, 'legtext':'er0'}
    styles['etarange_1'] = {'lc':root.kRed, 'ls':root.kSolid, 'fc':root.kRed, 'mc':root.kRed, 'ms':root.kOpenTriangleUp, 'legtext':'er1'}
    styles['etarange_2'] = {'lc':root.kBlue, 'ls':root.kSolid, 'fc':root.kBlue, 'mc':root.kBlue, 'ms':root.kOpenCircle, 'legtext':'er2'}
    styles['etarange_3'] = {'lc':root.kGreen, 'ls':root.kSolid, 'fc':root.kGreen, 'mc':root.kGreen, 'ms':root.kFullSquare, 'legtext':'er3'}
    styles['etarange_4'] = {'lc':root.kGreen+2, 'ls':root.kSolid, 'fc':root.kGreen+2, 'mc':root.kGreen+2, 'ms':root.kFullCircle, 'legtext':'er4'}
    styles['etarange_5'] = {'lc':root.kGreen+3, 'ls':root.kSolid, 'fc':root.kGreen+3, 'mc':root.kGreen+3, 'ms':root.kFullTriangleDown, 'legtext':'er5'}
    styles['etarange_6'] = {'lc':root.kGreen+4, 'ls':root.kSolid, 'fc':root.kGreen+4, 'mc':root.kGreen+4, 'ms':root.kFullTriangleUp, 'legtext':'er6'}
    styles['data_q12'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':root.kViolet, 'mc':root.kViolet, 'ms':root.kOpenTriangleUp, 'legtext':'data Q #geq 12'}
    styles['data_q8'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':root.kCyan+3, 'mc':root.kCyan+3, 'ms':root.kOpenCircle, 'legtext':'data Q #geq 8'}
    styles['data_q4'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':root.kYellow, 'mc':root.kYellow, 'ms':root.kFullSquare, 'legtext':'data Q #geq 4'}
    styles['data_qmin12'] = {'lc':root.kViolet, 'ls':root.kSolid, 'fc':root.kViolet, 'mc':root.kViolet, 'ms':root.kOpenTriangleUp, 'legtext':'data Q #geq 12'}
    styles['data_qmin8'] = {'lc':root.kCyan+3, 'ls':root.kSolid, 'fc':root.kCyan+3, 'mc':root.kCyan+3, 'ms':root.kOpenCircle, 'legtext':'data Q #geq 8'}
    styles['data_qmin4'] = {'lc':root.kOrange-2, 'ls':root.kSolid, 'fc':root.kYellow, 'mc':root.kYellow, 'ms':root.kFullSquare, 'legtext':'data Q #geq 4'}
    styles['emul'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':root.kBlue, 'mc':root.kBlue, 'ms':root.kFullCircle, 'legtext':'emul'}
    styles['emul_q12'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':root.kViolet, 'mc':root.kViolet, 'ms':root.kOpenTriangleUp, 'legtext':'emul Q #geq 12'}
    styles['emul_q8'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':root.kCyan+3, 'mc':root.kCyan+3, 'ms':root.kOpenCircle, 'legtext':'emul Q #geq 8'}
    styles['emul_q4'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':root.kYellow, 'mc':root.kYellow, 'ms':root.kFullSquare, 'legtext':'emul Q #geq 4'}
    styles['mu_1'] = {'lc':root.kRed, 'ls':root.kSolid, 'fc':root.kRed, 'mc':root.kRed, 'ms':root.kOpenTriangleUp, 'legtext':'1^{st} muon'}
    styles['mu_2'] = {'lc':root.kBlue, 'ls':root.kSolid, 'fc':root.kBlue, 'mc':root.kBlue, 'ms':root.kOpenCircle, 'legtext':'2^{nd} muon'}
    styles['mu_3'] = {'lc':root.kGreen, 'ls':root.kSolid, 'fc':root.kGreen, 'mc':root.kGreen, 'ms':root.kFullSquare, 'legtext':'3^{rd} muon'}
    return styles

def hist_style(key, filled=False, marker=False, lw=1):
    styles = define_styles()
    style = styles[key]
    if not filled:
        style['fc'] = None
    if not marker:
        style['mc'] = None
        style['ms'] = None
    style['lw'] = lw

    return style

