#!/usr/bin/env python
from ToolBox import parse_options_and_init_log
# have to do this first or ROOT masks the -h messages
opts, parser = parse_options_and_init_log()

from L1Analysis import L1Ana, L1Ntuple
from analysis_tools.plotting import HistManager, HistManager2d
from analysis_tools.plottools import define_styles, hist_style
from analysis_tools.selections import MuonSelections, Matcher
import ROOT as root
import re
import os

def parse_options_plotRates(parser):
    """
    Adds often used options to the OptionParser...
    """
    parsers = parser.add_subparsers()
    sub_parser = parsers.add_parser("plotIso")
    sub_parser.add_argument("-i", "--interactive", dest="interactive", action='store_false', help="Draw plots on screen.")
    sub_parser.add_argument("--towers", dest="towers", default=False, action='store_true', help="Plots for calo towers")
    sub_parser.add_argument("--2d", dest="twod", default=False, action='store_true', help="2D plots.")
    sub_parser.add_argument("--emul", dest="emul", default=False, action="store_true", help="Make emulator plots.")
    sub_parser.add_argument("--public", dest="public", default=False, action='store_true', help="Plot style for publication.")

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


def extract_notes_from_name(name, xBase, yBase, etaTxt=True, qualTxt=True, ptTxt=True):
    notes = []
    # extract eta range
    if etaTxt and name.rfind('.eta') == -1:
        eta_number_strs = re.findall(r'[\d\.\d]+', name[name.find('EtaMin')+6:name.find('EtaMax')+12])
        if len(eta_number_strs) > 1:
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
        ptminPos = name.find('l1ptmin')
        if ptminPos >= 0:
            l1_ptmin_strs = re.findall(r'\d+\.?\d*', name[ptminPos+5:ptminPos+9])
            if len(l1_ptmin_strs) > 0:
                if l1_ptmin_strs[0][-1] == '.':
                    l1_ptmin_strs[0] = l1_ptmin_strs[0][0:-1]
                notes.append([xBase, yBase+0.05, 'p_{T}^{L1} #geq '+l1_ptmin_strs[0]+' GeV', True])
    return notes

def plot_2dhist(hm2d, hName, drawDiag=True, data=False, xMin=None, xMax=None, yMin=None, yMax=None, scaleFactor=1.):
    canvas_name = hName
    if xMin:
        canvas_name += '_xMin{min}'.format(min=xMin)
    if xMax:
        canvas_name += '_xMax{max}'.format(max=xMax)
    if yMin:
        canvas_name += '_yMin{min}'.format(min=yMin)
    if yMax:
        canvas_name += '_yMax{max}'.format(max=yMax)

    # create canvas and draw on it
    canvas_title = canvas_name
    c = root.TCanvas(canvas_name, canvas_title, 100, 100, 600, 600)
    c.cd()
    set_root_style()
    root.gPad.SetRightMargin(0.14)

    if hName not in hm2d.get_varnames():
        return [c]
    h = hm2d.get(hName).Clone()

    if scaleFactor != 1.:
        h.Scale(scaleFactor)

    xAxis = h.GetXaxis()
    yAxis = h.GetYaxis()
    zAxis = h.GetZaxis()
    xAxis.SetTitleFont(font)
    xAxis.SetLabelFont(font)
    xAxis.SetLabelSize(fontSize)
    xAxis.SetNoExponent()
    if xMin and not xMax:
        xAxis.SetRangeUser(xMin, xAxis.GetBinUpEdge(xAxis.GetNbins()))
    elif xMax and not xMin:
        xAxis.SetRangeUser(xAxis.GetBinLowEdge(1), xMax)
    elif xMin and xMax:
        xAxis.SetRangeUser(xMin, xMax)
    yAxis.SetTitleOffset(1.5)
    yAxis.SetTitleFont(font)
    yAxis.SetLabelFont(font)
    yAxis.SetLabelSize(fontSize)
    if yMin and not yMax:
        yAxis.SetRangeUser(yMin, yAxis.GetBinUpEdge(yAxis.GetNbins()))
    elif yMax and not yMin:
        yAxis.SetRangeUser(yAxis.GetBinLowEdge(1), yMax)
    elif yMin and yMax:
        yAxis.SetRangeUser(yMin, yMax)
    zAxis.SetTitleFont(font)
    zAxis.SetLabelFont(font)
    zAxis.SetLabelSize(fontSize)
    h.Draw('colz')

    xBase = 0.56
    yBase = 0.13
    notes = extract_notes_from_name(hName, xBase, yBase, etaTxt=False, qualTxt=False, ptTxt=True)

    if data:
        notes.append([0.53, 0.93, 'CMS internal, 13 TeV', False])
    else:
        notes.append([0.53, 0.93, 'CMS Simulation, 13 TeV', False])
    tex = add_text(notes=notes)

    lines = draw_tf_eta_regions(hName=hName, xMin=h.GetXaxis().GetXmin(), yMin=h.GetYaxis().GetXmin(), xMax=h.GetXaxis().GetXmax(), yMax=h.GetYaxis().GetXmax(), twoD=True, drawDiag=drawDiag)

    c.Modified()
    c.Update()
    
    return [c, h, tex, lines]


def plot_hists(hm, hDefs, xTitle=None, yTitle='# muons', threshold=False, normToBinWidth=False, normalise=False, xMin=None, xMax=None, canvasPrefix='', notes=None, scaleFactor=1., logx=False, logy=False, data=False):
    den = hDefs[0]['den']
    if den:
        name = canvasPrefix+hDefs[0]['num']+'_over_'+den
    else:
        name = canvasPrefix+hDefs[0]['num']
    if normToBinWidth and not threshold and not den:
        name += '_normToBinWidth'
    if xMin:
        name += '_xMin{min}'.format(min=xMin)
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
        if threshold:
            h = hm.get_threshold_hist(hDef['num']).Clone()
            if den:
                hDen = hm.get_threshold_hist(den)
                h.Divide(h, hDen, 1, 1, "b")
        else:
            if den:
                h = hm.get_ratio(hDef['num'], den).Clone()
            else:
                h = hm.get(hDef['num']).Clone()

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
    yMax = 1.2*maxBinValue
    if c.GetLogy():
        yMin = 0.5 * minBinValueNonZero
        yMax = 2. * maxBinValue
    yAxis.SetRangeUser(yMin, yMax)
    if xMin and not xMax:
        xAxis.SetRangeUser(xMin, xAxis.GetBinUpEdge(xAxis.GetNbins()))
    elif xMax and not xMin:
        xAxis.SetRangeUser(xAxis.GetBinLowEdge(1), xMax)
    elif xMin and xMax:
        xAxis.SetRangeUser(xMin, xMax)

    # draw
    hs[0].Draw(hDefs[0]['drawopt'])
    for i, h in enumerate(hs[1:]):
        h.Draw(hDefs[i+1]['drawopt']+'same')
    hs[0].Draw('sameaxis')

    # draw vertical lines to mark TF boundaries
    lines = draw_tf_eta_regions(hName=hDefs[0]['num'], xMin=h.GetXaxis().GetXmin(), yMin=h.GetYaxis().GetXmin(), xMax=h.GetXaxis().GetXmax(), yMax=h.GetYaxis().GetXmax())

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
    

def plot_hist_standard(hm, hName, den=None, xTitle='', yTitle='# muons', threshold=False, stacked=False, normToBinWidth=False, xMin=None, xMax=None, reg='', scaleFactor=1., data=False):
    ugmt_dict = {'num':hName, 'den':den, 'drawopt':'', 'stacked':False, 'err':True}
    ugmt_dict.update(hist_style('data', lw=2))
    hDefs = []
    if reg == '':
        hDefs.append(ugmt_dict)

    logy = False
    if hName[-4:] == '.iet' and not den:
        logy = True

    xBase = 0.17
    yBase = 0.56
    notes = extract_notes_from_name(hName, xBase, yBase, etaTxt=False, qualTxt=False, ptTxt=True)
    if den:
        den_eta_number_strs = re.findall(r'[\d\.\d]+', den[den.find('EtaMin')+6:den.find('EtaMax')+12])
        if len(den_eta_number_strs) > 1 and eta_number_strs != den_eta_number_strs:
            den_note_str = den_eta_number_strs[0]+' < |#eta^{GMT}| < '+den_eta_number_strs[1]
            notes.append([xBase, yBase-0.05, den_note_str, True])

    prefix = 'c_'

    return plot_hists(hm, hDefs, xTitle, yTitle, threshold, normToBinWidth, False, xMin, xMax, prefix, notes, scaleFactor, False, logy, data)

def plot_hists_standard(hm, hNames, den=None, xTitle='', yTitle='# muons', legtexts=None, threshold=False, stacked=False, normToBinWidth=False, xMin=None, xMax=None, reg='', scaleFactor=1., data=False):
    dicts = []
    hDefs = []
    for i, hName in enumerate(hNames):
        dicts.append({'num':hName, 'den':den, 'drawopt':'', 'stacked':False, 'err':True})
        if legtexts:
            dicts[-1].update(hist_style('generic_{i}'.format(i=i), legtext=legtexts[i], lw=2))
        else:
            dicts[-1].update(hist_style('generic_{i}'.format(i=i), lw=2))
        if reg == '':
            hDefs.append(dicts[-1])

    logy = False
    #if (hName[-3:] == 'iet' or hName.find('.n') != -1) and not den:
    #    logy = True

    xBase = 0.17
    yBase = 0.56
    notes = extract_notes_from_name(hName, xBase, yBase, etaTxt=False, qualTxt=False, ptTxt=True)
    if den:
        den_eta_number_strs = re.findall(r'[\d\.\d]+', den[den.find('EtaMin')+6:den.find('EtaMax')+12])
        if len(den_eta_number_strs) > 1 and eta_number_strs != den_eta_number_strs:
            den_note_str = den_eta_number_strs[0]+' < |#eta^{GMT}| < '+den_eta_number_strs[1]
            notes.append([xBase, yBase-0.05, den_note_str, True])

    prefix = 'c_'

    return plot_hists(hm, hDefs, xTitle, yTitle, threshold, normToBinWidth, False, xMin, xMax, prefix, notes, scaleFactor, False, logy, data)

def main():
    opts = parse_options_plotRates(parser)
    batchRun = opts.interactive
    if batchRun:
        root.gROOT.SetBatch(True)

    global styles
    styles = define_styles()

    #l1PtMins = [0]
    l1PtMins = [0, 5, 12, 20]
    #l1PtMins = [0, 5, 12, 20, 25]

    namePrefix = ''
    if opts.emul:
        namePrefix += 'emu_'

    isData = True

    hm = HistManager(filename=opts.fname, subdir='all_runs')

    nEvents = hm.get(namePrefix+'l1_caloTower.n').GetEntries()
    print 'Found {n} processed events.'.format(n=nEvents)

    L1Ana.init_l1_analysis()
    print ""

    # holds the canvases, histograms, etc.
    objects = []

    ##########################################################################
    # L1 calo towers variables
    if opts.towers:
        objects.append(plot_hist_standard(hm, namePrefix+'l1_caloTower.n', data=isData))
        objects.append(plot_hist_standard(hm, namePrefix+'l1_caloTower.iet', xTitle='iE_{T}', yTitle='# towers', data=isData))
        objects.append(plot_hist_standard(hm, namePrefix+'l1_caloTower.ieta', xTitle='i#eta', data=isData))
        objects.append(plot_hist_standard(hm, namePrefix+'l1_caloTower.iphi', xTitle='i#phi', data=isData))
        objects.append(plot_hist_standard(hm, namePrefix+'l1_caloTower.iqual', xTitle='i qual', data=isData))
        for l1PtMin in l1PtMins:
            ptMinStr = '_l1ptmin{ptmin}'.format(ptmin=l1PtMin)

            hPrefix = namePrefix+'l1_caloTower'+ptMinStr
            objects.append(plot_hists_standard(hm, [hPrefix+'.area_1x1_iet', hPrefix+'.area_11x11-1x1_iet', hPrefix+'.area_11x11_iet'], xTitle='iE_{T}', legtexts=['1x1', '11x11-1x1', '11x11'], data=isData, xMax=50))
            objects.append(plot_hists_standard(hm, [hPrefix+'.area_1x3_iet', hPrefix+'.area_11x11-1x3_iet', hPrefix+'.area_11x11_iet'], xTitle='iE_{T}', legtexts=['1x3', '11x11-1x3', '11x11'], data=isData, xMax=50))
            objects.append(plot_hists_standard(hm, [hPrefix+'.area_1x5_iet', hPrefix+'.area_11x11-1x5_iet', hPrefix+'.area_11x11_iet'], xTitle='iE_{T}', legtexts=['1x5', '11x11-1x5', '11x11'], data=isData, xMax=50))
            objects.append(plot_hists_standard(hm, [hPrefix+'.area_3x3_iet', hPrefix+'.area_11x11-3x3_iet', hPrefix+'.area_11x11_iet'], xTitle='iE_{T}', legtexts=['3x3', '11x11-3x3', '11x11'], data=isData, xMax=50))
            objects.append(plot_hists_standard(hm, [hPrefix+'.area_1xm2to0_iet', hPrefix+'.area_1x0top2_iet'], xTitle='iE_{T}', legtexts=['1x-2to0', '1x0to+2'], data=isData, xMax=25))
            objects.append(plot_hists_standard(hm, [hPrefix+'.area_3xm3to0_iet', hPrefix+'.area_3x0top3_iet'], xTitle='iE_{T}', legtexts=['3x-3to0', '3x0to+3'], data=isData, xMax=50))
            objects.append(plot_hists_standard(hm, [hPrefix+'.area_3xm7to0_iet', hPrefix+'.area_3x0top7_iet'], xTitle='iE_{T}', legtexts=['3x-7to0', '3x0to+7'], data=isData, xMax=50))
            objects.append(plot_hists_standard(hm, [hPrefix+'.area_1xm2to0_iet_minus_area_1x0top2_iet',
                                                    hPrefix+'.area_1xm2to0_iet_minus_area_1x0top2_iet_mu_chg_pos',
                                                    hPrefix+'.area_1xm2to0_iet_minus_area_1x0top2_iet_mu_chg_neg'], xTitle='iE_{T}',
                                               legtexts=['1x-2to0 - 1x0to+2', '#mu^{+} 1x-2to0 - 1x0to+2', '#mu^{-} 1x-2to0 - 1x0to+2'], data=isData, xMin=-15, xMax=15))
            objects.append(plot_hists_standard(hm, [hPrefix+'.area_3xm3to0_iet_minus_area_3x0top3_iet',
                                                    hPrefix+'.area_3xm3to0_iet_minus_area_3x0top3_iet_mu_chg_pos',
                                                    hPrefix+'.area_3xm3to0_iet_minus_area_3x0top3_iet_mu_chg_neg'], xTitle='iE_{T}',
                                               legtexts=['3x-3to0 - 3x0to+3', '#mu^{+} 3x-3to0 - 3x0to+3', '#mu^{-} 3x-3to0 - 3x0to+3'], data=isData, xMin=-15, xMax=15))
            objects.append(plot_hists_standard(hm, [hPrefix+'.area_3xm7to0_iet_minus_area_3x0top7_iet',
                                                    hPrefix+'.area_3xm7to0_iet_minus_area_3x0top7_iet_mu_chg_pos',
                                                    hPrefix+'.area_3xm7to0_iet_minus_area_3x0top7_iet_mu_chg_neg'], xTitle='iE_{T}',
                                               legtexts=['3x-7to0 - 3x0to+7', '#mu^{+} 3x-7to0 - 3x0to+7', '#mu^{-} 3x-7to0 - 3x0to+7'], data=isData, xMin=-15, xMax=15))
            objects.append(plot_hists_standard(hm, [hPrefix+'.area_3xm3to0_iet_over_area_3x0top3_iet',
                                                    hPrefix+'.area_3xm3to0_iet_over_area_3x0top3_iet_mu_chg_pos',
                                                    hPrefix+'.area_3xm3to0_iet_over_area_3x0top3_iet_mu_chg_neg'],
                                               legtexts=['3x-3to0 / 3x0to+3', '#mu^{+} 3x-3to0 / 3x0to+3', '#mu^{-} 3x-3to0 / 3x0to+3'], data=isData))
            objects.append(plot_hists_standard(hm, [hPrefix+'.area_3xm3to0_iet_minus_area_3x0top3_iet_over_area_3xm3to0_iet_plus_area_3x0top3_iet',
                                                    hPrefix+'.area_3xm3to0_iet_minus_area_3x0top3_iet_over_area_3xm3to0_iet_plus_area_3x0top3_iet_mu_chg_pos',
                                                    hPrefix+'.area_3xm3to0_iet_minus_area_3x0top3_iet_over_area_3xm3to0_iet_plus_area_3x0top3_iet_mu_chg_neg'],
                                               legtexts=['(3x-3to0 - 3x0to+3)/(3x-3to0 + 3x0to+3)', '#mu^{+} (3x-3to0 - 3x0to+3)/(3x-3to0 + 3x0to+3)', '#mu^{-} (3x-3to0 - 3x0to+3)/(3x-3to0 + 3x0to+3)'], data=isData))
            objects.append(plot_hists_standard(hm, [hPrefix+'.area_3xm7to0_iet_minus_area_3x0top7_iet_over_area_3xm7to0_iet_plus_area_3x0top7_iet',
                                                    hPrefix+'.area_3xm7to0_iet_minus_area_3x0top7_iet_over_area_3xm7to0_iet_plus_area_3x0top7_iet_mu_chg_pos',
                                                    hPrefix+'.area_3xm7to0_iet_minus_area_3x0top7_iet_over_area_3xm7to0_iet_plus_area_3x0top7_iet_mu_chg_neg'],
                                               legtexts=['(3x-7to0 - 3x0to+7)/(3x-7to0 + 3x0to+7)', '#mu^{+} (3x-7to0 - 3x0to+7)/(3x-7to0 + 3x0to+7)', '#mu^{-} (3x-7to0 - 3x0to+7)/(3x-7to0 + 3x0to+7)'], data=isData))
            objects.append(plot_hist_standard(hm, hPrefix+'.area_3xm3to0_iet_minus_area_3x0top3_iet_over_area_3xm3to0_iet_plus_area_3x0top3_iet_times_mu_chg', xTitle='#mu charge * (3x-3to0 - 3x0to+3)/(3x-3to0 + 3x0to+3)', data=isData))
            objects.append(plot_hist_standard(hm, hPrefix+'.area_3xm7to0_iet_minus_area_3x0top7_iet_over_area_3xm7to0_iet_plus_area_3x0top7_iet_times_mu_chg', xTitle='#mu charge * (3x-7to0 - 3x0to+7)/(3x-7to0 + 3x0to+7)', data=isData))
            objects.append(plot_hist_standard(hm, hPrefix+'.area_11x11-1x1_over_area_11x11_iet', xTitle='iE_{T}^{11x11-1x1} / iE_{T}^{11x11}', data=isData))
            objects.append(plot_hist_standard(hm, hPrefix+'.area_11x11-1x3_over_area_11x11_iet', xTitle='iE_{T}^{11x11-1x3} / iE_{T}^{11x11}', data=isData))
            objects.append(plot_hist_standard(hm, hPrefix+'.area_11x11-1x5_over_area_11x11_iet', xTitle='iE_{T}^{11x11-1x5} / iE_{T}^{11x11}', data=isData))
            objects.append(plot_hist_standard(hm, hPrefix+'.area_11x11-3x3_over_area_11x11_iet', xTitle='iE_{T}^{11x11-3x3} / iE_{T}^{11x11}', data=isData))

            objects.append(plot_hists_standard(hm, [hPrefix+'.twobytwo_area_1x1_iet', hPrefix+'.twobytwo_area_5x5-1x1_iet', hPrefix+'.twobytwo_area_5x5_iet'], xTitle='iE_{T}^{2x2}', legtexts=['1x1', '5x5-1x1', '5x5'], data=isData, xMax=50))
            objects.append(plot_hists_standard(hm, [hPrefix+'.twobytwo_area_1xm1to0_iet', hPrefix+'.twobytwo_area_1x0top1_iet'], xTitle='2x2 area iE_{T}', legtexts=['1x-1to0', '1x0to+1'], data=isData, xMax=50))
            objects.append(plot_hists_standard(hm, [hPrefix+'.twobytwo_area_1xm1to0_iet_minus_twobytwo_area_1x0top1_iet',
                                                    hPrefix+'.twobytwo_area_1xm1to0_iet_minus_twobytwo_area_1x0top1_iet_mu_chg_pos',
                                                    hPrefix+'.twobytwo_area_1xm1to0_iet_minus_twobytwo_area_1x0top1_iet_mu_chg_neg'], xTitle='2x2 area iE_{T}',
                                               legtexts=['1x-1to0 - 1x0to+1', '#mu^{+} 1x-1to0 - 1x0to+1', '#mu^{-} 1x-1to0 - 1x0to+1'], data=isData, xMin=-15, xMax=15))
            objects.append(plot_hist_standard(hm, hPrefix+'.twobytwo_area_5x5-1x1_over_twobytwo_area_5x5_iet', xTitle='2x2 area iE_{T}^{5x5-1x1} / iE_{T}^{5x5}', data=isData))
            objects.append(plot_hist_standard(hm, hPrefix+'.twobytwo_area_5x5_iet_over_mu_ipt', xTitle='2x2 area iE_{T}^{5x5} / ip_{T}^{#mu}', data=isData))
            objects.append(plot_hist_standard(hm, hPrefix+'.twobytwo_area_5x5-1x1_iet_over_mu_ipt', xTitle='2x2 area iE_{T}^{5x5-1x1} / ip_{T}^{#mu}', data=isData))
            objects.append(plot_hist_standard(hm, hPrefix+'.twobytwo_area_5x5_5bit_iet_over_mu_ipt', xTitle='5 bit (iE_{T}^{max}=31) 2x2 area iE_{T}^{5x5} / ip_{T}^{#mu}', data=isData))
            objects.append(plot_hist_standard(hm, hPrefix+'.twobytwo_area_5x5-1x1_5bit_iet_over_mu_ipt', xTitle='5 bit (iE_{T}^{max}=31) 2x2 area iE_{T}^{5x5-1x1} / ip_{T}^{#mu}', data=isData))

    # 2d plots
    if opts.twod:
        hm2d = HistManager2d(filename=opts.fname, subdir='all_runs')

        histoprefix2d = namePrefix+'2d_caloTower'
        objects.append(plot_2dhist(hm2d, histoprefix2d+'.ieta_iphi', drawDiag=False, data=isData))
        objects.append(plot_2dhist(hm2d, histoprefix2d+'.iet_ieta_iet_iphi', drawDiag=False, data=isData, scaleFactor=1/nEvents))
        for l1PtMin in l1PtMins:
            ptMinStr = '_l1ptmin{ptmin}'.format(ptmin=l1PtMin)
            nMuEvents = hm.get(namePrefix+'l1_caloTower'+ptMinStr+'.n_mu').Integral(2, hm.get(namePrefix+'l1_caloTower'+ptMinStr+'.n_mu').GetNbinsX())
            print 'Found {nmu} L1 muons passing {cut} GeV cut in this event.'.format(nmu=nMuEvents, cut=l1PtMin)

            objects.append(plot_2dhist(hm2d, histoprefix2d+ptMinStr+'.iet_ietarel_iet_iphirel', drawDiag=False, data=isData, scaleFactor=1/nMuEvents))
            objects.append(plot_2dhist(hm2d, histoprefix2d+ptMinStr+'.iet_ietarel_red_iet_iphirel_red', drawDiag=False, data=isData, scaleFactor=1/nMuEvents))
            objects.append(plot_2dhist(hm2d, histoprefix2d+ptMinStr+'.area_11x11_iet_area_1x1_iet', drawDiag=False, data=isData, xMax=50, yMax=50))
            objects.append(plot_2dhist(hm2d, histoprefix2d+ptMinStr+'.area_11x11_iet_area_1x3_iet', drawDiag=False, data=isData, xMax=50, yMax=50))
            objects.append(plot_2dhist(hm2d, histoprefix2d+ptMinStr+'.area_11x11_iet_area_1x5_iet', drawDiag=False, data=isData, xMax=50, yMax=50))
            objects.append(plot_2dhist(hm2d, histoprefix2d+ptMinStr+'.area_11x11_iet_area_3x3_iet', drawDiag=False, data=isData, xMax=50, yMax=50))
            objects.append(plot_2dhist(hm2d, histoprefix2d+ptMinStr+'.area_11x11-1x1_iet_area_1x1_iet', drawDiag=False, data=isData, xMax=50, yMax=50))
            objects.append(plot_2dhist(hm2d, histoprefix2d+ptMinStr+'.area_11x11-1x3_iet_area_1x3_iet', drawDiag=False, data=isData, xMax=50, yMax=50))
            objects.append(plot_2dhist(hm2d, histoprefix2d+ptMinStr+'.area_11x11-1x5_iet_area_1x5_iet', drawDiag=False, data=isData, xMax=50, yMax=50))
            objects.append(plot_2dhist(hm2d, histoprefix2d+ptMinStr+'.area_11x11-3x3_iet_area_3x3_iet', drawDiag=False, data=isData, xMax=50, yMax=50))

            if l1PtMin == 0:
                objects.append(plot_2dhist(hm2d, histoprefix2d+ptMinStr+'.mu_pt_area_1x1_iet', drawDiag=False, data=isData, xMax=50, yMax=50))
                objects.append(plot_2dhist(hm2d, histoprefix2d+ptMinStr+'.mu_pt_area_1x3_iet', drawDiag=False, data=isData, xMax=50, yMax=50))
                objects.append(plot_2dhist(hm2d, histoprefix2d+ptMinStr+'.mu_pt_area_1x5_iet', drawDiag=False, data=isData, xMax=50, yMax=50))
                objects.append(plot_2dhist(hm2d, histoprefix2d+ptMinStr+'.mu_pt_area_3x3_iet', drawDiag=False, data=isData, xMax=50, yMax=50))
                objects.append(plot_2dhist(hm2d, histoprefix2d+ptMinStr+'.mu_pt_area_1xm2to0_iet_minus_area_1x0top2_iet', drawDiag=False, data=isData, xMax=50, yMin=-15, yMax=15))
                objects.append(plot_2dhist(hm2d, histoprefix2d+ptMinStr+'.mu_pt_area_1xm2to0_iet_minus_area_1x0top2_iet_mu_chg_pos', drawDiag=False, data=isData, xMax=50, yMin=-15, yMax=15))
                objects.append(plot_2dhist(hm2d, histoprefix2d+ptMinStr+'.mu_pt_area_1xm2to0_iet_minus_area_1x0top2_iet_mu_chg_neg', drawDiag=False, data=isData, xMax=50, yMin=-15, yMax=15))
                objects.append(plot_2dhist(hm2d, histoprefix2d+ptMinStr+'.mu_pt_twobytwo_area_1x1_iet', drawDiag=False, data=isData, xMax=50, yMax=50))
                objects.append(plot_2dhist(hm2d, histoprefix2d+ptMinStr+'.mu_pt_twobytwo_area_5x5_iet', drawDiag=False, data=isData, xMax=50, yMax=50))
                objects.append(plot_2dhist(hm2d, histoprefix2d+ptMinStr+'.mu_pt_twobytwo_area_5x5-1x1_iet', drawDiag=False, data=isData, xMax=50, yMax=50))

    ##########################################################################
    # save plots to root file
    if savePlots:
        plotdir = 'plots_'+opts.fname.replace('.root','').partition('/')[0]
        if opts.public:
            plotdir += '_public'
        if not os.path.exists(plotdir):
            os.makedirs(plotdir)
        output = root.TFile('./'+plotdir+'/l1_muon_iso_plots.root', 'recreate')
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
    styles = {}
    main()

