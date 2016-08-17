#!/usr/bin/env python
from ToolBox import parse_options_and_init_log
# have to do this first or ROOT masks the -h messages
opts, parser = parse_options_and_init_log()

from L1Analysis import L1Ana, L1Ntuple
from analysis_tools.plotting import HistManager
from analysis_tools.selections import MuonSelections, Matcher
import ROOT as root
import re

def parse_options_plotEff(parser):
    """
    Adds often used options to the OptionParser...
    """
    parsers = parser.add_subparsers()
    sub_parser = parsers.add_parser("plotEff")
    sub_parser.add_argument("-i", "--interactive", dest="interactive", action='store_false', help="Draw plots on screen.")

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

def plot_histos(hm):
    '''
    draw all histograms in hm
    cList: list with all canvases
    '''
    for varname in hm.get_varnames():
        c = root.TCanvas(varname, varname, 100, 100, 600, 600)
        c.cd()

        set_root_style()

        histo = hm.get(varname)
        histo.GetYaxis().SetTitleOffset(1.5)
        histo.SetLineWidth(2)
        histo.SetLineColor(root.kBlue)
        histo.Draw('hist')

        tex = root.TLatex()
        tex.SetNDC()
        tex.SetTextFont(font)
        tex.SetTextSize(0.04)
        tex.DrawLatex(0.555, 0.93, 'Simulation, 13 TeV')
        #tex.DrawLatex(0.555, 0.93, 'CMS Simulation, 13 TeV')

    return [c, histo, tex]

def plot_hists(hm, hDefs, xTitle=None, yTitle='# muons', normToBinWidth=False, canvasPrefix='', notes=None):
    if hDefs[0]['den']:
        name = canvasPrefix+hDefs[0]['num']+'_over_'+hDefs[0]['den']
    else:
        name = canvasPrefix+hDefs[0]['num']
    if normToBinWidth:
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
        if hDef['den']:
            h = hm.get_ratio(hDef['num'], hDef['den']).Clone()
        else:
            h = hm.get(hDef['num']).Clone()

        if normToBinWidth:
            for bin in range(1, h.GetNbinsX()+1):
               h.SetBinContent(bin, h.GetBinContent(bin) / h.GetBinWidth(bin))
               h.SetBinError(bin, h.GetBinError(bin) / h.GetBinWidth(bin))
        elif normToBinWidth:
            print 'Ignoring normToBinWidth flag for ratio plots'

        h.SetLineColor(hDef['lc'])
        h.SetLineStyle(hDef['ls'])
        h.SetLineWidth(2)
        legStyle = 'l'
        if hDef['fc']:
            h.SetFillColor(hDef['fc'])
            h.SetLineWidth(1)
            legStyle = 'f'
            # if a fill colour is defined stack this histogram with others
            hStack.Add(h)
        legEntries.append(legend.AddEntry(h, hDef['legtext'], legStyle))
        if len(hDef['legtext']) > 15:
            legend.SetX1(0.45)
            legend.SetX2(0.75)
        hs.append(h)

    # replace histograms to be stacked with stack histograms
    if hStack.GetNhists() > 0:
        canvas_name = 'c_eff_stacked_'+name
        stackHistos = hStack.GetStack()
        j = len(stackHistos)-1
        for i, hDef in enumerate(hDefs):
            if hDef['fc']:
                hs[i] = stackHistos[j].Clone()
                j -= 1
    else:
        canvas_name = 'c_eff_'+name

    # create canvas and draw on it
    canvas_title = canvas_name
    c = root.TCanvas(canvas_name, canvas_title, 100, 100, canvWidth, 600)
    c.cd()
    if name[-2:] == 'pt' and not hDefs[0]['den']:
        c.SetLogy(True)

    set_root_style()
    if legYmin < 0.6:
        root.gPad.SetRightMargin(0.2)

    if xTitle:
        hs[0].GetXaxis().SetTitle(xTitle)
    hs[0].GetYaxis().SetTitleOffset(1.5)
    hs[0].GetYaxis().SetTitle(yTitle)
    maxBinValue = hs[0].GetBinContent(hs[0].GetMaximumBin())
    if not c.GetLogy():
        yMax = 1.2*maxBinValue
        if maxBinValue <= 1.:
            yMax = 1.4
        hs[0].GetYaxis().SetRangeUser(0., yMax)
    # draw
    hs[0].SetLineWidth(2)
    legEntries[0].SetObject(hs[0])
    legEntries[0].SetOption(legEntries[0].GetOption()+'le')
    hs[0].Draw('hist')
    for h in hs[1:]:
        h.Draw('histsame')
    hs[0].Draw('same')
    hs[0].Draw('sameaxis')

    # draw vertical lines to mark TF boundaries
    lines = []
    if name[-3:] == 'eta':
        lines.append(root.TLine(0.83, 0., 0.83, yMax))
        lines[-1].SetLineStyle(root.kDotted)
        lines[-1].Draw('same')
        lines.append(root.TLine(-0.83, 0., -0.83, yMax))
        lines[-1].SetLineStyle(root.kDotted)
        lines[-1].Draw('same')
        lines.append(root.TLine(1.24, 0., 1.24, yMax))
        lines[-1].SetLineStyle(root.kDotted)
        lines[-1].Draw('same')
        lines.append(root.TLine(-1.24, 0., -1.24, yMax))
        lines[-1].SetLineStyle(root.kDotted)
        lines[-1].Draw('same')

    legend.Draw('same')

    tex = root.TLatex()
    tex.SetNDC()
    tex.SetTextFont(font)
    tex.SetTextSize(0.04)
    #tex.DrawLatex(0.555, 0.93, 'Simulation, 13 TeV')
    if canvWidth > 600:
        tex.DrawLatex(0.484, 0.93, 'CMS Simulation, 13 TeV')
    else:
        tex.DrawLatex(0.555, 0.93, 'CMS Simulation, 13 TeV')
    if notes:
        tex.SetTextSize(0.035)
        for note in notes:
            tex.DrawLatex(note[0], note[1], note[2])

    c.Modified()
    c.Update()

    return [c, hs, legend, lines, tex]

def hist_styles(stacked=False):
    styles = {}
    styles['gmt'] = {'lc':root.kCyan, 'ls':root.kSolid, 'fc':None, 'legtext':'GMT'}
    styles['ugmt'] = {'lc':root.kBlack, 'ls':root.kSolid, 'fc':None, 'legtext':'uGMT'}
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

def plot_hists_standard(hm, hName, den=None, hNamePrefix='', xTitle='', yTitle='# muons', stacked=False, normToBinWidth=False, tfMuonOrig='ugmt', reg=''):
    styles = hist_styles(stacked)

    if tfMuonOrig == 'ugmt':
        ugmt_str = '_ugmt'
        prefix = ''

    ugmt_dict = {'num':hNamePrefix+'ugmt_'+hName, 'den':den}
    ugmt_dict.update(styles['ugmt'])
    hDefs = []
    if reg == '':
        hDefs.append(ugmt_dict)

    # extract eta range from histogram name
    if hName.rfind('.eta') == -1:
        eta_number_strs = re.findall(r'[\d\.\d]+', hName[hName.find('EtaMin')+6:hName.find('EtaMax')+12])
        note_str = eta_number_strs[0]+' < |#eta^{reco}| < '+eta_number_strs[1]
        notes = [[0.17, 0.86, note_str]]
    else:
        notes = None

    return plot_hists(hm, hDefs, xTitle, yTitle, normToBinWidth, prefix, notes)

def plot_hists_qstack(hm, hName, den=None, hNamePrefix='', xTitle='', yTitle='# muons', normToBinWidth=False, tfMuonOrig='ugmt', reg=''):
    styles = hist_styles(stacked=False)

    if tfMuonOrig == 'ugmt':
        ugmt_str = '_ugmt'
        style_str = 'ugmt'
        prefix = 'q_'

    hDefs = []
    if reg == '':
        ugmt_dict = {'num':hNamePrefix+'ugmt_'+hName, 'den':den}
        ugmt_dict.update(styles['ugmt'])
        hDefs.append(ugmt_dict)
        for q in range(16):
            ugmt_q_dict = {'num':hNamePrefix+'ugmt_'+hName.replace('_matched', '_q{q}_matched'.format(q=q)), 'den':den}
            ugmt_q_dict.update(styles['ugmt_q{q}'.format(q=q)])
            hDefs.append(ugmt_q_dict)

    # extract eta range from histogram name
    if hName.rfind('.eta') == -1:
        eta_number_strs = re.findall(r'[\d\.\d]+', hName[hName.find('EtaMin')+6:hName.find('EtaMax')+12])
        note_str = eta_number_strs[0]+' < |#eta^{reco}| < '+eta_number_strs[1]
        notes = [[0.55, 0.86, note_str]]
    else:
        notes = None

    return plot_hists(hm, hDefs, xTitle, yTitle, normToBinWidth, prefix, notes)

def plot_hists_regional(hm, pt_mins, hNamePrefix='', xTitle='', yTitle='# muons', normToBinWidth=False, tfMuonOrig='ugmt', reg=''):
    styles = hist_styles(stacked=False)

    reco_cut_str = str(pt_mins[0])
    cut_str = str(pt_mins[1])

    if tfMuonOrig == 'ugmt':
        ugmt_str = '_ugmt'
        prefix = ''

    ugmt_dict = {'num':hNamePrefix+'ugmt_ptmin'+cut_str+'_matched_reco_muon_absEtaMin0_absEtaMax2.5_ptmin'+reco_cut_str+'.pt', 'den':'reco_muon_absEtaMin0_absEtaMax2.5_ptmin'+reco_cut_str+'.pt'}
    ugmt_dict.update(styles['ugmt'])
    hDefs = []
    if reg == '':
        hDefs.append(ugmt_dict)

    return plot_hists(hm, hDefs, xTitle, yTitle, normToBinWidth, prefix)

def plot_turn_on_curves(hm, eta_range = [0, 2.5], pt_cuts=[0.5], l1system='ugmt', ptSource='reco'):
    eta_min = eta_range[0]
    eta_max = eta_range[1]
    eta_min_str = '_absEtaMin'+str(eta_min)
    eta_max_str = '_absEtaMax'+str(eta_max)
    eta_title = '{eMin} < |#eta^{{reco}}| < {eMax}'.format(eMin=eta_min, eMax=eta_max)    

    colors = [root.kBlue, root.kRed, root.kGreen, root.kMagenta, root.kCyan, root.kOrange, root.kYellow, root.kPink]
    markers = [root.kFullCircle, root.kFullSquare, root.kFullTriangleUp, root.kFullTriangleDown, root.kOpenCircle, root.kOpenSquare, root.kOpenTriangleUp, root.kOpenTriangleDown]
    if ptSource == 'L1':
        canvas_name = 'turn_on_L1-pT_'+l1system+eta_min_str+eta_max_str
    else:
        canvas_name = 'turn_on_'+l1system+eta_min_str+eta_max_str
    canvas_title = canvas_name
    c = root.TCanvas(canvas_name, canvas_title, 100, 100, 600, 600)
    c.cd()
    
    labels = {'ugmt':'uGMT', 
             }

    set_root_style()

    legend = root.TLegend(0.68, 0.12, 0.95, 0.44, labels[l1system]+' p_{T} threshold')
    drawstring = ''
    hs = []
    for (marker, color, pt_min) in zip(markers, colors, pt_cuts):
        cut_str = '_ptmin'+str(pt_min)
        if ptSource == 'L1':
            hs.append(hm.get_ratio('best_reco'+eta_min_str+eta_max_str+'_ptmin0.5_matched_'+l1system+'_muon'+cut_str+'.pt', 'reco_muon'+eta_min_str+eta_max_str+'_ptmin0.5.pt').Clone(canvas_name+cut_str))
        else:
            hs.append(hm.get_ratio('best_'+l1system+cut_str+'_matched_reco_muon'+eta_min_str+eta_max_str+'_ptmin0.5.pt', 'reco_muon'+eta_min_str+eta_max_str+'_ptmin0.5.pt').Clone(canvas_name+cut_str))
        hs[-1].GetXaxis().SetTitleOffset(1.1)
        hs[-1].GetXaxis().SetTitle('p_{T}^{'+ptSource+'} (GeV/c)')
        hs[-1].GetYaxis().SetTitleOffset(1.5)
        hs[-1].GetYaxis().SetTitle(labels[l1system]+' efficiency')
        hs[-1].GetYaxis().SetRangeUser(0., 1.1)
        hs[-1].SetLineWidth(2)
        hs[-1].SetLineColor(color)
        #hs[-1].SetMarkerStyle(marker)
        hs[-1].SetMarkerColor(color)
        hs[-1].Draw(drawstring)
        drawstring = 'same'
        legend.AddEntry(hs[-1], str(pt_min)+' GeV/c', 'lep')
    
    legend.SetTextFont(font)
    legend.SetTextSize(0.03)
    legend.SetBorderSize(0)
    legend.SetFillColor(19)
    legend.SetFillStyle(0)
    #legend.SetNColumns(2)
    legend.Draw('same')

    tex = root.TLatex()
    tex.SetNDC()
    tex.SetTextFont(font)
    tex.SetTextSize(0.04)
    #tex.DrawLatex(0.555, 0.93, 'Simulation, 13 TeV')
    #tex.DrawLatex(0.68, 0.5, labels[l1system])
    tex.DrawLatex(0.555, 0.93, 'CMS Simulation, 13 TeV')
    tex.SetTextSize(0.035)
    tex.DrawLatex(0.17, 0.86, eta_title)

    c.Modified()
    c.Update()

    return [c, hs, legend, tex]

def main():
    opts = parse_options_plotEff(parser)
    batchRun = opts.interactive
    if batchRun:
        root.gROOT.SetBatch(True)

    # combinations of reco_pt_min and the corresponding pt_min values
    # combinations reco_pt_min == pt_min are necessary to fill non matched histograms
    # other combinations are optional
    ptmins_list = [[0.5, [0.5, 12, 16, 20, 24, 30]],
#    ptmins_list = [[0.5, [0.5, 16, 20, 30]],
#                  [12, [12]],
#                  [16, [12, 16]],
#                  [20, [12, 16, 20]],
#                  [24, [16, 20, 24]],
                  [24, [16]],
#                  [30, [20, 24, 30]],
                 ]

    eta_ranges = [[0, 2.5], [0, 0.83], [0.83, 1.24], [1.24, 2.5]]
    qualities = [12, 8]

    hm = HistManager(filename=opts.fname)

    # plot the histograms
    #cList = plot_histos(hm, cList)

    objects = []
    # plot the turn on curves vs threshold
    #objects.append(plot_turn_on_curves(hm, eta_ranges[0], ptmins_list[0][1], 'ugmt'))
    #objects.append(plot_turn_on_curves(hm, eta_ranges[1], ptmins_list[0][1], 'ugmt'))
    #objects.append(plot_turn_on_curves(hm, eta_ranges[2], ptmins_list[0][1], 'ugmt'))
    #objects.append(plot_turn_on_curves(hm, eta_ranges[3], ptmins_list[0][1], 'ugmt'))

    reco_0to2p5 = 'reco_muon_absEtaMin0_absEtaMax2.5_ptmin'
    reco_0to0p83 = 'reco_muon_absEtaMin0_absEtaMax0.83_ptmin'
    reco_0p83to1p24 = 'reco_muon_absEtaMin0.83_absEtaMax1.24_ptmin'
    reco_1p24to2p5 = 'reco_muon_absEtaMin1.24_absEtaMax2.5_ptmin'
    yTitle_best = '# best matched reco muons / # reco muons'

    # turn on plot for ugmt muons per region
    objects.append(plot_hists_standard(hm, 'ptmin16_matched_'+reco_0to2p5+'0.5.pt', reco_0to2p5+'0.5.pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_best, stacked=True))
    #objects.append(plot_hists_qstack(hm, 'ptmin16_matched_'+reco_0to2p5+'0.5.pt', reco_0to2p5+'0.5.pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_best))

    objects.append(plot_hists_standard(hm, 'ptmin16_matched_'+reco_0to0p83+'0.5.pt', reco_0to0p83+'0.5.pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_best, stacked=True))
    #objects.append(plot_hists_qstack(hm, 'ptmin16_matched_'+reco_0to0p83+'0.5.pt', reco_0to0p83+'0.5.pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_best))

    objects.append(plot_hists_standard(hm, 'ptmin16_matched_'+reco_0p83to1p24+'0.5.pt', reco_0p83to1p24+'0.5.pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_best, stacked=True))
    #objects.append(plot_hists_qstack(hm, 'ptmin16_matched_'+reco_0p83to1p24+'0.5.pt', reco_0p83to1p24+'0.5.pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_best))

    objects.append(plot_hists_standard(hm, 'ptmin16_matched_'+reco_1p24to2p5+'0.5.pt', reco_1p24to2p5+'0.5.pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_best, stacked=True))
    #objects.append(plot_hists_qstack(hm, 'ptmin16_matched_'+reco_1p24to2p5+'0.5.pt', reco_1p24to2p5+'0.5.pt', 'best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_best))

    # joint turn on plot for all subsystems
    #objects.append(plot_hists_regional(hm, [0.5, 16], hNamePrefix='best_', xTitle='p_{T}^{reco} (GeV/c)', yTitle=yTitle_best))

    ## eta plots
    objects.append(plot_hists_standard(hm, 'ptmin16_matched_'+reco_0to2p5+'24.eta', reco_0to2p5+'24.eta', 'best_', xTitle='#eta^{reco}', yTitle=yTitle_best, stacked=True))
    objects.append(plot_hists_standard(hm, 'ptmin16_matched_'+reco_0to2p5+'24.eta', hNamePrefix='best_', xTitle='#eta^{reco}', yTitle='# best matched reco muons', stacked=True))
    #objects.append(plot_hists_standard(hm, 'ptmin16_matched_'+reco_0to2p5+'24.eta', reco_0to2p5+'24.eta', 'best_', xTitle='#eta^{reco}', yTitle=yTitle_best, stacked=False))
    #objects.append(plot_hists_qstack(hm, 'ptmin16_matched_'+reco_0to2p5+'24.eta', reco_0to2p5+'24.eta', 'best_', xTitle='#eta^{reco}', yTitle=yTitle_best))

    # phi plots
    objects.append(plot_hists_standard(hm, 'ptmin16_matched_'+reco_0to2p5+'24.phi', reco_0to2p5+'24.phi', 'best_', xTitle='#phi^{reco}', yTitle=yTitle_best, stacked=True))
    #objects.append(plot_hists_standard(hm, 'ptmin16_matched_'+reco_0to2p5+'24.phi', reco_0to2p5+'24.phi', 'best_', xTitle='#phi^{reco}', yTitle=yTitle_best, stacked=False))
    #objects.append(plot_hists_qstack(hm, 'ptmin16_matched_'+reco_0to2p5+'24.phi', reco_0to2p5+'24.phi', 'best_', xTitle='#phi^{reco}', yTitle=yTitle_best))

    # delta R plots
    objects.append(plot_hists_standard(hm, 'ptmin16_matched_'+reco_0to2p5+'0.5.dr', hNamePrefix='best_', xTitle='#DeltaR', yTitle='# best matched reco muons', stacked=False))
    #objects.append(plot_hists_qstack(hm, 'ptmin16_matched_'+reco_0to2p5+'0.5.dr', hNamePrefix='best_', xTitle='#DeltaR', yTitle='# best matched reco muons'))

    # save canvases to root file
    if savePlots:
        output = root.TFile('./ugmt_eff_plots.root', 'recreate')
        output.cd()
        print objects
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
    savePlots = True
    batchRun = True
    best_only = False
    font = 42
    main()

