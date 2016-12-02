import ROOT as root
from array import array
from math import sqrt


class HistManager(object):
    """Class that manages and holds histograms"""
    def __init__(self, varnames=[], binning_dict={}, ytitle="# Muons", prefix="", filename=None, subdir=None):
        super(HistManager, self).__init__()
        self.varnames = varnames
        self.binnings = binning_dict
        self.prefix = prefix
        root.TGaxis().SetMaxDigits(3)

        self.hists = {}

        self._stackcache = {}
        self._effcache = {}
        self._ratiocache = {}
        self._thresholdcache = {}

        if filename is None:
            for vname in varnames:
                have_unit = type(self.binnings[vname][-2]) is str
                # variable binning when nBins == -1
                if self.binnings[vname][0] < 0:
                    if have_unit:
                        self.hists[vname] = root.TH1D(prefix+vname, "", len(self.binnings[vname])-4, array('d', self.binnings[vname][1:-2]))
                    else:
                        self.hists[vname] = root.TH1D(prefix+vname, "", len(self.binnings[vname])-3, array('d', self.binnings[vname][1:-1]))
                # fixed binning
                else:
                    self.hists[vname] = root.TH1D(prefix+vname, "", self.binnings[vname][0], self.binnings[vname][1], self.binnings[vname][2])
                self.hists[vname].Sumw2()
                if not have_unit:
                    xtitle = self.binnings[vname][-1]
                elif self.binnings[vname][-1] is None:
                    xtitle = self.binnings[vname][-2]
                else:
                    xtitle = "{title} ({unit})".format(title=self.binnings[vname][-2], unit=self.binnings[vname][-1])

                self.hists[vname].GetXaxis().SetTitle(xtitle)
                self.hists[vname].GetYaxis().SetTitle(ytitle)
        else:
            self.varnames = []
            input = root.TFile(filename)
            if subdir:
                directory = input.GetDirectory(subdir)
            else:
                directory = input.GetDirectory('')
            keyList = directory.GetListOfKeys()
            for key in keyList:
                hName = key.GetName()
                if not directory.Get(hName).InheritsFrom('TH1'):
                    continue
                self.varnames.append(hName)
                self.hists[hName] = directory.Get(hName)
                self.hists[hName].SetDirectory(0)
            input.Close()

    def fill(self, varname, val):
        self.hists[varname].Fill(val)

    def get(self, varname, addunderflow=False, addoverflow=False):
        h = self.hists[varname]
        if addunderflow:
            err = root.Double(0)
            integral = h.IntegralAndError(0, 1, err)
            h.SetBinContent(1, integral)
            h.SetBinError(1, err)
        if addoverflow:
            err = root.Double(0)
            integral = h.IntegralAndError(h.GetNbinsX(), h.GetNbinsX()+1, err)
            h.SetBinContent(h.GetNbinsX(), integral)
            h.SetBinError(h.GetNbinsX(), err)
        return h

    def get_varnames(self):
        return self.varnames

    def get_stack(self, varnames):
        keyname = "_".join(varnames)
        if keyname in self._stackcache.keys():
            return self._stackcache[keyname]
        stack = root.THStack()
        hs = []
        for vname in varnames:
            h = self.get(vname)
            stack.Add(h)
            hs.append(h)
        self._stackcache[keyname] = [hs, stack]
        return [hs, stack]

    def get_binning(self, varname):
        return self.binnings[varname]

    def get_threshold_hist(self, varname):
        if varname in self._thresholdcache.keys():
            return self._thresholdcache[varname]

        h_thr = self.hists[varname].Clone()
        h = self.hists[varname]
        bins = range(h_thr.GetNbinsX()+2)
        bins.reverse()
        bmax = h_thr.GetNbinsX()
        binsum = 0
        binerr2 = 0
        for ib in bins:
            binsum += h.GetBinContent(ib)
            binerr2 += h.GetBinError(ib)**2
            if ib > 0 and ib < bmax:
                h_thr.SetBinContent(ib, binsum)
                h_thr.SetBinError(ib, sqrt(binerr2))
        h_thr.GetYaxis().SetTitle("Integrated "+h_thr.GetYaxis().GetTitle())
        self._thresholdcache[varname] = h_thr
        return h_thr

    def get_threshold_stack(self, varnames):
        keyname = "thr_".join(varnames)
        if keyname in self._stackcache.keys():
            return self._stackcache[keyname]
        stack = root.THStack()
        hs = []
        for vname in varnames:
            h = self.get_threshold_hist(vname)
            stack.Add(h)
            hs.append(h)
        self._stackcache[keyname] = [hs, stack]
        return [hs, stack]

    def get_ratio(self, varname_nom, varname_denom, addunderflow=False, addoverflow=False):
        name = "{nom}_o_{denom}".format(nom=varname_nom, denom=varname_denom)
        if name in self._ratiocache.keys():
            return self._ratiocache[name]

        h_denom = self.hists[varname_denom]
        h_ratio = self.hists[varname_nom].Clone()
        if addunderflow:
            err = root.Double(0)
            integral = h_denom.IntegralAndError(0, 1, err)
            h_denom.SetBinContent(1, integral)
            h_denom.SetBinError(1, err)
            integral = h_ratio.IntegralAndError(0, 1, err)
            h_ratio.SetBinContent(1, integral)
            h_ratio.SetBinError(1, err)
        if addoverflow:
            err = root.Double(0)
            integral = h_denom.IntegralAndError(h_denom.GetNbinsX(), h_denom.GetNbinsX()+1, err)
            h_denom.SetBinContent(h_denom.GetNbinsX(), integral)
            h_denom.SetBinError(h_denom.GetNbinsX(), err)
            integral = h_ratio.IntegralAndError(h_ratio.GetNbinsX(), h_ratio.GetNbinsX()+1, err)
            h_ratio.SetBinContent(h_ratio.GetNbinsX(), integral)
            h_ratio.SetBinError(h_ratio.GetNbinsX(), err)
        h_ratio.Divide(h_ratio, h_denom, 1, 1, "b")

        if not addunderflow and not addoverflow:
            self._ratiocache[name] = h_ratio
        return h_ratio

    def get_ratio_stack(self, varnames_nom, varname_denom):
        name = "{nom}_o_{denom}".format(nom="_".join(varnames_nom), denom=varname_denom)
        if name in self._stackcache.keys():
            return self._stackcache[name]

        stack = root.THStack()
        hs = []
        for vname_nom in varnames_nom:
            h_ratio = self.get_ratio(vname_nom, varname_denom)
            hs.append(h_ratio)
            stack.Add(h_ratio)

        self._stackcache[name] = [hs, stack]
        return [hs, stack]

    def get_efficiency(self, varname_nom, varname_denom, addunderflow=False, addoverflow=False, rebin=1, removeProbeBinsNEntriesBelow=0):
        name = "{nom}_o_{denom}".format(nom=varname_nom, denom=varname_denom)
        if name in self._effcache.keys():
            return self._effcache[name]

        h_denom = self.hists[varname_denom].Clone()
        h_nom = self.hists[varname_nom].Clone()

        if rebin > 1:
            h_denom.Rebin(rebin)
            h_nom.Rebin(rebin)

        if addunderflow:
            err = root.Double(0)
            integral = h_denom.IntegralAndError(0, 1, err)
            h_denom.SetBinContent(1, integral)
            h_denom.SetBinError(1, err)
            integral = h_nom.IntegralAndError(0, 1, err)
            h_nom.SetBinContent(1, integral)
            h_nom.SetBinError(1, err)
        if addoverflow:
            err = root.Double(0)
            integral = h_denom.IntegralAndError(h_denom.GetNbinsX(), h_denom.GetNbinsX()+1, err)
            h_denom.SetBinContent(h_denom.GetNbinsX(), integral)
            h_denom.SetBinError(h_denom.GetNbinsX(), err)
            integral = h_nom.IntegralAndError(h_nom.GetNbinsX(), h_nom.GetNbinsX()+1, err)
            h_nom.SetBinContent(h_nom.GetNbinsX(), integral)
            h_nom.SetBinError(h_nom.GetNbinsX(), err)

        if removeProbeBinsNEntriesBelow > 0:
            for b in range(1, h_nom.GetNbinsX()):
                if h_denom.GetBinContent(b) < removeProbeBinsNEntriesBelow:
                    h_nom.SetBinContent(b, 0)
                    h_nom.SetBinError(b, 0)
                    h_denom.SetBinContent(b, 0)
                    h_denom.SetBinError(b, 0)

        eff = root.TEfficiency(h_nom, h_denom)
        if not addunderflow and not addoverflow:
            self._effcache[name] = eff
        return eff


    def get_efficiency_int(self, varname_nom, varname_denom, integrateToFromLeft=None, integrateToFromRight=None, rebin=1, removeProbeBinsNEntriesBelow=0):
        name = "{nom}_o_{denom}".format(nom=varname_nom, denom=varname_denom)
        if name in self._effcache.keys():
            return self._effcache[name]

        h_denom = self.hists[varname_denom].Clone()
        h_nom = self.hists[varname_nom].Clone()

        if rebin > 1:
            h_denom.Rebin(rebin)
            h_nom.Rebin(rebin)

        if integrateToFromLeft != None:
            err = root.Double(0)
            integral = h_denom.IntegralAndError(0, h_denom.FindBin(integrateToFromLeft), err)
            for b in range(h_denom.FindBin(integrateToFromLeft)):
                h_denom.SetBinContent(b+1, integral)
                h_denom.SetBinError(b+1, err)
            integral = h_nom.IntegralAndError(0, h_nom.FindBin(integrateToFromLeft), err)
            for b in range(h_nom.FindBin(integrateToFromLeft)):
                h_nom.SetBinContent(b+1, integral)
                h_nom.SetBinError(b+1, err)
        if integrateToFromRight != None:
            err = root.Double(0)
            # include the bin just below the axis cut if it falls on a bin edge
            rCutBin = h_denom.FindBin(integrateToFromRight)
            if h_denom.GetXaxis().GetBinLowEdge(rCutBin) == integrateToFromRight:
                rCutBin -= 1
            integral = h_denom.IntegralAndError(rCutBin, h_denom.GetNbinsX()+1, err)
            for b in range(rCutBin, h_denom.GetNbinsX()):
                h_denom.SetBinContent(b, integral)
                h_denom.SetBinError(b, err)
            rCutBin = h_nom.FindBin(integrateToFromRight)
            if h_nom.GetXaxis().GetBinLowEdge(rCutBin) == integrateToFromRight:
                rCutBin -= 1
            integral = h_nom.IntegralAndError(rCutBin, h_nom.GetNbinsX()+1, err)
            for b in range(rCutBin, h_nom.GetNbinsX()):
                h_nom.SetBinContent(b, integral)
                h_nom.SetBinError(b, err)

        if removeProbeBinsNEntriesBelow > 0:
            for b in range(1, h_nom.GetNbinsX()):
                if h_denom.GetBinContent(b) < removeProbeBinsNEntriesBelow:
                    h_nom.SetBinContent(b, 0)
                    h_nom.SetBinError(b, 0)
                    h_denom.SetBinContent(b, 0)
                    h_denom.SetBinError(b, 0)

        eff = root.TEfficiency(h_nom, h_denom)
        if not integrateToFromLeft and not integrateToFromRight:
            self._effcache[name] = eff
        return eff


class HistManager2d(object):
    """Class that manages and holds 2D histograms"""
    def __init__(self, varnames=[], binning_dict={}, profile_dict={}, ytitle="# Muons", prefix="", filename=None, subdir=None):
        super(HistManager2d, self).__init__()
        self.varnames = varnames
        self.binnings = binning_dict
        self.profiles = profile_dict
        self.prefix = prefix
        root.TGaxis().SetMaxDigits(3)

        self.hists = {}

        self._stackcache = {}
        self._effcache = {}
        self._ratiocache = {}
        self._thresholdcache = {}

        if filename is None:
            for vname in varnames:
                binning_x = self.binnings[vname][0]
                binning_y = self.binnings[vname][1]
                have_unit_x = type(binning_x[-2]) is str
                have_unit_y = type(binning_y[-2]) is str
                if have_unit_x:
                    nbinsx = len(binning_x)-4
                    binsx = array('d', binning_x[1:-2])
                else:
                    nbinsx = len(binning_x)-3
                    binsx = array('d', binning_x[1:-1])
                if have_unit_y:
                    nbinsy = len(binning_y)-4
                    binsy = array('d', binning_y[1:-2])
                else:
                    nbinsy = len(binning_y)-3
                    binsy = array('d', binning_y[1:-1])
                # variable binning when nBins == -1
                if binning_x[0] < 0 and binning_y[0] < 0:
                    if vname in self.profiles and self.profiles[vname] == True:
                        self.hists[vname] = root.TProfile2D(prefix+vname, "", nbinsx, binsx, nbinsy, binsy)
                    else:
                        self.hists[vname] = root.TH2D(prefix+vname, "", nbinsx, binsx, nbinsy, binsy)
                elif binning_x[0] < 0: # fixed binning on y axis
                    if vname in self.profiles and self.profiles[vname] == True:
                        self.hists[vname] = root.TProfile2D(prefix+vname, "", nbinsx, binsx, binning_y[0], binning_y[1], binning_y[2])
                    else:
                        self.hists[vname] = root.TH2D(prefix+vname, "", nbinsx, binsx, binning_y[0], binning_y[1], binning_y[2])
                elif binning_y[0] < 0: # fixed binning on x axis
                    if vname in self.profiles and self.profiles[vname] == True:
                        self.hists[vname] = root.TProfile2D(prefix+vname, "", binning_x[0], binning_x[1], binning_x[2], nbinsy, binsy)
                    else:
                        self.hists[vname] = root.TH2D(prefix+vname, "", binning_x[0], binning_x[1], binning_x[2], nbinsy, binsy)
                else: # fixed binning
                    if vname in self.profiles and self.profiles[vname] == True:
                        self.hists[vname] = root.TProfile2D(prefix+vname, "", binning_x[0], binning_x[1], binning_x[2], binning_y[0], binning_y[1], binning_y[2])
                    else:
                        self.hists[vname] = root.TH2D(prefix+vname, "", binning_x[0], binning_x[1], binning_x[2], binning_y[0], binning_y[1], binning_y[2])
                self.hists[vname].Sumw2()

                if not have_unit_x:
                    xtitle = binning_x[-1]
                elif binning_x[-1] is None:
                    xtitle = binning_x[-2]
                else:
                    xtitle = "{title} ({unit})".format(title=binning_x[-2], unit=binning_x[-1])

                if not have_unit_y:
                    ytitle = binning_y[-1]
                elif binning_y[-1] is None:
                    ytitle = binning_y[-2]
                else:
                    ytitle = "{title} ({unit})".format(title=binning_y[-2], unit=binning_y[-1])

                self.hists[vname].GetXaxis().SetTitle(xtitle)
                self.hists[vname].GetYaxis().SetTitle(ytitle)
        else:
            self.varnames = []
            input = root.TFile(filename)
            directory = input.GetDirectory(subdir)
            keyList = directory.GetListOfKeys()
            for key in keyList:
                hName = key.GetName()
                if not directory.Get(hName).InheritsFrom('TH2'):
                    continue
                self.varnames.append(hName)
                self.hists[hName] = directory.Get(hName)
                self.hists[hName].SetDirectory(0)
            input.Close()

    def fill(self, varname, valx, valy, valz=1.):
        if varname in self.profiles and self.profiles[varname] == True:
            self.hists[varname].Fill(valx, valy, valz)
        else:
            self.hists[varname].Fill(valx, valy)

    def get(self, varname):
        return self.hists[varname]

    def get_varnames(self):
        return self.varnames

    def get_binning(self, varname):
        return self.binnings[varname]

    def get_profile(self, varname):
        if varname in self.profiles:
            return self.profiles[varname]
        else:
            return False


class L1AnalysisHistManager(HistManager):
    """Class that manages and holds histograms"""
    def __init__(self, varnames, binning_dict, prefix=""):
        super(L1AnalysisHistManager, self).__init__(varnames, binning_dict, prefix=prefix)

    def fill(self, varname, val):
        self.hists[varname].Fill(val)


class VarExp(object):
    def __init__(self, name, varexp, cutexp):
        super(VarExp, self).__init__()
        self.name = name
        self.varexp = varexp
        self.cutexp = cutexp


class FlatHistManager(HistManager):
    def __init__(self, varnames, binning_dict, varexps, cutexp, prefix=""):
        super(FlatHistManager, self).__init__(varnames, binning_dict, prefix=prefix)
        self.varexps = varexps
        self.cutexp = cutexp

    def project(self, ntuple):
        for vname in self.varnames:
            ntuple.Project(self.prefix+vname, self.varexps[vname].varexp, "({cut1})*({cut2})".format(cut1=self.varexps[vname].cutexp, cut2=self.cutexp))
