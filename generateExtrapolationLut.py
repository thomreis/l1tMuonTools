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
import math

def parse_options_plotRates(parser):
    """
    Adds often used options to the OptionParser...
    """
    parsers = parser.add_subparsers()
    sub_parser = parsers.add_parser("generateExtrapolationLut")
    sub_parser.add_argument("--coordinate", dest="coordinate", type=str, default='phi', help="Coordinate (eta or phi) to generate the LUT for.")
    sub_parser.add_argument("--pt-bits", dest="ptbits", type=int, default=6, help="Number of pT input bits.")
    sub_parser.add_argument("--eta-bits", dest="etabits", type=int, default=6, help="Number of eta input bits.")
    sub_parser.add_argument("--out-bits", dest="outbits", type=int, default=3, help="Number of output bits.")
    sub_parser.add_argument("--out-shift", dest="outshift", type=int, default=3, help="Number of left shifts of output bits.")

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

def fit_extrapolation_hists(hm, coordinate, eta_ranges, fit_range):
    functions = []
    for eta_range in eta_ranges:
        eta_min = eta_range[0]
        eta_max = eta_range[1]
        histo_name = 'l1_muon_absEtaMin{etaMin}_absEtaMax{etaMax}.pt_d{coord}'.format(etaMin=eta_min, etaMax=eta_max, coord=coordinate)
        h = hm.get(histo_name).Clone()

        function = root.TF1('func', '[0] * x^(-1*[1]) + [2]', fit_range[0], fit_range[1])
        function.SetParameters(1., 1., 0.)

        # first try
        print histo_name
        fit_res_ptr = h.Fit(function, 'RS')
        fit_res = fit_res_ptr.Get()

        # increase lower bound for fit if the previous fit failed
        if fit_res and not fit_res.IsValid():
            for red in range(1, 6):
                function.SetParameters(1., 1., 0.)
                #function.FixParameter(1, 1.)
                fit_res_ptr = h.Fit(function, 'RS', '', fit_range[0] + red, fit_range[1])
                fit_res = fit_res_ptr.Get()
                if fit_res and fit_res.IsValid():
                    break

        # fix parameter if no fit was successful so far
        if fit_res and not fit_res.IsValid():
            function.SetParameters(1., 1., 0.)
            function.FixParameter(1, 1.)
            fit_res_ptr = h.Fit(function, 'RS')
            fit_res = fit_res_ptr.Get()

        if fit_res:
            print fit_res.Chi2()/fit_res.Ndf()

        # in case the fit did not work set the function parameters to be the same as for the previous range
        if not fit_res or (fit_res and not fit_res.IsValid()):
            prev_function = functions[-1]
            function.SetParameters(prev_function.GetParameter(0), prev_function.GetParameter(1), prev_function.GetParameter(2))

        functions.append(function)

    return functions

def main():
    opts = parse_options_plotRates(parser)
    coord = opts.coordinate
    pt_bits = opts.ptbits
    red_eta_bits = opts.etabits
    lut_out_bits = opts.outbits
    lut_scale_factor = 2**opts.outshift

    root.gROOT.SetBatch(True)

    hm = HistManager(filename=opts.fname, subdir='all_runs')

    L1Ana.init_l1_analysis()
    print ""

    pt_scale = 0.5
    eta_scale = 0.010875
    phi_scale = 0.010908
    lut_pt_values = 2**pt_bits
    eta_bits = 8
    red_eta_scale = 2**(eta_bits - red_eta_bits) * eta_scale
    if coord == 'eta':
        lut_scale = eta_scale * lut_scale_factor
    else:
        lut_scale = phi_scale * lut_scale_factor

    # calculate eta ranges
    # The LUT uses a reduced eta coordinate with the two LSBs removed and the MSB masked.
    eta_ranges = []
    for red_hw_eta in range(2**red_eta_bits):
        eta_ranges.append((red_hw_eta*red_eta_scale, (red_hw_eta+1)*red_eta_scale))

    functions = fit_extrapolation_hists(hm, coord, eta_ranges, (2., lut_pt_values*pt_scale))

    lut_header = '# '+coord+' extrapolation LUT\n'
    lut_header += '# anything after # is ignored with the exception of the header\n'
    lut_header += '#<header> V1 {i} {o} </header>\n'.format(i=pt_bits+red_eta_bits, o=lut_out_bits)
    lut_payload = ''
    lut_str = ''
    lut_entry = 0
    for i, eta_range in enumerate(eta_ranges):
        lut_str += '{etaMin:1.4f}-{etaMax:1.4f}: '.format(etaMin=eta_range[0], etaMax=eta_range[1])
        lut_str += '0 '
        # take special care of 0 hwPt value as it gives infinity with the fit function
        lut_payload += '{i} 0\n'.format(i=lut_entry)
        lut_entry += 1
        for hwPt in range(1, lut_pt_values):
            func_val = functions[i].Eval(hwPt * pt_scale)
            lut_val = int(math.floor(func_val / lut_scale))
            if lut_val > 2**lut_out_bits - 1:
                lut_val = 2**lut_out_bits - 1
            elif lut_val < 0:
                lut_val = 0
            lut_payload += '{i} {val}\n'.format(i=lut_entry, val=lut_val)
            lut_str += '{val:2d} '.format(val=lut_val)
            lut_entry += 1
        lut_str += '\n'

    print lut_str

    out_file_name = 'lut.txt'
    with open(out_file_name, 'w') as out_file:
        out_file.write(lut_header + lut_payload)


if __name__ == "__main__":
    savePlots = True
    main()

