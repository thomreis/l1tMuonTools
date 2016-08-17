import ROOT as root
from array import array
import argparse


def parse_options_ntupling(parser):
    """
    Adds often used options to the OptionParser...
    """
    parsers = parser.add_subparsers()
    sub_parser = parsers.add_parser("ntuple")
    sub_parser.add_argument("-o", "--outname", dest="outname", default="ntuple", type=str, help="A root file name where to save the ntuple.")

    opts, unknown = parser.parse_known_args()
    return opts


class VarDescriptor(object):
    def __init__(self, name, dimname=None, dimmax=1, iscntr=False):
        self.name = name
        self.dimname = dimname
        self.dimmax = dimmax
        self.is_counter = iscntr


class Ntupler(object):
    """docstring for Ntupler"""
    def __init__(self, fname, ntuplename, vardescriptions):
        super(Ntupler, self).__init__()
        self.fname = fname
        self.ntuplename = ntuplename
        self.vardescs = vardescriptions

        self._fobj = root.TFile(fname, 'recreate')
        self._open = True
        self._fobj.cd()
        self._tree = root.TTree(ntuplename, '')

        self._contents = {}
        for var in self.vardescs:
            if not var.is_counter:
                self._contents[var.name] = array('f', [-9999]*var.dimmax)
                if var.dimmax != 1:
                    self._tree.Branch(var.name, self._contents[var.name], '{v}/F'.format(v=var.name))
                else:
                    self._tree.Branch(var.name, self._contents[var.name], '{v}/F'.format(v=var.name))
            else:
                self._contents[var.name] = array('i', [0])
                self._tree.Branch(var.name, self._contents[var.name], '{v}/I'.format(v=var.name))

    def reset(self):
        for var in self.vardescs:
            if var.is_counter:
                self._contents[var.name][0] = 0
            else:
                for i in range(var.dimmax):
                    self._contents[var.name][i] = -9999

    def set(self, vname, value):
        self._contents[vname][0] = value

    def save_to_tree(self):
        self._tree.Fill()

    def save_file(self):
        if self._open:
            self._fobj.Write()

    def close_file(self):
        if self._open:
            self._fobj.Close()
            self._open = False
