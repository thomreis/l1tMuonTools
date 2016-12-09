import ROOT as root

from sys import exit
import logging


class L1Ana(object):

    """
    Static class that manages all the pyROOT / pyFWLite and initialization part
    It also manages logging for those who like colourful prompts
    """
    error_pre = '\x1b[31;01m'
    info_pre = '\x1b[32;01m'
    warning_pre = '\x1b[33;01m'
    fatal_pre = '\x1b[31;31m'
    debug_pre = '\x1b[36;01m'
    reset = '\x1b[39;49;00m'
    log = None

    loginit = False
    l1init = False

    @staticmethod
    def init_l1_analysis():
        """
        Initialize all the tools that are already present for analysing L1Ntuples
        Through ROOT we can just use whatever is already implemented
        """
        if not L1Ana.loginit:
            L1Ana.init_logging()

        # Import FWLite
        L1Ana.log.info("FWCoreFWLite library being loaded.")
        root.gSystem.Load("libFWCoreFWLite")
        root.gROOT.ProcessLine('AutoLibraryLoader::enable();')
        # root.gSystem.Load("libCintex")
        # root.gROOT.ProcessLine('ROOT::Cintex::Cintex::Enable();')

        # append the paths to the macro classes
        curr_path = root.gEnv.GetValue("Unix.*.Root.MacroPath", "")
        work_path = root.gROOT.GetMacroPath()

        # append the include directories to get access to
        # L1Analysis-DataFormats:
        L1Ana.log.info("Adding L1Ntuple include directories")
        root.gSystem.AddIncludePath(
            " -I$CMSSW_BASE/src/L1Trigger/L1NTtuples/interface")
        root.gROOT.ProcessLine(
            ".include $CMSSW_BASE/src/L1Trigger/L1TNtuples/interface")

        L1Ana.log.info("--- Initialization done ---")
        L1Ana.l1init = True

    @staticmethod
    def init_logging(name="L1Ana", level=None):
        """
        Initialize a logger with different colors for importance-levels
        """
        if not L1Ana.loginit or level != None:
            L1Ana.log = logging.getLogger(name)
            if level != None:
                L1Ana.log.setLevel(level)
            else:
                L1Ana.log.setLevel(logging.DEBUG)
            logging.addLevelName(
                logging.FATAL,   L1Ana.fatal_pre + logging.getLevelName(logging.FATAL) + L1Ana.reset)
            logging.addLevelName(
                logging.ERROR,   L1Ana.error_pre + logging.getLevelName(logging.ERROR) + L1Ana.reset)
            logging.addLevelName(
                logging.WARNING, L1Ana.warning_pre + logging.getLevelName(logging.WARNING) + L1Ana.reset)
            logging.addLevelName(
                logging.INFO,    L1Ana.info_pre + logging.getLevelName(logging.INFO) + L1Ana.reset)
            logging.addLevelName(
                logging.DEBUG,   L1Ana.debug_pre + logging.getLevelName(logging.DEBUG) + L1Ana.reset)

            logging.basicConfig(
                level=level, format='%(asctime)s (%(name)s) [%(levelname)s]: %(message)s', datefmt='%H:%M:%S')
            L1Ana.loginit = True


class L1Data(object):

    """
    This is the container that is returned by the iterator:
    The user will basically work on this container only and have access
    to all L1Ntuple-DataFormats through this.
    """

    def __init__(self):
        super(L1Data, self).__init__()
        self.event = None
        self.simulation = None
        self.gct = None
        self.gmt = None
        self.upgrade = None
        self.upgradeEmu = None
        self.upgradeBmtf = None
        self.upgradeBmtfEmu = None
        self.upgradeOmtf = None
        self.upgradeOmtfEmu = None
        self.upgradeEmtf = None
        self.upgradeEmtfEmu = None
        self.ugmt = None
        self.towers2x2 = None
        self.towers = None
        self.caloTowers = None
        self.caloTowersEmu = None
        self.gt = None
        self.rct = None
        self.dttf = None
        self.csctf = None
        self.recoMet = None
        self.recoMuon = None
        self.recoRpcHit = None
        self.recoJet = None
        self.recoBasicCluster = None
        self.recoSuperCluster = None
        self.l1extra = None
        self.l1emuextra = None
        self.recoVertex = None
        self.recoTrack = None
        self.l1menu = None
        self.gen = None


class L1Ntuple(object):

    """
    The interface to the user, it is based on the L1NTuple c++ class
    """

    def __init__(self, nevents=-1):
        super(L1Ntuple, self).__init__()
        self.data = L1Data()
        self.do_upgrade = False
        self.do_upgradeEmu = False
        self.do_upgradeTf = False
        self.do_upgradeTfEmu = False
        self.do_caloTowers = False
        self.do_caloTowersEmu = False
        self.do_reco = False
        self.do_muonreco = False
        self.do_l1extra = False
        self.do_l1emuextra = False
        self.do_l1menu = False
        self.do_muonup = False

        self.tree_main = None
        self.tree_upgrade = None
        self.tree_upgradeEmu = None
        self.tree_upgradeTf = None
        self.tree_upgradeTfEmu = None
        self.tree_caloTowers = None
        self.tree_caloTowersEmu = None
        self.tree_muon = None
        self.tree_muonupgrade = None
        self.tree_reco = None
        self.tree_extra = None
        self.tree_menu = None
        self.tree_emu_extra = None

        self.file_list = []
        self.nentries = -1
        self.nevents = nevents
        self.current = 0
        self.curr_file = None
        self.init = False

    def open_with_file_list(self, fname_list):
        """
        Initilize with a text file containing all root-files with
        L1Ntuples (one per line)
        TAKES: The file name pointing to the txt-file
        """
        if not L1Ana.l1init:
            L1Ana.init_l1_analysis()
        self.open_file_list(fname_list)
        self.check_first_file()
        self.open_no_init()
        self.init_branches()

        L1Ana.log.info("Ready to analyse.")
        self.init = True

    def open_with_file(self, fname):
        """
        Initialize with only one root-file containing a L1Ntuple
        TAKES: The file name pointing to the root-file
        """
        if not L1Ana.l1init:
            L1Ana.init_l1_analysis()
        self.file_list = [fname]
        self.check_first_file()
        self.open_no_init()
        self.init_branches()
        L1Ana.log.info("Ready to analyse.")
        self.init = True

    def open_no_init(self):
        """
        Initializes the TChains and adds present Trees as friends to the main tree
        this is needed so the tree's GetEntry stays in synch.
        """
        self.tree_main = root.TChain("l1EventTree/L1EventTree")
        self.tree_upgrade = root.TChain("l1UpgradeTree/L1UpgradeTree")
        self.tree_upgradeEmu = root.TChain("l1UpgradeEmuTree/L1UpgradeTree")
        self.tree_upgradeTf = root.TChain("l1UpgradeTfMuonTree/L1UpgradeTfMuonTree")
        self.tree_upgradeTfEmu = root.TChain("l1UpgradeTfMuonEmuTree/L1UpgradeTfMuonTree")
        self.tree_caloTowers = root.TChain("l1CaloTowerTree/L1CaloTowerTree")
        self.tree_caloTowersEmu = root.TChain("l1CaloTowerEmuTree/L1CaloTowerTree")
        self.tree_muon = root.TChain("l1MuonRecoTree/Muon2RecoTree")
        self.tree_reco = root.TChain("l1RecoTree/RecoTree")
        self.tree_extra = root.TChain("l1ExtraTreeProducer/L1ExtraTree")
        self.tree_emu_extra = root.TChain("l1EmulatorExtraTree/L1ExtraTree")
        self.tree_menu = root.TChain("l1MenuTreeProducer/L1MenuTree")
        self.tree_muonupgrade = root.TChain("l1MuonUpgradeTreeProducer/L1MuonUpgradeTree")

        for fname in self.file_list:
            self.tree_main.Add(fname)
            if self.do_upgrade:
                self.tree_upgrade.Add(fname)
            if self.do_upgradeEmu:
                self.tree_upgradeEmu.Add(fname)
            if self.do_upgradeTf:
                self.tree_upgradeTf.Add(fname)
            if self.do_upgradeTfEmu:
                self.tree_upgradeTfEmu.Add(fname)
            if self.do_caloTowers:
                self.tree_caloTowers.Add(fname)
            if self.do_caloTowersEmu:
                self.tree_caloTowersEmu.Add(fname)
            if self.do_reco:
                self.tree_reco.Add(fname)
            if self.do_muonreco:
                self.tree_muon.Add(fname)
            if self.do_l1emuextra:
                self.tree_emu_extra.Add(fname)
            if self.do_l1extra:
                self.tree_extra.Add(fname)
            if self.do_l1menu:
                self.tree_menu.Add(fname)
            if self.do_muonup:
                self.tree_muonupgrade.Add(fname)

        if self.do_upgrade:
            self.tree_main.AddFriend(self.tree_upgrade)
        if self.do_upgradeEmu:
            self.tree_main.AddFriend(self.tree_upgradeEmu)
        if self.do_upgradeTf:
            self.tree_main.AddFriend(self.tree_upgradeTf)
        if self.do_upgradeTfEmu:
            self.tree_main.AddFriend(self.tree_upgradeTfEmu)
        if self.do_caloTowers:
            self.tree_main.AddFriend(self.tree_caloTowers)
        if self.do_caloTowersEmu:
            self.tree_main.AddFriend(self.tree_caloTowersEmu)
        if self.do_reco:
            self.tree_main.AddFriend(self.tree_reco)
        if self.do_muonreco:
            self.tree_main.AddFriend(self.tree_muon)
        if self.do_l1emuextra:
            self.tree_main.AddFriend(self.tree_emu_extra)
        if self.do_l1extra:
            self.tree_main.AddFriend(self.tree_extra)
        if self.do_l1menu:
            self.tree_main.AddFriend(self.tree_menu)
        if self.do_muonup:
            self.tree_main.AddFriend(self.tree_muonupgrade)
        L1Ana.log.info("Files added to TChains.")

    def open_file_list(self, fname_list):
        """
        Open the txt-file and add all file-names to the list of files
        TAKES: file-name pointing to txt-file with one L1Ntuple-file per line
        """
        self.file_list = []
        L1Ana.log.info("Reading txt file with root-file list.")
        cntr = 0
        try:
            with open(fname_list) as fobj:
                for line in fobj:
                    fname = line.strip()
                    if fname == "":
                        continue
                    self.file_list.append(fname)
                    cntr += 1
        except EnvironmentError:
            L1Ana.log.fatal(
                "While reading file (probably it does not exist): {fname}".format(fname=fname_list))
            exit(0)

        L1Ana.log.info("Found list of {n} files:".format(n=cntr))
        for name in self.file_list:
            L1Ana.log.info("-- {fname}".format(fname=name))

    def check_first_file(self):
        """
        Checks which branches and trees are present in the first root-file.
        Sets flags which ones to add accordingly.
        """
        if not self.file_list:
            L1Ana.log.fatal("No root-files specified")
            exit(0)

        self.curr_file = root.TFile.Open(self.file_list[0])

        if not self.curr_file:
            L1Ana.log.fatal(
                "Could not open file: {fname}".format(fname=self.file_list[0]))
            exit(0)
        if self.curr_file.IsOpen() == 0:
            L1Ana.log.fatal(
                "Could not open file: {fname}".format(fname=self.file_list[0]))
            exit(0)

        my_chain = self.curr_file.Get("l1EventTree/L1EventTree")
        upgrade = self.curr_file.Get("l1UpgradeTree/L1UpgradeTree")
        upgradeEmu = self.curr_file.Get("l1UpgradeEmuTree/L1UpgradeTree")
        upgradeTf = self.curr_file.Get("l1UpgradeTfMuonTree/L1UpgradeTfMuonTree")
        upgradeTfEmu = self.curr_file.Get("l1UpgradeTfMuonEmuTree/L1UpgradeTfMuonTree")
        caloTowers = self.curr_file.Get("l1CaloTowerTree/L1CaloTowerTree")
        caloTowersEmu = self.curr_file.Get("l1CaloTowerEmuTree/L1CaloTowerTree")
        muon = self.curr_file.Get("l1MuonRecoTree/Muon2RecoTree")
        jets = self.curr_file.Get("l1RecoTree/RecoTree")
        extra = self.curr_file.Get("l1ExtraTreeProducer/L1ExtraTree")
        emuextra = self.curr_file.Get("l1EmulatorExtraTree/L1ExtraTree")
        menu = self.curr_file.Get("l1MenuTreeProducer/L1MenuTree")
        muonup = self.curr_file.Get("l1MuonUpgradeTreeProducer/L1MuonUpgradeTree")

        if my_chain:
            L1Ana.log.info("Found L1EventTree...")
        else:
            L1Ana.log.fatal("Could not find the main L1EventTree.")
            exit(0)

        if upgrade:
            L1Ana.log.info("Found L1UpgradeTree... Will add access to it.")
            self.do_upgrade = True
        else:
            L1Ana.log.warning(
                "Could not find L1UpgradeTree... It will be skipped.")

        if upgradeEmu:
            L1Ana.log.info("Found L1UpgradeEmuTree... Will add access to it.")
            self.do_upgradeEmu = True
        else:
            L1Ana.log.warning(
                "Could not find L1UpgradeEmuTree... It will be skipped.")

        if upgradeTf:
            L1Ana.log.info("Found L1UpgradeTfMuonTree... Will add access to it.")
            self.do_upgradeTf = True
        else:
            L1Ana.log.warning(
                "Could not find L1UpgradeTfMuonTree... It will be skipped.")

        if upgradeTfEmu:
            L1Ana.log.info("Found L1UpgradeTfMuonEmuTree... Will add access to it.")
            self.do_upgradeTfEmu = True
        else:
            L1Ana.log.warning(
                "Could not find L1UpgradeTfMuonEmuTree... It will be skipped.")

        if caloTowers:
            L1Ana.log.info("Found L1CaloTowerTree... Will add access to it.")
            self.do_caloTowers = True
        else:
            L1Ana.log.warning(
                "Could not find L1CaloTowerTree... It will be skipped.")

        if caloTowersEmu:
            L1Ana.log.info("Found L1CaloTowerEmuTree... Will add access to it.")
            self.do_caloTowersEmu = True
        else:
            L1Ana.log.warning(
                "Could not find L1CaloTowerEmuTree... It will be skipped.")

        if muon:
            L1Ana.log.info("Found MuonRecoTree... Will add access to it.")
            self.do_muonreco = True
        else:
            L1Ana.log.warning(
                "Could not find MuonRecoTree... It will be skipped.")

        if jets:
            L1Ana.log.info("Found RecoTree... Will add access to it.")
            self.do_reco = True
        else:
            L1Ana.log.warning("Could not find RecoTree... It will be skipped.")

        if extra:
            L1Ana.log.info("Found L1Extra... Will add access to it.")
            self.do_l1extra = True
        else:
            L1Ana.log.warning(
                "Could not find L1ExtraTree... It will be skipped.")

        if emuextra:
            L1Ana.log.info("Found L1EmulatorExtra... Will add access to it.")
            self.do_l1emuextra = True
        else:
            L1Ana.log.warning(
                "Could not find L1EmulatorExtra... It will be skipped.")

        if menu:
            L1Ana.log.info("Found L1MenuTree... Will add access to it.")
            self.do_l1menu = True
        else:
            L1Ana.log.warning(
                "Could not find L1MenuTree... It will be skipped.")

        if muonup:
            L1Ana.log.info("Found MuonUpgradeTree... Will add access to it.")
            self.do_muonup = True
        else:
            L1Ana.log.warning(
                "Could not find MuonUpgradeTree... It will be skipped.")

    def init_branches(self):
        """
        Connect the branches of the Trees with the corresponding members
        in the L1Data container.
        """
        if not self.tree_main:
            L1Ana.log.fatal(
                "There is no main L1Tree -- aborting initialization of branches")
            exit(0)

        self.nentries = self.tree_main.GetEntries()
        if self.nevents < 0 or self.nevents > self.nentries:
            self.nevents = self.nentries

        L1Ana.log.info("Approximate number of entries: {n}, running over: {n2}".format(
            n=self.nentries, n2=self.nevents))

        self.data.event = root.L1Analysis.L1AnalysisEventDataFormat()

        L1Ana.log.info("Setting branch addresses for main L1Tree.")

        self.tree_main.SetBranchAddress("Event", root.AddressOf(self.data.event))

        if self.tree_main.GetBranch("GCT"):
            self.data.gct = root.L1Analysis.L1AnalysisGCTDataFormat()
            self.tree_main.SetBranchAddress("GCT", root.AddressOf(self.data.gct))
        else:
            L1Ana.log.warning("GCT branch not present...")

        if self.tree_main.GetBranch("GMT"):
            self.data.gmt = root.L1Analysis.L1AnalysisGMTDataFormat()
            self.tree_main.SetBranchAddress("GMT", root.AddressOf(self.data.gmt))
        else:
            L1Ana.log.warning("GMT branch not present...")

        if self.tree_main.GetBranch("GT"):
            self.data.gt = root.L1Analysis.L1AnalysisGTDataFormat()
            self.tree_main.SetBranchAddress("GT", root.AddressOf(self.data.gt))
        else:
            L1Ana.log.warning("GT branch not present...")

        if self.tree_main.GetBranch("RCT"):
            self.data.rct = root.L1Analysis.L1AnalysisRCTDataFormat()
            self.tree_main.SetBranchAddress("RCT", root.AddressOf(self.data.rct))
        else:
            L1Ana.log.warning("RCT branch not present...")

        if self.tree_main.GetBranch("CSCTF"):
            self.data.csctf = root.L1Analysis.L1AnalysisCSCTFDataFormat()
            self.tree_main.SetBranchAddress(
                "CSCTF", root.AddressOf(self.data.csctf))
        else:
            L1Ana.log.warning("CSCTF branch not present...")

        if self.tree_main.GetBranch("DTTF"):
            self.data.dttf = root.L1Analysis.L1AnalysisDTTFDataFormat()
            self.tree_main.SetBranchAddress("DTTF", root.AddressOf(self.data.dttf))
        else:
            L1Ana.log.warning("DTTF branch not present...")

        if self.tree_main.GetBranch("Simulation"):
            self.data.simulation = root.L1Analysis.L1AnalysisSimulationDataFormat()
            self.tree_main.SetBranchAddress(
                "Simulation", root.AddressOf(self.data.simulation))
        else:
            L1Ana.log.warning("Simulation branch not present...")

        if self.tree_main.GetBranch("Generator"):
            self.data.gen = root.L1Analysis.L1AnalysisGeneratorDataFormat()
            self.tree_main.SetBranchAddress(
                "Generator", root.AddressOf(self.data.gen))
        else:
            L1Ana.log.warning("Generator branch not present...")

        if self.do_upgrade:
            L1Ana.log.info("Setting branch addresses for L1UpgradeTree")
            self.data.upgrade = root.L1Analysis.L1AnalysisL1UpgradeDataFormat()
            self.tree_upgrade.SetBranchAddress("L1Upgrade", root.AddressOf(self.data.upgrade))

        if self.do_upgradeEmu:
            L1Ana.log.info("Setting branch addresses for L1UpgradeEmuTree")
            self.data.upgradeEmu = root.L1Analysis.L1AnalysisL1UpgradeDataFormat()
            self.tree_upgradeEmu.SetBranchAddress("L1Upgrade", root.AddressOf(self.data.upgradeEmu))

        if self.do_upgradeTf:
            L1Ana.log.info("Setting branch addresses for L1UpgradeTfMuonTree")
            self.data.upgradeBmtf = root.L1Analysis.L1AnalysisL1UpgradeTfMuonDataFormat()
            self.data.upgradeOmtf = root.L1Analysis.L1AnalysisL1UpgradeTfMuonDataFormat()
            self.data.upgradeEmtf = root.L1Analysis.L1AnalysisL1UpgradeTfMuonDataFormat()
            self.tree_upgradeTf.SetBranchAddress("L1UpgradeBmtfMuon", root.AddressOf(self.data.upgradeBmtf))
            self.tree_upgradeTf.SetBranchAddress("L1UpgradeOmtfMuon", root.AddressOf(self.data.upgradeOmtf))
            self.tree_upgradeTf.SetBranchAddress("L1UpgradeEmtfMuon", root.AddressOf(self.data.upgradeEmtf))

        if self.do_upgradeTfEmu:
            L1Ana.log.info("Setting branch addresses for L1UpgradeTfMuonEmuTree")
            self.data.upgradeBmtfEmu = root.L1Analysis.L1AnalysisL1UpgradeTfMuonDataFormat()
            self.data.upgradeOmtfEmu = root.L1Analysis.L1AnalysisL1UpgradeTfMuonDataFormat()
            self.data.upgradeEmtfEmu = root.L1Analysis.L1AnalysisL1UpgradeTfMuonDataFormat()
            self.tree_upgradeTfEmu.SetBranchAddress("L1UpgradeBmtfMuon", root.AddressOf(self.data.upgradeBmtfEmu))
            self.tree_upgradeTfEmu.SetBranchAddress("L1UpgradeOmtfMuon", root.AddressOf(self.data.upgradeOmtfEmu))
            self.tree_upgradeTfEmu.SetBranchAddress("L1UpgradeEmtfMuon", root.AddressOf(self.data.upgradeEmtfEmu))

        if self.do_caloTowers:
            L1Ana.log.info("Setting branch addresses for L1CaloTowerTree")
            self.data.caloTowers = root.L1Analysis.L1AnalysisL1CaloTowerDataFormat()
            self.tree_caloTowers.SetBranchAddress("L1CaloTower", root.AddressOf(self.data.caloTowers))

        if self.do_caloTowersEmu:
            L1Ana.log.info("Setting branch addresses for L1CaloTowerEmuTree")
            self.data.caloTowersEmu = root.L1Analysis.L1AnalysisL1CaloTowerDataFormat()
            self.tree_caloTowersEmu.SetBranchAddress("L1CaloTower", root.AddressOf(self.data.caloTowersEmu))

        if self.do_reco:
            L1Ana.log.info("Setting branch addresses for RecoTree.")
            #self.data.recoJet = root.L1Analysis.L1AnalysisRecoJetDataFormat()
            #self.data.recoBasicCluster = root.L1Analysis.L1AnalysisRecoClusterDataFormat()
            #self.data.recoSuperCluster = root.L1Analysis.L1AnalysisRecoClusterDataFormat()
            #self.data.recoMet = root.L1Analysis.L1AnalysisRecoMetDataFormat()
            #self.data.recoTrack = root.L1Analysis.L1AnalysisRecoTrackDataFormat()
            self.data.recoVertex = root.L1Analysis.L1AnalysisRecoVertexDataFormat()

            #self.tree_reco.SetBranchAddress(
            #    "Jet", root.AddressOf(self.data.recoJet))
            #self.tree_reco.SetBranchAddress(
            #    "BasicClusters", root.AddressOf(self.data.recoBasicCluster))
            #self.tree_reco.SetBranchAddress(
            #    "SuperClusters", root.AddressOf(self.data.recoSuperCluster))
            #self.tree_reco.SetBranchAddress(
            #    "Met", root.AddressOf(self.data.recoMet))
            #self.tree_reco.SetBranchAddress(
            #    "Tracks", root.AddressOf(self.data.recoTrack))
            self.tree_reco.SetBranchAddress(
                "Vertex", root.AddressOf(self.data.recoVertex))

        if self.do_muonreco:
            L1Ana.log.info("Setting branch addresses for MuonRecoTree.")
            self.data.recoMuon = root.L1Analysis.L1AnalysisRecoMuon2DataFormat()
            self.data.recoRpcHit = root.L1Analysis.L1AnalysisRecoRpcHitDataFormat()

            self.tree_muon.SetBranchAddress(
                "Muon", root.AddressOf(self.data.recoMuon))
            if self.tree_muon.GetBranch("RpcHit"):
                self.tree_muon.SetBranchAddress(
                    "RpcHit", root.AddressOf(self.data.recoRpcHit))
            else:
                L1Ana.log.warning("RpcHit branch not present...")

        if self.do_l1menu:
            L1Ana.log.info("Setting branch addresses for L1Menu.")
            self.data.l1menu = root.L1Analysis.L1AnalysisL1MenuDataFormat()
            self.tree_menu.SetBranchAddress(
                "L1Menu", root.AddressOf(self.data.l1menu))
        if self.do_l1extra:
            L1Ana.log.info("Setting branch addresses for L1Extra.")
            self.data.l1extra = root.L1Analysis.L1AnalysisL1ExtraDataFormat()
            self.tree_extra.SetBranchAddress(
                "L1Extra", root.AddressOf(self.data.l1emuextra))
        if self.do_l1emuextra:
            L1Ana.log.info("Setting branch addresses for L1EmuExtra.")
            self.data.l1emuextra = root.L1Analysis.L1AnalysisL1ExtraDataFormat()
            self.tree_emu_extra.SetBranchAddress(
                "L1Extra", root.AddressOf(self.data.l1emuextra))
        if self.do_muonup:
            L1Ana.log.info("Setting branch addresses for MuonUpgradeTree")
            self.data.ugmt = root.L1Analysis.L1AnalysisUGMTDataFormat(8)
            self.data.towers2x2 = root.L1Analysis.L1AnalysisMuTwrDataFormat()
            self.data.towers = root.L1Analysis.L1AnalysisMuTwrDataFormat()
            self.tree_muonupgrade.SetBranchAddress(
                "L1TMuon", root.AddressOf(self.data.ugmt))
            self.tree_muonupgrade.SetBranchAddress(
                "L1TMuonCalo2x2", root.AddressOf(self.data.towers2x2))
            self.tree_muonupgrade.SetBranchAddress(
                "L1TMuonCalo", root.AddressOf(self.data.towers))

    def __len__(self):
        """
        RETURNS: number of entries
        """
        if not self.init:
            L1Ana.log.error("No estimate for elements, yet!")
            return -1
        return self.nentries

    def __getitem__(self, index):
        """
        This is the iterator, it will get the next entry and return the updated L1Data container
        TAKES: index of next event
        RETURNS: L1Data object with connected with L1Ntuple content
        """
        if not self.init:
            L1Ana.log.error(
                "L1Ntuple is not yet initialized! Aborting iteration.")
            raise IndexError("L1Ntuple is not yet initialized!")
        if not index < self.nentries:
            raise IndexError("Reached the end")

        self.tree_main.GetEntry(index)
        return self.data
