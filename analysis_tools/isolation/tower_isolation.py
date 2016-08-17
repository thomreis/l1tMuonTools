import math
import ROOT as root

class TowerIsolator(object):
    eta_tower = [0, 0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696, 0.783, 0.870, 0.957, 1.044, 1.131, 1.218,
                 1.305, 1.392, 1.479, 1.566, 1.653, 1.740, 1.830, 1.930, 2.043, 2.172, 2.322, 2.5, 2.650, 3.000, 3.5, 4.0, 4.5, 5.0]
    pi2 = 2*math.pi
    functions = {}
    ieta_to_ext_bin = {}
    calibration_factor = []

    @staticmethod
    def init_extrapolation(fname):
        fobj = root.TFile.Open(fname)
        for i, eta_bin in enumerate([0.5, 0.8, 1.1, 1.3, 1.5, 1.7, 1.8, 1.9, 2.0, 2.1]):
            TowerIsolator.functions[i] = fobj.Get("Reco_absdPhiME2Vertex_recoMaxEta"+str(eta_bin))

    @staticmethod
    def init_calibration():
        for i in range(100):
            if i in [0, 1]:
                TowerIsolator.calibration_factor.append(1.5)
            elif i in [2]:
                TowerIsolator.calibration_factor.append(1)
            elif i in [3, 4, 5, 6, 7, 8, 9, 10]:
                TowerIsolator.calibration_factor.append(0.75)
            else:
                TowerIsolator.calibration_factor.append(0.5)

    @staticmethod
    def extrapolate(eta, pt):
        if math.fabs(eta) > 2.1:
            return TowerIsolator.functions[9].Eval(pt)
        for i, eta_bin in enumerate([0.5, 0.8, 1.1, 1.3, 1.5, 1.7, 1.8, 1.9, 2.0, 2.1]):
            if math.fabs(eta) < eta_bin:
                f = TowerIsolator.functions[i]
                dphi = f.Eval(pt)
                return dphi

    @staticmethod
    def get_twr_ieta(eta):
        aeta = math.fabs(eta)
        for ieta, eta_max in enumerate(TowerIsolator.eta_tower):
            if eta_max > aeta:
                if eta > 0:
                    return ieta+1
                return -ieta

    @staticmethod
    def get_twr_iphi(phi):
        if phi < 0:
            phi += TowerIsolator.pi2
        return int(phi / 0.08726646259971647) + 1
        # iphi = int(phi / 0.08726646259971647)
        # if phi < 0:
        #     iphi = 72 + iphi
        # return iphi

    @staticmethod
    def get_tower_index(ieta, iphi):
        iEtaNoZero = ieta
        if (ieta > 0):
            iEtaNoZero -= 1
        return (iEtaNoZero+28)*72+iphi-1

    @staticmethod
    def get_2x2energy_sums(eta, phi, twrs, radius=2, extrapolate=False, pt=-99, chrg=0):
        ext_phi = phi
        if pt != -99 and extrapolate:
            dphi = TowerIsolator.extrapolate(eta, pt)
            ext_phi += chrg*dphi
            if ext_phi > math.pi:
                ext_phi -= 2*math.pi
            if ext_phi < -math.pi:
                ext_phi += 2*math.pi

        iphi = TowerIsolator.get_twr_iphi(ext_phi)
        ieta = TowerIsolator.get_twr_ieta(eta)

        iphi2x2 = iphi / 2
        ieta2x2 = (ieta + 27) / 2

        foot_print = 0
        et_sum = 0
        for i_twr in range(twrs.n):
            d_ieta = abs(twrs.packedEta[i_twr] - ieta2x2)
            if d_ieta > radius:
                continue

            d_iphi = twrs.packedPhi[i_twr] - iphi2x2
            if d_iphi > 18:
                d_iphi -= 36
            if d_iphi < -18:
                d_iphi += 36
            d_iphi = abs(d_iphi)

            if d_ieta == 0 and d_iphi == 0:
                foot_print = twrs.packedPt[i_twr]

            if d_ieta <= radius and d_iphi <= radius*2:
                et_sum += twrs.packedPt[i_twr]

        return et_sum, foot_print

    @staticmethod
    def get_twr_energy_sums(eta, phi, twrs, radius=5, extrapolate=False, pt=-99, chrg=0):
        ext_phi = phi
        if pt != -99 and extrapolate:
            dphi = TowerIsolator.extrapolate(eta, pt)
            ext_phi += chrg*dphi
            if ext_phi > math.pi:
                ext_phi -= 2*math.pi
            if ext_phi < -math.pi:
                ext_phi += 2*math.pi

        iphi = TowerIsolator.get_twr_iphi(ext_phi)
        ieta = TowerIsolator.get_twr_ieta(eta)

        et_sum = 0
        foot_print = 0

        ieta_min = ieta - radius
        ieta_max = ieta + radius - 1
        if ieta_max > 28:
            ieta_max = 27
        if ieta_min < -28:
            ieta_min = -28

        ietas = []
        for i in range(ieta_min, ieta_max):
            if i >= 0:
                ietas.append(i+1)
            else:
                ietas.append(i)

        iphi_min = iphi - radius*2 + 1
        iphi_max = iphi + radius*2
        iphis = []
        if iphi_min < 1:
            iphis = range(iphi, iphi_max+1)
            for i in range(iphi_min, iphi):
                if i < 1:
                    iphis.append(i + 72)
                else:
                    iphis.append(i)
        if iphi_max > 72:
            iphis = range(iphi_min, iphi+1)
            for i in range(iphi, iphi_max):
                if i > 72:
                    iphis.append(i-72)
                else:
                    iphis.append(i)
        if iphis == []:
            iphis = range(iphi_min, iphi_max+1)

        n_over_zero = 0

        for i_eta in ietas:
            for i_phi in iphis:
                itwr = TowerIsolator.get_tower_index(i_eta, i_phi)
                if not (twrs.packedEta[itwr] == i_eta and twrs.packedPhi[itwr] == i_phi):
                    print "shit", twrs.packedEta[itwr], i_eta, twrs.packedPhi[itwr], i_phi
                d_ieta = twrs.packedEta[itwr] - ieta
                d_iphi = twrs.packedPhi[itwr] - iphi
                if d_iphi > 36:
                    d_iphi -= 72
                if d_iphi < -36:
                    d_iphi += 72

                if (abs(d_iphi) <= radius*2 and abs(d_ieta) <= radius):
                    et_sum += twrs.packedPt[itwr]
                    if twrs.packedPt[itwr] > 0:
                        n_over_zero += 1
                else:
                    print "whoops..."

                if (abs(d_ieta) <= 1 and abs(d_iphi) <= 1):
                    foot_print += twrs.packedPt[itwr]
        # for i_twr in range(twrs.n):
        #     if abs(twrs.packedEta[i_twr]) > 28:
        #         continue
        #     d_ieta = twrs.packedEta[i_twr] - ieta
        #     if abs(d_ieta) > radius:
        #         continue
        #     d_iphi = twrs.packedPhi[i_twr] - iphi
        #     if d_iphi > 36:
        #         d_iphi -= 72
        #     if d_iphi < -36:
        #         d_iphi += 72

        #     if (abs(d_iphi) <= radius*2 and abs(d_ieta) <= radius):
        #         et_sum += twrs.packedPt[i_twr]

        #     if (abs(d_iphi) <= 1 and abs(d_ieta) == 0) or (abs(d_ieta) <= 1 and abs(d_iphi) == 0):
        #         foot_print += twrs.packedPt[i_twr]

        return et_sum, foot_print, n_over_zero
