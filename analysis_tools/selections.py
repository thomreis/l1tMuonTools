from analysis_tools.isolation.caloTowerIso import CaloTowerIsolator
import math
import ROOT as root

class Matcher(object):
    """Class containing static functions for matching L1Analysis collections"""
    twopi = math.pi*2

    @staticmethod
    def norm_phi(phi):
        """ Makes sure that phi is in -pi, pi, implementation stolen from TVector2"""
        nphi = phi
        while (nphi >= math.pi):
            nphi -= Matcher.twopi
        while (nphi < -math.pi):
            nphi += Matcher.twopi
        return nphi

    @staticmethod
    def delta_phi(phi1, phi2):
        return Matcher.norm_phi(phi1 - phi2)

    @staticmethod
    def delta_r(phi1, eta1, phi2, eta2, phi_normalize=True):
        deta = math.fabs(eta1 - eta2)
        if phi_normalize:
            phi1 = Matcher.norm_phi(phi1)
            phi2 = Matcher.norm_phi(phi2)
        dphi = Matcher.delta_phi(phi1, phi2)
        return math.sqrt(deta*deta + dphi*dphi)

    @staticmethod
    def match_dr(eta_coll1, phi_coll1, eta_coll2, phi_coll2, cut=0.5, phi_normalize=True, idcs1=None, idcs2=None):
        """
        Matching based on delta R (dR = sqrt(dEta^2 + dPhi^2)):
        Returns a sorted list of index pairs and corr. dR (first in tuple idx for collection1; second for collection2; third is dR)
        List is sorted by increasing dR
        Parameters eta_coll1/2 & phi_coll1/2 are the eta/phi vectors from the two collections
        Parameter cut specifies the maximum allowed dR
        Parameter phi_normalize specifies whether the phi scales are different (-pi, pi) vs (0, 2pi)
        Parameters idcs1/2 can specify if only a subset of the collections should be considered
        """
        index_tuples = []
        if idcs1 is None:
            idcs1 = range(eta_coll1.size())
        if idcs2 is None:
            idcs2 = range(eta_coll2.size())

        for i in idcs1:
            for j in idcs2:
                deta = eta_coll1[i] - eta_coll2[j]
                abs_deta = math.fabs(deta)
                if abs_deta > cut:  # don't bother with the rest
                    continue
                phi1 = phi_coll1[i]
                phi2 = phi_coll2[j]
                if phi_normalize:
                    phi1 = Matcher.norm_phi(phi1)
                    phi2 = Matcher.norm_phi(phi2)
                dphi = Matcher.delta_phi(phi1, phi2)
                dr = math.sqrt(abs_deta*abs_deta + dphi*dphi)
                if dr < cut:
                    index_tuples.append([i, j, dr, deta, dphi])

        return sorted(index_tuples, key=lambda idx_dr: idx_dr[2])


class MuonSelections(object):
    """Class containing functions for commonly used muon selections"""

    @staticmethod
    def getTfTypeFromTfMuonIdx(idx):
        if idx > 35 and idx < 72:
            return 0
        elif idx > 17 and idx < 90:
            return 1
        elif idx >= 0 and idx < 108:
            return 2
        else:
            return 3


    @staticmethod
    def select_ugmt_muons(ugmt, pt_min=0.5, qual_min=0, qual_max=15, abs_eta_min=0, abs_eta_max=4, bx_min=-1e6, bx_max=1e6, pos_eta=True, neg_eta=True, pos_charge=True, neg_charge=True, tftype=None, idcs=None, useVtxExtraCoord=False):
        type_acc = tftype
        if isinstance(tftype, int):
            type_acc.append(tftype)
        elif tftype is None:
            type_acc = [0, 1, 2]

        indices = []
        if idcs is None:
            idcs = range(ugmt.nMuons)
        for i in idcs:
            if useVtxExtraCoord:
                eta = ugmt.muonEtaAtVtx[i]
            else:
                eta = ugmt.muonEta[i]

            if ugmt.muonEt[i] < pt_min:
                continue
            if not pos_eta and eta >= 0:
                continue
            if not neg_eta and eta < 0:
                continue
            if not pos_charge and ugmt.muonChg[i] > 0:
                continue
            if not neg_charge and ugmt.muonChg[i] < 0:
                continue
            if math.fabs(eta) < abs_eta_min or math.fabs(eta) > abs_eta_max:
                continue
            if ugmt.muonQual[i] < qual_min:
                continue
            if ugmt.muonQual[i] > qual_max:
                continue
            if ugmt.muonBx[i] < bx_min or ugmt.muonBx[i] > bx_max:
                continue
            # find out which TF sent this muon
            if not (MuonSelections.getTfTypeFromTfMuonIdx(ugmt.muonTfMuonIdx[i]) in type_acc):
                continue

            indices.append(i)
        return indices

    @staticmethod
    def select_iso_ugmt_muons(ugmt, caloTowers, iso_min=0., iso_max=1., iso_eta_max=3.0, idcs=None, useVtxExtraCoord=False, iso_type=0):
        indices = []
        if idcs is None:
            idcs = range(ugmt.nMuons)
        for i in idcs:
            if useVtxExtraCoord:
                eta = ugmt.muonEtaAtVtx[i]
            else:
                eta = ugmt.muonEta[i]
            if abs(eta) > iso_eta_max:
                indices.append(i)
            else:
                if iso_type == 0 or iso_type == 1:
                    iso = CaloTowerIsolator.calc_calo_tower_2x2_sum(ugmt, caloTowers, i, (2, 2), maxSum=31) #absolute isolation
                    if iso_type == 1 and ugmt.muonEt[i] > 0.:
                        iso /= float(ugmt.muonEt[i]) #relative isolation
                elif iso_type == 2: #inner cone over threshold
                    iso = CaloTowerIsolator.calc_calo_tower_sum(ugmt, caloTowers, i, (1, 1))
                elif iso_type == 3: #outer cone / total cone below threshold
                    iso = CaloTowerIsolator.calc_out_over_tot_iso(ugmt, caloTowers, i)
                elif iso_type == 4: #inner cone over threshold 2x2
                    iso = CaloTowerIsolator.calc_calo_tower_2x2_sum(ugmt, caloTowers, i, (0, 0), maxSum=31)
                elif iso_type == 5: #outer cone / total cone below threshold 2x2
                    iso = CaloTowerIsolator.calc_out_over_tot_2x2_iso(ugmt, caloTowers, i, maxSum=31)
                #print '{imin} {i} {imax}'.format(imin=iso_min, i=iso, imax=iso_max)
                if iso >= iso_min and iso <= iso_max:
                    indices.append(i)

        return indices

    @staticmethod
    def select_tf_muons(tf, pt_min=0.5, qual_min=0, abs_eta_min=0, abs_eta_max=4, bx_min=-1e6, bx_max=1e6, pos_eta=True, neg_eta=True, pos_charge=True, neg_charge=True, idcs=None):
        ptScale = 0.5
        etaScale = 0.010875
        indices = []
        if idcs is None:
            idcs = range(tf.nTfMuons)
        for i in idcs:
            if tf.tfMuonHwPt[i] * ptScale < pt_min:
                continue
            if not pos_eta and tf.tfMuonHwEta[i] >= 0:
                continue
            if not neg_eta and tf.tfMuonHwEta[i] < 0:
                continue
            if not pos_charge and tf.sign[i] > 0:
                continue
            if not neg_charge and tf.sign[i] < 0:
                continue
            if math.fabs(tf.tfMuonHwEta[i] * etaScale) < abs_eta_min or math.fabs(tf.tfMuonHwEta[i] * etaScale) > abs_eta_max:
                continue
            if tf.tfMuonHwQual[i] < qual_min:
                continue
            if tf.tfMuonBx[i] < bx_min or tf.tfMuonBx[i] > bx_max:
                continue
            indices.append(i)
        return indices

    @staticmethod
    def select_gmt_muons(gmt, pt_min=0.5, qual_min=0, abs_eta_min=0, abs_eta_max=4, bx_min=-1e6, bx_max=1e6, pos_eta=True, neg_eta=True, pos_charge=True, neg_charge=True, idcs=None):
        # qual_min = 8 is interpreted as 2012 running conditions, i.e.:
        # take qualities 6, 7 and 5 if BX == 0
        indices = []

        if idcs is None:
            idcs = range(gmt.N)
        for i in idcs:
            if gmt.Pt[i] < pt_min:
                continue
            if not pos_eta and gmt.Eta[i] >= 0:
                continue
            if not neg_eta and gmt.Eta[i] < 0:
                continue
            if not pos_charge and gmt.Charge[i] > 0:
                continue
            if not neg_charge and gmt.Charge[i] < 0:
                continue
            if math.fabs(gmt.Eta[i]) < abs_eta_min or math.fabs(gmt.Eta[i]) > abs_eta_max:
                continue
            if gmt.Qual[i] < qual_min and qual_min < 8:
                continue
            elif qual_min == 8:
                if (gmt.Qual[i] < 5) or (gmt.Qual[i] == 5 and gmt.CandBx[i] != 0):
                    continue
            if gmt.CandBx[i] < bx_min or gmt.CandBx[i] > bx_max:
                continue
            indices.append(i)
        return indices

    @staticmethod
    def select_reco_muons(reco, pt_min=0.5, pt_max=1.e99, abs_eta_min=0, abs_eta_max=4, pos_eta=True, neg_eta=True, pos_charge=True, neg_charge=True, extrapolated=0, idcs=None):
        if idcs is None:
            idcs = range(reco.nMuons)
        indices = []
        for i in idcs:
            if reco.pt[i] < pt_min:
                continue
            if reco.pt[i] > pt_max:
                continue
            # select eta at vertex or extrapolaed eta at 1st or 2nd muon station
            eta = reco.eta[i]
            phi = reco.phi[i]
            if extrapolated == 1:
                eta = reco.etaSt1[i]
                phi = reco.phiSt1[i]
            elif extrapolated == 2:
                eta = reco.etaSt1[i]
                phi = reco.phiSt1[i]
            if math.fabs(eta) < abs_eta_min or math.fabs(eta) > abs_eta_max:
                continue
            if not pos_eta and eta >= 0:
                continue
            if not neg_eta and eta < 0:
                continue
            if not pos_charge and reco.charge[i] > 0:
                continue
            if not neg_charge and reco.charge[i] < 0:
                continue
            indices.append(i)
        return indices

    @staticmethod
    def select_tag_muons(reco, pt_min=0.5, pt_max=1.e99, abs_eta_min=0, abs_eta_max=4, pos_eta=True, neg_eta=True, pos_charge=True, neg_charge=True, pp_run=True, extrapolated=0, idcs=None):
        indices = []
        reco_indices = MuonSelections.select_reco_muons(reco, pt_min, pt_max, abs_eta_min, abs_eta_max, pos_eta, neg_eta, pos_charge, neg_charge, extrapolated, idcs)
        for i in reco_indices:
            if not reco.isTightMuon[i]:
                continue
            if reco.iso[i] >= 0.15:
                continue
            if pp_run:
                if reco.hlt_isomu[i] != 1:
                    continue
                if reco.hlt_isoDeltaR[i] >= 0.3:
                    continue
            else:
                if reco.hlt_mu[i] != 1:
                    continue
                if reco.hlt_deltaR[i] >= 0.3:
                    continue
            indices.append(i)
        return indices

    @staticmethod
    def select_probe_muons(reco, pt_min=0.5, pt_max=1.e99, abs_eta_min=0, abs_eta_max=4, pos_eta=True, neg_eta=True, pos_charge=True, neg_charge=True, extrapolated=0, idcs=None):
        indices = []
        reco_indices = MuonSelections.select_reco_muons(reco, pt_min, pt_max, abs_eta_min, abs_eta_max, pos_eta, neg_eta, pos_charge, neg_charge, extrapolated, idcs)
        for i in reco_indices:
            if not reco.isTightMuon[i]:
                continue
            if reco.iso[i] >= 0.15:
                continue
            indices.append(i)
        return indices

    @staticmethod
    def select_gen_muons(gen, pt_min=0.5, abs_eta_min=0, abs_eta_max=4, pos_eta=True, neg_eta=True, pos_charge=True, neg_charge=True, idcs=None):
        if idcs is None:
            idcs = range(gen.nPart)
        indices = []
        for i in idcs:
            if abs(gen.partId[i]) != 13:  # Select muons only!
                continue
            if gen.partPt[i] < pt_min:
                continue
            if math.fabs(gen.partEta[i]) < abs_eta_min or math.fabs(gen.partEta[i]) > abs_eta_max:
                continue
            if not pos_eta and gen.partEta[i] >= 0:
                continue
            if not neg_eta and gen.partEta[i] < 0:
                continue
            if not pos_charge and gen.partCh[i] > 0:
                continue
            if not neg_charge and gen.partCh[i] < 0:
                continue
            indices.append(i)
        return indices
