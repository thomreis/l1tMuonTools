from math import floor

class CaloTowerIsolator(object):
    caloTwrEtas = [x*0.087 for x in range(21)]
    caloTwrEtas += [1.83, 1.93, 2.043, 2.172, 2.322, 2.5, 2.65, 3.]
    etaScale = 0.010875
    phiScale = 0.010908
    caloTwrPhiScale = 0.087

    @staticmethod
    def calc_muon_calo_tower_ieta(muIEta):
        # calculate in which caloTower the muon falls in eta
        muEta = CaloTowerIsolator.etaScale * muIEta
        muInCaloTowerIEta = 0
        if abs(muEta) < CaloTowerIsolator.caloTwrEtas[21]:
            muInCaloTowerIEta = int(floor(muIEta / 8.))
            if muIEta >= 0:
                muInCaloTowerIEta += 1
        else:
            for i in range(21,29):
                if abs(muEta) < CaloTowerIsolator.caloTwrEtas[i]:
                    muInCaloTowerIEta = i
                    break
            if muIEta < 0:
                muInCaloTowerIEta *= -1
        if muInCaloTowerIEta > 0:
            muInCaloTowerIEta -= 1
        #print 'muon eta={eta}, ieta={ieta}, caloIEta={cieta}, caloTowerBounds=[{lower}, {upper}]'.format(eta=muEta, ieta=muIEta, cieta=muInCaloTowerIEta, lower=CaloTowerIsolator.caloTwrEtas[abs(muInCaloTowerIEta)-1], upper=CaloTowerIsolator.caloTwrEtas[abs(muInCaloTowerIEta)])
        return muInCaloTowerIEta

    @staticmethod
    def calc_muon_calo_tower_iphi(muIPhi):
        # calculate in which caloTower the muon falls in phi
        muInCaloTowerIPhi = floor(muIPhi * CaloTowerIsolator.phiScale / CaloTowerIsolator.caloTwrPhiScale)
        muInCaloTowerIPhi += 1
        return muInCaloTowerIPhi

    @staticmethod
    def calc_tower_delta_ieta(ieta, twrIEta):
        dIEta = twrIEta - ieta
        if twrIEta > 0:
            dIEta -= 1
        return dIEta

    @staticmethod
    def calc_tower_delta_iphi(iphi, twrIPhi):
        dIPhi = iphi - twrIPhi
        # wrap around
        if dIPhi < -36:
            dIPhi += 72
        elif dIPhi > 35:
            dIPhi -= 72
        return dIPhi

    @staticmethod
    def get_rel_calo_towers(caloTowers, ieta, iphi):
        relTwrs = []
        nCaloTwr = caloTowers.nTower
        for i in range(nCaloTwr):
            dIEta = CaloTowerIsolator.calc_tower_delta_ieta(ieta, caloTowers.ieta[i])
            dIPhi = CaloTowerIsolator.calc_tower_delta_iphi(iphi, caloTowers.iphi[i])
            iEt = caloTowers.iet[i]
            relTwrs.append((dIEta, dIPhi, iEt))
        return relTwrs

    @staticmethod
    def calc_calo_tower_sums(caloTowers, ieta, iphi, radii):
        #calc sums around ieta and iphi
        iEtSums = [0] * len(radii)
        relTwrs = CaloTowerIsolator.get_rel_calo_towers(caloTowers, ieta, iphi)

        for relTwr in relTwrs:
            for i, xyRadii in enumerate(radii):
                if abs(relTwr[0]) <= xyRadii[0] and abs(relTwr[1]) <= xyRadii[1]:
                    iEtSums[i] += relTwr[2]
        return relTwrs, iEtSums
