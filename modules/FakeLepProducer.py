import ROOT
from ROOT import TLorentzVector
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection 
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

import math
import os
import numpy as np
from numpy import sign

class FakeLepProducer(Module):
  def __init__(self , year):
    self.year = year
  def beginJob(self):
    pass
  def endJob(self):
    pass
  def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
    self.out = wrappedOutputTree
    self.out.branch("HLT_passEle32WPTight", "I")
    self.out.branch("n_tight_muon", "I")
    self.out.branch("n_fakeable_muon", "I")
    self.out.branch("n_loose_muon", "I")
    self.out.branch("n_tight_ele", "I")
    self.out.branch("n_fakeable_ele", "I")
    self.out.branch("n_loose_ele", "I")
    self.out.branch("n_tight_jet", "I")
    self.out.branch("n_bjet_DeepB", "I")
    self.out.branch("btag_SFall", "F")
    self.out.branch("n_cjet_DeepB_medium", "I")
    self.out.branch("HT", "F")
    self.out.branch("j1_pt", "F")
    self.out.branch("j1_eta", "F")
    self.out.branch("j1_phi", "F")
    self.out.branch("j1_mass", "F")
    self.out.branch("j2_pt", "F")
    self.out.branch("j2_eta", "F")
    self.out.branch("j2_phi", "F")
    self.out.branch("j2_mass", "F")
    self.out.branch("j3_pt", "F")
    self.out.branch("j3_eta", "F")
    self.out.branch("j3_phi", "F")
    self.out.branch("j3_mass", "F")
    self.out.branch("j4_pt", "F")
    self.out.branch("j4_eta", "F")
    self.out.branch("j4_phi", "F")
    self.out.branch("j4_mass", "F")
    self.out.branch("DeepB_j1_pt", "F")
    self.out.branch("DeepB_j1_eta", "F")
    self.out.branch("DeepB_j1_phi", "F")
    self.out.branch("DeepB_j1_mass", "F")
    self.out.branch("DeepB_j2_pt", "F")
    self.out.branch("DeepB_j2_eta", "F")
    self.out.branch("DeepB_j2_phi", "F")
    self.out.branch("DeepB_j2_mass", "F")
    self.out.branch("DeepC_medium_j1_pt", "F")
    self.out.branch("DeepC_medium_j1_eta", "F")
    self.out.branch("DeepC_medium_j1_phi", "F")
    self.out.branch("DeepC_medium_j1_mass", "F")
    self.out.branch("DeepC_medium_j2_pt", "F")
    self.out.branch("DeepC_medium_j2_eta", "F")
    self.out.branch("DeepC_medium_j2_phi", "F")
    self.out.branch("DeepC_medium_j2_mass", "F")
    self.out.branch("ttc_nl", "B")
    self.out.branch("ttc_jets", "B")
    self.out.branch("ttc_1P1F", "B")
    self.out.branch("ttc_0P2F", "B")
    self.out.branch("ttc_lep1_faketag", "B")
    self.out.branch("ttc_lep2_faketag", "B")
    self.out.branch("ttc_region", "I")
    self.out.branch("ttc_l1_id", "I")
    self.out.branch("ttc_l2_id", "I")
    self.out.branch("ttc_l1_pdgid", "I")
    self.out.branch("ttc_l2_pdgid", "I")
    self.out.branch("ttc_l1_pt", "F")
    self.out.branch("ttc_l1_eta", "F")
    self.out.branch("ttc_l1_phi", "F")
    self.out.branch("ttc_l1_mass", "F")
    self.out.branch("ttc_l2_pt", "F")
    self.out.branch("ttc_l2_eta", "F")
    self.out.branch("ttc_l2_phi", "F")
    self.out.branch("ttc_l2_mass", "F")
    self.out.branch("ttc_mll", "F")
    self.out.branch("ttc_drll", "F")
    self.out.branch("ttc_dphill", "F")
    self.out.branch("ttc_met", "F")
    self.out.branch("ttc_met_phi", "F")
    self.out.branch("ttc_dr_l1j1", "F")
    self.out.branch("ttc_dr_l1j2", "F")
    self.out.branch("ttc_dr_l1j3", "F")
    self.out.branch("ttc_dr_l1j4", "F")
    self.out.branch("ttc_dr_l2j1", "F")
    self.out.branch("ttc_dr_l2j2", "F")
    self.out.branch("ttc_dr_l2j3", "F")
    self.out.branch("ttc_dr_l2j4", "F")
    self.out.branch("ttc_mllj1", "F")
    self.out.branch("ttc_mllj2", "F")
    self.out.branch("ttc_mllj3", "F")
    self.out.branch("ttc_mllj4", "F")
    self.out.branch("WZ_region", "I")
    self.out.branch("WZ_zl1_id", "I")
    self.out.branch("WZ_zl2_id", "I")
    self.out.branch("WZ_wl_id", "I")
    self.out.branch("WZ_zl1_pdgid", "I")
    self.out.branch("WZ_zl2_pdgid", "I")
    self.out.branch("WZ_wl_pdgid", "I")
    self.out.branch("WZ_zl1_pt", "F")
    self.out.branch("WZ_zl1_eta", "F")
    self.out.branch("WZ_zl1_phi", "F")
    self.out.branch("WZ_zl1_mass", "F")
    self.out.branch("WZ_zl2_pt", "F")
    self.out.branch("WZ_zl2_eta", "F")
    self.out.branch("WZ_zl2_phi", "F")
    self.out.branch("WZ_zl2_mass", "F")
    self.out.branch("WZ_l3_pt", "F")
    self.out.branch("WZ_l3_eta", "F")
    self.out.branch("WZ_l3_phi", "F")
    self.out.branch("WZ_l3_mass", "F")
    self.out.branch("WZ_Z_mass", "F")
    self.out.branch("WZ_Z_pt", "F")
    self.out.branch("WZ_Z_eta", "F")
    self.out.branch("WZ_Z_phi", "F")
    self.out.branch("WZ_met", "F")
    self.out.branch("DY_region", "I")
    self.out.branch("DY_1P1F", "B")
    self.out.branch("DY_0P2F", "B")
    self.out.branch("DY_lep1_faketag", "B")
    self.out.branch("DY_lep2_faketag", "B")
    self.out.branch("DY_l1_id", "I")
    self.out.branch("DY_l2_id", "I")
    self.out.branch("DY_l1_pdgid", "I")
    self.out.branch("DY_l2_pdgid", "I")
    self.out.branch("DY_l1_pt", "F")
    self.out.branch("DY_l1_pt_raw", "F")
    self.out.branch("DY_l1_eta", "F")
    self.out.branch("DY_l1_phi", "F")
    self.out.branch("DY_l1_mass", "F")
    self.out.branch("DY_l2_pt", "F")
    self.out.branch("DY_l2_pt_raw", "F")
    self.out.branch("DY_l2_eta", "F")
    self.out.branch("DY_l2_phi", "F")
    self.out.branch("DY_l2_mass", "F")
    self.out.branch("DY_z_mass", "F")
    self.out.branch("DY_z_mass_raw", "F")
    self.out.branch("DY_z_pt", "F")
    self.out.branch("DY_z_pt_raw", "F")
    self.out.branch("DY_z_eta", "F")
    self.out.branch("DY_z_eta_raw", "F")
    self.out.branch("DY_z_phi", "F")
    self.out.branch("DY_z_phi_raw", "F")
    self.out.branch("DY_drll", "F")
    self.out.branch("DY_met", "F")
    self.out.branch("DY_met_phi", "F")
    self.out.branch("tightJets_id_in24","I",lenVar="nJet")
    self.out.branch("tightJets_id_in47","I",lenVar="nJet")
    self.out.branch("tightJets_b_DeepCSVmedium_id","I",lenVar="nJet")
    self.out.branch("tightJets_c_DeepCSVmedium_id","I",lenVar="nJet")
    self.out.branch("tightElectrons_id","I",lenVar="nElectron")
    self.out.branch("fakeable_Electrons_id","I",lenVar="nElectron")
    self.out.branch("additional_vetoElectrons_id","I",lenVar="nElectron")
    self.out.branch("tightMuons_id","I",lenVar="nMuon")
    self.out.branch("fakeable_Muons_id","I",lenVar="nMuon")
    self.out.branch("additional_looseMuons_id","I",lenVar="nMuon")
    self.is_mc = bool(inputTree.GetBranch("GenJet_pt"))

  def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
    pass

  def analyze(self, event):
    
    # PV selection
    if (event.PV_npvsGood<1): return False

    # trigger selection
    # special action for 2017 single ele HLT, https://twiki.cern.ch/twiki/bin/viewauth/CMS/Egamma2017DataRecommendations#Single_Electron_Triggers
    HLT_passEle32WPTight=0
    if self.year=="2017":
      trgobjs=Collection(event, 'TrigObj')
      if event.HLT_Ele32_WPTight_Gsf_L1DoubleEG==1:
	for iobj in range(0,event.nTrigObj):
	  if trgobjs[iobj].id==11 and (trgobjs[iobj].filterBits & (1<<10))== (1<<10):
	    HLT_passEle32WPTight=1

    self.out.fillBranch("HLT_passEle32WPTight",HLT_passEle32WPTight)

    # total number of ele+muon, currently require at least 1 leptons
    if ((event.nMuon + event.nElectron) < 2): return False

    # Muon selection: tight cut-based ID + tight PF iso, or loose cut-based ID + loose PF iso, with pt > 20 GeV
    muons = Collection(event, 'Muon')
    muon_v4_temp=TLorentzVector()
    muon_v4_temp_raw=TLorentzVector()
    tightMuons = []
    tightMuons_raw = []
    tightMuons_pdgid = []
    tightMuons_id = []
    fakeable_Muons = []
    fakeable_Muons_pdgid = []
    fakeable_Muons_id = []
    additional_looseMuons = []
    additional_looseMuons_pdgid = []
    additional_looseMuons_id = []
    for imu in range(0, event.nMuon):
      if (muons[imu].tightId):
        if (muons[imu].pfRelIso04_all<0.15 and abs(muons[imu].eta)<2.4 and muons[imu].tightCharge==2 and event.Muon_corrected_pt[imu]>15):
          muon_v4_temp.SetPtEtaPhiM(event.Muon_corrected_pt[imu], muons[imu].eta, muons[imu].phi, muons[imu].mass)
          muon_v4_temp_raw.SetPtEtaPhiM(muons[imu].pt, muons[imu].eta, muons[imu].phi, muons[imu].mass)
          tightMuons.append(muon_v4_temp.Clone())
          tightMuons_raw.append(muon_v4_temp_raw.Clone())
          tightMuons_pdgid.append(muons[imu].pdgId)
          tightMuons_id.append(imu)
	if (muons[imu].pfRelIso04_all<0.4 and muons[imu].pfRelIso04_all>0.2 and abs(muons[imu].eta)<2.4 and event.Muon_corrected_pt[imu]>15):
	  if self.is_mc and (muons[imu].genPartFlav==1 or muons[imu].genPartFlav==15):
            muon_v4_temp.SetPtEtaPhiM(event.Muon_corrected_pt[imu], muons[imu].eta, muons[imu].phi, muons[imu].mass)
            muon_v4_temp_raw.SetPtEtaPhiM(muons[imu].pt, muons[imu].eta, muons[imu].phi, muons[imu].mass)
            fakeable_Muons.append(muon_v4_temp.Clone())
            fakeable_Muons_pdgid.append(muons[imu].pdgId)
            fakeable_Muons_id.append(imu)
	  if not self.is_mc:
	    muon_v4_temp.SetPtEtaPhiM(event.Muon_corrected_pt[imu], muons[imu].eta, muons[imu].phi, muons[imu].mass)
            muon_v4_temp_raw.SetPtEtaPhiM(muons[imu].pt, muons[imu].eta, muons[imu].phi, muons[imu].mass)
            fakeable_Muons.append(muon_v4_temp.Clone())
            fakeable_Muons_pdgid.append(muons[imu].pdgId)
            fakeable_Muons_id.append(imu)
      elif (muons[imu].looseId):
        if (muons[imu].pfRelIso04_all<0.25 and abs(muons[imu].eta)<2.4 and event.Muon_corrected_pt[imu]>15):
          muon_v4_temp.SetPtEtaPhiM(event.Muon_corrected_pt[imu], muons[imu].eta, muons[imu].phi, muons[imu].mass)
          additional_looseMuons.append(muon_v4_temp.Clone())
          additional_looseMuons_pdgid.append(muons[imu].pdgId)
          additional_looseMuons_id.append(imu)

    n_tight_muon = len(tightMuons)
    n_fakeable_muon = len(fakeable_Muons)
    n_loose_muon = len(additional_looseMuons)

    self.out.fillBranch("n_tight_muon", n_tight_muon)
    self.out.fillBranch("n_fakeable_muon", n_fakeable_muon)
    self.out.fillBranch("n_loose_muon", n_loose_muon)
    if event.nMuon>0:
      tightMuons_id.extend(np.zeros(event.nMuon-len(tightMuons_id),int)-1)
      fakeable_Muons_id.extend(np.zeros(event.nMuon-len(fakeable_Muons_id),int)-1)
      additional_looseMuons_id.extend(np.zeros(event.nMuon-len(additional_looseMuons_id),int)-1)
      self.out.fillBranch("tightMuons_id", tightMuons_id)
      self.out.fillBranch("fakeable_Muons_id", fakeable_Muons_id)
      self.out.fillBranch("additional_looseMuons_id", additional_looseMuons_id)

    # electron selection: tight (veto) cut-based ID + impact parameter cut, with pt > 15 GeV
    electrons = Collection(event, 'Electron')
    electron_v4_temp=TLorentzVector()
    electron_v4_temp_raw=TLorentzVector()
    tightElectrons = []
    tightElectrons_raw = []
    tightElectrons_pdgid = []
    tightElectrons_id = []
    fakeable_Electrons = []
    fakeable_Electrons_pdgid = []
    fakeable_Electrons_id = []
    additional_vetoElectrons = []
    additional_vetoElectrons_pdgid = []
    additional_vetoElectrons_id = []
    for iele in range(0, event.nElectron):
      if not (((abs(electrons[iele].eta+electrons[iele].deltaEtaSC) <1.4442 and abs(electrons[iele].dxy)<0.05 and abs(electrons[iele].dz)<0.1) or (abs(electrons[iele].eta + electrons[iele].deltaEtaSC)>1.566 and abs(electrons[iele].eta + electrons[iele].deltaEtaSC)<2.4 and abs(electrons[iele].dxy)<0.1 and abs(electrons[iele].dz)<0.2)) and electrons[iele].pt>15): continue
      if (electrons[iele].cutBased==4 and electrons[iele].tightCharge==2):
        electron_v4_temp.SetPtEtaPhiM(electrons[iele].pt, electrons[iele].eta, electrons[iele].phi, electrons[iele].mass)
        electron_v4_temp_raw.SetPtEtaPhiM(electrons[iele].pt/electrons[iele].eCorr, electrons[iele].eta, electrons[iele].phi, electrons[iele].mass/electrons[iele].eCorr)
        tightElectrons.append(electron_v4_temp.Clone())
        tightElectrons_raw.append(electron_v4_temp_raw.Clone())
        tightElectrons_pdgid.append(electrons[iele].pdgId)
        tightElectrons_id.append(iele)
      if (electrons[iele].cutBased>1 and electrons[iele].cutBased<4):
        electron_v4_temp.SetPtEtaPhiM(electrons[iele].pt, electrons[iele].eta, electrons[iele].phi, electrons[iele].mass)
        additional_vetoElectrons.append(electron_v4_temp.Clone())
        additional_vetoElectrons_pdgid.append(electrons[iele].pdgId)
        additional_vetoElectrons_id.append(iele)
      if (electrons[iele].cutBased==1):
	if self.is_mc and (electrons[iele].genPartFlav==1 or electrons[iele].genPartFlav==15):
          electron_v4_temp.SetPtEtaPhiM(electrons[iele].pt, electrons[iele].eta, electrons[iele].phi, electrons[iele].mass)
          fakeable_Electrons.append(electron_v4_temp.Clone())
          fakeable_Electrons_pdgid.append(electrons[iele].pdgId)
          fakeable_Electrons_id.append(iele)
	if not self.is_mc:
	  electron_v4_temp.SetPtEtaPhiM(electrons[iele].pt, electrons[iele].eta, electrons[iele].phi, electrons[iele].mass)
          fakeable_Electrons.append(electron_v4_temp.Clone())
          fakeable_Electrons_pdgid.append(electrons[iele].pdgId)
          fakeable_Electrons_id.append(iele)

    n_tight_ele = len(tightElectrons)
    n_fakeable_ele = len(fakeable_Electrons)
    n_loose_ele = len(additional_vetoElectrons)
    self.out.fillBranch("n_tight_ele", n_tight_ele)
    self.out.fillBranch("n_fakeable_ele", n_fakeable_ele)
    self.out.fillBranch("n_loose_ele", n_loose_ele)
    if event.nElectron>0:
      tightElectrons_id.extend(np.zeros(event.nElectron-len(tightElectrons_id),int)-1)
      fakeable_Electrons_id.extend(np.zeros(event.nElectron-len(fakeable_Electrons_id),int)-1)
      additional_vetoElectrons_id.extend(np.zeros(event.nElectron-len(additional_vetoElectrons_id),int)-1)
      self.out.fillBranch("tightElectrons_id", tightElectrons_id)
      self.out.fillBranch("fakeable_Electrons_id", fakeable_Electrons_id)
      self.out.fillBranch("additional_vetoElectrons_id", additional_vetoElectrons_id)

    # tight leptons and additional loose leptons collection
    tightLeptons = tightMuons + tightElectrons
    tightLeptons.sort(key=lambda x: x.Pt(), reverse=True)
    tightLeptons_raw = tightMuons_raw + tightElectrons_raw
    tightLeptons_raw.sort(key=lambda x: x.Pt(), reverse=True)
    fakeableLeptons = fakeable_Muons + fakeable_Electrons
    fakeableLeptons.sort(key=lambda x: x.Pt(), reverse=True)
    looseLeptons = additional_looseMuons + additional_vetoElectrons
    looseLeptons.sort(key=lambda x: x.Pt(), reverse=True)

    # https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL17
    # https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL18
    # https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD#Jets
    # tightLepVeto PF jets (ak4), UL 2016/2017/2018 (jetId 110=6), medium B-tag WP
    # UL17 DeepCSV=(nanoaod btagDeepB) loose: 0.1355, medium: 0.4506, tight: 0.7738
    # UL18 DeepCSV=(nanoaod btagDeepB) loose: 0.1208, medium: 0.4168, tight: 0.7665
    # UL17 DeepFlavor=(nanoaod btagDeepFlavB) loose: 0.0532, medium: 0.3040, tight: 0.7476
    # UL18 DeepFlavor=(nanoaod btagDeepFlavB) loose: 0.0490, medium: 0.2783, tight: 0.7100

    # c-jet tag is based on two-D cuts, medium DeepJet WP:
    # UL17 CvsL=btagDeepFlavCvL: 0.085, CvsB=btagDeepFlavCvB: 0.34
    # UL18 CvsL=btagDeepFlavCvL: 0.099, CvsB=btagDeepFlavCvB: 0.325
    # c-tag not available in NANOAOD yet

    jets = Collection(event, 'Jet')

    j1_pt=-99
    j1_eta=-99
    j1_phi=-99
    j1_mass=-99
    j2_pt=-99
    j2_eta=-99
    j2_phi=-99
    j2_mass=-99
    j3_pt=-99
    j3_eta=-99
    j3_phi=-99
    j3_mass=-99
    j4_pt=-99
    j4_eta=-99
    j4_phi=-99
    j4_mass=-99
    DeepB_j1_pt=-99
    DeepB_j1_eta=-99
    DeepB_j1_phi=-99
    DeepB_j1_mass=-99
    DeepB_j2_pt=-99
    DeepB_j2_eta=-99
    DeepB_j2_phi=-99
    DeepB_j2_mass=-99
    DeepC_medium_j1_pt=-99
    DeepC_medium_j1_eta=-99
    DeepC_medium_j1_phi=-99
    DeepC_medium_j1_mass=-99
    DeepC_medium_j2_pt=-99
    DeepC_medium_j2_eta=-99
    DeepC_medium_j2_phi=-99
    DeepC_medium_j2_mass=-99

    tightJets_id_in24 = []
    tightJets_id_in47 = []

    tightJets_b_DeepCSVmedium_id = []

    tightJets_c_DeepCSVmedium_id = []

    # require DeltaR between Jets and tight leptons greater than 0.4
    jet_v4_temp=TLorentzVector()
    for ijet in range(0, event.nJet):
      pass_jet_lep_Dr=1
      jet_v4_temp.SetPtEtaPhiM(event.Jet_pt_nom[ijet],event.Jet_eta[ijet],event.Jet_phi[ijet],event.Jet_mass_nom[ijet])
      for ilep in range(0,len(tightLeptons)):
	if jet_v4_temp.DeltaR(tightLeptons[ilep])<0.4:pass_jet_lep_Dr=0
      for ilep in range(0,len(fakeableLeptons)):
	if jet_v4_temp.DeltaR(fakeableLeptons[ilep])<0.4:pass_jet_lep_Dr=0

      if not (pass_jet_lep_Dr>0):continue
      if self.year=="2016":
        if jets[ijet].jetId==6 and event.Jet_pt_nom[ijet]>30:
	  if abs(jets[ijet].eta)<4.7 and abs(jets[ijet].eta)>=2.4: 
	    tightJets_id_in47.append(ijet)
	  if abs(jets[ijet].eta)<2.4:
	    tightJets_id_in24.append(ijet)

            if (jets[ijet].btagDeepB > 0.4941):
              tightJets_b_DeepCSVmedium_id.append(ijet)

      if (self.year=="2017"):
	if jets[ijet].jetId==6 and event.Jet_pt_nom[ijet]>30:
	  if abs(jets[ijet].eta)<4.7 and abs(jets[ijet].eta)>=2.4:
	    tightJets_id_in47.append(ijet)
	  if abs(jets[ijet].eta)<2.4:
            tightJets_id_in24.append(ijet)

            if (jets[ijet].btagDeepFlavB > 0.3040):
              tightJets_b_DeepCSVmedium_id.append(ijet)

      if (self.year=="2018"):
	if jets[ijet].jetId==6 and event.Jet_pt_nom[ijet]>30:
	  if abs(jets[ijet].eta)<4.7 and abs(jets[ijet].eta)>=2.4:
	    tightJets_id_in47.append(ijet)
	  if abs(jets[ijet].eta)<2.4:
            tightJets_id_in24.append(ijet)

            if (jets[ijet].btagDeepFlavB > 0.2783):
              tightJets_b_DeepCSVmedium_id.append(ijet)

    HT=0
    for ijet in range(0,len(tightJets_id_in24)):
      HT=HT+event.Jet_pt_nom[tightJets_id_in24[ijet]]
    self.out.fillBranch("HT",HT)

    for ijet in range(0, event.nJet):
      if not (ijet in tightJets_id_in24):continue
      if (ijet in tightJets_b_DeepCSVmedium_id):continue
      #mudium WP
      if self.year=="2017" and jets[ijet].btagDeepFlavCvL>0.085 and jets[ijet].btagDeepFlavCvB<0.34:
	tightJets_c_DeepCSVmedium_id.append(ijet)
      if self.year=="2018" and jets[ijet].btagDeepFlavCvL>0.099 and jets[ijet].btagDeepFlavCvB<0.325:
	tightJets_c_DeepCSVmedium_id.append(ijet)
      
    n_tight_jet = len(tightJets_id_in24)
    n_bjet_DeepB = len(tightJets_b_DeepCSVmedium_id)

    btag_SFall=1.
    if n_bjet_DeepB>0 and self.is_mc:
      for ib in range(0,n_bjet_DeepB):
        btag_SFall=btag_SFall*event.Jet_btagSF_deepjet_M[tightJets_b_DeepCSVmedium_id[ib]]

    n_cjet_DeepB_medium = len(tightJets_c_DeepCSVmedium_id)
    self.out.fillBranch("n_tight_jet",n_tight_jet)
    self.out.fillBranch("n_bjet_DeepB",n_bjet_DeepB)
    self.out.fillBranch("btag_SFall",btag_SFall)
    self.out.fillBranch("n_cjet_DeepB_medium",n_cjet_DeepB_medium)

    if n_tight_jet>3:
      j4_pt=event.Jet_pt_nom[tightJets_id_in24[3]]
      j4_eta=event.Jet_eta[tightJets_id_in24[3]]
      j4_phi=event.Jet_phi[tightJets_id_in24[3]]
      j4_mass=event.Jet_mass_nom[tightJets_id_in24[3]]
      j3_pt=event.Jet_pt_nom[tightJets_id_in24[2]]
      j3_eta=event.Jet_eta[tightJets_id_in24[2]]
      j3_phi=event.Jet_phi[tightJets_id_in24[2]]
      j3_mass=event.Jet_mass_nom[tightJets_id_in24[2]]
      j2_pt=event.Jet_pt_nom[tightJets_id_in24[1]]
      j2_eta=event.Jet_eta[tightJets_id_in24[1]]
      j2_phi=event.Jet_phi[tightJets_id_in24[1]]
      j2_mass=event.Jet_mass_nom[tightJets_id_in24[1]]
      j1_pt=event.Jet_pt_nom[tightJets_id_in24[0]]
      j1_eta=event.Jet_eta[tightJets_id_in24[0]]
      j1_phi=event.Jet_phi[tightJets_id_in24[0]]
      j1_mass=event.Jet_mass_nom[tightJets_id_in24[0]]
    if n_tight_jet==3:
      j3_pt=event.Jet_pt_nom[tightJets_id_in24[2]]
      j3_eta=event.Jet_eta[tightJets_id_in24[2]]
      j3_phi=event.Jet_phi[tightJets_id_in24[2]]
      j3_mass=event.Jet_mass_nom[tightJets_id_in24[2]]
      j2_pt=event.Jet_pt_nom[tightJets_id_in24[1]]
      j2_eta=event.Jet_eta[tightJets_id_in24[1]]
      j2_phi=event.Jet_phi[tightJets_id_in24[1]]
      j2_mass=event.Jet_mass_nom[tightJets_id_in24[1]]
      j1_pt=event.Jet_pt_nom[tightJets_id_in24[0]]
      j1_eta=event.Jet_eta[tightJets_id_in24[0]]
      j1_phi=event.Jet_phi[tightJets_id_in24[0]]
      j1_mass=event.Jet_mass_nom[tightJets_id_in24[0]]
    if n_tight_jet==2:
      j2_pt=event.Jet_pt_nom[tightJets_id_in24[1]]
      j2_eta=event.Jet_eta[tightJets_id_in24[1]]
      j2_phi=event.Jet_phi[tightJets_id_in24[1]]
      j2_mass=event.Jet_mass_nom[tightJets_id_in24[1]]
      j1_pt=event.Jet_pt_nom[tightJets_id_in24[0]]
      j1_eta=event.Jet_eta[tightJets_id_in24[0]]
      j1_phi=event.Jet_phi[tightJets_id_in24[0]]
      j1_mass=event.Jet_mass_nom[tightJets_id_in24[0]]
    if n_tight_jet==1:
      j1_pt=event.Jet_pt_nom[tightJets_id_in24[0]]
      j1_eta=event.Jet_eta[tightJets_id_in24[0]]
      j1_phi=event.Jet_phi[tightJets_id_in24[0]]
      j1_mass=event.Jet_mass_nom[tightJets_id_in24[0]]

    if n_bjet_DeepB>1:
      DeepB_j1_pt=event.Jet_pt_nom[tightJets_b_DeepCSVmedium_id[0]]
      DeepB_j1_eta=event.Jet_eta[tightJets_b_DeepCSVmedium_id[0]]
      DeepB_j1_phi=event.Jet_phi[tightJets_b_DeepCSVmedium_id[0]]
      DeepB_j1_mass=event.Jet_mass_nom[tightJets_b_DeepCSVmedium_id[0]]
      DeepB_j2_pt=event.Jet_pt_nom[tightJets_b_DeepCSVmedium_id[1]]
      DeepB_j2_eta=event.Jet_eta[tightJets_b_DeepCSVmedium_id[1]]
      DeepB_j2_phi=event.Jet_phi[tightJets_b_DeepCSVmedium_id[1]]
      DeepB_j2_mass=event.Jet_mass_nom[tightJets_b_DeepCSVmedium_id[1]]
    if n_bjet_DeepB==1:
      DeepB_j1_pt=event.Jet_pt_nom[tightJets_b_DeepCSVmedium_id[0]]
      DeepB_j1_eta=event.Jet_eta[tightJets_b_DeepCSVmedium_id[0]]
      DeepB_j1_phi=event.Jet_phi[tightJets_b_DeepCSVmedium_id[0]]
      DeepB_j1_mass=event.Jet_mass_nom[tightJets_b_DeepCSVmedium_id[0]]

    if n_cjet_DeepB_medium>1:
      DeepC_medium_j1_pt=event.Jet_pt_nom[tightJets_c_DeepCSVmedium_id[0]]
      DeepC_medium_j1_eta=event.Jet_eta[tightJets_c_DeepCSVmedium_id[0]]
      DeepC_medium_j1_phi=event.Jet_phi[tightJets_c_DeepCSVmedium_id[0]]
      DeepC_medium_j1_mass=event.Jet_mass_nom[tightJets_c_DeepCSVmedium_id[0]]
      DeepC_medium_j2_pt=event.Jet_pt_nom[tightJets_c_DeepCSVmedium_id[1]]
      DeepC_medium_j2_eta=event.Jet_eta[tightJets_c_DeepCSVmedium_id[1]]
      DeepC_medium_j2_phi=event.Jet_phi[tightJets_c_DeepCSVmedium_id[1]]
      DeepC_medium_j2_mass=event.Jet_mass_nom[tightJets_c_DeepCSVmedium_id[1]]
    if n_cjet_DeepB_medium==1:
      DeepC_medium_j1_pt=event.Jet_pt_nom[tightJets_c_DeepCSVmedium_id[0]]
      DeepC_medium_j1_eta=event.Jet_eta[tightJets_c_DeepCSVmedium_id[0]]
      DeepC_medium_j1_phi=event.Jet_phi[tightJets_c_DeepCSVmedium_id[0]]
      DeepC_medium_j1_mass=event.Jet_mass_nom[tightJets_c_DeepCSVmedium_id[0]]

    self.out.fillBranch("j1_pt",j1_pt)
    self.out.fillBranch("j1_eta",j1_eta)
    self.out.fillBranch("j1_phi",j1_phi)
    self.out.fillBranch("j1_mass",j1_mass)
    self.out.fillBranch("j2_pt",j2_pt)
    self.out.fillBranch("j2_eta",j2_eta)
    self.out.fillBranch("j2_phi",j2_phi)
    self.out.fillBranch("j2_mass",j2_mass)
    self.out.fillBranch("j3_pt",j3_pt)
    self.out.fillBranch("j3_eta",j3_eta)
    self.out.fillBranch("j3_phi",j3_phi)
    self.out.fillBranch("j3_mass",j3_mass)
    self.out.fillBranch("j4_pt",j4_pt)
    self.out.fillBranch("j4_eta",j4_eta)
    self.out.fillBranch("j4_phi",j4_phi)
    self.out.fillBranch("j4_mass",j4_mass)
    self.out.fillBranch("DeepB_j1_pt",DeepB_j1_pt)
    self.out.fillBranch("DeepB_j1_eta",DeepB_j1_eta)
    self.out.fillBranch("DeepB_j1_phi",DeepB_j1_phi)
    self.out.fillBranch("DeepB_j1_mass",DeepB_j1_mass)
    self.out.fillBranch("DeepB_j2_pt",DeepB_j2_pt)
    self.out.fillBranch("DeepB_j2_eta",DeepB_j2_eta)
    self.out.fillBranch("DeepB_j2_phi",DeepB_j2_phi)
    self.out.fillBranch("DeepB_j2_mass",DeepB_j2_mass)
    self.out.fillBranch("DeepC_medium_j1_pt",DeepC_medium_j1_pt)
    self.out.fillBranch("DeepC_medium_j1_eta",DeepC_medium_j1_eta)
    self.out.fillBranch("DeepC_medium_j1_phi",DeepC_medium_j1_phi)
    self.out.fillBranch("DeepC_medium_j1_mass",DeepC_medium_j1_mass)
    self.out.fillBranch("DeepC_medium_j2_pt",DeepC_medium_j2_pt)
    self.out.fillBranch("DeepC_medium_j2_eta",DeepC_medium_j2_eta)
    self.out.fillBranch("DeepC_medium_j2_phi",DeepC_medium_j2_phi)
    self.out.fillBranch("DeepC_medium_j2_mass",DeepC_medium_j2_mass)

    tightJets_id_in24.extend(np.zeros(event.nJet-len(tightJets_id_in24),int)-1)
    tightJets_id_in47.extend(np.zeros(event.nJet-len(tightJets_id_in47),int)-1)
    tightJets_b_DeepCSVmedium_id.extend(np.zeros(event.nJet-len(tightJets_b_DeepCSVmedium_id),int)-1)
    tightJets_c_DeepCSVmedium_id.extend(np.zeros(event.nJet-len(tightJets_c_DeepCSVmedium_id),int)-1)

    self.out.fillBranch("tightJets_id_in24",tightJets_id_in24)
    self.out.fillBranch("tightJets_id_in47",tightJets_id_in47)
    self.out.fillBranch("tightJets_b_DeepCSVmedium_id",tightJets_b_DeepCSVmedium_id)
    self.out.fillBranch("tightJets_c_DeepCSVmedium_id",tightJets_c_DeepCSVmedium_id)

    if (len(tightLeptons)+len(fakeableLeptons))<2:return False

    #    t     t      cccc
    #  ttttt ttttt   cc
    #    t     t    cc
    #    t     t     cc
    #    ttt   ttt    cccc
    #region: only 2 tight leptons, at least three jets, 2 b-jet and 1 cjet (c-tag not available yet), mll>20, |Z-91.1876|>15 in 2 ele case
    #ttc region lepton number selections
    ttc_nl=False
    #ttc region jet and bjet selection
    ttc_jets=False
    ttc_1P1F=False
    ttc_0P2F=False
    ttc_lep1_faketag=False
    ttc_lep2_faketag=False
    #ttc region tag, 1:2 muon, 2:1 muon, 3:0 muon
    ttc_region=0
    ttc_l1_id=-1
    ttc_l2_id=-1
    ttc_l1_pdgid=-99
    ttc_l2_pdgid=-99
    ttc_l1_pt=-99
    ttc_l1_eta=-99
    ttc_l1_phi=-99
    ttc_l1_mass=-99
    ttc_l2_pt=-99
    ttc_l2_eta=-99
    ttc_l2_phi=-99
    ttc_l2_mass=-99
    ttc_mll=-99
    ttc_drll=-99
    ttc_dphill=-99
    ttc_met=-99
    ttc_met_phi=-99
    ttc_dr_l1j1=-99
    ttc_dr_l1j2=-99
    ttc_dr_l1j3=-99
    ttc_dr_l1j4=-99
    ttc_dr_l2j1=-99
    ttc_dr_l2j2=-99
    ttc_dr_l2j3=-99
    ttc_dr_l2j4=-99
    ttc_mllj1=-99
    ttc_mllj2=-99
    ttc_mllj3=-99
    ttc_mllj4=-99
    
    # the two leptons with pt 20, 3th lepton veto
    if ((len(fakeableLeptons)==2 and fakeableLeptons[1].Pt()>20) or (len(tightLeptons)==1 and len(fakeableLeptons)==1 and tightLeptons[0].Pt()>20 and fakeableLeptons[0].Pt()>20)) and len(looseLeptons)==0:
      ttc_nl=True
    # at least three jets
    if ttc_nl and n_tight_jet>2:
      ttc_jets=True

    if ttc_nl:
      if self.is_mc:
        ttc_met=event.MET_T1Smear_pt
	ttc_met_phi=event.MET_T1Smear_phi
      else:
        ttc_met=event.MET_T1_pt
        ttc_met_phi=event.MET_T1_phi

      if len(fakeableLeptons)==2:
	ttc_0P2F=True
	ttc_lep1_faketag=True
	ttc_lep2_faketag=True

        if len(fakeable_Electrons)==0 and abs(fakeable_Muons_pdgid[0]+fakeable_Muons_pdgid[1])==26:
  	  ttc_region=1
  	  ttc_l1_id=fakeable_Muons_id[0]
  	  ttc_l2_id=fakeable_Muons_id[1]
  	  ttc_l1_pdgid=fakeable_Muons_pdgid[0]
  	  ttc_l2_pdgid=fakeable_Muons_pdgid[1]
  	  ttc_l1_pt=fakeable_Muons[0].Pt()
  	  ttc_l1_eta=fakeable_Muons[0].Eta()
  	  ttc_l1_phi=fakeable_Muons[0].Phi()
  	  ttc_l1_mass=fakeable_Muons[0].M()
  	  ttc_l2_pt=fakeable_Muons[1].Pt()
  	  ttc_l2_eta=fakeable_Muons[1].Eta()
  	  ttc_l2_phi=fakeable_Muons[1].Phi()
  	  ttc_l2_mass=fakeable_Muons[1].M()
  	  ttc_mll=(fakeable_Muons[0]+fakeable_Muons[1]).M()
  	  ttc_drll=fakeable_Muons[0].DeltaR(fakeable_Muons[1])
  	  ttc_dphill=fakeable_Muons[0].DeltaPhi(fakeable_Muons[1])

        if len(fakeable_Electrons)==1 and abs(fakeable_Muons_pdgid[0]+fakeable_Electrons_pdgid[0])==24:
  	  ttc_region=2
  	  ttc_l1_id=fakeable_Muons_id[0]
  	  ttc_l1_pdgid=fakeable_Muons_pdgid[0]
  	  ttc_l2_id=fakeable_Electrons_id[0]
  	  ttc_l2_pdgid=fakeable_Electrons_pdgid[0]
  	  ttc_l1_pt=fakeable_Muons[0].Pt()
          ttc_l1_eta=fakeable_Muons[0].Eta()
          ttc_l1_phi=fakeable_Muons[0].Phi()
          ttc_l1_mass=fakeable_Muons[0].M()
          ttc_l2_pt=fakeable_Electrons[0].Pt()
          ttc_l2_eta=fakeable_Electrons[0].Eta()
          ttc_l2_phi=fakeable_Electrons[0].Phi()
          ttc_l2_mass=fakeable_Electrons[0].M()
          ttc_mll=(fakeable_Muons[0]+fakeable_Electrons[0]).M()
          ttc_drll=fakeable_Muons[0].DeltaR(fakeable_Electrons[0])
          ttc_dphill=fakeable_Muons[0].DeltaPhi(fakeable_Electrons[0])

        if len(fakeable_Electrons)==2 and abs(fakeable_Electrons_pdgid[0]+fakeable_Electrons_pdgid[1])==22:
  	  ttc_region=3
  	  ttc_l1_id=fakeable_Electrons_id[0]
  	  ttc_l2_id=fakeable_Electrons_id[1]
  	  ttc_l1_pdgid=fakeable_Electrons_pdgid[0]
  	  ttc_l2_pdgid=fakeable_Electrons_pdgid[1]
  	  ttc_l1_pt=fakeable_Electrons[0].Pt()
          ttc_l1_eta=fakeable_Electrons[0].Eta()
          ttc_l1_phi=fakeable_Electrons[0].Phi()
          ttc_l1_mass=fakeable_Electrons[0].M()
          ttc_l2_pt=fakeable_Electrons[1].Pt()
          ttc_l2_eta=fakeable_Electrons[1].Eta()
          ttc_l2_phi=fakeable_Electrons[1].Phi()
          ttc_l2_mass=fakeable_Electrons[1].M()
          ttc_mll=(fakeable_Electrons[0]+fakeable_Electrons[1]).M()
          ttc_drll=fakeable_Electrons[0].DeltaR(fakeable_Electrons[1])
          ttc_dphill=fakeable_Electrons[0].DeltaPhi(fakeable_Electrons[1])
	
      # one real and one fake lepton
      if len(fakeableLeptons)==1:
	ttc_1P1F=True

	#one real muon
	if len(tightElectrons)==0:
	  # one fake muon
	  if len(fakeable_Electrons)==0 and abs(tightMuons_pdgid[0]+fakeable_Muons_pdgid[0])==26:
  	    ttc_region=1
	    if tightMuons[0].Pt()>fakeable_Muons[0].Pt():
              ttc_lep1_faketag=False
              ttc_lep2_faketag=True
  	      ttc_l1_id=tightMuons_id[0]
      	      ttc_l2_id=fakeable_Muons_id[0]
	      ttc_l1_pdgid=tightMuons_pdgid[0]
	      ttc_l2_pdgid=fakeable_Muons_pdgid[0]
	      ttc_l1_pt=tightMuons[0].Pt()
	      ttc_l1_eta=tightMuons[0].Eta()
	      ttc_l1_phi=tightMuons[0].Phi()
	      ttc_l1_mass=tightMuons[0].M()
	      ttc_l2_pt=fakeable_Muons[0].Pt()
	      ttc_l2_eta=fakeable_Muons[0].Eta()
	      ttc_l2_phi=fakeable_Muons[0].Phi()
	      ttc_l2_mass=fakeable_Muons[0].M()
	    else:
              ttc_lep1_faketag=True
              ttc_lep2_faketag=False
  	      ttc_l2_id=tightMuons_id[0]
      	      ttc_l1_id=fakeable_Muons_id[0]
	      ttc_l2_pdgid=tightMuons_pdgid[0]
	      ttc_l1_pdgid=fakeable_Muons_pdgid[0]
	      ttc_l2_pt=tightMuons[0].Pt()
	      ttc_l2_eta=tightMuons[0].Eta()
	      ttc_l2_phi=tightMuons[0].Phi()
	      ttc_l2_mass=tightMuons[0].M()
	      ttc_l1_pt=fakeable_Muons[0].Pt()
	      ttc_l1_eta=fakeable_Muons[0].Eta()
	      ttc_l1_phi=fakeable_Muons[0].Phi()
	      ttc_l1_mass=fakeable_Muons[0].M()
      	    ttc_mll=(fakeable_Muons[0]+tightMuons[0]).M()
      	    ttc_drll=fakeable_Muons[0].DeltaR(tightMuons[0])
      	    ttc_dphill=fakeable_Muons[0].DeltaPhi(tightMuons[0])

	  # one fake electron
	  if len(fakeable_Electrons)==1 and abs(tightMuons_pdgid[0]+fakeable_Electrons_pdgid[0])==24:
  	    ttc_region=2
            ttc_lep1_faketag=False
            ttc_lep2_faketag=True
  	    ttc_l1_id=tightMuons_id[0]
      	    ttc_l2_id=fakeable_Electrons_id[0]
	    ttc_l1_pdgid=tightMuons_pdgid[0]
	    ttc_l2_pdgid=fakeable_Electrons_pdgid[0]
	    ttc_l1_pt=tightMuons[0].Pt()
	    ttc_l1_eta=tightMuons[0].Eta()
	    ttc_l1_phi=tightMuons[0].Phi()
	    ttc_l1_mass=tightMuons[0].M()
	    ttc_l2_pt=fakeable_Electrons[0].Pt()
	    ttc_l2_eta=fakeable_Electrons[0].Eta()
	    ttc_l2_phi=fakeable_Electrons[0].Phi()
	    ttc_l2_mass=fakeable_Electrons[0].M()
      	    ttc_mll=(fakeable_Electrons[0]+tightMuons[0]).M()
      	    ttc_drll=fakeable_Electrons[0].DeltaR(tightMuons[0])
      	    ttc_dphill=fakeable_Electrons[0].DeltaPhi(tightMuons[0])

	#one real electron
	if len(tightElectrons)==1:
	  # one fake muon
	  if len(fakeable_Electrons)==0 and abs(tightElectrons_pdgid[0]+fakeable_Muons_pdgid[0])==24:
  	    ttc_region=2
            ttc_lep1_faketag=True
            ttc_lep2_faketag=False
  	    ttc_l1_id=fakeable_Muons_id[0]
      	    ttc_l2_id=tightElectrons_id[0]
	    ttc_l1_pdgid=fakeable_Muons_pdgid[0]
	    ttc_l2_pdgid=tightElectrons_pdgid[0]
	    ttc_l1_pt=fakeable_Muons[0].Pt()
	    ttc_l1_eta=fakeable_Muons[0].Eta()
	    ttc_l1_phi=fakeable_Muons[0].Phi()
	    ttc_l1_mass=fakeable_Muons[0].M()
	    ttc_l2_pt=tightElectrons[0].Pt()
	    ttc_l2_eta=tightElectrons[0].Eta()
	    ttc_l2_phi=tightElectrons[0].Phi()
	    ttc_l2_mass=tightElectrons[0].M()
      	    ttc_mll=(fakeable_Muons[0]+tightElectrons[0]).M()
      	    ttc_drll=fakeable_Muons[0].DeltaR(tightElectrons[0])
      	    ttc_dphill=fakeable_Muons[0].DeltaPhi(tightElectrons[0])

	  # one fake electron
	  if len(fakeable_Electrons)==1 and abs(tightElectrons_pdgid[0]+fakeable_Electrons_pdgid[0])==22:
  	    ttc_region=3
	    if tightElectrons[0].Pt()>fakeable_Electrons[0].Pt():
              ttc_lep1_faketag=False
              ttc_lep2_faketag=True
  	      ttc_l1_id=tightElectrons_id[0]
      	      ttc_l2_id=fakeable_Electrons_id[0]
	      ttc_l1_pdgid=tightElectrongs_pdgid[0]
	      ttc_l2_pdgid=fakeable_Electrons_pdgid[0]
	      ttc_l1_pt=tightElectrongs[0].Pt()
	      ttc_l1_eta=tightElectrongs[0].Eta()
	      ttc_l1_phi=tightElectrongs[0].Phi()
	      ttc_l1_mass=tightElectrongs[0].M()
	      ttc_l2_pt=fakeable_Electrons[0].Pt()
	      ttc_l2_eta=fakeable_Electrons[0].Eta()
	      ttc_l2_phi=fakeable_Electrons[0].Phi()
	      ttc_l2_mass=fakeable_Electrons[0].M()
	    else:
              ttc_lep1_faketag=True
              ttc_lep2_faketag=False
  	      ttc_l2_id=tightElectrons_id[0]
      	      ttc_l1_id=fakeable_Electrons_id[0]
	      ttc_l2_pdgid=tightElectrongs_pdgid[0]
	      ttc_l1_pdgid=fakeable_Electrons_pdgid[0]
	      ttc_l2_pt=tightElectrongs[0].Pt()
	      ttc_l2_eta=tightElectrongs[0].Eta()
	      ttc_l2_phi=tightElectrongs[0].Phi()
	      ttc_l2_mass=tightElectrongs[0].M()
	      ttc_l1_pt=fakeable_Electrons[0].Pt()
	      ttc_l1_eta=fakeable_Electrons[0].Eta()
	      ttc_l1_phi=fakeable_Electrons[0].Phi()
	      ttc_l1_mass=fakeable_Electrons[0].M()
      	    ttc_mll=(fakeable_Electrons[0]+tightElectrongs[0]).M()
      	    ttc_drll=fakeable_Electrons[0].DeltaR(tightElectrongs[0])
      	    ttc_dphill=fakeable_Electrons[0].DeltaPhi(tightElectrongs[0])

    if ttc_region>0:
      l1_v4_temp=TLorentzVector()
      l2_v4_temp=TLorentzVector()
      l1_v4_temp.SetPtEtaPhiM(ttc_l1_pt,ttc_l1_eta,ttc_l1_phi,ttc_l1_mass)
      l2_v4_temp.SetPtEtaPhiM(ttc_l2_pt,ttc_l2_eta,ttc_l2_phi,ttc_l2_mass)
      j1_v4_temp=TLorentzVector()
      j2_v4_temp=TLorentzVector()
      j3_v4_temp=TLorentzVector()
      j4_v4_temp=TLorentzVector()
      if n_tight_jet>3:
	j1_v4_temp.SetPtEtaPhiM(j1_pt,j1_eta,j1_phi,j1_mass)
	j2_v4_temp.SetPtEtaPhiM(j2_pt,j2_eta,j2_phi,j2_mass)
	j3_v4_temp.SetPtEtaPhiM(j3_pt,j3_eta,j3_phi,j3_mass)
	j4_v4_temp.SetPtEtaPhiM(j4_pt,j4_eta,j4_phi,j4_mass)
	ttc_dr_l1j1=l1_v4_temp.DeltaR(j1_v4_temp)
	ttc_dr_l1j2=l1_v4_temp.DeltaR(j2_v4_temp)
	ttc_dr_l1j3=l1_v4_temp.DeltaR(j3_v4_temp)
	ttc_dr_l1j4=l1_v4_temp.DeltaR(j4_v4_temp)
	ttc_dr_l2j1=l2_v4_temp.DeltaR(j1_v4_temp)
	ttc_dr_l2j2=l2_v4_temp.DeltaR(j2_v4_temp)
	ttc_dr_l2j3=l2_v4_temp.DeltaR(j3_v4_temp)
	ttc_dr_l2j4=l2_v4_temp.DeltaR(j4_v4_temp)
	ttc_mllj1=(l1_v4_temp+l2_v4_temp+j1_v4_temp).M()
	ttc_mllj2=(l1_v4_temp+l2_v4_temp+j2_v4_temp).M()
	ttc_mllj3=(l1_v4_temp+l2_v4_temp+j3_v4_temp).M()
	ttc_mllj4=(l1_v4_temp+l2_v4_temp+j4_v4_temp).M()
      if n_tight_jet==3:
	j1_v4_temp.SetPtEtaPhiM(j1_pt,j1_eta,j1_phi,j1_mass)
        j2_v4_temp.SetPtEtaPhiM(j2_pt,j2_eta,j2_phi,j2_mass)
        j3_v4_temp.SetPtEtaPhiM(j3_pt,j3_eta,j3_phi,j3_mass)
        ttc_dr_l1j1=l1_v4_temp.DeltaR(j1_v4_temp)
        ttc_dr_l1j2=l1_v4_temp.DeltaR(j2_v4_temp)
        ttc_dr_l1j3=l1_v4_temp.DeltaR(j3_v4_temp)
        ttc_dr_l2j1=l2_v4_temp.DeltaR(j1_v4_temp)
        ttc_dr_l2j2=l2_v4_temp.DeltaR(j2_v4_temp)
        ttc_dr_l2j3=l2_v4_temp.DeltaR(j3_v4_temp)
        ttc_mllj1=(l1_v4_temp+l2_v4_temp+j1_v4_temp).M()
        ttc_mllj2=(l1_v4_temp+l2_v4_temp+j2_v4_temp).M()
        ttc_mllj3=(l1_v4_temp+l2_v4_temp+j3_v4_temp).M()
      if n_tight_jet==2:
	j1_v4_temp.SetPtEtaPhiM(j1_pt,j1_eta,j1_phi,j1_mass)
        j2_v4_temp.SetPtEtaPhiM(j2_pt,j2_eta,j2_phi,j2_mass)
        ttc_dr_l1j1=l1_v4_temp.DeltaR(j1_v4_temp)
        ttc_dr_l1j2=l1_v4_temp.DeltaR(j2_v4_temp)
        ttc_dr_l2j1=l2_v4_temp.DeltaR(j1_v4_temp)
        ttc_dr_l2j2=l2_v4_temp.DeltaR(j2_v4_temp)
        ttc_mllj1=(l1_v4_temp+l2_v4_temp+j1_v4_temp).M()
        ttc_mllj2=(l1_v4_temp+l2_v4_temp+j2_v4_temp).M()
      if n_tight_jet==1:
	j1_v4_temp.SetPtEtaPhiM(j1_pt,j1_eta,j1_phi,j1_mass)
        ttc_dr_l1j1=l1_v4_temp.DeltaR(j1_v4_temp)
        ttc_dr_l2j1=l2_v4_temp.DeltaR(j1_v4_temp)
        ttc_mllj1=(l1_v4_temp+l2_v4_temp+j1_v4_temp).M()

    self.out.fillBranch("ttc_nl", ttc_nl)
    self.out.fillBranch("ttc_jets", ttc_jets)
    self.out.fillBranch("ttc_1P1F", ttc_1P1F)
    self.out.fillBranch("ttc_0P2F", ttc_0P2F)
    self.out.fillBranch("ttc_lep1_faketag", ttc_lep1_faketag)
    self.out.fillBranch("ttc_lep2_faketag", ttc_lep2_faketag)
    self.out.fillBranch("ttc_region", ttc_region)
    self.out.fillBranch("ttc_l1_id", ttc_l1_id)
    self.out.fillBranch("ttc_l2_id", ttc_l2_id)
    self.out.fillBranch("ttc_l1_pdgid", ttc_l1_pdgid)
    self.out.fillBranch("ttc_l2_pdgid", ttc_l2_pdgid)
    self.out.fillBranch("ttc_l1_pt", ttc_l1_pt)
    self.out.fillBranch("ttc_l1_eta", ttc_l1_eta)
    self.out.fillBranch("ttc_l1_phi", ttc_l1_phi)
    self.out.fillBranch("ttc_l1_mass", ttc_l1_mass)
    self.out.fillBranch("ttc_l2_pt", ttc_l2_pt) 
    self.out.fillBranch("ttc_l2_eta", ttc_l2_eta) 
    self.out.fillBranch("ttc_l2_phi", ttc_l2_phi) 
    self.out.fillBranch("ttc_l2_mass", ttc_l2_mass)
    self.out.fillBranch("ttc_mll",ttc_mll) 
    self.out.fillBranch("ttc_drll",ttc_drll) 
    self.out.fillBranch("ttc_dphill",ttc_dphill) 
    self.out.fillBranch("ttc_met", ttc_met)
    self.out.fillBranch("ttc_met_phi", ttc_met_phi)
    self.out.fillBranch("ttc_dr_l1j1", ttc_dr_l1j1)
    self.out.fillBranch("ttc_dr_l1j2", ttc_dr_l1j2)
    self.out.fillBranch("ttc_dr_l1j3", ttc_dr_l1j3)
    self.out.fillBranch("ttc_dr_l1j4", ttc_dr_l1j4)
    self.out.fillBranch("ttc_dr_l2j1", ttc_dr_l2j1)
    self.out.fillBranch("ttc_dr_l2j2", ttc_dr_l2j2)
    self.out.fillBranch("ttc_dr_l2j3", ttc_dr_l2j3)
    self.out.fillBranch("ttc_dr_l2j4", ttc_dr_l2j4)
    self.out.fillBranch("ttc_mllj1", ttc_mllj1)
    self.out.fillBranch("ttc_mllj2", ttc_mllj2)
    self.out.fillBranch("ttc_mllj3", ttc_mllj3)
    self.out.fillBranch("ttc_mllj4", ttc_mllj4)

    # WW     WWW     WW ZZZZZZZZZ   region: only 3 tight leptons, no b-jet, mll>4, |Z-91.1876|<15
    # WW     WWW     WW       ZZ            MET>30
    #  WW   WW WW   WW      ZZ
    #   WW WW   WW WW     ZZ
    #    WWW     WWW    ZZZZZZZZZ

    #WZ region lepton number selections
    WZ_nl=False
    #WZ region b-jet selection
    WZ_nb=False
    #WZ region lepton kinematics selctions
    WZ_leptons=False
    #WZ region MET selection
    WZ_MET=False
    #WZ region tag, 0: fail to pass the WZ selection, 1:3 muon, 2:2muon, 3:1muon, 4:0 muon
    WZ_region=0
    WZ_zl1_id=-1
    WZ_zl2_id=-1
    WZ_wl_id=-1
    WZ_zl1_pdgid=-99
    WZ_zl2_pdgid=-99
    WZ_wl_pdgid=-99
    WZ_zl1_pt=-99
    WZ_zl1_eta=-99
    WZ_zl1_phi=-99
    WZ_zl1_mass=-99
    WZ_zl2_pt=-99
    WZ_zl2_eta=-99
    WZ_zl2_phi=-99
    WZ_zl2_mass=-99
    WZ_l3_pt=-99
    WZ_l3_eta=-99
    WZ_l3_phi=-99
    WZ_l3_mass=-99
    WZ_Z_mass=-99
    WZ_Z_pt=-99
    WZ_Z_eta=-99
    WZ_Z_phi=-99
    WZ_met=-99
    # the first two leading leptons with pt >20, 3rd lepton pt >15, 4th lepton veto
    if len(tightLeptons)==3 and tightLeptons[1].Pt()>20 and len(looseLeptons)==0:
      WZ_nl=True
    # no bjet
    if WZ_nl and tightJets_b_DeepCSVmedium_id[0]==-1:
      WZ_nb=True

    # mll>4 regardless the flavor and charge sign
    if WZ_nb and (tightLeptons[0]+tightLeptons[1]).M()>4 and (tightLeptons[2]+tightLeptons[1]).M()>4 and (tightLeptons[0]+tightLeptons[2]).M()>4:
      WZ_leptons=True
    
    if WZ_leptons and ((self.is_mc and event.MET_T1Smear_pt>30) or (event.MET_T1_pt>30 and (not (self.is_mc)))):
      WZ_MET=True

    if WZ_MET:
      if self.is_mc:
        WZ_met=event.MET_T1Smear_pt
      else:
        WZ_met=event.MET_T1_pt
      # 3 muons case
      if len(tightElectrons)==0 and abs(tightMuons_pdgid[0]+tightMuons_pdgid[1]+tightMuons_pdgid[2])==13:
	#two combination 0+2 or 1+2
        if (tightMuons_pdgid[0]-tightMuons_pdgid[1])==0:
          if abs((tightMuons[0]+tightMuons[2]).M()-91.1876)<abs((tightMuons[1]+tightMuons[2]).M()-91.1876) and abs((tightMuons[0]+tightMuons[2]).M()-91.1876)<15:
	    WZ_region=1
            WZ_zl1_id=tightMuons_id[0]
            WZ_zl2_id=tightMuons_id[2]
            WZ_wl_id=tightMuons_id[1]
            WZ_zl1_pdgid=tightMuons_pdgid[0]
            WZ_zl2_pdgid=tightMuons_pdgid[2]
            WZ_wl_pdgid=tightMuons_pdgid[1]
	    WZ_zl1_pt=tightMuons[0].Pt()
	    WZ_zl1_eta=tightMuons[0].Eta()
	    WZ_zl1_phi=tightMuons[0].Phi()
	    WZ_zl1_mass=tightMuons[0].M()
	    WZ_zl2_pt=tightMuons[2].Pt()
	    WZ_zl2_eta=tightMuons[2].Eta()
	    WZ_zl2_phi=tightMuons[2].Phi()
	    WZ_zl2_mass=tightMuons[2].M()
	    WZ_l3_pt=tightMuons[1].Pt()
	    WZ_l3_eta=tightMuons[1].Eta()
	    WZ_l3_phi=tightMuons[1].Phi()
	    WZ_l3_mass=tightMuons[1].M()
	    WZ_Z_mass=(tightMuons[0]+tightMuons[2]).M()
	    WZ_Z_pt=(tightMuons[0]+tightMuons[2]).Pt()
	    WZ_Z_eta=(tightMuons[0]+tightMuons[2]).Eta()
	    WZ_Z_phi=(tightMuons[0]+tightMuons[2]).Phi()

	  if abs((tightMuons[0]+tightMuons[2]).M()-91.1876)>abs((tightMuons[1]+tightMuons[2]).M()-91.1876) and abs((tightMuons[1]+tightMuons[2]).M()-91.1876)<15:
	    WZ_region=1
            WZ_zl1_id=tightMuons_id[1]
            WZ_zl2_id=tightMuons_id[2]
            WZ_wl_id=tightMuons_id[0]
            WZ_zl1_pdgid=tightMuons_pdgid[1]
            WZ_zl2_pdgid=tightMuons_pdgid[2]
            WZ_wl_pdgid=tightMuons_pdgid[0]
	    WZ_zl1_pt=tightMuons[1].Pt()
	    WZ_zl1_eta=tightMuons[1].Eta()
	    WZ_zl1_phi=tightMuons[1].Phi()
	    WZ_zl1_mass=tightMuons[1].M()
	    WZ_zl2_pt=tightMuons[2].Pt()
	    WZ_zl2_eta=tightMuons[2].Eta()
	    WZ_zl2_phi=tightMuons[2].Phi()
	    WZ_zl2_mass=tightMuons[2].M()
	    WZ_l3_pt=tightMuons[0].Pt()
	    WZ_l3_eta=tightMuons[0].Eta()
	    WZ_l3_phi=tightMuons[0].Phi()
	    WZ_l3_mass=tightMuons[0].M()
	    WZ_Z_mass=(tightMuons[1]+tightMuons[2]).M()
	    WZ_Z_pt=(tightMuons[1]+tightMuons[2]).Pt()
	    WZ_Z_eta=(tightMuons[1]+tightMuons[2]).Eta()
	    WZ_Z_phi=(tightMuons[1]+tightMuons[2]).Phi()
	#two combination 0+1 or 1+2
	elif (tightMuons_pdgid[0]-tightMuons_pdgid[2])==0:
	  if abs((tightMuons[0]+tightMuons[1]).M()-91.1876)<abs((tightMuons[1]+tightMuons[2]).M()-91.1876) and abs((tightMuons[0]+tightMuons[1]).M()-91.1876)<15:
	    WZ_region=1
            WZ_zl1_id=tightMuons_id[0]
            WZ_zl2_id=tightMuons_id[1]
            WZ_wl_id=tightMuons_id[2]
            WZ_zl1_pdgid=tightMuons_pdgid[0]
            WZ_zl2_pdgid=tightMuons_pdgid[1]
            WZ_wl_pdgid=tightMuons_pdgid[2]
	    WZ_zl1_pt=tightMuons[0].Pt()
	    WZ_zl1_eta=tightMuons[0].Eta()
	    WZ_zl1_phi=tightMuons[0].Phi()
	    WZ_zl1_mass=tightMuons[0].M()
	    WZ_zl2_pt=tightMuons[1].Pt()
	    WZ_zl2_eta=tightMuons[1].Eta()
	    WZ_zl2_phi=tightMuons[1].Phi()
	    WZ_zl2_mass=tightMuons[1].M()
	    WZ_l3_pt=tightMuons[2].Pt()
	    WZ_l3_eta=tightMuons[2].Eta()
	    WZ_l3_phi=tightMuons[2].Phi()
	    WZ_l3_mass=tightMuons[2].M()
	    WZ_Z_mass=(tightMuons[0]+tightMuons[1]).M()
	    WZ_Z_pt=(tightMuons[0]+tightMuons[1]).Pt()
	    WZ_Z_eta=(tightMuons[0]+tightMuons[1]).Eta()
	    WZ_Z_phi=(tightMuons[0]+tightMuons[1]).Phi()

	  if abs((tightMuons[0]+tightMuons[1]).M()-91.1876)>abs((tightMuons[1]+tightMuons[2]).M()-91.1876) and abs((tightMuons[1]+tightMuons[2]).M()-91.1876)<15:
	    WZ_region=1
            WZ_zl1_id=tightMuons_id[1]
            WZ_zl2_id=tightMuons_id[2]
            WZ_wl_id=tightMuons_id[0]
            WZ_zl1_pdgid=tightMuons_pdgid[1]
            WZ_zl2_pdgid=tightMuons_pdgid[2]
            WZ_wl_pdgid=tightMuons_pdgid[0]
	    WZ_zl1_pt=tightMuons[1].Pt()
	    WZ_zl1_eta=tightMuons[1].Eta()
	    WZ_zl1_phi=tightMuons[1].Phi()
	    WZ_zl1_mass=tightMuons[1].M()
	    WZ_zl2_pt=tightMuons[2].Pt()
	    WZ_zl2_eta=tightMuons[2].Eta()
	    WZ_zl2_phi=tightMuons[2].Phi()
	    WZ_zl2_mass=tightMuons[2].M()
	    WZ_l3_pt=tightMuons[0].Pt()
	    WZ_l3_eta=tightMuons[0].Eta()
	    WZ_l3_phi=tightMuons[0].Phi()
	    WZ_l3_mass=tightMuons[0].M()
	    WZ_Z_mass=(tightMuons[1]+tightMuons[2]).M()
	    WZ_Z_pt=(tightMuons[1]+tightMuons[2]).Pt()
	    WZ_Z_eta=(tightMuons[1]+tightMuons[2]).Eta()
	    WZ_Z_phi=(tightMuons[1]+tightMuons[2]).Phi()
	#two combination 0+1 or 0+2
	else:
	  if abs((tightMuons[0]+tightMuons[1]).M()-91.1876)<abs((tightMuons[0]+tightMuons[2]).M()-91.1876) and abs((tightMuons[0]+tightMuons[1]).M()-91.1876)<15:
	    WZ_region=1
            WZ_zl1_id=tightMuons_id[0]
            WZ_zl2_id=tightMuons_id[1]
            WZ_wl_id=tightMuons_id[2]
            WZ_zl1_pdgid=tightMuons_pdgid[0]
            WZ_zl2_pdgid=tightMuons_pdgid[1]
            WZ_wl_pdgid=tightMuons_pdgid[2]
	    WZ_zl1_pt=tightMuons[0].Pt()
	    WZ_zl1_eta=tightMuons[0].Eta()
	    WZ_zl1_phi=tightMuons[0].Phi()
	    WZ_zl1_mass=tightMuons[0].M()
	    WZ_zl2_pt=tightMuons[1].Pt()
	    WZ_zl2_eta=tightMuons[1].Eta()
	    WZ_zl2_phi=tightMuons[1].Phi()
	    WZ_zl2_mass=tightMuons[1].M()
	    WZ_l3_pt=tightMuons[2].Pt()
	    WZ_l3_eta=tightMuons[2].Eta()
	    WZ_l3_phi=tightMuons[2].Phi()
	    WZ_l3_mass=tightMuons[2].M()
	    WZ_Z_mass=(tightMuons[0]+tightMuons[1]).M()
	    WZ_Z_pt=(tightMuons[0]+tightMuons[1]).Pt()
	    WZ_Z_eta=(tightMuons[0]+tightMuons[1]).Eta()
	    WZ_Z_phi=(tightMuons[0]+tightMuons[1]).Phi()
	  if abs((tightMuons[0]+tightMuons[1]).M()-91.1876)>abs((tightMuons[0]+tightMuons[2]).M()-91.1876) and abs((tightMuons[0]+tightMuons[2]).M()-91.1876)<15:
	    WZ_region=1
            WZ_zl1_id=tightMuons_id[0]
            WZ_zl2_id=tightMuons_id[2]
            WZ_wl_id=tightMuons_id[1]
            WZ_zl1_pdgid=tightMuons_pdgid[0]
            WZ_zl2_pdgid=tightMuons_pdgid[2]
            WZ_wl_pdgid=tightMuons_pdgid[1]
	    WZ_zl1_pt=tightMuons[0].Pt()
	    WZ_zl1_eta=tightMuons[0].Eta()
	    WZ_zl1_phi=tightMuons[0].Phi()
	    WZ_zl1_mass=tightMuons[0].M()
	    WZ_zl2_pt=tightMuons[2].Pt()
	    WZ_zl2_eta=tightMuons[2].Eta()
	    WZ_zl2_phi=tightMuons[2].Phi()
	    WZ_zl2_mass=tightMuons[2].M()
	    WZ_l3_pt=tightMuons[1].Pt()
	    WZ_l3_eta=tightMuons[1].Eta()
	    WZ_l3_phi=tightMuons[1].Phi()
	    WZ_l3_mass=tightMuons[1].M()
	    WZ_Z_mass=(tightMuons[0]+tightMuons[2]).M()
	    WZ_Z_pt=(tightMuons[0]+tightMuons[2]).Pt()
	    WZ_Z_eta=(tightMuons[0]+tightMuons[2]).Eta()
	    WZ_Z_phi=(tightMuons[0]+tightMuons[2]).Phi()

      # 2 muons case
      if len(tightElectrons)==1 and (tightMuons_pdgid[0]-tightMuons_pdgid[1])==0:
	if abs((tightMuons[0]+tightMuons[1]).M()-91.1876)<15:
	  WZ_region=2
	  WZ_zl1_id=tightMuons_id[0]
	  WZ_zl2_id=tightMuons_id[1]
	  WZ_wl_id=tightElectrons_id[0]
	  WZ_zl1_pdgid=tightMuons_pdgid[0]
	  WZ_zl2_pdgid=tightMuons_pdgid[1]
	  WZ_wl_pdgid=tightElectrons_pdgid[0]
	  WZ_zl1_pt=tightMuons[0].Pt()
	  WZ_zl1_eta=tightMuons[0].Eta()
	  WZ_zl1_phi=tightMuons[0].Phi()
	  WZ_zl1_mass=tightMuons[0].M()
	  WZ_zl2_pt=tightMuons[1].Pt()
	  WZ_zl2_eta=tightMuons[1].Eta()
	  WZ_zl2_phi=tightMuons[1].Phi()
	  WZ_zl2_mass=tightMuons[1].M()
	  WZ_l3_pt=tightElectrons[0].Pt()
	  WZ_l3_eta=tightElectrons[0].Eta()
	  WZ_l3_phi=tightElectrons[0].Phi()
	  WZ_l3_mass=tightElectrons[0].M()
	  WZ_Z_mass=(tightMuons[0]+tightMuons[1]).M()
	  WZ_Z_pt=(tightMuons[0]+tightMuons[1]).Pt()
	  WZ_Z_eta=(tightMuons[0]+tightMuons[1]).Eta()
	  WZ_Z_phi=(tightMuons[0]+tightMuons[1]).Phi()

      # 1 muon case
      if len(tightElectrons)==2 and (tightElectrons_pdgid[0]-tightElectrons_pdgid[1])==0:
	if abs((tightElectrons[0]+tightElectrons[1]).M()-91.1876)<15:
	  WZ_region=3
	  WZ_zl1_id=tightElectrons_id[0]
	  WZ_zl2_id=tightElectrons_id[1]
	  WZ_wl_id=tightMuons_id[0]
	  WZ_zl1_pdgid=tightElectrons_pdgid[0]
	  WZ_zl2_pdgid=tightElectrons_pdgid[1]
	  WZ_wl_pdgid=tightMuons_pdgid[0]
	  WZ_zl1_pt=tightElectrons[0].Pt()
	  WZ_zl1_eta=tightElectrons[0].Eta()
	  WZ_zl1_phi=tightElectrons[0].Phi()
	  WZ_zl1_mass=tightElectrons[0].M()
	  WZ_zl2_pt=tightElectrons[1].Pt()
	  WZ_zl2_eta=tightElectrons[1].Eta()
	  WZ_zl2_phi=tightElectrons[1].Phi()
	  WZ_zl2_mass=tightElectrons[1].M()
	  WZ_l3_pt=tightMuons[0].Pt()
	  WZ_l3_eta=tightMuons[0].Eta()
	  WZ_l3_phi=tightMuons[0].Phi()
	  WZ_l3_mass=tightMuons[0].M()
	  WZ_Z_mass=(tightElectrons[0]+tightElectrons[1]).M()
	  WZ_Z_pt=(tightElectrons[0]+tightElectrons[1]).Pt()
	  WZ_Z_eta=(tightElectrons[0]+tightElectrons[1]).Eta()
	  WZ_Z_phi=(tightElectrons[0]+tightElectrons[1]).Phi()

      # 0 muon case
      if len(tightElectrons)==3 and abs(tightElectrons_pdgid[0]+tightElectrons_pdgid[1]+tightElectrons_pdgid[2])==11:
	#two combination 0+2 or 1+2
        if (tightElectrons_pdgid[0]-tightElectrons_pdgid[1])==0:
          if abs((tightElectrons[0]+tightElectrons[2]).M()-91.1876)<abs((tightElectrons[1]+tightElectrons[2]).M()-91.1876) and abs((tightElectrons[0]+tightElectrons[2]).M()-91.1876)<15:
	    WZ_region=4
            WZ_zl1_id=tightElectrons_id[0]
            WZ_zl2_id=tightElectrons_id[2]
            WZ_wl_id=tightElectrons_id[1]
            WZ_zl1_pdgid=tightElectrons_pdgid[0]
            WZ_zl2_pdgid=tightElectrons_pdgid[2]
            WZ_wl_pdgid=tightElectrons_pdgid[1]
	    WZ_zl1_pt=tightElectrons[0].Pt()
	    WZ_zl1_eta=tightElectrons[0].Eta()
	    WZ_zl1_phi=tightElectrons[0].Phi()
	    WZ_zl1_mass=tightElectrons[0].M()
	    WZ_zl2_pt=tightElectrons[2].Pt()
	    WZ_zl2_eta=tightElectrons[2].Eta()
	    WZ_zl2_phi=tightElectrons[2].Phi()
	    WZ_zl2_mass=tightElectrons[2].M()
	    WZ_l3_pt=tightElectrons[1].Pt()
	    WZ_l3_eta=tightElectrons[1].Eta()
	    WZ_l3_phi=tightElectrons[1].Phi()
	    WZ_l3_mass=tightElectrons[1].M()
	    WZ_Z_mass=(tightElectrons[0]+tightElectrons[2]).M()
	    WZ_Z_pt=(tightElectrons[0]+tightElectrons[2]).Pt()
	    WZ_Z_eta=(tightElectrons[0]+tightElectrons[2]).Eta()
	    WZ_Z_phi=(tightElectrons[0]+tightElectrons[2]).Phi()
	  if abs((tightElectrons[0]+tightElectrons[2]).M()-91.1876)>abs((tightElectrons[1]+tightElectrons[2]).M()-91.1876) and abs((tightElectrons[1]+tightElectrons[2]).M()-91.1876)<15:
	    WZ_region=4
            WZ_zl1_id=tightElectrons_id[1]
            WZ_zl2_id=tightElectrons_id[2]
            WZ_wl_id=tightElectrons_id[0]
            WZ_zl1_pdgid=tightElectrons_pdgid[1]
            WZ_zl2_pdgid=tightElectrons_pdgid[2]
            WZ_wl_pdgid=tightElectrons_pdgid[0]
	    WZ_zl1_pt=tightElectrons[1].Pt()
	    WZ_zl1_eta=tightElectrons[1].Eta()
	    WZ_zl1_phi=tightElectrons[1].Phi()
	    WZ_zl1_mass=tightElectrons[1].M()
	    WZ_zl2_pt=tightElectrons[2].Pt()
	    WZ_zl2_eta=tightElectrons[2].Eta()
	    WZ_zl2_phi=tightElectrons[2].Phi()
	    WZ_zl2_mass=tightElectrons[2].M()
	    WZ_l3_pt=tightElectrons[0].Pt()
	    WZ_l3_eta=tightElectrons[0].Eta()
	    WZ_l3_phi=tightElectrons[0].Phi()
	    WZ_l3_mass=tightElectrons[0].M()
	    WZ_Z_mass=(tightElectrons[1]+tightElectrons[2]).M()
	    WZ_Z_pt=(tightElectrons[1]+tightElectrons[2]).Pt()
	    WZ_Z_eta=(tightElectrons[1]+tightElectrons[2]).Eta()
	    WZ_Z_phi=(tightElectrons[1]+tightElectrons[2]).Phi()
	#two combination 0+1 or 1+2
	elif (tightElectrons_pdgid[0]-tightElectrons_pdgid[2])==0:
	  if abs((tightElectrons[0]+tightElectrons[1]).M()-91.1876)<abs((tightElectrons[1]+tightElectrons[2]).M()-91.1876) and abs((tightElectrons[0]+tightElectrons[1]).M()-91.1876)<15:
	    WZ_region=4
            WZ_zl1_id=tightElectrons_id[0]
            WZ_zl2_id=tightElectrons_id[1]
            WZ_wl_id=tightElectrons_id[2]
            WZ_zl1_pdgid=tightElectrons_pdgid[0]
            WZ_zl2_pdgid=tightElectrons_pdgid[1]
            WZ_wl_pdgid=tightElectrons_pdgid[2]
	    WZ_zl1_pt=tightElectrons[0].Pt()
	    WZ_zl1_eta=tightElectrons[0].Eta()
	    WZ_zl1_phi=tightElectrons[0].Phi()
	    WZ_zl1_mass=tightElectrons[0].M()
	    WZ_zl2_pt=tightElectrons[1].Pt()
	    WZ_zl2_eta=tightElectrons[1].Eta()
	    WZ_zl2_phi=tightElectrons[1].Phi()
	    WZ_zl2_mass=tightElectrons[1].M()
	    WZ_l3_pt=tightElectrons[2].Pt()
	    WZ_l3_eta=tightElectrons[2].Eta()
	    WZ_l3_phi=tightElectrons[2].Phi()
	    WZ_l3_mass=tightElectrons[2].M()
	    WZ_Z_mass=(tightElectrons[0]+tightElectrons[1]).M()
	    WZ_Z_pt=(tightElectrons[0]+tightElectrons[1]).Pt()
	    WZ_Z_eta=(tightElectrons[0]+tightElectrons[1]).Eta()
	    WZ_Z_phi=(tightElectrons[0]+tightElectrons[1]).Phi()
	  if abs((tightElectrons[0]+tightElectrons[1]).M()-91.1876)>abs((tightElectrons[1]+tightElectrons[2]).M()-91.1876) and abs((tightElectrons[1]+tightElectrons[2]).M()-91.1876)<15:
	    WZ_region=4
            WZ_zl1_id=tightElectrons_id[1]
            WZ_zl2_id=tightElectrons_id[2]
            WZ_wl_id=tightElectrons_id[0]
            WZ_zl1_pdgid=tightElectrons_pdgid[1]
            WZ_zl2_pdgid=tightElectrons_pdgid[2]
            WZ_wl_pdgid=tightElectrons_pdgid[0]
	    WZ_zl1_pt=tightElectrons[1].Pt()
	    WZ_zl1_eta=tightElectrons[1].Eta()
	    WZ_zl1_phi=tightElectrons[1].Phi()
	    WZ_zl1_mass=tightElectrons[1].M()
	    WZ_zl2_pt=tightElectrons[2].Pt()
	    WZ_zl2_eta=tightElectrons[2].Eta()
	    WZ_zl2_phi=tightElectrons[2].Phi()
	    WZ_zl2_mass=tightElectrons[2].M()
	    WZ_l3_pt=tightElectrons[0].Pt()
	    WZ_l3_eta=tightElectrons[0].Eta()
	    WZ_l3_phi=tightElectrons[0].Phi()
	    WZ_l3_mass=tightElectrons[0].M()
	    WZ_Z_mass=(tightElectrons[1]+tightElectrons[2]).M()
	    WZ_Z_pt=(tightElectrons[1]+tightElectrons[2]).Pt()
	    WZ_Z_eta=(tightElectrons[1]+tightElectrons[2]).Eta()
	    WZ_Z_phi=(tightElectrons[1]+tightElectrons[2]).Phi()
	#two combination 0+1 or 0+2
	else:
	  if abs((tightElectrons[0]+tightElectrons[1]).M()-91.1876)<abs((tightElectrons[0]+tightElectrons[2]).M()-91.1876) and abs((tightElectrons[0]+tightElectrons[1]).M()-91.1876)<15:
	    WZ_region=4
            WZ_zl1_id=tightElectrons_id[0]
            WZ_zl2_id=tightElectrons_id[1]
            WZ_wl_id=tightElectrons_id[2]
            WZ_zl1_pdgid=tightElectrons_pdgid[0]
            WZ_zl2_pdgid=tightElectrons_pdgid[1]
            WZ_wl_pdgid=tightElectrons_pdgid[2]
	    WZ_zl1_pt=tightElectrons[0].Pt()
	    WZ_zl1_eta=tightElectrons[0].Eta()
	    WZ_zl1_phi=tightElectrons[0].Phi()
	    WZ_zl1_mass=tightElectrons[0].M()
	    WZ_zl2_pt=tightElectrons[1].Pt()
	    WZ_zl2_eta=tightElectrons[1].Eta()
	    WZ_zl2_phi=tightElectrons[1].Phi()
	    WZ_zl2_mass=tightElectrons[1].M()
	    WZ_l3_pt=tightElectrons[2].Pt()
	    WZ_l3_eta=tightElectrons[2].Eta()
	    WZ_l3_phi=tightElectrons[2].Phi()
	    WZ_l3_mass=tightElectrons[2].M()
	    WZ_Z_mass=(tightElectrons[0]+tightElectrons[1]).M()
	    WZ_Z_pt=(tightElectrons[0]+tightElectrons[1]).Pt()
	    WZ_Z_eta=(tightElectrons[0]+tightElectrons[1]).Eta()
	    WZ_Z_phi=(tightElectrons[0]+tightElectrons[1]).Phi()
	  if abs((tightElectrons[0]+tightElectrons[1]).M()-91.1876)>abs((tightElectrons[0]+tightElectrons[2]).M()-91.1876) and abs((tightElectrons[0]+tightElectrons[2]).M()-91.1876)<15:
	    WZ_region=4
            WZ_zl1_id=tightElectrons_id[0]
            WZ_zl2_id=tightElectrons_id[2]
            WZ_wl_id=tightElectrons_id[1]
            WZ_zl1_pdgid=tightElectrons_pdgid[0]
            WZ_zl2_pdgid=tightElectrons_pdgid[2]
            WZ_wl_pdgid=tightElectrons_pdgid[1]
	    WZ_zl1_pt=tightElectrons[0].Pt()
	    WZ_zl1_eta=tightElectrons[0].Eta()
	    WZ_zl1_phi=tightElectrons[0].Phi()
	    WZ_zl1_mass=tightElectrons[0].M()
	    WZ_zl2_pt=tightElectrons[2].Pt()
	    WZ_zl2_eta=tightElectrons[2].Eta()
	    WZ_zl2_phi=tightElectrons[2].Phi()
	    WZ_zl2_mass=tightElectrons[2].M()
	    WZ_l3_pt=tightElectrons[1].Pt()
	    WZ_l3_eta=tightElectrons[1].Eta()
	    WZ_l3_phi=tightElectrons[1].Phi()
	    WZ_l3_mass=tightElectrons[1].M()
	    WZ_Z_mass=(tightElectrons[0]+tightElectrons[2]).M()
	    WZ_Z_pt=(tightElectrons[0]+tightElectrons[2]).Pt()
	    WZ_Z_eta=(tightElectrons[0]+tightElectrons[2]).Eta()
	    WZ_Z_phi=(tightElectrons[0]+tightElectrons[2]).Phi()
	
    self.out.fillBranch("WZ_region", WZ_region)
    self.out.fillBranch("WZ_zl1_id", WZ_zl1_id)
    self.out.fillBranch("WZ_zl2_id", WZ_zl2_id)
    self.out.fillBranch("WZ_wl_id", WZ_wl_id)
    self.out.fillBranch("WZ_zl1_pdgid", WZ_zl1_pdgid)
    self.out.fillBranch("WZ_zl2_pdgid", WZ_zl2_pdgid)
    self.out.fillBranch("WZ_wl_pdgid", WZ_wl_pdgid)
    self.out.fillBranch("WZ_zl1_pt", WZ_zl1_pt)
    self.out.fillBranch("WZ_zl1_eta", WZ_zl1_eta)
    self.out.fillBranch("WZ_zl1_phi", WZ_zl1_phi)
    self.out.fillBranch("WZ_zl1_mass", WZ_zl1_mass)
    self.out.fillBranch("WZ_zl2_pt", WZ_zl2_pt)
    self.out.fillBranch("WZ_zl2_eta", WZ_zl2_eta)
    self.out.fillBranch("WZ_zl2_phi", WZ_zl2_phi)
    self.out.fillBranch("WZ_zl2_mass", WZ_zl2_mass)
    self.out.fillBranch("WZ_l3_pt", WZ_l3_pt)
    self.out.fillBranch("WZ_l3_eta", WZ_l3_eta)
    self.out.fillBranch("WZ_l3_phi", WZ_l3_phi)
    self.out.fillBranch("WZ_l3_mass", WZ_l3_mass)
    self.out.fillBranch("WZ_Z_mass", WZ_Z_mass)
    self.out.fillBranch("WZ_Z_pt", WZ_Z_pt)
    self.out.fillBranch("WZ_Z_eta", WZ_Z_eta)
    self.out.fillBranch("WZ_Z_phi", WZ_Z_phi)
    self.out.fillBranch("WZ_met", WZ_met)
    

    #  DDDD   YY      YY  (opposite sign) region: two opposite sign lepton, with |mll-91.1876|<15
    #  D   D    YY  YY   
    #  D    D     YY
    #  D   D      YY
    #  DDDD       YY

    #DY region lepton number selections
    DY_nl=False
    DY_1P1F=False
    DY_0P2F=False
    DY_lep1_faketag=False
    DY_lep2_faketag=False
    #DY region tag, 0: fail to pass the DY selection, 1:2 muon, 2:1 muon, 3:0 muon
    DY_region=0
    DY_l1_id=-1
    DY_l2_id=-1
    DY_l1_pdgid=-99
    DY_l2_pdgid=-99
    DY_l1_pt=-99
    DY_l1_pt_raw=-99
    DY_l1_eta=-99
    DY_l1_phi=-99
    DY_l1_mass=-99
    DY_l2_pt=-99
    DY_l2_pt_raw=-99
    DY_l2_eta=-99
    DY_l2_phi=-99
    DY_l2_mass=-99
    DY_z_mass=-99
    DY_z_mass_raw=-99
    DY_z_pt=-99
    DY_z_pt_raw=-99
    DY_z_eta=-99
    DY_z_eta_raw=-99
    DY_z_phi=-99
    DY_z_phi_raw=-99
    DY_drll=-99
    DY_met=-99
    DY_met_phi=-99

    # the two leptons with pt >20, 3th lepton veto
    if ((len(fakeableLeptons)==2 and fakeableLeptons[1].Pt()>20) or (len(tightLeptons)==1 and len(fakeableLeptons)==1 and tightLeptons[0].Pt()>20 and fakeableLeptons[0].Pt()>20)) and len(looseLeptons)==0:
      DY_nl=True

    if DY_nl:

      if self.is_mc:
        DY_met=event.MET_T1Smear_pt
        DY_met_phi=event.MET_T1Smear_phi
      else:
        DY_met=event.MET_T1_pt
        DY_met_phi=event.MET_T1_phi

      if len(fakeableLeptons)==2:
        DY_0P2F=True
        DY_lep1_faketag=True
        DY_lep2_faketag=True

        # 2 muons case
	if len(fakeable_Electrons)==0 and abs(fakeable_Muons_pdgid[0]+fakeable_Muons_pdgid[1])==0:
	  DY_region=1
	  DY_l1_id=fakeable_Muons_id[0]
	  DY_l2_id=fakeable_Muons_id[1]
	  DY_l1_pdgid=fakeable_Muons_pdgid[0]
	  DY_l2_pdgid=fakeable_Muons_pdgid[1]
	  DY_l1_pt=fakeable_Muons[0].Pt()
	  DY_l1_eta=fakeable_Muons[0].Eta()
	  DY_l1_phi=fakeable_Muons[0].Phi()
	  DY_l1_mass=fakeable_Muons[0].M()
	  DY_l2_pt=fakeable_Muons[1].Pt()
	  DY_l2_eta=fakeable_Muons[1].Eta()
	  DY_l2_phi=fakeable_Muons[1].Phi()
	  DY_l2_mass=fakeable_Muons[1].M()
	  DY_z_mass=(fakeableLeptons[0]+fakeableLeptons[1]).M()
	  DY_z_pt=(fakeableLeptons[0]+fakeableLeptons[1]).Pt()
	  DY_z_eta=(fakeableLeptons[0]+fakeableLeptons[1]).Eta()
	  DY_z_phi=(fakeableLeptons[0]+fakeableLeptons[1]).Phi()
	  DY_drll=fakeableLeptons[0].DeltaR(fakeableLeptons[1])
        # 2 eles case
        if len(fakeable_Electrons)==2 and abs(fakeable_Electrons_pdgid[0]+fakeable_Electrons_pdgid[1])==0:
	  DY_region=3
          DY_l1_id=fakeable_Electrons_id[0]
          DY_l2_id=fakeable_Electrons_id[1]
          DY_l1_pdgid=fakeable_Electrons_pdgid[0]
          DY_l2_pdgid=fakeable_Electrons_pdgid[1]
	  DY_l1_pt=fakeable_Electrons[0].Pt()
	  DY_l1_eta=fakeable_Electrons[0].Eta()
	  DY_l1_phi=fakeable_Electrons[0].Phi()
	  DY_l1_mass=fakeable_Electrons[0].M()
	  DY_l2_pt=fakeable_Electrons[1].Pt()
	  DY_l2_eta=fakeable_Electrons[1].Eta()
	  DY_l2_phi=fakeable_Electrons[1].Phi()
	  DY_l2_mass=fakeable_Electrons[1].M()
	  DY_z_mass=(fakeableLeptons[0]+fakeableLeptons[1]).M()
	  DY_z_pt=(fakeableLeptons[0]+fakeableLeptons[1]).Pt()
	  DY_z_eta=(fakeableLeptons[0]+fakeableLeptons[1]).Eta()
	  DY_z_phi=(fakeableLeptons[0]+fakeableLeptons[1]).Phi()
	  DY_drll=fakeableLeptons[0].DeltaR(fakeableLeptons[1])
        # 1 ele case
        if len(fakeable_Electrons)==1 and (sign(fakeable_Muons_pdgid[0])+sign(fakeable_Electrons_pdgid[0]))==0:
	  DY_region=2
          DY_l1_id=fakeable_Muons_id[0]
          DY_l2_id=fakeable_Electrons_id[0]
          DY_l1_pdgid=fakeable_Muons_pdgid[0]
          DY_l2_pdgid=fakeable_Electrons_pdgid[0]
	  DY_l1_pt=fakeable_Muons[0].Pt()
	  DY_l1_eta=fakeable_Muons[0].Eta()
	  DY_l1_phi=fakeable_Muons[0].Phi()
	  DY_l1_mass=fakeable_Muons[0].M()
	  DY_l2_pt=fakeable_Electrons[0].Pt()
	  DY_l2_eta=fakeable_Electrons[0].Eta()
	  DY_l2_phi=fakeable_Electrons[0].Phi()
	  DY_l2_mass=fakeable_Electrons[0].M()
	  DY_z_mass=(fakeableLeptons[0]+fakeableLeptons[1]).M()
	  DY_z_pt=(fakeableLeptons[0]+fakeableLeptons[1]).Pt()
	  DY_z_eta=(fakeableLeptons[0]+fakeableLeptons[1]).Eta()
	  DY_z_phi=(fakeableLeptons[0]+fakeableLeptons[1]).Phi()
	  DY_drll=fakeableLeptons[0].DeltaR(fakeableLeptons[1])

      if len(fakeableLeptons)==1:
        DY_1P1F=True

        # one real muons
        if len(tightElectrons)==0:
          # one fake muon
          if len(fakeable_Electrons)==0 and abs(tightMuons_pdgid[0]+fakeable_Muons_pdgid[0])==0:
            DY_region=1
            if tightMuons[0].Pt()>fakeable_Muons[0].Pt():
              DY_lep1_faketag=False
              DY_lep2_faketag=True
	      DY_l1_id=tightMuons_id[0]
	      DY_l2_id=fakeable_Muons_id[0]
	      DY_l1_pdgid=tightMuons_pdgid[0]
	      DY_l2_pdgid=fakeable_Muons_pdgid[0]
	      DY_l1_pt=tightMuons[0].Pt()
	      DY_l1_eta=tightMuons[0].Eta()
	      DY_l1_phi=tightMuons[0].Phi()
	      DY_l1_mass=tightMuons[0].M()
	      DY_l2_pt=fakeable_Muons[0].Pt()
	      DY_l2_eta=fakeable_Muons[0].Eta()
	      DY_l2_phi=fakeable_Muons[0].Phi()
	      DY_l2_mass=fakeable_Muons[0].M()
            else:
              DY_lep1_faketag=True
              DY_lep2_faketag=False
	      DY_l2_id=tightMuons_id[0]
	      DY_l1_id=fakeable_Muons_id[0]
	      DY_l2_pdgid=tightMuons_pdgid[0]
	      DY_l1_pdgid=fakeable_Muons_pdgid[0]
	      DY_l2_pt=tightMuons[0].Pt()
	      DY_l2_eta=tightMuons[0].Eta()
	      DY_l2_phi=tightMuons[0].Phi()
	      DY_l2_mass=tightMuons[0].M()
	      DY_l1_pt=fakeable_Muons[0].Pt()
	      DY_l1_eta=fakeable_Muons[0].Eta()
	      DY_l1_phi=fakeable_Muons[0].Phi()
	      DY_l1_mass=fakeable_Muons[0].M()
	    DY_z_mass=(fakeable_Muons[0]+tightMuons[0]).M()
	    DY_z_pt=(fakeable_Muons[0]+tightMuons[0]).Pt()
	    DY_z_eta=(fakeable_Muons[0]+tightMuons[0]).Eta()
	    DY_z_phi=(fakeable_Muons[0]+tightMuons[0]).Phi()
	    DY_drll=fakeable_Muons[0].DeltaR(tightMuons[0])

          # one fake electron
          if len(fakeable_Electrons)==1 and (sign(tightMuons_pdgid[0])+sign(fakeable_Electrons_pdgid[0]))==2:
            DY_region=2
            DY_lep1_faketag=False
            DY_lep2_faketag=True
	    DY_l1_id=tightMuons_id[0]
	    DY_l2_id=fakeable_Electrons_id[0]
	    DY_l1_pdgid=tightMuons_pdgid[0]
	    DY_l2_pdgid=fakeable_Electrons_pdgid[0]
	    DY_l1_pt=tightMuons[0].Pt()
	    DY_l1_eta=tightMuons[0].Eta()
	    DY_l1_phi=tightMuons[0].Phi()
	    DY_l1_mass=tightMuons[0].M()
	    DY_l2_pt=fakeable_Electrons[0].Pt()
	    DY_l2_eta=fakeable_Electrons[0].Eta()
	    DY_l2_phi=fakeable_Electrons[0].Phi()
	    DY_l2_mass=fakeable_Electrons[0].M()
	    DY_z_mass=(fakeable_Electrons[0]+tightMuons[0]).M()
	    DY_z_pt=(fakeable_Electrons[0]+tightMuons[0]).Pt()
	    DY_z_eta=(fakeable_Electrons[0]+tightMuons[0]).Eta()
	    DY_z_phi=(fakeable_Electrons[0]+tightMuons[0]).Phi()
	    DY_drll=fakeable_Electrons[0].DeltaR(tightMuons[0])

        # one real electrons
        if len(tightElectrons)==1:
          # one fake muon
          if len(fakeable_Electrons)==0 and abs(tightElectrons_pdgid[0]+fakeable_Muons_pdgid[0])==0:
            DY_region=2
            DY_lep1_faketag=False
            DY_lep2_faketag=True
	    DY_l1_id=tightElectrons_id[0]
	    DY_l2_id=fakeable_Muons_id[0]
	    DY_l1_pdgid=tightElectrons_pdgid[0]
	    DY_l2_pdgid=fakeable_Muons_pdgid[0]
	    DY_l1_pt=tightElectrons[0].Pt()
	    DY_l1_eta=tightElectrons[0].Eta()
	    DY_l1_phi=tightElectrons[0].Phi()
	    DY_l1_mass=tightElectrons[0].M()
	    DY_l2_pt=fakeable_Muons[0].Pt()
	    DY_l2_eta=fakeable_Muons[0].Eta()
	    DY_l2_phi=fakeable_Muons[0].Phi()
	    DY_l2_mass=fakeable_Muons[0].M()
	    DY_z_mass=(fakeable_Muons[0]+tightElectrons[0]).M()
	    DY_z_pt=(fakeable_Muons[0]+tightElectrons[0]).Pt()
	    DY_z_eta=(fakeable_Muons[0]+tightElectrons[0]).Eta()
	    DY_z_phi=(fakeable_Muons[0]+tightElectrons[0]).Phi()
	    DY_drll=fakeable_Muons[0].DeltaR(tightElectrons[0])

          # one fake electron
          if len(fakeable_Electrons)==1 and abs(tightElectrons_pdgid[0]+fakeable_Electrons_pdgid[0])==0:
            DY_region=3
	    if tightElectrons[0].Pt()>fakeable_Electrons[0].Pt():
              DY_lep1_faketag=False
              DY_lep2_faketag=True
  	      DY_l1_id=tightElectrons_id[0]
  	      DY_l2_id=fakeable_Electrons_id[0]
  	      DY_l1_pdgid=tightElectrons_pdgid[0]
  	      DY_l2_pdgid=fakeable_Electrons_pdgid[0]
  	      DY_l1_pt=tightElectrons[0].Pt()
  	      DY_l1_eta=tightElectrons[0].Eta()
  	      DY_l1_phi=tightElectrons[0].Phi()
  	      DY_l1_mass=tightElectrons[0].M()
  	      DY_l2_pt=fakeable_Electrons[0].Pt()
  	      DY_l2_eta=fakeable_Electrons[0].Eta()
  	      DY_l2_phi=fakeable_Electrons[0].Phi()
  	      DY_l2_mass=fakeable_Electrons[0].M()
	    else:
              DY_lep1_faketag=True
              DY_lep2_faketag=False
  	      DY_l2_id=tightElectrons_id[0]
  	      DY_l1_id=fakeable_Electrons_id[0]
  	      DY_l2_pdgid=tightElectrons_pdgid[0]
  	      DY_l1_pdgid=fakeable_Electrons_pdgid[0]
  	      DY_l2_pt=tightElectrons[0].Pt()
  	      DY_l2_eta=tightElectrons[0].Eta()
  	      DY_l2_phi=tightElectrons[0].Phi()
  	      DY_l2_mass=tightElectrons[0].M()
  	      DY_l1_pt=fakeable_Electrons[0].Pt()
  	      DY_l1_eta=fakeable_Electrons[0].Eta()
  	      DY_l1_phi=fakeable_Electrons[0].Phi()
  	      DY_l1_mass=fakeable_Electrons[0].M()
	    DY_z_mass=(fakeable_Electrons[0]+tightElectrons[0]).M()
	    DY_z_pt=(fakeable_Electrons[0]+tightElectrons[0]).Pt()
	    DY_z_eta=(fakeable_Electrons[0]+tightElectrons[0]).Eta()
	    DY_z_phi=(fakeable_Electrons[0]+tightElectrons[0]).Phi()
	    DY_drll=fakeable_Electrons[0].DeltaR(tightElectrons[0])

    self.out.fillBranch("DY_region", DY_region)
    self.out.fillBranch("DY_1P1F", DY_1P1F)
    self.out.fillBranch("DY_0P2F", DY_0P2F)
    self.out.fillBranch("DY_lep1_faketag", DY_lep1_faketag)
    self.out.fillBranch("DY_lep2_faketag", DY_lep2_faketag)
    self.out.fillBranch("DY_l1_id", DY_l1_id)
    self.out.fillBranch("DY_l2_id", DY_l2_id)
    self.out.fillBranch("DY_l1_pdgid", DY_l1_pdgid)
    self.out.fillBranch("DY_l2_pdgid", DY_l2_pdgid)
    self.out.fillBranch("DY_l1_pt", DY_l1_pt)
    self.out.fillBranch("DY_l1_pt_raw", DY_l1_pt_raw)
    self.out.fillBranch("DY_l1_eta", DY_l1_eta)
    self.out.fillBranch("DY_l1_phi", DY_l1_phi)
    self.out.fillBranch("DY_l1_mass", DY_l1_mass)
    self.out.fillBranch("DY_l2_pt", DY_l2_pt)
    self.out.fillBranch("DY_l2_pt_raw", DY_l2_pt_raw)
    self.out.fillBranch("DY_l2_eta", DY_l2_eta)
    self.out.fillBranch("DY_l2_phi", DY_l2_phi)
    self.out.fillBranch("DY_l2_mass", DY_l2_mass)
    self.out.fillBranch("DY_z_mass", DY_z_mass)
    self.out.fillBranch("DY_z_mass_raw", DY_z_mass_raw)
    self.out.fillBranch("DY_z_pt", DY_z_pt)
    self.out.fillBranch("DY_z_pt_raw", DY_z_pt_raw)
    self.out.fillBranch("DY_z_eta", DY_z_eta)
    self.out.fillBranch("DY_z_eta_raw", DY_z_eta_raw)
    self.out.fillBranch("DY_z_phi", DY_z_phi)
    self.out.fillBranch("DY_z_phi_raw", DY_z_phi_raw)
    self.out.fillBranch("DY_drll", DY_drll)
    self.out.fillBranch("DY_met", DY_met)
    self.out.fillBranch("DY_met_phi", DY_met_phi)

    if not (ttc_nl or WZ_region >0 or DY_region>0):
      return False

    return True

FakeLep2016 = lambda: FakeLepProducer("2016")
FakeLep2017 = lambda: FakeLepProducer("2017")
FakeLep2018 = lambda: FakeLepProducer("2018")
