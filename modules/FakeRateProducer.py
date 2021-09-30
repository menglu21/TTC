import ROOT
from ROOT import TLorentzVector
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection 
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

import math
import os
import numpy as np
from numpy import sign
from numpy import sqrt, cos

class FakeRateProducer(Module):
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
    self.out.branch("mt", "F")
    self.out.branch("met", "F")
    self.out.branch("met_phi", "F")
    self.out.branch("l1_pt", "F")
    self.out.branch("l1_eta", "F")
    self.out.branch("l1_phi", "F")
    self.out.branch("l1_mass", "F")
    self.out.branch("jet_selection_20", "B")
    self.out.branch("jet_selection_25", "B")
    self.out.branch("jet_selection_30", "B")
    self.out.branch("jet_selection_35", "B")
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
    if ((event.nMuon + event.nElectron) < 1): return False

    # Muon selection: tight cut-based ID + tight PF iso, or loose cut-based ID + loose PF iso, with pt > 20 GeV
    muons = Collection(event, 'Muon')
    muon_v4_temp=TLorentzVector()
    tightMuons = []
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
        if (muons[imu].pfRelIso04_all<0.15 and abs(muons[imu].eta)<2.4 and muons[imu].tightCharge==2 and event.Muon_corrected_pt[imu]>20):
          muon_v4_temp.SetPtEtaPhiM(event.Muon_corrected_pt[imu], muons[imu].eta, muons[imu].phi, muons[imu].mass)
          tightMuons.append(muon_v4_temp.Clone())
          tightMuons_pdgid.append(muons[imu].pdgId)
          tightMuons_id.append(imu)
	if (muons[imu].pfRelIso04_all<0.4 and muons[imu].pfRelIso04_all>0.2 and abs(muons[imu].eta)<2.4 and event.Muon_corrected_pt[imu]>20):
	  if self.is_mc and (muons[imu].genPartFlav==1 or muons[imu].genPartFlav==15):
            muon_v4_temp.SetPtEtaPhiM(event.Muon_corrected_pt[imu], muons[imu].eta, muons[imu].phi, muons[imu].mass)
            fakeable_Muons.append(muon_v4_temp.Clone())
            fakeable_Muons_pdgid.append(muons[imu].pdgId)
            fakeable_Muons_id.append(imu)
	  if not self.is_mc:
	    muon_v4_temp.SetPtEtaPhiM(event.Muon_corrected_pt[imu], muons[imu].eta, muons[imu].phi, muons[imu].mass)
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
    tightElectrons = []
    tightElectrons_pdgid = []
    tightElectrons_id = []
    fakeable_Electrons = []
    fakeable_Electrons_pdgid = []
    fakeable_Electrons_id = []
    additional_vetoElectrons = []
    additional_vetoElectrons_pdgid = []
    additional_vetoElectrons_id = []
    for iele in range(0, event.nElectron):
      if not ((abs(electrons[iele].eta+electrons[iele].deltaEtaSC) <1.4442 and abs(electrons[iele].dxy)<0.05 and abs(electrons[iele].dz)<0.1) or (abs(electrons[iele].eta + electrons[iele].deltaEtaSC)>1.566 and abs(electrons[iele].eta + electrons[iele].deltaEtaSC)<2.4 and abs(electrons[iele].dxy)<0.1 and abs(electrons[iele].dz)<0.2)): continue
      if (electrons[iele].cutBased==4 and electrons[iele].tightCharge==2 and electrons[iele].pt>20):
        electron_v4_temp.SetPtEtaPhiM(electrons[iele].pt, electrons[iele].eta, electrons[iele].phi, electrons[iele].mass)
        tightElectrons.append(electron_v4_temp.Clone())
        tightElectrons_pdgid.append(electrons[iele].pdgId)
        tightElectrons_id.append(iele)
      if (electrons[iele].cutBased>1 and electrons[iele].cutBased<4 and electrons[iele].pt>15):
        electron_v4_temp.SetPtEtaPhiM(electrons[iele].pt, electrons[iele].eta, electrons[iele].phi, electrons[iele].mass)
        additional_vetoElectrons.append(electron_v4_temp.Clone())
        additional_vetoElectrons_pdgid.append(electrons[iele].pdgId)
        additional_vetoElectrons_id.append(iele)
      if (electrons[iele].cutBased==1):
	if self.is_mc and (electrons[iele].genPartFlav==1 or electrons[iele].genPartFlav==15 and electrons[iele].pt>20):
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

    mt=-99
    met=-99
    met_phi=-99
    l1_pt=-99
    l1_eta=-99
    l1_phi=-99
    l1_mass=-99

    if not n_tight_muon+n_fakeable_muon + n_tight_ele+n_fakeable_ele == 1: return False
    if not n_loose_muon + n_loose_ele == 0: return False

    # tight or fakeable leptons collection
    Leptons = tightMuons + tightElectrons + fakeable_Muons + fakeable_Electrons
    l1_pt=Leptons[0].Pt()
    l1_eta=Leptons[0].Eta()
    l1_phi=Leptons[0].Phi()
    l1_mass=Leptons[0].M()

    if self.is_mc:
      met=event.MET_T1Smear_pt
      met_phi=event.MET_T1Smear_phi
    else:
      met=event.MET_T1_pt
      met_phi=event.MET_T1_phi
    if met>40: return False

    mt = sqrt(2*l1_pt*met*(1 - cos(met_phi - l1_phi)))
    if mt>40: return False

    jet_selection_20=False
    jet_selection_25=False
    jet_selection_30=False
    jet_selection_35=False
    jets = Collection(event, 'Jet')

    # require DeltaR between Jets and tight leptons greater than 0.4
    jet_v4_temp=TLorentzVector()
    for ijet in range(0, event.nJet):
      if abs(jets[ijet].eta)>4.7: 
	continue
      jet_v4_temp.SetPtEtaPhiM(event.Jet_pt_nom[ijet],event.Jet_eta[ijet],event.Jet_phi[ijet],event.Jet_mass_nom[ijet])
      if jet_v4_temp.DeltaR(Leptons[0])>0.4 and event.Jet_pt_nom[ijet]>20:jet_selection_20=True
      if jet_v4_temp.DeltaR(Leptons[0])>0.4 and event.Jet_pt_nom[ijet]>25:jet_selection_25=True
      if jet_v4_temp.DeltaR(Leptons[0])>0.4 and event.Jet_pt_nom[ijet]>30:jet_selection_30=True
      if jet_v4_temp.DeltaR(Leptons[0])>0.4 and event.Jet_pt_nom[ijet]>35:jet_selection_35=True

    if not (jet_selection_20 or jet_selection_25 or jet_selection_30 or jet_selection_35):
      return False

    self.out.fillBranch("mt", mt)
    self.out.fillBranch("met", met)
    self.out.fillBranch("met_phi", met_phi)
    self.out.fillBranch("l1_pt", l1_pt)
    self.out.fillBranch("l1_eta", l1_eta)
    self.out.fillBranch("l1_phi", l1_phi)
    self.out.fillBranch("l1_mass", l1_mass)
    self.out.fillBranch("jet_selection_20", jet_selection_20)
    self.out.fillBranch("jet_selection_25", jet_selection_25)
    self.out.fillBranch("jet_selection_30", jet_selection_30)
    self.out.fillBranch("jet_selection_35", jet_selection_35)

    return True

FakeRate2016 = lambda: FakeRateProducer("2016")
FakeRate2017 = lambda: FakeRateProducer("2017")
FakeRate2018 = lambda: FakeRateProducer("2018")
