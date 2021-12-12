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

def flav_cut(cone_pt):
  if cone_pt<30:
    return 0.304
  if cone_pt>30 and cone_pt<50:
    return (0.304-0.01254*(cone_pt - 30))
  if cone_pt>50:
    return 0.0542

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
    self.out.branch("mt_fixed", "F")
    self.out.branch("met", "F")
    self.out.branch("met_phi", "F")
    self.out.branch("l1_pt", "F")
    self.out.branch("l1_eta", "F")
    self.out.branch("l1_phi", "F")
    self.out.branch("l1_mass", "F")
    self.out.branch("njet_20", "I")
    self.out.branch("njet_25", "I")
    self.out.branch("njet_30", "I")
    self.out.branch("njet_35", "I")
    self.out.branch("njet_20_dr07", "I")
    self.out.branch("njet_25_dr07", "I")
    self.out.branch("njet_30_dr07", "I")
    self.out.branch("njet_35_dr07", "I")
#    self.out.branch("nHad_tau", "I")
    self.out.branch("n_looseB", "I")
    self.out.branch("jet_selection_20", "B")
    self.out.branch("jet_selection_25", "B")
    self.out.branch("jet_selection_30", "B")
    self.out.branch("jet_selection_35", "B")
    self.out.branch("n_looseB_dr07", "I")
    self.out.branch("jet_selection_20_dr07", "B")
    self.out.branch("jet_selection_25_dr07", "B")
    self.out.branch("jet_selection_30_dr07", "B")
    self.out.branch("jet_selection_35_dr07", "B")
    self.out.branch("tightElectrons_id","I",lenVar="nElectron")
    self.out.branch("fakeable_Electrons_id","I",lenVar="nElectron")
    self.out.branch("additional_vetoElectrons_id","I",lenVar="nElectron")
    self.out.branch("tightMuons_id","I",lenVar="nMuon")
    self.out.branch("fakeable_Muons_id","I",lenVar="nMuon")
    self.out.branch("additional_looseMuons_id","I",lenVar="nMuon")
    self.out.branch("muon_genmatch","I",lenVar="nMuon")
    self.out.branch("ele_genmatch","I",lenVar="nElectron")
    self.out.branch("lep_genmatch","I",)
    self.out.branch("muon_conePt","F",lenVar="nMuon")
    self.out.branch("muon_jet_Ptratio","F",lenVar="nMuon")
    self.out.branch("muon_closest_jetid","I",lenVar="nMuon")
    self.out.branch("electron_conePt","F",lenVar="nElectron")
    self.out.branch("electron_jet_Ptratio","F",lenVar="nElectron")
    self.out.branch("electron_closest_jetid","I",lenVar="nElectron")
    self.out.branch("EleHLT1","B")
    self.out.branch("EleHLT2","B")
    self.out.branch("EleHLT3","B")
    self.out.branch("EleHLT4","B")
    self.out.branch("EleHLT5","B")
    self.out.branch("EleHLT6","B")
    self.out.branch("EleHLT7","B")
    self.out.branch("MuHLT1","B")
    self.out.branch("MuHLT2","B")
    self.out.branch("MuHLT3","B")
    self.out.branch("MuHLT4","B")
    self.out.branch("MuHLT5","B")
    self.out.branch("MuHLT6","B")
    self.out.branch("MuHLT7","B")
    self.out.branch("MuHLT8","B")
    self.out.branch("MuHLT9","B")
    self.out.branch("MuHLT10","B")
    self.is_mc = bool(inputTree.GetBranch("GenJet_pt"))
    self.has_EleHLT1 = bool(inputTree.GetBranch("HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30"))
    self.has_EleHLT2 = bool(inputTree.GetBranch("HLT_Ele8_CaloIdM_TrackIdM_PFJet30"))
    self.has_EleHLT3 = bool(inputTree.GetBranch("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30"))
    self.has_EleHLT4 = bool(inputTree.GetBranch("HLT_Ele15_CaloIdL_TrackIdL_IsoVL_PFJet30"))
    self.has_EleHLT5 = bool(inputTree.GetBranch("HLT_Ele17_CaloIdM_TrackIdM_PFJet30"))
    self.has_EleHLT6 = bool(inputTree.GetBranch("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30"))
    self.has_EleHLT7 = bool(inputTree.GetBranch("HLT_Ele23_CaloIdM_TrackIdM_PFJet30"))
    self.has_MuHLT1 = bool(inputTree.GetBranch("HLT_Mu3_PFJet40"))
    self.has_MuHLT2 = bool(inputTree.GetBranch("HLT_Mu8"))
    self.has_MuHLT3 = bool(inputTree.GetBranch("HLT_Mu12"))
    self.has_MuHLT4 = bool(inputTree.GetBranch("HLT_Mu15"))
    self.has_MuHLT5 = bool(inputTree.GetBranch("HLT_Mu17"))
    self.has_MuHLT6 = bool(inputTree.GetBranch("HLT_Mu17_TrkIsoVVL"))
    self.has_MuHLT7 = bool(inputTree.GetBranch("HLT_Mu19"))
    self.has_MuHLT8 = bool(inputTree.GetBranch("HLT_Mu19_TrkIsoVVL"))
    self.has_MuHLT9 = bool(inputTree.GetBranch("HLT_Mu20"))
    self.has_MuHLT10 = bool(inputTree.GetBranch("HLT_Mu27"))

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
    if not event.nJet>0: return False

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
    muon_genmatch = []

    muon_conePt = []
    muon_jet_Ptratio = []
    muon_closest_jetid = []
    jet_v4_temp=TLorentzVector()

    for imu in range(0, event.nMuon):
      dr_mu_jet=99.
      muon_closest_jetid_temp=-1
      muon_v4_temp.SetPtEtaPhiM(event.Muon_corrected_pt[imu], muons[imu].eta, muons[imu].phi, muons[imu].mass)
      for ijet in range(0, event.nJet):
	jet_v4_temp.SetPtEtaPhiM(event.Jet_pt_nom[ijet],event.Jet_eta[ijet],event.Jet_phi[ijet],event.Jet_mass_nom[ijet])
	if muon_v4_temp.DeltaR(jet_v4_temp)<dr_mu_jet:
	  dr_mu_jet=muon_v4_temp.DeltaR(jet_v4_temp)
	  muon_closest_jetid_temp=ijet

      if dr_mu_jet<0.4:
	muon_conePt.append(0.85*event.Jet_pt_nom[muon_closest_jetid_temp])
	muon_jet_Ptratio.append(event.Muon_corrected_pt[imu]/(0.85*event.Jet_pt_nom[muon_closest_jetid_temp]))
	muon_closest_jetid.append(muon_closest_jetid_temp)
      else:
	muon_conePt.append(event.Muon_corrected_pt[imu]*(1+event.Muon_miniPFRelIso_all[imu]))
	muon_jet_Ptratio.append(1./(1+event.Muon_miniPFRelIso_all[imu]))
	muon_closest_jetid.append(muon_closest_jetid_temp)


    for imu in range(0, event.nMuon):
      if not self.is_mc: muon_genmatch.append(0)
      else:
	if (muons[imu].genPartFlav==1 or muons[imu].genPartFlav==15): muon_genmatch.append(1)
        else: muon_genmatch.append(0)
      # following cuts are preseletion for MVA muon ID
      if abs(muons[imu].eta)>2.4 or muons[imu].sip3d>8 or abs(muons[imu].dxy)>0.05 or abs(muons[imu].dz)>0.1 or muons[imu].miniPFRelIso_all>0.4: continue
      if not (muons[imu].looseId and event.Muon_corrected_pt[imu]>10 and muon_conePt[imu]>10):continue
      if muons[imu].mediumId and muons[imu].tightCharge==2:
        if muons[imu].mvaTTH>-0.2:
 	  if (muon_conePt[imu]>20):
 	    muon_v4_temp.SetPtEtaPhiM(muon_conePt[imu], muons[imu].eta, muons[imu].phi, muons[imu].mass)
 	    tightMuons.append(muon_v4_temp.Clone())
 	    tightMuons_pdgid.append(muons[imu].pdgId)
 	    tightMuons_id.append(imu)
 	  if (muon_conePt[imu]<20 and event.Muon_corrected_pt[imu]>10):
 	    muon_v4_temp.SetPtEtaPhiM(event.Muon_corrected_pt[imu], muons[imu].eta, muons[imu].phi, muons[imu].mass)
  	    additional_looseMuons.append(muon_v4_temp.Clone())
  	    additional_looseMuons_pdgid.append(muons[imu].pdgId)
  	    additional_looseMuons_id.append(imu)
        else:
	  if muon_jet_Ptratio[imu]>0.5 and event.Jet_btagDeepFlavB[muon_closest_jetid[imu]]<flav_cut(muon_conePt[imu]):
	    if (muon_conePt[imu]>20):
              muon_v4_temp.SetPtEtaPhiM(muon_conePt[imu], muons[imu].eta, muons[imu].phi, muons[imu].mass)
              fakeable_Muons.append(muon_v4_temp.Clone())
              fakeable_Muons_pdgid.append(muons[imu].pdgId)
              fakeable_Muons_id.append(imu)
            else:
	      if event.Muon_corrected_pt[imu]>10:
                muon_v4_temp.SetPtEtaPhiM(event.Muon_corrected_pt[imu], muons[imu].eta, muons[imu].phi, muons[imu].mass)
                additional_looseMuons.append(muon_v4_temp.Clone())
                additional_looseMuons_pdgid.append(muons[imu].pdgId)
                additional_looseMuons_id.append(imu)
	  else:
	    if event.Muon_corrected_pt[imu]>10:
              muon_v4_temp.SetPtEtaPhiM(event.Muon_corrected_pt[imu], muons[imu].eta, muons[imu].phi, muons[imu].mass)
              additional_looseMuons.append(muon_v4_temp.Clone())
              additional_looseMuons_pdgid.append(muons[imu].pdgId)
              additional_looseMuons_id.append(imu)
      else:
	if event.Muon_corrected_pt[imu]>10:
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
    tightMuons_id.extend(np.zeros(event.nMuon-len(tightMuons_id),int)-1)
    fakeable_Muons_id.extend(np.zeros(event.nMuon-len(fakeable_Muons_id),int)-1)
    additional_looseMuons_id.extend(np.zeros(event.nMuon-len(additional_looseMuons_id),int)-1)
    self.out.fillBranch("tightMuons_id", tightMuons_id)
    self.out.fillBranch("fakeable_Muons_id", fakeable_Muons_id)
    self.out.fillBranch("additional_looseMuons_id", additional_looseMuons_id)
    self.out.fillBranch("muon_genmatch", muon_genmatch)
    self.out.fillBranch("muon_conePt", muon_conePt)
    self.out.fillBranch("muon_jet_Ptratio", muon_jet_Ptratio)
    self.out.fillBranch("muon_closest_jetid", muon_closest_jetid)

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
    ele_genmatch = []

    electron_conePt = []
    electron_jet_Ptratio = []
    electron_closest_jetid = []
    electron_closest_jetid_temp=-1

    for iele in range(0, event.nElectron):
      dr_ele_jet=99.
      electron_v4_temp.SetPtEtaPhiM(electrons[iele].pt, electrons[iele].eta, electrons[iele].phi, electrons[iele].mass)
      for ijet in range(0, event.nJet):
	jet_v4_temp.SetPtEtaPhiM(event.Jet_pt_nom[ijet],event.Jet_eta[ijet],event.Jet_phi[ijet],event.Jet_mass_nom[ijet])
	if electron_v4_temp.DeltaR(jet_v4_temp)<dr_ele_jet:
	  dr_ele_jet=electron_v4_temp.DeltaR(jet_v4_temp)
	  electron_closest_jetid_temp=ijet

      if dr_ele_jet<0.4:
	electron_conePt.append(0.85*event.Jet_pt_nom[electron_closest_jetid_temp])
	electron_jet_Ptratio.append(electrons[iele].pt/(0.85*event.Jet_pt_nom[electron_closest_jetid_temp]))
	electron_closest_jetid.append(electron_closest_jetid_temp)
      else:
	electron_conePt.append(electrons[iele].pt*(1+event.Electron_miniPFRelIso_all[iele]))
	electron_jet_Ptratio.append(1./(1+event.Electron_miniPFRelIso_all[iele]))
	electron_closest_jetid.append(electron_closest_jetid_temp)

    for iele in range(0, event.nElectron):
      if not self.is_mc: ele_genmatch.append(0)
      else:
	if (electrons[iele].genPartFlav==1 or electrons[iele].genPartFlav==15): ele_genmatch.append(1)
        else: ele_genmatch.append(0)
      # following cuts are preseletion for MVA electron ID
      if abs(electrons[iele].eta)>2.5 or electrons[iele].sip3d>8 or abs(electrons[iele].dxy)>0.05 or abs(electrons[iele].dz)>0.1 or electrons[iele].miniPFRelIso_all>0.4: continue
      if not (electrons[iele].mvaFall17V2noIso_WPL and electrons[iele].lostHits<2 and electrons[iele].pt>10 and electron_conePt[iele]>10 and electrons[iele].convVeto):continue
      if ((abs(electrons[iele].deltaEtaSC+electrons[iele].eta)<1.479 and electrons[iele].sieie<0.011) or (abs(electrons[iele].deltaEtaSC+electrons[iele].eta)>1.479 and electrons[iele].sieie<0.03)) and electrons[iele].hoe<0.1 and electrons[iele].eInvMinusPInv>-0.04:
	if electrons[iele].mvaTTH>0.25 and electrons[iele].tightCharge==2:
	  if (electron_conePt[iele]>20):
	    electron_v4_temp.SetPtEtaPhiM(electron_conePt[iele], electrons[iele].eta, electrons[iele].phi, electrons[iele].mass)
            tightElectrons.append(electron_v4_temp.Clone())
            tightElectrons_pdgid.append(electrons[iele].pdgId)
            tightElectrons_id.append(iele)
	  if (electron_conePt[iele]<20 and electrons[iele].pt>10):
	    electron_v4_temp.SetPtEtaPhiM(electrons[iele].pt, electrons[iele].eta, electrons[iele].phi, electrons[iele].mass)
            additional_vetoElectrons.append(electron_v4_temp.Clone())
            additional_vetoElectrons_pdgid.append(electrons[iele].pdgId)
            additional_vetoElectrons_id.append(iele)
	else:
	  if electrons[iele].mvaFall17V2noIso_WP90 and electron_jet_Ptratio[iele]>0.6 and event.Jet_btagDeepFlavB[electron_closest_jetid[iele]]<flav_cut(electron_conePt[iele]):
	    if (electron_conePt[iele]>20):
              electron_v4_temp.SetPtEtaPhiM(electron_conePt[iele], electrons[iele].eta, electrons[iele].phi, electrons[iele].mass)
              fakeable_Electrons.append(electron_v4_temp.Clone())
              fakeable_Electrons_pdgid.append(electrons[iele].pdgId)
              fakeable_Electrons_id.append(iele)
	    else:
	      if electrons[iele].pt>10:
                electron_v4_temp.SetPtEtaPhiM(electrons[iele].pt, electrons[iele].eta, electrons[iele].phi, electrons[iele].mass)
                additional_vetoElectrons.append(electron_v4_temp.Clone())
                additional_vetoElectrons_pdgid.append(electrons[iele].pdgId)
                additional_vetoElectrons_id.append(iele)
	  else:
	    if electrons[iele].pt>10:
              electron_v4_temp.SetPtEtaPhiM(electrons[iele].pt, electrons[iele].eta, electrons[iele].phi, electrons[iele].mass)
              additional_vetoElectrons.append(electron_v4_temp.Clone())
              additional_vetoElectrons_pdgid.append(electrons[iele].pdgId)
              additional_vetoElectrons_id.append(iele)
      else:
	if electrons[iele].pt>10:
	  electron_v4_temp.SetPtEtaPhiM(electrons[iele].pt, electrons[iele].eta, electrons[iele].phi, electrons[iele].mass)
          additional_vetoElectrons.append(electron_v4_temp.Clone())
          additional_vetoElectrons_pdgid.append(electrons[iele].pdgId)
          additional_vetoElectrons_id.append(iele)

    n_tight_ele = len(tightElectrons)
    n_fakeable_ele = len(fakeable_Electrons)
    n_loose_ele = len(additional_vetoElectrons)
    self.out.fillBranch("n_tight_ele", n_tight_ele)
    self.out.fillBranch("n_fakeable_ele", n_fakeable_ele)
    self.out.fillBranch("n_loose_ele", n_loose_ele)
    tightElectrons_id.extend(np.zeros(event.nElectron-len(tightElectrons_id),int)-1)
    fakeable_Electrons_id.extend(np.zeros(event.nElectron-len(fakeable_Electrons_id),int)-1)
    additional_vetoElectrons_id.extend(np.zeros(event.nElectron-len(additional_vetoElectrons_id),int)-1)
    self.out.fillBranch("tightElectrons_id", tightElectrons_id)
    self.out.fillBranch("fakeable_Electrons_id", fakeable_Electrons_id)
    self.out.fillBranch("additional_vetoElectrons_id", additional_vetoElectrons_id)
    self.out.fillBranch("ele_genmatch", ele_genmatch)
    self.out.fillBranch("electron_conePt", electron_conePt)
    self.out.fillBranch("electron_jet_Ptratio", electron_jet_Ptratio)
    self.out.fillBranch("electron_closest_jetid", electron_closest_jetid)

#    tightLeptons = tightMuons + tightElectrons
#    tightLeptons.sort(key=lambda x: x.Pt(), reverse=True)
#    fakeableLeptons = fakeable_Muons + fakeable_Electrons
#    fakeableLeptons.sort(key=lambda x: x.Pt(), reverse=True)
#
#    tau_v4_temp=TLorentzVector()
#    taus = Collection(event, 'Tau')
#    nHad_tau=0
#    Had_tau_id=[]
#    for itau in range(0, event.nTau):
#      tau_v4_temp.SetPtEtaPhiM(taus[itau].pt, taus[itau].eta, taus[itau].phi, taus[itau].mass)
#      pass_tau_lep_Dr=1
#      if taus[itau].pt>20 and abs(taus[itau].eta)<2.3 and taus[itau].idDecayModeNewDMs and taus[itau].idDeepTau2017v2p1VSe>=4 and taus[itau].idDeepTau2017v2p1VSjet>=4 and taus[itau].idDeepTau2017v2p1VSmu>=1:
#        for ilep in range(0,len(tightLeptons)):
#          if tau_v4_temp.DeltaR(tightLeptons[ilep])<0.4:pass_tau_lep_Dr=0
#        for ilep in range(0,len(fakeableLeptons)):
#          if tau_v4_temp.DeltaR(fakeableLeptons[ilep])<0.4:pass_tau_lep_Dr=0
#        if pass_tau_lep_Dr:
#          nHad_tau=nHad_tau+1
#          Had_tau_id.append(itau)
#    self.out.fillBranch("nHad_tau", nHad_tau)

    mt=-99
    mt_fixed=-99
    met=-99
    met_phi=-99
    l1_pt=-99
    l1_eta=-99
    l1_phi=-99
    l1_mass=-99

    if not n_tight_muon+n_fakeable_muon + n_tight_ele+n_fakeable_ele == 1: return False
    if not n_loose_muon + n_loose_ele == 0: return False
    if n_tight_muon==1:lep_genmatch=muon_genmatch[tightMuons_id[0]]
    if n_fakeable_muon==1:lep_genmatch=muon_genmatch[fakeable_Muons_id[0]]
    if n_tight_ele==1:lep_genmatch=ele_genmatch[tightElectrons_id[0]]
    if n_fakeable_ele==1:lep_genmatch=ele_genmatch[fakeable_Electrons_id[0]]
    self.out.fillBranch("lep_genmatch", lep_genmatch)

    # tight or fakeable leptons collection
    Leptons = tightMuons + tightElectrons + fakeable_Muons + fakeable_Electrons
    l1_pt=Leptons[0].Pt()
    l1_eta=Leptons[0].Eta()
    l1_phi=Leptons[0].Phi()
    l1_mass=Leptons[0].M()

    EleHLT1=False
    EleHLT2=False
    EleHLT3=False
    EleHLT4=False
    EleHLT5=False
    EleHLT6=False
    EleHLT7=False
    MuHLT1=False
    MuHLT2=False
    MuHLT3=False
    MuHLT4=False
    MuHLT5=False
    MuHLT6=False
    MuHLT7=False
    MuHLT8=False
    MuHLT9=False
    MuHLT10=False

    if (n_tight_muon+n_fakeable_muon)==1:
      if self.has_MuHLT1:
	if event.HLT_Mu3_PFJet40:MuHLT1=True
      if self.has_MuHLT2:
	if event.HLT_Mu8:MuHLT2=True
      if self.has_MuHLT3:
	if event.HLT_Mu12:MuHLT3=True
      if self.has_MuHLT4:
	if event.HLT_Mu15:MuHLT4=True
      if self.has_MuHLT5:
	if event.HLT_Mu17:MuHLT5=True
      if self.has_MuHLT6:
	if event.HLT_Mu17_TrkIsoVVL:MuHLT6=True
      if self.has_MuHLT7:
	if event.HLT_Mu19:MuHLT7=True
      if self.has_MuHLT8:
	if event.HLT_Mu19_TrkIsoVVL:MuHLT8=True
      if self.has_MuHLT9:
	if event.HLT_Mu20:MuHLT9=True
      if self.has_MuHLT10:
	if event.HLT_Mu27:MuHLT10=True

      if not (MuHLT1 or MuHLT2 or MuHLT3 or MuHLT4 or MuHLT5 or MuHLT6 or MuHLT7 or MuHLT8 or MuHLT9 or MuHLT10): return False

    if (n_tight_ele+n_fakeable_ele)==1:
      if self.has_EleHLT1:
	if event.HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30:EleHLT1=True
      if self.has_EleHLT2:
	if event.HLT_Ele8_CaloIdM_TrackIdM_PFJet30:EleHLT2=True
      if self.has_EleHLT3:
	if event.HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30:EleHLT3=True
      if self.has_EleHLT4:
	if event.HLT_Ele15_CaloIdL_TrackIdL_IsoVL_PFJet30:EleHLT4=True
      if self.has_EleHLT5:
	if event.HLT_Ele17_CaloIdM_TrackIdM_PFJet30:EleHLT5=True
      if self.has_EleHLT6:
	if event.HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30:EleHLT6=True
      if self.has_EleHLT7:
	if event.HLT_Ele23_CaloIdM_TrackIdM_PFJet30:EleHLT7=True
      if not (EleHLT1 or EleHLT2 or EleHLT3 or EleHLT4 or EleHLT5 or EleHLT6 or EleHLT7): return False

    if self.is_mc:
      met=event.MET_T1Smear_pt
      met_phi=event.MET_T1Smear_phi
    else:
      met=event.MET_T1_pt
      met_phi=event.MET_T1_phi
#    if met>40: return False

    mt = sqrt(2*l1_pt*met*(1 - cos(met_phi - l1_phi)))
    mt_fixed = sqrt(2*35*met*(1 - cos(met_phi - l1_phi)))
#    if mt>40: return False

    jet_selection_20=False
    jet_selection_25=False
    jet_selection_30=False
    jet_selection_35=False
    njet_20=0
    njet_25=0
    njet_30=0
    njet_35=0
    n_looseB=0
    jet_selection_20_dr07=False
    jet_selection_25_dr07=False
    jet_selection_30_dr07=False
    jet_selection_35_dr07=False
    njet_20_dr07=0
    njet_25_dr07=0
    njet_30_dr07=0
    njet_35_dr07=0
    n_looseB_dr07=0
    jets = Collection(event, 'Jet')

    # require DeltaR between Jets and tight leptons greater than 0.4
    for ijet in range(0, event.nJet):

#      jet_is_tau=0
#      if nHad_tau>0:
#        for ita in Had_tau_id:
#          if ijet==event.Tau_jetIdx[ita]:jet_is_tau=1
#      if jet_is_tau:continue

      #if abs(jets[ijet].eta)>4.7 or jets[ijet].jetId <2: 
      if abs(jets[ijet].eta)>4.7: 
	continue
      jet_v4_temp.SetPtEtaPhiM(event.Jet_pt_nom[ijet],event.Jet_eta[ijet],event.Jet_phi[ijet],event.Jet_mass_nom[ijet])
      if jet_v4_temp.DeltaR(Leptons[0])>0.4:
        if event.Jet_pt_nom[ijet]>35:njet_35=njet_35+1
        if event.Jet_pt_nom[ijet]>30:njet_30=njet_30+1
        if event.Jet_pt_nom[ijet]>25:njet_25=njet_25+1
        if event.Jet_pt_nom[ijet]>20:njet_20=njet_20+1
        if abs(jets[ijet].eta)<2.4 and jets[ijet].btagDeepB > 0.0532:n_looseB=n_looseB+1
      if jet_v4_temp.DeltaR(Leptons[0])>0.7:
        if event.Jet_pt_nom[ijet]>35:njet_35_dr07=njet_35_dr07+1
        if event.Jet_pt_nom[ijet]>30:njet_30_dr07=njet_30_dr07+1
        if event.Jet_pt_nom[ijet]>25:njet_25_dr07=njet_25_dr07+1
        if event.Jet_pt_nom[ijet]>20:njet_20_dr07=njet_20_dr07+1
        if abs(jets[ijet].eta)<2.4 and jets[ijet].btagDeepB > 0.0532:n_looseB_dr07=n_looseB_dr07+1

    if njet_35>0: jet_selection_35=True
    if njet_30>0: jet_selection_30=True
    if njet_25>0: jet_selection_25=True
    if njet_20>0: jet_selection_20=True
    if njet_35_dr07>0: jet_selection_35_dr07=True
    if njet_30_dr07>0: jet_selection_30_dr07=True
    if njet_25_dr07>0: jet_selection_25_dr07=True
    if njet_20_dr07>0: jet_selection_20_dr07=True

    if not (jet_selection_20 or jet_selection_25 or jet_selection_30 or jet_selection_35 or jet_selection_20_dr07 or jet_selection_25_dr07 or jet_selection_30_dr07 or jet_selection_35_dr07):
      return False

    self.out.fillBranch("mt", mt)
    self.out.fillBranch("met", met)
    self.out.fillBranch("met_phi", met_phi)
    self.out.fillBranch("l1_pt", l1_pt)
    self.out.fillBranch("l1_eta", l1_eta)
    self.out.fillBranch("l1_phi", l1_phi)
    self.out.fillBranch("l1_mass", l1_mass)
    self.out.fillBranch("njet_20", njet_20)
    self.out.fillBranch("njet_25", njet_25)
    self.out.fillBranch("njet_30", njet_30)
    self.out.fillBranch("njet_35", njet_35)
    self.out.fillBranch("n_looseB", n_looseB)
    self.out.fillBranch("jet_selection_20", jet_selection_20)
    self.out.fillBranch("jet_selection_25", jet_selection_25)
    self.out.fillBranch("jet_selection_30", jet_selection_30)
    self.out.fillBranch("jet_selection_35", jet_selection_35)
    self.out.fillBranch("njet_20_dr07", njet_20_dr07)
    self.out.fillBranch("njet_25_dr07", njet_25_dr07)
    self.out.fillBranch("njet_30_dr07", njet_30_dr07)
    self.out.fillBranch("njet_35_dr07", njet_35_dr07)
    self.out.fillBranch("n_looseB_dr07", n_looseB_dr07)
    self.out.fillBranch("jet_selection_20_dr07", jet_selection_20_dr07)
    self.out.fillBranch("jet_selection_25_dr07", jet_selection_25_dr07)
    self.out.fillBranch("jet_selection_30_dr07", jet_selection_30_dr07)
    self.out.fillBranch("jet_selection_35_dr07", jet_selection_35_dr07)
    self.out.fillBranch("EleHLT1",EleHLT1)
    self.out.fillBranch("EleHLT2",EleHLT2)
    self.out.fillBranch("EleHLT3",EleHLT3)
    self.out.fillBranch("EleHLT4",EleHLT4)
    self.out.fillBranch("EleHLT5",EleHLT5)
    self.out.fillBranch("EleHLT6",EleHLT6)
    self.out.fillBranch("EleHLT7",EleHLT7)
    self.out.fillBranch("MuHLT1",MuHLT1)
    self.out.fillBranch("MuHLT2",MuHLT2)
    self.out.fillBranch("MuHLT3",MuHLT3)
    self.out.fillBranch("MuHLT4",MuHLT4)
    self.out.fillBranch("MuHLT5",MuHLT5)
    self.out.fillBranch("MuHLT6",MuHLT6)
    self.out.fillBranch("MuHLT7",MuHLT7)
    self.out.fillBranch("MuHLT8",MuHLT8)
    self.out.fillBranch("MuHLT9",MuHLT9)
    self.out.fillBranch("MuHLT10",MuHLT10)

    return True

FakeRate2016 = lambda: FakeRateProducer("2016")
FakeRate2017 = lambda: FakeRateProducer("2017")
FakeRate2018 = lambda: FakeRateProducer("2018")
