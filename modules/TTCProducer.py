import ROOT
from ROOT import TLorentzVector
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection 
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

import math
import os
import numpy as np
from numpy import sign

def flav_cut(cone_pt):
  if cone_pt<30:
    return 0.304
  if cone_pt>30 and cone_pt<50:
    return (0.304-0.01254*(cone_pt - 30))
  if cone_pt>50:
    return 0.0542

class TTCProducer(Module):
  def __init__(self , year):
    self.year = year
  def beginJob(self):
    pass
  def endJob(self):
    pass
  def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
    self.out = wrappedOutputTree
    self.out.branch("HLT_passEle32WPTight", "I")
    self.out.branch("lhe_nlepton", "I")
    self.out.branch("n_tight_muon", "I")
    self.out.branch("n_fakeable_muon", "I")
    self.out.branch("n_loose_muon", "I")
    self.out.branch("n_tight_ele", "I")
    self.out.branch("n_fakeable_ele", "I")
    self.out.branch("n_loose_ele", "I")
    self.out.branch("n_tight_jet", "I")
    self.out.branch("n_bjet_DeepB", "I")
    self.out.branch("n_cjet_DeepB_medium", "I")
    self.out.branch("n_cjet_DeepB_loose", "I")
    self.out.branch("HT", "F")
    self.out.branch("nHad_tau", "I")
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
    self.out.branch("mj1j2", "F")
    self.out.branch("mj1j3", "F")
    self.out.branch("mj1j4", "F")
    self.out.branch("mj2j3", "F")
    self.out.branch("mj2j4", "F")
    self.out.branch("mj3j4", "F")
    self.out.branch("mj1j2j3", "F")
    self.out.branch("mj1j2j4", "F")
    self.out.branch("mj2j3j4", "F")
    self.out.branch("mj1j2j3j4", "F")
    self.out.branch("drj1j2", "F")
    self.out.branch("drj1j3", "F")
    self.out.branch("drj1j4", "F")
    self.out.branch("drj2j3", "F")
    self.out.branch("drj2j4", "F")
    self.out.branch("drj3j4", "F")
    self.out.branch("DeepB_j1_pt", "F")
    self.out.branch("DeepB_j1_eta", "F")
    self.out.branch("DeepB_j1_phi", "F")
    self.out.branch("DeepB_j1_mass", "F")
    self.out.branch("DeepB_j2_pt", "F")
    self.out.branch("DeepB_j2_eta", "F")
    self.out.branch("DeepB_j2_phi", "F")
    self.out.branch("DeepB_j2_mass", "F")
    self.out.branch("DeepC_loose_j1_pt", "F")
    self.out.branch("DeepC_loose_j1_eta", "F")
    self.out.branch("DeepC_loose_j1_phi", "F")
    self.out.branch("DeepC_loose_j1_mass", "F")
    self.out.branch("DeepC_loose_j2_pt", "F")
    self.out.branch("DeepC_loose_j2_eta", "F")
    self.out.branch("DeepC_loose_j2_phi", "F")
    self.out.branch("DeepC_loose_j2_mass", "F")
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
    self.out.branch("ttc_2P0F", "B")
    self.out.branch("ttc_1P1F", "B")
    self.out.branch("ttc_0P2F", "B")
    self.out.branch("ttc_lep1_faketag", "B")
    self.out.branch("ttc_lep2_faketag", "B")
    self.out.branch("ttc_region", "I")
    self.out.branch("ttc_l1_id", "I")
    self.out.branch("ttc_l2_id", "I")
    self.out.branch("ttc_l1_isprompt", "I")
    self.out.branch("ttc_l2_isprompt", "I")
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
    self.out.branch("OPS_region", "I")
    self.out.branch("OPS_2P0F", "B")
    self.out.branch("OPS_1P1F", "B")
    self.out.branch("OPS_0P2F", "B")
    self.out.branch("OPS_lep1_faketag", "B")
    self.out.branch("OPS_lep2_faketag", "B")
    self.out.branch("OPS_l1_id", "I")
    self.out.branch("OPS_l2_id", "I")
    self.out.branch("OPS_l1_isprompt", "I")
    self.out.branch("OPS_l2_isprompt", "I")
    self.out.branch("OPS_l1_pdgid", "I")
    self.out.branch("OPS_l2_pdgid", "I")
    self.out.branch("OPS_l1_pt", "F")
    self.out.branch("OPS_l1_pt_raw", "F")
    self.out.branch("OPS_l1_eta", "F")
    self.out.branch("OPS_l1_phi", "F")
    self.out.branch("OPS_l1_mass", "F")
    self.out.branch("OPS_l2_pt", "F")
    self.out.branch("OPS_l2_pt_raw", "F")
    self.out.branch("OPS_l2_eta", "F")
    self.out.branch("OPS_l2_phi", "F")
    self.out.branch("OPS_l2_mass", "F")
    self.out.branch("OPS_z_mass", "F")
    self.out.branch("OPS_z_mass_raw", "F")
    self.out.branch("OPS_z_pt", "F")
    self.out.branch("OPS_z_pt_raw", "F")
    self.out.branch("OPS_z_eta", "F")
    self.out.branch("OPS_z_eta_raw", "F")
    self.out.branch("OPS_z_phi", "F")
    self.out.branch("OPS_z_phi_raw", "F")
    self.out.branch("OPS_drll", "F")
    self.out.branch("tightJets_id_in24","I",lenVar="nJet")
    self.out.branch("tightJets_id_in47","I",lenVar="nJet")
    self.out.branch("tightJets_b_DeepCSVmedium_id","I",lenVar="nJet")
    self.out.branch("tightJets_c_DeepCSVmedium_id","I",lenVar="nJet")
    self.out.branch("tightJets_c_DeepCSVloose_id","I",lenVar="nJet")
    self.out.branch("tightElectrons_id","I",lenVar="nElectron")
    self.out.branch("fakeable_Electrons_id","I",lenVar="nElectron")
    self.out.branch("additional_vetoElectrons_id","I",lenVar="nElectron")
    self.out.branch("tightMuons_id","I",lenVar="nMuon")
    self.out.branch("fakeable_Muons_id","I",lenVar="nMuon")
    self.out.branch("additional_looseMuons_id","I",lenVar="nMuon")
    self.out.branch("Had_tau_id","I",lenVar="nTau")
    self.out.branch("muon_conePt","F",lenVar="nMuon")
    self.out.branch("muon_jet_Ptratio","F",lenVar="nMuon")
    self.out.branch("muon_closest_jetid","I",lenVar="nMuon")
    self.out.branch("electron_conePt","F",lenVar="nElectron")
    self.out.branch("electron_jet_Ptratio","F",lenVar="nElectron")
    self.out.branch("electron_closest_jetid","I",lenVar="nElectron")
    self.is_mc = bool(inputTree.GetBranch("GenJet_pt"))
    self.is_lhe = bool(inputTree.GetBranch("nLHEPart"))
    self.has_cjet_tag = bool(inputTree.GetBranch("Jet_btagDeepFlavCvL"))

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

    lhe_nlepton=0
    if self.is_lhe:
      lheparticle = Collection(event, 'LHEPart')
      for ilhe in range(0, event.nLHEPart):
        if lheparticle[ilhe].status==1 and (abs(lheparticle[ilhe].pdgId)==11 or abs(lheparticle[ilhe].pdgId)==13 or abs(lheparticle[ilhe].pdgId)==15):
          lhe_nlepton=lhe_nlepton+1

    self.out.fillBranch("lhe_nlepton", lhe_nlepton)

    # total number of ele+muon, currently require at least 1 leptons
    if ((event.nMuon + event.nElectron) < 2): return False
    if not event.nJet>0: return False

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
      # following cuts are preseletion for MVA muon ID
      if abs(muons[imu].eta)>2.4 or muons[imu].sip3d>8 or abs(muons[imu].dxy)>0.05 or abs(muons[imu].dz)>0.1 or muons[imu].miniPFRelIso_all>0.4: continue
      if not (muons[imu].looseId and event.Muon_corrected_pt[imu]>10):continue
      if muons[imu].mediumId and muons[imu].tightCharge==2:
        if muons[imu].mvaTTH>-0.2:
          if (event.Muon_corrected_pt[imu]>20):
            muon_v4_temp.SetPtEtaPhiM(event.Muon_corrected_pt[imu], muons[imu].eta, muons[imu].phi, muons[imu].mass)
            muon_v4_temp_raw.SetPtEtaPhiM(muons[imu].pt, muons[imu].eta, muons[imu].phi, muons[imu].mass)
            tightMuons.append(muon_v4_temp.Clone())
            tightMuons_raw.append(muon_v4_temp_raw.Clone())
            tightMuons_pdgid.append(muons[imu].pdgId)
            tightMuons_id.append(imu)
	  else:
            muon_v4_temp.SetPtEtaPhiM(event.Muon_corrected_pt[imu], muons[imu].eta, muons[imu].phi, muons[imu].mass)
            additional_looseMuons.append(muon_v4_temp.Clone())
            additional_looseMuons_pdgid.append(muons[imu].pdgId)
            additional_looseMuons_id.append(imu)
        else:
	  if muon_jet_Ptratio[imu]>0.5 and event.Jet_btagDeepFlavB[muon_closest_jetid[imu]]<flav_cut(muon_conePt[imu]):
            if (event.Muon_corrected_pt[imu]>20):
              muon_v4_temp.SetPtEtaPhiM(event.Muon_corrected_pt[imu], muons[imu].eta, muons[imu].phi, muons[imu].mass)
              fakeable_Muons.append(muon_v4_temp.Clone())
              fakeable_Muons_pdgid.append(muons[imu].pdgId)
              fakeable_Muons_id.append(imu)
            else:
              muon_v4_temp.SetPtEtaPhiM(event.Muon_corrected_pt[imu], muons[imu].eta, muons[imu].phi, muons[imu].mass)
              additional_looseMuons.append(muon_v4_temp.Clone())
              additional_looseMuons_pdgid.append(muons[imu].pdgId)
              additional_looseMuons_id.append(imu)
          else:
            muon_v4_temp.SetPtEtaPhiM(event.Muon_corrected_pt[imu], muons[imu].eta, muons[imu].phi, muons[imu].mass)
            additional_looseMuons.append(muon_v4_temp.Clone())
            additional_looseMuons_pdgid.append(muons[imu].pdgId)
            additional_looseMuons_id.append(imu)
      else:
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
    self.out.fillBranch("muon_conePt", muon_conePt)
    self.out.fillBranch("muon_jet_Ptratio", muon_jet_Ptratio)
    self.out.fillBranch("muon_closest_jetid", muon_closest_jetid)

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
      # following cuts are preseletion for MVA electron ID
      if abs(electrons[iele].eta)>2.5 or electrons[iele].sip3d>8 or abs(electrons[iele].dxy)>0.05 or abs(electrons[iele].dz)>0.1 or electrons[iele].miniPFRelIso_all>0.4: continue
      if not (electrons[iele].mvaFall17V2noIso_WPL and electrons[iele].lostHits<2 and electrons[iele].pt>10 and electrons[iele].convVeto):continue
      if ((abs(electrons[iele].deltaEtaSC+electrons[iele].eta)<1.479 and electrons[iele].sieie<0.011) or (abs(electrons[iele].deltaEtaSC+electrons[iele].eta)>1.479 and electrons[iele].sieie<0.03)) and electrons[iele].hoe<0.1 and electrons[iele].eInvMinusPInv>-0.04:
        if electrons[iele].mvaTTH>0.25 and electrons[iele].tightCharge==2:
          if (electrons[iele].pt>20):
            electron_v4_temp.SetPtEtaPhiM(electrons[iele].pt, electrons[iele].eta, electrons[iele].phi, electrons[iele].mass)
	    electron_v4_temp_raw.SetPtEtaPhiM(electrons[iele].pt/electrons[iele].eCorr, electrons[iele].eta, electrons[iele].phi, electrons[iele].mass/electrons[iele].eCorr)
            tightElectrons.append(electron_v4_temp.Clone())
	    tightElectrons_raw.append(electron_v4_temp_raw.Clone())
            tightElectrons_pdgid.append(electrons[iele].pdgId)
            tightElectrons_id.append(iele)
	  else:
            electron_v4_temp.SetPtEtaPhiM(electrons[iele].pt, electrons[iele].eta, electrons[iele].phi, electrons[iele].mass)
            additional_vetoElectrons.append(electron_v4_temp.Clone())
            additional_vetoElectrons_pdgid.append(electrons[iele].pdgId)
            additional_vetoElectrons_id.append(iele)
        else:
          if electrons[iele].mvaFall17V2noIso_WP90 and electron_jet_Ptratio[iele]>0.6 and event.Jet_btagDeepFlavB[electron_closest_jetid[iele]]<flav_cut(electron_conePt[iele]):
            if (electrons[iele].pt>20):
              electron_v4_temp.SetPtEtaPhiM(electrons[iele].pt, electrons[iele].eta, electrons[iele].phi, electrons[iele].mass)
              fakeable_Electrons.append(electron_v4_temp.Clone())
              fakeable_Electrons_pdgid.append(electrons[iele].pdgId)
              fakeable_Electrons_id.append(iele)
            else:
              electron_v4_temp.SetPtEtaPhiM(electrons[iele].pt, electrons[iele].eta, electrons[iele].phi, electrons[iele].mass)
              additional_vetoElectrons.append(electron_v4_temp.Clone())
              additional_vetoElectrons_pdgid.append(electrons[iele].pdgId)
              additional_vetoElectrons_id.append(iele)
          else:
            electron_v4_temp.SetPtEtaPhiM(electrons[iele].pt, electrons[iele].eta, electrons[iele].phi, electrons[iele].mass)
            additional_vetoElectrons.append(electron_v4_temp.Clone())
            additional_vetoElectrons_pdgid.append(electrons[iele].pdgId)
            additional_vetoElectrons_id.append(iele)
      else:
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
    self.out.fillBranch("electron_conePt", electron_conePt)
    self.out.fillBranch("electron_jet_Ptratio", electron_jet_Ptratio)
    self.out.fillBranch("electron_closest_jetid", electron_closest_jetid)

    # tight leptons and additional loose leptons collection
    tightLeptons = tightMuons + tightElectrons
    tightLeptons.sort(key=lambda x: x.Pt(), reverse=True)
    tightLeptons_raw = tightMuons_raw + tightElectrons_raw
    tightLeptons_raw.sort(key=lambda x: x.Pt(), reverse=True)
    fakeableLeptons = fakeable_Muons + fakeable_Electrons
    fakeableLeptons.sort(key=lambda x: x.Pt(), reverse=True)
    looseLeptons = additional_looseMuons + additional_vetoElectrons
    looseLeptons.sort(key=lambda x: x.Pt(), reverse=True)

    tau_v4_temp=TLorentzVector()
    taus = Collection(event, 'Tau')
    nHad_tau=0
    Had_tau_id=[]
    for itau in range(0, event.nTau):
      tau_v4_temp.SetPtEtaPhiM(taus[itau].pt, taus[itau].eta, taus[itau].phi, taus[itau].mass)
      pass_tau_lep_Dr=1
      if taus[itau].pt>20 and abs(taus[itau].eta)<2.3 and taus[itau].idDecayModeNewDMs and taus[itau].idDeepTau2017v2p1VSe>=4 and taus[itau].idDeepTau2017v2p1VSjet>=4 and taus[itau].idDeepTau2017v2p1VSmu>=1:
	for ilep in range(0,len(tightLeptons)):
          if tau_v4_temp.DeltaR(tightLeptons[ilep])<0.4:pass_tau_lep_Dr=0
	for ilep in range(0,len(fakeableLeptons)):
          if tau_v4_temp.DeltaR(fakeableLeptons[ilep])<0.4:pass_tau_lep_Dr=0
	if pass_tau_lep_Dr:
	  nHad_tau=nHad_tau+1
	  Had_tau_id.append(itau)
    self.out.fillBranch("nHad_tau", nHad_tau)

    # https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL17
    # tightLepVeto PF jets (ak4), 2016 (111=7), 2017/2018 (110=6), medium B-tag WP
    # DeepCSV=(nanoaod btagDeepB) loose: 0.1355, medium: 0.4506, tight: 0.7738
    # DeepFlavor=(nanoaod btagDeepFlavB) loose: 0.0532, medium: 0.3040, tight: 0.7476

    # c-jet tag is based on two-D cuts, medium DeepJet WP:
    # CvsL=btagDeepFlavCvL: 0.085, CvsB=btagDeepFlavCvB: 0.34
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
    mj1j2=-99
    mj1j3=-99
    mj1j4=-99
    mj2j3=-99
    mj2j4=-99
    mj3j4=-99
    mj1j2j3=-99
    mj1j2j4=-99
    mj2j3j4=-99
    mj1j2j3j4=-99
    drj1j2=-99
    drj1j3=-99
    drj1j4=-99
    drj2j3=-99
    drj2j4=-99
    drj3j4=-99
    DeepB_j1_pt=-99
    DeepB_j1_eta=-99
    DeepB_j1_phi=-99
    DeepB_j1_mass=-99
    DeepB_j2_pt=-99
    DeepB_j2_eta=-99
    DeepB_j2_phi=-99
    DeepB_j2_mass=-99
    DeepC_loose_j1_pt=-99
    DeepC_loose_j1_eta=-99
    DeepC_loose_j1_phi=-99
    DeepC_loose_j1_mass=-99
    DeepC_loose_j2_pt=-99
    DeepC_loose_j2_eta=-99
    DeepC_loose_j2_phi=-99
    DeepC_loose_j2_mass=-99
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
    tightJets_c_DeepCSVloose_id = []

    # require DeltaR between Jets and tight leptons greater than 0.4
    jet_v4_all = []
    for ijet in range(0, event.nJet):

      jet_is_tau=0
      if nHad_tau>0:
        for ita in Had_tau_id:
          if ijet==event.Tau_jetIdx[ita]:jet_is_tau=1
      if jet_is_tau:continue

      pass_jet_lep_Dr=1
      jet_v4_temp.SetPtEtaPhiM(event.Jet_pt_nom[ijet],event.Jet_eta[ijet],event.Jet_phi[ijet],event.Jet_mass_nom[ijet])
      for ilep in range(0,len(tightLeptons)):
	if jet_v4_temp.DeltaR(tightLeptons[ilep])<0.4:pass_jet_lep_Dr=0
      for ilep in range(0,len(fakeableLeptons)):
	if jet_v4_temp.DeltaR(fakeableLeptons[ilep])<0.4:pass_jet_lep_Dr=0

      if not (pass_jet_lep_Dr>0):continue
      if not (jets[ijet].jetId==6 and event.Jet_pt_nom[ijet]>30):continue

      if abs(jets[ijet].eta)<4.7 and abs(jets[ijet].eta)>=2.4: 
        tightJets_id_in47.append(ijet)
      if abs(jets[ijet].eta)<2.4:
        tightJets_id_in24.append(ijet)
        jet_v4_all.append(jet_v4_temp.Clone())

      if self.year=="2016" and jets[ijet].btagDeepB > 0.4941:
        tightJets_b_DeepCSVmedium_id.append(ijet)

      if self.year=="2017" and jets[ijet].btagDeepFlavB > 0.3040:
        tightJets_b_DeepCSVmedium_id.append(ijet)

      if self.year=="2018" and jets[ijet].btagDeepFlavB > 0.2783:
        tightJets_b_DeepCSVmedium_id.append(ijet)

    HT=0
    for ijet in range(0,len(tightJets_id_in24)):
      HT=HT+event.Jet_pt_nom[tightJets_id_in24[ijet]]
    self.out.fillBranch("HT",HT)

    for ijet in range(0, event.nJet):
      if not (ijet in tightJets_id_in24):continue
      if (ijet in tightJets_b_DeepCSVmedium_id):continue
      #loose/mudium WP
      if self.has_cjet_tag:
        if self.year=="2017" and jets[ijet].btagDeepFlavCvL>0.085 and jets[ijet].btagDeepFlavCvB<0.34:
  	  tightJets_c_DeepCSVmedium_id.append(ijet)
        if self.year=="2017" and jets[ijet].btagDeepFlavCvL>0.03 and jets[ijet].btagDeepFlavCvB<0.4:
  	  tightJets_c_DeepCSVloose_id.append(ijet)
        if self.year=="2018" and jets[ijet].btagDeepFlavCvL>0.099 and jets[ijet].btagDeepFlavCvB<0.325:
  	  tightJets_c_DeepCSVmedium_id.append(ijet)
        if self.year=="2018" and jets[ijet].btagDeepFlavCvL>0.038 and jets[ijet].btagDeepFlavCvB<0.246:
  	  tightJets_c_DeepCSVloose_id.append(ijet)
      
    n_tight_jet = len(tightJets_id_in24)
    n_bjet_DeepB = len(tightJets_b_DeepCSVmedium_id)
    n_cjet_DeepB_medium = len(tightJets_c_DeepCSVmedium_id)
    n_cjet_DeepB_loose = len(tightJets_c_DeepCSVloose_id)
    self.out.fillBranch("n_tight_jet",n_tight_jet)
    self.out.fillBranch("n_bjet_DeepB",n_bjet_DeepB)
    self.out.fillBranch("n_cjet_DeepB_medium",n_cjet_DeepB_medium)
    self.out.fillBranch("n_cjet_DeepB_loose",n_cjet_DeepB_loose)

    Had_tau_id.extend(np.zeros(event.nTau-len(Had_tau_id),int)-1)
    self.out.fillBranch("Had_tau_id", Had_tau_id)
    

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
      mj1j2=(jet_v4_all[0]+jet_v4_all[1]).M()
      mj1j3=(jet_v4_all[0]+jet_v4_all[2]).M()
      mj1j4=(jet_v4_all[0]+jet_v4_all[3]).M()
      mj2j3=(jet_v4_all[1]+jet_v4_all[2]).M()
      mj2j4=(jet_v4_all[1]+jet_v4_all[3]).M()
      mj3j4=(jet_v4_all[2]+jet_v4_all[3]).M()
      mj1j2j3=(jet_v4_all[0]+jet_v4_all[1]+jet_v4_all[2]).M()
      mj1j2j4=(jet_v4_all[0]+jet_v4_all[1]+jet_v4_all[3]).M()
      mj2j3j4=(jet_v4_all[1]+jet_v4_all[2]+jet_v4_all[3]).M()
      mj1j2j3j4=(jet_v4_all[0]+jet_v4_all[1]+jet_v4_all[2]+jet_v4_all[3]).M()
      drj1j2=jet_v4_all[0].DeltaR(jet_v4_all[1])
      drj1j3=jet_v4_all[0].DeltaR(jet_v4_all[2])
      drj1j4=jet_v4_all[0].DeltaR(jet_v4_all[3])
      drj2j3=jet_v4_all[1].DeltaR(jet_v4_all[2])
      drj2j4=jet_v4_all[1].DeltaR(jet_v4_all[3])
      drj3j4=jet_v4_all[2].DeltaR(jet_v4_all[3])
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
      mj1j2=(jet_v4_all[0]+jet_v4_all[1]).M()
      mj1j3=(jet_v4_all[0]+jet_v4_all[2]).M()
      mj2j3=(jet_v4_all[1]+jet_v4_all[2]).M()
      mj1j2j3=(jet_v4_all[0]+jet_v4_all[1]+jet_v4_all[2]).M()
      drj1j2=jet_v4_all[0].DeltaR(jet_v4_all[1])
      drj1j3=jet_v4_all[0].DeltaR(jet_v4_all[2])
      drj2j3=jet_v4_all[1].DeltaR(jet_v4_all[2])
    if n_tight_jet==2:
      j2_pt=event.Jet_pt_nom[tightJets_id_in24[1]]
      j2_eta=event.Jet_eta[tightJets_id_in24[1]]
      j2_phi=event.Jet_phi[tightJets_id_in24[1]]
      j2_mass=event.Jet_mass_nom[tightJets_id_in24[1]]
      j1_pt=event.Jet_pt_nom[tightJets_id_in24[0]]
      j1_eta=event.Jet_eta[tightJets_id_in24[0]]
      j1_phi=event.Jet_phi[tightJets_id_in24[0]]
      j1_mass=event.Jet_mass_nom[tightJets_id_in24[0]]
      mj1j2=(jet_v4_all[0]+jet_v4_all[1]).M()
      drj1j2=jet_v4_all[0].DeltaR(jet_v4_all[1])
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

    if n_cjet_DeepB_loose>1:
      DeepC_loose_j1_pt=event.Jet_pt_nom[tightJets_c_DeepCSVloose_id[0]]
      DeepC_loose_j1_eta=event.Jet_eta[tightJets_c_DeepCSVloose_id[0]]
      DeepC_loose_j1_phi=event.Jet_phi[tightJets_c_DeepCSVloose_id[0]]
      DeepC_loose_j1_mass=event.Jet_mass_nom[tightJets_c_DeepCSVloose_id[0]]
      DeepC_loose_j2_pt=event.Jet_pt_nom[tightJets_c_DeepCSVloose_id[1]]
      DeepC_loose_j2_eta=event.Jet_eta[tightJets_c_DeepCSVloose_id[1]]
      DeepC_loose_j2_phi=event.Jet_phi[tightJets_c_DeepCSVloose_id[1]]
      DeepC_loose_j2_mass=event.Jet_mass_nom[tightJets_c_DeepCSVloose_id[1]]
    if n_cjet_DeepB_loose==1:
      DeepC_loose_j1_pt=event.Jet_pt_nom[tightJets_c_DeepCSVloose_id[0]]
      DeepC_loose_j1_eta=event.Jet_eta[tightJets_c_DeepCSVloose_id[0]]
      DeepC_loose_j1_phi=event.Jet_phi[tightJets_c_DeepCSVloose_id[0]]
      DeepC_loose_j1_mass=event.Jet_mass_nom[tightJets_c_DeepCSVloose_id[0]]

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
    self.out.fillBranch("mj1j2", mj1j2)
    self.out.fillBranch("mj1j3", mj1j3)
    self.out.fillBranch("mj1j4", mj1j4)
    self.out.fillBranch("mj2j3", mj2j3)
    self.out.fillBranch("mj2j4", mj2j4)
    self.out.fillBranch("mj3j4", mj3j4)
    self.out.fillBranch("mj1j2j3", mj1j2j3)
    self.out.fillBranch("mj1j2j4", mj1j2j4)
    self.out.fillBranch("mj2j3j4", mj2j3j4)
    self.out.fillBranch("mj1j2j3j4", mj1j2j3j4)
    self.out.fillBranch("drj1j2", drj1j2)
    self.out.fillBranch("drj1j3", drj1j3)
    self.out.fillBranch("drj1j4", drj1j4)
    self.out.fillBranch("drj2j3", drj2j3)
    self.out.fillBranch("drj2j4", drj2j4)
    self.out.fillBranch("drj3j4", drj3j4)
    self.out.fillBranch("DeepB_j1_pt",DeepB_j1_pt)
    self.out.fillBranch("DeepB_j1_eta",DeepB_j1_eta)
    self.out.fillBranch("DeepB_j1_phi",DeepB_j1_phi)
    self.out.fillBranch("DeepB_j1_mass",DeepB_j1_mass)
    self.out.fillBranch("DeepB_j2_pt",DeepB_j2_pt)
    self.out.fillBranch("DeepB_j2_eta",DeepB_j2_eta)
    self.out.fillBranch("DeepB_j2_phi",DeepB_j2_phi)
    self.out.fillBranch("DeepB_j2_mass",DeepB_j2_mass)
    self.out.fillBranch("DeepC_loose_j1_pt",DeepC_loose_j1_pt)
    self.out.fillBranch("DeepC_loose_j1_eta",DeepC_loose_j1_eta)
    self.out.fillBranch("DeepC_loose_j1_phi",DeepC_loose_j1_phi)
    self.out.fillBranch("DeepC_loose_j1_mass",DeepC_loose_j1_mass)
    self.out.fillBranch("DeepC_loose_j2_pt",DeepC_loose_j2_pt)
    self.out.fillBranch("DeepC_loose_j2_eta",DeepC_loose_j2_eta)
    self.out.fillBranch("DeepC_loose_j2_phi",DeepC_loose_j2_phi)
    self.out.fillBranch("DeepC_loose_j2_mass",DeepC_loose_j2_mass)
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
    tightJets_c_DeepCSVloose_id.extend(np.zeros(event.nJet-len(tightJets_c_DeepCSVloose_id),int)-1)

    self.out.fillBranch("tightJets_id_in24",tightJets_id_in24)
    self.out.fillBranch("tightJets_id_in47",tightJets_id_in47)
    self.out.fillBranch("tightJets_b_DeepCSVmedium_id",tightJets_b_DeepCSVmedium_id)
    self.out.fillBranch("tightJets_c_DeepCSVmedium_id",tightJets_c_DeepCSVmedium_id)
    self.out.fillBranch("tightJets_c_DeepCSVloose_id",tightJets_c_DeepCSVloose_id)

    if len(tightLeptons)+len(fakeableLeptons)<2:return False

    #    t     t      cccc
    #  ttttt ttttt   cc
    #    t     t    cc
    #    t     t     cc
    #    ttt   ttt    cccc
    #region: only 2 tight leptons, at least three jets, (b/c-tag are not appplied here) 
    #ttc region lepton number selections
    ttc_nl=False
    #ttc region jet and bjet selection
    ttc_jets=False
    ttc_2P0F=False
    ttc_1P1F=False
    ttc_0P2F=False
    ttc_lep1_faketag=False
    ttc_lep2_faketag=False
    #ttc region tag, 1:2 muon, 2:1 muon, 3:0 muon
    ttc_region=0
    ttc_l1_id=-1
    ttc_l2_id=-1
    ttc_l1_isprompt=0
    ttc_l2_isprompt=0
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
    if len(tightLeptons)+len(fakeableLeptons)==2 and len(looseLeptons)==0:
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

      if len(tightLeptons)==2:
	ttc_2P0F=True
        ttc_1P1F=False
        ttc_0P2F=False
        ttc_lep1_faketag=False
        ttc_lep2_faketag=False
        if len(tightElectrons)==0 and abs(tightMuons_pdgid[0]+tightMuons_pdgid[1])==26:
          ttc_region=1
          ttc_l1_id=tightMuons_id[0]
          ttc_l2_id=tightMuons_id[1]
          ttc_l1_pdgid=tightMuons_pdgid[0]
          ttc_l2_pdgid=tightMuons_pdgid[1]
          ttc_l1_pt=tightMuons[0].Pt()
          ttc_l1_eta=tightMuons[0].Eta()
          ttc_l1_phi=tightMuons[0].Phi()
          ttc_l1_mass=tightMuons[0].M()
          ttc_l2_pt=tightMuons[1].Pt()
          ttc_l2_eta=tightMuons[1].Eta()
          ttc_l2_phi=tightMuons[1].Phi()
          ttc_l2_mass=tightMuons[1].M()
          ttc_mll=(tightMuons[0]+tightMuons[1]).M()
          ttc_drll=tightMuons[0].DeltaR(tightMuons[1])
          ttc_dphill=tightMuons[0].DeltaPhi(tightMuons[1])
	  if self.is_mc:
	    if event.Muon_genPartFlav[ttc_l1_id]==1 or event.Muon_genPartFlav[ttc_l1_id]==15:ttc_l1_isprompt=1
	    if event.Muon_genPartFlav[ttc_l2_id]==1 or event.Muon_genPartFlav[ttc_l2_id]==15:ttc_l2_isprompt=1
        if len(tightElectrons)==1 and abs(tightMuons_pdgid[0]+tightElectrons_pdgid[0])==24:
          ttc_region=2
          ttc_l1_id=tightMuons_id[0]
          ttc_l1_pdgid=tightMuons_pdgid[0]
          ttc_l2_id=tightElectrons_id[0]
          ttc_l2_pdgid=tightElectrons_pdgid[0]
          ttc_l1_pt=tightMuons[0].Pt()
          ttc_l1_eta=tightMuons[0].Eta()
          ttc_l1_phi=tightMuons[0].Phi()
          ttc_l1_mass=tightMuons[0].M()
          ttc_l2_pt=tightElectrons[0].Pt()
          ttc_l2_eta=tightElectrons[0].Eta()
          ttc_l2_phi=tightElectrons[0].Phi()
          ttc_l2_mass=tightElectrons[0].M()
          ttc_mll=(tightMuons[0]+tightElectrons[0]).M()
          ttc_drll=tightMuons[0].DeltaR(tightElectrons[0])
          ttc_dphill=tightMuons[0].DeltaPhi(tightElectrons[0])
	  if self.is_mc:
	    if event.Muon_genPartFlav[ttc_l1_id]==1 or event.Muon_genPartFlav[ttc_l1_id]==15:ttc_l1_isprompt=1
	    if event.Electron_genPartFlav[ttc_l2_id]==1 or event.Electron_genPartFlav[ttc_l2_id]==15:ttc_l2_isprompt=1
        if len(tightElectrons)==2 and abs(tightElectrons_pdgid[0]+tightElectrons_pdgid[1])==22:
          ttc_region=3
          ttc_l1_id=tightElectrons_id[0]
          ttc_l2_id=tightElectrons_id[1]
          ttc_l1_pdgid=tightElectrons_pdgid[0]
          ttc_l2_pdgid=tightElectrons_pdgid[1]
          ttc_l1_pt=tightElectrons[0].Pt()
          ttc_l1_eta=tightElectrons[0].Eta()
          ttc_l1_phi=tightElectrons[0].Phi()
          ttc_l1_mass=tightElectrons[0].M()
          ttc_l2_pt=tightElectrons[1].Pt()
          ttc_l2_eta=tightElectrons[1].Eta()
          ttc_l2_phi=tightElectrons[1].Phi()
          ttc_l2_mass=tightElectrons[1].M()
          ttc_mll=(tightElectrons[0]+tightElectrons[1]).M()
          ttc_drll=tightElectrons[0].DeltaR(tightElectrons[1])
          ttc_dphill=tightElectrons[0].DeltaPhi(tightElectrons[1])
	  if self.is_mc:
	    if event.Electron_genPartFlav[ttc_l1_id]==1 or event.Electron_genPartFlav[ttc_l1_id]==15:ttc_l1_isprompt=1
	    if event.Electron_genPartFlav[ttc_l2_id]==1 or event.Electron_genPartFlav[ttc_l2_id]==15:ttc_l2_isprompt=1

      # one real and one fake lepton
      if len(tightLeptons)==1:
	ttc_2P0F=False
        ttc_1P1F=True
        ttc_0P2F=False
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
	      if self.is_mc:
	        if event.Muon_genPartFlav[ttc_l1_id]==1 or event.Muon_genPartFlav[ttc_l1_id]==15:ttc_l1_isprompt=1
	        if event.Muon_genPartFlav[ttc_l2_id]==1 or event.Muon_genPartFlav[ttc_l2_id]==15:ttc_l2_isprompt=1
	    else:
              ttc_lep1_faketag=True
              ttc_lep2_faketag=False
              ttc_l1_id=fakeable_Muons_id[0]
              ttc_l2_id=tightMuons_id[0]
              ttc_l1_pdgid=fakeable_Muons_pdgid[0]
              ttc_l2_pdgid=tightMuons_pdgid[0]
              ttc_l1_pt=fakeable_Muons[0].Pt()
              ttc_l1_eta=fakeable_Muons[0].Eta()
              ttc_l1_phi=fakeable_Muons[0].Phi()
              ttc_l1_mass=fakeable_Muons[0].M()
              ttc_l2_pt=tightMuons[0].Pt()
              ttc_l2_eta=tightMuons[0].Eta()
              ttc_l2_phi=tightMuons[0].Phi()
              ttc_l2_mass=tightMuons[0].M()
	      if self.is_mc:
	        if event.Muon_genPartFlav[ttc_l1_id]==1 or event.Muon_genPartFlav[ttc_l1_id]==15:ttc_l1_isprompt=1
	        if event.Muon_genPartFlav[ttc_l2_id]==1 or event.Muon_genPartFlav[ttc_l2_id]==15:ttc_l2_isprompt=1
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
	    if self.is_mc:
	      if event.Muon_genPartFlav[ttc_l1_id]==1 or event.Muon_genPartFlav[ttc_l1_id]==15:ttc_l1_isprompt=1
	      if event.Electron_genPartFlav[ttc_l2_id]==1 or event.Electron_genPartFlav[ttc_l2_id]==15:ttc_l2_isprompt=1
	
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
	    if self.is_mc:
	      if event.Muon_genPartFlav[ttc_l1_id]==1 or event.Muon_genPartFlav[ttc_l1_id]==15:ttc_l1_isprompt=1
	      if event.Electron_genPartFlav[ttc_l2_id]==1 or event.Electron_genPartFlav[ttc_l2_id]==15:ttc_l2_isprompt=1

	  # one fake electron
          if len(fakeable_Electrons)==1 and abs(tightElectrons_pdgid[0]+fakeable_Electrons_pdgid[0])==22:
            ttc_region=3
            if tightElectrons[0].Pt()>fakeable_Electrons[0].Pt():
              ttc_lep1_faketag=False
              ttc_lep2_faketag=True
              ttc_l1_id=tightElectrons_id[0]
              ttc_l2_id=fakeable_Electrons_id[0]
              ttc_l1_pdgid=tightElectrons_pdgid[0]
              ttc_l2_pdgid=fakeable_Electrons_pdgid[0]
              ttc_l1_pt=tightElectrons[0].Pt()
              ttc_l1_eta=tightElectrons[0].Eta()
              ttc_l1_phi=tightElectrons[0].Phi()
              ttc_l1_mass=tightElectrons[0].M()
              ttc_l2_pt=fakeable_Electrons[0].Pt()
              ttc_l2_eta=fakeable_Electrons[0].Eta()
              ttc_l2_phi=fakeable_Electrons[0].Phi()
              ttc_l2_mass=fakeable_Electrons[0].M()
	      if self.is_mc:
	        if event.Electron_genPartFlav[ttc_l1_id]==1 or event.Electron_genPartFlav[ttc_l1_id]==15:ttc_l1_isprompt=1
	        if event.Electron_genPartFlav[ttc_l2_id]==1 or event.Electron_genPartFlav[ttc_l2_id]==15:ttc_l2_isprompt=1
	    else:
              ttc_lep1_faketag=True
              ttc_lep2_faketag=False
              ttc_l1_id=fakeable_Electrons_id[0]
              ttc_l2_id=tightElectrons_id[0]
              ttc_l1_pdgid=fakeable_Electrons_pdgid[0]
              ttc_l2_pdgid=tightElectrons_pdgid[0]
              ttc_l1_pt=fakeable_Electrons[0].Pt()
              ttc_l1_eta=fakeable_Electrons[0].Eta()
              ttc_l1_phi=fakeable_Electrons[0].Phi()
              ttc_l1_mass=fakeable_Electrons[0].M()
              ttc_l2_pt=tightElectrons[0].Pt()
              ttc_l2_eta=tightElectrons[0].Eta()
              ttc_l2_phi=tightElectrons[0].Phi()
              ttc_l2_mass=tightElectrons[0].M()
	      if self.is_mc:
	        if event.Electron_genPartFlav[ttc_l1_id]==1 or event.Electron_genPartFlav[ttc_l1_id]==15:ttc_l1_isprompt=1
	        if event.Electron_genPartFlav[ttc_l2_id]==1 or event.Electron_genPartFlav[ttc_l2_id]==15:ttc_l2_isprompt=1
            ttc_mll=(fakeable_Electrons[0]+tightElectrons[0]).M()
            ttc_drll=fakeable_Electrons[0].DeltaR(tightElectrons[0])
            ttc_dphill=fakeable_Electrons[0].DeltaPhi(tightElectrons[0])

      # two fake leptons
      if len(tightLeptons)==0:
	ttc_2P0F=False
        ttc_1P1F=False
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
	  if self.is_mc:
	    if event.Muon_genPartFlav[ttc_l1_id]==1 or event.Muon_genPartFlav[ttc_l1_id]==15:ttc_l1_isprompt=1
	    if event.Muon_genPartFlav[ttc_l2_id]==1 or event.Muon_genPartFlav[ttc_l2_id]==15:ttc_l2_isprompt=1

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
	  if self.is_mc:
	    if event.Muon_genPartFlav[ttc_l1_id]==1 or event.Muon_genPartFlav[ttc_l1_id]==15:ttc_l1_isprompt=1
	    if event.Electron_genPartFlav[ttc_l2_id]==1 or event.Electron_genPartFlav[ttc_l2_id]==15:ttc_l2_isprompt=1

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
	  if self.is_mc:
	    if event.Electron_genPartFlav[ttc_l1_id]==1 or event.Electron_genPartFlav[ttc_l1_id]==15:ttc_l1_isprompt=1
	    if event.Electron_genPartFlav[ttc_l2_id]==1 or event.Electron_genPartFlav[ttc_l2_id]==15:ttc_l2_isprompt=1

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

    self.out.fillBranch("ttc_2P0F", ttc_2P0F)
    self.out.fillBranch("ttc_1P1F", ttc_1P1F)
    self.out.fillBranch("ttc_0P2F", ttc_0P2F)
    self.out.fillBranch("ttc_lep1_faketag", ttc_lep1_faketag)
    self.out.fillBranch("ttc_lep2_faketag", ttc_lep2_faketag)
    self.out.fillBranch("ttc_nl", ttc_nl)
    self.out.fillBranch("ttc_jets", ttc_jets)
    self.out.fillBranch("ttc_region", ttc_region)
    self.out.fillBranch("ttc_l1_id", ttc_l1_id)
    self.out.fillBranch("ttc_l2_id", ttc_l2_id)
    self.out.fillBranch("ttc_l1_isprompt", ttc_l1_isprompt)
    self.out.fillBranch("ttc_l2_isprompt", ttc_l2_isprompt)
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

    #Opposite sign region lepton number selections
    OPS_nl=False
    #OPS region tag, 0: fail to pass the OPS selection, 1:2 muon, 2:1 muon, 3:0 muon
    OPS_2P0F=False
    OPS_1P1F=False
    OPS_0P2F=False
    OPS_lep1_faketag=False
    OPS_lep2_faketag=False
    OPS_region=0
    OPS_l1_id=-1
    OPS_l2_id=-1
    OPS_l1_isprompt=0
    OPS_l2_isprompt=0
    OPS_l1_pdgid=-99
    OPS_l2_pdgid=-99
    OPS_l1_pt=-99
    OPS_l1_pt_raw=-99
    OPS_l1_eta=-99
    OPS_l1_phi=-99
    OPS_l1_mass=-99
    OPS_l2_pt=-99
    OPS_l2_pt_raw=-99
    OPS_l2_eta=-99
    OPS_l2_phi=-99
    OPS_l2_mass=-99
    OPS_z_mass=-99
    OPS_z_mass_raw=-99
    OPS_z_pt=-99
    OPS_z_pt_raw=-99
    OPS_z_eta=-99
    OPS_z_eta_raw=-99
    OPS_z_phi=-99
    OPS_z_phi_raw=-99
    OPS_drll=-99

    # the two leptons with pt >20, 3th lepton veto
    if len(tightLeptons)+len(fakeableLeptons)==2 and len(looseLeptons)==0:
      OPS_nl=True
    if OPS_nl:
      if len(tightLeptons)==2:
	OPS_2P0F=True
        OPS_1P1F=False
        OPS_0P2F=False
        OPS_lep1_faketag=False
        OPS_lep2_faketag=False
        # 2 muons case
        if len(tightElectrons)==0 and abs(tightMuons_pdgid[0]+tightMuons_pdgid[1])==0:
	  OPS_region=1
	  OPS_l1_id=tightMuons_id[0]
	  OPS_l2_id=tightMuons_id[1]
	  OPS_l1_pdgid=tightMuons_pdgid[0]
	  OPS_l2_pdgid=tightMuons_pdgid[1]
	  OPS_l1_pt=tightMuons[0].Pt()
	  OPS_l1_pt_raw=tightMuons_raw[0].Pt()
	  OPS_l1_eta=tightMuons[0].Eta()
	  OPS_l1_phi=tightMuons[0].Phi()
	  OPS_l1_mass=tightMuons[0].M()
	  OPS_l2_pt=tightMuons[1].Pt()
	  OPS_l2_pt_raw=tightMuons_raw[1].Pt()
	  OPS_l2_eta=tightMuons[1].Eta()
	  OPS_l2_phi=tightMuons[1].Phi()
	  OPS_l2_mass=tightMuons[1].M()
	  OPS_z_mass=(tightLeptons[0]+tightLeptons[1]).M()
	  OPS_z_pt=(tightLeptons[0]+tightLeptons[1]).Pt()
	  OPS_z_eta=(tightLeptons[0]+tightLeptons[1]).Eta()
	  OPS_z_phi=(tightLeptons[0]+tightLeptons[1]).Phi()
	  OPS_z_mass_raw=(tightLeptons_raw[0]+tightLeptons_raw[1]).M()
	  OPS_z_pt_raw=(tightLeptons_raw[0]+tightLeptons_raw[1]).Pt()
	  OPS_z_eta_raw=(tightLeptons_raw[0]+tightLeptons_raw[1]).Eta()
	  OPS_z_phi_raw=(tightLeptons_raw[0]+tightLeptons_raw[1]).Phi()
	  OPS_drll=tightLeptons[0].DeltaR(tightLeptons[1])
	  if self.is_mc:
	    if event.Muon_genPartFlav[OPS_l1_id]==1 or event.Muon_genPartFlav[OPS_l1_id]==15:OPS_l1_isprompt=1
	    if event.Muon_genPartFlav[OPS_l2_id]==1 or event.Muon_genPartFlav[OPS_l2_id]==15:OPS_l2_isprompt=1
        # 2 eles case
        if len(tightElectrons)==2 and abs(tightElectrons_pdgid[0]+tightElectrons_pdgid[1])==0:
    	  OPS_region=3
          OPS_l1_id=tightElectrons_id[0]
          OPS_l2_id=tightElectrons_id[1]
          OPS_l1_pdgid=tightElectrons_pdgid[0]
          OPS_l2_pdgid=tightElectrons_pdgid[1]
  	  OPS_l1_pt=tightElectrons[0].Pt()
  	  OPS_l1_pt_raw=tightElectrons_raw[0].Pt()
  	  OPS_l1_eta=tightElectrons[0].Eta()
  	  OPS_l1_phi=tightElectrons[0].Phi()
  	  OPS_l1_mass=tightElectrons[0].M()
  	  OPS_l2_pt=tightElectrons[1].Pt()
  	  OPS_l2_pt_raw=tightElectrons_raw[1].Pt()
  	  OPS_l2_eta=tightElectrons[1].Eta()
  	  OPS_l2_phi=tightElectrons[1].Phi()
  	  OPS_l2_mass=tightElectrons[1].M()
  	  OPS_z_mass=(tightLeptons[0]+tightLeptons[1]).M()
  	  OPS_z_pt=(tightLeptons[0]+tightLeptons[1]).Pt()
  	  OPS_z_eta=(tightLeptons[0]+tightLeptons[1]).Eta()
  	  OPS_z_phi=(tightLeptons[0]+tightLeptons[1]).Phi()
  	  OPS_z_mass_raw=(tightLeptons_raw[0]+tightLeptons_raw[1]).M()
  	  OPS_z_pt_raw=(tightLeptons_raw[0]+tightLeptons_raw[1]).Pt()
  	  OPS_z_eta_raw=(tightLeptons_raw[0]+tightLeptons_raw[1]).Eta()
  	  OPS_z_phi_raw=(tightLeptons_raw[0]+tightLeptons_raw[1]).Phi()
  	  OPS_drll=tightLeptons[0].DeltaR(tightLeptons[1])
	  if self.is_mc:
	    if event.Electron_genPartFlav[OPS_l1_id]==1 or event.Electron_genPartFlav[OPS_l1_id]==15:OPS_l1_isprompt=1
	    if event.Electron_genPartFlav[OPS_l2_id]==1 or event.Electron_genPartFlav[OPS_l2_id]==15:OPS_l2_isprompt=1
        # 1 ele case
        if len(tightElectrons)==1 and (sign(tightMuons_pdgid[0])+sign(tightElectrons_pdgid[0]))==0:
	  OPS_region=2
          OPS_l1_id=tightMuons_id[0]
          OPS_l2_id=tightElectrons_id[0]
          OPS_l1_pdgid=tightMuons_pdgid[0]
          OPS_l2_pdgid=tightElectrons_pdgid[0]
	  OPS_l1_pt=tightMuons[0].Pt()
	  OPS_l1_pt_raw=tightMuons_raw[0].Pt()
	  OPS_l1_eta=tightMuons[0].Eta()
	  OPS_l1_phi=tightMuons[0].Phi()
	  OPS_l1_mass=tightMuons[0].M()
	  OPS_l2_pt=tightElectrons[0].Pt()
	  OPS_l2_pt_raw=tightElectrons_raw[0].Pt()
	  OPS_l2_eta=tightElectrons[0].Eta()
	  OPS_l2_phi=tightElectrons[0].Phi()
	  OPS_l2_mass=tightElectrons[0].M()
	  OPS_z_mass=(tightLeptons[0]+tightLeptons[1]).M()
	  OPS_z_pt=(tightLeptons[0]+tightLeptons[1]).Pt()
	  OPS_z_eta=(tightLeptons[0]+tightLeptons[1]).Eta()
	  OPS_z_phi=(tightLeptons[0]+tightLeptons[1]).Phi()
  	  OPS_z_mass_raw=(tightLeptons_raw[0]+tightLeptons_raw[1]).M()
  	  OPS_z_pt_raw=(tightLeptons_raw[0]+tightLeptons_raw[1]).Pt()
  	  OPS_z_eta_raw=(tightLeptons_raw[0]+tightLeptons_raw[1]).Eta()
  	  OPS_z_phi_raw=(tightLeptons_raw[0]+tightLeptons_raw[1]).Phi()
	  OPS_drll=tightLeptons[0].DeltaR(tightLeptons[1])
	  if self.is_mc:
	    if event.Muon_genPartFlav[OPS_l1_id]==1 or event.Muon_genPartFlav[OPS_l1_id]==15:OPS_l1_isprompt=1
	    if event.Electron_genPartFlav[OPS_l2_id]==1 or event.Electron_genPartFlav[OPS_l2_id]==15:OPS_l2_isprompt=1

      if len(tightLeptons)==1:
	OPS_2P0F=False
        OPS_1P1F=True
        OPS_0P2F=False
	# one real muons
        if len(tightElectrons)==0:
          # one fake muon
          if len(fakeable_Electrons)==0 and abs(tightMuons_pdgid[0]+fakeable_Muons_pdgid[0])==0:
            OPS_region=1
            if tightMuons[0].Pt()>fakeable_Muons[0].Pt():
              OPS_lep1_faketag=False
              OPS_lep2_faketag=True
              OPS_l1_id=tightMuons_id[0]
              OPS_l2_id=fakeable_Muons_id[0]
              OPS_l1_pdgid=tightMuons_pdgid[0]
              OPS_l2_pdgid=fakeable_Muons_pdgid[0]
              OPS_l1_pt=tightMuons[0].Pt()
              OPS_l1_eta=tightMuons[0].Eta()
              OPS_l1_phi=tightMuons[0].Phi()
              OPS_l1_mass=tightMuons[0].M()
              OPS_l2_pt=fakeable_Muons[0].Pt()
              OPS_l2_eta=fakeable_Muons[0].Eta()
              OPS_l2_phi=fakeable_Muons[0].Phi()
              OPS_l2_mass=fakeable_Muons[0].M()
	      if self.is_mc:
	        if event.Muon_genPartFlav[OPS_l1_id]==1 or event.Muon_genPartFlav[OPS_l1_id]==15:OPS_l1_isprompt=1
	        if event.Muon_genPartFlav[OPS_l2_id]==1 or event.Muon_genPartFlav[OPS_l2_id]==15:OPS_l2_isprompt=1
	    else:
              OPS_lep1_faketag=True
              OPS_lep2_faketag=False
              OPS_l1_id=fakeable_Muons_id[0]
              OPS_l2_id=tightMuons_id[0]
              OPS_l1_pdgid=fakeable_Muons_pdgid[0]
              OPS_l2_pdgid=tightMuons_pdgid[0]
              OPS_l1_pt=fakeable_Muons[0].Pt()
              OPS_l1_eta=fakeable_Muons[0].Eta()
              OPS_l1_phi=fakeable_Muons[0].Phi()
              OPS_l1_mass=fakeable_Muons[0].M()
              OPS_l2_pt=tightMuons[0].Pt()
              OPS_l2_eta=tightMuons[0].Eta()
              OPS_l2_phi=tightMuons[0].Phi()
              OPS_l2_mass=tightMuons[0].M()
	      if self.is_mc:
	        if event.Muon_genPartFlav[OPS_l1_id]==1 or event.Muon_genPartFlav[OPS_l1_id]==15:OPS_l1_isprompt=1
	        if event.Muon_genPartFlav[OPS_l2_id]==1 or event.Muon_genPartFlav[OPS_l2_id]==15:OPS_l2_isprompt=1
            OPS_z_mass=(fakeable_Muons[0]+tightMuons[0]).M()
            OPS_z_pt=(fakeable_Muons[0]+tightMuons[0]).Pt()
            OPS_z_eta=(fakeable_Muons[0]+tightMuons[0]).Eta()
            OPS_z_phi=(fakeable_Muons[0]+tightMuons[0]).Phi()
            OPS_drll=fakeable_Muons[0].DeltaR(tightMuons[0])
	  
	  # one fake electron
          if len(fakeable_Electrons)==1 and (sign(tightMuons_pdgid[0])+sign(fakeable_Electrons_pdgid[0]))==0:
            OPS_region=2
            OPS_lep1_faketag=False
            OPS_lep2_faketag=True
            OPS_l1_id=tightMuons_id[0]
            OPS_l2_id=fakeable_Electrons_id[0]
            OPS_l1_pdgid=tightMuons_pdgid[0]
            OPS_l2_pdgid=fakeable_Electrons_pdgid[0]
            OPS_l1_pt=tightMuons[0].Pt()
            OPS_l1_eta=tightMuons[0].Eta()
            OPS_l1_phi=tightMuons[0].Phi()
            OPS_l1_mass=tightMuons[0].M()
            OPS_l2_pt=fakeable_Electrons[0].Pt()
            OPS_l2_eta=fakeable_Electrons[0].Eta()
            OPS_l2_phi=fakeable_Electrons[0].Phi()
            OPS_l2_mass=fakeable_Electrons[0].M()
            OPS_z_mass=(fakeable_Electrons[0]+tightMuons[0]).M()
            OPS_z_pt=(fakeable_Electrons[0]+tightMuons[0]).Pt()
            OPS_z_eta=(fakeable_Electrons[0]+tightMuons[0]).Eta()
            OPS_z_phi=(fakeable_Electrons[0]+tightMuons[0]).Phi()
            OPS_drll=fakeable_Electrons[0].DeltaR(tightMuons[0])
	    if self.is_mc:
	      if event.Muon_genPartFlav[OPS_l1_id]==1 or event.Muon_genPartFlav[OPS_l1_id]==15:OPS_l1_isprompt=1
	      if event.Electron_genPartFlav[OPS_l2_id]==1 or event.Electron_genPartFlav[OPS_l2_id]==15:OPS_l2_isprompt=1

	# one real electrons
        if len(tightElectrons)==1:
          # one fake muon
          if len(fakeable_Electrons)==0 and abs(sign(tightElectrons_pdgid[0])+sign(fakeable_Muons_pdgid[0]))==0:
            OPS_region=2
            OPS_lep1_faketag=True
            OPS_lep2_faketag=False
            OPS_l1_id=fakeable_Muons_id[0]
            OPS_l2_id=tightElectrons_id[0]
            OPS_l1_pdgid=fakeable_Muons_pdgid[0]
            OPS_l2_pdgid=tightElectrons_pdgid[0]
            OPS_l1_pt=fakeable_Muons[0].Pt()
            OPS_l1_eta=fakeable_Muons[0].Eta()
            OPS_l1_phi=fakeable_Muons[0].Phi()
            OPS_l1_mass=fakeable_Muons[0].M()
            OPS_l2_pt=tightElectrons[0].Pt()
            OPS_l2_eta=tightElectrons[0].Eta()
            OPS_l2_phi=tightElectrons[0].Phi()
            OPS_l2_mass=tightElectrons[0].M()
            OPS_z_mass=(fakeable_Muons[0]+tightElectrons[0]).M()
            OPS_z_pt=(fakeable_Muons[0]+tightElectrons[0]).Pt()
            OPS_z_eta=(fakeable_Muons[0]+tightElectrons[0]).Eta()
            OPS_z_phi=(fakeable_Muons[0]+tightElectrons[0]).Phi()
            OPS_drll=fakeable_Muons[0].DeltaR(tightElectrons[0])    
	    if self.is_mc:
	      if event.Muon_genPartFlav[OPS_l1_id]==1 or event.Muon_genPartFlav[OPS_l1_id]==15:OPS_l1_isprompt=1
	      if event.Electron_genPartFlav[OPS_l2_id]==1 or event.Electron_genPartFlav[OPS_l2_id]==15:OPS_l2_isprompt=1
	  # one fake electron
          if len(fakeable_Electrons)==1 and abs(tightElectrons_pdgid[0]+fakeable_Electrons_pdgid[0])==0:
            OPS_region=3
            if tightElectrons[0].Pt()>fakeable_Electrons[0].Pt():
              OPS_lep1_faketag=False
              OPS_lep2_faketag=True
              OPS_l1_id=tightElectrons_id[0]
              OPS_l2_id=fakeable_Electrons_id[0]
              OPS_l1_pdgid=tightElectrons_pdgid[0]
              OPS_l2_pdgid=fakeable_Electrons_pdgid[0]
              OPS_l1_pt=tightElectrons[0].Pt()
              OPS_l1_eta=tightElectrons[0].Eta()
              OPS_l1_phi=tightElectrons[0].Phi()
              OPS_l1_mass=tightElectrons[0].M()
              OPS_l2_pt=fakeable_Electrons[0].Pt()
              OPS_l2_eta=fakeable_Electrons[0].Eta()
              OPS_l2_phi=fakeable_Electrons[0].Phi()
              OPS_l2_mass=fakeable_Electrons[0].M()
	      if self.is_mc:
	        if event.Electron_genPartFlav[OPS_l1_id]==1 or event.Electron_genPartFlav[OPS_l1_id]==15:OPS_l1_isprompt=1
	        if event.Electron_genPartFlav[OPS_l2_id]==1 or event.Electron_genPartFlav[OPS_l2_id]==15:OPS_l2_isprompt=1
            else:
	      OPS_lep1_faketag=True
              OPS_lep2_faketag=False
              OPS_l2_id=tightElectrons_id[0]
              OPS_l1_id=fakeable_Electrons_id[0]
              OPS_l2_pdgid=tightElectrons_pdgid[0]
              OPS_l1_pdgid=fakeable_Electrons_pdgid[0]
              OPS_l2_pt=tightElectrons[0].Pt()
              OPS_l2_eta=tightElectrons[0].Eta()
              OPS_l2_phi=tightElectrons[0].Phi()
              OPS_l2_mass=tightElectrons[0].M()
              OPS_l1_pt=fakeable_Electrons[0].Pt()
              OPS_l1_eta=fakeable_Electrons[0].Eta()
              OPS_l1_phi=fakeable_Electrons[0].Phi()
              OPS_l1_mass=fakeable_Electrons[0].M()
	      if self.is_mc:
	        if event.Electron_genPartFlav[OPS_l1_id]==1 or event.Electron_genPartFlav[OPS_l1_id]==15:OPS_l1_isprompt=1
	        if event.Electron_genPartFlav[OPS_l2_id]==1 or event.Electron_genPartFlav[OPS_l2_id]==15:OPS_l2_isprompt=1
            OPS_z_mass=(fakeable_Electrons[0]+tightElectrons[0]).M()
            OPS_z_pt=(fakeable_Electrons[0]+tightElectrons[0]).Pt()
            OPS_z_eta=(fakeable_Electrons[0]+tightElectrons[0]).Eta()
            OPS_z_phi=(fakeable_Electrons[0]+tightElectrons[0]).Phi()
            OPS_drll=fakeable_Electrons[0].DeltaR(tightElectrons[0])

      if len(tightLeptons)==0:
	OPS_2P0F=False
        OPS_1P1F=False
        OPS_0P2F=True
	# 2 muons case
        if len(fakeable_Electrons)==0 and abs(fakeable_Muons_pdgid[0]+fakeable_Muons_pdgid[1])==0:
          OPS_region=1
          OPS_l1_id=fakeable_Muons_id[0]
          OPS_l2_id=fakeable_Muons_id[1]
          OPS_l1_pdgid=fakeable_Muons_pdgid[0]
          OPS_l2_pdgid=fakeable_Muons_pdgid[1]
          OPS_l1_pt=fakeable_Muons[0].Pt()
          OPS_l1_eta=fakeable_Muons[0].Eta()
          OPS_l1_phi=fakeable_Muons[0].Phi()
          OPS_l1_mass=fakeable_Muons[0].M()
          OPS_l2_pt=fakeable_Muons[1].Pt()
          OPS_l2_eta=fakeable_Muons[1].Eta()
          OPS_l2_phi=fakeable_Muons[1].Phi()
          OPS_l2_mass=fakeable_Muons[1].M()
          OPS_z_mass=(fakeableLeptons[0]+fakeableLeptons[1]).M()
          OPS_z_pt=(fakeableLeptons[0]+fakeableLeptons[1]).Pt()
          OPS_z_eta=(fakeableLeptons[0]+fakeableLeptons[1]).Eta()
          OPS_z_phi=(fakeableLeptons[0]+fakeableLeptons[1]).Phi()
          OPS_drll=fakeableLeptons[0].DeltaR(fakeableLeptons[1])
	  if self.is_mc:
	    if event.Muon_genPartFlav[OPS_l1_id]==1 or event.Muon_genPartFlav[OPS_l1_id]==15:OPS_l1_isprompt=1
	    if event.Muon_genPartFlav[OPS_l2_id]==1 or event.Muon_genPartFlav[OPS_l2_id]==15:OPS_l2_isprompt=1
        # 2 eles case
        if len(fakeable_Electrons)==2 and abs(fakeable_Electrons_pdgid[0]+fakeable_Electrons_pdgid[1])==0:
	  OPS_region=3
          OPS_l1_id=fakeable_Electrons_id[0]
          OPS_l2_id=fakeable_Electrons_id[1]
          OPS_l1_pdgid=fakeable_Electrons_pdgid[0]
          OPS_l2_pdgid=fakeable_Electrons_pdgid[1]
          OPS_l1_pt=fakeable_Electrons[0].Pt()
          OPS_l1_eta=fakeable_Electrons[0].Eta()
          OPS_l1_phi=fakeable_Electrons[0].Phi()
          OPS_l1_mass=fakeable_Electrons[0].M()
          OPS_l2_pt=fakeable_Electrons[1].Pt()
          OPS_l2_eta=fakeable_Electrons[1].Eta()
          OPS_l2_phi=fakeable_Electrons[1].Phi()
          OPS_l2_mass=fakeable_Electrons[1].M()
          OPS_z_mass=(fakeableLeptons[0]+fakeableLeptons[1]).M()
          OPS_z_pt=(fakeableLeptons[0]+fakeableLeptons[1]).Pt()
          OPS_z_eta=(fakeableLeptons[0]+fakeableLeptons[1]).Eta()
          OPS_z_phi=(fakeableLeptons[0]+fakeableLeptons[1]).Phi()
          OPS_drll=fakeableLeptons[0].DeltaR(fakeableLeptons[1])
	  if self.is_mc:
	    if event.Electron_genPartFlav[OPS_l1_id]==1 or event.Electron_genPartFlav[OPS_l1_id]==15:OPS_l1_isprompt=1
	    if event.Electron_genPartFlav[OPS_l2_id]==1 or event.Electron_genPartFlav[OPS_l2_id]==15:OPS_l2_isprompt=1
	# 1 ele case
        if len(fakeable_Electrons)==1 and (sign(fakeable_Muons_pdgid[0])+sign(fakeable_Electrons_pdgid[0]))==0:
          OPS_region=2
          OPS_l1_id=fakeable_Muons_id[0]
          OPS_l2_id=fakeable_Electrons_id[0]
          OPS_l1_pdgid=fakeable_Muons_pdgid[0]
          OPS_l2_pdgid=fakeable_Electrons_pdgid[0]
          OPS_l1_pt=fakeable_Muons[0].Pt()
          OPS_l1_eta=fakeable_Muons[0].Eta()
          OPS_l1_phi=fakeable_Muons[0].Phi()
          OPS_l1_mass=fakeable_Muons[0].M()
          OPS_l2_pt=fakeable_Electrons[0].Pt()
          OPS_l2_eta=fakeable_Electrons[0].Eta()
          OPS_l2_phi=fakeable_Electrons[0].Phi()
          OPS_l2_mass=fakeable_Electrons[0].M()
          OPS_z_mass=(fakeableLeptons[0]+fakeableLeptons[1]).M()
          OPS_z_pt=(fakeableLeptons[0]+fakeableLeptons[1]).Pt()
          OPS_z_eta=(fakeableLeptons[0]+fakeableLeptons[1]).Eta()
          OPS_z_phi=(fakeableLeptons[0]+fakeableLeptons[1]).Phi()
          OPS_drll=fakeableLeptons[0].DeltaR(fakeableLeptons[1])
	  if self.is_mc:
	    if event.Muon_genPartFlav[OPS_l1_id]==1 or event.Muon_genPartFlav[OPS_l1_id]==15:OPS_l1_isprompt=1
	    if event.Electron_genPartFlav[OPS_l2_id]==1 or event.Electron_genPartFlav[OPS_l2_id]==15:OPS_l2_isprompt=1

    self.out.fillBranch("OPS_region", OPS_region)
    self.out.fillBranch("OPS_2P0F", OPS_2P0F)
    self.out.fillBranch("OPS_1P1F", OPS_1P1F)
    self.out.fillBranch("OPS_0P2F", OPS_0P2F)
    self.out.fillBranch("OPS_lep1_faketag", OPS_lep1_faketag)
    self.out.fillBranch("OPS_lep2_faketag", OPS_lep2_faketag)
    self.out.fillBranch("OPS_l1_id", OPS_l1_id)
    self.out.fillBranch("OPS_l2_id", OPS_l2_id)
    self.out.fillBranch("OPS_l1_isprompt", OPS_l1_isprompt)
    self.out.fillBranch("OPS_l2_isprompt", OPS_l2_isprompt)
    self.out.fillBranch("OPS_l1_pdgid", OPS_l1_pdgid)
    self.out.fillBranch("OPS_l2_pdgid", OPS_l2_pdgid)
    self.out.fillBranch("OPS_l1_pt", OPS_l1_pt)
    self.out.fillBranch("OPS_l1_pt_raw", OPS_l1_pt_raw)
    self.out.fillBranch("OPS_l1_eta", OPS_l1_eta)
    self.out.fillBranch("OPS_l1_phi", OPS_l1_phi)
    self.out.fillBranch("OPS_l1_mass", OPS_l1_mass)
    self.out.fillBranch("OPS_l2_pt", OPS_l2_pt)
    self.out.fillBranch("OPS_l2_pt_raw", OPS_l2_pt_raw)
    self.out.fillBranch("OPS_l2_eta", OPS_l2_eta)
    self.out.fillBranch("OPS_l2_phi", OPS_l2_phi)
    self.out.fillBranch("OPS_l2_mass", OPS_l2_mass)
    self.out.fillBranch("OPS_z_mass", OPS_z_mass)
    self.out.fillBranch("OPS_z_pt", OPS_z_pt)
    self.out.fillBranch("OPS_z_eta", OPS_z_eta)
    self.out.fillBranch("OPS_z_phi", OPS_z_phi)
    self.out.fillBranch("OPS_z_mass_raw", OPS_z_mass_raw)
    self.out.fillBranch("OPS_z_pt_raw", OPS_z_pt_raw)
    self.out.fillBranch("OPS_z_eta_raw", OPS_z_eta_raw)
    self.out.fillBranch("OPS_z_phi_raw", OPS_z_phi_raw)
    self.out.fillBranch("OPS_drll", OPS_drll)

    if not (ttc_nl or WZ_region >0 or OPS_region>0):
      return False

    return True

TTC2016 = lambda: TTCProducer("2016")
TTC2017 = lambda: TTCProducer("2017")
TTC2018 = lambda: TTCProducer("2018")
