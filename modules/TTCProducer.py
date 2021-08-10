import ROOT
from ROOT import TLorentzVector
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection 
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

import math
import os
import numpy as np
from numpy import sign

class TTCProducer(Module):
  def __init__( self , year ):
    self.year = year
  def beginJob(self):
    pass
  def endJob(self):
    pass
  def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
    self.out = wrappedOutputTree
    self.out.branch("HLT_passEle32WPTight", "I")
    self.out.branch("n_tight_muon", "I")
    self.out.branch("n_loose_muon", "I")
    self.out.branch("n_tight_ele", "I")
    self.out.branch("n_loose_ele", "I")
    self.out.branch("ttc_nl", "B")
    self.out.branch("ttc_jets", "B")
    self.out.branch("n_bjet", "I")
    self.out.branch("n_nobjet", "I")
    self.out.branch("ttc_region", "I")
    self.out.branch("ttc_l1_id", "I")
    self.out.branch("ttc_l2_id", "I")
    self.out.branch("ttc_l1_pt", "F")
    self.out.branch("ttc_l1_eta", "F")
    self.out.branch("ttc_l1_phi", "F")
    self.out.branch("ttc_l2_pt", "F")
    self.out.branch("ttc_l2_eta", "F")
    self.out.branch("ttc_l2_phi", "F")
    self.out.branch("ttc_mll", "F")
    self.out.branch("ttc_drll", "F")
    self.out.branch("ttc_met", "F")
    self.out.branch("WZ_region", "I")
    self.out.branch("WZ_zl1_id", "I")
    self.out.branch("WZ_zl2_id", "I")
    self.out.branch("WZ_wl_id", "I")
    self.out.branch("WZ_zl1_pt", "F")
    self.out.branch("WZ_zl1_eta", "F")
    self.out.branch("WZ_zl1_phi", "F")
    self.out.branch("WZ_zl2_pt", "F")
    self.out.branch("WZ_zl2_eta", "F")
    self.out.branch("WZ_zl2_phi", "F")
    self.out.branch("WZ_l3_pt", "F")
    self.out.branch("WZ_l3_eta", "F")
    self.out.branch("WZ_l3_phi", "F")
    self.out.branch("WZ_Z_mass", "F")
    self.out.branch("WZ_Z_pt", "F")
    self.out.branch("WZ_Z_eta", "F")
    self.out.branch("WZ_Z_phi", "F")
    self.out.branch("WZ_met", "F")
    self.out.branch("DY_region", "I")
    self.out.branch("DY_l1_id", "I")
    self.out.branch("DY_l2_id", "I")
    self.out.branch("DY_l1_pt", "F")
    self.out.branch("DY_l1_eta", "F")
    self.out.branch("DY_l1_phi", "F")
    self.out.branch("DY_l2_pt", "F")
    self.out.branch("DY_l2_eta", "F")
    self.out.branch("DY_l2_phi", "F")
    self.out.branch("DY_z_mass", "F")
    self.out.branch("DY_z_pt", "F")
    self.out.branch("DY_z_eta", "F")
    self.out.branch("DY_z_phi", "F")
    self.out.branch("DY_drll", "F")
    self.out.branch("tightJets_id_in24","I",lenVar="nJet")
    self.out.branch("tightJets_id_in47","I",lenVar="nJet")
    self.out.branch("tightJets_nob_CSVmedium_id","I",lenVar="nJet")
    self.out.branch("tightJets_b_CSVmedium_id","I",lenVar="nJet")
    self.out.branch("tightJets_nob_DeepCSVmedium_id","I",lenVar="nJet")
    self.out.branch("tightJets_b_DeepCSVmedium_id","I",lenVar="nJet")
    if event.nElectron>0:
      self.out.branch("tightElectrons_id","I",lenVar="nElectron")
      self.out.branch("additional_vetoElectrons_id","I",lenVar="nElectron")
    if event.nMuon>0:
      self.out.branch("tightMuons_id","I",lenVar="nMuon")
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
    tightMuons = []
    tightMuons_pdgid = []
    tightMuons_id = []
    additional_looseMuons = []
    additional_looseMuons_pdgid = []
    additional_looseMuons_id = []
    for imu in range(0, event.nMuon):
      if (muons[imu].tightId):
        if (muons[imu].pfRelIso04_all<0.15 and abs(muons[imu].eta)<2.4 and muons[imu].tightCharge==2 and event.Muon_corrected_pt[imu]>15):
          muon_v4_temp.SetPtEtaPhiM(event.Muon_corrected_pt[imu], muons[imu].eta, muons[imu].phi, muons[imu].mass)
          tightMuons.append(muon_v4_temp.Clone())
          tightMuons_pdgid.append(muons[imu].pdgId)
          tightMuons_id.append(imu)
      elif (muons[imu].looseId):
        if (muons[imu].pfRelIso04_all<0.25 and abs(muons[imu].eta)<2.4 and event.Muon_corrected_pt[imu]>15):
          muon_v4_temp.SetPtEtaPhiM(event.Muon_corrected_pt[imu], muons[imu].eta, muons[imu].phi, muons[imu].mass)
          additional_looseMuons.append(muon_v4_temp.Clone())
          additional_looseMuons_pdgid.append(muons[imu].pdgId)
          additional_looseMuons_id.append(imu)

    n_tight_muon = len(tightMuons)
    n_loose_muon = len(additional_looseMuons)
    tightMuons_id.extend(np.zeros(event.nMuon-len(tightMuons_id),int)-1)
    additional_looseMuons_id.extend(np.zeros(event.nMuon-len(additional_looseMuons_id),int)-1)
    self.out.fillBranch("n_tight_muon", n_tight_muon)
    self.out.fillBranch("n_loose_muon", n_loose_muon)
    self.out.fillBranch("tightMuons_id", tightMuons_id)
    self.out.fillBranch("additional_looseMuons_id", additional_looseMuons_id)

    # electron selection: tight (veto) cut-based ID + impact parameter cut, with pt > 15 GeV
    electrons = Collection(event, 'Electron')
    electron_v4_temp=TLorentzVector()
    tightElectrons = []
    tightElectrons_pdgid = []
    tightElectrons_id = []
    additional_vetoElectrons = []
    additional_vetoElectrons_pdgid = []
    additional_vetoElectrons_id = []
    for iele in range(0, event.nElectron):
      if (electrons[iele].cutBased==4):
        if (electrons[iele].tightCharge==2 and ((abs(electrons[iele].eta+electrons[iele].deltaEtaSC) <1.4442 and abs(electrons[iele].dxy)<0.05 and abs(electrons[iele].dz<0.1)) or (abs(electrons[iele].eta + electrons[iele].deltaEtaSC)>1.566 and abs(electrons[iele].eta + electrons[iele].deltaEtaSC)<2.4 and abs(electrons[iele].dxy)<0.1 and abs(electrons[iele].dz)<0.2)) and electrons[iele].pt>15):
          electron_v4_temp.SetPtEtaPhiM(electrons[iele].pt, electrons[iele].eta, electrons[iele].phi, electrons[iele].mass)
          tightElectrons.append(electron_v4_temp.Clone())
          tightElectrons_pdgid.append(electrons[iele].pdgId)
          tightElectrons_id.append(iele)
      elif (electrons[iele].cutBased>0):
        if (((abs(electrons[iele].eta+electrons[iele].deltaEtaSC) <1.4442 and abs(electrons[iele].dxy)<0.05 and abs(electrons[iele].dz<0.1)) or (abs(electrons[iele].eta + electrons[iele].deltaEtaSC)>1.566 and abs(electrons[iele].eta + electrons[iele].deltaEtaSC)<2.4 and abs(electrons[iele].dxy)<0.1 and abs(electrons[iele].dz)<0.2)) and electrons[iele].pt>15):
          electron_v4_temp.SetPtEtaPhiM(electrons[iele].pt, electrons[iele].eta, electrons[iele].phi, electrons[iele].mass)
          additional_vetoElectrons.append(electron_v4_temp.Clone())
          additional_vetoElectrons_pdgid.append(electrons[iele].pdgId)
          additional_vetoElectrons_id.append(iele)

    n_tight_ele = len(tightElectrons)
    n_loose_ele = len(additional_vetoElectrons)
    tightElectrons_id.extend(np.zeros(event.nElectron-len(tightElectrons_id),int)-1)
    additional_vetoElectrons_id.extend(np.zeros(event.nElectron-len(additional_vetoElectrons_id),int)-1)
    self.out.fillBranch("n_tight_ele", n_tight_ele)
    self.out.fillBranch("n_loose_ele", n_loose_ele)
    self.out.fillBranch("tightElectrons_id", tightElectrons_id)
    self.out.fillBranch("additional_vetoElectrons_id", additional_vetoElectrons_id)

    # https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL17
    # tight ak4 jets, 2016 (111=7), 2017/2018 (110=6), medium B-tag WP
    # DeepCSV=(nanoaod btagDeepB) loose: 0.1355, medium: 0.4506, tight: 0.7738
    # DeepFlavor=(nanoaod btagDeepFlavB) loose: 0.0532, medium: 0.3040, tight: 0.7476

    # c-jet tag is based on two-D cuts, medium DeepJet WP:
    # CvsL=btagDeepFlavCvL: 0.085, CvsB=btagDeepFlavCvB: 0.34
    # c-tag not available in NANOAOD yet

    jets = Collection(event, 'Jet')

    tightJets_id_in24 = []
    tightJets_id_in47 = []

    tightJets_nob_CSVmedium_id = []
    tightJets_b_CSVmedium_id = []

    tightJets_nob_DeepCSVmedium_id = []
    tightJets_b_DeepCSVmedium_id = []

    for ijet in range(0, event.nJet):
      if self.year=="2016":
        if jets[ijet].jetId==7 and jets[ijet].pt>30:
	  if abs(jets[ijet].eta)<4.7 and abs(jets[ijet].eta)>=2.4: 
	    tightJets_id_in47.append(ijet)
	  if abs(jets[ijet].eta)<2.4:
	    tightJets_id_in24.append(ijet)

            if (jets[ijet].btagCSVV2 > 0.8838):
              tightJets_b_CSVmedium_id.append(ijet)
            else:tightJets_nob_CSVmedium_id.append(ijet)
            if (jets[ijet].btagDeepB > 0.4941):
              tightJets_b_DeepCSVmedium_id.append(ijet)
            else:tightJets_nob_DeepCSVmedium_id.append(ijet)

      elif (self.year=="2017" or self.year=="2018"):
	if jets[ijet].jetId==6 and jets[ijet].pt>30:
	  if abs(jets[ijet].eta)<4.7 and abs(jets[ijet].eta)>=2.4:
	    tightJets_id_in47.append(ijet)
	  if abs(jets[ijet].eta)<2.4:
            tightJets_id_in24.append(ijet)

            if (jets[ijet].btagDeepB > 0.4506):
              tightJets_b_CSVmedium_id.append(ijet)
            else:tightJets_nob_CSVmedium_id.append(ijet)
            if (jets[ijet].btagDeepFlavB > 0.3040):
              tightJets_b_DeepCSVmedium_id.append(ijet)
            else:tightJets_nob_DeepCSVmedium_id.append(ijet)

    tightJets_id_in24.extend(np.zeros(event.nJet-len(tightJets_id_in24),int)-1)
    tightJets_id_in47.extend(np.zeros(event.nJet-len(tightJets_id_in47),int)-1)
    tightJets_b_CSVmedium_id.extend(np.zeros(event.nJet-len(tightJets_b_CSVmedium_id),int)-1)
    tightJets_nob_CSVmedium_id.extend(np.zeros(event.nJet-len(tightJets_nob_CSVmedium_id),int)-1)
    tightJets_b_DeepCSVmedium_id.extend(np.zeros(event.nJet-len(tightJets_b_DeepCSVmedium_id),int)-1)
    tightJets_nob_DeepCSVmedium_id.extend(np.zeros(event.nJet-len(tightJets_nob_DeepCSVmedium_id),int)-1)

    self.out.fillBranch("tightJets_id_in24",tightJets_id_in24)
    self.out.fillBranch("tightJets_id_in47",tightJets_id_in47)
    self.out.fillBranch("tightJets_nob_CSVmedium_id",tightJets_nob_CSVmedium_id)
    self.out.fillBranch("tightJets_b_CSVmedium_id",tightJets_b_CSVmedium_id)
    self.out.fillBranch("tightJets_nob_DeepCSVmedium_id",tightJets_nob_DeepCSVmedium_id)
    self.out.fillBranch("tightJets_b_DeepCSVmedium_id",tightJets_b_DeepCSVmedium_id)

    # tight leptons and additional loose leptons collection
    tightLeptons = tightMuons + tightElectrons
    tightLeptons.sort(key=lambda x: x.Pt(), reverse=True)
    looseLeptons = additional_looseMuons + additional_vetoElectrons
    looseLeptons.sort(key=lambda x: x.Pt(), reverse=True)

    if len(tightLeptons)<2:return False

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
    #ttc region tag, 1:2 muon, 2:1 muon, 3:0 muon
    ttc_region=0
    ttc_l1_id=-1
    ttc_l2_id=-1
    ttc_l1_pt=-99
    ttc_l1_eta=-99
    ttc_l1_phi=-99
    ttc_l2_pt=-99
    ttc_l2_eta=-99
    ttc_l2_phi=-99
    ttc_l3_pt=-99
    ttc_l3_eta=-99
    ttc_l3_phi=-99
    ttc_mll=-99
    ttc_drll=-99
    ttc_met=-99
    
    # the two leptons with pt >30/20, 3th lepton veto
    if len(tightLeptons)==2 and tightLeptons[0].Pt()>30 and tightLeptons[1].Pt()>20 and len(looseLeptons)==0:
      ttc_nl=True
    # at least three jets and one bjets
    n_bjet=0
    n_nobjet=0
    for ijet in range(0,len(tightJets_b_DeepCSVmedium_id)):
      if tightJets_b_DeepCSVmedium_id[ijet]>-1:
	n_bjet=n_bjet+1
    for ijet in range(0,len(tightJets_nob_DeepCSVmedium_id)):
      if tightJets_nob_DeepCSVmedium_id[ijet]>-1:
        n_nobjet=n_nobjet+1

    if ttc_nl and (n_nobjet+n_bjet)>2 and n_bjet>0:
      ttc_jets=True

    if ttc_nl:
      if self.is_mc:
        ttc_met=event.MET_T1Smear_pt
      else:
        ttc_met=event.MET_T1_pt
      if len(tightElectrons)==0 and abs(tightMuons_pdgid[0]+tightMuons_pdgid[1])==26:
	ttc_region=1
	ttc_l1_id=tightMuons_id[0]
	ttc_l2_id=tightMuons_id[1]
	ttc_l1_pt=tightMuons[0].Pt()
	ttc_l1_eta=tightMuons[0].Eta()
	ttc_l1_phi=tightMuons[0].Phi()
	ttc_l2_pt=tightMuons[1].Pt()
	ttc_l2_eta=tightMuons[1].Eta()
	ttc_l2_phi=tightMuons[1].Phi()
	ttc_mll=(tightMuons[0]+tightMuons[1]).M()
	ttc_drll=tightMuons[0].DeltaR(tightMuons[1])
      if len(tightElectrons)==1 and abs(tightMuons_pdgid[0]+tightElectrons_pdgid[0])==24:
	ttc_region=2
	ttc_l1_id=tightMuons_id[0]
	ttc_l2_id=tightElectrons_id[0]
	ttc_l1_pt=tightMuons[0].Pt()
        ttc_l1_eta=tightMuons[0].Eta()
        ttc_l1_phi=tightMuons[0].Phi()
        ttc_l2_pt=tightElectrons[0].Pt()
        ttc_l2_eta=tightElectrons[0].Eta()
        ttc_l2_phi=tightElectrons[0].Phi()
        ttc_mll=(tightMuons[0]+tightElectrons[0]).M()
        ttc_drll=tightMuons[0].DeltaR(tightElectrons[0])
      if len(tightElectrons)==2 and abs(tightElectrons_pdgid[0]+tightElectrons_pdgid[1])==22:
	ttc_region=3
	ttc_l1_id=tightElectrons_id[0]
	ttc_l2_id=tightElectrons_id[1]
	ttc_l1_pt=tightElectrons[0].Pt()
        ttc_l1_eta=tightElectrons[0].Eta()
        ttc_l1_phi=tightElectrons[0].Phi()
        ttc_l2_pt=tightElectrons[1].Pt()
        ttc_l2_eta=tightElectrons[1].Eta()
        ttc_l2_phi=tightElectrons[1].Phi()
        ttc_mll=(tightElectrons[0]+tightElectrons[1]).M()
        ttc_drll=tightElectrons[0].DeltaR(tightElectrons[1])

    self.out.fillBranch("ttc_nl", ttc_nl)
    self.out.fillBranch("ttc_jets", ttc_jets)
    self.out.fillBranch("n_bjet", n_bjet)
    self.out.fillBranch("n_nobjet", n_nobjet)
    self.out.fillBranch("ttc_region", ttc_region)
    self.out.fillBranch("ttc_l1_id", ttc_l1_id)
    self.out.fillBranch("ttc_l2_id", ttc_l2_id)
    self.out.fillBranch("ttc_l1_pt", ttc_l1_pt) 
    self.out.fillBranch("ttc_l1_eta", ttc_l1_eta) 
    self.out.fillBranch("ttc_l1_phi", ttc_l1_phi) 
    self.out.fillBranch("ttc_l2_pt", ttc_l2_pt) 
    self.out.fillBranch("ttc_l2_eta", ttc_l2_eta) 
    self.out.fillBranch("ttc_l2_phi", ttc_l2_phi) 
    self.out.fillBranch("ttc_mll",ttc_mll) 
    self.out.fillBranch("ttc_drll",ttc_drll) 
    self.out.fillBranch("ttc_met", ttc_met)

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
    WZ_zl1_pt=-99
    WZ_zl1_eta=-99
    WZ_zl1_phi=-99
    WZ_zl2_pt=-99
    WZ_zl2_eta=-99
    WZ_zl2_phi=-99
    WZ_l3_pt=-99
    WZ_l3_eta=-99
    WZ_l3_phi=-99
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
	    WZ_zl1_pt=tightMuons[0].Pt()
	    WZ_zl1_eta=tightMuons[0].Eta()
	    WZ_zl1_phi=tightMuons[0].Phi()
	    WZ_zl2_pt=tightMuons[2].Pt()
	    WZ_zl2_eta=tightMuons[2].Eta()
	    WZ_zl2_phi=tightMuons[2].Phi()
	    WZ_l3_pt=tightMuons[1].Pt()
	    WZ_l3_eta=tightMuons[1].Eta()
	    WZ_l3_phi=tightMuons[1].Phi()
	    WZ_Z_mass=(tightMuons[0]+tightMuons[2]).M()
	    WZ_Z_pt=(tightMuons[0]+tightMuons[2]).Pt()
	    WZ_Z_eta=(tightMuons[0]+tightMuons[2]).Eta()
	    WZ_Z_phi=(tightMuons[0]+tightMuons[2]).Phi()

	  if abs((tightMuons[0]+tightMuons[2]).M()-91.1876)>abs((tightMuons[1]+tightMuons[2]).M()-91.1876) and abs((tightMuons[1]+tightMuons[2]).M()-91.1876)<15:
	    WZ_region=1
            WZ_zl1_id=tightMuons_id[1]
            WZ_zl2_id=tightMuons_id[2]
            WZ_wl_id=tightMuons_id[0]
	    WZ_zl1_pt=tightMuons[1].Pt()
	    WZ_zl1_eta=tightMuons[1].Eta()
	    WZ_zl1_phi=tightMuons[1].Phi()
	    WZ_zl2_pt=tightMuons[2].Pt()
	    WZ_zl2_eta=tightMuons[2].Eta()
	    WZ_zl2_phi=tightMuons[2].Phi()
	    WZ_l3_pt=tightMuons[0].Pt()
	    WZ_l3_eta=tightMuons[0].Eta()
	    WZ_l3_phi=tightMuons[0].Phi()
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
	    WZ_zl1_pt=tightMuons[0].Pt()
	    WZ_zl1_eta=tightMuons[0].Eta()
	    WZ_zl1_phi=tightMuons[0].Phi()
	    WZ_zl2_pt=tightMuons[1].Pt()
	    WZ_zl2_eta=tightMuons[1].Eta()
	    WZ_zl2_phi=tightMuons[1].Phi()
	    WZ_l3_pt=tightMuons[2].Pt()
	    WZ_l3_eta=tightMuons[2].Eta()
	    WZ_l3_phi=tightMuons[2].Phi()
	    WZ_Z_mass=(tightMuons[0]+tightMuons[1]).M()
	    WZ_Z_pt=(tightMuons[0]+tightMuons[1]).Pt()
	    WZ_Z_eta=(tightMuons[0]+tightMuons[1]).Eta()
	    WZ_Z_phi=(tightMuons[0]+tightMuons[1]).Phi()

	  if abs((tightMuons[0]+tightMuons[1]).M()-91.1876)>abs((tightMuons[1]+tightMuons[2]).M()-91.1876) and abs((tightMuons[1]+tightMuons[2]).M()-91.1876)<15:
	    WZ_region=1
            WZ_zl1_id=tightMuons_id[1]
            WZ_zl2_id=tightMuons_id[2]
            WZ_wl_id=tightMuons_id[0]
	    WZ_zl1_pt=tightMuons[1].Pt()
	    WZ_zl1_eta=tightMuons[1].Eta()
	    WZ_zl1_phi=tightMuons[1].Phi()
	    WZ_zl2_pt=tightMuons[2].Pt()
	    WZ_zl2_eta=tightMuons[2].Eta()
	    WZ_zl2_phi=tightMuons[2].Phi()
	    WZ_l3_pt=tightMuons[0].Pt()
	    WZ_l3_eta=tightMuons[0].Eta()
	    WZ_l3_phi=tightMuons[0].Phi()
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
	    WZ_zl1_pt=tightMuons[0].Pt()
	    WZ_zl1_eta=tightMuons[0].Eta()
	    WZ_zl1_phi=tightMuons[0].Phi()
	    WZ_zl2_pt=tightMuons[1].Pt()
	    WZ_zl2_eta=tightMuons[1].Eta()
	    WZ_zl2_phi=tightMuons[1].Phi()
	    WZ_l3_pt=tightMuons[2].Pt()
	    WZ_l3_eta=tightMuons[2].Eta()
	    WZ_l3_phi=tightMuons[2].Phi()
	    WZ_Z_mass=(tightMuons[0]+tightMuons[1]).M()
	    WZ_Z_pt=(tightMuons[0]+tightMuons[1]).Pt()
	    WZ_Z_eta=(tightMuons[0]+tightMuons[1]).Eta()
	    WZ_Z_phi=(tightMuons[0]+tightMuons[1]).Phi()
	  if abs((tightMuons[0]+tightMuons[1]).M()-91.1876)>abs((tightMuons[0]+tightMuons[2]).M()-91.1876) and abs((tightMuons[0]+tightMuons[2]).M()-91.1876)<15:
	    WZ_region=1
            WZ_zl1_id=tightMuons_id[0]
            WZ_zl2_id=tightMuons_id[2]
            WZ_wl_id=tightMuons_id[1]
	    WZ_zl1_pt=tightMuons[0].Pt()
	    WZ_zl1_eta=tightMuons[0].Eta()
	    WZ_zl1_phi=tightMuons[0].Phi()
	    WZ_zl2_pt=tightMuons[2].Pt()
	    WZ_zl2_eta=tightMuons[2].Eta()
	    WZ_zl2_phi=tightMuons[2].Phi()
	    WZ_l3_pt=tightMuons[1].Pt()
	    WZ_l3_eta=tightMuons[1].Eta()
	    WZ_l3_phi=tightMuons[1].Phi()
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
	  WZ_zl1_pt=tightMuons[0].Pt()
	  WZ_zl1_eta=tightMuons[0].Eta()
	  WZ_zl1_phi=tightMuons[0].Phi()
	  WZ_zl2_pt=tightMuons[1].Pt()
	  WZ_zl2_eta=tightMuons[1].Eta()
	  WZ_zl2_phi=tightMuons[1].Phi()
	  WZ_l3_pt=tightElectrons[0].Pt()
	  WZ_l3_eta=tightElectrons[0].Eta()
	  WZ_l3_phi=tightElectrons[0].Phi()
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
	  WZ_zl1_pt=tightElectrons[0].Pt()
	  WZ_zl1_eta=tightElectrons[0].Eta()
	  WZ_zl1_phi=tightElectrons[0].Phi()
	  WZ_zl2_pt=tightElectrons[1].Pt()
	  WZ_zl2_eta=tightElectrons[1].Eta()
	  WZ_zl2_phi=tightElectrons[1].Phi()
	  WZ_l3_pt=tightMuons[0].Pt()
	  WZ_l3_eta=tightMuons[0].Eta()
	  WZ_l3_phi=tightMuons[0].Phi()
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
	    WZ_zl1_pt=tightElectrons[0].Pt()
	    WZ_zl1_eta=tightElectrons[0].Eta()
	    WZ_zl1_phi=tightElectrons[0].Phi()
	    WZ_zl2_pt=tightElectrons[2].Pt()
	    WZ_zl2_eta=tightElectrons[2].Eta()
	    WZ_zl2_phi=tightElectrons[2].Phi()
	    WZ_l3_pt=tightElectrons[1].Pt()
	    WZ_l3_eta=tightElectrons[1].Eta()
	    WZ_l3_phi=tightElectrons[1].Phi()
	    WZ_Z_mass=(tightElectrons[0]+tightElectrons[2]).M()
	    WZ_Z_pt=(tightElectrons[0]+tightElectrons[2]).Pt()
	    WZ_Z_eta=(tightElectrons[0]+tightElectrons[2]).Eta()
	    WZ_Z_phi=(tightElectrons[0]+tightElectrons[2]).Phi()
	  if abs((tightElectrons[0]+tightElectrons[2]).M()-91.1876)>abs((tightElectrons[1]+tightElectrons[2]).M()-91.1876) and abs((tightElectrons[1]+tightElectrons[2]).M()-91.1876)<15:
	    WZ_region=4
            WZ_zl1_id=tightElectrons_id[1]
            WZ_zl2_id=tightElectrons_id[2]
            WZ_wl_id=tightElectrons_id[0]
	    WZ_zl1_pt=tightElectrons[1].Pt()
	    WZ_zl1_eta=tightElectrons[1].Eta()
	    WZ_zl1_phi=tightElectrons[1].Phi()
	    WZ_zl2_pt=tightElectrons[2].Pt()
	    WZ_zl2_eta=tightElectrons[2].Eta()
	    WZ_zl2_phi=tightElectrons[2].Phi()
	    WZ_l3_pt=tightElectrons[0].Pt()
	    WZ_l3_eta=tightElectrons[0].Eta()
	    WZ_l3_phi=tightElectrons[0].Phi()
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
	    WZ_zl1_pt=tightElectrons[0].Pt()
	    WZ_zl1_eta=tightElectrons[0].Eta()
	    WZ_zl1_phi=tightElectrons[0].Phi()
	    WZ_zl2_pt=tightElectrons[1].Pt()
	    WZ_zl2_eta=tightElectrons[1].Eta()
	    WZ_zl2_phi=tightElectrons[1].Phi()
	    WZ_l3_pt=tightElectrons[2].Pt()
	    WZ_l3_eta=tightElectrons[2].Eta()
	    WZ_l3_phi=tightElectrons[2].Phi()
	    WZ_Z_mass=(tightElectrons[0]+tightElectrons[1]).M()
	    WZ_Z_pt=(tightElectrons[0]+tightElectrons[1]).Pt()
	    WZ_Z_eta=(tightElectrons[0]+tightElectrons[1]).Eta()
	    WZ_Z_phi=(tightElectrons[0]+tightElectrons[1]).Phi()
	  if abs((tightElectrons[0]+tightElectrons[1]).M()-91.1876)>abs((tightElectrons[1]+tightElectrons[2]).M()-91.1876) and abs((tightElectrons[1]+tightElectrons[2]).M()-91.1876)<15:
	    WZ_region=4
            WZ_zl1_id=tightElectrons_id[1]
            WZ_zl2_id=tightElectrons_id[2]
            WZ_wl_id=tightElectrons_id[0]
	    WZ_zl1_pt=tightElectrons[1].Pt()
	    WZ_zl1_eta=tightElectrons[1].Eta()
	    WZ_zl1_phi=tightElectrons[1].Phi()
	    WZ_zl2_pt=tightElectrons[2].Pt()
	    WZ_zl2_eta=tightElectrons[2].Eta()
	    WZ_zl2_phi=tightElectrons[2].Phi()
	    WZ_l3_pt=tightElectrons[0].Pt()
	    WZ_l3_eta=tightElectrons[0].Eta()
	    WZ_l3_phi=tightElectrons[0].Phi()
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
	    WZ_zl1_pt=tightElectrons[0].Pt()
	    WZ_zl1_eta=tightElectrons[0].Eta()
	    WZ_zl1_phi=tightElectrons[0].Phi()
	    WZ_zl2_pt=tightElectrons[1].Pt()
	    WZ_zl2_eta=tightElectrons[1].Eta()
	    WZ_zl2_phi=tightElectrons[1].Phi()
	    WZ_l3_pt=tightElectrons[2].Pt()
	    WZ_l3_eta=tightElectrons[2].Eta()
	    WZ_l3_phi=tightElectrons[2].Phi()
	    WZ_Z_mass=(tightElectrons[0]+tightElectrons[1]).M()
	    WZ_Z_pt=(tightElectrons[0]+tightElectrons[1]).Pt()
	    WZ_Z_eta=(tightElectrons[0]+tightElectrons[1]).Eta()
	    WZ_Z_phi=(tightElectrons[0]+tightElectrons[1]).Phi()
	  if abs((tightElectrons[0]+tightElectrons[1]).M()-91.1876)>abs((tightElectrons[0]+tightElectrons[2]).M()-91.1876) and abs((tightElectrons[0]+tightElectrons[2]).M()-91.1876)<15:
	    WZ_region=4
            WZ_zl1_id=tightElectrons_id[0]
            WZ_zl2_id=tightElectrons_id[2]
            WZ_wl_id=tightElectrons_id[1]
	    WZ_zl1_pt=tightElectrons[0].Pt()
	    WZ_zl1_eta=tightElectrons[0].Eta()
	    WZ_zl1_phi=tightElectrons[0].Phi()
	    WZ_zl2_pt=tightElectrons[2].Pt()
	    WZ_zl2_eta=tightElectrons[2].Eta()
	    WZ_zl2_phi=tightElectrons[2].Phi()
	    WZ_l3_pt=tightElectrons[1].Pt()
	    WZ_l3_eta=tightElectrons[1].Eta()
	    WZ_l3_phi=tightElectrons[1].Phi()
	    WZ_Z_mass=(tightElectrons[0]+tightElectrons[2]).M()
	    WZ_Z_pt=(tightElectrons[0]+tightElectrons[2]).Pt()
	    WZ_Z_eta=(tightElectrons[0]+tightElectrons[2]).Eta()
	    WZ_Z_phi=(tightElectrons[0]+tightElectrons[2]).Phi()
	
    self.out.fillBranch("WZ_region", WZ_region)
    self.out.fillBranch("WZ_zl1_id", WZ_zl1_id)
    self.out.fillBranch("WZ_zl2_id", WZ_zl2_id)
    self.out.fillBranch("WZ_wl_id", WZ_wl_id)
    self.out.fillBranch("WZ_zl1_pt", WZ_zl1_pt)
    self.out.fillBranch("WZ_zl1_eta", WZ_zl1_eta)
    self.out.fillBranch("WZ_zl1_phi", WZ_zl1_phi)
    self.out.fillBranch("WZ_zl2_pt", WZ_zl2_pt)
    self.out.fillBranch("WZ_zl2_eta", WZ_zl2_eta)
    self.out.fillBranch("WZ_zl2_phi", WZ_zl2_phi)
    self.out.fillBranch("WZ_l3_pt", WZ_l3_pt)
    self.out.fillBranch("WZ_l3_eta", WZ_l3_eta)
    self.out.fillBranch("WZ_l3_phi", WZ_l3_phi)
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
    #DY region b-jet selection
    DY_nb=False
    #DY region tag, 0: fail to pass the DY selection, 1:2 muon, 2:1 muon, 3:0 muon
    DY_region=0
    DY_l1_id=-1
    DY_l2_id=-1
    DY_l1_pt=-99
    DY_l1_eta=-99
    DY_l1_phi=-99
    DY_l2_pt=-99
    DY_l2_eta=-99
    DY_l2_phi=-99
    DY_z_mass=-99
    DY_z_pt=-99
    DY_z_eta=-99
    DY_z_phi=-99
    DY_drll=-99

    # the two leptons with pt >20, 3th lepton veto
    if len(tightLeptons)==2 and tightLeptons[1].Pt()>20 and len(looseLeptons)==0:
      DY_nl=True
    # no bjet
    if DY_nl and tightJets_b_DeepCSVmedium_id[0]==-1:
      DY_nb=True
    if DY_nb:
      # 2 muons case
      if len(tightElectrons)==0 and abs(tightMuons_pdgid[0]+tightMuons_pdgid[1])==0:
	DY_region=1
	DY_l1_id=tightMuons_id[0]
	DY_l2_id=tightMuons_id[1]
	DY_l1_pt=tightMuons[0].Pt()
	DY_l1_eta=tightMuons[0].Eta()
	DY_l1_phi=tightMuons[0].Phi()
	DY_l2_pt=tightMuons[1].Pt()
	DY_l2_eta=tightMuons[1].Eta()
	DY_l2_phi=tightMuons[1].Phi()
	DY_z_mass=(tightLeptons[0]+tightLeptons[1]).M()
	DY_z_pt=(tightLeptons[0]+tightLeptons[1]).Pt()
	DY_z_eta=(tightLeptons[0]+tightLeptons[1]).Eta()
	DY_z_phi=(tightLeptons[0]+tightLeptons[1]).Phi()
	DY_drll=tightLeptons[0].DeltaR(tightLeptons[1])
      # 2 eles case
      if len(tightElectrons)==2 and abs(tightElectrons_pdgid[0]+tightElectrons_pdgid[1])==0:
	DY_region=3
        DY_l1_id=tightElectrons_id[0]
        DY_l2_id=tightElectrons_id[1]
	DY_l1_pt=tightElectrons[0].Pt()
	DY_l1_eta=tightElectrons[0].Eta()
	DY_l1_phi=tightElectrons[0].Phi()
	DY_l2_pt=tightElectrons[1].Pt()
	DY_l2_eta=tightElectrons[1].Eta()
	DY_l2_phi=tightElectrons[1].Phi()
	DY_z_mass=(tightLeptons[0]+tightLeptons[1]).M()
	DY_z_pt=(tightLeptons[0]+tightLeptons[1]).Pt()
	DY_z_eta=(tightLeptons[0]+tightLeptons[1]).Eta()
	DY_z_phi=(tightLeptons[0]+tightLeptons[1]).Phi()
	DY_drll=tightLeptons[0].DeltaR(tightLeptons[1])
      # 1 ele case
      if len(tightElectrons)==1 and (sign(tightMuons_pdgid[0])+sign(tightElectrons_pdgid[0]))==0:
	DY_region=2
        DY_l1_id=tightMuons_id[0]
        DY_l2_id=tightElectrons_id[0]
	DY_l1_pt=tightLeptons[0].Pt()
	DY_l1_eta=tightLeptons[0].Eta()
	DY_l1_phi=tightLeptons[0].Phi()
	DY_l2_pt=tightLeptons[1].Pt()
	DY_l2_eta=tightLeptons[1].Eta()
	DY_l2_phi=tightLeptons[1].Phi()
	DY_z_mass=(tightLeptons[0]+tightLeptons[1]).M()
	DY_z_pt=(tightLeptons[0]+tightLeptons[1]).Pt()
	DY_z_eta=(tightLeptons[0]+tightLeptons[1]).Eta()
	DY_z_phi=(tightLeptons[0]+tightLeptons[1]).Phi()
	DY_drll=tightLeptons[0].DeltaR(tightLeptons[1])

    self.out.fillBranch("DY_region", DY_region)
    self.out.fillBranch("DY_l1_id", DY_l1_id)
    self.out.fillBranch("DY_l2_id", DY_l2_id)
    self.out.fillBranch("DY_l1_pt", DY_l1_pt)
    self.out.fillBranch("DY_l1_eta", DY_l1_eta)
    self.out.fillBranch("DY_l1_phi", DY_l1_phi)
    self.out.fillBranch("DY_l2_pt", DY_l2_pt)
    self.out.fillBranch("DY_l2_eta", DY_l2_eta)
    self.out.fillBranch("DY_l2_phi", DY_l2_phi)
    self.out.fillBranch("DY_z_mass", DY_z_mass)
    self.out.fillBranch("DY_z_pt", DY_z_pt)
    self.out.fillBranch("DY_z_eta", DY_z_eta)
    self.out.fillBranch("DY_z_phi", DY_z_phi)
    self.out.fillBranch("DY_drll", DY_drll)

    if not (ttc_nl or WZ_region >0 or DY_region>0):
      return False

    return True

TTC2016 = lambda: TTCProducer("2016")
TTC2017 = lambda: TTCProducer("2017")
TTC2018 = lambda: TTCProducer("2018")
