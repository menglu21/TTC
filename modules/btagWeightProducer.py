#!/usr/bin/env python
# Refer to 1d) method from https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods#1d_Event_reweighting_using_discr
# details: https://twiki.cern.ch/twiki/bin/view/CMS/BTagShapeCalibration
# Need to process after btv official module: https://github.com/cms-nanoAOD/nanoAOD-tools/blob/master/python/postprocessing/modules/btv/btagSFProducer.py
# Also refer to Latino framework: https://github.com/latinos/LatinoAnalysis/blob/master/NanoGardener/python/modules/BTagEventWeightProducer.py
# Noticed that this module should performed before the b tag cut and measure the sum of event weights before and after applying b-tag event weights.

import os, sys
import math
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from importlib import import_module
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor

from PhysicsTools.NanoAODTools.postprocessing.tools import deltaR

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

import math
import os

class btagWeightProduce(Module):
    def __init__(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.systs_shape_corr = []
        for syst in [ 'jes',
                      'lf', 'hf',
                      'hfstats1', 'hfstats2',
                      'lfstats1', 'lfstats2',
                      'cferr1', 'cferr2' ]:
            self.systs_shape_corr.append("up_%s" % syst)
            self.systs_shape_corr.append("down_%s" % syst)
            self.central_and_systs_shape_corr = [ "central" ]
            self.central_and_systs_shape_corr.extend(self.systs_shape_corr)
            self.branchNames_central_and_systs_shape_corr={}
            for central_or_syst in self.central_and_systs_shape_corr:
                if central_or_syst == "central":
                    self.branchNames_central_and_systs_shape_corr[central_or_syst] = "btagWeight"
                else:
                    self.branchNames_central_and_systs_shape_corr[central_or_syst] = "btagWeight_%s" % central_or_syst
                self.out.branch(self.branchNames_central_and_systs_shape_corr[central_or_syst],'F')     

        self.h_neventsgenweighted_in24_btag = ROOT.TH1D('nEventsGenWeighted_in24_btag','nEventsGenWeighted_in24_btag', 10, 0, 10)
        self.h_neventsgenweighted_in47_btag = ROOT.TH1D('nEventsGenWeighted_in47_btag','nEventsGenWeighted_in47_btag', 10, 0, 10)
        self.h_neventsgenweighted_in24 = ROOT.TH1D('nEventsGenWeighted_in24','nEventsGenWeighted_in24', 10, 0, 10)
        self.h_neventsgenweighted_in47 = ROOT.TH1D('nEventsGenWeighted_in47','nEventsGenWeighted_in47', 10, 0, 10)

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        prevdir = ROOT.gDirectory
        outputFile.cd()
        self.h_neventsgenweighted_in24_btag.Write()
        self.h_neventsgenweighted_in47_btag.Write()
        self.h_neventsgenweighted_in24.Write()
        self.h_neventsgenweighted_in47.Write()
        self.h_neventsgenweighted_in24_ratio = self.h_neventsgenweighted_in24.Clone("nEventsGenWeighted_in24_ratio")
        self.h_neventsgenweighted_in24_ratio.Divide(self.h_neventsgenweighted_in24_btag)
        self.h_neventsgenweighted_in47_ratio = self.h_neventsgenweighted_in24.Clone("nEventsGenWeighted_in47_ratio")
        self.h_neventsgenweighted_in47_ratio.Divide(self.h_neventsgenweighted_in47_btag)
        self.h_neventsgenweighted_in24_ratio.Write()
        self.h_neventsgenweighted_in47_ratio.Write()
        prevdir.cd()
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        btag_weight = 1.
        for central_or_syst in self.central_and_systs_shape_corr:
            btag_weight = 1.
            if central_or_syst == "central":
                for i in range(len(event.tightJets_b_DeepCSVmedium_id)):
                    idx = event.tightJets_b_DeepCSVmedium_id[i]
                    if idx < 0: continue
                    btag_weight *= event.Jet_btagSF_deepjet_shape[idx]
            else:
                for i in range(len(event.tightJets_b_DeepCSVmedium_id)):
                    idx = event.tightJets_b_DeepCSVmedium_id[i]
                    if idx < 0: continue
                    btag_weight *= getattr(event, "Jet_btagSF_deepjet_shape_%s" %central_or_syst)[idx]
            self.out.fillBranch(self.branchNames_central_and_systs_shape_corr[central_or_syst], btag_weight)

        if hasattr(event, 'Generator_weight') and event.Generator_weight < 0:
            self.h_neventsgenweighted_in24_btag.Fill(len(event.tightJets_id_in24), -1 * btag_weight)
            self.h_neventsgenweighted_in47_btag.Fill(len(event.tightJets_id_in47), -1 * btag_weight)
            self.h_neventsgenweighted_in24.Fill(len(event.tightJets_id_in24), -1)
            self.h_neventsgenweighted_in47.Fill(len(event.tightJets_id_in47), -1)
        else:
            self.h_neventsgenweighted_in24_btag.Fill(len(event.tightJets_id_in24), 1 * btag_weight)
            self.h_neventsgenweighted_in47_btag.Fill(len(event.tightJets_id_in47), 1 * btag_weight)
            self.h_neventsgenweighted_in24.Fill(len(event.tightJets_id_in24), 1)
            self.h_neventsgenweighted_in47.Fill(len(event.tightJets_id_in47), 1)

        return True

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

btagWeightModule = lambda : btagWeightProduce()