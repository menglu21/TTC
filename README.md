# TTC

##Steps for setup:
1. release CMSSW_10_6_27
2. cd $CMSSW_BASE/src
   git clone https://github.com/cms-nanoAOD/nanoAOD-tools.git PhysicsTools/NanoAODTools
   cd PhysicsTools/NanoAODTools
   cmsenv
   scram b

3. cd python/postprocessing
   ##clone this repository
   git clone https://github.com/menglu21/TTC.git analysis
   cd $CMSSW_BASE/src
   scram b

## submit jobs
cd analysis/crab
using the configure files under 'configs', namely,
crab submit -c configs/DoubleEGB_cfg.py

##corrections
the modules (most of them are corrections) used can be seen from analysis/crab/crab_script.py, 
N.B. the egamma correction is already applied default in NanoAOD
####for MC:
countHistogramsModule(): store the opsitive and negative events number for weight apply
puWeight_2017(): pileup reweight
PrefCorr(): L1-prefiring correction
muonIDISOSF2017(): muon ID/ISO SF
muonScaleRes2017(): muon momentum correction, i.e., the Rochester correction
eleRECOSF2017(): electron RECO SF
eleIDSF2017(): electron IS SF
jmeCorrections_UL2017MC(): JetMET correction
btagSF2017UL(): b tag SF
####for Data:
muonScaleRes2017(): muon momentum correction, i.e., the Rochester correction
jmeCorrections_UL2017*(): JetMET correction

1. pileup reweight (this correction is applied using the official module, so we need to update the rootfiles for pileup and do some modification on the official module. The files under others/for_pileup/ can be used directly)
###data
according to https://twiki.cern.ch/twiki/bin/view/CMS/PileupJSONFileforData#Centrally_produced_ROOT_histogra, use histograms under /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PileUp/UltraLegacy/, combine three histograms to a single one with name “pileup, pileup_plus, pileup_minus”
###MC
https://twiki.cern.ch/twiki/bin/view/CMS/PileupScenariosRun2

move "mcPileupUL2017.root" and "PileupHistogram-goldenJSON-13tev-UL2017-99bins_withVar.root" to python/postprocessing/data/pileup/, and move "puWeightProducer.py" to python/postprocessing/modules/common/

2. prefiring correction (needed files are in others/for_prefiring, can be used directly)
details are here: Pre-firing: https://twiki.cern.ch/twiki/bin/viewauth/CMS/L1ECALPrefiringWeightRecipe#Accessing_the_UL2017_maps, in order to use the current NanoAOD module, extract separate rootfiles from https://github.com/cms-data/PhysicsTools-PatUtils/raw/master/L1PrefiringMaps.root

### data & MC
move "others/for_prefiring/*.root" to NanoAODTools/data/prefire_maps/, and move "others/for_prefiring/PrefireCorr.py" to postprocessing/modules/common/

   
