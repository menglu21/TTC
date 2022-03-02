#!/bin/bash

year=${1}
if [ -z "$year" ]; then
  echo "Please specify which year! 2016apv; 2016; 2017; 2018";
  exit 1;
fi

# design to substitute some outdated files for NanoAOD-tools
# Remove some outdated files to save IO for crab jobs
echo Initing

export WORKING_PATH="$CMSSW_BASE/src/PhysicsTools/NanoAODTools/python/postprocessing/"
echo $WORKING_PATH

if [ ${year}==2016apv ]; then
  echo "Initiating setup for 2016apv......";

  echo "Updating pileup"
  cp $WORKING_PATH/analysis/others/for_pileup/mcPileupUL2016.root $WORKING_PATH/data/pileup/
  cp $WORKING_PATH/analysis/others/for_pileup/PileupHistogram-goldenJSON-13tev-UL2016-preVFP-99bins_withVar.root $WORKING_PATH/data/pileup/
  cp $WORKING_PATH/analysis/others/for_pileup/puWeightProducer.py $WORKING_PATH/modules/common/

  echo Updating prefiring correction
  cp $WORKING_PATH/analysis/others/for_prefiring/L1prefiring_jetempt_UL2016BtoH.root $CMSSW_BASE/src/PhysicsTools/NanoAODTools/data/prefire_maps/
  cp $WORKING_PATH/analysis/others/for_prefiring/L1prefiring_jetpt_UL2016BtoH.root $CMSSW_BASE/src/PhysicsTools/NanoAODTools/data/prefire_maps/
  cp $WORKING_PATH/analysis/others/for_prefiring/L1prefiring_photonpt_UL2016BtoH.root $CMSSW_BASE/src/PhysicsTools/NanoAODTools/data/prefire_maps/
  cp $WORKING_PATH/analysis/others/for_prefiring/PrefireCorr.py $WORKING_PATH/modules/common/
  
  echo Updateing JME correction
  echo cleaning unused files
  rm $CMSSW_BASE/src/PhysicsTools/NanoAODTools/data/jme/*.tgz
  cp $WORKING_PATH/analysis/others/for_jme/Summer19UL16APV* $CMSSW_BASE/src/PhysicsTools/NanoAODTools/data/jme/
  cp $WORKING_PATH/analysis/others/for_jme/Summer20UL16APV* $CMSSW_BASE/src/PhysicsTools/NanoAODTools/data/jme/
  cp $WORKING_PATH/analysis/others/for_jme/jetmetHelperRun2.py $WORKING_PATH/modules/jme
  
  echo Updating BJet related
  cp $WORKING_PATH/analysis/others/for_btv/btagSFProducer.py $WORKING_PATH/modules/btv
  cp $WORKING_PATH/analysis/others/for_btv/DeepJet_106XUL16preVFPSF_v1.csv $CMSSW_BASE/src/PhysicsTools/NanoAODTools/data/btagSF/

fi

if [ ${year}==2016 ]; then
  echo "Initiating setup for 2016......";

  echo "Updating pileup"
  cp $WORKING_PATH/analysis/others/for_pileup/mcPileupUL2016.root $WORKING_PATH/data/pileup/
  cp $WORKING_PATH/analysis/others/for_pileup/PileupHistogram-goldenJSON-13tev-UL2016-postVFP-99bins_withVar.root $WORKING_PATH/data/pileup/
  cp $WORKING_PATH/analysis/others/for_pileup/puWeightProducer.py $WORKING_PATH/modules/common/

  echo Updating prefiring correction
  cp $WORKING_PATH/analysis/others/for_prefiring/L1prefiring_jetempt_UL2016BtoH.root $CMSSW_BASE/src/PhysicsTools/NanoAODTools/data/prefire_maps/
  cp $WORKING_PATH/analysis/others/for_prefiring/L1prefiring_jetpt_UL2016BtoH.root $CMSSW_BASE/src/PhysicsTools/NanoAODTools/data/prefire_maps/
  cp $WORKING_PATH/analysis/others/for_prefiring/L1prefiring_photonpt_UL2016BtoH.root $CMSSW_BASE/src/PhysicsTools/NanoAODTools/data/prefire_maps/
  cp $WORKING_PATH/analysis/others/for_prefiring/PrefireCorr.py $WORKING_PATH/modules/common/
  
  echo Updateing JME correction
  echo cleaning unused files
  rm $CMSSW_BASE/src/PhysicsTools/NanoAODTools/data/jme/*.tgz
  cp $WORKING_PATH/analysis/others/for_jme/Summer19UL16_* $CMSSW_BASE/src/PhysicsTools/NanoAODTools/data/jme/
  cp $WORKING_PATH/analysis/others/for_jme/Summer20UL16_* $CMSSW_BASE/src/PhysicsTools/NanoAODTools/data/jme/
  cp $WORKING_PATH/analysis/others/for_jme/jetmetHelperRun2.py $WORKING_PATH/modules/jme
  
  echo Updating BJet related
  cp $WORKING_PATH/analysis/others/for_btv/btagSFProducer.py $WORKING_PATH/modules/btv
  cp $WORKING_PATH/analysis/others/for_btv/DeepJet_106XUL16postVFPSF_v2.csv $CMSSW_BASE/src/PhysicsTools/NanoAODTools/data/btagSF/

fi

if [ ${year}==2017 ]; then
  echo "Initiating setup for 2017......";

  echo "Updating pileup"
  cp $WORKING_PATH/analysis/others/for_pileup/mcPileupUL2017.root $WORKING_PATH/data/pileup/
  cp $WORKING_PATH/analysis/others/for_pileup/PileupHistogram-goldenJSON-13tev-UL2017-99bins_withVar.root $WORKING_PATH/data/pileup/
  cp $WORKING_PATH/analysis/others/for_pileup/puWeightProducer.py $WORKING_PATH/modules/common/

  echo Updating prefiring correction
  cp $WORKING_PATH/analysis/others/for_prefiring/L1prefiring_jetempt_UL2017BtoF.root $CMSSW_BASE/src/PhysicsTools/NanoAODTools/data/prefire_maps/
  cp $WORKING_PATH/analysis/others/for_prefiring/L1prefiring_jetpt_UL2017BtoF.root $CMSSW_BASE/src/PhysicsTools/NanoAODTools/data/prefire_maps/
  cp $WORKING_PATH/analysis/others/for_prefiring/L1prefiring_photonpt_UL2017BtoF.root $CMSSW_BASE/src/PhysicsTools/NanoAODTools/data/prefire_maps/
  cp $WORKING_PATH/analysis/others/for_prefiring/PrefireCorr.py $WORKING_PATH/modules/common/
  
  echo Updateing JME correction
  echo cleaning unused files
  rm $CMSSW_BASE/src/PhysicsTools/NanoAODTools/data/jme/*.tgz
  cp $WORKING_PATH/analysis/others/for_jme/Summer19UL17* $CMSSW_BASE/src/PhysicsTools/NanoAODTools/data/jme/
  cp $WORKING_PATH/analysis/others/for_jme/jetmetHelperRun2.py $WORKING_PATH/modules/jme
  
  echo Updating BJet related
  cp $WORKING_PATH/analysis/others/for_btv/btagSFProducer.py $WORKING_PATH/modules/btv
  cp $WORKING_PATH/analysis/others/for_btv/DeepJet_106XUL17SF_V2p1.csv $CMSSW_BASE/src/PhysicsTools/NanoAODTools/data/btagSF/

fi

if [ ${year}==2018 ]; then
  echo "Initiating setup for 2018......";

  echo "Updating pileup"
  cp $WORKING_PATH/analysis/others/for_pileup/mcPileupUL2018.root $WORKING_PATH/data/pileup/
  cp $WORKING_PATH/analysis/others/for_pileup/PileupHistogram-goldenJSON-13tev-UL2018-99bins_withVar.root $WORKING_PATH/data/pileup/
  cp $WORKING_PATH/analysis/others/for_pileup/puWeightProducer.py $WORKING_PATH/modules/common/

  echo Updateing JME correction
  echo cleaning unused files
  rm $CMSSW_BASE/src/PhysicsTools/NanoAODTools/data/jme/*.tgz
  cp $WORKING_PATH/analysis/others/for_jme/Summer19UL18* $CMSSW_BASE/src/PhysicsTools/NanoAODTools/data/jme/
  cp $WORKING_PATH/analysis/others/for_jme/jetmetHelperRun2.py $WORKING_PATH/modules/jme
  
  echo Updating BJet related
  cp $WORKING_PATH/analysis/others/for_btv/btagSFProducer.py $WORKING_PATH/modules/btv
  cp $WORKING_PATH/analysis/others/for_btv/DeepCSV_106XUL18SF_V1p1.csv $CMSSW_BASE/src/PhysicsTools/NanoAODTools/data/btagSF/

fi

echo Cleaning
echo Cleaning
rm -r $WORKING_PATH/analysis/others/
rm -r $WORKING_PATH/data/roccor*

echo redo scram
cd $CMSSW_BASE/src
scram b

echo "Initing Done \(ᵔᵕᵔ)/"
