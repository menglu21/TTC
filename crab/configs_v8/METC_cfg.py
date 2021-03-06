from WMCore.Configuration import Configuration

config = Configuration()

config.section_("General")
config.General.requestName = 'MET_C'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.allowUndistributedCMSSW = True
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'PSet.py'
config.JobType.scriptExe = 'crab_script_dataC.sh'
# hadd nano will not be needed once nano tools are in cmssw
config.JobType.inputFiles = ['crab_script.py', '../scripts/haddnano.py','keep_and_drop.txt','Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt']
config.JobType.sendPythonFolder = True

config.section_("Data")
config.Data.inputDataset = '/MET/Run2017C-UL2017_MiniAODv1_NanoAODv2-v1/NANOAOD'
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 80
config.Data.lumiMask = 'Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt'
config.Data.publication = False
config.Data.outputDatasetTag = 'MET_C'

config.section_("Site")
config.Site.storageSite = "T3_CH_CERNBOX"
