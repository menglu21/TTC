from WMCore.Configuration import Configuration

config = Configuration()

config.section_("General")
config.General.requestName = 'ZZ'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.allowUndistributedCMSSW = True
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'PSet.py'
config.JobType.scriptExe = 'crab_script.sh'
# hadd nano will not be needed once nano tools are in cmssw
config.JobType.inputFiles = ['crab_script.py', '../scripts/haddnano.py','keep_and_drop.txt']
config.JobType.sendPythonFolder = True

config.section_("Data")
config.Data.inputDataset = '/ZZ_TuneCP5_13TeV-pythia8/RunIISummer19UL17NanoAOD-106X_mc2017_realistic_v6-v3/NANOAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
#config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 1
config.Data.totalUnits = -1
config.Data.publication = False
config.Data.outputDatasetTag = 'ZZ'

config.section_("Site")
config.Site.storageSite = "T3_CH_CERNBOX"
