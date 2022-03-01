from WMCore.Configuration import Configuration

config = Configuration()

config.section_("General")
config.General.requestName = 'dummy'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.allowUndistributedCMSSW = True
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'PSet.py'
config.JobType.scriptExe = 'dummy'
# hadd nano will not be needed once nano tools are in cmssw
config.JobType.inputFiles = ['crab_script.py', '../scripts/haddnano.py','keep_and_drop.txt','dummy']
config.JobType.sendPythonFolder = True

config.section_("Data")
config.Data.inputDataset = 'dummy'
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 80
config.Data.lumiMask = 'dummy'
config.Data.outLFNDirBase = '/store/group/phys_top/ExtraYukawa/TTC_version9/'
config.Data.publication = False
config.Data.outputDatasetTag = 'dummy'

config.section_("Site")
config.Site.storageSite = "T2_CH_CERN"
