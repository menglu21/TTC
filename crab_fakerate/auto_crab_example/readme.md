# Auto crab help for TTC analysis
The `crab_help.py` is designed to handle crab jobs automatically, which include `prepare` `submit` `kill` `status` `resubmit` and `hadd` working mode. The `hadd` mode is not stable and still in progress currently.

This is just a preliminary example. In order not to package the environment several times, it would be better to submit crab job from the outside of `python/postprocessing` folder.

Since Nanoaod-tools is a compatible framework, one can easily modify the modules in `TTC_postproc.py` for different purposes.

```bash
cd ../
cp -r auto_crab_example $CMSSW_BASE/src/PhysicsTools/NanoAODTools
cd $CMSSW_BASE/src/PhysicsTools/NanoAODTools/auto_crab_example
```

Noticed that the `crab_help.py` is written in python3, hence the `scram b` in CMSSW would leave some error message. Since this crab helper normally would not be included by other codes, you can ignore these errors.

## Input Json
A json file is needed to store the dataset name and related information. It is interpreted like what `input.json` list. You can specify the input json by using `-f` option. e.g.
```bash
python3 crab_help.py -f input.json -m prepare
```

## Several working modes
So several working modes can be activated by using `-m` option as the example above.

### Prepare mode
The `prepare` mode is designed to prepare the crab submit code automatically. It will create folders which contain `*_cfg.py` codes for crab submitting.

Different users need to modify the corresponding parts to make it work properly.

First the store site and the outputdir for crab job should be modified for different users.
e.g.
```python
outLFNDirBase = '/store/user/sdeng/TTC/test/'

f.write('config.Site.storageSite = "T2_CN_Beijing"\n')
```

Second the input files and execution files could be changed for different usage.

### Submit mode
e.g.
```bash
python3 crab_help -f input.json -m submit
```
After the preparation, one can use `submit` mode directly. Crab working directories would be created automatically. You can check these directories manually as how you use crab normally.

### Status mode
e.g.
```bash
python3 crab_help -f input.json -m status
```
After the submitting, one can use `status` mode to check job status.

### Resubmit mode
e.g.
```bash
python3 crab_help -f input.json -m resubmit
```
After the submitting, one can use `resubmit` mode to resubmit jobs.

### Kill mode
e.g.
```bash
python3 crab_help -f input.json -m kill
```
After the submitting, one can use `kill` mode to kill jobs and remove related working directories.