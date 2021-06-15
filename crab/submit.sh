crab submit -c configs/METB_cfg.py 
rm crab_MET_B/inputs/*.tgz
crab submit -c configs/METC_cfg.py 
rm crab_MET_C/inputs/*.tgz
crab submit -c configs/METD_cfg.py 
rm crab_MET_D/inputs/*.tgz
crab submit -c configs/METE_cfg.py 
rm crab_MET_E/inputs/*.tgz
crab submit -c configs/METF_cfg.py 
rm crab_MET_F/inputs/*.tgz
crab submit -c configs/ttjets_cfg.py 
rm crab_ttJets/inputs/*.tgz
