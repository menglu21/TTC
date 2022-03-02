import os, sys, json
from shutil import copyfile

# DeepJet_106XUL16preVFPSF_v2_old.csv is combined from wp_deepJet_106XUL16preVFP_v2.csv and reshaping_deepJet_106XUL16preVFP_v2.csv
# DeepJet_106XUL16postVFPSF_v3_old.csv is combined from wp_deepJet_106XUL16postVFP_v3.csv and reshaping_deepJet_106XUL16postVFP_v3.csv
# DeepJet_106XUL17SF_V3_old.csv is combined from wp_deepJet_106XUL17_v3.csv and reshaping_deepJet_106XUL17_v3.csv
# DeepJet_106XUL18SF_V2_old.csv is combined from wp_deepJet_106XUL18_v2.csv and reshaping_deepJet_106XUL18_v2.csv

year = sys.argv[1]

if year=='2016apv':
  fin="DeepJet_106XUL16preVFPSF_v2_old.csv"
  fout='DeepJet_106XUL16preVFPSF_v2.csv'
if year=='2016':
  fin="DeepJet_106XUL16postVFPSF_v3_old.csv"
  fout='DeepJet_106XUL16postVFPSF_v3.csv'
if year=='2017':
  fin="DeepJet_106XUL17SF_V3_old.csv"
  fout='DeepJet_106XUL17SF_V3.csv'
if year=='2018':
  fin="DeepJet_106XUL18SF_V2_old.csv"
  fout='DeepJet_106XUL18SF_V2.csv'

copyfile(fin,fout)

with open(fout, 'r') as f:
  line=f.read()
  length=len(line.split('\n'))
  for il in range(2,length):
    fla_temp=line.split('\n')[il-1].split(',')[3]
    if fla_temp=='0':
      os.system("sed -i '%ss/0/2/' %s" %(il, fout))
    if fla_temp=='4':
      os.system("sed -i '%ss/4/1/' %s" %(il, fout))
    if fla_temp=='5':
      os.system("sed -i '%ss/5/0/' %s" %(il, fout))
    wp_temp=line.split('\n')[il-1].split(',')[0]
    if wp_temp=='L':
      os.system("sed -i '%ss/L/0/' %s" %(il, fout))
    if wp_temp=='M':
      os.system("sed -i '%ss/M/1/' %s" %(il, fout))
    if wp_temp=='T':
      os.system("sed -i '%ss/T/2/' %s" %(il, fout))
    if wp_temp=='shape':
      os.system("sed -i '%ss/shape/3/' %s" %(il, fout))
