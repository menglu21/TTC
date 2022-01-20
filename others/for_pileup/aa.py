import ROOT
from ROOT import TH1D, TFile

f1=TFile.Open('center.root')
f2=TFile.Open('up.root')
f3=TFile.Open('down.root')

h1=TH1D()
h2=TH1D()
h3=TH1D()

f1.GetObject("pileup",h1)
f2.GetObject("pileup",h2)
f3.GetObject("pileup",h3)

h2.SetName('pileup_plus')
h3.SetName('pileup_minus')

fout=TFile.Open('aa.root','RECREATE')
fout.cd()
h1.Write()
h2.Write()
h3.Write()
fout.Close()
