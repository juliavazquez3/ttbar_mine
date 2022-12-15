import ROOT, os, sys
from ROOT import *
import json
import argparse
import numpy as np

#if not sys.flags.interactive: ROOT.EnableImplicitMT()

# Some defaults
gROOT.SetStyle("Plain")
gStyle.SetOptStat(1111111)
gStyle.SetPadGridX(True)
gStyle.SetPadGridY(True)
gStyle.SetGridStyle(3)
gStyle.SetCanvasDefW(1600)
gStyle.SetCanvasDefH(800)

# Argument parsing
parser = argparse.ArgumentParser()
parser.add_argument("--ssos", action="store_true", default=False,
                    help="Are these ssos plots?")
parser.add_argument("--png", action="store_true", default=False,
                    help="png format")

# Use like:
# python arg.py --data="No"
# python hist_plot.py --data="No" --stack --ratio

args = parser.parse_args()

#if (args.data == "No" or args.data == "2016" or args.data == "2017" or args.data == "2018"): data_op = str(args.data)
#else: raise NameError('Incorrect data option')

plotdir = '/nfs/cms/vazqueze/ttbaranalisis/plotspng/'

if not os.path.exists(plotdir):
    os.makedirs(plotdir)

## Open hists files

if args.ssos: filePath = "/nfs/cms/vazqueze/ttbaranalisis/hists/SF/ssos/"
else: filePath = "/nfs/cms/vazqueze/ttbaranalisis/hists/SF/"

if args.ssos: ssos_add = "SSOS"
else: ssos_add = ""

term = "SFnj_histstt_v1v2"

datayears = ["2016","2016B","2017","2018"]

samplesHT = ["WW_semi_charm","WW_semi_nocharm","WW_hadronic","WW_leptonic","Wjets1_light","Wjets1_bottom","Wjets1_charm","Wjets1_doublecharm","Wjets2_light","Wjets2_bottom","Wjets2_charm","Wjets2_doublecharm",
        "Wjets3_light","Wjets3_bottom","Wjets3_charm","Wjets3_doublecharm","Wjets4_light","Wjets4_bottom","Wjets4_charm","Wjets4_doublecharm",
        "Wjets5_light","Wjets5_bottom","Wjets5_charm","Wjets5_doublecharm","Wjets6_light","Wjets6_bottom","Wjets6_charm","Wjets6_doublecharm",
        "Wjets7_light","Wjets7_bottom","Wjets7_charm","Wjets7_doublecharm","Wjets8_light","Wjets8_bottom","Wjets8_charm","Wjets8_doublecharm",
        "ttbar_charm","ttbar_charmbottom","ttbar_nocharm","ttbarlep_charm","ttbarlep_nocharm","ttbarhad_charm","ttbarhad_nocharm","DY1",
        "DY2","DY3","DY4","DY5","DY6","DY7","DY8","WZ","ZZ","ST1_charm","ST2_charm","ST3_charm","ST4_charm","ST1_nocharm","ST2_nocharm","ST3_nocharm","ST4_nocharm",
        "Wjets0J_charm","Wjets0J_doublecharm","Wjets0J_light","Wjets0J_bottom",
        "Wjets1J_charm","Wjets1J_doublecharm","Wjets1J_light","Wjets1J_bottom","Wjets2J_charm","Wjets2J_doublecharm","Wjets2J_light","Wjets2J_bottom",
        "DY0J","DY1J","DY2J"]

## Adding QCD

histFile = {}

for s in samplesHT:
        ## mc files
        histFile[s] = {}
        for data_op in datayears:
                if s[0:2] == "WW":
                        histFile[s][data_op] = TFile.Open(filePath + term+ssos_add+"_"+s[0:2]+data_op+s[2:]+".root","READ")
                elif (s[0:5] == "Wjets" and s[6]=="_"):
                        histFile[s][data_op] = TFile.Open(filePath + term+ssos_add+"_"+s[0:6]+data_op+s[6:]+".root","READ")
                elif (s[0:5] == "Wjets" and s[6]=="J"):
                        histFile[s][data_op] = TFile.Open(filePath + term+ssos_add+"_"+s[0:7]+data_op+s[7:]+".root","READ")
                elif s[0:6] == "ttbar_":
                        histFile[s][data_op] = TFile.Open(filePath + term+ssos_add+"_"+s[0:5]+data_op+s[5:]+".root","READ")
                elif s[0:8] == "ttbarlep":
                        histFile[s][data_op] = TFile.Open(filePath + term+ssos_add+"_"+s[0:8]+data_op+s[8:]+".root","READ")
                elif s[0:8] == "ttbarhad":
                        histFile[s][data_op] = TFile.Open(filePath + term+ssos_add+"_"+s[0:8]+data_op+s[8:]+".root","READ")
                elif s[0:2] == "ST":
                        histFile[s][data_op] = TFile.Open(filePath + term+ssos_add+"_"+s[0:3]+data_op+s[3:]+".root","READ")
                else:
                        histFile[s][data_op] = TFile.Open(filePath + term+ssos_add+"_"+s+data_op+".root","READ")


samples = ["WW","Wjets1","Wjets2","Wjets3","Wjets4","Wjets5","Wjets6","Wjets7","Wjets8","ttbar","ttbarlep","ttbarhad","DY1",
        "DY2","DY3","DY4","DY5","DY6","DY7","DY8","WZ","ZZ","ST1","ST2","ST3","ST4","QCD","Wjets0J","Wjets1J","Wjets2J","DY0J","DY1J","DY2J"]
samples_d = ["2016","2016B","2017","2018"]
#samples = ["WW","Wjets","ttbar"]

samples_mc = []

for s in samples:
  for y in samples_d:
    samples_mc.append(s+y)

samples_year = []

for y in samples_d:
  samples_year.append(y+"M")
  samples_year.append(y+"E")

## lumi info

lumi = {}
xsecs = {}
nevents = {}

for data_op in datayears:
        files = json.load(open("/nfs/cms/vazqueze/analisisWW/mc_info_v9"+data_op+".json"))
        processes = files.keys()
        lumi[data_op] = {}
        for p in processes:
                num_events = files[p]["events"] # Number of events
                num_files = files[p]["files"] # Number of files
                luminosity = files[p]["lumi"] # Luminosity
                #print(files[p]["type"])
                for s in samples:
                        if files[p]["type"]==s+data_op:
                                #print(s)
                                lumi[data_op][s] = luminosity
#print(lumi)

for data_op in datayears:
        lumi[data_op]["WW_semi_charm"] = lumi[data_op]["WW"]
        lumi[data_op]["WW_semi_nocharm"] = lumi[data_op]["WW"]
        lumi[data_op]["WW_hadronic"] = lumi[data_op]["WW"]
        lumi[data_op]["WW_leptonic"] = lumi[data_op]["WW"]
        lumi[data_op]["ttbar_charm"] = lumi[data_op]["ttbar"]
        lumi[data_op]["ttbar_charmbottom"] = lumi[data_op]["ttbar"]
        lumi[data_op]["ttbar_nocharm"] = lumi[data_op]["ttbar"]
        lumi[data_op]["ttbarlep_charm"] = lumi[data_op]["ttbarlep"]
        lumi[data_op]["ttbarlep_nocharm"] = lumi[data_op]["ttbarlep"]
        lumi[data_op]["ttbarhad_charm"] = lumi[data_op]["ttbarhad"]
        lumi[data_op]["ttbarhad_nocharm"] = lumi[data_op]["ttbarhad"]

        for i in np.arange(8)+1:
                lumi[data_op]["Wjets"+str(i)+"_charm"] = lumi[data_op]["Wjets"+str(i)]
                lumi[data_op]["Wjets"+str(i)+"_doublecharm"] = lumi[data_op]["Wjets"+str(i)]
                lumi[data_op]["Wjets"+str(i)+"_bottom"] = lumi[data_op]["Wjets"+str(i)]
                lumi[data_op]["Wjets"+str(i)+"_light"] = lumi[data_op]["Wjets"+str(i)]

        for i in np.arange(4)+1:
                lumi[data_op]["ST"+str(i)+"_charm"] = lumi[data_op]["ST"+str(i)]
                lumi[data_op]["ST"+str(i)+"_nocharm"] = lumi[data_op]["ST"+str(i)]

## Correcting lumis for Wjets and DY
lumi["2018"]["DY0J"] = 13.55
lumi["2018"]["DY1J"] = 45.56
lumi["2018"]["DY2J"] = 37.95
lumi["2018"]["Wjets0J"] = 2.53
lumi["2018"]["Wjets1J"] = 9.78
lumi["2018"]["Wjets2J"] = 8.78

lumi["2017"]["DY0J"] = 12.27
lumi["2017"]["DY1J"] = 41.64
lumi["2017"]["DY2J"] = 39.62
lumi["2017"]["Wjets0J"] = 2.48
lumi["2017"]["Wjets1J"] = 9.75
lumi["2017"]["Wjets2J"] = 9.06

lumi["2016"]["DY0J"] = 12.5
lumi["2016"]["DY1J"] = 43.57
lumi["2016"]["DY2J"] = 36.65
lumi["2016"]["Wjets0J"] = 2.46
lumi["2016"]["Wjets1J"] = 10.17
lumi["2016"]["Wjets2J"] = 8.81

lumi["2016B"]["DY0J"] = 11.58
lumi["2016B"]["DY1J"] = 40.61
lumi["2016B"]["DY2J"] = 35.70
lumi["2016B"]["Wjets0J"] = 2.60
lumi["2016B"]["Wjets1J"] = 9.81
lumi["2016B"]["Wjets2J"] = 8.60

for data_op in datayears:
        for i in np.arange(3):
                lumi[data_op]["Wjets"+str(i)+"J_charm"] = lumi[data_op]["Wjets"+str(i)+"J"]
                lumi[data_op]["Wjets"+str(i)+"J_doublecharm"] = lumi[data_op]["Wjets"+str(i)+"J"]
                lumi[data_op]["Wjets"+str(i)+"J_bottom"] = lumi[data_op]["Wjets"+str(i)+"J"]
                lumi[data_op]["Wjets"+str(i)+"J_light"] = lumi[data_op]["Wjets"+str(i)+"J"]


#print(lumi)

files_d = json.load(open("/nfs/cms/vazqueze/analisisWW/data_info_v9.json"))
processes_d = files_d.keys()

lumi_d = {}
xsecs_d = {}
nevents_d = {}

for p in processes_d:
    num_events = files_d[p]["events"] # Number of events
    num_files = files_d[p]["files"] # Number of files
    luminosity = files_d[p]["lumi"] # Luminosity
    #print(len(list_files))
    for s in samples_d:
      if files_d[p]["type"]==s+"M":
        lumi_d[s] = luminosity
        nevents_d[s] = num_events

print(lumi_d)
print(nevents_d)

#histNames = [h for h in histNames if (h[-1]=="M")]

#######################################################
########### Start of plot creation ####################
#######################################################

histNames = ["jet_muon_pt_coef_M","jet_muon_pt_coef_E"]

for name in histNames:

  if args.png: c1 = TCanvas("c1","",1200,800)
  else: c1 = TCanvas("c1","",600,400)

  samples = ["ttbar_charm","ttbar_charmbottom","ttbarlep_charm","ttbarhad_charm","ttbar_nocharm","ttbarlep_nocharm","ttbarhad_nocharm"]

  ## HISTS
  histss = {}
  for data_op in datayears:
    histss[data_op] = {}
    for s in samples:
      histss[data_op][s] = histFile[s][data_op].Get(name)

  #samples = ["WW_semi_charm","WW_semi_nocharm","WW_hadronic","WW_leptonic","ttbar_nocharm","ttbar_charm","ttbar_charmbottom","DY","WZ","ZZ","ST_charm","ST_nocharm","Wjets_charm","Wjets_doublecharm","Wjets_bottom","Wjets_light","QCD"]
  samples = ["ttbar_nocharm","ttbar_charm","ttbarlep_charm","ttbar_charmbottom","ttbarlep_nocharm","ttbarhad_charm","ttbarhad_nocharm"]

  histT = {}
  for s in samples:
    histT[s] = histss["2018"][s]
    histT[s].Add(histss["2017"][s])
    histT[s].Add(histss["2016B"][s])
    histT[s].Add(histss["2016"][s])

  if (name=="jet_muon_pt_coef_E" or name=="jet_muon_pt_coef_M"):
    for s in samples:
      print("Number of events for "+name[-1]+" channel in the sample "+s+" is "+str(histT[s].Integral()))

  gStyle.SetOptStat(kFALSE);  ## remove statistics box in histos

  ymax = 0
  for s in samples:
      ## Axis
      histT[s].GetYaxis().SetTitle("coef")
      histT[s].GetXaxis().SetTitle("pt")
      #histT[s].GetXaxis().SetRange(histT[s].GetXaxis().GetFirst(),histT[s].GetXaxis().GetLast()+1)

      y = histT[s].GetMaximum()
      if y>ymax: ymax=y

  histT["ttbar_charm"].Draw()

  ## Plot settings

  term= "profile_ttcharm_v1v2"

  if args.ssos: ssos_add="ssos_"
  else: ssos_add=""

  if args.png: c1.Print(plotdir+ssos_add+term+ name + ".png")
  else: c1.Print(plotdir+ssos_add+term + name + ".pdf")

for s in samplesHT:
        for data_op in datayears:
                        histFile[s][data_op].Close()


