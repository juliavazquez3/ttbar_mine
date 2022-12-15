import sys
import ROOT
import os
from ROOT import *
import json
import argparse
import numpy as np
from os.path import isfile, join, isdir

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
parser.add_argument("--stack", action="store_true", default=False,
                    help="Stack simulation or not")
parser.add_argument("--ssos", action="store_true", default=False,
                    help="Are these ssos plots?")
parser.add_argument("--ratio", action="store_true", default=False,
                    help="Plot ratio or not")
parser.add_argument("--linear", action="store_true", default=False,
                    help="Plot linearly")
parser.add_argument("--png", action="store_true", default=False,
                    help="png format")
parser.add_argument("--qcd", action="store_true", default=False,
                    help="include qcd samples")
parser.add_argument("--nodata", action="store_true", default=False,
                    help="Do not plot data")

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

if args.ssos: filePath = "/nfs/cms/vazqueze/hists_ttbar/hists/fromJF/ssos/"
else: filePath = "/nfs/cms/vazqueze/hists_ttbar/hists/fromJF/"

if args.ssos: ssos_add = "SSOS"
else: ssos_add = "" 

term = "histfromJF_v1v2vBWP"

datayears = ["2016","2016B","2017","2018"]
#datayears = ["2018","2016","2016B"]

samplesHT = ["ww","wjets_2_light","wjets_2_bottom","wjets_1_charm","wjets_1_doublecharm","wjets_2_light","wjets_2_bottom","wjets_2_charm","wjets_2_doublecharm",
        "wjets_3_light","wjets_3_bottom","wjets_3_charm","wjets_3_doublecharm","wjets_4_light","wjets_4_bottom","wjets_4_charm","wjets_4_doublecharm",
        "wjets_5_light","wjets_5_bottom","wjets_5_charm","wjets_5_doublecharm","wjets_6_light","wjets_6_bottom","wjets_6_charm","wjets_6_doublecharm",
        "wjets_7_light","wjets_7_bottom","wjets_7_charm","wjets_7_doublecharm","wjets_8_light","wjets_8_bottom","wjets_8_charm","wjets_8_doublecharm",
        "ttbar_sl_charm","ttbar_sl_nocharm","ttbar_dl_charm","ttbar_dl_nocharm","ttbar_dh_charm","ttbar_dh_nocharm","zjets_1",
        "zjets_2","zjets_3","zjets_4","zjets_5","zjets_6","zjets_7","zjets_8","wz","zz","st_1_charm","st_2_charm","st_3_charm","st_4_charm","st_1_nocharm",
        "st_2_nocharm","st_3_nocharm","st_4_nocharm"]

## Adding QCD

histFile = {}

for data_op in datayears:
	## mc files
	histFile[data_op] = {}
	for s in samplesHT:
		if (s[0:6] == "wjets_") and isfile(filePath + term+ssos_add+"_"+s[0:7]+data_op+s[7:]+".root"):
			histFile[data_op][s] = TFile.Open(filePath + term+ssos_add+"_"+s[0:7]+data_op+s[7:]+".root","READ")
		elif s[0:6] == "ttbar_" and isfile(filePath + term+ssos_add+"_"+s[0:8]+data_op+s[8:]+".root"):
			histFile[data_op][s] = TFile.Open(filePath + term+ssos_add+"_"+s[0:8]+data_op+s[8:]+".root","READ")
		elif s[0:2] == "st" and isfile(filePath + term+ssos_add+"_"+s[0:4]+data_op+s[4:]+".root"):
			histFile[data_op][s] = TFile.Open(filePath + term+ssos_add+"_"+s[0:4]+data_op+s[4:]+".root","READ")
		elif isfile(filePath + term+ssos_add+"_"+s+data_op+".root"):
			histFile[data_op][s] = TFile.Open(filePath + term+ssos_add+"_"+s+data_op+".root","READ")
	#print(histFile[data_op].keys())
#print(histFile.keys())

histFileD = {}

for data_op in datayears:
	histFileD[data_op] = {}
	# data files
	histFileD[data_op]["M"] = TFile.Open(filePath + term+ssos_add+"_"+data_op+"M.root","READ")
	histFileD[data_op]["E"] = TFile.Open(filePath + term+ssos_add+"_"+data_op+"E.root","READ")

histNames = []

histNames.append("nJetGood_M")
histNames.append("nJetGood_E")
histNames.append("nMuoninJet_M")
histNames.append("nMuoninJet_E")
histNames.append("nTau_M")
histNames.append("nTau_E")
histNames.append("second_muon_pt_M")
histNames.append("second_muon_pt_E")
histNames.append("second_electron_pt_M")
histNames.append("second_electron_pt_E")
histNames.append("jet_muon_pt_M")
histNames.append("jet_muon_nmu_M")
histNames.append("jet_muon_mass_M")
histNames.append("jet_not_muon_pt_M")
histNames.append("jet_muon_eta_M")
histNames.append("jet_not_muon_eta_M")
histNames.append("jet_notmuon_qgl_M")
histNames.append("jet_notmuon_nmu_M")
histNames.append("jet_notmuon_mass_M")
histNames.append("jet_muon_pt_E")
histNames.append("jet_muon_nmu_E")
histNames.append("jet_muon_mass_E")
histNames.append("jet_not_muon_pt_E")
histNames.append("jet_muon_eta_E")
histNames.append("jet_not_muon_eta_E")
histNames.append("jet_notmuon_qgl_E")
histNames.append("jet_notmuon_nmu_E")
histNames.append("jet_notmuon_mass_E")
#histNames.append("jet_muon_nmu_good_M")
#histNames.append("jet_muon_nmu_good_E")
histNames.append("lepton_pt_M")
histNames.append("lepton_eta_M")
histNames.append("lepton_pt_E")
histNames.append("lepton_eta_E")
histNames.append("muon_jet_pt_M")
histNames.append("muon_jet_eta_M")
histNames.append("muon_jet_pt_E")
histNames.append("muon_jet_eta_E")
histNames.append("InvM_2jets_M")
histNames.append("InvM_2jets_E")
histNames.append("InvM_jetM_lepM")
histNames.append("InvM_jetM_lepE")
histNames.append("deltaR_jetM_jetNM_M")
histNames.append("deltaR_jetM_jetNM_E")
histNames.append("deltapt_jetM_jetNM_M")
histNames.append("deltapt_jetM_jetNM_E")
histNames.append("deltaeta_jetM_jetNM_M")
histNames.append("deltaeta_jetM_jetNM_E")
histNames.append("deltaphi_jetM_jetNM_M")
histNames.append("deltaphi_jetM_jetNM_E")
histNames.append("MET_pt_M")
histNames.append("MET_pt_E")
histNames.append("transverse_massM")
histNames.append("transverse_massE")
histNames.append("tracks_jetM_M")
histNames.append("tracks_jetNM_M")
histNames.append("tracks_jetM_E")
histNames.append("tracks_jetNM_E")
histNames.append("EMN_jetM_M")
histNames.append("EMC_jetM_M")
histNames.append("EMN_jetM_E")
histNames.append("EMC_jetM_E")
histNames.append("EMtotal_jetM_M")
histNames.append("EMtotal_jetM_E")
histNames.append("InvM_muon_jet_M")
histNames.append("muon_jet_mva_M")
histNames.append("muon_jet_mva_E")
histNames.append("muon_jet_relpt_M")
histNames.append("muon_jet_relpt_E")
histNames.append("muon_jet_sigxy_M")
histNames.append("muon_jet_sigxy_E")
histNames.append("muon_jet_sigz_M")
histNames.append("muon_jet_sigz_E")
histNames.append("muon_jet_xy_M")
histNames.append("muon_jet_xy_E")
histNames.append("muon_jet_z_M")
histNames.append("muon_jet_z_E")
histNames.append("muon_jet_sigr_M")
histNames.append("muon_jet_sigr_E")
histNames.append("muon_jet_iso_M")
histNames.append("muon_jet_iso_E")
histNames.append("my_tau_pt_M")
histNames.append("my_tau_pt_E")
histNames.append("my_tau_mass_M")
histNames.append("my_tau_mass_E")
histNames.append("SSOS_M")
histNames.append("SSOS_E")
histNames.append("pT_sum_M")
histNames.append("pT_sum_E")
histNames.append("pT_product_M")
histNames.append("pT_product_E")
histNames.append("deltaR_lep_2jets_M")
histNames.append("deltaR_lep_2jets_E")
histNames.append("eta_2jets_M")
histNames.append("eta_2jets_E")
histNames.append("pt_2jets_M")
histNames.append("pt_2jets_E")
histNames.append("pt_Wlep_M")
histNames.append("pt_Wlep_E")
histNames.append("pT_proy_M")
histNames.append("pT_proy_E")
histNames.append("deltaphi_MET_2jets_M")
histNames.append("deltaphi_MET_2jets_E")
histNames.append("deltaphi_lephad_M")
histNames.append("deltaphi_lephad_E")
histNames.append("deltaR_lep_jet1_M")
histNames.append("deltaR_lep_jet2_M")
histNames.append("deltaPhi_lep_jet1_M")
histNames.append("deltaPhi_lep_jet2_M")
histNames.append("deltaEta_lep_jet1_M")
histNames.append("deltaEta_lep_jet2_M")
histNames.append("deltaR_lep_jet1_E")
histNames.append("deltaR_lep_jet2_E")
histNames.append("deltaPhi_lep_jet1_E")
histNames.append("deltaPhi_lep_jet2_E")
histNames.append("deltaEta_lep_jet1_E")
histNames.append("deltaEta_lep_jet2_E")
histNames.append("jet_bot1_btag_M")
histNames.append("jet_bot1_btag_E")
histNames.append("jet_bot2_btag_M")
histNames.append("jet_bot2_btag_E")
histNames.append("jet_muon_btag_M")
histNames.append("jet_muon_btag_E")
histNames.append("jet_notmuon_btag_M")
histNames.append("jet_notmuon_btag_E")
histNames.append("jet_notmuon_deeptagG_M")
histNames.append("jet_notmuon_deeptagG_E")
histNames.append("jet_notmuon_deeptagC_M")
histNames.append("jet_notmuon_deeptagC_E")
histNames.append("nMuon_nobot_M")
histNames.append("nMuon_nobot_E")
histNames.append("InvM_bot_farther_M")
histNames.append("InvM_bot_farther_E")
histNames.append("InvM_bot_closer_M")
histNames.append("InvM_bot_closer_E")

#histNames = ["nMuon_nobot_M","nMuon_nobot_E"]
#histNames = ["jet_muon_coef_M","jet_muon_coef_E","jet_muon_ptresol_E","jet_muon_ptresol_M",
#	"jet_muon_scalefactor_M","jet_muon_scalefactor_E"]

histNames = [
       "jet_muon_flavourH_M","jet_muon_flavourH_E","jet_notmuon_flavourH_M","jet_notmuon_flavourH_E",
       "jet_muon_flavourP_M","jet_muon_flavourP_E","jet_notmuon_flavourP_M","jet_notmuon_flavourP_E",
       "btag_sf_M","btag_sf_E","Br_weight_sl_M","Br_weight_sl_E","Frag_weight_sl_M","Frag_weight_sl_E",
       "lep_id_sf_M","lep_id_sf_E","lep_id_lowpt_sf_M","lep_id_lowpt_sf_E","puWeight_M","puWeight_E",
       "PUjetID_SF_M","PUjetID_SF_E"
]

not_rebin = ["nJetGood_M","nJetGood_E","nMuoninJet_M","nMuoninJet_E","jet_muon_nmu_M","jet_muon_nmu_E","SSOS_M","SSOS_E",
       "nLooseLepton_M","nLooseLepton_E", "muon_jet_tight_M","muon_jet_tight_E","nTau_M","nTau_E",
       "jet_muon_flavourH_M","jet_muon_flavourH_E","jet_notmuon_flavourH_M","jet_notmuon_flavourH_E",
       "jet_muon_flavourP_M","jet_muon_flavourP_E","jet_notmuon_flavourP_M","jet_notmuon_flavourP_E","nMuon_nobot_M","nMuon_nobot_E",
       "second_muon_pt_M","second_muon_pt_E","second_electron_pt_M","second_electron_pt_E"]

samples = ["ww","wjets_1","wjets_2","wjets_3","wjets_4","wjets_5","wjets_6","wjets_7","wjets_8","ttbar_sl","ttbar_dl","ttbar_dh","zjets_1",
        "zjets_2","zjets_3","zjets_4","zjets_5","zjets_6","zjets_7","zjets_8","wz","zz","st_1","st_2","st_3","st_4"]

#samples_d = ["2018"]
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

for data_op in samples_d:
	files = json.load(open("/nfs/cms/vazqueze/ttbaranalisis/mcinfo"+data_op+".json"))
	lumi[data_op] = {}
	for p in samples:
		#num_events = files[p]["events"] # Number of events
		num_files = files[p]["files"] # Number of files
		luminosity = files[p]["lumi"] # Luminosity
		#print(files[p]["type"])
		lumi[data_op][p] = luminosity

for data_op in samples_d:
	lumi[data_op]["ttbar_sl_charm"] = lumi[data_op]["ttbar_sl"]
	lumi[data_op]["ttbar_sl_nocharm"] = lumi[data_op]["ttbar_sl"]
	lumi[data_op]["ttbar_dl_charm"] = lumi[data_op]["ttbar_dl"]
	lumi[data_op]["ttbar_dl_nocharm"] = lumi[data_op]["ttbar_dl"]
	lumi[data_op]["ttbar_dh_charm"] = lumi[data_op]["ttbar_dh"]
	lumi[data_op]["ttbar_dh_nocharm"] = lumi[data_op]["ttbar_dh"]

	for i in np.arange(8)+1:
		lumi[data_op]["wjets_"+str(i)+"_charm"] = lumi[data_op]["wjets_"+str(i)]
		lumi[data_op]["wjets_"+str(i)+"_doublecharm"] = lumi[data_op]["wjets_"+str(i)]
		lumi[data_op]["wjets_"+str(i)+"_bottom"] = lumi[data_op]["wjets_"+str(i)]
		lumi[data_op]["wjets_"+str(i)+"_light"] = lumi[data_op]["wjets_"+str(i)]

	for i in np.arange(4)+1:
		lumi[data_op]["st_"+str(i)+"_charm"] = lumi[data_op]["st_"+str(i)]
		lumi[data_op]["st_"+str(i)+"_nocharm"] = lumi[data_op]["st_"+str(i)]

#print("mc lumis",lumi)

files_d = json.load(open("/nfs/cms/vazqueze/ttbaranalisis/datainfo.json"))

lumi_d = {}
xsecs_d = {}
nevents_d = {}

for p in samples_year:
    num_events = files_d[p]["events"] # Number of events
    num_files = files_d[p]["files"] # Number of files
    luminosity = files_d[p]["lumi"] # Luminosity
    #print(len(list_files))
    lumi_d[p[:-1]] = luminosity
    nevents_d[p[:-1]] = num_events

#print("data lumis are",lumi_d)
#print(nevents_d)

#histNames = [h for h in histNames if (h[-1]=="M")]

#######################################################
########### Start of plot creation ####################
#######################################################

for name in histNames:

  if args.png: c1 = TCanvas("c1","",1200,800)
  else: c1 = TCanvas("c1","",600,400)

  if args.ratio and not args.nodata:
    ## In case of ratio plot
    upper_pad = ROOT.TPad("upper_pad", "", 0, 0.35, 1, 1)
    lower_pad = ROOT.TPad("lower_pad", "", 0, 0, 1, 0.35)
    for p in [upper_pad, lower_pad]:
        p.SetLeftMargin(0.14)
        p.SetRightMargin(0.05)
        p.SetTickx(False)
        p.SetTicky(False)
    upper_pad.SetBottomMargin(0)
    lower_pad.SetTopMargin(0)
    lower_pad.SetBottomMargin(0.3)
 
    if not args.linear: upper_pad.SetLogy()
    upper_pad.Draw()
    lower_pad.Draw()

  samples = ["ww","wjets_1_charm","wjets_1_doublecharm","wjets_1_bottom","wjets_1_light","wjets_2_charm","wjets_2_doublecharm","wjets_2_bottom","wjets_2_light",
        "wjets_3_charm","wjets_3_light","wjets_4_charm","wjets_4_light","wjets_5_charm","wjets_5_light",
        "wjets_6_charm","wjets_6_light","wjets_7_charm","wjets_7_light","wjets_8_charm","wjets_8_light",
        "ttbar_sl_charm","ttbar_sl_nocharm","ttbar_dl_charm","ttbar_dl_nocharm","ttbar_dh_charm","ttbar_dh_nocharm",
        "zjets_1","zjets_2","zjets_3","zjets_4","zjets_5","zjets_6","zjets_7","zjets_8","wz","zz","st_1_charm","st_2_charm","st_3_charm","st_4_charm","st_1_nocharm","st_2_nocharm","st_3_nocharm","st_4_nocharm"]

  ## HISTS
  samples_foryear = {}
  histss = {}
  hdataM = {}
  hdataE = {}
  for data_op in datayears:
    samples_foryear[data_op] = [s for s in samples if s in histFile[data_op].keys()]
    histss[data_op] = {}
    for s in samples_foryear[data_op]:
      histss[data_op][s] = histFile[data_op][s].Get(name)
    if not args.nodata: hdataM[data_op] = histFileD[data_op]["M"].Get(name)
    if not args.nodata: hdataE[data_op] = histFileD[data_op]["E"].Get(name)
    if (args.qcd and name[-1] == "M") : histss[data_op]["QCD"] = histFile["QCD"][data_op].Get(name)
  
  samples_stc = {}
  samples_stnc = {}
  samples_wjetsc = {}
  samples_wjetsdc = {}
  samples_wjetsl = {}
  samples_wjetsb = {}
  samples_zjets = {}
  ## Scaling to lumi
  #print(name) 
  for data_op in datayears:
    lumi_data = lumi_d[data_op]
    for s in samples_foryear[data_op]:
      #print(s)
      #print(data_op)
      if s != "QCD": histss[data_op][s].Scale(lumi_data/lumi[data_op][s])
      ## Fixing single top
    samples_stc[data_op] = [s for s in samples_foryear[data_op] if (s[0:2]=="st" and s[-6:]=="_charm")]
    samples_stnc[data_op] = [s for s in samples_foryear[data_op] if (s[0:2]=="st" and s[-6:]=="ocharm")]
    samples_wjetsc[data_op] = [s for s in samples_foryear[data_op] if (s[0:2]=="wj" and s[-6:]=="_charm")] 
    samples_wjetsdc[data_op] = [s for s in samples_foryear[data_op] if (s[0:2]=="wj" and s[-6:]=="echarm")]
    samples_wjetsl[data_op] = [s for s in samples_foryear[data_op] if (s[0:2]=="wj" and s[-6:]=="_light")]
    samples_wjetsb[data_op] = [s for s in samples_foryear[data_op] if (s[0:2]=="wj" and s[-6:]=="bottom")]
    samples_zjets[data_op] = [s for s in samples_foryear[data_op] if s[0:2]=="zj"]
    #print(samples_stc[data_op],samples_stnc[data_op],samples_wjetsc[data_op],samples_wjetsdc[data_op],samples_wjetsl[data_op],samples_wjetsb[data_op],samples_zjets[data_op]) 
    histss[data_op]["st_charm"] = histss[data_op][samples_stc[data_op][0]]
    histss[data_op]["st_nocharm"] = histss[data_op][samples_stnc[data_op][0]]
    histss[data_op]["wjets_charm"] = histss[data_op][samples_wjetsc[data_op][0]]
    histss[data_op]["wjets_doublecharm"] = histss[data_op][samples_wjetsdc[data_op][0]]
    histss[data_op]["wjets_light"] = histss[data_op][samples_wjetsl[data_op][0]]
    histss[data_op]["wjets_bottom"] = histss[data_op][samples_wjetsb[data_op][0]]
    histss[data_op]["zjets"] = histss[data_op][samples_zjets[data_op][0]]
    for s in samples_stc[data_op][1:]:
      histss[data_op]["st_charm"].Add(histss[data_op][s])
    for s in samples_stnc[data_op][1:]:
      histss[data_op]["st_nocharm"].Add(histss[data_op][s])
    for s in samples_wjetsc[data_op][1:]:
      histss[data_op]["wjets_charm"].Add(histss[data_op][s])
    for s in samples_wjetsdc[data_op][1:]:
      histss[data_op]["wjets_doublecharm"].Add(histss[data_op][s])
    for s in samples_wjetsl[data_op][1:]:
      histss[data_op]["wjets_light"].Add(histss[data_op][s])
    for s in samples_wjetsb[data_op][1:]:
      histss[data_op]["wjets_bottom"].Add(histss[data_op][s])
    for s in samples_zjets[data_op][1:]:
      histss[data_op]["zjets"].Add(histss[data_op][s])

  samples = ["ww","ttbar_sl_nocharm","ttbar_sl_charm","ttbar_dl_charm","ttbar_dl_nocharm","ttbar_dh_charm","ttbar_dh_nocharm","zjets","wz","zz",
    "st_charm","st_nocharm","wjets_charm","wjets_doublecharm","wjets_bottom","wjets_light"]
  samples_foryear["2018"] = ["ww","ttbar_sl_nocharm","ttbar_sl_charm","ttbar_dl_charm","ttbar_dl_nocharm","ttbar_dh_charm","ttbar_dh_nocharm","zjets","wz",
    "st_charm","st_nocharm","wjets_charm","wjets_doublecharm","wjets_bottom","wjets_light"]
  samples_foryear["2017"] = ["ww","ttbar_sl_nocharm","ttbar_sl_charm","ttbar_dl_charm","ttbar_dl_nocharm","ttbar_dh_charm","ttbar_dh_nocharm","zjets","wz",
    "st_charm","st_nocharm","wjets_charm","wjets_doublecharm","wjets_bottom","wjets_light"]
  samples_foryear["2016"] = ["ww","ttbar_sl_nocharm","ttbar_sl_charm","ttbar_dl_charm","ttbar_dl_nocharm","ttbar_dh_charm","ttbar_dh_nocharm","zjets","wz",
    "st_charm","st_nocharm","wjets_charm","wjets_doublecharm","wjets_bottom","wjets_light"]
  samples_foryear["2016B"] = ["ww","ttbar_sl_nocharm","ttbar_sl_charm","ttbar_dl_charm","ttbar_dl_nocharm","ttbar_dh_charm","ttbar_dh_nocharm","zjets","wz",
    "st_charm","st_nocharm","wjets_charm","wjets_doublecharm","wjets_bottom","wjets_light"]

  histT = {}
  for s in samples:
    if s in samples_foryear[datayears[0]]:
       histT[s] = histss[datayears[0]][s]
    for d in datayears[1:]:
       if s in samples_foryear[d]: histT[s].Add(histss[d][s])

  #print(histT)
  ###### Rearrange according to years used
  samples = ["ww","ttbar_dl_charm","ttbar_dl_nocharm","ttbar_dh_charm","ttbar_dh_nocharm","zjets","wz",
    "wjets_charm","wjets_doublecharm","wjets_bottom","wjets_light","st_nocharm","ttbar_sl_nocharm","st_charm","ttbar_sl_charm"]

  if not args.nodata:
    histD = {}
    histD["E"] = hdataE[datayears[0]]
    for d in datayears[1:]:
       histD["E"].Add(hdataE[d])
    histD["M"] = hdataM[datayears[0]]
    for d in datayears[1:]:
       histD["M"].Add(hdataM[d])

  if (name=="nJetGood_E" or name=="nJetGood_M"):
    for s in samples:
      print("Number of events for "+name[-1]+" channel in the sample "+s+" is "+str(histT[s].Integral()))
    if (args.qcd and name[-1] == "M"): print("Number of events for "+name[-1]+" channel in the sample QCD is "+str(histT["QCD"].Integral()))

  gStyle.SetOptStat(kFALSE);  ## remove statistics box in histos

  colors = {}
  colors["ww"] = (222,90,106)
  colors["WW_semi_nocharm"] = (246,165,42)
  colors["WW_hadronic"] = (183,2,2)
  colors["WW_leptonic"] = (153,76,0)
  colors["wjets_doublecharm"] = (51,51,255)
  colors["wjets_charm"] = (155,152,204)
  colors["wjets_bottom"] = (255,0,127)
  colors["wjets_light"] = (255,153,204)
  colors["ttbar_sl_charm"] = (204,255,153)
  colors["ttbar_sl_nocharm"] = (120,154,86)
  colors["ttbar_dl_charm"] = (255,255,0)
  colors["ttbar_dl_nocharm"] = (255,153,51)
  colors["ttbar_dh_charm"] = (204,204,0)
  colors["ttbar_dh_nocharm"] = (153,153,0)
  colors["zjets"] = (153,255,255)
  colors["zz"] = (255,255,102)
  colors["wz"] = (153,153,0)
  colors["st_charm"] = (153,51,255)
  colors["st_nocharm"] = (190,153,228)
  colors["QCD"] = (0,153,76)

  if args.stack:
    ymax = 0
    ymin = 0
    for s in samples:
      histT[s].SetLineWidth(1)
      histT[s].SetFillColor(ROOT.TColor.GetColor(*colors[s]))
      ## Axis
      histT[s].GetYaxis().SetTitle("Number of events")
      histT[s].GetXaxis().SetTitle(name)
      if args.ssos:
         if name not in not_rebin: histT[s].Rebin(2)
      #histT[s].GetXaxis().SetRange(histT[s].GetXaxis().GetFirst(),histT[s].GetXaxis().GetLast()+1)

      y = histT[s].GetMaximum()
      ym = histT[s].GetMinimum()
      if y>ymax: ymax=y
      if ymin>ym: ymin=ym
      if name[-1] == "M" and not args.nodata: y = histD["M"].GetMaximum()
      if name[-1] == "E" and not args.nodata: y = histD["E"].GetMaximum()
      if y>ymax: ymax=y

      histT["ww"].SetMinimum(ymin)
      histT["ww"].SetMaximum(5*ymax)
      if args.linear: histT["ww"].SetMaximum(1.3*ymax)
      if args.linear: histT["ww"].SetMinimum(1.3*ymin)

    ## Stack creation

    if args.ratio and not args.nodata: upper_pad.cd()
    stack = ROOT.THStack()
    for s in samples:
      stack.Add(histT[s])
    if (args.qcd and name[-1] == "M"): stack.Add(histT["QCD"])
    y = stack.GetMaximum()
    if y>ymax: ymax=y
    stack.SetMinimum(1.)
    stack.SetMaximum(5*ymax)
    if args.linear: stack.SetMaximum(1.3*ymax)
    stack.Draw("HIST")

    if args.ssos:
      ## auxiliar line
      hAux = histT["ww"].Clone('hAux')
      for s in samples[1:]:
          hAux.Add(histT[s])
      hAux.SetMarkerStyle(20)
      hAux.SetMarkerSize(0.3)
      hAux.SetLineWidth(1)
      hAux.SetLineColor(ROOT.kGreen)
      hAux.Draw("E SAME")
    
  if args.stack and not args.nodata:
    # Draw data
    if name[-1] == "M": data = histD["M"]
    if name[-1] == "E": data = histD["E"]
    data.SetMarkerStyle(20)
    data.SetMarkerSize(0.3)
    data.SetLineWidth(1)
    data.SetLineColor(ROOT.kBlack)
    if args.ssos:
         if name not in not_rebin: data.Rebin(2)
    data.Draw("E SAME")

    if args.ratio and not args.nodata:
      lower_pad.cd()
      ratio = data.Clone("ratio")
      ratio.SetLineColor(kBlack)
      ratio.SetMarkerStyle(21)
      ratio.SetTitle("")
      ratio.SetMinimum(0)
      ratio.SetMaximum(2)
      ratio.GetYaxis().SetTitle("Data/MC")
      ratio.GetXaxis().SetTitle(name)
      ratio.GetXaxis().SetLabelSize(0.08)
      ratio.GetXaxis().SetTitleSize(0.12)
      ratio.GetXaxis().SetTitleOffset(1.0)
      ratio.GetYaxis().SetLabelSize(0.05)
      ratio.GetYaxis().SetTitleSize(0.09)
      ratio.GetYaxis().CenterTitle()
      ratio.GetYaxis().SetTitleOffset(0.5)
      # Set up plot for markers and errors
      ratio.Sumw2()
      ratio.SetStats(0)
      hTotal = histT["ww"].Clone('hTotal')
      for s in samples[1:]:
        hTotal.Add(histT[s])
      if (args.qcd and name[-1] == "M"): hTotal.Add(histT["QCD"])
      ratio.Divide(hTotal)
      ratio.Draw("ep")

  ## Legends
  if args.ratio and not args.nodata: upper_pad.cd()
  leg = TLegend(0.77,0.7,0.89,0.89)
  leg.SetBorderSize(1)
  for s in samples:
    leg.AddEntry(histT[s],s,"f")
  if (args.qcd and name[-1] == "M"): leg.AddEntry(histT["QCD"],"QCD","f")
  if args.stack and not args.nodata: leg.AddEntry(data, "Data" ,"lep")
  leg.Draw()

  ## Plot settings

  ## print image to a gif file
  term= "totalHT_v1v2vBWP_aftercorr"

  if args.ratio: 
    notation = "_ratio_"
    if args.linear:
      notation = "_linratio_"
  else: 
    notation = "_normed_"

  if args.ssos: ssos_add="ssos_"
  else: ssos_add=""

  if args.png: c1.Print(plotdir+ssos_add+term+notation+ name + ".png")
  else: c1.Print(plotdir+ssos_add+term+notation + name + ".pdf")

for s in samplesHT:
        for data_op in datayears:
                        histFile[s][data_op].Close()

for data_op in datayears:
	histFileD[data_op]["M"].Close()
	histFileD[data_op]["E"].Close()

