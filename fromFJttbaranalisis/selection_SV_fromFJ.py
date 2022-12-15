######################################                                        
######################################                                        
#######      GOOD  VERSION     #######
######################################                                     
######################################                                     

print('SV CHANNEL')

## Modified version, meant to be run locally, not sent to condor 

## Selection for noth MC and data samples, created to distinguish between years
## Different histogram files are produced for each situation 
## It includes an option for SSOS substraction

import ROOT, os, sys
from ROOT import *
from os import listdir
from os.path import isfile, join, isdir

import json
import argparse

sys.path.append('/nfs/cms/vazqueze/ttbaranalisis/last_corrections/')

# Some defaults
gROOT.SetStyle("Plain")
gStyle.SetOptStat(1111111)
gStyle.SetPadGridX(True)
gStyle.SetPadGridY(True)
gStyle.SetGridStyle(3)
gStyle.SetCanvasDefW(1600)
gStyle.SetCanvasDefH(800)

ROOT.EnableImplicitMT()

# Argument parsing
parser = argparse.ArgumentParser()
parser.add_argument("--process", type=string, default="WW",
                    help="Select type of process to run")
parser.add_argument("--year", type=string, default="2016",
                    help="Select year of process to run")
parser.add_argument("--type", type=string, default="data",
                    help="Selec type of data to run")
parser.add_argument("--notfull", action="store_true", default=False,
                    help="Use just a range of the sample")
parser.add_argument("--ssos", action="store_true", default=False,
                    help="Perform ssos substraction")
parser.add_argument('-l','--list', nargs='+', help='range of sample to use')
# Use like:
# python arg.py -l 1234 2345 3456 4567

args = parser.parse_args()

if args.process == "allMC": proc = ["ww","wjets_1","wjets_2","wjets_3","wjets_4","wjets_5","wjets_6","wjets_7","wjets_8",
        "zjets_1","zjets_2","zjets_3","zjets_4","zjets_5","zjets_6","zjets_7","zjets_8",
        "ttbar_sl","ttbar_dl","ttbar_dh","zz","wz","st_1","st_2","st_3","st_4"]
#if args.process == "allMC": proc = ["ww","wjets_1","wjets_2","wjets_3","wjets_5","wjets_6","wjets_7","wjets_8",
#        "zjets_1","zjets_2","zjets_3","zjets_4","zjets_5","zjets_6","zjets_7","zjets_8",
#        "ttbar_sl","ttbar_dl","ttbar_dh","zz","wz","st_1","st_2","st_3","st_4"]
elif (args.process == "ww" or args.process == "wjets_1" or args.process == "wjets_2" or args.process == "wjets_3" or args.process == "wjets_4" 
        or args.process == "wjets_5" or args.process == "wjets_6" or args.process == "wjets_7" or args.process == "wjets_8"
        or args.process == "zjets_1" or args.process == "zjets_2" or args.process == "zjets_3" or args.process == "zjets_4"
        or args.process == "zjets_5" or args.process == "zjets_6" or args.process == "zjets_7" or args.process == "zjets_8"
        or args.process == "ttbar_sl" or args.process == "zjets" or args.process == "wz"or args.process == "zz" or args.process == "st_1"
        or args.process == "st_2" or args.process == "st_3" or args.process == "st_4"  or args.process == "zjets"  or args.process == "ttbar_dl"
        or args.process == "ttbar_dh" or args.process == "M" or args.process == "E"): proc = [str(args.process)]
else: raise NameError('Incorrect process name')

if (args.year == "2016" or args.year == "2016B" or args.year == "2017" or args.year == "2018"): year = str(args.year)
else: raise NameError('Incorrect year')

if args.type == "data": mode = "data"
elif args.type == "mc": mode = "mc"
else: raise NameError('Incorrect type')

samples = []

for p in proc:
    if mode == "mc": samples.append(p+year)
    else: samples.append(year+p)

samples_test = samples[:]

print("the samples treated are",samples)

if args.notfull:
	if len(args.list) != 2: raise NameError('List has to have 2 elements')
	print("the range of files is",args.list)
 
if mode == "mc": files = json.load(open("/nfs/cms/vazqueze/ttbaranalisis/mcinfo"+year+".json"))
else: files = json.load(open("/nfs/cms/vazqueze/ttbaranalisis/datainfo.json"))

infos = files.keys()

df = {}
xsecs = {}
sumws = {}
archives = {}
event_test = {}

for s in samples:
  archives[s]=[]

#term1 = "/nfs/cms/vazqueze/cmt/PreprocessRDF/"
#if year == "2018": term1 = "/pnfs/ciemat.es/data/cms/store/user/juvazque/PreprocessRDF/"

term1 =	"/pnfs/ciemat.es/data/cms/store/user/juvazque/PreprocessRDF/"

for p in proc:
    # Construct the dataframes
    folder = term1+"myconfig"+year+"/{sample}/cat_base/prod_test/" # Folder name
    if args.type == "mc":
        folder = folder.format(sample = p) # Sample name
        num_files = files[p]["files"] # Number of files
        list_files = [f for f in listdir(folder) if isfile(join(folder, f))] # Lista de archivos
        if (num_files == len(list_files)):
           for f in list_files:
               file_r = TFile(join(folder, f))
               if file_r.GetListOfKeys().Contains("Events"):
                  archives[p+year].append(join(folder,f))
    else:
        if p == "M":
          folder1 = term1+"myconfig"+year+"/"
          for f in [f1 for f1 in listdir(folder1) if (isdir(join(folder1, f1)) and f1[0:7] == "data_mu")]:
             folder = folder1 + f + "/cat_base/prod_test/"
             #num_files = files[f]["files"] # Number of files
             list_files = [f1 for f1 in listdir(folder) if isfile(join(folder, f1))] # Lista de archivos
             #if (num_files == len(list_files)):
             for fil in list_files:
                 file_r = TFile(join(folder, fil))
                 if file_r.GetListOfKeys().Contains("Events"):
                   archives[year+p].append(join(folder,fil))
       	if p ==	"E":
          folder1 = term1+"myconfig"+year+"/"
       	  for f in [f1 for f1 in listdir(folder1) if (isdir(join(folder1, f1)) and f1[0:7] == "data_el")]:
             folder = folder1 + f + "/cat_base/prod_test/"
             #num_files = files[f]["files"] # Number of files
             list_files = [f1 for f1 in listdir(folder) if isfile(join(folder, f1))] # Lista de archivos
             #if (num_files == len(list_files)):
             for fil in list_files:
                 file_r = TFile(join(folder, fil))
                 if file_r.GetListOfKeys().Contains("Events"):
                   archives[year+p].append(join(folder,fil))

for s in samples:
  if args.notfull: archives[s]=archives[s][int(args.list[0]):int(args.list[1])]
  df[s] = ROOT.RDataFrame("Events",set(archives[s]))
  print("Number of files for",s,len(archives[s]))

samples = [s for s in samples if len(archives[s]) > 0]

df_muon = {}
df_electron = {}

## Cuts per year

cuts_btag = {}

cuts_btag["2016"]=[0.0480, 0.2489, 0.6377]
cuts_btag["2016B"]=[0.0480, 0.2489, 0.6377]
cuts_btag["2017"]=[0.0532, 0.3040, 0.7476]
cuts_btag["2018"]=[0.0490,0.2783,0.7100]

## Muon and Electron pT cuts per year

muon_pt = {}

muon_pt["2016"]=30
muon_pt["2016B"]=30
muon_pt["2017"]=30
muon_pt["2018"]=30

el_pt = {}

el_pt["2016"]=35
el_pt["2016B"]=35
el_pt["2017"]=35
el_pt["2018"]=35

#################################
###### Extra definitions ########
#################################

## Numero de SV dentro de un jet
gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto muoninjet(UInt_t nmu, Vint mu_jetid, Vint mu_good) {
            vector<int> vb;
            for (unsigned int i=0; i<nmu; ++i) {
                if (mu_good.size()>0){
                	if (i!=mu_good[0] && mu_jetid[i]!=-1){
                                vb.push_back(i);
                	}
                } else {
                        if (mu_jetid[i]!=-1){
                                vb.push_back(i);
                        }
                }

            }
            return vb;
      };
""")

#######   JETS   #######

## Funciones para seleccionar JETS
gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto JetInd(UInt_t njet, Vfloat pt, Vfloat eta, Vfloat phi, Vint mu_good, Vfloat mu_eta, Vfloat mu_phi, Vint el_good, Vfloat el_eta, Vfloat el_phi, Vint puID, Vint jetID) {
            vector<int> vb;
	    bool cond = false;
	    bool cond1 = false; 
            for (unsigned int i=0; i<njet; ++i) {
		if(mu_good.size()>0){
			pt[i]<50. ? cond1 = puID[i]>=4 : cond1 = true;
			cond = ROOT::VecOps::DeltaR(mu_eta[mu_good[0]],eta[i],mu_phi[mu_good[0]],phi[i]) > 0.4;
                        if (pt[i]>30. && fabs(eta[i])<2.4 && cond && cond1 && jetID[i]>1){
                                vb.push_back(i);
                        }
		}
                if(el_good.size()>0){
			pt[i]<50. ? cond1 = puID[i]>=4 : cond1 = true;
                        cond = ROOT::VecOps::DeltaR(el_eta[el_good[0]],eta[i],el_phi[el_good[0]],phi[i]) > 0.4;
                        if (pt[i]>30. && fabs(eta[i])<2.4 && cond && cond1 && jetID[i]>1){
                                vb.push_back(i);
                        }
                }
            }
            return vb;
      };

      auto InvariantM(const float pt, const float eta, const float phi, const float mass, const float pt1, const float eta1, const float phi1, const float mass1) {
            auto x = pt*std::cos(phi);
	    auto x1 = pt1*std::cos(phi1);
            auto y = pt*std::sin(phi);
            auto y1 = pt1*std::sin(phi1);
            auto z = pt*std::sinh(eta);
	    auto z1 = pt1*std::sinh(eta1);
            auto e = std::sqrt(x*x+y*y+z*z+mass*mass);
	    auto e1 = std::sqrt(x1*x1+y1*y1+z1*z1+mass1*mass1);

            auto mJet = std::sqrt((e+e1)*(e+e1)-(x+x1)*(x+x1)-(y+y1)*(y+y1)-(z+z1)*(z+z1));
            return mJet;
      };
      auto InvariantM3(const float pt, const float eta, const float phi, const float mass, const float pt1, const float eta1, const float phi1, const float mass1, const float pt2, const float eta2, const float phi2, const float mass2) {
            auto x = pt*std::cos(phi);
            auto x1 = pt1*std::cos(phi1);
            auto x2 = pt2*std::cos(phi2);
            auto y = pt*std::sin(phi);
            auto y1 = pt1*std::sin(phi1);
            auto y2 = pt2*std::sin(phi2);
            auto z = pt*std::sinh(eta);
            auto z1 = pt1*std::sinh(eta1);
            auto z2 = pt2*std::sinh(eta2);
            auto e = std::sqrt(x*x+y*y+z*z+mass*mass);
            auto e1 = std::sqrt(x1*x1+y1*y1+z1*z1+mass1*mass1);
            auto e2 = std::sqrt(x2*x2+y2*y2+z2*z2+mass2*mass2);

            auto mJet = std::sqrt((e+e1+e2)*(e+e1+e2)-(x+x1+x2)*(x+x1+x2)-(y+y1+y2)*(y+y1+y2)-(z+z1+z2)*(z+z1+z2));
            return mJet;
      };
""")

######### bottom jets ordering

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      template <typename T>
      vector<size_t> sort_indexes(const vector<T> &v) {

          // initialize original index locations
          vector<size_t> idx(v.size());
          iota(idx.begin(), idx.end(), 0);

          // sort indexes based on comparing values in v
          // using std::stable_sort instead of std::sort
          // to avoid unnecessary index re-orderings
          // when v contains elements of equal values
          stable_sort(idx.begin(), idx.end(),[&v](size_t i1, size_t i2) {return v[i1] > v[i2];});

          return idx;
      };
      auto bottomjets(UInt_t njetgood, Vint jetgood, UInt_t njet, Vfloat jet_btag) {
          vector<float> vb;
          vector<int> fin;
          if (njetgood>0) {
                for (unsigned int i=0; i<njetgood; ++i){
                    vb.push_back(jet_btag[jetgood[i]]);
                }
          }
          for (auto i: sort_indexes(vb)){
                fin.push_back(i);
          }
          return fin;      
      };
""")

## Muon dentro del jet

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto JetMuonIndJet(UInt_t njet, Vint good, Vfloat pt, Vfloat eta, Vfloat phi, UInt_t nmu, Vint mu_good, Vfloat mu_eta, Vfloat mu_phi, Vfloat mu_pt, Vfloat mu_iso, Vint el_good, Vint mu_jetid, Vint jetbotind) {
            vector<int> vb;
            bool cond = false;
            bool cond1 = true;
            bool condb = true;
            int ind=-1;
            float ptM{-10.};
            float ptJ{-10.};
            if (good.size() == 1){
                for (unsigned int i=0; i<nmu; ++i){
                        cond = ROOT::VecOps::DeltaR(mu_eta[i],eta[good[0]],mu_phi[i],phi[good[0]]) < 0.4;
                        if (mu_good.size() > 0) cond1 = mu_good[0] != i;
                        //if(cond && mu_pt[i] > ptM && cond1 && mu_jetid[i]==good[0] && mu_iso[i] > 0.2){
                        if(cond && mu_pt[i] > ptM && cond1 && mu_jetid[i]==good[0]){
                                ind = good[0];
                                ptM = mu_pt[i];
                        }
                }
            }
            if (good.size() > 1){
                for (unsigned int j=0; j<good.size(); ++j){
                        for (unsigned int i=0; i<nmu; ++i){
                                cond = ROOT::VecOps::DeltaR(mu_eta[i],eta[good[j]],mu_phi[i],phi[good[j]]) < 0.4;
                                if (mu_good.size() > 0) cond1 = mu_good[0] != i;
                                if (good.size() > 2) condb = (good[j] != good[jetbotind[0]] && good[j] != good[jetbotind[1]]);
                                //if(cond && mu_pt[i] > ptM && cond1 && mu_jetid[i]==good[j] && mu_iso[i] > 0.2 && condb){
                                if(cond && mu_pt[i] > ptM && cond1 && mu_jetid[i]==good[j] && condb){
                                        ind = good[j];
                                        ptM = mu_pt[i];
                                        ptJ = pt[good[j]];
                                }
                        }
                }
            }
            if (ind>-1) {
		vb.push_back(ind);
	    }
            return vb;
      };
      auto JetMuonInd(UInt_t njet, Vint good, Vfloat pt, Vfloat eta, Vfloat phi, UInt_t nmu, Vint mu_good, Vfloat mu_eta, Vfloat mu_phi, Vfloat mu_pt, Vfloat mu_iso, Vint el_good, Vint mu_jetid, Vint jetbotind) {
            vector<int> vb;
            bool cond = false;
	    bool cond1 = true;
            bool condb = true;
            int ind=-1;
            float ptM{-10.};
            float ptJ{-10.};
	    if (good.size() == 1){
		for (unsigned int i=0; i<nmu; ++i){
			cond = ROOT::VecOps::DeltaR(mu_eta[i],eta[good[0]],mu_phi[i],phi[good[0]]) < 0.4;
			if (mu_good.size() > 0) cond1 = mu_good[0] != i;
                        //if(cond && mu_pt[i] > ptM && cond1 && mu_jetid[i]==good[0] && mu_iso[i] > 0.2){
			if(cond && mu_pt[i] > ptM && cond1 && mu_jetid[i]==good[0]){
                                ind = i;
				ptM = mu_pt[i];
			}
		}
	    }
	    if (good.size() > 1){
		for (unsigned int j=0; j<good.size(); ++j){
                	for (unsigned int i=0; i<nmu; ++i){
                        	cond = ROOT::VecOps::DeltaR(mu_eta[i],eta[good[j]],mu_phi[i],phi[good[j]]) < 0.4;
				if (mu_good.size() > 0) cond1 = mu_good[0] != i;
                                if (good.size() > 2) condb = (good[j] != good[jetbotind[0]] && good[j] != good[jetbotind[1]]);
                                //if(cond && mu_pt[i] > ptM && cond1 && mu_jetid[i]==good[j] && mu_iso[i] > 0.2 && condb){
                        	if(cond && mu_pt[i] > ptM && cond1 && mu_jetid[i]==good[j] && condb){
                                        ind = i;
                                	ptM = mu_pt[i];
                                        ptJ = pt[good[j]];
                        	}
                	}
		}
	    }
            if (ind>-1) {
		vb.push_back(ind);
            }
	    return vb;
      };
""")

######### pedimos algunas condiciones al muon seleccionado

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto muonCond( Vint mu_jet,Vfloat mu_pt, Vfloat mu_eta, Vfloat mu_iso, Vbool mu_softid) {
	    bool cond = false;
	    if (mu_jet.size()>0){ 
	    	if (mu_pt[mu_jet[0]]<25. && fabs(mu_eta[mu_jet[0]])<2.4 && mu_softid[mu_jet[0]]) {
			cond = true;
	    	}
	    }
            return cond;
      };
""")

## Secondary vertex

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto JetSVIndJet(UInt_t njet, Vint good, Vfloat pt, Vfloat eta, Vfloat phi, UInt_t nSV, Vfloat sv_pt, Vfloat sv_phi, Vfloat sv_eta, Vint jetbotind) {
            vector<int> vb;
            bool cond = false;
            bool condb = true;
            float ptM{-10.};
            float ptJ{-10.};
            int ind = -1;
            for (unsigned int j=0; j<good.size(); ++j){
                        for (unsigned int i=0; i<nSV; ++i){
                                cond = ROOT::VecOps::DeltaR(sv_eta[i],eta[good[j]],sv_phi[i],phi[good[j]]) < 0.4;
                                if (good.size() > 2) condb = (good[j] != good[jetbotind[0]] && good[j] != good[jetbotind[1]]);
                                if(cond && pt[good[j]] > ptJ && condb){
                                        ind = good[j];
                                        ptJ = pt[good[j]];
                                }
                        }
            }
            if (ind>-1) {
		vb.push_back(ind);
            }
            return vb;
      };
      auto JetSVIndSV(UInt_t njet, Vint good, Vint jet_sv, Vfloat pt, Vfloat eta, Vfloat phi, UInt_t nSV, Vfloat sv_pt, Vfloat sv_phi, Vfloat sv_eta) {
            vector<int> vb;
            bool cond = false;
            bool condb = true;
            float ptM{-10.};
            float ptJ{-10.};
            int ind = -1;
            for (unsigned int i=0; i<nSV; ++i){
                  if (jet_sv.size() > 0) cond = ROOT::VecOps::DeltaR(sv_eta[i],eta[jet_sv[0]],sv_phi[i],phi[jet_sv[0]]) < 0.4;
                  if(cond && sv_pt[i] > ptM){
                               ind = i;
                               ptM = sv_pt[i];
                  }
            }
            if (ind>-1) {
                vb.push_back(ind);
            }
            return vb;
      };
""")

######### pedimos algunas condiciones al SV seleccionado

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto svCond( Vint sv_jet,Vfloat sv_pt, Vfloat sv_eta, Vint sv_charge) {
	    bool cond = false;
	    if (sv_jet.size()>0){ 
	    	if (sv_charge[sv_jet[0]] != 0) {
			cond = true;
	    	}
	    }
            return cond;
      };
""")

####### Cantidad de svs total

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto nSVtotal(UInt_t njet, Vint jetbot, Vint jetgood, Vfloat pt, Vfloat eta, Vfloat phi, UInt_t nSV, Vfloat sv_pt, Vfloat sv_phi, Vfloat sv_eta, Vint sv_charge) {
            vector<int> vb;
            bool condb = false;
            bool cond = false;
            int ind = -1;
            for (unsigned int i=0; i<nSV; ++i){
                  ind = -1;
                  for (unsigned int j=0; j<jetgood.size(); ++j){
                             if (jetbot.size() > 1) condb = (jetgood[j] != jetgood[jetbot[0]] && jetgood[j] != jetgood[jetbot[1]]);
                             cond = ROOT::VecOps::DeltaR(sv_eta[i],eta[jetgood[j]],sv_phi[i],phi[jetgood[j]]) < 0.4;
                             if(cond && sv_charge[i] != 0 && condb){
                                     ind = i;
                             }

                  }
                  if(ind>-1){
                             vb.push_back(ind);
                  }

            }
            return vb;
      };
""")


####### Cantidad de svs apropiados

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto JetSVMult(UInt_t njet, Vint jetsv, Vfloat pt, Vfloat eta, Vfloat phi, UInt_t nSV, Vfloat sv_pt, Vfloat sv_phi, Vfloat sv_eta, Vint sv_charge, int n) {
            vector<int> vb;
            bool cond = false;
            int ind = -1;
            for (unsigned int i=0; i<nSV; ++i){
                  if (jetsv.size() > 0) cond = ROOT::VecOps::DeltaR(sv_eta[i],eta[jetsv[n]],sv_phi[i],phi[jetsv[n]]) < 0.4;
                  if(cond && sv_charge[i] != 0){
                          vb.push_back(i);
                  }
            }
            return vb;
      };
""")

######### Segundo jet

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto secondjet( UInt_t njetgood,Vint good, Vint jetsv, Vint jetbotind, Vfloat jet_pt) {
            int ind = -1;
            float pT = -10.;
            for (unsigned int i=0; i<njetgood; ++i){
              if (jet_pt[good[i]] > pT && good[i]!=jetsv[0] && good[i]!=good[jetbotind[0]] && good[i]!=good[jetbotind[1]]) {
                    ind = good[i];
                    pT = jet_pt[good[i]];
              }
            }
            return ind;
      };
""")

######### Masa invariante con los bottom

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto InvMassBot(Vint jetgood, Vint jetbot, Vint jetmuon, Vfloat Jet_pt, Vfloat Jet_eta, Vfloat Jet_phi, UInt_t jet_notmuon, Vfloat Jet_mass) {
            int ind = -1;
            vector<float> vb;
            float dR1 = ROOT::VecOps::DeltaR(Jet_eta[jetmuon[0]],Jet_eta[jetgood[jetbot[0]]],Jet_phi[jetmuon[0]],Jet_phi[jetgood[jetbot[0]]]);
            float dR2 = ROOT::VecOps::DeltaR(Jet_eta[jetmuon[0]],Jet_eta[jetgood[jetbot[1]]],Jet_phi[jetmuon[0]],Jet_phi[jetgood[jetbot[1]]]);
            if (dR1 > dR2){
                   vb.push_back(InvariantM3(Jet_pt[jetgood[jetbot[1]]],Jet_eta[jetgood[jetbot[1]]],Jet_phi[jetgood[jetbot[1]]],Jet_mass[jetgood[jetbot[1]]],Jet_pt[jetmuon[0]],Jet_eta[jetmuon[0]],Jet_phi[jetmuon[0]],Jet_mass[jetmuon[0]],Jet_pt[jet_notmuon],Jet_eta[jet_notmuon],Jet_phi[jet_notmuon],Jet_mass[jet_notmuon]));
                   vb.push_back(InvariantM3(Jet_pt[jetgood[jetbot[0]]],Jet_eta[jetgood[jetbot[0]]],Jet_phi[jetgood[jetbot[0]]],Jet_mass[jetgood[jetbot[0]]],Jet_pt[jetmuon[0]],Jet_eta[jetmuon[0]],Jet_phi[jetmuon[0]],Jet_mass[jetmuon[0]],Jet_pt[jet_notmuon],Jet_eta[jet_notmuon],Jet_phi[jet_notmuon],Jet_mass[jet_notmuon]));
            }else{
                   vb.push_back(InvariantM3(Jet_pt[jetgood[jetbot[0]]],Jet_eta[jetgood[jetbot[0]]],Jet_phi[jetgood[jetbot[0]]],Jet_mass[jetgood[jetbot[0]]],Jet_pt[jetmuon[0]],Jet_eta[jetmuon[0]],Jet_phi[jetmuon[0]],Jet_mass[jetmuon[0]],Jet_pt[jet_notmuon],Jet_eta[jet_notmuon],Jet_phi[jet_notmuon],Jet_mass[jet_notmuon]));
                   vb.push_back(InvariantM3(Jet_pt[jetgood[jetbot[1]]],Jet_eta[jetgood[jetbot[1]]],Jet_phi[jetgood[jetbot[1]]],Jet_mass[jetgood[jetbot[1]]],Jet_pt[jetmuon[0]],Jet_eta[jetmuon[0]],Jet_phi[jetmuon[0]],Jet_mass[jetmuon[0]],Jet_pt[jet_notmuon],Jet_eta[jet_notmuon],Jet_phi[jet_notmuon],Jet_mass[jet_notmuon]));
            }
            return vb;
      };
""")

####   SSOS   ####

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto SSOS(Vint mu_good, Vint el_good, Vint sv_jet,Vint el_charge, Vint mu_charge, Vint sv_charge) {
	    int ind = 0;
	    if(sv_jet.size()>0){
            	if(mu_good.size()>0){
			ind = mu_charge[mu_good[0]]*sv_charge[sv_jet[0]];
            	}
            	if(el_good.size()>0){
                	ind = el_charge[el_good[0]]*sv_charge[sv_jet[0]];
            	}
	    }
            return ind;
      };
""")

## pT component calculations

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto pTsum(Vint mu_good, Vint el_good, Vint jet_muon, UInt_t jet_notmuon, Vfloat mu_pt, Vfloat mu_eta, Vfloat mu_phi, Vfloat mu_mass,Vfloat el_pt, Vfloat el_eta, Vfloat el_phi, Vfloat el_mass, Float_t met_pt, Float_t met_phi, Vfloat jet_pt, Vfloat jet_eta, Vfloat jet_phi, Vfloat jet_mass) {
            ROOT::Math::PtEtaPhiMVector plep1;
            if(mu_good.size()>0){
                  plep1.SetPt(mu_pt[mu_good[0]]); plep1.SetEta(mu_eta[mu_good[0]]); plep1.SetPhi(mu_phi[mu_good[0]]); plep1.SetM(mu_mass[mu_good[0]]);
            } else if(el_good.size()>0){
                  plep1.SetPt(el_pt[el_good[0]]); plep1.SetEta(el_eta[el_good[0]]); plep1.SetPhi(el_phi[el_good[0]]); plep1.SetM(el_mass[el_good[0]]);
            }
            ROOT::Math::PtEtaPhiMVector plep2(met_pt, 0., met_phi, 0.);
            ROOT::Math::PtEtaPhiMVector phad1(jet_pt[jet_muon[0]], jet_eta[jet_muon[0]], jet_phi[jet_muon[0]], jet_mass[jet_muon[0]]);
            ROOT::Math::PtEtaPhiMVector phad2(jet_pt[jet_notmuon], jet_eta[jet_notmuon], jet_phi[jet_notmuon], jet_mass[jet_notmuon]);
            auto plep = plep1+plep2;
            auto phad = phad1+phad2;
            auto ptot = plep1+plep2+phad1+phad2;
            auto plepX = plep.Px(); auto plepY = plep.Py(); auto phadX = phad.Px(); auto phadY = phad.Py(); auto ptotX = ptot.Px(); auto ptotY = ptot.Py();
            float plepmod; float phadmod; float ptotmod;
            plepmod = std::sqrt(plepX*plepX + plepY*plepY);
            phadmod = std::sqrt(phadX*phadX + phadY*phadY);
            ptotmod = std::sqrt(ptotX*ptotX + ptotY*ptotY);
            return ptotmod/(plepmod+phadmod);
      };
      auto pTprod(Vint mu_good, Vint el_good, Vint jet_muon, UInt_t jet_notmuon, Vfloat mu_pt, Vfloat mu_eta, Vfloat mu_phi, Vfloat mu_mass,Vfloat el_pt, Vfloat el_eta, Vfloat el_phi, Vfloat el_mass, Float_t met_pt, Float_t met_phi, Vfloat jet_pt, Vfloat jet_eta, Vfloat jet_phi, Vfloat jet_mass) {
            ROOT::Math::PtEtaPhiMVector plep1;
            if(mu_good.size()>0){
                  plep1.SetPt(mu_pt[mu_good[0]]); plep1.SetEta(mu_eta[mu_good[0]]); plep1.SetPhi(mu_phi[mu_good[0]]); plep1.SetM(mu_mass[mu_good[0]]);
            } else if(el_good.size()>0){
                  plep1.SetPt(el_pt[el_good[0]]); plep1.SetEta(el_eta[el_good[0]]); plep1.SetPhi(el_phi[el_good[0]]); plep1.SetM(el_mass[el_good[0]]);
            }
            ROOT::Math::PtEtaPhiMVector plep2(met_pt, 0., met_phi, 0.);
            ROOT::Math::PtEtaPhiMVector phad1(jet_pt[jet_muon[0]], jet_eta[jet_muon[0]], jet_phi[jet_muon[0]], jet_mass[jet_muon[0]]);
            ROOT::Math::PtEtaPhiMVector phad2(jet_pt[jet_notmuon], jet_eta[jet_notmuon], jet_phi[jet_notmuon], jet_mass[jet_notmuon]);
            auto plep = plep1+plep2;
            auto phad = phad1+phad2;
            auto plepX = plep.Px(); auto plepY = plep.Py(); auto phadX = phad.Px(); auto phadY = phad.Py();
            float plepmod; float phadmod; float ptotmod;
            plepmod = std::sqrt(plepX*plepX + plepY*plepY);
            phadmod = std::sqrt(phadX*phadX + phadY*phadY);
            ptotmod = plepX*phadX + plepY*phadY;
            return ptotmod/(plepmod*phadmod);
      };
      auto variousSUM(Vint mu_good, Vint el_good, Vint jet_muon, UInt_t jet_notmuon, Vfloat mu_pt, Vfloat mu_eta, Vfloat mu_phi, Vfloat mu_mass,Vfloat el_pt, Vfloat el_eta, Vfloat el_phi, Vfloat el_mass, Float_t met_pt, Float_t met_phi, Vfloat jet_pt, Vfloat jet_eta, Vfloat jet_phi, Vfloat jet_mass) {
            vector<float> vb;
            ROOT::Math::PtEtaPhiMVector plep1;
            float deltaRlepjet1;
            float deltaRlepjet2;
            float deltaPhilepjet1;
            float deltaPhilepjet2;
            float deltaEtalepjet1;
            float deltaEtalepjet2;
            if(mu_good.size()>0){
                  plep1.SetPt(mu_pt[mu_good[0]]); plep1.SetEta(mu_eta[mu_good[0]]); plep1.SetPhi(mu_phi[mu_good[0]]); plep1.SetM(mu_mass[mu_good[0]]);
                  deltaRlepjet1 = ROOT::VecOps::DeltaR(mu_eta[mu_good[0]],jet_eta[jet_muon[0]],mu_phi[mu_good[0]],jet_phi[jet_muon[0]]);
                  deltaRlepjet2 = ROOT::VecOps::DeltaR(mu_eta[mu_good[0]],jet_eta[jet_notmuon],mu_phi[mu_good[0]],jet_phi[jet_notmuon]);
                  deltaPhilepjet1 = fabs(mu_phi[mu_good[0]]-jet_phi[jet_muon[0]]);
                  deltaPhilepjet2 = fabs(mu_phi[mu_good[0]]-jet_phi[jet_notmuon]);
                  deltaEtalepjet1 = fabs(mu_eta[mu_good[0]]-jet_eta[jet_muon[0]]);
                  deltaEtalepjet2 = fabs(mu_eta[mu_good[0]]-jet_eta[jet_notmuon]);
            } else if(el_good.size()>0){
                  plep1.SetPt(el_pt[el_good[0]]); plep1.SetEta(el_eta[el_good[0]]); plep1.SetPhi(el_phi[el_good[0]]); plep1.SetM(el_mass[el_good[0]]);
                  deltaRlepjet1 = ROOT::VecOps::DeltaR(el_eta[el_good[0]],jet_eta[jet_muon[0]],el_phi[el_good[0]],jet_phi[jet_muon[0]]);
                  deltaRlepjet2 = ROOT::VecOps::DeltaR(el_eta[el_good[0]],jet_eta[jet_notmuon],el_phi[el_good[0]],jet_phi[jet_notmuon]);
                  deltaPhilepjet1 = fabs(el_phi[el_good[0]]-jet_phi[jet_muon[0]]);
                  deltaPhilepjet2 = fabs(el_phi[el_good[0]]-jet_phi[jet_notmuon]);
                  deltaEtalepjet1 = fabs(el_eta[el_good[0]]-jet_eta[jet_muon[0]]);
                  deltaEtalepjet2 = fabs(el_eta[el_good[0]]-jet_eta[jet_notmuon]);
            }
            ROOT::Math::PtEtaPhiMVector plep2(met_pt, 0., met_phi, 0.);
            ROOT::Math::PtEtaPhiMVector phad1(jet_pt[jet_muon[0]], jet_eta[jet_muon[0]], jet_phi[jet_muon[0]], jet_mass[jet_muon[0]]);
            ROOT::Math::PtEtaPhiMVector phad2(jet_pt[jet_notmuon], jet_eta[jet_notmuon], jet_phi[jet_notmuon], jet_mass[jet_notmuon]);
            auto plep = plep1+plep2;
            auto phad = phad1+phad2;
            vb.push_back(ROOT::VecOps::DeltaR(plep1.Eta(),phad.Eta(),plep1.Phi(),phad.Phi()));
            vb.push_back(fabs(plep2.Phi()-phad.Phi()));
            vb.push_back(phad.Eta());
            vb.push_back(phad.Pt());
            vb.push_back(fabs(plep.Phi()-phad.Phi()));
            vb.push_back(ROOT::VecOps::DeltaR(plep.Eta(),phad.Eta(),plep.Phi(),phad.Phi()));
            vb.push_back(fabs(plep1.Phi()-phad.Phi()));
            vb.push_back(fabs(plep.Eta()-phad.Eta()));
            vb.push_back(fabs(plep1.Eta()-phad.Eta()));
            vb.push_back(fabs(plep.Pt()-phad.Pt()));
            vb.push_back(fabs(plep1.Pt()-phad.Pt()));
            vb.push_back(phad2.Pt()*std::sin(phad1.Phi()-phad2.Phi()));
            vb.push_back(phad.Pt()/(phad1.Pt()+phad2.Pt()));
            vb.push_back(plep.Pt());
            vb.push_back(deltaRlepjet1);
            vb.push_back(deltaRlepjet2);
            vb.push_back(deltaPhilepjet1);
            vb.push_back(deltaPhilepjet2);
            vb.push_back(deltaEtalepjet1);
            vb.push_back(deltaEtalepjet2);
            return vb;
      };
""")

## Trigger function for 2017 electron data

gInterpreter.Declare("""
      #include <iomanip>
      #include <math.h>
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto calDeltaR(float eta1, float phi1, float eta2, float phi2) {
          float dphi = phi1-phi2;
          float deta = eta1-eta2;
          if (dphi <= - M_PI) {
             dphi = dphi + M_PI;
          }
          if (dphi >= M_PI) {
             dphi = dphi - M_PI;
          }
          float prod = deta*deta + dphi*dphi;
          return prod;
      };
      auto triggeremulator(UInt_t nel, UInt_t ntrig, Vfloat el_eta, Vfloat el_deltaeta, Vfloat el_phi, Vfloat trig_eta, Vfloat trig_phi, Vint trig_bits, Vint trig_id) {
          bool cond = false;
          bool cond1 = false;
          bool cond2 = false;
          for (unsigned int j=0; j<nel; ++j) {
            for (unsigned int i=0; i<ntrig; ++i) {
                cond1 = calDeltaR(el_eta[j]+el_deltaeta[j],el_phi[j], trig_eta[i], trig_phi[i]) < 0.01;
                cond2 = trig_bits[i] & (0x1 << 10);
                if( trig_id[i]==11 && cond1 && cond2) {
                    cond = true;
                }
            }
          }
          return cond;
      };
""")

## Funciones para seleccionar muones y electrones
gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto muonSecInd(UInt_t nmu, Vfloat pt, Vfloat eta, Vfloat iso, Vbool tID, Float_t cutpt) {
            vector<int> vb;
            float ptest = -10.;
            int ind= -1;
            for (unsigned int i=0; i<nmu; ++i) {
                if (fabs(eta[i])<2.4 && iso[i]<0.15 && tID[i] && pt[i] < cutpt && pt[i]>ptest){
                        ptest = pt[i];
                        ind = i;
                }
            }
            if (ind > -1) vb.push_back(ind);
            return vb;
      };
      auto elSecInd(UInt_t nmu, Vfloat pt, Vfloat eta, Vfloat iso, Vint cutB, Vbool mva80 , Vbool mva90, Float_t cutpt) {
            vector<int> vb;
            float ptest = -10.;
            int ind= -1;
            for (unsigned int i=0; i<nmu; ++i) {
                //if (pt[i]>cutpt && fabs(eta[i])<2.5 && iso[i]<0.15 && mva80[i]){
                if (fabs(eta[i])<2.5 && mva80[i] && pt[i] < cutpt && pt[i]>ptest){
                        ptest = pt[i];
                        ind = i;
                }
            }
            if (ind > -1) vb.push_back(ind);
            return vb;
      };
""")

######### Taus study

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto tau_best(Vint jetmuon, Vfloat tau_pt, Vint tau_idx) {
            int ind = -1;
            float pT = -10.;
            for (unsigned int i=0; i<tau_pt.size(); ++i){
              if (tau_pt[i] > pT && jetmuon[0] == tau_idx[i]) {
                    ind = i;
                    pT = tau_pt[i];
              }
            }
            return ind;
      };
""")

###############################################################################################################################################################################################################

if mode == "data":
        for s in samples:
           df[s] = df[s].Define('Jet_pt_nom','Jet_pt')
           df[s] = df[s].Define('Jet_mass_nom','Jet_mass')
           df[s] = df[s].Define('MET_smeared_phi','MET_phi')
           df[s] = df[s].Define('MET_smeared_pt','MET_pt')

for s in samples:
	df[s] = df[s].Define('weightSSOS','-1*SVLepSign/std::abs(SVLepSign)')
	df[s] = df[s].Filter('nMuonGood<2 && nElectronGood<2').Filter('!(nMuonGood==1) != !(nElectronGood==1)').Filter('nJetGood >= 4').Filter('SVJetInd.size() >= 1').Filter('SVJetGood')
	df[s] = df[s].Filter('!(MuonJetInd.size() >= 1 && MuonJetGood)')
	df[s] = df[s].Define('JetnotSVInd','secondjet(nJetGood, JetGoodInd, JetSVInd, JetBotInd, Jet_pt_nom)')
	df[s] = df[s].Define('SVmaxpTInd','std::distance(SV_pt.begin(),std::max_element(SV_pt.begin(),SV_pt.end()))')
	### hists definitions
	df[s] = df[s].Define('jet_muon_pt','Jet_pt_nom[JetSVInd[0]]')
	df[s] = df[s].Define('jet_muon_nmu','Jet_nMuons[JetSVInd[0]]')
	df[s] = df[s].Define('jet_muon_mass','Jet_mass_nom[JetSVInd[0]]')
	df[s] = df[s].Define('jet_notmuon_pt','nJetGood>1 ? Jet_pt_nom[JetnotSVInd] : 0')
	df[s] = df[s].Define('jet_muon_eta','Jet_eta[JetSVInd[0]]')
	df[s] = df[s].Define('jet_notmuon_eta','nJetGood>1 ? Jet_eta[JetnotSVInd] : 0')
	df[s] = df[s].Define('jet_notmuon_mass','nJetGood>1 ? Jet_mass_nom[JetnotSVInd] : 0')
	df[s] = df[s].Define('jet_notmuon_qgl','nJetGood>1 ? Jet_qgl[JetnotSVInd] : 0')
	df[s] = df[s].Define('jet_notmuon_nmu','nJetGood>1 ? Jet_nMuons[JetnotSVInd] : 0')
	df[s] = df[s].Define('sv_jet_pt','SV_pt[SVJetInd[0]]')
	df[s] = df[s].Define('sv_jet_eta','SV_eta[SVJetInd[0]]')
	df[s] = df[s].Define('sv_jet_ntracks','SV_ntracks[SVJetInd[0]]')
	df[s] = df[s].Define('sv_jet_mass','SV_mass[SVJetInd[0]]')
	df[s] = df[s].Define('sv_jet_pangle','SV_pAngle[SVJetInd[0]]')
	df[s] = df[s].Define('sv_jet_chi','SV_chi2[SVJetInd[0]]')
	df[s] = df[s].Define('sv_jet_ndof','SV_ndof[SVJetInd[0]]')
	df[s] = df[s].Define('sv_jet_charge','SV_charge[SVJetInd[0]]')
	df[s] = df[s].Define('sv_jet_relpt','SV_pt[SVJetInd[0]]/Jet_pt[JetSVInd[0]]')
	df[s] = df[s].Define('sv_max_pt','SV_pt[SVmaxpTInd]')
	df[s] = df[s].Define('sv_max_eta','SV_eta[SVmaxpTInd]')
	df[s] = df[s].Define('sv_max_ntracks','SV_ntracks[SVmaxpTInd]')
	df[s] = df[s].Define('sv_max_mass','SV_mass[SVmaxpTInd]')
	df[s] = df[s].Define('sv_max_pangle','SV_pAngle[SVmaxpTInd]')
	df[s] = df[s].Define('sv_max_chi','SV_chi2[SVmaxpTInd]')
	df[s] = df[s].Define('sv_max_ndof','SV_ndof[SVmaxpTInd]')
	df[s] = df[s].Define('sv_max_charge','SV_charge[SVmaxpTInd]')
	df[s] = df[s].Define('InvM_2jets','nJetGood>1 ? InvariantM(Jet_pt_nom[JetSVInd[0]],Jet_eta[JetSVInd[0]],Jet_phi[JetSVInd[0]],Jet_mass_nom[JetSVInd[0]],Jet_pt_nom[JetnotSVInd],Jet_eta[JetnotSVInd],Jet_phi[JetnotSVInd],Jet_mass_nom[JetnotSVInd]) : 0')
	df[s] = df[s].Define('deltaR_jetM_jetNM','nJetGood>1 ? ROOT::VecOps::DeltaR(Jet_eta[JetnotSVInd], Jet_eta[JetSVInd[0]] , Jet_phi[JetnotSVInd], Jet_phi[JetSVInd[0]])  : 10')
	df[s] = df[s].Define('deltaphi_jetM_jetNM','fabs(Jet_phi[JetSVInd[0]]-Jet_phi[JetnotSVInd])')
	df[s] = df[s].Define('deltaeta_jetM_jetNM','fabs(Jet_eta[JetSVInd[0]]-Jet_eta[JetnotSVInd])')
	df[s] = df[s].Define('deltapt_jetM_jetNM','fabs(Jet_pt_nom[JetSVInd[0]]-Jet_pt_nom[JetnotSVInd])')
	df[s] = df[s].Define('tracks_jetM','Jet_nConstituents[JetSVInd[0]]')
	df[s] = df[s].Define('tracks_jetNM','nJetGood>1 ? Jet_nConstituents[JetnotSVInd] : 0')
	df[s] = df[s].Define('EMN_jetM','Jet_neEmEF[JetSVInd[0]]')
	df[s] = df[s].Define('EMC_jetM','Jet_chEmEF[JetSVInd[0]]')
	df[s] = df[s].Define('EMtotal_jetM','Jet_chEmEF[JetSVInd[0]]+Jet_neEmEF[JetSVInd[0]]')
	df[s] = df[s].Define('MuoninJetAux','muoninjet(nMuon, Muon_jetIdx, MuonGoodInd)')
	df[s] = df[s].Define('nMuoninJet','MuoninJetAux.size()')
	df[s] = df[s].Define('sv_jet_xy','SV_dxy[SVJetInd[0]]')
	df[s] = df[s].Define('sv_jet_sigxy','SV_dxySig[SVJetInd[0]]')
	df[s] = df[s].Define('sv_jet_r','SV_dlen[SVJetInd[0]]')
	df[s] = df[s].Define('sv_jet_sigr','SV_dlenSig[SVJetInd[0]]')
	df[s] = df[s].Define('sv_max_xy','SV_dxy[SVmaxpTInd]')
	df[s] = df[s].Define('sv_max_sigxy','SV_dxySig[SVmaxpTInd]')
	df[s] = df[s].Define('sv_max_r','SV_dlen[SVmaxpTInd]')
	df[s] = df[s].Define('sv_max_sigr','SV_dlenSig[SVmaxpTInd]')
	df[s] = df[s].Define('pT_sum','pTsum(MuonGoodInd, ElectronGoodInd, JetSVInd, JetnotSVInd, Muon_pt, Muon_eta, Muon_phi, Muon_mass, Electron_pt, Electron_eta, Electron_phi, Electron_mass, MET_smeared_pt, MET_smeared_phi, Jet_pt_nom, Jet_eta, Jet_phi, Jet_mass_nom)')
	df[s] = df[s].Define('pT_product','pTprod(MuonGoodInd, ElectronGoodInd, JetSVInd, JetnotSVInd, Muon_pt, Muon_eta, Muon_phi, Muon_mass, Electron_pt, Electron_eta, Electron_phi, Electron_mass, MET_smeared_pt, MET_smeared_phi, Jet_pt_nom, Jet_eta, Jet_phi, Jet_mass_nom)')
	df[s] = df[s].Define('aux_various','variousSUM(MuonGoodInd, ElectronGoodInd, JetSVInd, JetnotSVInd, Muon_pt, Muon_eta, Muon_phi, Muon_mass, Electron_pt, Electron_eta, Electron_phi, Electron_mass, MET_smeared_pt, MET_smeared_phi, Jet_pt_nom, Jet_eta, Jet_phi, Jet_mass_nom)')
	df[s] = df[s].Define('deltaR_lep_2jets','aux_various[0]')
	df[s] = df[s].Define('deltaphi_MET_2jets','aux_various[1]')
	df[s] = df[s].Define('eta_2jets','aux_various[2]')
	df[s] = df[s].Define('pt_2jets','aux_various[3]')
	df[s] = df[s].Define('deltaphi_lephad','aux_various[4]')
	df[s] = df[s].Define('deltaR_lephad','aux_various[5]')
	df[s] = df[s].Define('deltaphi_lep_2jets','aux_various[6]')
	df[s] = df[s].Define('deltaeta_lephad','aux_various[7]')
	df[s] = df[s].Define('deltaeta_lep_2jets','aux_various[8]')
	df[s] = df[s].Define('deltapt_lephad','aux_various[9]')
	df[s] = df[s].Define('deltapt_lep_2jets','aux_various[10]')
	df[s] = df[s].Define('pT_proy','aux_various[11]')
	df[s] = df[s].Define('pT_sum_2J','aux_various[12]')
	df[s] = df[s].Define('pT_Wlep','aux_various[13]')
	df[s] = df[s].Define('deltaR_lep_jet1','aux_various[14]')
	df[s] = df[s].Define('deltaR_lep_jet2','aux_various[15]')
	df[s] = df[s].Define('deltaPhi_lep_jet1','aux_various[16]')
	df[s] = df[s].Define('deltaPhi_lep_jet2','aux_various[17]')
	df[s] = df[s].Define('deltaEta_lep_jet1','aux_various[18]')
	df[s] = df[s].Define('deltaEta_lep_jet2','aux_various[19]')
	df[s] = df[s].Define('jet_bot1_btag','Jet_btagDeepFlavB[JetGoodInd[JetBotInd[0]]]')
	df[s] = df[s].Define('jet_bot2_btag','Jet_btagDeepFlavB[JetGoodInd[JetBotInd[1]]]')
	df[s] = df[s].Define('jet_muon_btag','Jet_btagDeepFlavB[JetSVInd[0]]')
	df[s] = df[s].Define('jet_notmuon_btag','Jet_btagDeepFlavB[JetnotSVInd]')
	df[s] = df[s].Define('nLooseLepton','nMuon+nElectron-1')
	df[s] = df[s].Define('jet_bot1_nsv_good_aux','JetSVMult(nJet, JetBotInd, Jet_pt_nom, Jet_eta, Jet_phi, nSV, SV_pt, SV_phi, SV_eta, SV_charge, 0)')
	df[s] = df[s].Define('jet_bot2_nsv_good_aux','JetSVMult(nJet, JetBotInd, Jet_pt_nom, Jet_eta, Jet_phi, nSV, SV_pt, SV_phi, SV_eta, SV_charge, 1)')
	df[s] = df[s].Define('jet_bot1_nsv_good','jet_bot1_nsv_good_aux.size()')
	df[s] = df[s].Define('jet_bot2_nsv_good','jet_bot2_nsv_good_aux.size()')
	df[s] = df[s].Define('jet_muon_nsv_good_aux','JetSVMult(nJet, JetSVInd, Jet_pt_nom, Jet_eta, Jet_phi, nSV, SV_pt, SV_phi, SV_eta, SV_charge, 0)')
	df[s] = df[s].Define('jet_muon_nsv_good','jet_muon_nsv_good_aux.size()')
	df[s] = df[s].Define('nSV_nobot_aux','nSVtotal(nJet, JetBotInd, JetGoodInd, Jet_pt_nom, Jet_eta, Jet_phi, nSV, SV_pt, SV_phi, SV_eta, SV_charge)')
	df[s] = df[s].Define('nSV_nobot','nSV_nobot_aux.size()')
	df[s] = df[s].Define('InvM3_aux','InvMassBot(JetGoodInd, JetBotInd, JetSVInd, Jet_pt_nom, Jet_eta, Jet_phi, JetnotSVInd, Jet_mass_nom)')
	df[s] = df[s].Define('InvM_bot_closer','InvM3_aux[0]')
	df[s] = df[s].Define('InvM_bot_farther','InvM3_aux[1]')
	df[s] = df[s].Define('second_muon_aux','muonSecInd(nMuon,Muon_pt,Muon_eta,Muon_pfRelIso04_all, Muon_tightId,'+str(muon_pt[year])+')')
	df[s] = df[s].Define('second_electron_aux','elSecInd(nElectron, Electron_pt, Electron_eta, Electron_pfRelIso03_all, Electron_cutBased, Electron_mvaFall17V2Iso_WP80, Electron_mvaFall17V2Iso_WP90,'+str(el_pt[year])+')')
	df[s] = df[s].Define('second_muon_pt','second_muon_aux.size() > 0 ? Muon_pt[second_muon_aux[0]] : -1')
	df[s] = df[s].Define('second_electron_pt','second_electron_aux.size() > 0 ? Electron_pt[second_electron_aux[0]] : -1')
	df[s] = df[s].Define('my_tau','tau_best(JetMuonInd, Tau_pt, Tau_jetIdx)')
	df[s] = df[s].Define('my_tau_pt','my_tau > -1 ? Tau_pt[my_tau] : -10.')
	df[s] = df[s].Define('my_tau_mass','my_tau > -1 ? Tau_mass[my_tau] : -1.')

############ Gen level definitions

if mode == "mc":
       for s in samples:
                df[s] = df[s].Define('jet_muon_flavourH','Jet_hadronFlavour[JetSVInd[0]]')
                df[s] = df[s].Define('jet_notmuon_flavourH','Jet_hadronFlavour[JetnotSVInd]')
                df[s] = df[s].Define('jet_muon_flavourP','Jet_partonFlavour[JetSVInd[0]]')
                df[s] = df[s].Define('jet_notmuon_flavourP','Jet_partonFlavour[JetnotSVInd]')
                df[s] = df[s].Define('jet_bot1_flavourP','Jet_partonFlavour[JetGoodInd[JetBotInd[0]]]')
                df[s] = df[s].Define('jet_bot2_flavourP','Jet_partonFlavour[JetGoodInd[JetBotInd[1]]]')
                df[s] = df[s].Define('btag_sf_aux1','(fabs(jet_bot1_flavourP) == 4 || fabs(jet_bot1_flavourP) == 5) ? btag_MED_sf[JetGoodInd[JetBotInd[0]]] : btag_MED_incl_sf[JetGoodInd[JetBotInd[0]]]')
                df[s] = df[s].Define('btag_sf_aux2','(fabs(jet_bot2_flavourP) == 4 || fabs(jet_bot2_flavourP) == 5) ? btag_LOO_sf[JetGoodInd[JetBotInd[1]]] : btag_LOO_incl_sf[JetGoodInd[JetBotInd[1]]]')
                df[s] = df[s].Define('btag_sf','btag_sf_aux1*btag_sf_aux2')

#########################################
######## Further corrections ############
#########################################

############ BR corrections ##############

from branchingfractions_corrections import *

if mode == "mc":
        for s in samples:
                if s[-1]=="B":
                       kwargs = {"year":s[-5:],"isMC":True, "isUL":True}
                else:
                       kwargs = {"year":s[-4:],"isMC":True, "isUL":True}
                       #print(kwargs)
                b= branchingfractions_SV_corRDF(**kwargs)
                df[s] = b().run(df[s])

#################################
#######    GEN FILTER    ########
#################################

########## Wjets sectioning for W plus c discrimination

if mode == "mc":
        for s in samples:
                if (s[0]+s[1]+s[2]+s[3]+s[4] == "wjets" and (s[-1]=="6" or s[-1]=="7" or s[-1]=="8" or s[-1]=="B")):
                        df[s+"_charm"] = df[s].Filter('isWplusc == 2')
                        df[s+"_bottom"] = df[s].Filter('isWplusc == 3')
                        df[s+"_doublecharm"] = df[s].Filter('isWplusc == 1')
                        df[s+"_light"] = df[s].Filter('isWplusc == 0')
                        ## Samples correction
                        samples.append(s+"_charm")
                        samples.append(s+"_bottom")
                        samples.append(s+"_doublecharm")
                        samples.append(s+"_light")

samples = [s for s in samples if not (s[0]+s[1]+s[2]+s[3]+s[4] == "wjets" and (s[-1]=="6" or s[-1]=="7" or s[-1]=="8" or s[-1]=="B"))]

########## ttbar sectioning for charm discrimination

if mode == "mc":
        for s in samples:
                if (s[0]+s[1]+s[2]+s[3]+s[4] == "ttbar" and (s[-1]=="6" or s[-1]=="7" or s[-1]=="8" or s[-1]=="B")):
                        df[s+"_charm"] = df[s].Filter('isttbarC')
                        df[s+"_nocharm"] = df[s].Filter('!isttbarC')
                        ## Samples correction
                        samples.append(s+"_charm")
                        samples.append(s+"_nocharm")

samples = [s for s in samples if not (s[0]+s[1]+s[2]+s[3]+s[4] == "ttbar"  and (s[-1]=="6" or s[-1]=="7" or s[-1]=="8" or s[-1]=="B"))]

########## ST sectioning for charm discrimination

if mode == "mc":
        for s in samples:
                if (s[0]+s[1] == "st" and (s[-1]=="6" or s[-1]=="7" or s[-1]=="8" or s[-1]=="B")):
                        df[s+"_charm"] = df[s].Filter('isttbarC')
                        df[s+"_nocharm"] = df[s].Filter('!isttbarC')
                        ## Samples correction
                        samples.append(s+"_charm")
                        samples.append(s+"_nocharm")

samples = [s for s in samples if not (s[0]+s[1] == "st"  and (s[-1]=="6" or s[-1]=="7" or s[-1]=="8" or s[-1]=="B"))]

print(samples)

#################### Final definitions and filters ###############################

for s in samples:
        if args.ssos: df[s] = df[s].Define('weight_aux','weightSSOS')
        else: df[s] = df[s].Define('weight_aux','1')
        df[s] = df[s].Define('deltaR_jetM_lep','nMuonGood>0 ? ROOT::VecOps::DeltaR(Muon_eta[MuonGoodInd[0]],Jet_eta[JetSVInd[0]] , Muon_phi[MuonGoodInd[0]], Jet_phi[JetSVInd[0]]) : ROOT::VecOps::DeltaR(Electron_eta[ElectronGoodInd[0]],Jet_eta[JetSVInd[0]] , Electron_phi[ElectronGoodInd[0]], Jet_phi[JetSVInd[0]])')
        df[s] = df[s].Define('InvM_jetM_lep', 'nMuonGood>0 ? InvariantM(Jet_pt_nom[JetSVInd[0]],Jet_eta[JetSVInd[0]],Jet_phi[JetSVInd[0]],Jet_mass_nom[JetSVInd[0]],Muon_pt[MuonGoodInd[0]],Muon_eta[MuonGoodInd[0]],Muon_phi[MuonGoodInd[0]],Muon_mass[MuonGoodInd[0]]) : InvariantM(Jet_pt_nom[JetSVInd[0]],Jet_eta[JetSVInd[0]],Jet_phi[JetSVInd[0]],Jet_mass_nom[JetSVInd[0]],Electron_pt[ElectronGoodInd[0]],Electron_eta[ElectronGoodInd[0]],Electron_phi[ElectronGoodInd[0]], Electron_mass[ElectronGoodInd[0]])')
        #### Filters
        df[s] = df[s].Filter('transverse_mass > 50')
        df[s] = df[s].Filter('sv_jet_ntracks > 2')
        df[s] = df[s].Filter('jet_bot1_btag >'+str(cuts_btag[year][1]))
        df[s] = df[s].Filter('jet_bot2_btag >'+str(cuts_btag[year][0]))
        df_muon[s] = df[s].Filter('nMuonGood>0')
        df_electron[s] = df[s].Filter('nElectronGood >0')
        if args.type == "mc":
               df_muon[s] = df_muon[s].Define('weightSSOS_final','weight_aux*btag_sf*musf_tight_id[MuonGoodInd[0]]*musf_tight_reliso[MuonGoodInd[0]]*puWeight*PUjetID_SF')
               df_electron[s] = df_electron[s].Define('weightSSOS_final','weight_aux*elesf_wp80iso[ElectronGoodInd[0]]*btag_sf*puWeight*PUjetID_SF')
        else:
               df_muon[s] = df_muon[s].Define('weightSSOS_final','weight_aux')
               df_electron[s] = df_electron[s].Define('weightSSOS_final','weight_aux')

############################################################
####################     HISTS    ##########################
############################################################

hist_nJetGood_M = {}
hist_nJetGood_E = {}
hist_nSV_M = {}
hist_nSV_E = {}
hist_nSV_nobot_M = {}
hist_nSV_nobot_E = {}
hist_nLooseLepton_M = {}
hist_nLooseLepton_E = {}
hist_nTau_M = {}
hist_nTau_E = {}
hist_second_muon_pt_M = {}
hist_second_muon_pt_E = {}
hist_second_electron_pt_M = {}
hist_second_electron_pt_E = {}
hist_jet_muon_pt_M = {}
hist_jet_muon_nmu_M = {}
hist_jet_muon_mass_M = {}
hist_jet_notmuon_pt_M = {}
hist_jet_muon_eta_M = {}
hist_jet_notmuon_eta_M = {}
hist_jet_notmuon_nmu_M = {}
hist_jet_notmuon_qgl_M = {}
hist_jet_notmuon_mass_M = {}
hist_jet_muon_pt_E = {}
hist_jet_muon_nmu_E = {}
hist_jet_muon_mass_E = {}
hist_jet_notmuon_pt_E = {}
hist_jet_muon_eta_E = {}
hist_jet_notmuon_eta_E = {}
hist_jet_notmuon_nmu_E = {}
hist_jet_notmuon_qgl_E = {}
hist_jet_notmuon_mass_E = {}
hist_jet_bot1_nsv_good_M = {}
hist_jet_bot1_nsv_good_E = {}
hist_jet_bot2_nsv_good_M = {}
hist_jet_bot2_nsv_good_E = {}
hist_jet_muon_nsv_good_M = {}
hist_jet_muon_nsv_good_E = {}
hist_lepton_pt_M = {}
hist_lepton_eta_M = {}
hist_lepton_pt_E = {}
hist_lepton_eta_E = {}
hist_sv_jet_pt_M = {}
hist_sv_jet_eta_M = {}
hist_sv_jet_mass_M = {}
hist_sv_jet_pangle_M = {}
hist_sv_jet_chi_M = {}
hist_sv_jet_ndof_M = {}
hist_sv_jet_ntracks_M = {}
hist_sv_jet_charge_M = {}
hist_sv_jet_pt_E = {}
hist_sv_jet_eta_E = {}
hist_sv_jet_mass_E = {}
hist_sv_jet_pangle_E = {}
hist_sv_jet_chi_E = {}
hist_sv_jet_ndof_E = {}
hist_sv_jet_ntracks_E = {}
hist_sv_jet_charge_E = {}
hist_sv_jet_relpt_M = {}
hist_sv_jet_relpt_E = {}
hist_sv_max_pt_M = {}
hist_sv_max_eta_M = {}
hist_sv_max_mass_M = {}
hist_sv_max_pangle_M = {}
hist_sv_max_chi_M = {}
hist_sv_max_ndof_M = {}
hist_sv_max_ntracks_M = {}
hist_sv_max_charge_M = {}
hist_sv_max_pt_E = {}
hist_sv_max_eta_E = {}
hist_sv_max_mass_E = {}
hist_sv_max_pangle_E = {}
hist_sv_max_chi_E = {}
hist_sv_max_ndof_E = {}
hist_sv_max_ntracks_E = {}
hist_sv_max_charge_E = {}
hist_my_tau_pt_M = {}
hist_my_tau_pt_E = {}
hist_my_tau_mass_M = {}
hist_my_tau_mass_E = {}
hist_InvM_2jets_M = {}
hist_InvM_2jets_E = {}
hist_InvM_jetM_lepM = {}
hist_InvM_jetM_lepE = {}
hist_InvM_bot_closer_M = {}
hist_InvM_bot_closer_E = {}
hist_InvM_bot_farther_M = {}
hist_InvM_bot_farther_E = {}
hist_deltaR_jetM_jetNM_M = {}
hist_deltaR_jetM_jetNM_E = {}
hist_deltaphi_jetM_jetNM_M = {}
hist_deltaphi_jetM_jetNM_E = {}
hist_deltaeta_jetM_jetNM_M = {}
hist_deltaeta_jetM_jetNM_E = {}
hist_deltapt_jetM_jetNM_M = {}
hist_deltapt_jetM_jetNM_E = {}
hist_MET_M = {}
hist_MET_E = {}
hist_tranverse_massM = {}
hist_tranverse_massE = {}
hist_tracks_jetM_M = {}
hist_tracks_jetNM_M = {}
hist_tracks_jetM_E = {}
hist_tracks_jetNM_E = {}
hist_EMN_jetM_M = {}
hist_EMC_jetM_M = {}
hist_EMN_jetM_E = {}
hist_EMC_jetM_E = {}
hist_EMtotal_jetM_M = {}
hist_EMtotal_jetM_E = {}
hist_SSOS_M = {}
hist_SSOS_E = {}
hist_sv_jet_sigxy_M = {}
hist_sv_jet_sigxy_E = {}
hist_sv_jet_sigr_M = {}
hist_sv_jet_sigr_E = {}
hist_sv_jet_xy_M = {}
hist_sv_jet_xy_E = {}
hist_sv_jet_r_M = {}
hist_sv_jet_r_E = {}
hist_sv_max_sigxy_M = {}
hist_sv_max_sigxy_E = {}
hist_sv_max_sigr_M = {}
hist_sv_max_sigr_E = {}
hist_sv_max_xy_M = {}
hist_sv_max_xy_E = {}
hist_sv_max_r_M = {}
hist_sv_max_r_E = {}
hist_pT_sum_M = {}
hist_pT_sum_E = {}
hist_pT_product_M = {}
hist_pT_product_E = {}
hist_deltaR_lep_2jets_M = {}
hist_deltaR_lep_2jets_E = {}
hist_deltaphi_MET_2jets_M = {}
hist_deltaphi_MET_2jets_E = {}
hist_deltaphi_lephad_M = {}
hist_deltaphi_lephad_E = {}
hist_eta_2jets_M = {}
hist_eta_2jets_E = {}
hist_pt_2jets_M = {}
hist_pt_2jets_E = {}
hist_pt_Wlep_M = {}
hist_pt_Wlep_E = {}
hist_deltaR_lephad_M = {}
hist_deltaR_lephad_E = {}
hist_deltaphi_lep_2jets_M = {}
hist_deltaphi_lep_2jets_E = {}
hist_deltaeta_lephad_M = {}
hist_deltaeta_lephad_E = {}
hist_deltaeta_lep_2jets_M = {}
hist_deltaeta_lep_2jets_E = {}
hist_deltapt_lephad_M = {}
hist_deltapt_lephad_E = {}
hist_deltapt_lep_2jets_M = {}
hist_deltapt_lep_2jets_E = {}
hist_jet_bot1_btag_M = {}
hist_jet_bot1_btag_E = {}
hist_jet_bot2_btag_M = {}
hist_jet_bot2_btag_E = {}
hist_jet_muon_btag_M = {}
hist_jet_muon_btag_E = {}
hist_jet_notmuon_btag_M = {}
hist_jet_notmuon_btag_E = {}
hist_jet_muon_flavourH_M = {}
hist_jet_muon_flavourH_E = {}
hist_jet_notmuon_flavourH_M = {}
hist_jet_notmuon_flavourH_E = {}
hist_jet_muon_flavourP_M = {}
hist_jet_muon_flavourP_E = {}
hist_jet_notmuon_flavourP_M = {}
hist_jet_notmuon_flavourP_E = {}
hist_deltaR_lep_jet1_M = {}
hist_deltaR_lep_jet2_M = {}
hist_deltaPhi_lep_jet1_M = {}
hist_deltaPhi_lep_jet2_M = {}
hist_deltaEta_lep_jet1_M = {}
hist_deltaEta_lep_jet2_M = {}
hist_deltaR_lep_jet1_E = {}
hist_deltaR_lep_jet2_E = {}
hist_deltaPhi_lep_jet1_E = {}
hist_deltaPhi_lep_jet2_E = {}
hist_deltaEta_lep_jet1_E = {}
hist_deltaEta_lep_jet2_E = {}
hist_pT_proy_M = {}
hist_pT_proy_E = {}
hist_pT_sum_2J_M = {}
hist_pT_sum_2J_E = {}
#### hists sf
hist_br_sf_M = {}
hist_br_sf_E = {}
hist_frag_sf_M = {}
hist_frag_sf_E = {}

for s in samples:
        hist_nJetGood_M[s] = df_muon[s].Histo1D(("nJetGood_M","",10,0,10),"nJetGood","weightSSOS_final")
        hist_nJetGood_E[s] = df_electron[s].Histo1D(("nJetGood_E","",10,0,10),"nJetGood","weightSSOS_final")

        hist_nSV_M[s] = df_muon[s].Histo1D(("nSV_M","",10,0,10),"nSV","weightSSOS_final")
        hist_nSV_E[s] = df_electron[s].Histo1D(("nSV_E","",10,0,10),"nSV","weightSSOS_final")

        hist_nSV_nobot_M[s] = df_muon[s].Histo1D(("nSV_nobot_M","",10,0,10),"nSV_nobot","weightSSOS_final")
        hist_nSV_nobot_E[s] = df_electron[s].Histo1D(("nSV_nobot_E","",10,0,10),"nSV_nobot","weightSSOS_final")

        hist_nLooseLepton_M[s] = df_muon[s].Histo1D(("nLooseLepton_M","",10,0,10),"nLooseLepton","weightSSOS_final")
        hist_nLooseLepton_E[s] = df_electron[s].Histo1D(("nLooseLepton_E","",10,0,10),"nLooseLepton","weightSSOS_final")

        hist_nTau_M[s] = df_muon[s].Histo1D(("nTau_M","",10,0,10),"nTau","weightSSOS_final")
        hist_nTau_E[s] = df_electron[s].Histo1D(("nTau_E","",10,0,10),"nTau","weightSSOS_final")

        hist_second_muon_pt_M[s] = df_muon[s].Histo1D(("second_muon_pt_M","",100,-1,30),"second_muon_pt","weightSSOS_final")
        hist_second_muon_pt_E[s] = df_electron[s].Histo1D(("second_muon_pt_E","",100,-1,30),"second_muon_pt","weightSSOS_final")
        hist_second_electron_pt_M[s] = df_muon[s].Histo1D(("second_electron_pt_M","",100,-1,30),"second_electron_pt","weightSSOS_final")
        hist_second_electron_pt_E[s] = df_electron[s].Histo1D(("second_electron_pt_E","",100,-1,30),"second_electron_pt","weightSSOS_final")

        hist_jet_muon_pt_M[s] = df_muon[s].Histo1D(("jet_sv_pt_M","",50,20,120),"jet_muon_pt","weightSSOS_final")
        hist_jet_muon_nmu_M[s] = df_muon[s].Histo1D(("jet_sv_nmu_M","",10,0,10),"jet_muon_nmu","weightSSOS_final")
        hist_jet_muon_mass_M[s] = df_muon[s].Histo1D(("jet_sv_mass_M","",40,0,40),"jet_muon_mass","weightSSOS_final")
        hist_jet_notmuon_pt_M[s] = df_muon[s].Histo1D(("jet_not_sv_pt_M","",50,20,120),"jet_notmuon_pt","weightSSOS_final")
        hist_jet_muon_eta_M[s] = df_muon[s].Histo1D(("jet_sv_eta_M","",80,-4,4),"jet_muon_eta","weightSSOS_final")
        hist_jet_notmuon_eta_M[s] = df_muon[s].Histo1D(("jet_not_sv_eta_M","",80,-4,4),"jet_notmuon_eta","weightSSOS_final")
        hist_jet_notmuon_nmu_M[s] = df_muon[s].Histo1D(("jet_notsv_nmu_M","",10,0,10),"jet_notmuon_nmu","weightSSOS_final")
        hist_jet_notmuon_mass_M[s] = df_muon[s].Histo1D(("jet_notsv_mass_M","",40,0,40),"jet_notmuon_mass","weightSSOS_final")
        hist_jet_notmuon_qgl_M[s] = df_muon[s].Histo1D(("jet_notsv_qgl_M","",100,0,1),"jet_notmuon_qgl","weightSSOS_final")

        hist_jet_muon_pt_E[s] = df_electron[s].Histo1D(("jet_sv_pt_E","",50,20,120),"jet_muon_pt","weightSSOS_final")
        hist_jet_muon_nmu_E[s] = df_electron[s].Histo1D(("jet_sv_nmu_E","",10,0,10),"jet_muon_nmu","weightSSOS_final")
        hist_jet_muon_mass_E[s] = df_electron[s].Histo1D(("jet_sv_mass_E","",40,0,40),"jet_muon_mass","weightSSOS_final")
        hist_jet_notmuon_pt_E[s] = df_electron[s].Histo1D(("jet_not_sv_pt_E","",50,20,120),"jet_notmuon_pt","weightSSOS_final")
        hist_jet_muon_eta_E[s] = df_electron[s].Histo1D(("jet_sv_eta_E","",80,-4,4),"jet_muon_eta","weightSSOS_final")
        hist_jet_notmuon_eta_E[s] = df_electron[s].Histo1D(("jet_not_sv_eta_E","",80,-4,4),"jet_notmuon_eta","weightSSOS_final")
        hist_jet_notmuon_nmu_E[s] = df_electron[s].Histo1D(("jet_notsv_nmu_E","",10,0,10),"jet_notmuon_nmu","weightSSOS_final")
        hist_jet_notmuon_mass_E[s] = df_electron[s].Histo1D(("jet_notsv_mass_E","",40,0,40),"jet_notmuon_mass","weightSSOS_final")
        hist_jet_notmuon_qgl_E[s] = df_electron[s].Histo1D(("jet_notsv_qgl_E","",100,0,1),"jet_notmuon_qgl","weightSSOS_final")

        hist_jet_bot1_nsv_good_M[s] = df_muon[s].Histo1D(("jet_bot1_nsv_good_M","",10,0,10),"jet_bot1_nsv_good","weightSSOS_final")
        hist_jet_bot1_nsv_good_E[s] = df_electron[s].Histo1D(("jet_bot1_nsv_good_E","",10,0,10),"jet_bot1_nsv_good","weightSSOS_final")
        hist_jet_bot2_nsv_good_M[s] = df_muon[s].Histo1D(("jet_bot2_nsv_good_M","",10,0,10),"jet_bot2_nsv_good","weightSSOS_final")
        hist_jet_bot2_nsv_good_E[s] = df_electron[s].Histo1D(("jet_bot2_nsv_good_E","",10,0,10),"jet_bot2_nsv_good","weightSSOS_final")
        hist_jet_muon_nsv_good_M[s] = df_muon[s].Histo1D(("jet_sv_nsv_good_M","",10,0,10),"jet_muon_nsv_good","weightSSOS_final")
        hist_jet_muon_nsv_good_E[s] = df_electron[s].Histo1D(("jet_sv_nsv_good_E","",10,0,10),"jet_muon_nsv_good","weightSSOS_final")

        hist_lepton_pt_M[s] = df_muon[s].Histo1D(("lepton_pt_M","",50,20,120),"lepton_pt","weightSSOS_final")
        hist_lepton_eta_M[s] = df_muon[s].Histo1D(("lepton_eta_M","",80,-4,4),"lepton_eta","weightSSOS_final")

        hist_lepton_pt_E[s] = df_electron[s].Histo1D(("lepton_pt_E","",50,20,120),"lepton_pt","weightSSOS_final")
        hist_lepton_eta_E[s] = df_electron[s].Histo1D(("lepton_eta_E","",80,-4,4),"lepton_eta","weightSSOS_final")

        hist_sv_jet_pt_M[s] = df_muon[s].Histo1D(("sv_jet_pt_M","",50,0,25),"sv_jet_pt","weightSSOS_final")
        hist_sv_jet_eta_M[s] = df_muon[s].Histo1D(("sv_jet_eta_M","",80,-4,4),"sv_jet_eta","weightSSOS_final")
        hist_sv_jet_mass_M[s] = df_muon[s].Histo1D(("sv_jet_mass_M","",50,0,5),"sv_jet_mass","weightSSOS_final")
        hist_sv_jet_pangle_M[s] = df_muon[s].Histo1D(("sv_jet_pangle_M","",50,1,3.5),"sv_jet_pangle","weightSSOS_final")
        hist_sv_jet_chi_M[s] = df_muon[s].Histo1D(("sv_jet_chi_M","",50,0,20),"sv_jet_chi","weightSSOS_final")
        hist_sv_jet_ndof_M[s] = df_muon[s].Histo1D(("sv_jet_ndof_M","",10,0,10),"sv_jet_ndof","weightSSOS_final")
        hist_sv_jet_charge_M[s] = df_muon[s].Histo1D(("sv_jet_charge_M","",10,-5,5),"sv_jet_charge","weightSSOS_final")

        hist_sv_jet_pt_E[s] = df_electron[s].Histo1D(("sv_jet_pt_E","",50,0,25),"sv_jet_pt","weightSSOS_final")
        hist_sv_jet_eta_E[s] = df_electron[s].Histo1D(("sv_jet_eta_E","",80,-4,4),"sv_jet_eta","weightSSOS_final")
        hist_sv_jet_mass_E[s] = df_electron[s].Histo1D(("sv_jet_mass_E","",50,0,5),"sv_jet_mass","weightSSOS_final")
        hist_sv_jet_pangle_E[s] = df_electron[s].Histo1D(("sv_jet_pangle_E","",50,1,3.5),"sv_jet_pangle","weightSSOS_final")
        hist_sv_jet_chi_E[s] = df_electron[s].Histo1D(("sv_jet_chi_E","",50,0,20),"sv_jet_chi","weightSSOS_final")
        hist_sv_jet_ndof_E[s] = df_electron[s].Histo1D(("sv_jet_ndof_E","",10,0,10),"sv_jet_ndof","weightSSOS_final")
        hist_sv_jet_charge_E[s] = df_electron[s].Histo1D(("sv_jet_charge_E","",10,-5,5),"sv_jet_charge","weightSSOS_final")

        hist_sv_jet_ntracks_M[s] = df_muon[s].Histo1D(("sv_jet_ntracks_M","",15,0,15),"sv_jet_ntracks","weightSSOS_final")
        hist_sv_jet_ntracks_E[s] = df_electron[s].Histo1D(("sv_jet_ntracks_E","",15,0,15),"sv_jet_ntracks","weightSSOS_final")

        hist_sv_jet_relpt_M[s] = df_muon[s].Histo1D(("sv_jet_relpt_M","",50,0,1),"sv_jet_relpt","weightSSOS_final")
        hist_sv_jet_relpt_E[s] = df_electron[s].Histo1D(("sv_jet_relpt_E","",50,0,1),"sv_jet_relpt","weightSSOS_final")

        hist_sv_max_pt_M[s] = df_muon[s].Histo1D(("sv_max_pt_M","",50,0,25),"sv_max_pt","weightSSOS_final")
        hist_sv_max_eta_M[s] = df_muon[s].Histo1D(("sv_max_eta_M","",80,-4,4),"sv_max_eta","weightSSOS_final")
        hist_sv_max_mass_M[s] = df_muon[s].Histo1D(("sv_max_mass_M","",50,0,5),"sv_max_mass","weightSSOS_final")
        hist_sv_max_pangle_M[s] = df_muon[s].Histo1D(("sv_max_pangle_M","",50,1,3.5),"sv_max_pangle","weightSSOS_final")
        hist_sv_max_chi_M[s] = df_muon[s].Histo1D(("sv_max_chi_M","",50,0,20),"sv_max_chi","weightSSOS_final")
        hist_sv_max_ndof_M[s] = df_muon[s].Histo1D(("sv_max_ndof_M","",10,0,10),"sv_max_ndof","weightSSOS_final")
        hist_sv_max_charge_M[s] = df_muon[s].Histo1D(("sv_max_charge_M","",10,-5,5),"sv_max_charge","weightSSOS_final")

        hist_sv_max_pt_E[s] = df_electron[s].Histo1D(("sv_max_pt_E","",50,0,25),"sv_max_pt","weightSSOS_final")
        hist_sv_max_eta_E[s] = df_electron[s].Histo1D(("sv_max_eta_E","",80,-4,4),"sv_max_eta","weightSSOS_final")
        hist_sv_max_mass_E[s] = df_electron[s].Histo1D(("sv_max_mass_E","",50,0,5),"sv_max_mass","weightSSOS_final")
        hist_sv_max_pangle_E[s] = df_electron[s].Histo1D(("sv_max_pangle_E","",50,1,3.5),"sv_max_pangle","weightSSOS_final")
        hist_sv_max_chi_E[s] = df_electron[s].Histo1D(("sv_max_chi_E","",50,0,20),"sv_max_chi","weightSSOS_final")
        hist_sv_max_ndof_E[s] = df_electron[s].Histo1D(("sv_max_ndof_E","",10,0,10),"sv_max_ndof","weightSSOS_final")
        hist_sv_max_charge_E[s] = df_electron[s].Histo1D(("sv_max_charge_E","",10,-5,5),"sv_max_charge","weightSSOS_final")

        hist_sv_max_ntracks_M[s] = df_muon[s].Histo1D(("sv_max_ntracks_M","",15,0,15),"sv_max_ntracks","weightSSOS_final")
        hist_sv_max_ntracks_E[s] = df_electron[s].Histo1D(("sv_max_ntracks_E","",15,0,15),"sv_max_ntracks","weightSSOS_final")

        hist_InvM_2jets_M[s] = df_muon[s].Histo1D(("InvM_2jets_M","",100,0,300),"InvM_2jets","weightSSOS_final")
        hist_InvM_2jets_E[s] = df_electron[s].Histo1D(("InvM_2jets_E","",100,0,300),"InvM_2jets","weightSSOS_final")
        hist_InvM_jetM_lepM[s] = df_muon[s].Histo1D(("InvM_jetM_lepM","",100,0,300),"InvM_jetM_lep","weightSSOS_final")
        hist_InvM_jetM_lepE[s] = df_electron[s].Histo1D(("InvM_jetM_lepE","",100,0,300),"InvM_jetM_lep","weightSSOS_final")

        hist_InvM_bot_closer_M[s] = df_muon[s].Histo1D(("InvM_bot_closer_M","",100,0,300),"InvM_bot_closer","weightSSOS_final")
        hist_InvM_bot_closer_E[s] = df_electron[s].Histo1D(("InvM_bot_closer_E","",100,0,300),"InvM_bot_closer","weightSSOS_final")
        hist_InvM_bot_farther_M[s] = df_muon[s].Histo1D(("InvM_bot_farther_M","",100,0,300),"InvM_bot_farther","weightSSOS_final")
        hist_InvM_bot_farther_E[s] = df_electron[s].Histo1D(("InvM_bot_farther_E","",100,0,300),"InvM_bot_farther","weightSSOS_final")

        hist_MET_M[s] = df_muon[s].Histo1D(("MET_pt_M","",100,0,150),"MET_pt","weightSSOS_final")
        hist_MET_E[s] = df_electron[s].Histo1D(("MET_pt_E","",100,0,150),"MET_pt","weightSSOS_final")

        hist_deltaR_jetM_jetNM_M[s] = df_muon[s].Histo1D(("deltaR_jetM_jetNM_M","",100,0,5),"deltaR_jetM_jetNM","weightSSOS_final")
        hist_deltaR_jetM_jetNM_E[s] = df_electron[s].Histo1D(("deltaR_jetM_jetNM_E","",100,0,5),"deltaR_jetM_jetNM","weightSSOS_final")

        hist_deltaphi_jetM_jetNM_M[s] = df_muon[s].Histo1D(("deltaphi_jetM_jetNM_M","",100,0,5),"deltaphi_jetM_jetNM","weightSSOS_final")
        hist_deltaphi_jetM_jetNM_E[s] = df_electron[s].Histo1D(("deltaphi_jetM_jetNM_E","",100,0,5),"deltaphi_jetM_jetNM","weightSSOS_final")

        hist_deltaeta_jetM_jetNM_M[s] = df_muon[s].Histo1D(("deltaeta_jetM_jetNM_M","",100,0,5),"deltaeta_jetM_jetNM","weightSSOS_final")
        hist_deltaeta_jetM_jetNM_E[s] = df_electron[s].Histo1D(("deltaeta_jetM_jetNM_E","",100,0,5),"deltaeta_jetM_jetNM","weightSSOS_final")

        hist_deltapt_jetM_jetNM_M[s] = df_muon[s].Histo1D(("deltapt_jetM_jetNM_M","",100,0,50),"deltapt_jetM_jetNM","weightSSOS_final")
        hist_deltapt_jetM_jetNM_E[s] = df_electron[s].Histo1D(("deltapt_jetM_jetNM_E","",100,0,50),"deltapt_jetM_jetNM","weightSSOS_final")

        hist_tranverse_massM[s] = df_muon[s].Histo1D(("transverse_massM","",50,0,150),"transverse_mass","weightSSOS_final")
        hist_tranverse_massE[s] = df_electron[s].Histo1D(("transverse_massE","",50,0,150),"transverse_mass","weightSSOS_final")

        hist_tracks_jetM_M[s] = df_muon[s].Histo1D(("tracks_jetM_M","",60,0,60),"tracks_jetM","weightSSOS_final")
        hist_tracks_jetNM_M[s] = df_muon[s].Histo1D(("tracks_jetNM_M","",60,0,60),"tracks_jetNM","weightSSOS_final")

        hist_tracks_jetM_E[s] = df_electron[s].Histo1D(("tracks_jetM_E","",60,0,60),"tracks_jetM","weightSSOS_final")
        hist_tracks_jetNM_E[s] = df_electron[s].Histo1D(("tracks_jetNM_E","",60,0,60),"tracks_jetNM","weightSSOS_final")

        hist_EMN_jetM_M[s] = df_muon[s].Histo1D(("EMN_jetM_M","",60,0,1),"EMN_jetM","weightSSOS_final")
        hist_EMC_jetM_M[s] = df_muon[s].Histo1D(("EMC_jetM_M","",60,0,1),"EMC_jetM","weightSSOS_final")
        hist_EMN_jetM_E[s] = df_electron[s].Histo1D(("EMN_jetM_E","",60,0,1),"EMN_jetM","weightSSOS_final")
        hist_EMC_jetM_E[s] = df_electron[s].Histo1D(("EMC_jetM_E","",60,0,1),"EMC_jetM","weightSSOS_final")

        hist_EMtotal_jetM_M[s] = df_muon[s].Histo1D(("EMtotal_jetM_M","",60,0,1),"EMtotal_jetM","weightSSOS_final")
        hist_EMtotal_jetM_E[s] = df_electron[s].Histo1D(("EMtotal_jetM_E","",60,0,1),"EMtotal_jetM","weightSSOS_final")

        hist_sv_jet_sigxy_M[s] = df_muon[s].Histo1D(("sv_jet_sigxy_M","",100,-4,4),"sv_jet_sigxy","weightSSOS_final")
        hist_sv_jet_sigxy_E[s] = df_electron[s].Histo1D(("sv_jet_sigxy_E","",100,-4,4),"sv_jet_sigxy","weightSSOS_final")

        hist_sv_jet_sigr_M[s] = df_muon[s].Histo1D(("sv_jet_sigr_M","",100,0,30),"sv_jet_sigr","weightSSOS_final")
        hist_sv_jet_sigr_E[s] = df_electron[s].Histo1D(("sv_jet_sigr_E","",100,0,30),"sv_jet_sigr","weightSSOS_final")

        hist_sv_jet_xy_M[s] = df_muon[s].Histo1D(("sv_jet_xy_M","",100,-4,4),"sv_jet_xy","weightSSOS_final")
        hist_sv_jet_xy_E[s] = df_electron[s].Histo1D(("sv_jet_xy_E","",100,-4,4),"sv_jet_xy","weightSSOS_final")

        hist_sv_jet_r_M[s] = df_muon[s].Histo1D(("sv_jet_r_M","",100,0,6),"sv_jet_r","weightSSOS_final")
        hist_sv_jet_r_E[s] = df_electron[s].Histo1D(("sv_jet_r_E","",100,0,6),"sv_jet_r","weightSSOS_final")

        hist_sv_max_sigxy_M[s] = df_muon[s].Histo1D(("sv_max_sigxy_M","",100,-4,4),"sv_max_sigxy","weightSSOS_final")
        hist_sv_max_sigxy_E[s] = df_electron[s].Histo1D(("sv_max_sigxy_E","",100,-4,4),"sv_max_sigxy","weightSSOS_final")

        hist_sv_max_sigr_M[s] = df_muon[s].Histo1D(("sv_max_sigr_M","",100,3,30),"sv_max_sigr","weightSSOS_final")
        hist_sv_max_sigr_E[s] = df_electron[s].Histo1D(("sv_max_sigr_E","",100,3,30),"sv_max_sigr","weightSSOS_final")

        hist_sv_max_xy_M[s] = df_muon[s].Histo1D(("sv_max_xy_M","",100,-4,4),"sv_max_xy","weightSSOS_final")
        hist_sv_max_xy_E[s] = df_electron[s].Histo1D(("sv_max_xy_E","",100,-4,4),"sv_max_xy","weightSSOS_final")

        hist_sv_max_r_M[s] = df_muon[s].Histo1D(("sv_max_r_M","",100,0,6),"sv_max_r","weightSSOS_final")
        hist_sv_max_r_E[s] = df_electron[s].Histo1D(("sv_max_r_E","",100,0,6),"sv_max_r","weightSSOS_final")

        hist_my_tau_pt_M[s] = df_muon[s].Histo1D(("my_tau_pt_M","",100,-10,90),"my_tau_pt","weightSSOS_final")
        hist_my_tau_pt_E[s] = df_electron[s].Histo1D(("my_tau_pt_E","",100,-10,90),"my_tau_pt","weightSSOS_final")
        hist_my_tau_mass_M[s] = df_muon[s].Histo1D(("my_tau_mass_M","",100,-1,20),"my_tau_mass","weightSSOS_final")
        hist_my_tau_mass_E[s] = df_electron[s].Histo1D(("my_tau_mass_E","",100,-1,20),"my_tau_mass","weightSSOS_final")

        hist_SSOS_M[s] = df_muon[s].Histo1D(("SSOS_M","",6,-3,3),"SVLepSign","weightSSOS_final")
        hist_SSOS_E[s] = df_electron[s].Histo1D(("SSOS_E","",6,-3,3),"SVLepSign","weightSSOS_final")

        hist_pT_sum_M[s] = df_muon[s].Histo1D(("pT_sum_M","",100,0,1),"pT_sum","weightSSOS_final")
        hist_pT_sum_E[s] = df_electron[s].Histo1D(("pT_sum_E","",100,0,1),"pT_sum","weightSSOS_final")

        hist_pT_product_M[s] = df_muon[s].Histo1D(("pT_product_M","",100,-1,1),"pT_product","weightSSOS_final")
        hist_pT_product_E[s] = df_electron[s].Histo1D(("pT_product_E","",100,-1,1),"pT_product","weightSSOS_final")

        hist_deltaR_lep_2jets_M[s] = df_muon[s].Histo1D(("deltaR_lep_2jets_M","",100,0,5),"deltaR_lep_2jets","weightSSOS_final")
        hist_deltaR_lep_2jets_E[s] = df_electron[s].Histo1D(("deltaR_lep_2jets_E","",100,0,5),"deltaR_lep_2jets","weightSSOS_final")
        hist_deltaphi_MET_2jets_M[s] = df_muon[s].Histo1D(("deltaphi_MET_2jets_M","",100,0,5),"deltaphi_MET_2jets","weightSSOS_final")
        hist_deltaphi_MET_2jets_E[s] = df_electron[s].Histo1D(("deltaphi_MET_2jets_E","",100,0,5),"deltaphi_MET_2jets","weightSSOS_final")
        hist_deltaphi_lephad_M[s] = df_muon[s].Histo1D(("deltaphi_lephad_M","",100,0,5),"deltaphi_lephad","weightSSOS_final")
        hist_deltaphi_lephad_E[s] = df_electron[s].Histo1D(("deltaphi_lephad_E","",100,0,5),"deltaphi_lephad","weightSSOS_final")

        hist_eta_2jets_M[s] = df_muon[s].Histo1D(("eta_2jets_M","",50,-5,5),"eta_2jets","weightSSOS_final")
        hist_eta_2jets_E[s] = df_electron[s].Histo1D(("eta_2jets_E","",50,-5,5),"eta_2jets","weightSSOS_final")
        hist_pt_2jets_M[s] = df_muon[s].Histo1D(("pt_2jets_M","",100,0,200),"pt_2jets","weightSSOS_final")
        hist_pt_2jets_E[s] = df_electron[s].Histo1D(("pt_2jets_E","",100,0,200),"pt_2jets","weightSSOS_final")

        hist_pt_Wlep_M[s] = df_muon[s].Histo1D(("pt_Wlep_M","",100,0,200),"pT_Wlep","weightSSOS_final")
        hist_pt_Wlep_E[s] = df_electron[s].Histo1D(("pt_Wlep_E","",100,0,200),"pT_Wlep","weightSSOS_final")

        hist_deltaR_lephad_M[s] = df_muon[s].Histo1D(("deltaR_lephad_M","",100,0,5),"deltaR_lephad","weightSSOS_final")
        hist_deltaR_lephad_E[s] = df_electron[s].Histo1D(("deltaR_lephad_E","",100,0,5),"deltaR_lephad","weightSSOS_final")
        hist_deltaphi_lep_2jets_M[s] = df_muon[s].Histo1D(("deltaphi_lep_2jets_M","",100,0,5),"deltaphi_lep_2jets","weightSSOS_final")
        hist_deltaphi_lep_2jets_E[s] = df_electron[s].Histo1D(("deltaphi_lep_2jets_E","",100,0,5),"deltaphi_lep_2jets","weightSSOS_final")
        hist_deltaeta_lephad_M[s] = df_muon[s].Histo1D(("deltaeta_lephad_M","",100,0,5),"deltaeta_lephad","weightSSOS_final")
        hist_deltaeta_lephad_E[s] = df_electron[s].Histo1D(("deltaeta_lephad_E","",100,0,5),"deltaeta_lephad","weightSSOS_final")
        hist_deltaeta_lep_2jets_M[s] = df_muon[s].Histo1D(("deltaeta_lep_2jets_M","",100,0,5),"deltaeta_lep_2jets","weightSSOS_final")
        hist_deltaeta_lep_2jets_E[s] = df_electron[s].Histo1D(("deltaeta_lep_2jets_E","",100,0,5),"deltaeta_lep_2jets","weightSSOS_final")
        hist_deltapt_lephad_M[s] = df_muon[s].Histo1D(("deltapt_lephad_M","",100,0,5),"deltapt_lephad","weightSSOS_final")
        hist_deltapt_lephad_E[s] = df_electron[s].Histo1D(("deltapt_lephad_E","",100,0,5),"deltapt_lephad","weightSSOS_final")
        hist_deltapt_lep_2jets_M[s] = df_muon[s].Histo1D(("deltapt_lep_2jets_M","",100,0,5),"deltapt_lep_2jets","weightSSOS_final")
        hist_deltapt_lep_2jets_E[s] = df_electron[s].Histo1D(("deltapt_lep_2jets_E","",100,0,5),"deltapt_lep_2jets","weightSSOS_final")

        hist_deltaR_lep_jet1_M[s] = df_muon[s].Histo1D(("deltaR_lep_jet1_M","",100,0,5),"deltaR_lep_jet1","weightSSOS_final")
        hist_deltaR_lep_jet2_M[s] = df_muon[s].Histo1D(("deltaR_lep_jet2_M","",100,0,5),"deltaR_lep_jet2","weightSSOS_final")
        hist_deltaPhi_lep_jet1_M[s] = df_muon[s].Histo1D(("deltaPhi_lep_jet1_M","",100,0,5),"deltaPhi_lep_jet1","weightSSOS_final")
        hist_deltaPhi_lep_jet2_M[s] = df_muon[s].Histo1D(("deltaPhi_lep_jet2_M","",100,0,5),"deltaPhi_lep_jet2","weightSSOS_final")
        hist_deltaEta_lep_jet1_M[s] = df_muon[s].Histo1D(("deltaEta_lep_jet1_M","",100,0,5),"deltaEta_lep_jet1","weightSSOS_final")
        hist_deltaEta_lep_jet2_M[s] = df_muon[s].Histo1D(("deltaEta_lep_jet2_M","",100,0,5),"deltaEta_lep_jet2","weightSSOS_final")

        hist_deltaR_lep_jet1_E[s] = df_electron[s].Histo1D(("deltaR_lep_jet1_E","",100,0,5),"deltaR_lep_jet1","weightSSOS_final")
        hist_deltaR_lep_jet2_E[s] = df_electron[s].Histo1D(("deltaR_lep_jet2_E","",100,0,5),"deltaR_lep_jet2","weightSSOS_final")
        hist_deltaPhi_lep_jet1_E[s] = df_electron[s].Histo1D(("deltaPhi_lep_jet1_E","",100,0,5),"deltaPhi_lep_jet1","weightSSOS_final")
        hist_deltaPhi_lep_jet2_E[s] = df_electron[s].Histo1D(("deltaPhi_lep_jet2_E","",100,0,5),"deltaPhi_lep_jet2","weightSSOS_final")
        hist_deltaEta_lep_jet1_E[s] = df_electron[s].Histo1D(("deltaEta_lep_jet1_E","",100,0,5),"deltaEta_lep_jet1","weightSSOS_final")
        hist_deltaEta_lep_jet2_E[s] = df_electron[s].Histo1D(("deltaEta_lep_jet2_E","",100,0,5),"deltaEta_lep_jet2","weightSSOS_final")

        hist_jet_muon_btag_M[s] = df_muon[s].Histo1D(("jet_sv_btag_M","",50,0,1),"jet_muon_btag","weightSSOS_final")
        hist_jet_muon_btag_E[s] = df_electron[s].Histo1D(("jet_sv_btag_E","",50,0,1),"jet_muon_btag","weightSSOS_final")
        hist_jet_notmuon_btag_M[s] = df_muon[s].Histo1D(("jet_notsv_btag_M","",50,0,1),"jet_notmuon_btag","weightSSOS_final")
        hist_jet_notmuon_btag_E[s] = df_electron[s].Histo1D(("jet_notsv_btag_E","",50,0,1),"jet_notmuon_btag","weightSSOS_final")

        hist_jet_bot1_btag_M[s] = df_muon[s].Histo1D(("jet_bot1_btag_M","",50,0,1),"jet_bot1_btag","weightSSOS_final")
        hist_jet_bot1_btag_E[s] = df_electron[s].Histo1D(("jet_bot1_btag_E","",50,0,1),"jet_bot1_btag","weightSSOS_final")
        hist_jet_bot2_btag_M[s] = df_muon[s].Histo1D(("jet_bot2_btag_M","",50,0,1),"jet_bot2_btag","weightSSOS_final")
        hist_jet_bot2_btag_E[s] = df_electron[s].Histo1D(("jet_bot2_btag_E","",50,0,1),"jet_bot2_btag","weightSSOS_final")

        hist_pT_proy_M[s] = df_muon[s].Histo1D(("pT_proy_M","",100,-100,100),"pT_proy","weightSSOS_final")
        hist_pT_proy_E[s] = df_electron[s].Histo1D(("pT_proy_E","",100,-100,100),"pT_proy","weightSSOS_final")

        hist_pT_sum_2J_M[s] = df_muon[s].Histo1D(("pT_sum 2J_M","",100,0,1),"pT_sum_2J","weightSSOS_final")
        hist_pT_sum_2J_E[s] = df_electron[s].Histo1D(("pT_sum 2J_E","",100,0,1),"pT_sum_2J","weightSSOS_final")

########## gen hists

if mode == "mc":
        for s in samples:
               hist_jet_muon_flavourH_M[s] = df_muon[s].Histo1D(("jet_sv_flavourH_M","",6,0,6),"jet_muon_flavourH","weightSSOS_final")
               hist_jet_muon_flavourH_E[s] = df_electron[s].Histo1D(("jet_sv_flavourH_E","",6,0,6),"jet_muon_flavourH","weightSSOS_final")
               hist_jet_notmuon_flavourH_M[s] = df_muon[s].Histo1D(("jet_notsv_flavourH_M","",6,0,6),"jet_notmuon_flavourH","weightSSOS_final")
               hist_jet_notmuon_flavourH_E[s] = df_electron[s].Histo1D(("jet_notsv_flavourH_E","",6,0,6),"jet_notmuon_flavourH","weightSSOS_final")
               hist_jet_muon_flavourP_M[s] = df_muon[s].Histo1D(("jet_sv_flavourP_M","",28,-6,22),"jet_muon_flavourP","weightSSOS_final")
               hist_jet_muon_flavourP_E[s] = df_electron[s].Histo1D(("jet_sv_flavourP_E","",28,-6,22),"jet_muon_flavourP","weightSSOS_final")
               hist_jet_notmuon_flavourP_M[s] = df_muon[s].Histo1D(("jet_notsv_flavourP_M","",28,-6,22),"jet_notmuon_flavourP","weightSSOS_final")
               hist_jet_notmuon_flavourP_E[s] = df_electron[s].Histo1D(("jet_notsv_flavourP_E","",28,-6,22),"jet_notmuon_flavourP","weightSSOS_final")
               hist_br_sf_M[s] = df_muon[s].Histo1D(("Br_weight_sv_M","",100,0.5,1.5),"Br_weight_sv","weightSSOS_final")
               hist_br_sf_E[s] = df_electron[s].Histo1D(("Br_weight_sv_E","",100,0.5,1.5),"Br_weight_sv","weightSSOS_final")
               hist_frag_sf_M[s] = df_muon[s].Histo1D(("Frag_weight_sv_M","",100,0.5,1.5),"Frag_weight_sv","weightSSOS_final")
               hist_frag_sf_E[s] = df_electron[s].Histo1D(("Frag_weight_sv_E","",100,0.5,1.5),"Frag_weight_sv","weightSSOS_final")

#############################
####     DATA SAVING     ####
#############################

if args.notfull:
        for s in samples:
                if args.ssos : path_hist = '/nfs/cms/vazqueze/ttbaranalisis/hists/fromJF/SV/ssos/histfromJF_sv_v1v2SSOS_'+s+'_range_'+str(args.list[0])+'_'+str(args.list[1])+'.root'
                else: path_hist = '/nfs/cms/vazqueze/ttbaranalisis/hists/fromJF/SV/histfromJF_sv_v1v2_'+s+'_range_'+str(args.list[0])+'_'+str(args.list[1])+'.root'
                myfile = TFile( path_hist, 'RECREATE' )

                hist_nJetGood_M[s].Write()
                hist_nJetGood_E[s].Write()
                hist_nSV_M[s].Write()
                hist_nSV_E[s].Write()
                hist_nSV_nobot_M[s].Write()
                hist_nSV_nobot_E[s].Write()
                hist_nLooseLepton_M[s].Write()
                hist_nLooseLepton_E[s].Write()
                hist_nTau_M[s].Write()
                hist_nTau_E[s].Write()
                hist_second_muon_pt_M[s].Write()
                hist_second_muon_pt_E[s].Write()
                hist_second_electron_pt_M[s].Write()
                hist_second_electron_pt_E[s].Write()
                hist_jet_muon_pt_M[s].Write()
                hist_jet_muon_nmu_M[s].Write()
                hist_jet_muon_mass_M[s].Write()
                hist_jet_muon_eta_M[s].Write()
                hist_jet_muon_pt_E[s].Write()
                hist_jet_muon_eta_E[s].Write()
                hist_jet_muon_nmu_E[s].Write()
                hist_jet_muon_mass_E[s].Write()
                hist_jet_notmuon_pt_M[s].Write()
                hist_jet_notmuon_eta_M[s].Write()
                hist_jet_notmuon_nmu_M[s].Write()
                hist_jet_notmuon_mass_M[s].Write()
                hist_jet_notmuon_qgl_M[s].Write()
                hist_jet_notmuon_pt_E[s].Write()
                hist_jet_notmuon_eta_E[s].Write()
                hist_jet_notmuon_nmu_E[s].Write()
                hist_jet_notmuon_nmu_E[s].Write()
                hist_jet_notmuon_mass_E[s].Write()
                hist_jet_notmuon_qgl_E[s].Write()
                hist_jet_bot1_nsv_good_M[s].Write()
                hist_jet_bot1_nsv_good_E[s].Write()
                hist_jet_bot2_nsv_good_M[s].Write()
                hist_jet_bot2_nsv_good_E[s].Write()
                hist_jet_muon_nsv_good_M[s].Write()
                hist_jet_muon_nsv_good_E[s].Write()
                hist_lepton_pt_M[s].Write()
                hist_lepton_eta_M[s].Write()
                hist_lepton_pt_E[s].Write()
                hist_lepton_eta_E[s].Write()
                hist_my_tau_pt_M[s].Write()
                hist_my_tau_pt_E[s].Write()
                hist_my_tau_mass_M[s].Write()
                hist_my_tau_mass_E[s].Write()
                hist_sv_jet_pt_M[s].Write()
                hist_sv_jet_eta_M[s].Write()
                hist_sv_jet_mass_M[s].Write()
                hist_sv_jet_pangle_M[s].Write()
                hist_sv_jet_chi_M[s].Write()
                hist_sv_jet_ndof_M[s].Write()
                hist_sv_jet_pt_E[s].Write()
                hist_sv_jet_eta_E[s].Write()
                hist_sv_jet_mass_E[s].Write()
                hist_sv_jet_pangle_E[s].Write()
                hist_sv_jet_chi_E[s].Write()
                hist_sv_jet_ndof_E[s].Write()
                hist_sv_jet_ntracks_M[s].Write()
                hist_sv_jet_ntracks_E[s].Write()
                hist_sv_jet_charge_M[s].Write()
                hist_sv_jet_charge_E[s].Write()
                hist_sv_jet_relpt_M[s].Write()
                hist_sv_jet_relpt_E[s].Write()
                hist_sv_max_pt_M[s].Write()
                hist_sv_max_eta_M[s].Write()
                hist_sv_max_mass_M[s].Write()
                hist_sv_max_pangle_M[s].Write()
                hist_sv_max_chi_M[s].Write()
                hist_sv_max_ndof_M[s].Write()
                hist_sv_max_pt_E[s].Write()
                hist_sv_max_eta_E[s].Write()
                hist_sv_max_mass_E[s].Write()
                hist_sv_max_pangle_E[s].Write()
                hist_sv_max_chi_E[s].Write()
                hist_sv_max_ndof_E[s].Write()
                hist_sv_max_ntracks_M[s].Write()
                hist_sv_max_ntracks_E[s].Write()
                hist_sv_max_charge_M[s].Write()
                hist_sv_max_charge_E[s].Write()
                hist_InvM_2jets_M[s].Write()
                hist_InvM_2jets_E[s].Write()
                hist_InvM_jetM_lepM[s].Write()
                hist_InvM_jetM_lepE[s].Write()
                hist_InvM_bot_closer_M[s].Write()
                hist_InvM_bot_closer_E[s].Write()
                hist_InvM_bot_farther_M[s].Write()
                hist_InvM_bot_farther_E[s].Write()
                hist_deltaR_jetM_jetNM_M[s].Write()
                hist_deltaR_jetM_jetNM_E[s].Write()
                hist_deltaphi_jetM_jetNM_M[s].Write()
                hist_deltaphi_jetM_jetNM_E[s].Write()
                hist_deltaeta_jetM_jetNM_M[s].Write()
                hist_deltaeta_jetM_jetNM_E[s].Write()
                hist_deltapt_jetM_jetNM_M[s].Write()
                hist_deltapt_jetM_jetNM_E[s].Write()
                hist_MET_M[s].Write()
                hist_MET_E[s].Write()
                hist_tranverse_massM[s].Write()
                hist_tranverse_massE[s].Write()
                hist_tracks_jetM_M[s].Write()
                hist_tracks_jetNM_M[s].Write()
                hist_tracks_jetM_E[s].Write()
                hist_tracks_jetNM_E[s].Write()
                hist_EMN_jetM_M[s].Write()
                hist_EMC_jetM_M[s].Write()
                hist_EMN_jetM_E[s].Write()
                hist_EMC_jetM_E[s].Write()
                hist_EMtotal_jetM_M[s].Write()
                hist_EMtotal_jetM_E[s].Write()
                hist_sv_jet_sigxy_M[s].Write()
                hist_sv_jet_sigxy_E[s].Write()
                hist_sv_jet_sigr_M[s].Write()
                hist_sv_jet_sigr_E[s].Write()
                hist_sv_jet_xy_M[s].Write()
                hist_sv_jet_xy_E[s].Write()
                hist_sv_jet_r_M[s].Write()
                hist_sv_jet_r_E[s].Write()
                hist_sv_max_sigxy_M[s].Write()
                hist_sv_max_sigxy_E[s].Write()
                hist_sv_max_sigr_M[s].Write()
                hist_sv_max_sigr_E[s].Write()
                hist_sv_max_xy_M[s].Write()
                hist_sv_max_xy_E[s].Write()
                hist_sv_max_r_M[s].Write()
                hist_sv_max_r_E[s].Write()
                hist_SSOS_M[s].Write()
                hist_SSOS_E[s].Write()
                hist_pT_sum_M[s].Write()
                hist_pT_sum_E[s].Write()
                hist_pT_product_M[s].Write()
                hist_pT_product_E[s].Write()
                hist_deltaR_lep_2jets_M[s].Write()
                hist_deltaR_lep_2jets_E[s].Write()
                hist_deltaphi_MET_2jets_M[s].Write()
                hist_deltaphi_MET_2jets_E[s].Write()
                hist_deltaphi_lephad_M[s].Write()
                hist_deltaphi_lephad_E[s].Write()
                hist_eta_2jets_M[s].Write()
                hist_eta_2jets_E[s].Write()
                hist_pt_2jets_M[s].Write()
                hist_pt_2jets_E[s].Write()
                hist_pt_Wlep_M[s].Write()
                hist_pt_Wlep_E[s].Write()
                hist_deltaR_lephad_M[s].Write()
                hist_deltaR_lephad_E[s].Write()
                hist_deltaphi_lep_2jets_M[s].Write()
                hist_deltaphi_lep_2jets_E[s].Write()
                hist_deltaeta_lephad_M[s].Write()
                hist_deltaeta_lephad_E[s].Write()
                hist_deltaeta_lep_2jets_M[s].Write()
                hist_deltaeta_lep_2jets_E[s].Write()
                hist_deltapt_lephad_M[s].Write()
                hist_deltapt_lephad_E[s].Write()
                hist_deltapt_lep_2jets_M[s].Write()
                hist_deltapt_lep_2jets_E[s].Write()
                hist_deltaR_lep_jet1_M[s].Write()
                hist_deltaR_lep_jet2_M[s].Write()
                hist_deltaPhi_lep_jet1_M[s].Write()
                hist_deltaPhi_lep_jet2_M[s].Write()
                hist_deltaEta_lep_jet1_M[s].Write()
                hist_deltaEta_lep_jet2_M[s].Write()
                hist_deltaR_lep_jet1_E[s].Write()
                hist_deltaR_lep_jet2_E[s].Write()
                hist_deltaPhi_lep_jet1_E[s].Write()
                hist_deltaPhi_lep_jet2_E[s].Write()
                hist_deltaEta_lep_jet1_E[s].Write()
                hist_deltaEta_lep_jet2_E[s].Write()
                hist_jet_bot1_btag_M[s].Write()
                hist_jet_bot1_btag_E[s].Write()
                hist_jet_bot2_btag_M[s].Write()
                hist_jet_bot2_btag_E[s].Write()
                hist_jet_muon_btag_M[s].Write()
                hist_jet_muon_btag_E[s].Write()
                hist_jet_notmuon_btag_M[s].Write()
                hist_jet_notmuon_btag_E[s].Write()
                hist_pT_proy_M[s].Write()
                hist_pT_proy_E[s].Write()
                hist_pT_sum_2J_M[s].Write()
                hist_pT_sum_2J_E[s].Write()

                if mode=="mc":
                        hist_jet_muon_flavourH_M[s].Write()
                        hist_jet_muon_flavourH_E[s].Write()
                        hist_jet_notmuon_flavourH_M[s].Write()
                        hist_jet_notmuon_flavourH_E[s].Write()
                        hist_jet_muon_flavourP_M[s].Write()
                        hist_jet_muon_flavourP_E[s].Write()
                        hist_jet_notmuon_flavourP_M[s].Write()
                        hist_jet_notmuon_flavourP_E[s].Write()

                myfile.Close()

else:
	for s in samples:
                if args.ssos: path_hist = '/nfs/cms/vazqueze/hists_ttbar/hists/fromJF/SV/ssos/histfromJF_sv_v1v2vBWPSSOS_'+s+'.root'
                else: path_hist = '/nfs/cms/vazqueze/hists_ttbar/hists/fromJF/SV/histfromJF_sv_v1v2vBWP_'+s+'.root'
                myfile = TFile( path_hist, 'RECREATE' )

                hist_nJetGood_M[s].Write()
                hist_nJetGood_E[s].Write()
                hist_nSV_M[s].Write()
                hist_nSV_E[s].Write()
                hist_nSV_nobot_M[s].Write()
                hist_nSV_nobot_E[s].Write()
                hist_nLooseLepton_M[s].Write()
                hist_nLooseLepton_E[s].Write()
                hist_nTau_M[s].Write()
                hist_nTau_E[s].Write()
                hist_second_muon_pt_M[s].Write()
                hist_second_muon_pt_E[s].Write()
                hist_second_electron_pt_M[s].Write()
                hist_second_electron_pt_E[s].Write()
                hist_jet_muon_pt_M[s].Write()
                hist_jet_muon_nmu_M[s].Write()
                hist_jet_muon_mass_M[s].Write()
                hist_jet_muon_eta_M[s].Write()
                hist_jet_muon_pt_E[s].Write()
                hist_jet_muon_eta_E[s].Write()
                hist_jet_muon_nmu_E[s].Write()
                hist_jet_muon_mass_E[s].Write()
                hist_jet_notmuon_pt_M[s].Write()
                hist_jet_notmuon_eta_M[s].Write()
                hist_jet_notmuon_nmu_M[s].Write()
                hist_jet_notmuon_mass_M[s].Write()
                hist_jet_notmuon_qgl_M[s].Write()
                hist_jet_notmuon_pt_E[s].Write()
                hist_jet_notmuon_eta_E[s].Write()
                hist_jet_notmuon_nmu_E[s].Write()
                hist_jet_notmuon_nmu_E[s].Write()
                hist_jet_notmuon_mass_E[s].Write()
                hist_jet_notmuon_qgl_E[s].Write()
                hist_jet_bot1_nsv_good_M[s].Write()
                hist_jet_bot1_nsv_good_E[s].Write()
                hist_jet_bot2_nsv_good_M[s].Write()
                hist_jet_bot2_nsv_good_E[s].Write()
                hist_jet_muon_nsv_good_M[s].Write()
                hist_jet_muon_nsv_good_E[s].Write()
                hist_my_tau_pt_M[s].Write()
                hist_my_tau_pt_E[s].Write()
                hist_my_tau_mass_M[s].Write()
                hist_my_tau_mass_E[s].Write()
                hist_lepton_pt_M[s].Write()
                hist_lepton_eta_M[s].Write()
                hist_lepton_pt_E[s].Write()
                hist_lepton_eta_E[s].Write()
                hist_sv_jet_pt_M[s].Write()
                hist_sv_jet_eta_M[s].Write()
                hist_sv_jet_mass_M[s].Write()
                hist_sv_jet_pangle_M[s].Write()
                hist_sv_jet_chi_M[s].Write()
                hist_sv_jet_ndof_M[s].Write()
                hist_sv_jet_pt_E[s].Write()
                hist_sv_jet_eta_E[s].Write()
                hist_sv_jet_mass_E[s].Write()
                hist_sv_jet_pangle_E[s].Write()
                hist_sv_jet_chi_E[s].Write()
                hist_sv_jet_ndof_E[s].Write()
                hist_sv_jet_ntracks_M[s].Write()
                hist_sv_jet_ntracks_E[s].Write()
                hist_sv_jet_charge_M[s].Write()
                hist_sv_jet_charge_E[s].Write()
                hist_sv_jet_relpt_M[s].Write()
                hist_sv_jet_relpt_E[s].Write()
                hist_sv_max_pt_M[s].Write()
                hist_sv_max_eta_M[s].Write()
                hist_sv_max_mass_M[s].Write()
                hist_sv_max_pangle_M[s].Write()
                hist_sv_max_chi_M[s].Write()
                hist_sv_max_ndof_M[s].Write()
                hist_sv_max_pt_E[s].Write()
                hist_sv_max_eta_E[s].Write()
                hist_sv_max_mass_E[s].Write()
                hist_sv_max_pangle_E[s].Write()
                hist_sv_max_chi_E[s].Write()
                hist_sv_max_ndof_E[s].Write()
                hist_sv_max_ntracks_M[s].Write()
                hist_sv_max_ntracks_E[s].Write()
                hist_sv_max_charge_M[s].Write()
                hist_sv_max_charge_E[s].Write()
                hist_InvM_2jets_M[s].Write()
                hist_InvM_2jets_E[s].Write()
                hist_InvM_jetM_lepM[s].Write()
                hist_InvM_jetM_lepE[s].Write()
                hist_InvM_bot_closer_M[s].Write()
                hist_InvM_bot_closer_E[s].Write()
                hist_InvM_bot_farther_M[s].Write()
                hist_InvM_bot_farther_E[s].Write()
                hist_deltaR_jetM_jetNM_M[s].Write()
                hist_deltaR_jetM_jetNM_E[s].Write()
                hist_deltaphi_jetM_jetNM_M[s].Write()
                hist_deltaphi_jetM_jetNM_E[s].Write()
                hist_deltaeta_jetM_jetNM_M[s].Write()
                hist_deltaeta_jetM_jetNM_E[s].Write()
                hist_deltapt_jetM_jetNM_M[s].Write()
                hist_deltapt_jetM_jetNM_E[s].Write()
                hist_MET_M[s].Write()
                hist_MET_E[s].Write()
                hist_tranverse_massM[s].Write()
                hist_tranverse_massE[s].Write()
                hist_tracks_jetM_M[s].Write()
                hist_tracks_jetNM_M[s].Write()
                hist_tracks_jetM_E[s].Write()
                hist_tracks_jetNM_E[s].Write()
                hist_EMN_jetM_M[s].Write()
                hist_EMC_jetM_M[s].Write()
                hist_EMN_jetM_E[s].Write()
                hist_EMC_jetM_E[s].Write()
                hist_EMtotal_jetM_M[s].Write()
                hist_EMtotal_jetM_E[s].Write()
                hist_sv_jet_sigxy_M[s].Write()
                hist_sv_jet_sigxy_E[s].Write()
                hist_sv_jet_sigr_M[s].Write()
                hist_sv_jet_sigr_E[s].Write()
                hist_sv_jet_xy_M[s].Write()
                hist_sv_jet_xy_E[s].Write()
                hist_sv_jet_r_M[s].Write()
                hist_sv_jet_r_E[s].Write()
                hist_sv_max_sigxy_M[s].Write()
                hist_sv_max_sigxy_E[s].Write()
                hist_sv_max_sigr_M[s].Write()
                hist_sv_max_sigr_E[s].Write()
                hist_sv_max_xy_M[s].Write()
                hist_sv_max_xy_E[s].Write()
                hist_sv_max_r_M[s].Write()
                hist_sv_max_r_E[s].Write()
                hist_SSOS_M[s].Write()
                hist_SSOS_E[s].Write()
                hist_pT_sum_M[s].Write()
                hist_pT_sum_E[s].Write()
                hist_pT_product_M[s].Write()
                hist_pT_product_E[s].Write()
                hist_deltaR_lep_2jets_M[s].Write()                                                                   
                hist_deltaR_lep_2jets_E[s].Write()
                hist_deltaphi_MET_2jets_M[s].Write()
                hist_deltaphi_MET_2jets_E[s].Write()
                hist_deltaphi_lephad_M[s].Write()
                hist_deltaphi_lephad_E[s].Write()
                hist_eta_2jets_M[s].Write()
                hist_eta_2jets_E[s].Write()
                hist_pt_2jets_M[s].Write()
                hist_pt_2jets_E[s].Write()
                hist_pt_Wlep_M[s].Write()
                hist_pt_Wlep_E[s].Write()
                hist_deltaR_lephad_M[s].Write()
                hist_deltaR_lephad_E[s].Write()
                hist_deltaphi_lep_2jets_M[s].Write()
                hist_deltaphi_lep_2jets_E[s].Write()
                hist_deltaeta_lephad_M[s].Write()
                hist_deltaeta_lephad_E[s].Write()
                hist_deltaeta_lep_2jets_M[s].Write()
                hist_deltaeta_lep_2jets_E[s].Write()
                hist_deltapt_lephad_M[s].Write()
                hist_deltapt_lephad_E[s].Write()
                hist_deltapt_lep_2jets_M[s].Write()
                hist_deltapt_lep_2jets_E[s].Write()
                hist_deltaR_lep_jet1_M[s].Write()
                hist_deltaR_lep_jet2_M[s].Write()
                hist_deltaPhi_lep_jet1_M[s].Write()
                hist_deltaPhi_lep_jet2_M[s].Write()
                hist_deltaEta_lep_jet1_M[s].Write()
                hist_deltaEta_lep_jet2_M[s].Write()
                hist_deltaR_lep_jet1_E[s].Write()
                hist_deltaR_lep_jet2_E[s].Write()
                hist_deltaPhi_lep_jet1_E[s].Write()
                hist_deltaPhi_lep_jet2_E[s].Write()
                hist_deltaEta_lep_jet1_E[s].Write()
                hist_deltaEta_lep_jet2_E[s].Write()
                hist_jet_bot1_btag_M[s].Write()
                hist_jet_bot1_btag_E[s].Write()
                hist_jet_bot2_btag_M[s].Write()
                hist_jet_bot2_btag_E[s].Write()
                hist_jet_muon_btag_M[s].Write()
                hist_jet_muon_btag_E[s].Write()
                hist_jet_notmuon_btag_M[s].Write()
                hist_jet_notmuon_btag_E[s].Write()
                hist_pT_proy_M[s].Write()
                hist_pT_proy_E[s].Write()
                hist_pT_sum_2J_M[s].Write()
                hist_pT_sum_2J_E[s].Write()

                if mode=="mc":
                        hist_jet_muon_flavourH_M[s].Write()
                        hist_jet_muon_flavourH_E[s].Write()
                        hist_jet_notmuon_flavourH_M[s].Write()
                        hist_jet_notmuon_flavourH_E[s].Write()
                        hist_jet_muon_flavourP_M[s].Write()
                        hist_jet_muon_flavourP_E[s].Write()
                        hist_jet_notmuon_flavourP_M[s].Write()
                        hist_jet_notmuon_flavourP_E[s].Write()
                        hist_br_sf_M[s].Write()
                        hist_br_sf_E[s].Write()
                        hist_frag_sf_M[s].Write()
                        hist_frag_sf_E[s].Write()

                myfile.Close()

if args.ssos: print('SSOS version')

print('Ended succesfully')

