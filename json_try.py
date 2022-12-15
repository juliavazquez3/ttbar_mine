######### Muon trigger sfs

gInterpreter.Declare("""
      #include <nlohmann/json.hpp>
      using nlohmann::json;
      #include <fstream>
      #include <iostream>
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      using namespace std;
      auto muon_trigger_sf(Vint muon_good, Vfloat pt, Vint eta, string file_json) {
            vector<float> vb;
            std::ifstream fJson(file_json);
            stringstream buffer;
            buffer << fJson.rdbuf();
            auto my_sfs = nlohmann::json::parse(buffer.str());
            //auto my_sfs = my_json["TightISO_TightID_pt_eta"]; // Values for my requirements
            for (unsigned int i=0; i<pt.size(); ++i) {
              if (pt[i] > 20. && pt[i] <= 25.) {
                  if (fabs(eta[i]) <= 0.9) {
                     vb.push_back(my_sfs["pt_abseta_ratio"]["pt:[20.0,25.0]"]["abseta:[0.0,0.9]"]["value"]);
                  } else if (fabs(eta[i]) <= 1.2 && fabs(eta[i]) > 0.9) {
                     vb.push_back(my_sfs["pt_abseta_ratio"]["pt:[20.0,25.0]"]["abseta:[0.9,1.2]"]["value"]);
                  } else if (fabs(eta[i]) <= 2.1 && fabs(eta[i]) > 1.2) {
                     vb.push_back(my_sfs["pt_abseta_ratio"]["pt:[20.0,25.0]"]["abseta:[1.2,2.1]"]["value"]);
                  } else if (fabs(eta[i]) <= 2.4 && fabs(eta[i]) > 2.1) {
                     vb.push_back(my_sfs["pt_abseta_ratio"]["pt:[20.0,25.0]"]["abseta:[2.1,2.4]"]["value"]);
                  } else {
                     vb.push_back(1.);
                  }
              } else if (pt[i] > 25. && pt[i] <= 30.) {
                  if (fabs(eta[i]) <= 0.9) {
                     vb.push_back(my_sfs["pt_abseta_ratio"]["pt:[25.0,30.0]"]["abseta:[0.0,0.9]"]["value"]);
                  } else if (fabs(eta[i]) <= 1.2 && fabs(eta[i]) > 0.9) {
                     vb.push_back(my_sfs["pt_abseta_ratio"]["pt:[25.0,30.0]"]["abseta:[0.9,1.2]"]["value"]);
                  } else if (fabs(eta[i]) <= 2.1 && fabs(eta[i]) > 1.2) {
                     vb.push_back(my_sfs["pt_abseta_ratio"]["pt:[25.0,30.0]"]["abseta:[1.2,2.1]"]["value"]);
                  } else if (fabs(eta[i]) <= 2.4 && fabs(eta[i]) > 2.1) {
                     vb.push_back(my_sfs["pt_abseta_ratio"]["pt:[25.0,30.0]"]["abseta:[2.1,2.4]"]["value"]);
                  }  else {
                     vb.push_back(1.);
                  }
              } else if (pt[i] > 30. && pt[i] <= 40.) {
                  if (fabs(eta[i]) <= 0.9) {
                     vb.push_back(my_sfs["pt_abseta_ratio"]["pt:[30.0,40.0]"]["abseta:[0.0,0.9]"]["value"]);
                  } else if (fabs(eta[i]) <= 1.2 && fabs(eta[i]) > 0.9) {
                     vb.push_back(my_sfs["pt_abseta_ratio"]["pt:[30.0,40.0]"]["abseta:[0.9,1.2]"]["value"]);
                  } else if (fabs(eta[i]) <= 2.1 && fabs(eta[i]) > 1.2) {
                     vb.push_back(my_sfs["pt_abseta_ratio"]["pt:[30.0,40.0]"]["abseta:[1.2,2.1]"]["value"]);
                  } else if (fabs(eta[i]) <= 2.4 && fabs(eta[i]) > 2.1) {
                     vb.push_back(my_sfs["pt_abseta_ratio"]["pt:[30.0,40.0]"]["abseta:[2.1,2.4]"]["value"]);
                  }  else {
                     vb.push_back(1.);
                  }
              } else if (pt[i] > 40. && pt[i] <= 50.) {
                  if (fabs(eta[i]) <= 0.9) {
                     vb.push_back(my_sfs["pt_abseta_ratio"]["pt:[40.0,50.0]"]["abseta:[0.0,0.9]"]["value"]);
                  } else if (fabs(eta[i]) <= 1.2 && fabs(eta[i]) > 0.9) {
                     vb.push_back(my_sfs["pt_abseta_ratio"]["pt:[40.0,50.0]"]["abseta:[0.9,1.2]"]["value"]);
                  } else if (fabs(eta[i]) <= 2.1 && fabs(eta[i]) > 1.2) {
                     vb.push_back(my_sfs["pt_abseta_ratio"]["pt:[40.0,50.0]"]["abseta:[1.2,2.1]"]["value"]);
                  } else if (fabs(eta[i]) <= 2.4 && fabs(eta[i]) > 2.1) {
                     vb.push_back(my_sfs["pt_abseta_ratio"]["pt:[40.0,50.0]"]["abseta:[2.1,2.4]"]["value"]);
                  }  else {
                     vb.push_back(1.);
                  }
              } else if (pt[i] > 50. && pt[i] <= 60.) {
                  if (fabs(eta[i]) <= 0.9) {
                     vb.push_back(my_sfs["pt_abseta_ratio"]["pt:[50.0,60.0]"]["abseta:[0.0,0.9]"]["value"]);
                  } else if (fabs(eta[i]) <= 1.2 && fabs(eta[i]) > 0.9) {
                     vb.push_back(my_sfs["pt_abseta_ratio"]["pt:[50.0,60.0]"]["abseta:[0.9,1.2]"]["value"]);
                  } else if (fabs(eta[i]) <= 2.1 && fabs(eta[i]) > 1.2) {
                     vb.push_back(my_sfs["pt_abseta_ratio"]["pt:[50.0,60.0]"]["abseta:[1.2,2.1]"]["value"]);
                  } else if (fabs(eta[i]) <= 2.4 && fabs(eta[i]) > 2.1) {
                     vb.push_back(my_sfs["pt_abseta_ratio"]["pt:[50.0,60.0]"]["abseta:[2.1,2.4]"]["value"]);
                  }  else {
                     vb.push_back(1.);
                  }
              } else if (pt[i] > 60. && pt[i] <= 120.) {
                  if (fabs(eta[i]) <= 0.9) {
                     vb.push_back(my_sfs["pt_abseta_ratio"]["pt:[60.0,120.0]"]["abseta:[0.0,0.9]"]["value"]);
                  } else if (fabs(eta[i]) <= 1.2 && fabs(eta[i]) > 0.9) {
                     vb.push_back(my_sfs["pt_abseta_ratio"]["pt:[60.0,120.0]"]["abseta:[0.9,1.2]"]["value"]);
                  } else if (fabs(eta[i]) <= 2.1 && fabs(eta[i]) > 1.2) {
                     vb.push_back(my_sfs["pt_abseta_ratio"]["pt:[60.0,120.0]"]["abseta:[1.2,2.1]"]["value"]);
                  } else if (fabs(eta[i]) <= 2.4 && fabs(eta[i]) > 2.1) {
                     vb.push_back(my_sfs["pt_abseta_ratio"]["pt:[60.0,120.0]"]["abseta:[2.1,2.4]"]["value"]);
                  }  else {
                     vb.push_back(1.);
                  }
              } else {
                    vb.push_back(1.);
              }
            }
            return vb;
      };
""")

