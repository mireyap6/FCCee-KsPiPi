#ifndef ZHfunctions_H
#define ZHfunctions_H

#include <cmath>
#include <vector>
#include <math.h>
#include <unordered_map>
#include <iostream>

#include "TLorentzVector.h"
#include "ROOT/RVec.hxx"
#include "edm4hep/ReconstructedParticleData.h"
#include "edm4hep/MCParticleData.h"
#include "edm4hep/ParticleIDData.h"
#include "ReconstructedParticle2MC.h"
#include "ReconstructedTrack.h"

// edm4hep renamed TrackerHitData -> TrackerHit3DData around key4hep 2024-10.
// Older releases (e.g. 2023-11-23, matching the winter2023 Delphes
// production files) only have the old name; alias it so the code below
// compiles unchanged against either release.
#if !__has_include("edm4hep/TrackerHit3DData.h")
namespace edm4hep { using TrackerHit3DData = edm4hep::TrackerHitData; }
#endif


namespace FCCAnalyses { namespace ZHfunctions {


// build the Z resonance based on the available leptons. Returns the best lepton pair compatible with the Z mass and recoil at 125 GeV
// technically, it returns a ReconstructedParticleData object with index 0 the di-lepton system, index and 2 the leptons of the pair
struct resonanceBuilder_mass_recoil {
    float m_resonance_mass;
    float m_recoil_mass;
    float chi2_recoil_frac;
    float ecm;
    bool m_use_MC_Kinematics;
    resonanceBuilder_mass_recoil(float arg_resonance_mass, float arg_recoil_mass, float arg_chi2_recoil_frac, float arg_ecm, bool arg_use_MC_Kinematics);
    Vec_rp operator()(Vec_rp legs, Vec_i recind, Vec_i mcind, Vec_rp reco, Vec_mc mc, Vec_i parents, Vec_i daugthers) ;
};
// Function prototypes
ROOT::VecOps::RVec<double> constraints_Ks(bool tight);
ROOT::VecOps::RVec<double> constraints_Lambda0(bool tight);
ROOT::VecOps::RVec<double> constraints_Gamma(bool tight);
ROOT::VecOps::RVec<double> get_V0candidate(VertexingUtils::FCCAnalysesVertex &V0_vtx, ROOT::VecOps::RVec<edm4hep::TrackState> tr_pair,VertexingUtils::FCCAnalysesVertex PV,bool chi2,double chi2_cut);


bool debug_me = false;
// if particle masses defined in a dedicated file, call those rather than defining here
const double m_pi = 0.13957039; // pi+- mass [GeV]
const double m_p  = 0.93827208; // p+- mass [GeV]
const double m_e  = 0.00051099; // e+- mass [GeV]
const double m_mu = 0.10565837; // mu+- mass [GeV]

resonanceBuilder_mass_recoil::resonanceBuilder_mass_recoil(float arg_resonance_mass, float arg_recoil_mass, float arg_chi2_recoil_frac, float arg_ecm, bool arg_use_MC_Kinematics) {m_resonance_mass = arg_resonance_mass, m_recoil_mass = arg_recoil_mass, chi2_recoil_frac = arg_chi2_recoil_frac, ecm = arg_ecm, m_use_MC_Kinematics = arg_use_MC_Kinematics;}

Vec_rp resonanceBuilder_mass_recoil::resonanceBuilder_mass_recoil::operator()(Vec_rp legs, Vec_i recind, Vec_i mcind, Vec_rp reco, Vec_mc mc, Vec_i parents, Vec_i daugthers) {

    Vec_rp result;
    result.reserve(3);
    std::vector<std::vector<int>> pairs; // for each permutation, add the indices of the muons
    int n = legs.size();
  
    if(n > 1) {
        ROOT::VecOps::RVec<bool> v(n);
        std::fill(v.end() - 2, v.end(), true); // helper variable for permutations
        do {
            std::vector<int> pair;
            rp reso;
            reso.charge = 0;
            TLorentzVector reso_lv; 
            for(int i = 0; i < n; ++i) {
                if(v[i]) {
                    pair.push_back(i);
                    reso.charge += legs[i].charge;
                    TLorentzVector leg_lv;

                    if(m_use_MC_Kinematics) { // MC kinematics
                        int track_index = legs[i].tracks_begin;   // index in the Track array
                        int mc_index = ReconstructedParticle2MC::getTrack2MC_index(track_index, recind, mcind, reco);
                        if (mc_index >= 0 && mc_index < mc.size()) {
                            leg_lv.SetXYZM(mc.at(mc_index).momentum.x, mc.at(mc_index).momentum.y, mc.at(mc_index).momentum.z, mc.at(mc_index).mass);
                        }
                    }
                    else { // reco kinematics
                         leg_lv.SetXYZM(legs[i].momentum.x, legs[i].momentum.y, legs[i].momentum.z, legs[i].mass);
                    }

                    reso_lv += leg_lv;
                }
            }

            if(reso.charge != 0) continue; // neglect non-zero charge pairs
            reso.momentum.x = reso_lv.Px();
            reso.momentum.y = reso_lv.Py();
            reso.momentum.z = reso_lv.Pz();
            reso.mass = reso_lv.M();
            result.emplace_back(reso);
            pairs.push_back(pair);

        } while(std::next_permutation(v.begin(), v.end()));
    }
    else {
        std::cout << "ERROR: resonanceBuilder_mass_recoil, at least two leptons required." << std::endl;
        exit(1);
    }
  
    if(result.size() > 1) {
  
        Vec_rp bestReso;
        
        int idx_min = -1;
        float d_min = 9e9;
        for (int i = 0; i < result.size(); ++i) {
            
            // calculate recoil
            auto recoil_p4 = TLorentzVector(0, 0, 0, ecm);
            TLorentzVector tv1;
            tv1.SetXYZM(result.at(i).momentum.x, result.at(i).momentum.y, result.at(i).momentum.z, result.at(i).mass);
            recoil_p4 -= tv1;
      
            auto recoil_fcc = edm4hep::ReconstructedParticleData();
            recoil_fcc.momentum.x = recoil_p4.Px();
            recoil_fcc.momentum.y = recoil_p4.Py();
            recoil_fcc.momentum.z = recoil_p4.Pz();
            recoil_fcc.mass = recoil_p4.M();
            
            TLorentzVector tg;
            tg.SetXYZM(result.at(i).momentum.x, result.at(i).momentum.y, result.at(i).momentum.z, result.at(i).mass);
        
            float boost = tg.P();
            float mass = std::pow(result.at(i).mass - m_resonance_mass, 2); // mass
            float rec = std::pow(recoil_fcc.mass - m_recoil_mass, 2); // recoil
            float d = (1.0-chi2_recoil_frac)*mass + chi2_recoil_frac*rec;
            
            if(d < d_min) {
                d_min = d;
                idx_min = i;
            }

     
        }
        if(idx_min > -1) { 
            bestReso.push_back(result.at(idx_min));
            auto & l1 = legs[pairs[idx_min][0]];
            auto & l2 = legs[pairs[idx_min][1]];
            bestReso.emplace_back(l1);
            bestReso.emplace_back(l2);
        }
        else {
            std::cout << "ERROR: resonanceBuilder_mass_recoil, no mininum found." << std::endl;
            exit(1);
        }
        return bestReso;
    }
    else {
        auto & l1 = legs[0];
        auto & l2 = legs[1];
        result.emplace_back(l1);
        result.emplace_back(l2);
        return result;
    }
}    




struct sel_iso {
    sel_iso(float arg_max_iso);
    float m_max_iso = .25;
    Vec_rp operator() (Vec_rp in, Vec_f iso);
  };

sel_iso::sel_iso(float arg_max_iso) : m_max_iso(arg_max_iso) {};
ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>  sel_iso::operator() (Vec_rp in, Vec_f iso) {
    Vec_rp result;
    result.reserve(in.size());
    for (size_t i = 0; i < in.size(); ++i) {
        auto & p = in[i];
        if (iso[i] < m_max_iso) {
            result.emplace_back(p);
        }
    }
    return result;
}

 
// compute the cone isolation for reco particles
struct coneIsolation {

    coneIsolation(float arg_dr_min, float arg_dr_max);
    double deltaR(double eta1, double phi1, double eta2, double phi2) { return TMath::Sqrt(TMath::Power(eta1-eta2, 2) + (TMath::Power(phi1-phi2, 2))); };

    float dr_min = 0;
    float dr_max = 0.4;
    Vec_f operator() (Vec_rp in, Vec_rp rps) ;
};

coneIsolation::coneIsolation(float arg_dr_min, float arg_dr_max) : dr_min(arg_dr_min), dr_max( arg_dr_max ) { };
Vec_f coneIsolation::coneIsolation::operator() (Vec_rp in, Vec_rp rps) {
  
    Vec_f result;
    result.reserve(in.size());

    std::vector<ROOT::Math::PxPyPzEVector> lv_reco;
    std::vector<ROOT::Math::PxPyPzEVector> lv_charged;
    std::vector<ROOT::Math::PxPyPzEVector> lv_neutral;

    for(size_t i = 0; i < rps.size(); ++i) {

        ROOT::Math::PxPyPzEVector tlv;
        tlv.SetPxPyPzE(rps.at(i).momentum.x, rps.at(i).momentum.y, rps.at(i).momentum.z, rps.at(i).energy);
        
        if(rps.at(i).charge == 0) lv_neutral.push_back(tlv);
        else lv_charged.push_back(tlv);
    }
    
    for(size_t i = 0; i < in.size(); ++i) {

        ROOT::Math::PxPyPzEVector tlv;
        tlv.SetPxPyPzE(in.at(i).momentum.x, in.at(i).momentum.y, in.at(i).momentum.z, in.at(i).energy);
        lv_reco.push_back(tlv);
    }

    
    // compute the isolation (see https://github.com/delphes/delphes/blob/master/modules/Isolation.cc#L154) 
    for (auto & lv_reco_ : lv_reco) {
    
        double sumNeutral = 0.0;
        double sumCharged = 0.0;
    
        // charged
        for (auto & lv_charged_ : lv_charged) {
    
            double dr = coneIsolation::deltaR(lv_reco_.Eta(), lv_reco_.Phi(), lv_charged_.Eta(), lv_charged_.Phi());
            if(dr > dr_min && dr < dr_max) sumCharged += lv_charged_.P();
        }
        
        // neutral
        for (auto & lv_neutral_ : lv_neutral) {
    
            double dr = coneIsolation::deltaR(lv_reco_.Eta(), lv_reco_.Phi(), lv_neutral_.Eta(), lv_neutral_.Phi());
            if(dr > dr_min && dr < dr_max) sumNeutral += lv_neutral_.P();
        }
        
        double sum = sumCharged + sumNeutral;
        double ratio= sum / lv_reco_.P();
        result.emplace_back(ratio);
    }
    return result;
}
 
 
 

ROOT::VecOps::RVec<float> get_MC_Vertex_mass(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertexMC> vertex,
                                                         ROOT::VecOps::RVec<edm4hep::MCParticleData> mc){
  ROOT::VecOps::RVec<float> result;
  for (auto &p:vertex){
    TLorentzVector tlv; 
    for (size_t i = 0; i < p.mc_ind.size(); ++i){   
	TLorentzVector tmp_tlv; 
	tmp_tlv.SetXYZM(mc.at(p.mc_ind.at(i)).momentum.x, mc.at(p.mc_ind.at(i)).momentum.y, mc.at(p.mc_ind.at(i)).momentum.z, mc.at(p.mc_ind.at(i)).mass);
        tlv += tmp_tlv;
    }
    for (size_t i = 0; i < p.mc_indneutral.size(); ++i) {
	TLorentzVector tmp_tlv; 
	tmp_tlv.SetXYZM(mc.at(p.mc_indneutral.at(i)).momentum.x, mc.at(p.mc_indneutral.at(i)).momentum.y, mc.at(p.mc_indneutral.at(i)).momentum.z, mc.at(p.mc_indneutral.at(i)).mass);
        tlv += tmp_tlv;
    }	
    result.push_back(tlv.M());
  }
  return result;
}

ROOT::VecOps::RVec<float> get_MC_Vertex_p(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertexMC> vertex,
                                                         ROOT::VecOps::RVec<edm4hep::MCParticleData> mc){
  ROOT::VecOps::RVec<float> result;
  for (auto &p:vertex){
    TLorentzVector tlv;
    for (size_t i = 0; i < p.mc_ind.size(); ++i){
        TLorentzVector tmp_tlv;
        tmp_tlv.SetXYZM(mc.at(p.mc_ind.at(i)).momentum.x, mc.at(p.mc_ind.at(i)).momentum.y, mc.at(p.mc_ind.at(i)).momentum.z, mc.at(p.mc_ind.at(i)).mass);
        tlv += tmp_tlv;
    }
    for (size_t i = 0; i < p.mc_indneutral.size(); ++i) {
        TLorentzVector tmp_tlv;
        tmp_tlv.SetXYZM(mc.at(p.mc_indneutral.at(i)).momentum.x, mc.at(p.mc_indneutral.at(i)).momentum.y, mc.at(p.mc_indneutral.at(i)).momentum.z, mc.at(p.mc_indneutral.at(i)).mass);
        tlv += tmp_tlv;
    }
    result.push_back(tlv.P());
  }
  return result;
}

ROOT::VecOps::RVec<float> get_MC_Vertex_pt(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertexMC> vertex,
                                                         ROOT::VecOps::RVec<edm4hep::MCParticleData> mc){
  ROOT::VecOps::RVec<float> result;
  for (auto &p:vertex){
    TLorentzVector tlv;
    for (size_t i = 0; i < p.mc_ind.size(); ++i){
        TLorentzVector tmp_tlv;
        tmp_tlv.SetXYZM(mc.at(p.mc_ind.at(i)).momentum.x, mc.at(p.mc_ind.at(i)).momentum.y, mc.at(p.mc_ind.at(i)).momentum.z, mc.at(p.mc_ind.at(i)).mass);
        tlv += tmp_tlv;
    }
    for (size_t i = 0; i < p.mc_indneutral.size(); ++i) {
        TLorentzVector tmp_tlv;
        tmp_tlv.SetXYZM(mc.at(p.mc_indneutral.at(i)).momentum.x, mc.at(p.mc_indneutral.at(i)).momentum.y, mc.at(p.mc_indneutral.at(i)).momentum.z, mc.at(p.mc_indneutral.at(i)).mass);
        tlv += tmp_tlv;
    }
    result.push_back(tlv.Pt());
  }
  return result;
}

ROOT::VecOps::RVec<float> get_MC_Vertex_eta(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertexMC> vertex,
                                            ROOT::VecOps::RVec<edm4hep::MCParticleData> mc){
  ROOT::VecOps::RVec<float> result;
  for (auto &p : vertex){
    TLorentzVector tlv;

    for (size_t i = 0; i < p.mc_ind.size(); ++i){
        TLorentzVector tmp_tlv;
        auto& part = mc.at(p.mc_ind.at(i));
        tmp_tlv.SetXYZM(part.momentum.x, part.momentum.y, part.momentum.z, part.mass);
        tlv += tmp_tlv;
    }

    for (size_t i = 0; i < p.mc_indneutral.size(); ++i) {
        TLorentzVector tmp_tlv;
        auto& part = mc.at(p.mc_indneutral.at(i));
        tmp_tlv.SetXYZM(part.momentum.x, part.momentum.y, part.momentum.z, part.mass);
        tlv += tmp_tlv;
    }
    
    if (tlv.P() > 0) {
        result.push_back(tlv.Eta());
    } else {
        result.push_back(-999.0);
    }
  }
  return result;
}

ROOT::VecOps::RVec<float> get_MC_Vertex_phi(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertexMC> vertex,
                                            ROOT::VecOps::RVec<edm4hep::MCParticleData> mc) {
  ROOT::VecOps::RVec<float> result;
  for (auto &p : vertex) {
    TLorentzVector tlv;

    for (size_t i = 0; i < p.mc_ind.size(); ++i) {
        TLorentzVector tmp_tlv;
        auto& part = mc.at(p.mc_ind.at(i));
        tmp_tlv.SetXYZM(part.momentum.x, part.momentum.y, part.momentum.z, part.mass);
        tlv += tmp_tlv;
    }

    for (size_t i = 0; i < p.mc_indneutral.size(); ++i) {
        TLorentzVector tmp_tlv;
        auto& part = mc.at(p.mc_indneutral.at(i));
        tmp_tlv.SetXYZM(part.momentum.x, part.momentum.y, part.momentum.z, part.mass);
        tlv += tmp_tlv;
    }
    
    if (tlv.Pt() > 0) {
        result.push_back(tlv.Phi());
    } else {
        result.push_back(-999.0);
    }
  }
  return result;
};


ROOT::VecOps::RVec<float> get_Vertex_pt(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex,
                                                   ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco){

  ROOT::VecOps::RVec<float> result;
  for (auto &p:vertex){
    ROOT::VecOps::RVec<int> reco_ind = p.reco_ind;
    TLorentzVector tmp_tlv  = myUtils::build_tlv(reco, reco_ind);
    result.push_back(tmp_tlv.Pt());
  }
  return result;
}

ROOT::VecOps::RVec<float> get_Vertex_eta(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex,
                                                   ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco){

  ROOT::VecOps::RVec<float> result;
  for (auto &p:vertex){
    ROOT::VecOps::RVec<int> reco_ind = p.reco_ind;
    TLorentzVector tmp_tlv  = myUtils::build_tlv(reco, reco_ind);
    result.push_back(tmp_tlv.Eta());
  }
  return result;
}

ROOT::VecOps::RVec<float> get_Vertex_phi(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex,
                                         ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco) {
  ROOT::VecOps::RVec<float> result;
  for (auto &p : vertex) {
    ROOT::VecOps::RVec<int> reco_ind = p.reco_ind;
    TLorentzVector tmp_tlv = myUtils::build_tlv(reco, reco_ind);
    
    if (tmp_tlv.Pt() > 0) {
        result.push_back(tmp_tlv.Phi());
    } else {
        result.push_back(-999.0);
    }
  }
  return result;
};

ROOT::VecOps::RVec<float> get_Vertex_p(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex,
                                       ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco) {
  ROOT::VecOps::RVec<float> result;
  for (auto &p : vertex) {
    ROOT::VecOps::RVec<int> reco_ind = p.reco_ind;
    TLorentzVector tmp_tlv = myUtils::build_tlv(reco, reco_ind);
    result.push_back(tmp_tlv.P());
  }
  return result;
};

ROOT::VecOps::RVec<float> get_Vertex_px(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex,
                                       ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco) {
  ROOT::VecOps::RVec<float> result;
  for (auto &p : vertex) {
    ROOT::VecOps::RVec<int> reco_ind = p.reco_ind;
    TLorentzVector tmp_tlv = myUtils::build_tlv(reco, reco_ind);
    result.push_back(tmp_tlv.Px());
  }
  return result;
};

ROOT::VecOps::RVec<float> get_Vertex_py(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex,
                                        ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco) {
  ROOT::VecOps::RVec<float> result;
  for (auto &p : vertex) {
    ROOT::VecOps::RVec<int> reco_ind = p.reco_ind;
    TLorentzVector tmp_tlv = myUtils::build_tlv(reco, reco_ind);
    result.push_back(tmp_tlv.Py());
  }
  return result;
};

ROOT::VecOps::RVec<float> get_Vertex_pz(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex,
                                        ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco) {
  ROOT::VecOps::RVec<float> result;
  for (auto &p : vertex) {
    ROOT::VecOps::RVec<int> reco_ind = p.reco_ind;
    TLorentzVector tmp_tlv = myUtils::build_tlv(reco, reco_ind);
    result.push_back(tmp_tlv.Pz());
  }
  return result;
};

ROOT::VecOps::RVec<float> get_MC_KS_mass_mumu(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertexMC> vertex,
                                              ROOT::VecOps::RVec<edm4hep::MCParticleData> mc) {
  ROOT::VecOps::RVec<float> result;
  double m_mu = 0.105658; 

  for (auto &p : vertex) {
    TLorentzVector tlv;
    
    // For charged particles
    for (size_t i = 0; i < p.mc_ind.size(); ++i) {
      TLorentzVector tmp_tlv;
      auto& part = mc.at(p.mc_ind.at(i));
      double p2 = part.momentum.x * part.momentum.x + part.momentum.y * part.momentum.y + part.momentum.z * part.momentum.z;
      //Muon mass hypothesis: E = sqrt(p^2 + m_mu^2)
      tmp_tlv.SetXYZM(part.momentum.x, part.momentum.y, part.momentum.z, m_mu);
      tlv += tmp_tlv;
    }
    
    // For neutral particles
    for (size_t i = 0; i < p.mc_indneutral.size(); ++i) {
      TLorentzVector tmp_tlv;
      auto& part = mc.at(p.mc_indneutral.at(i));
      tmp_tlv.SetXYZM(part.momentum.x, part.momentum.y, part.momentum.z, m_mu);
      tlv += tmp_tlv;
    }
    
    result.push_back(tlv.M());
  }
  return result;
}

ROOT::VecOps::RVec<float> get_MC_KS_mass_mumu_smeared(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertexMC> vertex,
                                                     ROOT::VecOps::RVec<edm4hep::MCParticleData> mc) {
  ROOT::VecOps::RVec<float> result;
  double m_mu = 0.105658; 
  double resolution = 0.001; // 0.1% resolución

  // random number
  static TRandom3 gen(42); 

  for (auto &p : vertex) {
    TLorentzVector tlv;
    
    // 1. charged particles
    for (size_t i = 0; i < p.mc_ind.size(); ++i) {
      TLorentzVector tmp_tlv;
      auto& part = mc.at(p.mc_ind.at(i));
      
      double smearFactor = gen.Gaus(1.0, resolution);
      
      // apply smearing to momentum components
      tmp_tlv.SetXYZM(part.momentum.x * smearFactor, 
                      part.momentum.y * smearFactor, 
                      part.momentum.z * smearFactor, 
                      m_mu);
      tlv += tmp_tlv;
    }
    
    // 2. neutral particles
    for (size_t i = 0; i < p.mc_indneutral.size(); ++i) {
      TLorentzVector tmp_tlv;
      auto& part = mc.at(p.mc_indneutral.at(i));
      
      double smearFactor = gen.Gaus(1.0, resolution);
      
      tmp_tlv.SetXYZM(part.momentum.x * smearFactor, 
                      part.momentum.y * smearFactor, 
                      part.momentum.z * smearFactor, 
                      m_mu);
      tlv += tmp_tlv;
    }
    
    // reconstructed mass
    result.push_back(tlv.M());
  }
  return result;
}

ROOT::VecOps::RVec<float> get_KS_mass_mumu(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertices) {
    ROOT::VecOps::RVec<float> result;
    double m_mu = 0.105658; // muon mass

    for (auto &vertex : vertices) {
        // get momenta updated at the vertex point
        ROOT::VecOps::RVec<TVector3> p_tracks = vertex.updated_track_momentum_at_vertex;
        
        TLorentzVector k_short_tlv(0, 0, 0, 0);

        // loop over tracks associated to the vertex and build the K_S candidate 4-momentum using the muon mass hypothesis
        for (size_t i = 0; i < p_tracks.size(); ++i) {
            TLorentzVector muon_tlv;
            muon_tlv.SetXYZM(p_tracks[i].X(), p_tracks[i].Y(), p_tracks[i].Z(), m_mu);
            k_short_tlv += muon_tlv;
        }

        // store the invariant mass of the K_S candidate
        result.push_back(k_short_tlv.M());
    }
    return result;
}

// get the momentum of the tracks at the vertex point
ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> get_momentum_updated_at_vertex(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertices) {
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> result;

    for (auto &vertex : vertices) {
        ROOT::VecOps::RVec<TVector3> p_tracks = vertex.updated_track_momentum_at_vertex;
        ROOT::VecOps::RVec<float> magnitudes;
        
        for (size_t i = 0; i < p_tracks.size(); ++i) {
            magnitudes.push_back(p_tracks[i].Mag());
        }
        result.push_back(magnitudes);
    }
    return result;
}

// get the x component of the momentum of the tracks at the vertex point
ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> get_px_updated_at_vertex(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertices) {
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> result;

    for (auto &vertex : vertices) {
        ROOT::VecOps::RVec<TVector3> p_tracks = vertex.updated_track_momentum_at_vertex;
        ROOT::VecOps::RVec<float> px_values;
        
        for (size_t i = 0; i < p_tracks.size(); ++i) {
            px_values.push_back(p_tracks[i].X());
        }
        result.push_back(px_values);
    }
    return result;
}

// get the y component of the momentum of the tracks at the vertex point
ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> get_py_updated_at_vertex(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertices) {
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> result;

    for (auto &vertex : vertices) {
        ROOT::VecOps::RVec<TVector3> p_tracks = vertex.updated_track_momentum_at_vertex;
        ROOT::VecOps::RVec<float> py_values;
        
        for (size_t i = 0; i < p_tracks.size(); ++i) {
            py_values.push_back(p_tracks[i].Y());
        }
        result.push_back(py_values);
    }
    return result;
}

// get the z component of the momentum of the tracks at the vertex point
ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> get_pz_updated_at_vertex(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertices) {
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> result;

    for (auto &vertex : vertices) {
        ROOT::VecOps::RVec<TVector3> p_tracks = vertex.updated_track_momentum_at_vertex;
        ROOT::VecOps::RVec<float> pz_values;
        
        for (size_t i = 0; i < p_tracks.size(); ++i) {
            pz_values.push_back(p_tracks[i].Z());
        }
        result.push_back(pz_values);
    }
    return result;
}

// omega, tanLambda, and their diagonal covariance elements (covMatrix[5],
// covMatrix[14]) from the ORIGINAL (pre-vertex-fit) TrackState of each track
// that formed each V0 vertex, indexed via reco_ind -> tracks_begin, reading
// the raw EDM4hep track parameters and their covariance instead of the
// vertex-refit momentum. This is what
// Delphes' TrackCovariance module actually assigned to that specific track,
// so it can be compared directly against the analytic IDEA resolution grid
// (track_resolution/compute_track_resolution_grid.cpp) for muons vs pions.
ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> get_V0_tracks_omega(
        ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertices,
        ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recoparticles,
        ROOT::VecOps::RVec<edm4hep::TrackState> tracks) {
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> result;
    for (auto &vertex : vertices) {
        ROOT::VecOps::RVec<float> vals;
        for (int reco_idx : vertex.reco_ind) {
            int trk_idx = recoparticles[reco_idx].tracks_begin;
            vals.push_back(tracks[trk_idx].omega);
        }
        result.push_back(vals);
    }
    return result;
}

ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> get_V0_tracks_tanLambda(
        ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertices,
        ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recoparticles,
        ROOT::VecOps::RVec<edm4hep::TrackState> tracks) {
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> result;
    for (auto &vertex : vertices) {
        ROOT::VecOps::RVec<float> vals;
        for (int reco_idx : vertex.reco_ind) {
            int trk_idx = recoparticles[reco_idx].tracks_begin;
            vals.push_back(tracks[trk_idx].tanLambda);
        }
        result.push_back(vals);
    }
    return result;
}

ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> get_V0_tracks_omega_cov(
        ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertices,
        ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recoparticles,
        ROOT::VecOps::RVec<edm4hep::TrackState> tracks) {
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> result;
    for (auto &vertex : vertices) {
        ROOT::VecOps::RVec<float> vals;
        for (int reco_idx : vertex.reco_ind) {
            int trk_idx = recoparticles[reco_idx].tracks_begin;
            vals.push_back(tracks[trk_idx].covMatrix[5]);
        }
        result.push_back(vals);
    }
    return result;
}

ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> get_V0_tracks_tanLambda_cov(
        ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertices,
        ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recoparticles,
        ROOT::VecOps::RVec<edm4hep::TrackState> tracks) {
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> result;
    for (auto &vertex : vertices) {
        ROOT::VecOps::RVec<float> vals;
        for (int reco_idx : vertex.reco_ind) {
            int trk_idx = recoparticles[reco_idx].tracks_begin;
            vals.push_back(tracks[trk_idx].covMatrix[14]);
        }
        result.push_back(vals);
    }
    return result;
}

// D0 and Z0 (track impact parameters w.r.t. origin) for each V0 daughter leg.
// covMatrix uses 6×6 lower-triangle row-major packing (21 elements):
//   σ²(D0)  = covMatrix[0]   (row 0, col 0)
//   σ²(Z0)  = covMatrix[9]   (row 3, col 3)
// Consistent with covMatrix[5]=σ²(omega) and covMatrix[14]=σ²(tanLambda).
ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> get_V0_tracks_D0(
        ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertices,
        ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recoparticles,
        ROOT::VecOps::RVec<edm4hep::TrackState> tracks) {
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> result;
    for (auto &vertex : vertices) {
        ROOT::VecOps::RVec<float> vals;
        for (int reco_idx : vertex.reco_ind) {
            int trk_idx = recoparticles[reco_idx].tracks_begin;
            vals.push_back(tracks[trk_idx].D0);
        }
        result.push_back(vals);
    }
    return result;
}

ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> get_V0_tracks_D0_cov(
        ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertices,
        ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recoparticles,
        ROOT::VecOps::RVec<edm4hep::TrackState> tracks) {
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> result;
    for (auto &vertex : vertices) {
        ROOT::VecOps::RVec<float> vals;
        for (int reco_idx : vertex.reco_ind) {
            int trk_idx = recoparticles[reco_idx].tracks_begin;
            vals.push_back(tracks[trk_idx].covMatrix[0]);
        }
        result.push_back(vals);
    }
    return result;
}

ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> get_V0_tracks_Z0(
        ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertices,
        ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recoparticles,
        ROOT::VecOps::RVec<edm4hep::TrackState> tracks) {
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> result;
    for (auto &vertex : vertices) {
        ROOT::VecOps::RVec<float> vals;
        for (int reco_idx : vertex.reco_ind) {
            int trk_idx = recoparticles[reco_idx].tracks_begin;
            vals.push_back(tracks[trk_idx].Z0);
        }
        result.push_back(vals);
    }
    return result;
}

ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> get_V0_tracks_Z0_cov(
        ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertices,
        ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recoparticles,
        ROOT::VecOps::RVec<edm4hep::TrackState> tracks) {
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> result;
    for (auto &vertex : vertices) {
        ROOT::VecOps::RVec<float> vals;
        for (int reco_idx : vertex.reco_ind) {
            int trk_idx = recoparticles[reco_idx].tracks_begin;
            vals.push_back(tracks[trk_idx].covMatrix[9]);
        }
        result.push_back(vals);
    }
    return result;
}

// cov(omega, tanLambda) cross term, needed to propagate sigma(p) from the
// track's (omega, tanLambda) covariance block: row 4 (tanLambda), col 2
// (omega) of the 6x6 lower triangle -> index 4*5/2 + 2 = 12.
ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> get_V0_tracks_omega_tanLambda_cov(
        ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertices,
        ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recoparticles,
        ROOT::VecOps::RVec<edm4hep::TrackState> tracks) {
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> result;
    for (auto &vertex : vertices) {
        ROOT::VecOps::RVec<float> vals;
        for (int reco_idx : vertex.reco_ind) {
            int trk_idx = recoparticles[reco_idx].tracks_begin;
            vals.push_back(tracks[trk_idx].covMatrix[12]);
        }
        result.push_back(vals);
    }
    return result;
}

// Reconstructed impact parameter of each V0 (K_S) candidate w.r.t. the fitted
// primary vertex. The K_S is neutral and leaves no track of its own, so this
// is computed purely geometrically: the V0's flight line is a straight line
// through the fitted secondary vertex, along the direction of the summed
// daughter momentum (vertex.updated_track_momentum_at_vertex) -- no curvature,
// since a neutral particle isn't bent by the field. IP is the distance of
// closest approach of that line to the primary vertex, i.e. the component of
// (SV - PV) perpendicular to the flight direction. Same geometric definition
// as a charged track's D0, just applied to a straight line instead of a helix.
ROOT::VecOps::RVec<float> get_V0_impact_parameter(
        ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertices,
        edm4hep::VertexData primaryVertex) {
    ROOT::VecOps::RVec<float> result;
    TVector3 pv(primaryVertex.position.x, primaryVertex.position.y, primaryVertex.position.z);
    for (auto &vertex : vertices) {
        TVector3 sv(vertex.vertex.position.x, vertex.vertex.position.y, vertex.vertex.position.z);
        TVector3 p(0., 0., 0.);
        for (auto &leg_p : vertex.updated_track_momentum_at_vertex) p += leg_p;
        if (p.Mag() < 1e-9) { result.push_back(-1.); continue; }
        TVector3 phat = p.Unit();
        TVector3 d = sv - pv;
        TVector3 perp = d - d.Dot(phat) * phat;
        result.push_back(perp.Mag());
    }
    return result;
}

// Same as above, but only the transverse (x-y) component of the impact
// parameter -- the conventional d0-like definition, if you want the value
// directly comparable to a charged track's D0 rather than the full 3D miss
// distance.
ROOT::VecOps::RVec<float> get_V0_impact_parameter_xy(
        ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertices,
        edm4hep::VertexData primaryVertex) {
    ROOT::VecOps::RVec<float> result;
    TVector3 pv(primaryVertex.position.x, primaryVertex.position.y, 0.);
    for (auto &vertex : vertices) {
        TVector3 sv(vertex.vertex.position.x, vertex.vertex.position.y, 0.);
        TVector3 p(0., 0., 0.);
        for (auto &leg_p : vertex.updated_track_momentum_at_vertex) p += leg_p;
        double pt = std::sqrt(p.X()*p.X() + p.Y()*p.Y());
        if (pt < 1e-9) { result.push_back(-1.); continue; }
        TVector3 phat(p.X()/pt, p.Y()/pt, 0.);
        TVector3 d = sv - pv;
        TVector3 perp = d - d.Dot(phat) * phat;
        result.push_back(perp.Mag());
    }
    return result;
}

// TOF (ps, hit at the calo entrance) of the tracks that formed each V0 vertex
ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> get_V0_tracks_TOF(
        ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertices,
        ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recoparticles,
        ROOT::VecOps::RVec<edm4hep::TrackData> trackdata,          // EFlowTrack
        ROOT::VecOps::RVec<edm4hep::TrackerHit3DData> trackerhits // TrackerHits
        ) {
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> result;
    for (auto &vertex : vertices) {
        ROOT::VecOps::RVec<int> trk_idx;
        for (int reco_idx : vertex.reco_ind) {
            trk_idx.push_back(recoparticles[reco_idx].tracks_begin);
        }
        result.push_back(FCCAnalyses::ReconstructedTrack::tracks_TOF(trk_idx, trackdata, trackerhits));
    }
    return result;
}

// path length (mm) of the tracks that formed each V0 vertex
ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> get_V0_tracks_length(
        ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertices,
        ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recoparticles,
        ROOT::VecOps::RVec<float> track_length // EFlowTrack_L
        ) {
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> result;
    for (auto &vertex : vertices) {
        ROOT::VecOps::RVec<int> trk_idx;
        for (int reco_idx : vertex.reco_ind) {
            trk_idx.push_back(recoparticles[reco_idx].tracks_begin);
        }
        result.push_back(FCCAnalyses::ReconstructedTrack::tracks_length(trk_idx, track_length));
    }
    return result;
}

// Raw (pre-vertex-fit) momentum of the tracks that formed each V0 vertex --
// i.e. ReconstructedParticles' own momentum, NOT the vertex-constrained
// refit momentum (updated_track_momentum_at_vertex). This is the momentum
// consistent with the TOF and length above, which are also derived from the
// original, unrefit track/hits -- the vertex-fit-updated momentum comes from
// a joint 2-track constrained refit (VertexFitterSimple.cc, GetNewPar/
// GetMomentum) that can shift the momentum magnitude itself, not just
// direction, so it should NOT be paired with TOF/length for a mass
// calculation.
ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> get_V0_tracks_p_raw(
        ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertices,
        ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recoparticles
        ) {
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> result;
    for (auto &vertex : vertices) {
        ROOT::VecOps::RVec<float> p_values;
        for (int reco_idx : vertex.reco_ind) {
            auto &mom = recoparticles[reco_idx].momentum;
            p_values.push_back(std::sqrt(mom.x*mom.x + mom.y*mom.y + mom.z*mom.z));
        }
        result.push_back(p_values);
    }
    return result;
}

// Ks flight time, per vertex, from its production point to its decay vertex
// (= the production point of the two V0 daughters), in ps.
//
// TrackerHits.time (and therefore V0_tracks_TOF_ps) is measured from the
// collision itself (t=0 at the IP), NOT from when the daughter track was
// created. Since the Ks travels some distance before decaying, that raw TOF
// includes the Ks's own flight time before the daughter (e.g. the muon) even
// exists. Subtracting this t0 from V0_tracks_TOF_ps before computing beta
// corrects for that.
//
// Ks production point is assumed to be the origin (the IP): a good
// approximation for Ks produced directly in light-quark (here Z->ss)
// hadronization, with no long-lived heavy-flavour decay chain upstream to
// displace it further. Ks momentum is reconstructed as the vector sum of the
// two daughter momenta at the vertex (V0_tracks_px/py/pz_at_vertex).
ROOT::VecOps::RVec<float> get_V0_ks_t0(
        ROOT::VecOps::RVec<float> Vertex_x, // mm, one entry per vertex
        ROOT::VecOps::RVec<float> Vertex_y,
        ROOT::VecOps::RVec<float> Vertex_z,
        ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> px_at_vertex, // GeV, 2 entries per vertex
        ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> py_at_vertex,
        ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> pz_at_vertex
        ) {
    const float c_mm_per_ps = 0.299792458f;
    const float m_KS = 0.497611f; // GeV
    ROOT::VecOps::RVec<float> result;
    for (size_t iv = 0; iv < Vertex_x.size(); ++iv) {
        float d = std::sqrt(Vertex_x[iv]*Vertex_x[iv] + Vertex_y[iv]*Vertex_y[iv]
                             + Vertex_z[iv]*Vertex_z[iv]); // mm, IP -> decay vertex

        float px = 0.f, py = 0.f, pz = 0.f;
        for (size_t k = 0; k < px_at_vertex[iv].size(); ++k) {
            px += px_at_vertex[iv][k];
            py += py_at_vertex[iv][k];
            pz += pz_at_vertex[iv][k];
        }
        float p_KS = std::sqrt(px*px + py*py + pz*pz);
        float beta_KS = p_KS / std::sqrt(p_KS*p_KS + m_KS*m_KS);

        result.push_back(beta_KS > 0.f ? d / (c_mm_per_ps * beta_KS) : 0.f);
    }
    return result;
}

// Same as get_V0_tracks_mass_TOF, but subtracts the per-vertex Ks flight
// time (get_V0_ks_t0) from the raw TOF before computing beta -- i.e. uses
// the muon's own flight time (decay vertex -> Wrapper/calo hit), not the
// time since the collision.
ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> get_V0_tracks_mass_TOF(
        ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> momentum,   // V0_tracks_p_raw
        ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> length_mm,  // V0_tracks_length_mm
        ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> tof_ps,     // V0_tracks_TOF_ps
        ROOT::VecOps::RVec<float> ks_t0_ps                        // get_V0_ks_t0
        ) {
    const float c_mm_per_ps = 0.299792458f;
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> result;
    for (size_t iv = 0; iv < momentum.size(); ++iv) {
        ROOT::VecOps::RVec<float> masses;
        for (size_t k = 0; k < momentum[iv].size(); ++k) {
            float p = momentum[iv][k];
            float L = length_mm[iv][k];
            float T = tof_ps[iv][k] - ks_t0_ps[iv];
            if (L <= 0 || T <= 0) { masses.push_back(-1.f); continue; }
            float beta = L / (c_mm_per_ps * T);
            if (beta >= 1.f) { masses.push_back(-1.f); continue; }
            masses.push_back(p * std::sqrt(1.f / (beta * beta) - 1.f));
        }
        result.push_back(masses);
    }
    return result;
}



///// I am adding a function to get the indices of the MC particles associated to particles with reconstructed tracks.
ROOT::VecOps::RVec<int> get_MCindTRK(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,
                                     ROOT::VecOps::RVec<edm4hep::TrackState> tracks,
                                     ROOT::VecOps::RVec<int> recind,
                                     ROOT::VecOps::RVec<int> mcind) {
    ROOT::VecOps::RVec<int> result;
    for (auto i : ReconstructedParticle2Track::get_recoindTRK(reco, tracks)) {
        int mc_index = ReconstructedParticle2MC::getTrack2MC_index(i, recind, mcind, reco);
        if (mc_index >= 0) result.push_back(mc_index);
    }
    return result;
}

struct sel_indices {
  sel_indices() {};

  ROOT::VecOps::RVec<edm4hep::MCParticleData>
  operator()(ROOT::VecOps::RVec<edm4hep::MCParticleData> in,
             ROOT::VecOps::RVec<int> indices) {

    ROOT::VecOps::RVec<edm4hep::MCParticleData> result;
    result.reserve(indices.size());

    for (auto idx : indices) {
      if (idx >= 0 && idx < (int)in.size()) {
        result.emplace_back(in[idx]);
      }
    }

    return result;
  }
};

ROOT::VecOps::RVec<int> get_MC_Vertex_isKSpipi(
    const std::vector<std::vector<int>>& PDGmother,
    const std::vector<std::vector<int>>& PDG)
{
    ROOT::VecOps::RVec<int> result;

    for (size_t i = 0; i < PDGmother.size(); ++i) {
        int isKS = 0;
        
        for (size_t j = 0; j < PDGmother[i].size(); ++j) {
            int mom_pdg = std::abs(PDGmother[i][j]);
            
            if (mom_pdg == 310) {
                isKS = 10; 
                int charged_daughters = 0;
                
                for (size_t k = 0; k < PDG[i].size(); ++k) {
                    int pdg_id = std::abs(PDG[i][k]);
                    
                    if (pdg_id == 211 || pdg_id == 13) {
                        charged_daughters++;
                    }
                }
                isKS += charged_daughters;
                break; 
            }
        }
        result.push_back(isKS);
    }
    return result;
}

ROOT::VecOps::RVec<int> get_KS_decays(const ROOT::VecOps::RVec<edm4hep::MCParticleData>& particles) {
    ROOT::VecOps::RVec<int> counts;
    for (const auto& p : particles) {
        if (std::abs(p.PDG) == 310) {
            // 
            int n_daughters = p.daughters_end - p.daughters_begin;
            if (n_daughters == 2) {
                counts.push_back(1);
            }
        }
    }
    return counts;
}

//A function to get the indices of those pions that are Ks' daughters
ROOT::VecOps::RVec<int> get_Ks_pions_indices(
    const ROOT::VecOps::RVec<edm4hep::MCParticleData>& particles,
    const ROOT::VecOps::RVec<int>& daughter_indices) { 
    
    ROOT::VecOps::RVec<int> indices;
    
    for (const auto& p : particles) {
        if (p.PDG == 310) { // Ks

            for (unsigned int i = p.daughters_begin; i < p.daughters_end; ++i) {
                
                int daughter_idx = daughter_indices.at(i);
                const auto& daughter = particles.at(daughter_idx);
                
                if (std::abs(daughter.PDG) == 211 || std::abs(daughter.PDG) == 13) {
                    indices.push_back(daughter_idx); 
                }
            }
        }
    }
    return indices;
}

//A function that takes the daughter of the Ks with less momentum
ROOT::VecOps::RVec<int> get_Ks_pion_min_p_idx(
    const ROOT::VecOps::RVec<edm4hep::MCParticleData>& particles,
    const ROOT::VecOps::RVec<int>& daughter_indices) { 
    
    ROOT::VecOps::RVec<int> min_indices;
    
    for (const auto& p : particles) {
        if (p.PDG == 310) {
            
            int best_idx = -1;
            float min_p = 1e30; 

            for (unsigned int i = p.daughters_begin; i < p.daughters_end; ++i) {
                if (i >= daughter_indices.size()) continue;

                int d_idx = daughter_indices.at(i);
                const auto& d = particles.at(d_idx);
                
                if (std::abs(d.PDG) == 211 || std::abs(d.PDG) == 13) {
                   
                    float p_mag = std::sqrt(d.momentum.x*d.momentum.x + 
                                            d.momentum.y*d.momentum.y + 
                                            d.momentum.z*d.momentum.z);
                    
                    if (p_mag < min_p) {
                        min_p = p_mag;
                        best_idx = d_idx;
                    }
                }
            }
            
            if (best_idx != -1) {
                min_indices.push_back(best_idx);
            }
        }
    }
    return min_indices;
}

ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>> get_Ks_daughter_status(
    const ROOT::VecOps::RVec<edm4hep::MCParticleData>& in, 
    const ROOT::VecOps::RVec<int>& ind) {
    
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>> res;

    for (size_t i = 0; i < in.size(); ++i) {
        // Buscamos el Ks (PDG 310)
        if (std::abs(in[i].PDG) != 310) continue;

        ROOT::VecOps::RVec<int> current_status;
        
        unsigned start = in[i].daughters_begin;
        unsigned end = in[i].daughters_end;

        for (unsigned j = start; j < end; ++j) {
            if (j < ind.size()) {
                int daughter_idx = ind.at(j);
                if (daughter_idx >= 0 && static_cast<size_t>(daughter_idx) < in.size()) {
                    current_status.push_back(in[daughter_idx].generatorStatus);
                }
            }
        }
        res.push_back(current_status);
    }
    return res;
}


ROOT::VecOps::RVec<int> get_Ks_daughter_pdg(
    const ROOT::VecOps::RVec<edm4hep::MCParticleData>& in, 
    const ROOT::VecOps::RVec<int>& ind) {
    
    ROOT::VecOps::RVec<int> result;
    for (size_t i = 0; i < in.size(); ++i) {
        if (std::abs(in[i].PDG) == 310) {
            for (unsigned j = in[i].daughters_begin; j < in[i].daughters_end; ++j) {
                if (j < ind.size()) {
                    int idx = ind.at(j);
                    if (idx >= 0 && (size_t)idx < in.size()) result.push_back(in[idx].PDG);
                }
            }
        }
    }
    return result;
}



ROOT::VecOps::RVec<int> get_Vertex_isMCKSpipi(
    ROOT::VecOps::RVec<int> Vertex_mcind,
    const std::vector<std::vector<int>>& PDGmother,  // These are the names of the inputs
    const std::vector<std::vector<int>>& PDG) 
{
    ROOT::VecOps::RVec<int> result;

    for (size_t i = 0; i < Vertex_mcind.size(); ++i) {
        int isKS = 0;
        int mc_index = Vertex_mcind[i];

        // 1. Check if vertex has the MC association
        if (mc_index >= 0 && mc_index < (int)PDGmother.size()) {
            
            bool hasKSMother = false;
            // 2. Check if the mother is a KS
            for (size_t j = 0; j < PDGmother[mc_index].size(); ++j) {
                if (std::abs(PDGmother[mc_index][j]) == 310) {
                    hasKSMother = true;
                    isKS = 10; 
                    break;
                }
            }

            // 3. Count the daughter pions if the KS mother is found
            if (hasKSMother) {
                for (size_t l = 0; l < PDG[mc_index].size(); ++l) {
                    if (std::abs(PDG[mc_index][l]) == 211) {
                        isKS += 1;
                    }
                }
            }
        } else {
            isKS = -1; 
        }
        
        result.push_back(isKS);
    }
    return result;
};

ROOT::VecOps::RVec<int> get_MC_Vertex_isPimunu(
    const std::vector<std::vector<int>>& PDGmother, 
    const std::vector<std::vector<int>>& PDG)
{
    ROOT::VecOps::RVec<int> result;
    for (size_t i = 0; i < PDGmother.size(); ++i) {
        bool hasPionMother = false;
        bool hasMuon = false;
        bool hasNu = false;

        for (size_t j = 0; j < PDGmother[i].size(); ++j) {
            if (std::abs(PDGmother[i][j]) == 211) {
                hasPionMother = true;
                break;
            }
        }

        if (hasPionMother) {
            for (size_t k = 0; k < PDG[i].size(); ++k) {
                int pdg_id = std::abs(PDG[i][k]);
                if (pdg_id == 13) hasMuon = true;
                if (pdg_id == 14) hasNu = true;
            }
        }

        int code = 0;
        if (hasPionMother) {
            code = 10 + (hasMuon ? 1 : 0) + (hasNu ? 1 : 0);
        }
        result.push_back(code);
    }
    return result;
};

//Function to get primary tracks.
ROOT::VecOps::RVec<edm4hep::TrackState> getPrimaryTracks(
        const ROOT::VecOps::RVec<edm4hep::TrackState>& tracks,
        double chi2max = 25.,
        double beamspotX = 0., double beamspotY = 0., double beamspotZ = 0.){

        // convert beamspot position units from centimeter to 10 micrometer
        // note: the FCCAnalyses function expect these values in micrometer,
        //       corresponding to track parameters in mm;
        //       but since our track parameters are in cm instead of mm,
        //       we use beamspot units of 10 micrometers.
        beamspotX = beamspotX * 1e3;
        beamspotY = beamspotY * 1e3;
        beamspotZ = beamspotZ * 1e3;

        // define beamspot width
        double sigma_beamspotX = 20; // unit: 10 micrometer
        double sigma_beamspotY = 10; // unit: 10 micrometer
        double sigma_beamspotZ = 2000; // unit: 10 micrometer
        bool doBeamSpotConstraint = true;

        // intitialize output
        ROOT::VecOps::RVec<edm4hep::TrackState> primaryTracks;

        // do filtering
        // note: the input is assumed to have already passed the baseline selection;
        //       here we just apply an extra cut on D0 and Z0 to focus on primary tracks
        ROOT::VecOps::RVec<edm4hep::TrackState> tracksToUse;
        for (const edm4hep::TrackState& trk : tracks) {
            if (std::abs(trk.D0)>50 || std::abs(trk.Z0)>50) continue;
            tracksToUse.push_back(trk);
        }
        if( tracksToUse.size() < 2 ){ return primaryTracks; }

        // call primary track finder from FCCAnalyses
        primaryTracks = FCCAnalyses::VertexFitterSimple::get_PrimaryTracks(
            tracksToUse,
            doBeamSpotConstraint,
            sigma_beamspotX, sigma_beamspotY, sigma_beamspotZ,
            beamspotX, beamspotY, beamspotZ
        );
        return primaryTracks;
};

// helper function to find tracks passing some baseline selection
// note: this function uses both the track states and the track objects (the latter for chi2);
//       it is assumed that both collections are corresponding trivially, i.e. one-to-one by same index.
ROOT::VecOps::RVec<edm4hep::TrackState> getSelectedTracks(
        const ROOT::VecOps::RVec<edm4hep::TrackState>& tracks,
        const ROOT::VecOps::RVec<edm4hep::TrackData>& trackDatas,
        const float D0max, const float Z0max){

        // safety check
        if( tracks.size() != trackDatas.size() ){
            throw std::runtime_error("Vector sizes do not match (in getSelectedTracks plain).");
        }

        // do filtering
        ROOT::VecOps::RVec<edm4hep::TrackState> selectedTracks;
        for (int idx = 0; idx < tracks.size(); idx++) {
            const edm4hep::TrackState trk = tracks.at(idx);
            const edm4hep::TrackData trkobj = trackDatas.at(idx);
            const auto& c = trk.covMatrix;
            if (c[0] <= 0 || c[2] <= 0 || c[9] <= 0) continue;
            if (c[0] < 1e-12 || c[2] < 1e-12 || c[9] <= 1e-12) continue;
            if (!std::isfinite(c[0]) || !std::isfinite(c[2]) || !std::isfinite(c[9])) continue;
            if (D0max > 0 && std::abs(trk.D0)>D0max) continue;
            if (Z0max > 0 && std::abs(trk.Z0)>Z0max) continue;
            if (trkobj.ndf==0) continue;
            if (trkobj.chi2 / trkobj.ndf > 10.) continue;
            selectedTracks.push_back(trk);
        }
        return selectedTracks;
};

FCCAnalyses::VertexingUtils::FCCAnalysesVertex fitRecoPrimaryVertex(
        const ROOT::VecOps::RVec<edm4hep::TrackState>& tracks,
        double bsX = 0., double bsY = 0., double bsZ = 0.,            // beam-spot position, um
        double sigX = 5.96, double sigY = 0.0238, double sigZ = 397.){ // beam-spot size, um (Z pole)

        edm4hep::VertexData dummyVertex;
        dummyVertex.chi2 = -1;
        dummyVertex.position = edm4hep::Vector3f(bsX*1e-3, bsY*1e-3, bsZ*1e-3); // um -> mm
        FCCAnalyses::VertexingUtils::FCCAnalysesVertex dummyVertexObject;
        dummyVertexObject.vertex = dummyVertex;
        dummyVertexObject.ntracks = 0;
        dummyVertexObject.mc_ind = -1;
        if( tracks.size() < 2 ){ return dummyVertexObject; }

        return FCCAnalyses::VertexFitterSimple::VertexFitter_Tk(
            1, tracks, tracks, true, sigX, sigY, sigZ, bsX, bsY, bsZ);
};

// helper function to find tracks not compatible with primary vertex.
ROOT::VecOps::RVec<edm4hep::TrackState> getSecondaryTracks(
        const ROOT::VecOps::RVec<edm4hep::TrackState>& tracks,
        const ROOT::VecOps::RVec<edm4hep::TrackState>& primaryTracks){

        // skip tracks compatible with primary vertex
        ROOT::VecOps::RVec<edm4hep::TrackState> secondaryTracks;
        secondaryTracks = FCCAnalyses::VertexFitterSimple::get_NonPrimaryTracks(tracks, primaryTracks);
        return secondaryTracks;
};

struct MatchMCV0 {
    std::vector<float> distance;
    std::vector<int> mc_vtx_idx;

    std::vector<std::vector<int>> mc_daughter_pdg;
    std::vector<std::vector<int>> mc_mother_pdg;
    std::vector<std::vector<int>> mc_gmother_pdgs;

    std::vector<std::vector<float>> mc_vtx_x;
    std::vector<std::vector<float>> mc_vtx_y;
    std::vector<std::vector<float>> mc_vtx_z;

    ROOT::VecOps::RVec<int> is_valid; // 1 within tolerance, 0 for particles that formed V0 but have different origin 

    std::vector<std::vector<int>> mc_daughter_indices;
    std::vector<std::vector<int>> mc_mother_indices;
    std::vector<std::vector<int>> mc_gmother_indices;
};

MatchMCV0 match_MC2Reco_V0s(
    const ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex>& recoVertices,
    const ROOT::VecOps::RVec<int>& rp2mc, // Association reco
    const ROOT::VecOps::RVec<int>& mc2rp, // Association MC
    const ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertexMC>& mcvertices,
    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco_particles,
    const ROOT::VecOps::RVec<edm4hep::MCParticleData>& mc_particles,
    const ROOT::VecOps::RVec<int>& parents_ind,
    const ROOT::VecOps::RVec<float>& MC_Vertex_x, 
    const ROOT::VecOps::RVec<float>& MC_Vertex_y,
    const ROOT::VecOps::RVec<float>& MC_Vertex_z
) {
    MatchMCV0 result;

    // Iterate over each reco vertex
    for (const auto& vtx : recoVertices) {
        int matched_vtx_idx = -1;
        float dist = -999.0;

        bool valid_mc_origin = true;
        
        // this will store the reco particles that form such reco vertex
        std::vector<int> asso_vtx_mc_particles;
        std::vector<int> vtx_pdgs;
        std::vector<int> vtx_mother_pdgs;
        std::vector<int> vtx_gmother_pdgs;

        std::vector<float> vtx_xs;
        std::vector<float> vtx_ys;
        std::vector<float> vtx_zs;

        std::vector<int> mc_daughter_indices;
        std::vector<int> mc_mother_indices;
        std::vector<int> mc_gmother_indices;
        
        for (int rp_idx : vtx.reco_ind) {
            bool found_mc = false; 

            for (size_t j = 0; j < rp2mc.size(); ++j) {
                if (rp2mc.at(j) == rp_idx) {
                    int mc_p_id = mc2rp.at(j);
                    if (mc_p_id >= 0 && mc_p_id < (int)mc_particles.size()) {
                        found_mc = true;
                        asso_vtx_mc_particles.push_back(mc_p_id);
                        vtx_pdgs.push_back(mc_particles.at(mc_p_id).PDG);
                        vtx_xs.push_back(mc_particles.at(mc_p_id).vertex.x);
                        vtx_ys.push_back(mc_particles.at(mc_p_id).vertex.y);
                        vtx_zs.push_back(mc_particles.at(mc_p_id).vertex.z);

                        mc_daughter_indices.push_back(mc_p_id);
                        //store index of daughter MC particle: I will need later to access the mc particle collection
                        //and get momentum etc

                        int mother_pdg = -999;
                        int gmother_pdg = -999;
                        unsigned int beg = mc_particles.at(mc_p_id).parents_begin;
                        unsigned int end = mc_particles.at(mc_p_id).parents_end;
                        int real_gmother_idx = -1;

                        if (beg < end && beg < parents_ind.size()) {
                            
                            int real_mother_idx = parents_ind.at(beg);
                            
                            if (real_mother_idx >= 0 && real_mother_idx < (int)mc_particles.size()) {
                                mother_pdg = mc_particles.at(real_mother_idx).PDG;

                                mc_mother_indices.push_back(real_mother_idx);
                                //I need to also store the index of the mother particle to be able to check later if the mother is the same for all daughters.
                                
                                unsigned int g_beg = mc_particles.at(real_mother_idx).parents_begin;
                                unsigned int g_end = mc_particles.at(real_mother_idx).parents_end;

                                if (g_beg < g_end && g_beg < parents_ind.size()) {
                                    real_gmother_idx = parents_ind.at(g_beg);

                                    if (real_gmother_idx >= 0 && real_gmother_idx < (int)mc_particles.size()) {
                                        gmother_pdg = mc_particles.at(real_gmother_idx).PDG;
                                        vtx_gmother_pdgs.push_back(gmother_pdg);
                                        mc_gmother_indices.push_back(real_gmother_idx);
                                    } else {
                                        mc_gmother_indices.push_back(-1);
                                        vtx_gmother_pdgs.push_back(-999);
                                    }
                                } else {
                                    mc_gmother_indices.push_back(-1);
                                    vtx_gmother_pdgs.push_back(0);
                                }
                            } else { 
                                mc_mother_indices.push_back(-1);
                                mc_gmother_indices.push_back(-1);
                                vtx_gmother_pdgs.push_back(-999);
                            }
                        } else { 
                            mc_mother_indices.push_back(-1);
                            mc_gmother_indices.push_back(-1);
                            vtx_gmother_pdgs.push_back(-999);
                        }
                        vtx_mother_pdgs.push_back(mother_pdg);
                    }
                    break;
                }
            }

            if (!found_mc) {
                vtx_pdgs.push_back(0);
                vtx_xs.push_back(-999.0);
                vtx_ys.push_back(-999.0);
                vtx_zs.push_back(-999.0);
                mc_daughter_indices.push_back(-1);
                vtx_mother_pdgs.push_back(0);
                mc_mother_indices.push_back(-1);
                vtx_gmother_pdgs.push_back(0);
                mc_gmother_indices.push_back(-1);
            }
        } 

        // 2. check if daughters MC come from the same position
        //if (asso_vtx_mc_particles.size() >= 2) {
            
         //   float tolerance = 1.0; // Tolerance of 1 mm

            // Vertex coordinates of first daughter MC particle
         //   const auto& vtx0 = mc_particles.at(asso_vtx_mc_particles[0]).vertex;

        //    for (size_t i = 1; i < asso_vtx_mc_particles.size(); ++i) {
        //       const auto& vtx_i = mc_particles.at(asso_vtx_mc_particles[i]).vertex;
                
                // Calculate distance: sqrt(dx^2 + dy^2 + dz^2)
        //        float dx = vtx_i.x - vtx0.x;
        //      float dy = vtx_i.y - vtx0.y;
        //        float dz = vtx_i.z - vtx0.z;
                
        //       float distance = std::sqrt(dx*dx + dy*dy + dz*dz);

          //    if (distance > tolerance) {
            //    valid_mc_origin = false;
              //    break;
                //}
            //}
        //}

            // 3. Look for the corresponding MC vertex
         //   if (valid_mc_origin) {
                std::vector<int> matching_mc_vertices;
                
               // for (size_t i = 0; i < mcvertices.size(); ++i) {
                    //const auto& mc_vtx_indices = mcvertices.at(i).mc_ind;
                    
                    // Check that the mc particles are in the mc vertex
                    //Note that if I do this, the previous check does not make sense anymore 
                    //because mc daughters of an mc vertex have to be at least 10^-5 mm close.
                    //for (int mc_p_id : asso_vtx_mc_particles) {
                    //    if (std::find(mc_vtx_indices.begin(), mc_vtx_indices.end(), mc_p_id) != mc_vtx_indices.end()) {
                    //        matching_mc_vertices.push_back(i);
                    //        break;
                    //    }
                  //  }
                //}

                // calculate distance reco - mc vertex
            //    if (!matching_mc_vertices.empty()) {
              //      matched_vtx_idx = matching_mc_vertices[0];
                    
           //          float dx = vtx.vertex.position.x - MC_Vertex_x[matched_vtx_idx];
            //        float dy = vtx.vertex.position.y - MC_Vertex_y[matched_vtx_idx];
            //        float dz = vtx.vertex.position.z - MC_Vertex_z[matched_vtx_idx];
            //        dist = std::sqrt(dx*dx + dy*dy + dz*dz);
            //    }
        //}

        // CRITICAL: Validate that all daughters come from the same MC mother
        // This prevents creating fake V0s from muons originating from different Ks particles
        bool same_mc_mother = true;
        if (mc_daughter_indices.size() >= 2) {
            // Check if all mother indices are valid and identical
            for (size_t i = 0; i < mc_mother_indices.size(); ++i) {
                if (mc_mother_indices[i] < 0) {
                    same_mc_mother = false;
                    break;
                }
                if (i > 0 && mc_mother_indices[i] != mc_mother_indices[0]) {
                    same_mc_mother = false;
                    break;
                }
            }
        } else {
            same_mc_mother = false;
        }

        result.mc_vtx_idx.push_back(matched_vtx_idx);
        result.distance.push_back(dist);
        result.mc_daughter_pdg.push_back(vtx_pdgs);
        result.mc_mother_pdg.push_back(vtx_mother_pdgs);
        result.mc_gmother_pdgs.push_back(vtx_gmother_pdgs);

        result.mc_vtx_x.push_back(vtx_xs);
        result.mc_vtx_y.push_back(vtx_ys);
        result.mc_vtx_z.push_back(vtx_zs);

        result.mc_daughter_indices.push_back(mc_daughter_indices);
        result.mc_mother_indices.push_back(mc_mother_indices);
        result.mc_gmother_indices.push_back(mc_gmother_indices);
    //    result.is_valid.push_back(same_mc_mother ? 1 : 0);
    }

    return result;
}



std::vector<float> get_v0_dist(MatchMCV0 res) { return res.distance; }
std::vector<int> get_v0_idx(MatchMCV0 res) { return res.mc_vtx_idx; }

// function to calculate minimum distance between a Kaon trajectory (vertex + momentum) and another vertex
// to check compability between the KS decay vertex and other vertices. 
float get_distance_traj2vertex(edm4hep::Vector3f vPos, TVector3 pVec, edm4hep::Vector3f v2Pos) {
    
    // unitary vector alomg the Kaon momentum direction
    float pMag = sqrt(pVec.X()*pVec.X() + pVec.Y()*pVec.Y() + pVec.Z()*pVec.Z());
    if (pMag < 1e-6) return 999.0; // avoid errors if momentum is too small
    
    TVector3 u(pVec.X()/pMag, pVec.Y()/pMag, pVec.Z()/pMag);

    // vector from the kaon vertex to the other vertex
    TVector3 d(v2Pos.x - vPos.x, v2Pos.y - vPos.y, v2Pos.z - vPos.z);

    // one can measure the distance of the kaon trajectory to the other vertex as the magnitude of the
    // cross product between unitary vector and vector d.
    TVector3 crossProd = d.Cross(u);
    
    return crossProd.Mag();
}

//function to get the minimum distance between all KS candidates and other vertices in the event, 
// to check if there are other vertices compatible with the KS decay vertex.

struct KSVertexCompatibility {
    ROOT::VecOps::RVec<float> distances;
    ROOT::VecOps::RVec<int> indices;
    ROOT::VecOps::RVec<int> is_correct_origin;
};

KSVertexCompatibility get_ks2vertex_min_dist(
    ROOT::VecOps::RVec<int> ks_mask,
    ROOT::VecOps::RVec<float> px,
    ROOT::VecOps::RVec<float> py, 
    ROOT::VecOps::RVec<float> pz,
    ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vtx,
    ROOT::VecOps::RVec<int> mc_vtx_ind, // Vertex_mcind
    ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertexMC> mc_vtx
) {
    ROOT::VecOps::RVec<float> min_distances;
    ROOT::VecOps::RVec<int> best_vtx_idx;
    ROOT::VecOps::RVec<int> is_correct_origin(vtx.size(), 0);

    // first check if vertex is a KS candidate
    for (size_t i = 0; i < vtx.size(); ++i) {
        if (ks_mask[i] == 0) continue; 

        // get the position and momentum of the KS candidate vertex
        edm4hep::Vector3f vPos = vtx.at(i).vertex.position;
        TVector3 pVec(px[i], py[i], pz[i]);

        float min_d = 1e3;
        int min_vtx_idx = -1;

        // compare the ks candidate vertex with all other vertices in the event
        for (size_t j = 0; j < vtx.size(); ++j) {
            if (i == j) continue; // you don't compare the vertex with itself
            
            float d = get_distance_traj2vertex(vPos, pVec, vtx.at(j).vertex.position);
            if (d < min_d) {
                min_d = d;
                min_vtx_idx = j;
            }
        }
        
        // store the minimum distance and the index of the closest vertex
        if (min_vtx_idx != -1) {
            min_distances.push_back(min_d);
            best_vtx_idx.push_back(min_vtx_idx);
            is_correct_origin.push_back(0);
            
            int mc_idx_ks = mc_vtx_ind[i];
            int mc_idx_parent = mc_vtx_ind[min_vtx_idx]; // MC vertex of closest reco vertex

            if (mc_idx_ks >= 0 && mc_idx_parent >= 0) {
                // pdg of kaon's mother
                auto mothers_ks = mc_vtx.at(mc_idx_ks).gmother_ind;
                auto mothers_found = mc_vtx.at(mc_idx_parent).mother_ind; 
                
                for (const auto& mks : mothers_ks) {
                    for (const auto& mf : mothers_found) {
                        if (mks == mf) {
                            is_correct_origin[i] = 1;
                            break;
                        }
                    }
                    if (is_correct_origin[i] == 1) break;
                }
            }
        }
    }

    return {min_distances, best_vtx_idx, is_correct_origin};
}


//Defining a function to get the distance between close particles. I do this to check which pions and muons
//are close to each other at the vertex level 
struct CloseParticles {
        ROOT::VecOps::RVec<float> d;
        ROOT::VecOps::RVec<float> pix, piy, piz;
        ROOT::VecOps::RVec<float> mux, muy, muz;
    };

// Function to get the ID of the grid cell for a given 3D point, to organise particles in a grid for fast searching of close pairs.
    inline long long get_grid_ID(float x, float y, float z, float step) {
        long long ix = std::floor(x / step);
        long long iy = std::floor(y / step);
        long long iz = std::floor(z / step);
        return (ix + 10000) * 1000000 + (iy + 10000) * 1000 + (iz + 10000);
    }

inline CloseParticles get_close_pairs(
        ROOT::VecOps::RVec<float> pix, ROOT::VecOps::RVec<float> piy, ROOT::VecOps::RVec<float> piz,
        ROOT::VecOps::RVec<float> mux, ROOT::VecOps::RVec<float> muy, ROOT::VecOps::RVec<float> muz,
        float threshold) 
    {
        CloseParticles res;
        // 1. organise muons in a grid for fast searching
        std::unordered_map<long long, std::vector<size_t>> grid;
        for (size_t j = 0; j < mux.size(); ++j) {
            grid[get_grid_ID(mux[j], muy[j], muz[j], threshold)].push_back(j);
        }

        // 2. iterate over pions and search for close muons in the grid
        for (size_t i = 0; i < pix.size(); ++i) {
            long long ix = std::floor(pix[i] / threshold);
            long long iy = std::floor(piy[i] / threshold);
            long long iz = std::floor(piz[i] / threshold);

            for (int dx = -1; dx <= 1; ++dx) {
                for (int dy = -1; dy <= 1; ++dy) {
                    for (int dz = -1; dz <= 1; ++dz) {
                        long long key = (ix + dx + 10000) * 1000000 + (iy + dy + 10000) * 1000 + (iz + dz + 10000);
                        
                        if (grid.count(key)) {
                            for (size_t muon_idx : grid[key]) {
                                float dist = std::sqrt(std::pow(pix[i]-mux[muon_idx],2) + 
                                                       std::pow(piy[i]-muy[muon_idx],2) + 
                                                       std::pow(piz[i]-muz[muon_idx],2));
                                if (dist < threshold) {
                                    res.d.push_back(dist);
                                    res.pix.push_back(pix[i]); res.piy.push_back(piy[i]); res.piz.push_back(piz[i]);
                                    res.mux.push_back(mux[muon_idx]); res.muy.push_back(muy[muon_idx]); res.muz.push_back(muz[muon_idx]);
                                }
                            }
                        }
                    }
                }
            }
        }
        return res;
    }






//add all functions to reconstruct V0s

VertexingUtils::FCCAnalysesV0 get_V0s(ROOT::VecOps::RVec<edm4hep::TrackState> np_tracks,
                                      const ROOT::VecOps::RVec<edm4hep::TrackState>& alltracks,
                                      ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recoparticles,
                                      VertexingUtils::FCCAnalysesVertex PV,
                                      bool tight,
                                      double chi2_cut) {
  // V0 reconstruction
  // if(tight)  -> tight constraints
  // if(!tight) -> loose constraints

  VertexingUtils::FCCAnalysesV0 result;
  ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vtx; // FCCAnalyses vertex object
  ROOT::VecOps::RVec<int> pdgAbs;                            // absolute PDG ID
  ROOT::VecOps::RVec<double> invM;                           // invariant mass
  result.vtx = vtx;
  result.pdgAbs = pdgAbs;
  result.invM = invM;


  VertexingUtils::FCCAnalysesVertex V0_vtx;
  
  int nTr = np_tracks.size();
  if(nTr<2) return result;
  ROOT::VecOps::RVec<bool> isInV0(nTr, false);

  // set constraints (if(tight==true) tight_set)
  ROOT::VecOps::RVec<double> isKs      = constraints_Ks(tight);
  ROOT::VecOps::RVec<double> isLambda0 = constraints_Lambda0(tight);
  ROOT::VecOps::RVec<double> isGamma   = constraints_Gamma(tight);

  ROOT::VecOps::RVec<edm4hep::TrackState> tr_pair;
  // push empty tracks to make a size=2 vector
  edm4hep::TrackState tr_i, tr_j;
  tr_pair.push_back(tr_i);
  tr_pair.push_back(tr_j);
  //
  for(unsigned int i=0; i<nTr-1; i++) {
    if(isInV0[i] == true) continue; // don't pair a track if it already forms a V0
    tr_pair[0] = np_tracks[i];

    for(unsigned int j=i+1; j<nTr; j++) {
      if(isInV0[j] == true) continue; // don't pair a track if it already forms a V0
      if(tr_pair[0].omega * np_tracks[j].omega > 0) continue; // don't pair tracks with same charge (same sign curvature = same sign charge)
      tr_pair[1] = np_tracks[j];

      ROOT::VecOps::RVec<double> V0_cand = get_V0candidate(V0_vtx, tr_pair, PV, true, chi2_cut);
      if(V0_cand[0] == -1) continue;
      
      // Ks
      if(V0_cand[0]>isKs[0] && V0_cand[0]<isKs[1] && V0_cand[5]>isKs[2] && V0_cand[6]>isKs[3]) {
if(debug_me) std::cout<<"Found a Ks"<<std::endl;
        isInV0[i] = true;
        isInV0[j] = true;

        V0_vtx.reco_ind.clear(); 
        for (int k = 0; k < recoparticles.size(); k++) {
            auto &p = recoparticles[k];
            if (p.tracks_begin >= 0 && p.tracks_begin < alltracks.size()) {
                edm4hep::TrackState track_from_reco = alltracks[p.tracks_begin];
                
                for (auto &track_in_vertex : tr_pair) {
                    if (VertexingUtils::compare_Tracks(track_from_reco, track_in_vertex)) {
                        V0_vtx.reco_ind.push_back(k);
                        break;
                    }
                }
            }
        }

        vtx.push_back(V0_vtx);
        pdgAbs.push_back(310);
        invM.push_back(V0_cand[0]);
        break;
      }
      
      // Lambda0
      else if(V0_cand[1]>isLambda0[0] && V0_cand[1]<isLambda0[1] && V0_cand[5]>isLambda0[2] && V0_cand[6]>isLambda0[3]) {
if(debug_me) std::cout<<"Found a Lambda0"<<std::endl;
        isInV0[i] = true;
        isInV0[j] = true;

        V0_vtx.reco_ind.clear();
        for (int k = 0; k < recoparticles.size(); k++) {
            auto &p = recoparticles[k];
            if (p.tracks_begin >= 0 && p.tracks_begin < alltracks.size()) {
                edm4hep::TrackState track_from_reco = alltracks[p.tracks_begin];
                for (auto &track_in_vertex : tr_pair) {
                    if (VertexingUtils::compare_Tracks(track_from_reco, track_in_vertex)) {
                        V0_vtx.reco_ind.push_back(k);
                        break;
                    }
                }
            }
        }

        vtx.push_back(V0_vtx);
        pdgAbs.push_back(3122);
        invM.push_back(V0_cand[2]);
        break;
      }

      else if(V0_cand[2]>isLambda0[0] && V0_cand[2]<isLambda0[1] && V0_cand[5]>isLambda0[2] && V0_cand[6]>isLambda0[3]) {
if(debug_me) std::cout<<"Found a Lambda0"<<std::endl;
        isInV0[i] = true;
        isInV0[j] = true;

        V0_vtx.reco_ind.clear(); 
        for (int k = 0; k < recoparticles.size(); k++) {
            auto &p = recoparticles[k];
            if (p.tracks_begin >= 0 && p.tracks_begin < alltracks.size()) {
                edm4hep::TrackState track_from_reco = alltracks[p.tracks_begin];
                for (auto &track_in_vertex : tr_pair) {
                    if (VertexingUtils::compare_Tracks(track_from_reco, track_in_vertex)) {
                        V0_vtx.reco_ind.push_back(k);
                        break;
                    }
                }
            }
        }

        vtx.push_back(V0_vtx);
        pdgAbs.push_back(3122);
        invM.push_back(V0_cand[1]);
        break;
      }
	
      // photon conversion
      else if(V0_cand[3]<isGamma[1] && V0_cand[5]>isGamma[2] && V0_cand[6]>isGamma[3]) {
if(debug_me) std::cout<<"Found a Photon coversion"<<std::endl;
        isInV0[i] = true;
        isInV0[j] = true;

        V0_vtx.reco_ind.clear();
        for (int k = 0; k < recoparticles.size(); k++) {
            auto &p = recoparticles[k];
            if (p.tracks_begin >= 0 && p.tracks_begin < alltracks.size()) {
                edm4hep::TrackState track_from_reco = alltracks[p.tracks_begin];
                for (auto &track_in_vertex : tr_pair) {
                    if (VertexingUtils::compare_Tracks(track_from_reco, track_in_vertex)) {
                        V0_vtx.reco_ind.push_back(k);
                        break;
                    }
                }
            }
        }

        vtx.push_back(V0_vtx);
        pdgAbs.push_back(22);
        invM.push_back(V0_cand[3]);
        break;
      }
      //dimuon 

      else if(V0_cand[4]>isKs[0] && V0_cand[4]<isKs[1] && V0_cand[5]>isKs[2] && V0_cand[6]>isKs[3]) {
        if (debug_me) std::cout<<"Found a dimuon"<<std::endl;
        isInV0[i] = true;
        isInV0[j] = true;
        vtx.push_back(V0_vtx);
        pdgAbs.push_back(1313);
        invM.push_back(V0_cand[4]);
        break;
            }

    }
  }

  result.vtx = vtx;
  result.pdgAbs = pdgAbs;
  result.invM = invM;
  //
  return result;
}



ROOT::VecOps::RVec<double> get_V0candidate(VertexingUtils::FCCAnalysesVertex &V0_vtx,
					   ROOT::VecOps::RVec<edm4hep::TrackState> tr_pair,
					   VertexingUtils::FCCAnalysesVertex PV,
					   bool chi2,
					   double chi2_cut)
{
  // get invariant mass, distance from PV, and colliniarity variables for all V0 candidates
  
  // [0] -> invM_Ks [GeV]
  // [1] -> invM_Lambda1 [GeV]
  // [2] -> invM_Lambda2 [GeV]
  // [3] -> invM_Gamma [GeV]
  // [4] -> r (distance from PV) [mm]
  // [5] -> r.p (colinearity) [r & p - unit vectors]
  // skip the candidate with output entries = -1
  
  ROOT::VecOps::RVec<double> result(7, -1);

  edm4hep::Vector3f r_PV = PV.vertex.position; // in mm
  
  V0_vtx = VertexFitterSimple::VertexFitter_Tk(2, tr_pair);

  if(chi2) {
    // constraint on chi2: chi2 < cut (9)
    double chi2_V0 = V0_vtx.vertex.chi2; // normalised but nDOF=1
    if(chi2_V0 >= chi2_cut) return result;
  }
  
  // invariant masses for V0 candidates
  result[0] = VertexingUtils::get_invM_pairs(V0_vtx, m_pi, m_pi);
  result[1] = VertexingUtils::get_invM_pairs(V0_vtx, m_pi, m_p);
  result[2] = VertexingUtils::get_invM_pairs(V0_vtx, m_p, m_pi);
  result[3] = VertexingUtils::get_invM_pairs(V0_vtx, m_e, m_e);
  result[4] = VertexingUtils::get_invM_pairs(V0_vtx, m_mu, m_mu);

  // V0 candidate distance from PV
  edm4hep::Vector3f r_V0 = V0_vtx.vertex.position; // in mm
  TVector3 r_V0_PV(r_V0[0] - r_PV[0], r_V0[1] - r_PV[1], r_V0[2] - r_PV[2]);
  result[5] = r_V0_PV.Mag(); // in mm

  // angle b/n V0 candidate momentum & PV-V0 displacement vector
  result[6] = VertexingUtils::get_PV2V0angle(V0_vtx, PV);

  return result;
}

ROOT::VecOps::RVec<double> constraints_Ks(bool tight) {

  ROOT::VecOps::RVec<double> result(4, 0);

  if(tight) {
    result[0] = 0.493;
    result[1] = 0.503;
    result[2] = 0.5;
    result[3] = 0.999;
  }

  else {
    result[0] = 0.0;
    result[1] = 10;
    result[2] = 0.0;
    result[3] = 0.999;
  }
  //
  return result;
}

ROOT::VecOps::RVec<double> constraints_Lambda0(bool tight) {

  ROOT::VecOps::RVec<double> result(4, 0);

  if(tight) {
    result[0] = 1.111;
    result[1] = 1.121;
    result[2] = 0.5;
    result[3] = 0.99995;
  }

  else {
    result[0] = 1.106;
    result[1] = 1.126;
    result[2] = 0.3;
    result[3] = 0.999;
  }
  //
  return result;
}

ROOT::VecOps::RVec<double> constraints_Gamma(bool tight) {

  ROOT::VecOps::RVec<double> result(4, 0);

  if(tight) {
    result[1] = 0.005;
    result[2] = 9;
    result[3] = 0.99995;
  }

  else {
    result[1] = 0.01;
    result[2] = 9;
    result[3] = 0.999;
  }
  //
  return result;
}

// I will add a function to look for muons that come from pions and check if those pions come from Kaons.
// Loop over all MC particles and identify the decay chain:
//     K0S -> pi± -> mu±
//
// For every muon found in the event:
//   1. Check whether its direct mother is a charged pion (PDG = ±211)
//   2. If true, check whether that pion comes from a KS (PDG = 310)
//   3. Store the indices of the muon, pion, and KS
//
// The function returns a struct containing:
//   - muon indices
//   - pion mother indices
//   - KS grandmother indices
struct PiMuFromKs {
    ROOT::VecOps::RVec<int> mu_idx;   // muon MC index (from pi->mu nu)
    ROOT::VecOps::RVec<int> pi_idx;   // decaying pion MC index (no reco track in Delphes)
    ROOT::VecOps::RVec<int> ks_idx;   // Ks MC index
    ROOT::VecOps::RVec<int> pi2_idx;  // stable sibling pion MC index (the other Ks daughter)
};

// Particle0 = Particle#0.index (parent links), Particle1 = Particle#1.index (daughter links)
PiMuFromKs get_PiMuFromKs(
    ROOT::VecOps::RVec<edm4hep::MCParticleData> particles,
    ROOT::VecOps::RVec<int> Particle0,
    ROOT::VecOps::RVec<int> Particle1
) {
    PiMuFromKs result;
    for (size_t i = 0; i < particles.size(); ++i) {
        const auto& mu = particles.at(i);
        if (std::abs(mu.PDG) != 13) continue;
        if (mu.parents_begin >= mu.parents_end) continue;
        for (size_t j = mu.parents_begin; j < mu.parents_end; ++j) {
            if (j >= Particle0.size()) continue;
            int pi_idx = Particle0.at(j);
            if (pi_idx < 0 || pi_idx >= (int)particles.size()) continue;
            const auto& pion = particles.at(pi_idx);
            if (std::abs(pion.PDG) != 211) continue;
            if (pion.parents_begin >= pion.parents_end) continue;
            for (size_t k = pion.parents_begin; k < pion.parents_end; ++k) {
                if (k >= Particle0.size()) continue;
                int ks_idx = Particle0.at(k);
                if (ks_idx < 0 || ks_idx >= (int)particles.size()) continue;
                const auto& ks = particles.at(ks_idx);
                if (ks.PDG != 310) continue;

                // Find the stable sibling pion: the other charged daughter of this Ks
                int pi2_idx = -1;
                for (int d = ks.daughters_begin; d < ks.daughters_end; ++d) {
                    if (d < 0 || d >= (int)Particle1.size()) continue;
                    int dau = Particle1.at(d);
                    if (dau < 0 || dau >= (int)particles.size()) continue;
                    if (dau == pi_idx) continue;
                    if (std::abs(particles.at(dau).PDG) == 211) { pi2_idx = dau; break; }
                }

                result.mu_idx.push_back(i);
                result.pi_idx.push_back(pi_idx);
                result.ks_idx.push_back(ks_idx);
                result.pi2_idx.push_back(pi2_idx);
            }
        }
    }
    return result;
}

// For each Ks->pi->mu pair identified by get_PiMuFromKs, return the reconstructed
// track states of the pion and muon by walking the MC-Reco association table.
// For each Ks->pi->mu pair, retrieve the reco track states of:
//   - the muon (from pi->mu nu, has a reco track in Delphes)
//   - the stable pion pi2 (the other Ks daughter that did NOT decay in flight)
// The decaying pion itself has no reco track in Delphes parametric simulation.
// Results are stored as consecutive pairs: [mu_ts, pi2_ts, mu_ts, pi2_ts, ...]
ROOT::VecOps::RVec<edm4hep::TrackState> get_pimu_reco_tracks(
    const ROOT::VecOps::RVec<int>& mu_mc_idx,
    const ROOT::VecOps::RVec<int>& pi2_mc_idx,  // stable sibling pion from same Ks
    const ROOT::VecOps::RVec<int>& assoc_reco,   // MCRecoAssociations#0.index
    const ROOT::VecOps::RVec<int>& assoc_mc,     // MCRecoAssociations#1.index
    const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& RecoParticles,
    const ROOT::VecOps::RVec<edm4hep::TrackState>& AllTracks)
{
    ROOT::VecOps::RVec<edm4hep::TrackState> result;

    // Build MC index → reco index map
    std::unordered_map<int,int> mc2reco;
    mc2reco.reserve(assoc_mc.size());
    for (size_t j = 0; j < assoc_mc.size(); ++j)
        mc2reco[assoc_mc[j]] = assoc_reco[j];

    size_t n = std::min(mu_mc_idx.size(), pi2_mc_idx.size());
    for (size_t i = 0; i < n; ++i) {
        if (pi2_mc_idx[i] < 0) continue;  // no sibling pion found
        auto it_mu  = mc2reco.find(mu_mc_idx[i]);
        auto it_pi2 = mc2reco.find(pi2_mc_idx[i]);
        if (it_mu == mc2reco.end() || it_pi2 == mc2reco.end()) continue;

        int mu_reco  = it_mu->second;
        int pi2_reco = it_pi2->second;
        if (mu_reco  < 0 || mu_reco  >= (int)RecoParticles.size()) continue;
        if (pi2_reco < 0 || pi2_reco >= (int)RecoParticles.size()) continue;

        ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> mu_rp  = {RecoParticles[mu_reco]};
        ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> pi2_rp = {RecoParticles[pi2_reco]};
        auto mu_ts  = FCCAnalyses::ReconstructedParticle2Track::getRP2TRK(mu_rp,  AllTracks);
        auto pi2_ts = FCCAnalyses::ReconstructedParticle2Track::getRP2TRK(pi2_rp, AllTracks);

        if (mu_ts.empty() || pi2_ts.empty()) continue;
        result.push_back(mu_ts[0]);
        result.push_back(pi2_ts[0]);
    }
    return result;
}


// Returns the 0-based indices (into the MC pair list) of pairs where both tracks were found.
// Use with ROOT::VecOps::Take(pimu_mc_ksdec_r, pimu_found_idx) to get the "found" sub-distribution.
ROOT::VecOps::RVec<int> get_pimu_found_idx(
    const ROOT::VecOps::RVec<int>& mu_mc_idx,
    const ROOT::VecOps::RVec<int>& pi2_mc_idx,
    const ROOT::VecOps::RVec<int>& assoc_reco,
    const ROOT::VecOps::RVec<int>& assoc_mc,
    const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& RecoParticles,
    const ROOT::VecOps::RVec<edm4hep::TrackState>& AllTracks)
{
    ROOT::VecOps::RVec<int> result;

    std::unordered_map<int,int> mc2reco;
    mc2reco.reserve(assoc_mc.size());
    for (size_t j = 0; j < assoc_mc.size(); ++j)
        mc2reco[assoc_mc[j]] = assoc_reco[j];

    size_t n = std::min(mu_mc_idx.size(), pi2_mc_idx.size());
    for (size_t i = 0; i < n; ++i) {
        if (pi2_mc_idx[i] < 0) continue;
        auto it_mu  = mc2reco.find(mu_mc_idx[i]);
        auto it_pi2 = mc2reco.find(pi2_mc_idx[i]);
        if (it_mu == mc2reco.end() || it_pi2 == mc2reco.end()) continue;

        int mu_reco  = it_mu->second;
        int pi2_reco = it_pi2->second;
        if (mu_reco  < 0 || mu_reco  >= (int)RecoParticles.size()) continue;
        if (pi2_reco < 0 || pi2_reco >= (int)RecoParticles.size()) continue;

        ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> mu_rp  = {RecoParticles[mu_reco]};
        ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> pi2_rp = {RecoParticles[pi2_reco]};
        auto mu_ts  = FCCAnalyses::ReconstructedParticle2Track::getRP2TRK(mu_rp,  AllTracks);
        auto pi2_ts = FCCAnalyses::ReconstructedParticle2Track::getRP2TRK(pi2_rp, AllTracks);

        if (mu_ts.empty() || pi2_ts.empty()) continue;
        result.push_back((int)i);
    }
    return result;
}

// Invariant mass (GeV) of a two-body system under the muon mass hypothesis for both particles.
// idx1 and idx2 are MC particle indices; both are assigned m_mu regardless of their true identity.
// Used to compare: idx1=pi (decaying pion treated as mu) vs idx1=mu (true muon) — both paired with pi2.
ROOT::VecOps::RVec<float> get_pimu_invM_mumu_hyp(
    const ROOT::VecOps::RVec<int>& idx1,
    const ROOT::VecOps::RVec<int>& idx2,
    const ROOT::VecOps::RVec<edm4hep::MCParticleData>& particles)
{
    const float m_mu = 0.10566f; // GeV
    ROOT::VecOps::RVec<float> result;
    size_t n = std::min(idx1.size(), idx2.size());
    result.reserve(n);
    for (size_t i = 0; i < n; ++i) {
        int i1 = idx1[i], i2 = idx2[i];
        if (i1 < 0 || i2 < 0) { result.push_back(-1.f); continue; }
        auto& p1 = particles[i1].momentum;
        auto& p2 = particles[i2].momentum;
        float p1sq = p1.x*p1.x + p1.y*p1.y + p1.z*p1.z;
        float p2sq = p2.x*p2.x + p2.y*p2.y + p2.z*p2.z;
        float E1 = std::sqrt(p1sq + m_mu*m_mu);
        float E2 = std::sqrt(p2sq + m_mu*m_mu);
        float px = p1.x + p2.x, py = p1.y + p2.y, pz = p1.z + p2.z;
        float M2 = (E1+E2)*(E1+E2) - px*px - py*py - pz*pz;
        result.push_back(M2 > 0.f ? std::sqrt(M2) : 0.f);
    }
    return result;
}


// Step 1: for each V0, extract the Ks MC index from the nested mother-index table.
// Returns -1 for V0s not matched to a Ks (PDG 310) mother.
struct extract_V0_ks_idx {
    ROOT::VecOps::RVec<int> operator()(
        const std::vector<std::vector<int>>& V0_mother_idx,
        const std::vector<std::vector<int>>& V0_mother_pdg)
    {
        ROOT::VecOps::RVec<int> result;
        result.reserve(V0_mother_idx.size());
        for (size_t v = 0; v < V0_mother_idx.size(); ++v) {
            int idx = -1;
            if (v < V0_mother_pdg.size()) {
                const auto& mothers = V0_mother_idx[v];
                const auto& mpdgs   = V0_mother_pdg[v];
                for (size_t t = 0; t < mothers.size() && t < mpdgs.size(); ++t) {
                    if (std::abs(mpdgs[t]) == 310) { idx = mothers[t]; break; }
                }
            }
            result.push_back(idx);
        }
        return result;
    }
};

// Step 2: for each V0, look up the opening angle between the decaying pion and the muon
// using the Ks MC index extracted by extract_V0_ks_idx. Returns -1 if no match.
ROOT::VecOps::RVec<float> get_V0_pimu_opening_angle(
    const ROOT::VecOps::RVec<int>& V0_ks_idx,
    const PiMuFromKs& pimu_cands,
    const ROOT::VecOps::RVec<edm4hep::MCParticleData>& particles)
{
    // Build map: Ks MC index -> opening angle between decaying pi and muon
    std::unordered_map<int, float> ks_to_angle;
    size_t n = std::min({pimu_cands.ks_idx.size(), pimu_cands.pi_idx.size(), pimu_cands.mu_idx.size()});
    for (size_t i = 0; i < n; ++i) {
        int ks_i = pimu_cands.ks_idx[i], pi_i = pimu_cands.pi_idx[i], mu_i = pimu_cands.mu_idx[i];
        if (ks_i < 0 || pi_i < 0 || mu_i < 0) continue;
        if (pi_i >= (int)particles.size() || mu_i >= (int)particles.size()) continue;
        const auto& p_pi = particles[pi_i].momentum;
        const auto& p_mu = particles[mu_i].momentum;
        float mag_pi = std::sqrt(p_pi.x*p_pi.x + p_pi.y*p_pi.y + p_pi.z*p_pi.z);
        float mag_mu = std::sqrt(p_mu.x*p_mu.x + p_mu.y*p_mu.y + p_mu.z*p_mu.z);
        if (mag_pi <= 0.f || mag_mu <= 0.f) continue;
        float cosine = (p_pi.x*p_mu.x + p_pi.y*p_mu.y + p_pi.z*p_mu.z) / (mag_pi * mag_mu);
        cosine = std::max(-1.f, std::min(1.f, cosine));
        ks_to_angle[ks_i] = std::acos(cosine);
    }

    ROOT::VecOps::RVec<float> result;
    result.reserve(V0_ks_idx.size());
    for (int ks_i : V0_ks_idx) {
        if (ks_i < 0) { result.push_back(-1.f); continue; }
        auto it = ks_to_angle.find(ks_i);
        result.push_back(it != ks_to_angle.end() ? it->second : -1.f);
    }
    return result;
}


// Return the momentum magnitude of one track (by index 0 or 1) for each fitted pimu vertex.
ROOT::VecOps::RVec<float> get_pimu_vtx_track_p(
    const ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex>& vtx,
    int track_idx)
{
    ROOT::VecOps::RVec<float> result;
    for (const auto& v : vtx) {
        if (track_idx < (int)v.updated_track_momentum_at_vertex.size())
            result.push_back(v.updated_track_momentum_at_vertex[track_idx].Mag());
        else
            result.push_back(-1.f);
    }
    return result;
}


ROOT::VecOps::RVec<float> get_pimu_vtx_p_min(
    const ROOT::VecOps::RVec<float>& p_mu,
    const ROOT::VecOps::RVec<float>& p_pi)
{
    ROOT::VecOps::RVec<float> result;
    size_t n = std::min(p_mu.size(), p_pi.size());
    result.reserve(n);
    for (size_t i = 0; i < n; ++i)
        result.push_back(std::min(p_mu[i], p_pi[i]));
    return result;
}


ROOT::VecOps::RVec<float> get_MC_vertex_x_from_indices(
    const ROOT::VecOps::RVec<edm4hep::MCParticleData>& mcparticle,
    const ROOT::VecOps::RVec<int>& indices) {

    ROOT::VecOps::RVec<float> result;

    for (int idx : indices) {
        result.push_back(mcparticle[idx].vertex.x);
    }

    return result;
}

ROOT::VecOps::RVec<float> get_MC_vertex_y_from_indices(
    const ROOT::VecOps::RVec<edm4hep::MCParticleData>& mcparticle,
    const ROOT::VecOps::RVec<int>& indices) {

    ROOT::VecOps::RVec<float> result;

    for (int idx : indices) {
        result.push_back(mcparticle[idx].vertex.y);
    }

    return result;
}

ROOT::VecOps::RVec<float> get_MC_vertex_z_from_indices(
    const ROOT::VecOps::RVec<edm4hep::MCParticleData>& mcparticle,
    const ROOT::VecOps::RVec<int>& indices) {

    ROOT::VecOps::RVec<float> result;

    for (int idx : indices) {
        result.push_back(mcparticle[idx].vertex.z);
    }

    return result;
}

ROOT::VecOps::RVec<float> get_MC_pt_from_indices(
    const ROOT::VecOps::RVec<edm4hep::MCParticleData>& mcparticle,
    const ROOT::VecOps::RVec<int>& indices) {

    ROOT::VecOps::RVec<float> result;

    for (int idx : indices) {

        float px = mcparticle[idx].momentum.x;
        float py = mcparticle[idx].momentum.y;

        result.push_back(std::sqrt(px*px + py*py));
    }

    return result;
}

ROOT::VecOps::RVec<float> get_pi_mu_opening_angle(
    const ROOT::VecOps::RVec<edm4hep::MCParticleData>& mcparticle,
    const ROOT::VecOps::RVec<int>& mu_idx,
    const ROOT::VecOps::RVec<int>& pi_idx) {

    ROOT::VecOps::RVec<float> result;

    for (size_t i = 0; i < mu_idx.size(); ++i) {

        const auto& mu = mcparticle[mu_idx[i]];
        const auto& pi = mcparticle[pi_idx[i]];

        float mu_px = mu.momentum.x;
        float mu_py = mu.momentum.y;
        float mu_pz = mu.momentum.z;

        float pi_px = pi.momentum.x;
        float pi_py = pi.momentum.y;
        float pi_pz = pi.momentum.z;

        float dot =
            mu_px*pi_px +
            mu_py*pi_py +
            mu_pz*pi_pz;

        float mu_p =
            std::sqrt(mu_px*mu_px +
                      mu_py*mu_py +
                      mu_pz*mu_pz);

        float pi_p =
            std::sqrt(pi_px*pi_px +
                      pi_py*pi_py +
                      pi_pz*pi_pz);

        float cos_theta = dot/(mu_p*pi_p);

        // protection against numerical precision
        cos_theta = std::max(-1.f, std::min(1.f, cos_theta));

        result.push_back(std::acos(cos_theta));
    }

    return result;
}

ROOT::VecOps::RVec<float> get_MC_p_from_indices(
    const ROOT::VecOps::RVec<edm4hep::MCParticleData>& parts,
    const ROOT::VecOps::RVec<int>& indices) {

    ROOT::VecOps::RVec<float> result;
    result.reserve(indices.size());

    for (int idx : indices) {
        if (idx >= 0 && idx < (int)parts.size()) {

            const auto& p = parts[idx];

            float px = p.momentum.x;
            float py = p.momentum.y;
            float pz = p.momentum.z;

            result.push_back(std::sqrt(px*px + py*py + pz*pz));
        } else {
            result.push_back(-1.0); // safety flag
        }
    }

    return result;
}


// 1. Momentos de las hijas
struct get_nested_daughter_p {
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> operator()(
        const ROOT::VecOps::RVec<edm4hep::MCParticleData>& mc_particles,
        const std::vector<std::vector<int>>& daughter_indices)
    {
        ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> result;
        for (const auto& indices : daughter_indices) {
            ROOT::VecOps::RVec<float> vtx_p;
            for (int idx : indices) {
                if (idx >= 0 && idx < (int)mc_particles.size()) {
                    auto& p = mc_particles.at(idx).momentum;
                    float p_mag = std::sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
                    vtx_p.push_back(p_mag);
                } else {
                    vtx_p.push_back(-999.0);
                }
            }
            result.push_back(vtx_p);
        }
        return result;
    }
};

struct get_nested_daughter_px {
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> operator()(
        const ROOT::VecOps::RVec<edm4hep::MCParticleData>& mc_particles,
        const std::vector<std::vector<int>>& daughter_indices)
    {
        ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> result;
        for (const auto& indices : daughter_indices) {
            ROOT::VecOps::RVec<float> vals;
            for (int idx : indices) {
                vals.push_back(idx >= 0 && idx < (int)mc_particles.size() ? mc_particles.at(idx).momentum.x : -999.0f);
            }
            result.push_back(vals);
        }
        return result;
    }
};

struct get_nested_daughter_py {
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> operator()(
        const ROOT::VecOps::RVec<edm4hep::MCParticleData>& mc_particles,
        const std::vector<std::vector<int>>& daughter_indices)
    {
        ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> result;
        for (const auto& indices : daughter_indices) {
            ROOT::VecOps::RVec<float> vals;
            for (int idx : indices) {
                vals.push_back(idx >= 0 && idx < (int)mc_particles.size() ? mc_particles.at(idx).momentum.y : -999.0f);
            }
            result.push_back(vals);
        }
        return result;
    }
};

struct get_nested_daughter_pz {
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> operator()(
        const ROOT::VecOps::RVec<edm4hep::MCParticleData>& mc_particles,
        const std::vector<std::vector<int>>& daughter_indices)
    {
        ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> result;
        for (const auto& indices : daughter_indices) {
            ROOT::VecOps::RVec<float> vals;
            for (int idx : indices) {
                vals.push_back(idx >= 0 && idx < (int)mc_particles.size() ? mc_particles.at(idx).momentum.z : -999.0f);
            }
            result.push_back(vals);
        }
        return result;
    }
};

// Momentos de las madres
struct get_nested_mother_p {
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> operator()(
        const ROOT::VecOps::RVec<edm4hep::MCParticleData>& mc_particles,
        const std::vector<std::vector<int>>& mother_indices)
    {
        ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> result;
        for (const auto& indices : mother_indices) {
            ROOT::VecOps::RVec<float> vtx_p;
            for (int idx : indices) {
                if (idx >= 0 && idx < (int)mc_particles.size()) {
                    auto& p = mc_particles.at(idx).momentum;
                    float p_mag = std::sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
                    vtx_p.push_back(p_mag);
                } else {
                    vtx_p.push_back(-999.0);
                }
            }
            result.push_back(vtx_p);
        }
        return result;
    }
};

struct get_nested_mother_px {
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> operator()(
        const ROOT::VecOps::RVec<edm4hep::MCParticleData>& mc_particles,
        const std::vector<std::vector<int>>& mother_indices)
    {
        ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> result;
        for (const auto& indices : mother_indices) {
            ROOT::VecOps::RVec<float> vals;
            for (int idx : indices) {
                vals.push_back(idx >= 0 && idx < (int)mc_particles.size() ? mc_particles.at(idx).momentum.x : -999.0f);
            }
            result.push_back(vals);
        }
        return result;
    }
};

struct get_nested_mother_py {
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> operator()(
        const ROOT::VecOps::RVec<edm4hep::MCParticleData>& mc_particles,
        const std::vector<std::vector<int>>& mother_indices)
    {
        ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> result;
        for (const auto& indices : mother_indices) {
            ROOT::VecOps::RVec<float> vals;
            for (int idx : indices) {
                vals.push_back(idx >= 0 && idx < (int)mc_particles.size() ? mc_particles.at(idx).momentum.y : -999.0f);
            }
            result.push_back(vals);
        }
        return result;
    }
};

struct get_nested_mother_pz {
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> operator()(
        const ROOT::VecOps::RVec<edm4hep::MCParticleData>& mc_particles,
        const std::vector<std::vector<int>>& mother_indices)
    {
        ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> result;
        for (const auto& indices : mother_indices) {
            ROOT::VecOps::RVec<float> vals;
            for (int idx : indices) {
                vals.push_back(idx >= 0 && idx < (int)mc_particles.size() ? mc_particles.at(idx).momentum.z : -999.0f);
            }
            result.push_back(vals);
        }
        return result;
    }
};

// 3. PDG de las abuelas
struct get_nested_gmother_pdg {
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>> operator()(
        const ROOT::VecOps::RVec<edm4hep::MCParticleData>& mc_particles,
        const ROOT::VecOps::RVec<int>& parents_ind,
        const std::vector<std::vector<int>>& mother_indices) // <-- Cambiado aquí
    {
        ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>> result;
        for (const auto& indices : mother_indices) {
            ROOT::VecOps::RVec<int> vtx_gmother_pdg;
            for (int m_idx : indices) {
                if (m_idx >= 0 && m_idx < (int)mc_particles.size()) {
                    const auto& mother = mc_particles.at(m_idx);
                    if (mother.parents_end > mother.parents_begin && mother.parents_begin < parents_ind.size()) {
                        int gm_idx = parents_ind.at(mother.parents_begin);
                        if (gm_idx >= 0 && gm_idx < (int)mc_particles.size()) {
                            vtx_gmother_pdg.push_back(mc_particles.at(gm_idx).PDG);
                        } else {
                            vtx_gmother_pdg.push_back(-999);
                        }
                    } else {
                        vtx_gmother_pdg.push_back(0);
                    }
                } else {
                    vtx_gmother_pdg.push_back(-999);
                }
            }
            result.push_back(vtx_gmother_pdg);
        }
        return result;
    }
};

// Pure generator-truth kinematics for KS -> pi1 pi2, pi1 -> mu nu decays.
// For every KS with exactly two pion daughters where exactly one of them
// decays in-flight to mu+nu, stores the true muon 4-momentum, the sibling
// pion's true 4-momentum, its energy recomputed under the muon-mass
// hypothesis (3-momentum unchanged), the resulting mu-mu invariant mass
// under that hypothesis, and the true KS 4-momentum. No reco.
struct PimunuGenTruth {
    ROOT::VecOps::RVec<float> KS_p;
    ROOT::VecOps::RVec<float> m_mumu_muhyp;
};

PimunuGenTruth get_KS_pimunu_gentruth(
    const ROOT::VecOps::RVec<edm4hep::MCParticleData>& mc,
    const ROOT::VecOps::RVec<int>& daughter_ind) // Particle1 = Particle#1.index
{
    PimunuGenTruth res;
    constexpr double m_mu = 0.1056583755; // PDG muon mass [GeV]

    auto find_muon_daughter = [&](int idx) -> int {
        const auto& p = mc.at(idx);
        for (unsigned j = p.daughters_begin; j < p.daughters_end; ++j) {
            if (j >= daughter_ind.size()) continue;
            int dau = daughter_ind.at(j);
            if (dau >= 0 && (size_t)dau < mc.size() && std::abs(mc.at(dau).PDG) == 13) return dau;
        }
        return -1;
    };

    for (size_t i = 0; i < mc.size(); ++i) {
        if (mc[i].PDG != 310) continue; // KS
        const auto& ks = mc[i];
        if (ks.daughters_end - ks.daughters_begin != 2) continue;
        if (ks.daughters_begin + 1 >= daughter_ind.size()) continue;

        int d0 = daughter_ind.at(ks.daughters_begin);
        int d1 = daughter_ind.at(ks.daughters_begin + 1);
        if (d0 < 0 || d1 < 0 || (size_t)d0 >= mc.size() || (size_t)d1 >= mc.size()) continue;
        if (std::abs(mc[d0].PDG) != 211 || std::abs(mc[d1].PDG) != 211) continue; // both pions

        int mu_from_d0 = find_muon_daughter(d0);
        int mu_from_d1 = find_muon_daughter(d1);

        int pi2_idx = -1, mu_idx = -1;
        if (mu_from_d0 >= 0 && mu_from_d1 < 0)      { pi2_idx = d1; mu_idx = mu_from_d0; }
        else if (mu_from_d1 >= 0 && mu_from_d0 < 0) { pi2_idx = d0; mu_idx = mu_from_d1; }
        else continue; // neither or both decayed to mu nu -> not this channel / ambiguous

        const auto& mu  = mc.at(mu_idx);
        const auto& pi2 = mc.at(pi2_idx);

        double mu_p2 = mu.momentum.x*mu.momentum.x + mu.momentum.y*mu.momentum.y + mu.momentum.z*mu.momentum.z;
        double mu_E  = std::sqrt(mu_p2 + mu.mass*mu.mass);

        double pi2_p2      = pi2.momentum.x*pi2.momentum.x + pi2.momentum.y*pi2.momentum.y + pi2.momentum.z*pi2.momentum.z;
        double pi2_E_muhyp = std::sqrt(pi2_p2 + m_mu*m_mu);

        double Esum  = mu_E + pi2_E_muhyp;
        double pxsum = mu.momentum.x + pi2.momentum.x;
        double pysum = mu.momentum.y + pi2.momentum.y;
        double pzsum = mu.momentum.z + pi2.momentum.z;
        double m2    = Esum*Esum - (pxsum*pxsum + pysum*pysum + pzsum*pzsum);
        double m_mumu = (m2 > 0) ? std::sqrt(m2) : -std::sqrt(-m2);

        double ks_p2 = ks.momentum.x*ks.momentum.x + ks.momentum.y*ks.momentum.y + ks.momentum.z*ks.momentum.z;

        res.KS_p.push_back(std::sqrt(ks_p2));
        res.m_mumu_muhyp.push_back((float)m_mumu);
    }
    return res;
}

ROOT::VecOps::RVec<float> get_pimunu_KS_p(PimunuGenTruth t)        { return t.KS_p; }
ROOT::VecOps::RVec<float> get_pimunu_m_mumu_muhyp(PimunuGenTruth t){ return t.m_mumu_muhyp; }

// Pure generator-truth m_mumu for genuine (non-decayed) KS -> pi+ pi- decays,
// under the muon-mass hypothesis applied to BOTH true pion 3-momenta. This is
// the truth-level analogue of the "pipi misID" background (both legs
// misidentified as muons), for direct comparison against the pimu DIF edge.
struct PipiGenTruth {
    ROOT::VecOps::RVec<float> m_mumu_muhyp;
};

PipiGenTruth get_KS_pipi_gentruth(
    const ROOT::VecOps::RVec<edm4hep::MCParticleData>& mc,
    const ROOT::VecOps::RVec<int>& daughter_ind) // Particle1 = Particle#1.index
{
    PipiGenTruth res;
    constexpr double m_mu = 0.1056583755; // PDG muon mass [GeV]

    auto has_muon_daughter = [&](int idx) -> bool {
        const auto& p = mc.at(idx);
        for (unsigned j = p.daughters_begin; j < p.daughters_end; ++j) {
            if (j >= daughter_ind.size()) continue;
            int dau = daughter_ind.at(j);
            if (dau >= 0 && (size_t)dau < mc.size() && std::abs(mc.at(dau).PDG) == 13) return true;
        }
        return false;
    };

    for (size_t i = 0; i < mc.size(); ++i) {
        if (mc[i].PDG != 310) continue; // KS
        const auto& ks = mc[i];
        if (ks.daughters_end - ks.daughters_begin != 2) continue;
        if (ks.daughters_begin + 1 >= daughter_ind.size()) continue;

        int d0 = daughter_ind.at(ks.daughters_begin);
        int d1 = daughter_ind.at(ks.daughters_begin + 1);
        if (d0 < 0 || d1 < 0 || (size_t)d0 >= mc.size() || (size_t)d1 >= mc.size()) continue;
        if (std::abs(mc[d0].PDG) != 211 || std::abs(mc[d1].PDG) != 211) continue; // both pions

        if (has_muon_daughter(d0) || has_muon_daughter(d1)) continue; // require genuine (undecayed) pipi

        const auto& p0 = mc.at(d0);
        const auto& p1 = mc.at(d1);

        double p0_2 = p0.momentum.x*p0.momentum.x + p0.momentum.y*p0.momentum.y + p0.momentum.z*p0.momentum.z;
        double p1_2 = p1.momentum.x*p1.momentum.x + p1.momentum.y*p1.momentum.y + p1.momentum.z*p1.momentum.z;
        double E0 = std::sqrt(p0_2 + m_mu*m_mu);
        double E1 = std::sqrt(p1_2 + m_mu*m_mu);

        double Esum  = E0 + E1;
        double pxsum = p0.momentum.x + p1.momentum.x;
        double pysum = p0.momentum.y + p1.momentum.y;
        double pzsum = p0.momentum.z + p1.momentum.z;
        double m2    = Esum*Esum - (pxsum*pxsum + pysum*pysum + pzsum*pzsum);
        double m_mumu = (m2 > 0) ? std::sqrt(m2) : -std::sqrt(-m2);

        res.m_mumu_muhyp.push_back((float)m_mumu);
    }
    return res;
}

ROOT::VecOps::RVec<float> get_pipi_m_mumu_muhyp(PipiGenTruth t) { return t.m_mumu_muhyp; }

// ─── Pion decay-in-tracker: muon daughter track lookup ────────────────────
// For each V0 (indexed by V0_matched_daughters_idx), for each truth-matched
// pion leg, check whether that pion has a muon child whose production vertex
// is inside the tracker (vertex_r < 1500 mm, |vertex_z| < 2300 mm). If so,
// look up the muon's reconstructed track via MCRecoAssociations and return
// its raw (pre-vertex-fit) parameters.
//
// Layout of all output vectors: [V0 index][leg index (0 or 1)].
// Zero-filled when the pion did not decay in the tracker, or when the muon
// could not be matched to a reconstructed track.
//
// Parameter naming:
//   Particle1   = Particle#1.index  (parent-link table: for particle i,
//                 Particle[i].parents_begin..parents_end index into Particle1
//                 to retrieve parent indices in Particle)
//   MCRecoAssoc0 = MCRecoAssociations#0.index  (reco-particle index)
//   MCRecoAssoc1 = MCRecoAssociations#1.index  (MC-particle index)
struct V0PionDecayInfo {
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> muon_r;         // muon vertex_r [mm]
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> muon_omega;     // track curvature [1/mm]
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> muon_omega_cov; // covMatrix[5] = sigma^2(omega)
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> muon_tL;        // tanLambda
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> muon_tL_cov;    // covMatrix[14] = sigma^2(tanLambda)
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> muon_omega_tL_cov; // covMatrix[12] = cov(omega,tanLambda)
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> muon_D0;        // impact parameter [mm]
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> muon_D0_cov;    // covMatrix[0] = sigma^2(D0)
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> muon_Z0;        // longitudinal impact parameter [mm]
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> muon_Z0_cov;    // covMatrix[9] = sigma^2(Z0)
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> muon_px;        // reco momentum x [GeV]
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> muon_py;
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> muon_pz;
};

V0PionDecayInfo get_V0_pion_decay_info(
    const std::vector<std::vector<int>>&  daughters_idx,
    const ROOT::VecOps::RVec<edm4hep::MCParticleData>&  Particle,
    const ROOT::VecOps::RVec<int>&                      Particle1,
    const ROOT::VecOps::RVec<int>&                      MCRecoAssoc0,
    const ROOT::VecOps::RVec<int>&                      MCRecoAssoc1,
    const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& reco,
    const ROOT::VecOps::RVec<edm4hep::TrackState>&      tracks
) {
    // Build reverse map: children[i] = particles whose parent is i
    int N = (int)Particle.size();
    std::vector<std::vector<int>> children(N);
    for (int i = 0; i < N; ++i) {
        for (int j = (int)Particle[i].parents_begin; j < (int)Particle[i].parents_end; ++j) {
            if (j >= 0 && j < (int)Particle1.size()) {
                int par = Particle1[j];
                if (par >= 0 && par < N) children[par].push_back(i);
            }
        }
    }

    // Build MC-index → reco-particle-index map (last association wins)
    std::unordered_map<int,int> mc2reco;
    for (int k = 0; k < (int)MCRecoAssoc1.size(); ++k)
        mc2reco[MCRecoAssoc1[k]] = MCRecoAssoc0[k];

    V0PionDecayInfo out;

    for (const auto& v0_legs : daughters_idx) {
        ROOT::VecOps::RVec<float> r_v, om_v, omC_v, tL_v, tLC_v, omtLC_v;
        ROOT::VecOps::RVec<float> D0_v, D0C_v, Z0_v, Z0C_v, px_v, py_v, pz_v;
        for (int pi_idx : v0_legs) {
            float mu_r=0.f, mu_om=0.f, mu_omC=0.f, mu_tL=0.f, mu_tLC=0.f, mu_omtLC=0.f;
            float mu_D0=0.f, mu_D0C=0.f, mu_Z0=0.f, mu_Z0C=0.f;
            float mu_px=0.f, mu_py=0.f, mu_pz=0.f;

            if (pi_idx >= 0 && pi_idx < N) {
                int mu_mc = -1;
                for (int ch : children[pi_idx]) {
                    if (std::abs(Particle[ch].PDG) != 13) continue;
                    float vx = Particle[ch].vertex.x;
                    float vy = Particle[ch].vertex.y;
                    float vz = Particle[ch].vertex.z;
                    float vr = std::sqrt(vx*vx + vy*vy);
                    if (vr > 1.f && vr < 1500.f && std::abs(vz) < 2300.f) {
                        mu_r  = vr;
                        mu_mc = ch;
                        break;
                    }
                }
                if (mu_mc >= 0) {
                    auto it = mc2reco.find(mu_mc);
                    if (it != mc2reco.end()) {
                        int ri = it->second;
                        if (ri >= 0 && ri < (int)reco.size()) {
                            int ti = reco[ri].tracks_begin;
                            if (ti >= 0 && ti < (int)tracks.size()) {
                                mu_om    = tracks[ti].omega;
                                mu_omC   = tracks[ti].covMatrix[5];
                                mu_tL    = tracks[ti].tanLambda;
                                mu_tLC   = tracks[ti].covMatrix[14];
                                mu_omtLC = tracks[ti].covMatrix[12];
                                mu_D0    = tracks[ti].D0;
                                mu_D0C   = tracks[ti].covMatrix[0];
                                mu_Z0    = tracks[ti].Z0;
                                mu_Z0C   = tracks[ti].covMatrix[9];
                            }
                            const auto& m = reco[ri].momentum;
                            mu_px = m.x; mu_py = m.y; mu_pz = m.z;
                        }
                    }
                }
            }
            r_v.push_back(mu_r);       om_v.push_back(mu_om);     omC_v.push_back(mu_omC);
            tL_v.push_back(mu_tL);     tLC_v.push_back(mu_tLC);   omtLC_v.push_back(mu_omtLC);
            D0_v.push_back(mu_D0);     D0C_v.push_back(mu_D0C);
            Z0_v.push_back(mu_Z0);     Z0C_v.push_back(mu_Z0C);
            px_v.push_back(mu_px);     py_v.push_back(mu_py);     pz_v.push_back(mu_pz);
        }
        out.muon_r.push_back(r_v);        out.muon_omega.push_back(om_v);
        out.muon_omega_cov.push_back(omC_v);
        out.muon_tL.push_back(tL_v);      out.muon_tL_cov.push_back(tLC_v);
        out.muon_omega_tL_cov.push_back(omtLC_v);
        out.muon_D0.push_back(D0_v);      out.muon_D0_cov.push_back(D0C_v);
        out.muon_Z0.push_back(Z0_v);      out.muon_Z0_cov.push_back(Z0C_v);
        out.muon_px.push_back(px_v);      out.muon_py.push_back(py_v);
        out.muon_pz.push_back(pz_v);
    }
    return out;
}

// Flat accessors for the V0PionDecayInfo struct (one per field)
ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> get_V0_muon_r           (V0PionDecayInfo d) { return d.muon_r; }
ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> get_V0_muon_omega       (V0PionDecayInfo d) { return d.muon_omega; }
ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> get_V0_muon_omega_cov   (V0PionDecayInfo d) { return d.muon_omega_cov; }
ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> get_V0_muon_tL          (V0PionDecayInfo d) { return d.muon_tL; }
ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> get_V0_muon_tL_cov      (V0PionDecayInfo d) { return d.muon_tL_cov; }
ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> get_V0_muon_omega_tL_cov(V0PionDecayInfo d) { return d.muon_omega_tL_cov; }
ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> get_V0_muon_D0          (V0PionDecayInfo d) { return d.muon_D0; }
ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> get_V0_muon_D0_cov      (V0PionDecayInfo d) { return d.muon_D0_cov; }
ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> get_V0_muon_Z0          (V0PionDecayInfo d) { return d.muon_Z0; }
ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> get_V0_muon_Z0_cov      (V0PionDecayInfo d) { return d.muon_Z0_cov; }
ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> get_V0_muon_px          (V0PionDecayInfo d) { return d.muon_px; }
ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> get_V0_muon_py          (V0PionDecayInfo d) { return d.muon_py; }
ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> get_V0_muon_pz          (V0PionDecayInfo d) { return d.muon_pz; }

}}

#endif
