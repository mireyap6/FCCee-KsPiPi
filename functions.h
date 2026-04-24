#ifndef ZHfunctions_H
#define ZHfunctions_H

#include <cmath>
#include <vector>
#include <math.h>

#include "TLorentzVector.h"
#include "ROOT/RVec.hxx"
#include "edm4hep/ReconstructedParticleData.h"
#include "edm4hep/MCParticleData.h"
#include "edm4hep/ParticleIDData.h"
#include "ReconstructedParticle2MC.h"


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
 
 
 
// returns missing energy vector, based on reco particles
Vec_rp missingEnergy(float ecm, Vec_rp in, float p_cutoff = 0.0) {
    float px = 0, py = 0, pz = 0, e = 0;
    for(auto &p : in) {
        if (std::sqrt(p.momentum.x * p.momentum.x + p.momentum.y*p.momentum.y) < p_cutoff) continue;
        px += -p.momentum.x;
        py += -p.momentum.y;
        pz += -p.momentum.z;
        e += p.energy;
    }
    
    Vec_rp ret;
    rp res;
    res.momentum.x = px;
    res.momentum.y = py;
    res.momentum.z = pz;
    res.energy = ecm-e;
    ret.emplace_back(res);
    return ret;
}

// calculate the cosine(theta) of the missing energy vector
float get_cosTheta_miss(Vec_rp met){
    
    float costheta = 0.;
    if(met.size() > 0) {
        
        TLorentzVector lv_met;
        lv_met.SetPxPyPzE(met[0].momentum.x, met[0].momentum.y, met[0].momentum.z, met[0].energy);
        costheta = fabs(std::cos(lv_met.Theta()));
    }
    return costheta;
}

 
ROOT::VecOps::RVec<TLorentzVector> build_p4(ROOT::VecOps::RVec<float> px, ROOT::VecOps::RVec<float> py, ROOT::VecOps::RVec<float> pz, ROOT::VecOps::RVec<float> mass) {
    ROOT::VecOps::RVec<TLorentzVector> p4;
    for (size_t i = 0; i < px.size(); ++i) {
        TLorentzVector tlv;
        tlv.SetXYZM(px[i], py[i], pz[i], mass[i]);
        p4.push_back(tlv);
    }
    return p4;
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

ROOT::VecOps::RVec<int> get_MC_Vertex_isKSpipi(
    const std::vector<std::vector<int>>& PDGmother, 
    const std::vector<std::vector<int>>& PDG)
{
    ROOT::VecOps::RVec<int> result;
    for (size_t i = 0; i < PDGmother.size(); ++i) {
        int isKS = 0;
        for (size_t j = 0; j < PDGmother[i].size(); ++j) {
            if (std::abs(PDGmother[i][j]) == 310) {
                isKS = 10; 
                int charged_daughters = 0;
                for (size_t k = 0; k < PDG[i].size(); ++k) {
                    int pdg_id = std::abs(PDG[i][k]);
                    // accept both pions and muons
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
        double beamspotX = 0, double beamspotY = 0, double beamspotZ = 0){

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

        // define dummy vertex in case the fit cannot be performed
        edm4hep::VertexData dummyVertex;
        dummyVertex.chi2 = -1;
        //dummyVertex.ndf = 0;
        dummyVertex.position = edm4hep::Vector3f(beamspotX, beamspotY, beamspotZ);
        FCCAnalyses::VertexingUtils::FCCAnalysesVertex dummyVertexObject;
        dummyVertexObject.vertex = dummyVertex;
        dummyVertexObject.ntracks = 0;
        dummyVertexObject.mc_ind = -1;
        if( tracks.size() < 2 ){ return dummyVertexObject; }

        // call primary vertex finder from FCCAnalyses
        FCCAnalyses::VertexingUtils::FCCAnalysesVertex vertex;
        vertex = FCCAnalyses::VertexFitterSimple::VertexFitter_Tk(
            1, tracks, doBeamSpotConstraint,
            sigma_beamspotX, sigma_beamspotY, sigma_beamspotZ,
            beamspotX, beamspotY, beamspotZ
        );
        return vertex;
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
};

inline MatchMCV0 match_MC2Reco_V0s(
    const ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex>& recoVertices,
    const ROOT::VecOps::RVec<int>& rp2mc, // association indices that point to the mc particle collection
    const ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertexMC>& mcvertices,
    const ROOT::VecOps::RVec<float>& MC_Vertex_x, 
    const ROOT::VecOps::RVec<float>& MC_Vertex_y,
    const ROOT::VecOps::RVec<float>& MC_Vertex_z
) {
    MatchMCV0 result;

    //loop over all reconstructed vertices
    for (const auto& vtx : recoVertices) {
        int matched_vtx_idx = -1;
        float dist = -999.0;
        
        //initialise vector of associated MC particles for the vertex
        std::vector<int> asso_vtx_mc_particles;
        
        //loop over reco particles that form the reco vertices
        for (int rp_idx : vtx.reco_ind) {
            //get the MC particle index associated to the reco particle
            if (rp_idx >= 0 && rp_idx < rp2mc.size()) {
                int mc_p_id = rp2mc.at(rp_idx);
                if (mc_p_id >= 0) {
                    asso_vtx_mc_particles.push_back(mc_p_id);
                }
            }
        }
        //once two associated MC particles are found, look for them in the MC vertex object
        if (asso_vtx_mc_particles.size() >= 2) {
            for (size_t i = 0; i < mcvertices.size(); ++i) {
                int matches = 0;
                const auto& mc_ind = mcvertices.at(i).mc_ind; //this gives the indices of mc particles associated to the MC vertex i
                
                for (int mc_p_id : asso_vtx_mc_particles) {
                    if (std::find(mc_ind.begin(), mc_ind.end(), mc_p_id) != mc_ind.end()) {
                        matches++;
                    }
                }

                if (matches >= 2) {
                    matched_vtx_idx = i;
                    float dx = vtx.vertex.position.x - MC_Vertex_x[i];
                    float dy = vtx.vertex.position.y - MC_Vertex_y[i];
                    float dz = vtx.vertex.position.z - MC_Vertex_z[i];
                    dist = std::sqrt(dx*dx + dy*dy + dz*dz);
                    break;
                }
            }
        }
        result.mc_vtx_idx.push_back(matched_vtx_idx);
        result.distance.push_back(dist);
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

}}

#endif

