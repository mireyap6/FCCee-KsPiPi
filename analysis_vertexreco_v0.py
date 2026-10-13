#Mandatory: List of processes
#inputDir = "/eos/user/m/mpandopo/TestDataPythia"


processList = {
                'p8_ee_Zss_ecm91_Ksmumu':{'chunks':10, 'fraction':1},
                'p8_ee_Zcc_ecm91_Ksmumu':{'chunks':10, 'fraction':1},
                'p8_ee_Zbb_ecm91_Ksmumu':{'chunks':10, 'fraction':1},
                'p8_ee_Zud_ecm91_Ksmumu':{'chunks':10, 'fraction':1},
}

#Mandatory: Production tag when running over EDM4Hep centrally produced events, this points to the yaml files for getting sample statistics
prodTag     = "FCCee/winter2023/IDEA/"

#Optional: output directory, default is local running directory
outputDir = "/afs/cern.ch/work/m/mpandopo/Sensitivity_1M"

#EOS output directory for batch jobs
#outputDirEos = "/eos/user/m/mpandopo"

#Optional

nCPUS       = 8
runBatch    = True
batchQueue = "workday"
compGroup = "group_u_FCC.local_gen"

includePaths = ["functions.h"]

#Mandatory: RDFanalysis class where the use defines the operations on the TTree
class RDFanalysis():

    #__________________________________________________________
    #Mandatory: analysers funtion to define the analysers to process, please make sure you return the last dataframe, in this example it is df2
    def analysers(df):
        df2 = (df
               .Alias("Particle0", "Particle#0.index")  #parent particles 
               .Alias("Particle1", "Particle#1.index")  #daugther particles 
               .Alias("MCRecoAssociations0", "MCRecoAssociations#0.index") #points to the ReconstructedParticles collection
               .Alias("MCRecoAssociations1", "MCRecoAssociations#1.index") #points to the Particle collection
               .Alias("Track1", "EFlowTrack#1.index")
               .Alias("TrackState", "EFlowTrack_1")
               .Alias("Muon0", "Muon#0.index")
      #         .Alias("Particle0", "_Particle_parents.index")  #parent particles 
       #        .Alias("Particle1", "_Particle_daughters.index")  #daugther particles 
        #       .Alias("MCRecoAssociations0", "_MCRecoAssociations_rec.index") #points to the ReconstructedParticles collection
         #      .Alias("MCRecoAssociations1", "_MCRecoAssociations_sim.index") #points to the Particle collection
          #     .Alias("Track1", "_EFlowTrack_tracks.index")
           #    .Alias("TrackState", "_EFlowTrack_trackStates")
            #   .Alias("Muon0", "Muon_objIdx.index")

        #       .Define("ReconstructedParticles_filtered_pt", "FCCAnalyses::ReconstructedParticle::sel_pt(1.0)(ReconstructedParticles)") #selecting only reconstructed particles with pt > 1.5 GeV because otherwise muons don't reach the muon chambers
         #      .Define("ReconstructedParticles_filtered_p", "FCCAnalyses::ReconstructedParticle::sel_p(1.0)(ReconstructedParticles)") #selecting only reconstructed particles with p > 3 GeV because otherwise muons don't reach the muon chambers

               #get muons parameters
               .Define("muons", "FCCAnalyses::ReconstructedParticle::get(Muon0, ReconstructedParticles)")
               .Define("muons_p", "FCCAnalyses::ReconstructedParticle::get_p(muons)")
               .Define("muons_px", "FCCAnalyses::ReconstructedParticle::get_px(muons)")
               .Define("muons_py", "FCCAnalyses::ReconstructedParticle::get_py(muons)")
               .Define("muons_pz", "FCCAnalyses::ReconstructedParticle::get_pz(muons)")
               .Define("muon_eta", "FCCAnalyses::ReconstructedParticle::get_eta(muons)")
               .Define("muon_energy", "FCCAnalyses::ReconstructedParticle::get_e(muons)")
               .Define("muon_mass", "FCCAnalyses::ReconstructedParticle::get_mass(muons)")
               .Define("muon_phi", "FCCAnalyses::ReconstructedParticle::get_phi(muons)")
               .Define("muon_n", "FCCAnalyses::ReconstructedParticle::get_n(muons)")

               .Define("muon_TrackState", "FCCAnalyses::ReconstructedParticle2Track::getRP2TRK(muons, TrackState)")
               
               .Define("recopart_px", "FCCAnalyses::ReconstructedParticle::get_px(ReconstructedParticles)")
               .Define("recopart_py", "FCCAnalyses::ReconstructedParticle::get_py(ReconstructedParticles)")
               .Define("recopart_pz", "FCCAnalyses::ReconstructedParticle::get_pz(ReconstructedParticles)")
               .Define("recopart_charge", "FCCAnalyses::ReconstructedParticle::get_charge(ReconstructedParticles)")
               .Define("recopart_p", "FCCAnalyses::ReconstructedParticle::get_p(ReconstructedParticles)")
               .Define("recopart_eta", "FCCAnalyses::ReconstructedParticle::get_eta(ReconstructedParticles)")
               .Define("recopart_energy", "FCCAnalyses::ReconstructedParticle::get_e(ReconstructedParticles)")
               .Define("recopart_mass", "FCCAnalyses::ReconstructedParticle::get_mass(ReconstructedParticles)")
               .Define("recopart_phi", "FCCAnalyses::ReconstructedParticle::get_phi(ReconstructedParticles)")
               .Define("AllReco_isChargedPion", "recopart_mass > 0.139 and recopart_mass < 0.140")
               .Define("AllReco_ChargedPion_indices", "ROOT::VecOps::Nonzero(AllReco_isChargedPion)")
               .Define("chargedpion_particles", "FCCAnalyses::ReconstructedParticle::get(AllReco_ChargedPion_indices, ReconstructedParticles)")

               .Define("Pipm_TrackState", "FCCAnalyses::ReconstructedParticle2Track::getRP2TRK(chargedpion_particles, TrackState)")

              #.Define("TrackState_ptcut", "FCCAnalyses::ReconstructedParticle2Track::getRP2TRK(ReconstructedParticles_filtered_pt, TrackState)")
              #.Define("TrackState_pcut", "FCCAnalyses::ReconstructedParticle2Track::getRP2TRK(ReconstructedParticles_filtered_p, TrackState)")

               .Define("chargedpion_px",       "recopart_px[AllReco_isChargedPion]")
               .Define("chargedpion_py",       "recopart_py[AllReco_isChargedPion]")
               .Define("chargedpion_pz",       "recopart_pz[AllReco_isChargedPion]")
               .Define("chargedpion_p",        "recopart_p[AllReco_isChargedPion]")
               .Define("chargedpion_phi",      "recopart_phi[AllReco_isChargedPion]")
               .Define("chargedpion_eta",      "recopart_eta[AllReco_isChargedPion]")
               .Define("chargedpion_energy",   "recopart_energy[AllReco_isChargedPion]")
               .Define("chargedpion_mass",     "recopart_mass[AllReco_isChargedPion]")
               .Define("chargedpion_charge",   "recopart_charge[AllReco_isChargedPion]")

               #############################################################
               ##  these lines get truth info (generator info) of  kaons  ##
               ##  number of each kaon type per event, and their energy   ##
               #############################################################
               ##  plot1_multi.py will use the following branches         ##
               #############################################################
               .Define("genKS",   "FCCAnalyses::MCParticle::sel_pdgID(310, true)(Particle)")
               .Define("genmuons", "FCCAnalyses::MCParticle::sel_pdgID(abs(13), true)(Particle)")
               .Define("genpions", "FCCAnalyses::MCParticle::sel_pdgID(abs(211), true)(Particle)")
               #.Define("genmuons_from_pion", "ROOT::VecOps::RVec<int> result; for (size_t i = 0; i < Particle.size(); ++i) { if (abs(Particle.at(i).PDG) == 13) { int flag = 0; for (size_t j = Particle.at(i).parents_begin; j < Particle.at(i).parents_end; ++j) { int parent_index = Particle0.at(j); if (parent_index >= 0 && parent_index < Particle.size() && abs(Particle.at(parent_index).PDG) == 211) { flag = 1; break; } } result.push_back(flag); } } return result;")

               ## some variables of generated muons, kaons and pions
               .Define("KS_found", "FCCAnalyses::ZHfunctions::get_KS_decays(Particle)")
               .Define("n_KS_decays", "KS_found.size()")

               .Define("n_genKSs",      "FCCAnalyses::MCParticle::get_n(genKS)")
               .Define("n_genKSs_status", "FCCAnalyses::MCParticle::get_genStatus(genKS)")
               .Define("n_genmuons",    "FCCAnalyses::MCParticle::get_n(genmuons)")
               .Define("n_genpions",    "FCCAnalyses::MCParticle::get_n(genpions)")

               .Define("genKS_energy",     "FCCAnalyses::MCParticle::get_e(genKS)")
               .Define("genmuons_energy", "FCCAnalyses::MCParticle::get_e(genmuons)")
               .Define("genmuons_px", "FCCAnalyses::MCParticle::get_px(genmuons)")
               .Define("genmuons_py", "FCCAnalyses::MCParticle::get_py(genmuons)")
               .Define("genmuons_pz", "FCCAnalyses::MCParticle::get_pz(genmuons)")
               .Define("genmuons_p", "FCCAnalyses::MCParticle::get_p(genmuons)")
               .Define("genmuons_eta", "FCCAnalyses::MCParticle::get_eta(genmuons)")
               .Define("genmuons_phi", "FCCAnalyses::MCParticle::get_phi(genmuons)")
               .Define("genpions_energy", "FCCAnalyses::MCParticle::get_e(genpions)")
               .Define("genpions_px", "FCCAnalyses::MCParticle::get_px(genpions)")
               .Define("genpions_py", "FCCAnalyses::MCParticle::get_py(genpions)")
               .Define("genpions_pz", "FCCAnalyses::MCParticle::get_pz(genpions)")
               .Define("genpions_p", "FCCAnalyses::MCParticle::get_p(genpions)")
               .Define("genpions_eta", "FCCAnalyses::MCParticle::get_eta(genpions)")
               .Define("genpions_phi", "FCCAnalyses::MCParticle::get_phi(genpions)")
               ## DEBUG: Check pion decay information
               .Define("debug_n_particles", "int(Particle.size())")
               .Define("debug_muon_parents", "ROOT::VecOps::RVec<int> parents; for(size_t i=0; i<genmuons.size(); ++i) { int mu_idx = -1; for(size_t j=0; j<Particle.size(); ++j) { if(abs(Particle.at(j).PDG)==13 && Particle.at(j).momentum.x == genmuons.at(i).momentum.x) { mu_idx = j; break; } } if(mu_idx >= 0) { for(size_t k=Particle.at(mu_idx).parents_begin; k<Particle.at(mu_idx).parents_end; ++k) { if(k < Particle0.size()) parents.push_back(Particle.at(Particle0.at(k)).PDG); } } } return parents;")

               #some lines to check how many cases we have in which a muon and a pion are very close to each other
               .Define("genpions_vertex_x", "FCCAnalyses::MCParticle::get_vertex_x(genpions)")
               .Define("genpions_vertex_y", "FCCAnalyses::MCParticle::get_vertex_y(genpions)")
               .Define("genpions_vertex_z", "FCCAnalyses::MCParticle::get_vertex_z(genpions)")
               .Define("genmuons_vertex_x", "FCCAnalyses::MCParticle::get_vertex_x(genmuons)")
               .Define("genmuons_vertex_y", "FCCAnalyses::MCParticle::get_vertex_y(genmuons)")
               .Define("genmuons_vertex_z", "FCCAnalyses::MCParticle::get_vertex_z(genmuons)")
               .Define("genmuons_vertex_r", "sqrt(genmuons_vertex_x*genmuons_vertex_x + genmuons_vertex_y*genmuons_vertex_y)")
               .Define("genmuons_vertex_d", "sqrt(genmuons_vertex_x*genmuons_vertex_x + genmuons_vertex_y*genmuons_vertex_y + genmuons_vertex_z*genmuons_vertex_z)")
               .Define("pionsmuons_close", "FCCAnalyses::ZHfunctions::get_close_pairs(genpions_vertex_x, genpions_vertex_y, genpions_vertex_z, genmuons_vertex_x, genmuons_vertex_y, genmuons_vertex_z, 1)")
               .Define("pionsmuons_close_d", "pionsmuons_close.d")
               .Define("pionsmuons_close_pix", "pionsmuons_close.pix")
               .Define("pionsmuons_close_piy", "pionsmuons_close.piy")
               .Define("pionsmuons_close_piz", "pionsmuons_close.piz")
               .Define("pionsmuons_close_mux", "pionsmuons_close.mux")
               .Define("pionsmuons_close_muy", "pionsmuons_close.muy")
               .Define("pionsmuons_close_muz", "pionsmuons_close.muz")
               .Define("pionsmuons_close_pid", "sqrt(pionsmuons_close_pix*pionsmuons_close_pix + pionsmuons_close_piy*pionsmuons_close_piy + pionsmuons_close_piz*pionsmuons_close_piz)")
               .Define("pionsmuons_close_mud", "sqrt(pionsmuons_close_mux*pionsmuons_close_mux + pionsmuons_close_muy*pionsmuons_close_muy + pionsmuons_close_muz*pionsmuons_close_muz)")
               .Define("pionsmuons_close_pid_displaced", "pionsmuons_close_pid[pionsmuons_close_d > 1]")
               .Define("pionsmuons_close_mud_displaced", "pionsmuons_close_mud[pionsmuons_close_d > 1]")

               .Define("Ks_daughter_status", "FCCAnalyses::ZHfunctions::get_Ks_daughter_status(Particle, Particle0)")
               .Define("Ks_daughter_PDG", "FCCAnalyses::ZHfunctions::get_Ks_daughter_pdg(Particle, Particle0)")

               .Define("genKS_status", "FCCAnalyses::MCParticle::get_genStatus(genKS)")

               .Define("genmuons_status", "FCCAnalyses::MCParticle::get_genStatus(genmuons)")

               .Define("genpions_status", "FCCAnalyses::MCParticle::get_genStatus(genpions)")

               .Define("pions_from_Ks_idx", "FCCAnalyses::ZHfunctions::get_Ks_pions_indices(Particle, Particle1)")
               .Define("pions_from_Ks", "FCCAnalyses::ZHfunctions::sel_indices()(Particle, pions_from_Ks_idx)")
               .Define("pions_from_Ks_p", "FCCAnalyses::MCParticle::get_p(pions_from_Ks)")
               .Define("pions_from_Ks_pt", "FCCAnalyses::MCParticle::get_pt(pions_from_Ks)")

               #########################################################
               ##  these lines get true vertex properties             ##
               ##  "MCVertexObject" is a vector of MC vertex objects ##
               ##  the rest are vector of float (or int) variables    ##
               #########################################################
               
               
               .Define("MCVertexObject", "myUtils::get_MCVertexObject(Particle, Particle0)")
               .Define("MC_Vertex_x",    "myUtils::get_MCVertex_x(MCVertexObject)")
               .Define("MC_Vertex_y",    "myUtils::get_MCVertex_y(MCVertexObject)")
               .Define("MC_Vertex_r",     "sqrt(MC_Vertex_x*MC_Vertex_x+MC_Vertex_y*MC_Vertex_y)")
               .Define("MC_Vertex_z",    "myUtils::get_MCVertex_z(MCVertexObject)")
               .Define("MC_Vertex_d",    "sqrt(MC_Vertex_x*MC_Vertex_x+MC_Vertex_y*MC_Vertex_y+MC_Vertex_z*MC_Vertex_z)")
               .Define("MC_Vertex_ind",  "myUtils::get_MCindMCVertex(MCVertexObject)")
               .Define("MC_Vertex_ntrk", "myUtils::get_NTracksMCVertex(MCVertexObject)")
               .Define("MC_Vertex_n",    "int(MC_Vertex_x.size())")
               .Define("MC_Vertex_PDG",  "myUtils::get_MCpdgMCVertex(MCVertexObject, Particle)")
               .Define("MC_Vertex_PDGmother",  "myUtils::get_MCpdgMotherMCVertex(MCVertexObject, Particle)")
               .Define("MC_Vertex_PDGgmother", "myUtils::get_MCpdgGMotherMCVertex(MCVertexObject, Particle)")
               .Define("MC_Vertex_mass",  "FCCAnalyses::ZHfunctions::get_MC_Vertex_mass(MCVertexObject, Particle)")
               .Define("MC_Vertex_eta",   "FCCAnalyses::ZHfunctions::get_MC_Vertex_eta(MCVertexObject, Particle)")
               .Define("MC_Vertex_phi",   "FCCAnalyses::ZHfunctions::get_MC_Vertex_phi(MCVertexObject, Particle)")
               .Define("MC_Vertex_p",     "FCCAnalyses::ZHfunctions::get_MC_Vertex_p(MCVertexObject, Particle)")
               .Define("MC_Vertex_pt",     "FCCAnalyses::ZHfunctions::get_MC_Vertex_pt(MCVertexObject, Particle)")

               .Define("MC_Vertex_isKSpipi", "FCCAnalyses::ZHfunctions::get_MC_Vertex_isKSpipi(MC_Vertex_PDGmother, MC_Vertex_PDG)")
               .Define("MC_Vertex_isPimunu", "FCCAnalyses::ZHfunctions::get_MC_Vertex_isPimunu(MC_Vertex_PDGmother, MC_Vertex_PDG)")

               .Define("MC_Vertex_KSflag", "1.0*(MC_Vertex_isKSpipi >= 10)")
               .Define("MC_Vertex_PimunuFlag", "1.0*(MC_Vertex_isPimunu == 12)")
               .Define("n_total_vertices", "int(MC_Vertex_PDGmother.size())")
               .Define("genPimunu_Vertex_mass", "MC_Vertex_mass[MC_Vertex_PimunuFlag>0]")
               .Define("genPimunu_Vertex_r", "MC_Vertex_r[MC_Vertex_PimunuFlag>0]")
               .Define("genPimunu_Vertex_z", "MC_Vertex_z[MC_Vertex_PimunuFlag>0]")   
               .Define("genPimunu_Vertex_d", "MC_Vertex_d[MC_Vertex_PimunuFlag>0]")
               .Define("genPimunu_Vertex_d_displaced", "genPimunu_Vertex_d[genPimunu_Vertex_d>0.5]" )
               .Define("genPimunu_Vertex_eta", "MC_Vertex_eta[MC_Vertex_PimunuFlag>0]")
               .Define("genPimunu_Vertex_phi", "MC_Vertex_phi[MC_Vertex_PimunuFlag>0]")
               .Define("genPimunu_Vertex_p", "MC_Vertex_p[MC_Vertex_PimunuFlag>0]")
               .Define("genPimunu_Vertex_pt", "MC_Vertex_pt[MC_Vertex_PimunuFlag>0]")
               .Define("genKS_Vertex_mass", "MC_Vertex_mass[MC_Vertex_KSflag>0]")
               .Define("genKS_Vertex_r", "MC_Vertex_r[MC_Vertex_KSflag>0]")
               .Define("genKS_Vertex_z", "MC_Vertex_z[MC_Vertex_KSflag>0]")
               .Define("genKS_Vertex_d", "MC_Vertex_d[MC_Vertex_KSflag>0]")
               .Define("genKS_Vertex_eta", "MC_Vertex_eta[MC_Vertex_KSflag>0]")
               .Define("genKS_Vertex_phi", "MC_Vertex_phi[MC_Vertex_KSflag>0]")
               .Define("genKS_Vertex_p", "MC_Vertex_p[MC_Vertex_KSflag>0]")
               .Define("genKS_Vertex_pt", "MC_Vertex_pt[MC_Vertex_KSflag>0]")



               .Define("Kspimumu_candidates", "FCCAnalyses::ZHfunctions::get_PiMuFromKs(Particle, Particle0, Particle1)")
               .Define("Kspimumu_muon_indices", "Kspimumu_candidates.mu_idx")
               .Define("Kspimumu_pion_indices", "Kspimumu_candidates.pi_idx")
               .Define("Kspimumu_ks_indices", "Kspimumu_candidates.ks_idx")

               .Define("Kspimumu_muon_vertex_x", "FCCAnalyses::ZHfunctions::get_MC_vertex_x_from_indices(Particle, Kspimumu_muon_indices)")
               .Define("Kspimumu_muon_vertex_y", "FCCAnalyses::ZHfunctions::get_MC_vertex_y_from_indices(Particle, Kspimumu_muon_indices)")
               .Define("Kspimumu_muon_vertex_z", "FCCAnalyses::ZHfunctions::get_MC_vertex_z_from_indices(Particle, Kspimumu_muon_indices)")
               .Define("Kspimumu_muon_vertex_d", "sqrt(Kspimumu_muon_vertex_x*Kspimumu_muon_vertex_x+Kspimumu_muon_vertex_y*Kspimumu_muon_vertex_y+Kspimumu_muon_vertex_z*Kspimumu_muon_vertex_z)")
               .Define("Kspimumu_muon_vertex_r", "sqrt(Kspimumu_muon_vertex_x*Kspimumu_muon_vertex_x+Kspimumu_muon_vertex_y*Kspimumu_muon_vertex_y)")

               .Define("Kspimumu_muon_pt", "FCCAnalyses::ZHfunctions::get_MC_pt_from_indices(Particle, Kspimumu_muon_indices)")
               .Define("Kspimumu_pion_pt", "FCCAnalyses::ZHfunctions::get_MC_pt_from_indices(Particle, Kspimumu_pion_indices)")

               .Define("Kspimumu_opening_angle", "FCCAnalyses::ZHfunctions::get_pi_mu_opening_angle(Particle, Kspimumu_muon_indices, Kspimumu_pion_indices)")
               .Define("Kspimumu_pion_p", "FCCAnalyses::ZHfunctions::get_MC_p_from_indices(Particle, Kspimumu_pion_indices)")
               .Define("Kspimumu_muon_p", "FCCAnalyses::ZHfunctions::get_MC_p_from_indices(Particle, Kspimumu_muon_indices)")


                ##############################################
                # These lines are to define all the 
                # secondary vertices from reconstructed data.
                ##############################################

               .Define("Beamspot_x", "0.0")    # um
               .Define("Beamspot_y", "0.0")    # um
               .Define("Beamspot_z", "0.0")    # um
               .Define("sig_x", "5.96")        # um
               .Define("sig_y", "0.0238")      # um (23.8 nm)
               .Define("sig_z", "397.")        # um


               #.Define("SelectedTracks", "FCCAnalyses::ZHfunctions::getSelectedTracks(TrackState, TrackData, 100, 100)") #this is done to select some tracks that are
               #close to the PV so that we get the primary tracks faster.
          #     .Define("PrimaryTracks", "FCCAnalyses::ZHfunctions::getPrimaryTracks(TrackState, 3., Beamspot_x, Beamspot_y, Beamspot_z)")
               .Define("PrimaryTracks", "FCCAnalyses::VertexFitterSimple::get_PrimaryTracks(TrackState, true, sig_x, sig_y, sig_z, Beamspot_x, Beamspot_y, Beamspot_z)")
               .Define("PrimaryVertexObject", "FCCAnalyses::ZHfunctions::fitRecoPrimaryVertex(PrimaryTracks, Beamspot_x, Beamspot_y, Beamspot_z, sig_x, sig_y, sig_z)")
               .Define("PrimaryVertex", "FCCAnalyses::VertexingUtils::get_VertexData(PrimaryVertexObject)")
               .Define("SecondaryTracks", "FCCAnalyses::ZHfunctions::getSecondaryTracks(TrackState, PrimaryTracks)")
               .Define("EventV0s", "FCCAnalyses::ZHfunctions::get_V0s(SecondaryTracks, TrackState, ReconstructedParticles, PrimaryVertexObject, false, 10)")
               .Define("inv_mass", "EventV0s.invM")
               .Define("V0_pdg", "EventV0s.pdgAbs")
               .Define("inv_mass_pipi", "inv_mass[V0_pdg == 310]")
               .Define("inv_mass_mumu", "inv_mass[V0_pdg == 1313]")
               .Define("VertexObject", "EventV0s.vtx")
               .Define("VertexObject_mumu", "VertexObject[V0_pdg == 1313]")
               .Define("VertexObject_pipi", "VertexObject[V0_pdg == 310]")
               .Define("Vertex_x", "myUtils::get_Vertex_x(VertexObject)")
               .Define("Vertex_y", "myUtils::get_Vertex_y(VertexObject)")
               .Define("Vertex_z", "myUtils::get_Vertex_z(VertexObject)")
               .Define("Vertex_r", "sqrt(Vertex_x*Vertex_x+Vertex_y*Vertex_y)")
               .Define("recoKS_Vertex_r_mumu", "Vertex_r[V0_pdg == 1313]")
               .Define("inv_mass_before_VerDet", "inv_mass[Vertex_r<13.7]")  #mass of the V0 with two pions that decayed in flight to muons
               .Define("inv_mass_within_VerDet", "inv_mass[Vertex_r>12.7 && Vertex_r<34]")
               .Define("inv_mass_within_DC", "inv_mass[Vertex_r>34]")
               .Define("inv_mass_displaced", "inv_mass[Vertex_r>1]")
               .Define("Vertex_eta", "FCCAnalyses::ZHfunctions::get_Vertex_eta(VertexObject, ReconstructedParticles)")
               .Define("chi2", "FCCAnalyses::VertexingUtils::get_chi2_SV(VertexObject)")
               .Define("V0_tracks_momentum_at_vertex", "FCCAnalyses::ZHfunctions::get_momentum_updated_at_vertex(VertexObject)")
               .Define("V0_tracks_px_at_vertex", "FCCAnalyses::ZHfunctions::get_px_updated_at_vertex(VertexObject)")
               .Define("V0_tracks_py_at_vertex", "FCCAnalyses::ZHfunctions::get_py_updated_at_vertex(VertexObject)")
               .Define("V0_tracks_pz_at_vertex", "FCCAnalyses::ZHfunctions::get_pz_updated_at_vertex(VertexObject)")
               .Define("V0_tracks_TOF_ps", "FCCAnalyses::ZHfunctions::get_V0_tracks_TOF(VertexObject, ReconstructedParticles, EFlowTrack, TrackerHits)")
               .Define("V0_tracks_length_mm", "FCCAnalyses::ZHfunctions::get_V0_tracks_length(VertexObject, ReconstructedParticles, EFlowTrack_L)")
               .Define("V0_tracks_p_raw", "FCCAnalyses::ZHfunctions::get_V0_tracks_p_raw(VertexObject, ReconstructedParticles)")
               .Define("V0_tracks_mass_TOF", "FCCAnalyses::ZHfunctions::get_V0_tracks_mass_TOF(V0_tracks_p_raw, V0_tracks_length_mm, V0_tracks_TOF_ps)")
               .Define("V0_ks_t0_ps", "FCCAnalyses::ZHfunctions::get_V0_ks_t0(Vertex_x, Vertex_y, Vertex_z, V0_tracks_px_at_vertex, V0_tracks_py_at_vertex, V0_tracks_pz_at_vertex)")
               .Define("V0_tracks_mass_TOF_corrected", "FCCAnalyses::ZHfunctions::get_V0_tracks_mass_TOF_corrected(V0_tracks_p_raw, V0_tracks_length_mm, V0_tracks_TOF_ps, V0_ks_t0_ps)")
               
      #         .Define("V0_tracks_px_at_vertex", "FCCAnalyses::ZHfunctions::get_momentum_updated_at_vertex_x(VertexObject)")
       #        .Define("V0_tracks_py_at_vertex", "FCCAnalyses::ZHfunctions::get_momentum_updated_at_vertex_y(VertexObject)")
        #       .Define("V0_tracks_pz_at_vertex", "FCCAnalyses::ZHfunctions::get_momentum_updated_at_vertex_z(VertexObject)")
               .Define("KS_mass_mumuhyp_all", "FCCAnalyses::ZHfunctions::get_KS_mass_mumu(VertexObject)") #all V0s (pions and pions that decayed in flight) with muon mass hyp
         #      .Define("KS_mass_mumuhyp_all_before_VerDet", "KS_mass_mumuhyp_all[Vertex_r<13.7]")
         #      .Define("KS_mass_mumuhyp_all_within_VerDet", "KS_mass_mumuhyp_all[Vertex_r>12.7 && Vertex_r<34]")
         #      .Define("KS_mass_mumuhyp_all_within_DC", "KS_mass_mumuhyp_all[Vertex_r>34]")
          #     .Define("KS_mass_mumuhyp_all_displaced", "KS_mass_mumuhyp_all[Vertex_r>1000]")
           #    .Define("KS_mass_mumuhyp_pipi", "FCCAnalyses::ZHfunctions::get_KS_mass_mumu(VertexObject_pipi)") #muon mass hyp only for pions that didnt decay in flight to muons
           #    .Define("Vertex_x_pipi", "myUtils::get_Vertex_x(VertexObject_pipi)")
           #    .Define("Vertex_y_pipi", "myUtils::get_Vertex_y(VertexObject_pipi)")
           #    .Define("Vertex_z_pipi", "myUtils::get_Vertex_z(VertexObject_pipi)")
           #    .Define("Vertex_r_pipi", "sqrt(Vertex_x_pipi*Vertex_x_pipi+Vertex_y_pipi*Vertex_y_pipi)")
           #    .Define("KS_mass_mumuhyp_pipi_before_VerDet", "KS_mass_mumuhyp_pipi[Vertex_r_pipi<13.7]")
           #    .Define("KS_mass_mumuhyp_pipi_within_VerDet", "KS_mass_mumuhyp_pipi[Vertex_r_pipi>12.7 && Vertex_r_pipi<34]")
            #   .Define("KS_mass_mumuhyp_pipi_within_DC", "KS_mass_mumuhyp_pipi[Vertex_r_pipi>34]")
            #   .Define("KS_mass_mumuhyp_pipi_displaced", "KS_mass_mumuhyp_pipi[Vertex_r_pipi>1000]")
               .Define("n_vertices", "EventV0s.vtx.size()")

               .Define("VtxMatchResult", "FCCAnalyses::ZHfunctions::match_MC2Reco_V0s(VertexObject, MCRecoAssociations0, MCRecoAssociations1, MCVertexObject, ReconstructedParticles, Particle, Particle0, MC_Vertex_x, MC_Vertex_y, MC_Vertex_z)")
               .Define("V0_dist_to_MC", "VtxMatchResult.distance")
               .Define("V0_matched_MC_idx", "VtxMatchResult.mc_vtx_idx")
               .Define("V0_matched_MC_pdg", "VtxMatchResult.mc_daughter_pdg")
               .Define("V0_matched_MC_pdgmother", "VtxMatchResult.mc_mother_pdg")
               .Define("V0_matched_MC_x", "VtxMatchResult.mc_vtx_x")
               .Define("V0_matched_MC_y", "VtxMatchResult.mc_vtx_y")
               .Define("V0_matched_MC_z", "VtxMatchResult.mc_vtx_z")

                .Define("V0_matched_daughters_idx", "VtxMatchResult.mc_daughter_indices")
                .Define("V0_matched_mother_idx", "VtxMatchResult.mc_mother_indices")
                .Define("V0_matched_gmother_idx", "VtxMatchResult.mc_gmother_indices")


                .Define("V0_is_pure_match", """
                        ROOT::VecOps::RVec<int> res;
                        size_t idx = 0;
                        for (const auto& mothers : V0_matched_mother_idx) {
                            if (mothers.size() < 2) {
                                res.push_back(0);
                                idx++;
                                continue;
                            }
                            
                            // Filter fake -1 values
                            bool all_valid = true;
                            for (int m : mothers) if (m == -1) { all_valid = false; break; }
                            
                            if (!all_valid) {
                                res.push_back(0);
                            } else {
                                // Check if the indices are the same
                                bool all_same = true;
                                for (size_t i = 1; i < mothers.size(); ++i) {
                                    if (mothers[i] != mothers[0]) {
                                        all_same = false;
                                        break;
                                    }
                                }
                                
                                // Also check that both mother PDGs are Ks (310)
                                if (all_same && idx < V0_matched_MC_pdgmother.size()) {
                                    const auto& mother_pdgs = V0_matched_MC_pdgmother[idx];
                                    if (mother_pdgs.size() >= 2) {
                                        // Both mothers should be Ks (310)
                                        bool both_ks = (mother_pdgs[0] == 310 && mother_pdgs[1] == 310);
                                        res.push_back(both_ks ? 1 : 0);
                                    } else {
                                        res.push_back(0);
                                    }
                                } else {
                                    res.push_back(0);
                                }
                            }
                            idx++;
                        }
                        return res;
                    """)

                .Define("V0_matched_daughters_p", "FCCAnalyses::ZHfunctions::get_nested_daughter_p()(Particle, V0_matched_daughters_idx)")
                .Define("V0_matched_mothers_p", "FCCAnalyses::ZHfunctions::get_nested_mother_p()(Particle, V0_matched_mother_idx)")
                .Define("V0_matched_gmother_pdg", "FCCAnalyses::ZHfunctions::get_nested_gmother_pdg()(Particle, Particle0, V0_matched_mother_idx)")
                .Define("V0_ks_mother_idx",
                        "FCCAnalyses::ZHfunctions::extract_V0_ks_idx()(V0_matched_mother_idx, V0_matched_MC_pdgmother)")
                .Define("V0_pimu_opening_angle",
                        "FCCAnalyses::ZHfunctions::get_V0_pimu_opening_angle(V0_ks_mother_idx, Kspimumu_candidates, Particle)")

        )
        return df2

    #__________________________________________________________
    #Mandatory: output function, please make sure you return the branchlist as a python list
    def output():
        # Stripped branch list for the large-scale background run: only what
        # background_estimation.py needs to classify/count pipi misID and pimu DIF.
        branchList = [
                 "n_KS_decays",
                 "KS_mass_mumuhyp_all",
                 "chi2",
                 "genPimunu_Vertex_mass",
                 "V0_matched_MC_pdg", "V0_matched_MC_pdgmother", "V0_matched_gmother_pdg",
                 "V0_matched_mother_idx", "V0_matched_gmother_idx",
                 "V0_is_pure_match",
                 "V0_matched_MC_x", "V0_matched_MC_y", "V0_matched_MC_z",
                 "Vertex_x", "Vertex_y", "Vertex_z",
                 "V0_tracks_momentum_at_vertex",
                 "V0_tracks_px_at_vertex", "V0_tracks_py_at_vertex", "V0_tracks_pz_at_vertex",
                 "V0_tracks_TOF_ps", "V0_tracks_length_mm", "V0_tracks_p_raw", "V0_tracks_mass_TOF",
                 "V0_ks_t0_ps", "V0_tracks_mass_TOF_corrected",
                 "V0_pimu_opening_angle",
                ]
        return branchList