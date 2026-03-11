#Mandatory: List of processes
processList = {
             'p8_ee_Zss_ecm91':{'chunks':1, 'fraction':0.0000001},
#             'p8_ee_Zcc_ecm91':{'chunks':1, 'fraction':0.0000001},
#             'p8_ee_Zbb_ecm91':{'chunks':1, 'fraction':0.0000001},
#             'p8_ee_Zud_ecm91':{'chunks':1, 'fraction':0.0000001},

            }

#Mandatory: Production tag when running over EDM4Hep centrally produced events, this points to the yaml files for getting sample statistics
prodTag     = "FCCee/winter2023/IDEA/"

#Optional: output directory, default is local running directory
outputDir   = "outputs/strange/analysis_stage1_1118/"

#EOS output directory for batch jobs
outputDirEos = "/eos/experiment/fcc/ee/analyses/case-studies/strange/displaced_vertex/flatNtuples/winter2023"

#Optional

nCPUS       = 8
runBatch    = False
batchQueue = "workday"
compGroup = "group_u_FCC.local_gen"

includePaths = ["functions.h"]#, "functions_mireya.h"]

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

               #############################################################
               ##  these lines get truth info (generator info) of  kaons  ##
               ##  number of each kaon type per event, and their energy   ## 
               #############################################################
               ##  plot1_multi.py will use the following branches         ##
               #############################################################

               .Define("genK0",   "FCCAnalyses::MCParticle::sel_pdgID(311, true)(Particle)")
               .Define("genKL",   "FCCAnalyses::MCParticle::sel_pdgID(130, true)(Particle)")
               .Define("genKS",   "FCCAnalyses::MCParticle::sel_pdgID(310, true)(Particle)")
               .Define("genKpos", "FCCAnalyses::MCParticle::sel_pdgID(321, false)(Particle)")
               .Define("genKneg", "FCCAnalyses::MCParticle::sel_pdgID(-321, false)(Particle)")

               ## number of kaons in the event

               .Define("n_genK0s",      "FCCAnalyses::MCParticle::get_n(genK0)")
               .Define("n_genKLs",      "FCCAnalyses::MCParticle::get_n(genKL)")
               .Define("n_genKSs",      "FCCAnalyses::MCParticle::get_n(genKS)")
               .Define("n_genKposs",    "FCCAnalyses::MCParticle::get_n(genKpos)")
               .Define("n_genKnegs",    "FCCAnalyses::MCParticle::get_n(genKneg)")

               #pions
               .Define("genPipm",       "FCCAnalyses::MCParticle::sel_pdgID(211, true)(Particle)")
               .Define("genPipms_n",    "FCCAnalyses::MCParticle::get_n(genPipm)")
               .Define("genPipm_p", "FCCAnalyses::MCParticle::get_p(genPipm)")
               .Define("genPipm_indices", "ROOT::VecOps::RVec<int> idx; for (size_t i=0; i < Particle.size(); ++i){ if (abs(Particle[i].PDG)==211) idx.push_back(i); }return idx;")

               #get the position of the pions at the decay vertex
               .Define("genPipm_vertex_x", "FCCAnalyses::MCParticle::get_vertex_x(genPipm)")
               .Define("genPipm_vertex_y", "FCCAnalyses::MCParticle::get_vertex_y(genPipm)")
               .Define("genPipm_vertex_z", "FCCAnalyses::MCParticle::get_vertex_z(genPipm)")
               .Define("genPipm_vertex_r", "sqrt(genPipm_vertex_x*genPipm_vertex_x + genPipm_vertex_y*genPipm_vertex_y)")

               #get the generated pions that come from the decay of the kaon
                #I do not know why this doesnt work so I will comment it :(
                #.Define("KS_to_Pipm_indices_obj", "FCCAnalyses::MCParticle::get_indices(310, {211, -211}, true, false, false, false)")
                #.Define("genPipm_fromKS_indices", "KS_to_Pipm_indices_obj(Particle, genPipm_indices)")
                #.Define("n_genPipm_fromKS", "int(genPipm_fromKS_indices.size())")
                #.Define("genPipm_fromKS_p", "genPipm_p[genPipm_fromKS_indices]")
                #.Define("genPipm_fromKS_vertex_r", "genPipm_vertex_r[genPipm_fromKS_indices]")
                #Trying a different approach to get the pions from KS decays.
                .Define("parents_Pipm_indices", "FCCAnalyses::MCParticle::get_parentid(genPipm_indices, Particle, Particle0)")
                .Define("pdg_allmc", "FCCAnalyses::McParticle::get_pdg(Particle)")
                .Define("parents_KS_Pipm_indices", "ROOT::VecOps::RVec<int> result; for (int idx : parents_Pipm_indices) { if (pdg_allmc[idx] == 310) result.push_back(idx); } return result;")

                .Define("genPipm_fromKS_indices", "ROOT::VecOps::RVec<int> result; for (int idx : parents_KS_Pipm_indices) {auto daughters = FCCAnalyses::MCParticle::get_list_of_particles_from_decay(idx, Particle, Particle1); for (int v : daughters) { result.push_back(v); } } return result;")
                .Define("genPipm_fromKS_vertex_r", "genPipm_vertex_r[gen_Pipm_fromKS_indices]")
                .Define("genPipm_fromKS_vertex_p", "genPipm_p[gen_Pipm_fromKS_indices]")

               #get those pions with reconstructed tracks
               #I am commenting these lines since they are giving dimensional error when looking at MC_recotracks_indices and genPipm_indices)
               #.Define("MC_recotracks_indices", "FCCAnalyses::ZHfunctions::get_MCindTRK(ReconstructedParticles,EFlowTrack_1,MCRecoAssociations0,MCRecoAssociations1)")
               #.Define("genPipm_recotracks_mask", "ROOT::VecOps::RVec<int> mask(genPipm_indices.size(),0);for (size_t i=0;i<genPipm_indices.size();++i){ int idx = genPipm_indices[i];"
                #    "for (auto &trk : MC_recotracks_indices){if (trk == idx){ mask[i]=1; break; } } } return mask;")
               #.Define("genPipm_recotracks_indices", "ROOT::VecOps::RVec<int> result; for (size_t i=0; i < genPipm.size(); ++i) { int idx = MC_recotracks_indices[i]; if (idx>=0) result.push_back(idx); } return result;")
               #.Define("genPipm_recotracks_r", "genPipm_vertex_r[genPipm_recotracks_mask>0]")
               #.Define("genPipm_recotracks_p", "genPipm_p[genPipm_recotracks_mask>0]")

               ## energy of kaons

               .Define("genKS_energy",     "FCCAnalyses::MCParticle::get_e(genKS)")
               .Define("genKpos_energy",   "FCCAnalyses::MCParticle::get_e(genKpos)")
               .Define("genKneg_energy",   "FCCAnalyses::MCParticle::get_e(genKneg)")


               #########################################################
               ##  these lines get true vertex properties             ##
               ##  "MCVertexObject" is an vector of MC vertex objects ##
               ##  the rest are vector of float (or int) variables    ##
               #########################################################
               .Define("MCVertexObject", "myUtils::get_MCVertexObject(Particle, Particle0)")
               .Define("MC_Vertex_x",    "myUtils::get_MCVertex_x(MCVertexObject)")
               .Define("MC_Vertex_y",    "myUtils::get_MCVertex_y(MCVertexObject)")
               .Define("MC_Vertex_z",    "myUtils::get_MCVertex_z(MCVertexObject)")
               .Define("MC_Vertex_ind",  "myUtils::get_MCindMCVertex(MCVertexObject)")
               .Define("MC_Vertex_ntrk", "myUtils::get_NTracksMCVertex(MCVertexObject)")
               .Define("MC_Vertex_n",    "int(MC_Vertex_x.size())")
               .Define("MC_Vertex_PDG",  "myUtils::get_MCpdgMCVertex(MCVertexObject, Particle)")
               .Define("MC_Vertex_PDGmother",  "myUtils::get_MCpdgMotherMCVertex(MCVertexObject, Particle)")
               .Define("MC_Vertex_PDGgmother", "myUtils::get_MCpdgGMotherMCVertex(MCVertexObject, Particle)")
               .Define("MC_Vertex_mass",  "FCCAnalyses::ZHfunctions::get_MC_Vertex_mass(MCVertexObject, Particle)")
               .Define("MC_Vertex_p",     "FCCAnalyses::ZHfunctions::get_MC_Vertex_p(MCVertexObject, Particle)")
               .Define("MC_Vertex_pt",    "FCCAnalyses::ZHfunctions::get_MC_Vertex_pt(MCVertexObject, Particle)")

                #Now I will only select those pions associated with a reconstructed vertex.
                #I am commenting these lines since I am using variables that come from the part of the code that crashes.
               #.Define("genPipm_recotracks_fromKS", "ROOT::VecOps::RVec<int> result;"
                #"for (size_t i = 0; i < genPipm_recotracks_indices.size(); ++i) {int idx = genPipm_recotracks_indices[i];"
                 #   "bool fromKS = false;"
                  #  "for (size_t v = 0; v < MC_Vertex_PDGmother.size(); ++v) {"
                   #     "for (size_t d = 0; d < MC_Vertex_PDGmother[v].size(); ++d) {"
                    #        "if (MC_Vertex_PDGmother[v][d] == 310 && d == idx) {"
                     #           "fromKS = true; break;} }"
                      #  "if (fromKS) break; } result.push_back(fromKS ? 1 : 0);} return result;")
               #.Define("genPipm_KSrecotracks_r", "genPipm_vertex_r[genPipm_recotracks_fromKS > 0]")
               #.Define("genPipm_KSrecotracks_p", "genPipm_p[genPipm_recotracks_fromKS > 0]")
               #.Define("genPipm_KSrecotracks_n", "int(genPipm_KSrecotracks_r.size())")


                ########################################################
                ##           single out true Bs vertices              ##
                ##  The "MC_Vertex_isKSpipi" is a hacky way to check  ##
                ##  for each MC vertex if it is a K_S to pi pi decay  ##
                ##  A better way is to put this long line of for loop ##
                ##  into a function and call the function             ##
                ########################################################
                ##  plot2_decay.py will use the following branches    ##
                ########################################################


               #This was wrong: .Define("MC_Vertex_isKSpipi",   "ROOT::VecOps::RVec<int> result; for (size_t i=0; i < MC_Vertex_PDGmother.size(); ++i) {int isKS=0; for (size_t j=0; j < MC_Vertex_PDGmother[i].size(); ++j) {if (abs(MC_Vertex_PDGmother[i][j])==310) isKS+=10;} for (size_t j=0; j < MC_Vertex_PDG[i].size(); ++j) {if (abs(MC_Vertex_PDG[i][j])==211) isKS+=1;} result.push_back(isKS);} return result;")
               .Define("MC_Vertex_isKSpipi", "ROOT::VecOps::RVec<int> result; for (size_t i=0; i < MC_Vertex_PDGmother.size(); ++i) {int isKS=0; for (size_t j=0; j < MC_Vertex_PDGmother[i].size(); ++j) {if (abs(MC_Vertex_PDGmother[i][j])==310) {isKS+=10; for (size_t k=0; k < MC_Vertex_PDG[i].size(); ++k) {if (abs(MC_Vertex_PDG[i][k])==211) isKS+=1;} break;}} result.push_back(isKS);} return result;")

               .Define("MC_Vertex_KSflag", "1.0*(MC_Vertex_isKSpipi > 10 && MC_Vertex_isKSpipi % 10 !=0)")
               .Define("genKS_Vertex_x",   "MC_Vertex_x[MC_Vertex_KSflag>0]") #find Ks meson
               .Define("genKS_Vertex_y",   "MC_Vertex_y[MC_Vertex_KSflag>0]") #no implementation of abs() for vector
               .Define("genKS_Vertex_z",   "MC_Vertex_z[MC_Vertex_KSflag>0]")
               .Define("genKS_Vertex_n", "int(genKS_Vertex_x.size())")
               .Define("genKS_Vertex_p",   "MC_Vertex_p[MC_Vertex_KSflag>0]")
               .Define("genKS_Vertex_pt",  "MC_Vertex_pt[MC_Vertex_KSflag>0]")
               .Define("genKS_Vertex_r",   "sqrt(genKS_Vertex_x*genKS_Vertex_x + genKS_Vertex_y*genKS_Vertex_y)")
               .Define("genKS_Vertex_acceptance_r", "genKS_Vertex_r[(genKS_Vertex_r < 2000) && (abs(genKS_Vertex_z) < 2000)]")
               .Define("genKS_Vertex_acceptance_n", "int(genKS_Vertex_acceptance_r.size())")
               .Define("genKS_Vertex_d", "sqrt(genKS_Vertex_x*genKS_Vertex_x + genKS_Vertex_y*genKS_Vertex_y + genKS_Vertex_z*genKS_Vertex_z)")

               #############################################
               ##              Build Reco Vertex          ##
               #############################################
               .Define("VertexObject", "myUtils::get_VertexObject(MCVertexObject,ReconstructedParticles,EFlowTrack_1,MCRecoAssociations0,MCRecoAssociations1)")
               .Define("EVT_NVertex",   "float(VertexObject.size())")

               ####################################################
               ##    read RECO particles with PID at vertex      ##
               ####################################################
               .Define("RecoPartPID" ,"myUtils::PID(ReconstructedParticles, MCRecoAssociations0,MCRecoAssociations1,Particle)")
               .Define("RecoPartPIDAtVertex" ,"myUtils::get_RP_atVertex(RecoPartPID, VertexObject)")

               #############################################
               ##         Build vertex variables          ##
               #############################################

               ## these are the x,y,z positions of reconstructed vertices
               .Define("Vertex_x",        "myUtils::get_Vertex_x(VertexObject)")
               .Define("Vertex_y",        "myUtils::get_Vertex_y(VertexObject)")
               .Define("Vertex_z",        "myUtils::get_Vertex_z(VertexObject)")
               .Define("Vertex_xErr",     "myUtils::get_Vertex_xErr(VertexObject)")
               .Define("Vertex_yErr",     "myUtils::get_Vertex_yErr(VertexObject)")
               .Define("Vertex_zErr",     "myUtils::get_Vertex_zErr(VertexObject)")

               .Define("Vertex_chi2",     "myUtils::get_Vertex_chi2(VertexObject)")
               .Define("Vertex_mcind",    "myUtils::get_Vertex_indMC(VertexObject)")
               .Define("Vertex_ind",      "myUtils::get_Vertex_ind(VertexObject)")
               .Define("Vertex_isPV",     "myUtils::get_Vertex_isPV(VertexObject)")
               .Define("Vertex_ntrk",     "myUtils::get_Vertex_ntracks(VertexObject)")
               .Define("Vertex_n",        "int(Vertex_x.size())")
               .Define("Vertex_mass",     "myUtils::get_Vertex_mass(VertexObject,RecoPartPIDAtVertex)")
               .Define("Vertex_pt",       "FCCAnalyses::ZHfunctions::get_Vertex_pt(VertexObject,RecoPartPIDAtVertex)")
               .Define("Vertex_eta",      "FCCAnalyses::ZHfunctions::get_Vertex_eta(VertexObject,RecoPartPIDAtVertex)")

               #check for Kaons decaying to pions that have reconstructable tracks 
                .Define("genKS_tracks_indices",
                    "ROOT::VecOps::RVec<int> result;"
                    "auto trk_mc = FCCAnalyses::ZHfunctions::get_MCindTRK(ReconstructedParticles,EFlowTrack_1,MCRecoAssociations0,MCRecoAssociations1);"
                    "for (size_t v=0; v<MC_Vertex_PDGmother.size(); ++v) {bool isKS = false; for (size_t m=0; m<MC_Vertex_PDGmother[v].size(); ++m) {"
                            "if (abs(MC_Vertex_PDGmother[v][m]) == 310) { isKS = true; break; } } if (!isKS) continue; int nRecoPions = 0;"
                        "for (size_t d=0; d<MC_Vertex_PDG[v].size(); ++d) {if (abs(MC_Vertex_PDG[v][d])==211 && ROOT::VecOps::Any(trk_mc==MC_Vertex_ind[v][d])) nRecoPions++; }"
                        "if (nRecoPions==2) result.push_back(v);} return result;")

                .Define("genKS_tracks_n", "int(genKS_tracks_indices.size())")

                .Define("genKS_tracks_x", "ROOT::VecOps::RVec<float> result;for (auto i : genKS_tracks_indices) result.push_back(MC_Vertex_x[i]); return result;")
                .Define("genKS_tracks_y", "ROOT::VecOps::RVec<float> result; for (auto i : genKS_tracks_indices) result.push_back(MC_Vertex_y[i]); return result;")
                .Define("genKS_tracks_z",
                    "ROOT::VecOps::RVec<float> result;" "for (auto i : genKS_tracks_indices) result.push_back(MC_Vertex_z[i]);"
                    "return result;")
               .Define("MC_tracks_p", "ROOT::VecOps::RVec<float> result; for (auto i : genKS_tracks_indices) result.push_back(MC_Vertex_p[i]);"
                    "return result;")
               .Define("genKS_tracks_r", "sqrt(genKS_tracks_x*genKS_tracks_x + genKS_tracks_y*genKS_tracks_y)")
               .Define("genKS_tracks_d", "sqrt(genKS_tracks_x*genKS_tracks_x + genKS_tracks_y*genKS_tracks_y + genKS_tracks_z*genKS_tracks_z)")


               ######################################################################
               ## MATCHING RECO PIpm TO MC PIpm
               ######################################################################
               #getRP2MC_index creates a vector that maps the reconstructed particles and the MC particles: RP_MC_index[ ireco ] = imc
               .Define("RP_MC_index", "ReconstructedParticle2MC::getRP2MC_index(MCRecoAssociations0,MCRecoAssociations1,ReconstructedParticles)") 
               #now I need the indices of the reco particles that are pions
               .Define("Pip_reco_indices",  "myUtils::sel_PID(211)(ReconstructedParticles)")
               .Define("Pim_reco_indices", "myUtils::sel_PID(-211)(ReconstructedParticles)")
               .Define("Pipm_reco_indices", "ROOT::VecOps::Concatenate(Pip_reco_indices, Pim_reco_indices)")
               #now I can get the indices of the MC particles that are matched to the reconstructed pions
               .Define("Pipm_RP2MC_indices", "RP_MC_index[Pipm_reco_indices]")
               #momentum and production position of all reconstructed pions
                .Define("genPipm_RP2MC_p", "genPipm_p[Pipm_RP2MC_indices]")
                .Define("genPipm_RP2MC_vertex_r", "genPipm_vertex_r[Pipm_RP2MC_indices]")
                #Now I will select only those that come from a KS. I will do this at MC level, not reconstructed level... Maybe it's cheating.
               #.Define("genPipm_fromKSRP2MC_indices", "KS_to_Pipm_indices_obj(Particle, Pipm_RP2MC_indices)")
                #now I will select the p and r of the reconstructed pions that come from KS
               #.Define("genPipm_fromKSRP2MC_p", "genPipm_p[genPipm_fromKSRP2MC_indices]")
               #.Define("genPipm_fromKSRP2MC_vertex_r", "genPipm_vertex_r[genPipm_fromKSRP2MC_indices]")
                 #different approach (same as above)
                #already done: .Define("pdg_allmc", "FCCAnalyses::McParticle::get_pdg(Particle)")
                .Define("parents_recoPipm_indices", "FCCAnalyses::MCParticle::get_parentid(genPipm_RP2MC_vertex_r, Particle, Particle0)")
                .Define("parents_recoPipmfromKS_indices", "ROOT::VecOps::RVec<int> result; for (int idx : parents_recoPipm_indices) { "
                        "if (pdg_allmc[idx] == 310) result.push_back(idx); } return result;")

                .Define("genPipm_recofromKS_indices", "ROOT::VecOps::RVec<int> result; for (int idx : parents_recoPipmfromKS_indices) { 
                        "auto daughters = FCCAnalyses::MCParticle::get_list_of_particles_from_decay(idx, Particle, Particle1); "
                        "for (int v : daughters) { result.push_back(v); } } return result;")
                .Define("genPipm_recofromKS_vertex_r", "genPipm_vertex_r[gen_Pipm_fromKS_indices]")
                .Define("genPipm_fromKS_vertex_p", "genPipm_p[gen_Pipm_fromKS_indices]")
        

               ## these are the true x,y,z positions of the MC vertices matched to the reco vertices
               .Define("Vertex_MCx", "ROOT::VecOps::RVec<float> result; for (size_t i=0; i < Vertex_mcind.size(); ++i) result.push_back(MC_Vertex_x.at(Vertex_mcind[i])); return result;")
               .Define("Vertex_MCy",      "ROOT::VecOps::RVec<float> result; for (size_t i=0; i < Vertex_mcind.size(); ++i) result.push_back(MC_Vertex_y.at(Vertex_mcind[i])); return result;")
               .Define("Vertex_MCz",      "ROOT::VecOps::RVec<float> result; for (size_t i=0; i < Vertex_mcind.size(); ++i) result.push_back(MC_Vertex_z.at(Vertex_mcind[i])); return result;")
               .Define("Vertex_MCp",      "ROOT::VecOps::RVec<float> result; for (size_t i=0; i < Vertex_mcind.size(); ++i) result.push_back(MC_Vertex_p.at(Vertex_mcind[i])); return result;")
               .Define("Vertex_MCd",      "sqrt( (Vertex_MCx-Vertex_x)*(Vertex_MCx-Vertex_x) + (Vertex_MCy-Vertex_y)*(Vertex_MCy-Vertex_y) + (Vertex_MCz-Vertex_z)*(Vertex_MCz-Vertex_z)   )")
               .Define("Vertex_isMCKSpipi", "ROOT::VecOps::RVec<int> res; for (size_t i=0; i < Vertex_mcind.size(); ++i) { int isKS=0, idx=Vertex_mcind[i]; if (idx>=0) { for (size_t j=0; j < MC_Vertex_PDGmother[idx].size(); ++j) { if (abs(MC_Vertex_PDGmother[idx][j])==310) { isKS+=10; for (size_t k=0; k < MC_Vertex_PDG[idx].size(); ++k) { if (abs(MC_Vertex_PDG[idx][k])==211) isKS+=1; } break; } } } res.push_back(isKS); } return res;")
              # This was wrong:.Define("Vertex_isMCKSpipi",   "ROOT::VecOps::RVec<int> result; for (size_t i=0; i < Vertex_mcind.size(); ++i) {int isKS=0; for (size_t j=0; j < MC_Vertex_PDGmother[Vertex_mcind[i]].size(); ++j) {if (abs(MC_Vertex_PDGmother[i][j])==310) isKS+=10;} for (size_t l=0; l < MC_Vertex_PDG[Vertex_mcind[i]].size(); ++l) {if (abs(MC_Vertex_PDG[i][l])==211) isKS+=1;} result.push_back(isKS);} return result;")

               .Define("Vertex_KSmatch",  "1.0*(Vertex_isMCKSpipi > 10 && Vertex_isMCKSpipi % 10 !=0 && Vertex_chi2<10 && Vertex_MCd<10)") #I changed it from 2 mm
               #(too tight according to Michele) to 1 cm (maybe now it's too loose, need to check)

               ########################################################
               ##   check for mass reso for displaced vertices (KS)  ##
               ########################################################

               .Define("Vertex_rErr", "sqrt(Vertex_xErr*Vertex_xErr + Vertex_yErr*Vertex_yErr)")
               .Define("Vertex_dErr", "sqrt(Vertex_xErr*Vertex_xErr + Vertex_yErr*Vertex_yErr + Vertex_zErr*Vertex_zErr)")

               ## mass of all reconstructed vertices
               .Define("Vertex_r", "sqrt(Vertex_x*Vertex_x + Vertex_y*Vertex_y)")
               .Define("Vertex_d", "sqrt(Vertex_x*Vertex_x + Vertex_y*Vertex_y + Vertex_z*Vertex_z)")
               .Define("Vertex_mass_before_VerDet", "Vertex_mass[Vertex_chi2<10 && Vertex_r<13.7]") ## decays before first vertex detector layer
               .Define("Vertex_mass_within_VerDet", "Vertex_mass[Vertex_chi2<10 && Vertex_r>13.7 && Vertex_r<34]") ## decays within vertex detector volume
               .Define("Vertex_mass_within_DC",     "Vertex_mass[Vertex_chi2<10 && Vertex_r>35]")
               .Define("Vertex_mass_super_displaced", "Vertex_mass[Vertex_chi2<10 && Vertex_r>1000]")

               ## mass of reconstructed vertices that are matched as K_S decays
               .Define("recoKS_Vertex_r", "Vertex_r[Vertex_KSmatch>0]")
               .Define("recoKS_Vertex_z", "Vertex_z[Vertex_KSmatch>0]")
               .Define("recoKS_Vertex_mass", "Vertex_mass[Vertex_KSmatch>0]")
               .Define("recoKS_Vertex_mass_before_VerDet", "Vertex_mass[Vertex_KSmatch>0 && Vertex_r<13.7]")
               .Define("recoKS_Vertex_mass_within_VerDet", "Vertex_mass[Vertex_KSmatch>0 && Vertex_r>13.7 && Vertex_r<34]")
               .Define("recoKS_Vertex_mass_beyond_VerDet", "Vertex_mass[Vertex_KSmatch>0 && Vertex_r>35]")

               #the generated kaons that are associated with a reconstructed vertex and it's matched to the MC vertex
               .Define("genKS_recoVertex_n", "int(Sum(Vertex_KSmatch))")
               .Define("genKS_recoVertex_x", "ROOT::VecOps::RVec<float> result;"
                       "for (size_t i=0; i < Vertex_MCx.size(); ++i) {if (Vertex_KSmatch[i] == 1) result.push_back(Vertex_MCx[i]); }"
                       "return result;")
               .Define("genKS_recoVertex_y", "ROOT::VecOps::RVec<float> result;"
                       "for (size_t i=0; i < Vertex_MCy.size(); ++i) {if (Vertex_KSmatch[i] == 1) result.push_back(Vertex_MCy[i]); }"
                       "return result;")

               .Define("genKS_recoVertex_z", "ROOT::VecOps::RVec<float> result;"
                       "for (size_t i=0; i < Vertex_MCz.size(); ++i) {if (Vertex_KSmatch[i] == 1) result.push_back(Vertex_MCz[i]); }"
                       "return result;")
               .Define("genKS_recoVertex_p", "ROOT::VecOps::RVec<float> result;"
                        "for (size_t i=0; i < Vertex_MCp.size(); ++i) {if (Vertex_KSmatch[i] == 1) result.push_back(Vertex_MCp[i]); }"
                        "return result;")
               .Define("genKS_recoVertex_r", "sqrt(genKS_recoVertex_x*genKS_recoVertex_x + genKS_recoVertex_y*genKS_recoVertex_y)")
               .Define("genKS_recoVertex_d", "sqrt(genKS_recoVertex_x*genKS_recoVertex_x + genKS_recoVertex_y*genKS_recoVertex_y + genKS_recoVertex_z*genKS_recoVertex_z)")

                
        return df2)



    #__________________________________________________________
    #Mandatory: output function, please make sure you return the branchlist as a python list
    def output():
        branchList = [
                 "n_genK0s", "n_genKLs", "n_genKSs",
                 "n_genKposs", "n_genKnegs",
                 "genKS_energy", "genKpos_energy", "genKneg_energy",

                 "genPipm", "n_genPipms", "genPipm_p",
                 "genPipm_vertex_x", "genPipm_vertex_y", "genPipm_vertex_z", "genPipm_vertex_r",

                 "genKS_Vertex_x", "genKS_Vertex_y", "genKS_Vertex_z", "genKS_Vertex_n", "genKS_Vertex_r", "genKS_Vertex_acceptance_r",
                 "genKS_Vertex_acceptance_n", "genKS_Vertex_d", "genKS_Vertex_p", "genKS_Vertex_pt",
                 "EVT_NVertex",
                 "MC_Vertex_mass", "MC_Vertex_p",
                 "Vertex_x", "Vertex_y", "Vertex_z", "Vertex_r", "Vertex_d",
                 "Vertex_xErr", "Vertex_yErr", "Vertex_zErr", "Vertex_rErr", "Vertex_dErr", 
                 "Vertex_chi2", "Vertex_isPV", "Vertex_ntrk", "Vertex_mass",
                 "Vertex_n",
                 "Vertex_mass_before_VerDet", "Vertex_mass_within_VerDet", "Vertex_mass_within_DC", "Vertex_mass_super_displaced",
                 "Vertex_pt", "Vertex_eta",
                 "genKS_tracks_indices", "genKS_tracks_x", "genKS_tracks_y", "genKS_tracks_z", "MC_tracks_p","genKS_tracks_r", "genKS_tracks_d", 
                 "genKS_tracks_n",
                 "Vertex_MCx", "Vertex_MCy", "Vertex_MCz",
                 "Vertex_isMCKSpipi",
                 "recoKS_Vertex_r", "recoKS_Vertex_z", "recoKS_Vertex_mass",
                 "recoKS_Vertex_mass_before_VerDet", "recoKS_Vertex_mass_within_VerDet", "recoKS_Vertex_mass_beyond_VerDet", 
                 
                 "genKS_recoVertex_n", "genKS_recoVertex_x", "genKS_recoVertex_y", "genKS_recoVertex_z", "genKS_recoVertex_p","genKS_recoVertex_r", 
                 "genKS_recoVertex_d"
                ]
        return branchList
