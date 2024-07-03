#include <iostream>
#include <cmath>
#include <vector>
#include "Pythia8/Pythia.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/tools/Filter.hh"
#include "TH2F.h"
#include "TFile.h"
#include "TROOT.h"
#include "TCanvas.h"

using namespace Pythia8;
using namespace fastjet;
using namespace std;

double tau_calc(const fastjet::PseudoJet &jet, const fastjet::PseudoJet &subjet1, const fastjet::PseudoJet &subjet2);
double deltaR(const fastjet::PseudoJet &jet1, const fastjet::PseudoJet &jet2);

int main(int argc, char* argv[]) {
    Pythia pythia;
    pythia.readFile(argv[1]);
    pythia.init();

    double tau_uncluster;
    vector<double> tau_unc;

    double tau_parton;
    vector<double> tau_par;
    vector<double> outgoingParticles;

    // start reconstructing using anti-kt
    double R_akt = 0.5;
    JetDefinition recon_jet(antikt_algorithm, R_akt);

    double R_CA = 1.0;
    JetDefinition CA_jet(antikt_algorithm, R_CA);

    // defining for reclustering using kt algorithm
    double R_kt = 1.0;

    for (int iEvent = 0; iEvent <= 5000; ++iEvent) {  // event loop
        if (!pythia.next()) continue;

        vector<PseudoJet> particles;

        for (int i = 0; i < pythia.event.size(); ++i) {  // particle loop

            if(pythia.event[i].isFinal()){

                particles.push_back(PseudoJet(
                pythia.event[i].px(),
                pythia.event[i].py(),
                pythia.event[i].pz(),
                pythia.event[i].e()
                ));
                
                if(pythia.event[i].status() == 23){     // identify final state particle
                    outgoingParticles.push_back(i);
                }
            }

        }

        // cluster using anti-kt
        ClusterSequence cluster_antikt(particles, recon_jet);
        vector<PseudoJet> jet_antikt = sorted_by_pt(cluster_antikt.inclusive_jets());                              

        // choosing the leading jet
        PseudoJet leadingJet;

        vector<PseudoJet> leadingJet_values;
        for(const auto& jet : jet_antikt){
            if(jet.pt() > 300.0 && abs(jet.eta()) < 1.0){   // jets selection condition
                leadingJet_values.push_back(jet);
            }
        }

        vector<PseudoJet> leadingJet_sort = sorted_by_pt(leadingJet_values);
        vector<PseudoJet> leadingParts;
        if(!leadingJet_sort.empty()){
            leadingJet = leadingJet_sort[0];    // leading jet is the highest pT jet
            // start reclustering using kt algorithm
            leadingParts = leadingJet.constituents(); // new particles vector is comprised of the leading jet particles
        }

        double p = 1.0;
        JetDefinition genkt_jet(genkt_algorithm, R_kt, p);

        ClusterSequence cluster_kt(leadingParts, genkt_jet);
        vector<PseudoJet> jet_genkt = sorted_by_pt(cluster_kt.inclusive_jets());

        // new leading jet

        PseudoJet parent1, parent2;

        if(!jet_genkt.empty()){
            leadingJet = jet_genkt[0];    // new leading jet

            vector<PseudoJet> subjets = {leadingJet};   // desclustering vector

            while(!subjets.empty()){
            PseudoJet jet = subjets.back(); // jet = last contituint of subject
            subjets.pop_back(); // delete last constituint 
            
                if(cluster_kt.has_parents(jet, parent1, parent2)){
                    // add parents to the declustering vector
                    subjets.push_back(parent1);
                    subjets.push_back(parent2);

                    if(subjets.size() >= 2){  // garantee that there are at least 2 subjets for each calc
                        tau_uncluster = tau_calc(leadingJet, subjets[0], subjets[1]);
                        tau_unc.push_back(tau_uncluster);
                    }

                }

            }
        }

        // tag all final state produced by initial partons
        vector<PseudoJet> taggedParticles;

        for(int i : outgoingParticles){
            for(int j = 0; j < pythia.event.size(); ++j){
                if(pythia.event[j].mother1() == i || pythia.event[j].mother2() == i){
                    taggedParticles.push_back(PseudoJet(
                    pythia.event[i].px(),
                    pythia.event[i].py(),
                    pythia.event[i].pz(),
                    pythia.event[i].e()
                    ));
                    
                }
            }
        }

//        cout << taggedParticles.size() << endl;

        // R = 0.5
        // recluster anti-kt with R = 0.5
        double R_tau = 0.5;
        JetDefinition jetTau(antikt_algorithm, R_tau);

        ClusterSequence cs_tag(taggedParticles, jetTau);
//        vector<PseudoJet> jets_tag = sorted_by_pt(cs_tag.inclusive_jets());

//        PseudoJet leadingJetPar;
//        //cout << jets_tag.size() << endl;
//        for(const auto &jet : jets_tag){
//            if(jet.pt() > 300.0 && abs(jet.eta()) < 1.0){
//                leadingJetPar = jet;
//                break;
//            }
//        }

        // compare deltaR between the tagged jets and the leading jet
//        double minDeltaR = 1.0;
//        PseudoJet closestJet;
//        for(const auto &jet : jets_tag){
//            double dR = deltaR(leadingJetPar, jet);
//            if(dR < minDeltaR){
//                minDeltaR = dR;
//                closestJet = jet;
//                break;
//            }
//        }


        // R = 1.0
        // recluster anti-kt with R = 1.0
//        ClusterSequence cs_tagCA(taggedParticles, CA_jet);
//        vector<PseudoJet> jets_tagCA = sorted_by_pt(cs_tagCA.inclusive_jets());
//
//        PseudoJet leadingJetParCA;
//        for(const auto &jet : jets_tagCA){
//            if(jet.pt() > 300.0 && abs(jet.eta()) < 1.0){
//                leadingJetParCA = jet;
//                break;
//            }
//        }

        // compare deltaR between the tagged jets and the leading jet
//        PseudoJet closestJetCA;
//        for(const auto &jet : jets_tagCA){
//            double dR = deltaR(leadingJetParCA, jet);
//            if(dR < minDeltaR){
//                minDeltaR = dR;
//                closestJetCA = jet;
//            }
//        }

    }    

//    cout << "tau_unc size = " << tau_unc.size() << "; tau_par size = " << tau_par.size() << endl;

    pythia.stat();

    return 0;
}

double tau_calc(const fastjet::PseudoJet &jet, const fastjet::PseudoJet &subjet1, const fastjet::PseudoJet &subjet2){
    double E = jet.e();
    double z1 = subjet1.e()/jet.e();

    double cosTheta_12 = (subjet1.px()*subjet2.px()+subjet1.py()*subjet2.py()+subjet1.pz()*subjet2.pz())/sqrt((pow(subjet1.px(), 2)+pow(subjet1.py(), 2)+pow(subjet1.pz(), 2))*(pow(subjet2.px(), 2)+pow(subjet2.py(), 2)+pow(subjet2.pz(), 2)));

    double tau_form = 1/(2*E*z1*(1-z1)*(1-cosTheta_12));

    return tau_form;
}

double deltaR(const fastjet::PseudoJet &jet1, const fastjet::PseudoJet &jet2){
    double phi1 = jet1.phi();
    double phi2 = jet2.phi();
    double y1 = jet1.rap();
    double y2 = jet2.rap();

    double phi_dif = phi1 - phi2;
    double y_dif = y1 - y2;

    double delR = sqrt(pow(phi_dif, 2) + pow(y_dif, 2));

    return delR;
}