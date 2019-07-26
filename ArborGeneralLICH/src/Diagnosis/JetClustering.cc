#include <JetClustering.hh>
#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/LCFloatVec.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <IMPL/MCParticleImpl.h>
#include <EVENT/LCRelation.h>
#include <IMPL/LCRelationImpl.h>
#include <values.h>
#include <string>
#include <iostream>
#include <EVENT/LCFloatVec.h>
#include <EVENT/LCParameters.h>
#include <stdexcept>
#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>
#include <TRandom.h>
#include <Rtypes.h>
#include <sstream>
#include <cmath>
#include <vector>
#include <TMath.h>
#include "TLorentzVector.h"


#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/Config.h"

using namespace std;

JetClustering a_JetClustering_instance;

JetClustering::JetClustering()
: Processor("JetClustering"),
_output(0)
{
    _description = "Print MC Truth" ;
    
    _treeFileName="MCTruth.root";
    registerProcessorParameter( "TreeOutputFile" ,
                               "The name of the file to which the ROOT tree will be written" ,
                               _treeFileName ,
                               _treeFileName);
    
    _treeName="Tau";
    registerProcessorParameter( "TreeName" ,
                               "The name of the ROOT tree" ,
                               _treeName ,
                               _treeName);
    
    _overwrite=0;
    registerProcessorParameter( "OverwriteFile" ,
                               "If zero an already existing file will not be overwritten." ,
                               _overwrite ,
                               _overwrite);
    
}

void findCombine(int index, vector<int> &a, vector<int>& tmp, const int& n, vector<vector<int> >& res){
	if(tmp.size() == n/2){
		res.push_back(tmp);
		return;
	}
	for(int i=index; i<=n; i++){
		tmp.push_back(a[i-1]);
		findCombine(i+1, a, tmp, n, res);
		tmp.pop_back();
	}
	return;
}

vector<vector<int> > Combine(int n, vector<int> a){
	vector<vector<int> > res;
	if(n<1 || n%2!=0 || a.size()<1)
		return res;
	vector<int> tmp;
	findCombine(1, a, tmp, n, res);
	return res;
}

vector<int> DifferenceSet(vector<int> total, vector<int> original, vector<int> left){
    for(int i = 0; i<total.size(); i++)
    {
        int element = total[i];
        vector<int >::iterator it;
        it = find(original.begin(), original.end(), element);
        if(it == original.end()){
            left.push_back(element);
        }
    }
    return left;
}


void JetClustering::init() {
    
    printParameters();
    TFile *tree_file=new TFile(_treeFileName.c_str(),(_overwrite ? "RECREATE" : "UPDATE"));
    if (!tree_file->IsOpen()) {
        delete tree_file;
        tree_file=new TFile(_treeFileName.c_str(),"NEW");
    }
    
    _outputTree = new TTree(_treeName.c_str(),_treeName.c_str());
    _outputTree->SetAutoSave(32*1024*1024);  // autosave every 32MB
    _outputTree->Branch("EventNr", &eventNr, "EventNr/I");
    _outputTree->Branch("Num", &Num, "Num/I");
    _outputTree->Branch("eventType", &eventType, "eventType/I");
    _outputTree->Branch("GeventType", &GeventType, "GeventType/I");
    
    _outputTree->Branch("vP_Px", &vP_Px);
    _outputTree->Branch("vP_Py", &vP_Py);
    _outputTree->Branch("vP_Pz", &vP_Pz);
    _outputTree->Branch("vP_En", &vP_En);

    _outputTree->Branch("vB_Px", &vB_Px);
    _outputTree->Branch("vB_Py", &vB_Py);
    _outputTree->Branch("vB_Pz", &vB_Pz);
    _outputTree->Branch("vB_En", &vB_En);
    
    _outputTree->Branch("vG_Px", &vG_Px);
    _outputTree->Branch("vG_Py", &vG_Py);
    _outputTree->Branch("vG_Pz", &vG_Pz);
    _outputTree->Branch("vG_En", &vG_En);
    
    _outputTree->Branch("vR_Px", &vR_Px);
    _outputTree->Branch("vR_Py", &vR_Py);
    _outputTree->Branch("vR_Pz", &vR_Pz);
    _outputTree->Branch("vR_En", &vR_En);
    
    _outputTree->Branch("count_c", &count_c, "count_c/I");
    _outputTree->Branch("count_u", &count_u, "count_u/I");
    
    _outputTree->Branch("mass_B1", &mass_B1, "mass_B1/F");
    _outputTree->Branch("mass_B2", &mass_B2, "mass_B2/F");
    _outputTree->Branch("mass_G1", &mass_G1, "mass_G1/F");
    _outputTree->Branch("mass_G2", &mass_G2, "mass_G2/F");
    _outputTree->Branch("mass_R1", &mass_R1, "mass_R1/F");
    _outputTree->Branch("mass_R2", &mass_R2, "mass_R2/F");
    
    _outputTree->Branch("misMinEn", &misMinEn, "misMinEn/D");
    _outputTree->Branch("ISR_En", &ISR_En, "ISR_En/D");
    
    _outputTree->Branch("count_boson", &count_boson, "count_boson/I");
    _outputTree->Branch("Pselect_combi", &Pselect_combi, "Pselect_combi/I");
    _outputTree->Branch("En_minPart1", &En_minPart1, "En_minPart1/F");
    _outputTree->Branch("En_minPart2", &En_minPart2, "En_minPart2/F");
    
    _outputTree->Branch("CParameter", &CParameter, "CParameter/D");
    _outputTree->Branch("hemiMass1", &hemiMass1, "hemiMass1/D");
    _outputTree->Branch("hemiMass2", &hemiMass2, "hemiMass2/D");
    _outputTree->Branch("hemiBroadening1", &hemiBroadening1, "hemiBroadening1/D");
    _outputTree->Branch("hemiBroadening2", &hemiBroadening2, "hemiBroadening2/D");
    _outputTree->Branch("T", &T, "T/D");

    
    Num = 0;
}

void JetClustering::processEvent( LCEvent * evtP )
{
    
    if (evtP)
    {
        try{
            
            vB_Px.clear();
            vB_Py.clear();
            vB_Pz.clear();
            vB_En.clear();
            
            vG_Px.clear();
            vG_Py.clear();
            vG_Pz.clear();
            vG_En.clear();
            
            vR_Px.clear();
            vR_Py.clear();
            vR_Pz.clear();
            vR_En.clear();
            
            vP_Px.clear();
            vP_Py.clear();
            vP_Pz.clear();
            vP_En.clear();
            
            cout<<"Next Event ***********************"<<endl;
            eventNr = evtP->getEventNumber();
            cout<<"eventNr : "<<eventNr<<" Num : "<<Num<<endl;
            
            double k1 = 0.0;
            double k2 = 0.0;
            double wmass = 80.4, zmass = 91.2, hmass = 125;
            
            //ArborPFOs
            LCCollection* col_PFO = evtP->getCollection( "ArborPFOs" );
            int nPFO = col_PFO->getNumberOfElements();
//            cout<<"nPFO : "<<nPFO<<endl;
            TLorentzVector TLPFO(0,0,0,0);
            for(int i = 0; i<nPFO; i++){
                ReconstructedParticle* pfo = dynamic_cast<ReconstructedParticle*>(col_PFO->getElementAt(i));
                TLorentzVector temp(pfo->getMomentum(), pfo->getEnergy());
                TLPFO += temp;
            }
            cout<<"total energy of ArborPFOs is "<<TLPFO.E()<<endl;
            

            
            //RecoJet
            LCCollection* col_Jet = evtP->getCollection( "FastJets" );
            int num_jet = col_Jet->getNumberOfElements();
            cout<<"the number of recojet is "<<num_jet<<endl;
            vector<int > a(num_jet);
            for(int i = 0; i<num_jet; i++){
                a[i] = i;
            }
            vector<vector<int> > res;
            res = Combine(num_jet, a);
            
            vector<vector<int> > vect1;
            for(int i = 0; i<(num_jet-1); i++)
            {
                for(int j = i+1; j<num_jet; j++)
                {
                        vector<int> temp(2);
                        vector<int> result;
                        
                        temp[0] = a[i];
                        temp[1] = a[j];
                        vect1.push_back(temp);
                        result = DifferenceSet(a, temp, result);
                        vect1.push_back(result);
                }
            }
            
            
            //used for pairing of 4 partons
            vector<int> b(4);
            for(int i=0; i<4; i++){
                b[i] = i;
            }
            vector<vector<int>> bres;
            bres = Combine(4, b);
            vector<vector<int>> vectP;
            for(int i = 0; i<(4-1); i++){
                for(int j = i+1; j<4; j++){
                    vector<int> temp(2);
                    vector<int> result;
                    temp[0] = b[i];
                    temp[1] = b[j];
                    vectP.push_back(temp);
                    result = DifferenceSet(b, temp, result);
                    vectP.push_back(result);
                }
            }
            
            
            
            
            TLorentzVector RTL1(0,0,0,0);
            TLorentzVector RTL2(0,0,0,0);
            
            double minDif = 9999.0;
            int select_combi = 0;
            int eventType = 9;    // WW = 0, ZZ = 1
            for(int i=0; i<(vect1.size() - 1); i = i+2){
                vector<ReconstructedParticle* > group1;
                vector<ReconstructedParticle* > group2;
                vector<TLorentzVector > TLVGroup1;
                vector<TLorentzVector > TLVGroup2;
                group1.clear();     group2.clear();
                TLVGroup1.clear(); TLVGroup2.clear();
                int charge_one = 0;
                int charge_two = 0;
                int total_RecoParticles = 0;
                for(int j=0; j<vect1[i].size(); j++){
                    ReconstructedParticle* jet = dynamic_cast<ReconstructedParticle*>(col_Jet->getElementAt(vect1[i][j]));
                    group1.push_back(jet);
                    ReconstructedParticleVec components = jet->getParticles();
                    int ncomps = components.size();
                    total_RecoParticles += ncomps;
                    for(int icomp = 0; icomp<ncomps; icomp++)
                    {
                        ReconstructedParticle* compPar = components.at(icomp);
                        charge_one += compPar->getCharge();
                    }
                }
                for(int j=0; j<vect1[i+1].size(); j++){
                    ReconstructedParticle* jet = dynamic_cast<ReconstructedParticle*>(col_Jet->getElementAt(vect1[i+1][j]));
                    group2.push_back(jet);
                    ReconstructedParticleVec components = jet->getParticles();
                    int ncomps = components.size();
                    total_RecoParticles += ncomps;
                    for(int icomp = 0; icomp<ncomps; icomp++)
                    {
                        ReconstructedParticle* compPar = components.at(icomp);
                        charge_two += compPar->getCharge();
                    }
                }
//                cout<<group1.size()<<" : "<<group2.size()<<endl;
//                cout<<"In FastJets, total particles is "<<total_RecoParticles<<endl;
                TLorentzVector TLGroup1(0,0,0,0);
                TLorentzVector TLGroup2(0,0,0,0);
                
                for(int k=0; k<group1.size(); k++){
                    ReconstructedParticle* recojet = group1.at(k);
                    TLorentzVector TLRecoJet(recojet->getMomentum(), recojet->getEnergy());
                    TLGroup1 += TLRecoJet;
                    TLVGroup1.push_back(TLRecoJet);
                }
                for(int k=0; k<group2.size(); k++){
                    ReconstructedParticle* recojet = group2.at(k);
                    TLorentzVector TLRecoJet(recojet->getMomentum(), recojet->getEnergy());
                    TLGroup2 += TLRecoJet;
                    TLVGroup2.push_back(TLRecoJet);
                }

                double compareW = (TLGroup1.M() - wmass)*(TLGroup1.M() - wmass)/14.44 + (TLGroup2.M() - wmass)*(TLGroup2.M() - wmass)/14.44 + k1 * pow(abs(TLGroup1.M() - TLGroup2.M()),2) + k2 * pow(abs(charge_one) - 1 + abs(charge_two) - 1,2);
                double compareZ = (TLGroup1.M() - zmass)*(TLGroup1.M() - zmass)/19.36 + (TLGroup2.M() - zmass)*(TLGroup2.M() - zmass)/19.36 + k1 * pow(abs(TLGroup1.M() - TLGroup2.M()),2) + k2 * pow(abs(charge_one) + abs(charge_two),2);
                
                if(compareW < minDif){minDif = compareW; select_combi = i; eventType = 0; RTL1 = TLGroup1; RTL2 = TLGroup2;}
                if(compareZ < minDif){minDif = compareZ; select_combi = i; eventType = 1; RTL1 = TLGroup1; RTL2 = TLGroup2;}
                
            }
            vR_Px.push_back(RTL1.Px()); vR_Py.push_back(RTL1.Py()); vR_Pz.push_back(RTL1.Pz()); vR_En.push_back(RTL1.E());
            vR_Px.push_back(RTL2.Px()); vR_Py.push_back(RTL2.Py()); vR_Pz.push_back(RTL2.Pz()); vR_En.push_back(RTL2.E());
            mass_R1 = RTL1.M(); mass_R2 = RTL2.M();
            cout<<"the total energy if FastJet is "<<(RTL1 + RTL2).E()<<endl;
            
            
            // MCParticle
            std::vector<MCParticle* >MCPUsedForThrust;

            LCCollection* col_MCP = evtP->getCollection( "MCParticle" );
            std::vector<MCParticle* > MCTQuark;
            TLorentzVector TLMCP(0,0,0,0);
            TLorentzVector TLISR(0,0,0,0);
            ISR_En = 0;
            double mass_941 = 0;
            int n_MCP = col_MCP->getNumberOfElements();
            count_boson = 0;
            int total_MCPFLParticles = 0;
            TLorentzVector TLBoson1(0,0,0,0);    TLorentzVector TLBoson2(0,0,0,0);   TLorentzVector TLBosonTest(0,0,0,0);
            mass_B1 = 0; mass_B2 = 0;
            MCParticle* MCP941;
            MCParticle* MCP942;
            
            for(int i = 0; i<n_MCP; i++)
            {
                MCParticle* a_MCP = dynamic_cast<MCParticle*>(col_MCP->getElementAt(i));
                int NParents = a_MCP->getParents().size();
                int NDaughters = a_MCP->getDaughters().size();
                int PDG = a_MCP->getPDG();
                TLorentzVector temp(a_MCP->getMomentum(), a_MCP->getEnergy());
                TVector3 temp_TV = a_MCP->getMomentum();
                //find the original two bosons, energy and costheta
                if(PDG == 94)
                {
                    count_boson += 1;
                    MCParticle* Par941 = a_MCP->getParents()[0]; MCParticle* Par942 = a_MCP->getParents()[1];
                    MCParticle* PPar941 = Par941; MCParticle* PPar942 = Par942;
                    do{PPar941 = PPar941->getParents()[0];}
                    while(PPar941->getParents().size() != 0);
                    do{PPar942 = PPar942->getParents()[0];}
                    while(PPar942->getParents().size() != 0);
                    if(PPar941->getParents().size() == 0 && PPar942->getParents().size() == 0)
                    {
                        cout<<"the PDG of two quarks merge to 94 is "<<PPar941->getPDG()<<" : "<<PPar942->getPDG()<<endl;
                        MCTQuark.push_back(PPar941); MCTQuark.push_back(PPar942); //MCTQuark store the four original quarks
                        TLorentzVector boson1(PPar941->getMomentum(), PPar941->getEnergy());
                        TLorentzVector boson2(PPar942->getMomentum(), PPar942->getEnergy());
                        TLBosonTest = boson1 + boson2;
                    }
                    TLorentzVector TL_B11(PPar941->getMomentum(), PPar941->getEnergy());
                    TLorentzVector TL_B12(PPar942->getMomentum(), PPar942->getEnergy());
                    TLorentzVector TL_B1(0,0,0,0);
                    TL_B1 = TL_B11 + TL_B12;
                    if(count_boson == 1)
                    {
                        MCP941 = a_MCP;
                        TLBoson1 = TLBosonTest;
                        vB_Px.push_back(TL_B1.Px()); vB_Py.push_back(TL_B1.Py()); vB_Pz.push_back(TL_B1.Pz()); vB_En.push_back(TL_B1.E());
                        mass_B1 = TL_B1.M();
                        mass_941 = a_MCP->getMass();
                        cout<<"mass_B1 : "<<mass_B1<<endl;
                    }
                    else if(count_boson == 2)
                    {
                        MCP942 = a_MCP;
                        TLBoson2 = TLBosonTest;
                        vB_Px.push_back(TL_B1.Px()); vB_Py.push_back(TL_B1.Py()); vB_Pz.push_back(TL_B1.Pz()); vB_En.push_back(TL_B1.E());
                        mass_B2 = TL_B1.M();
                        cout<<"mass_B2 : "<<mass_B2<<endl;
                    }
                    
                }
                
//                if(NParents == 0 && (abs(PDG) == 1 || abs(PDG) == 2 || abs(PDG) == 3 || abs(PDG) == 4 || abs(PDG) == 5 || abs(PDG) == 6) )
//                {
//                    MCTQuark.push_back(a_MCP);
//                }
                
                if(a_MCP->getGeneratorStatus() == 1)
                {
                    if(abs(PDG) != 12 || abs(PDG) != 14 || abs(PDG) != 16){
                        total_MCPFLParticles += 1;
                        TLMCP += temp;
                    }
                    
                    if(NDaughters == 0 && PDG != 22){
                        MCPUsedForThrust.push_back(a_MCP);
                    }
                    else if(NDaughters == 0 && PDG == 22){
                        MCParticle* a_parent = a_MCP;
                        do{a_parent = a_parent->getParents()[0];}
                        while(a_parent->getParents().size() != 0);
                        if(a_parent->getParents().size() == 0 && a_parent->getPDG() != 22){
                            MCPUsedForThrust.push_back(a_MCP);
                        }
                    }
                }
                    
            
                
                if(NParents == 0 && PDG == 22)
                {
                    TLISR += temp;
                }
            }
            ISR_En = TLISR.E();
//            cout<<"the ISR energy in MCParticles is "<<TLISR.E()<<endl;
//            cout<<"the total energy in MCParticles, except neutrinos, is "<<TLMCP.E()<<endl;
            TVector3 TVBoson1 = TLBoson1.Vect();
            TVector3 TVBoson2 = TLBoson2.Vect();
//            cout<<"the angle betwoon two boson is "<<TVBoson1.Angle(TVBoson2)<<endl;
            
            
            vector<TLorentzVector > TLMCTQuark;
            count_u = 0; count_c = 0;
            if(MCTQuark.size() == 4)
            {
                for(int i = 0; i< 4; i++)
                {
                    MCParticle* quark = MCTQuark.at(i);
//                    TLorentzVector temp(quark->getMomentum(), quark->getEnergy());
                    cout<<"the PDG of original quark is "<<quark->getPDG()<<endl;
                    if(abs(quark->getPDG()) == 2 ){count_u += 1;}
                    else if(abs(quark->getPDG()) == 4 ){count_c += 1;}
                    
                    double noise1 = gRandom->Gaus(0, 0.05);
                    double noise2 = gRandom->Gaus(0, 0.05);
                    double noise3 = gRandom->Gaus(0, 0.05);
                    double noise4 = gRandom->Gaus(0, 0.05);
                    TLorentzVector temp((quark->getMomentum()[0])*(noise1 + 1), (quark->getMomentum()[1])*(noise2 + 1), (quark->getMomentum()[2])*(noise3 + 1), (quark->getEnergy())*(1 + noise4));
                    TLMCTQuark.push_back(temp);
                    
                }
                
                cout<<"test2******"<<endl;
                
                
                int num_quark = MCTQuark.size();
                double PminDif = 9999.0;
                Pselect_combi = 99;
                int PeventType = 9;
                TLorentzVector TLParton12(0,0,0,0);
                TLorentzVector TLParton34(0,0,0,0);
                for(int i = 0; i<(vectP.size() - 1); i = i + 2){
                    TLorentzVector TLparton12(0,0,0,0);
                    TLorentzVector TLparton34(0,0,0,0);
                    
                    for(int j = 0; j<vectP[i].size(); j++){
//                        MCParticle* quark = MCTQuark.at(vect1[i][j]);
//                        TLorentzVector temp(quark->getMomentum(), quark->getEnergy());
//                        TLorentzVector temp
//                        TLparton12 += temp;
                        TLparton12 += TLMCTQuark.at(vectP[i][j]);
                    }
                    for(int j = 0; j<vectP[i+1].size(); j++){
//                        MCParticle* quark = MCTQuark.at(vect1[i+1][j]);
//                        TLorentzVector temp(quark->getMomentum(), quark->getEnergy());
//                        TLparton34 += temp;
                        TLparton34 += TLMCTQuark.at(vectP[i+1][j]);
                    }
                    
                    double compareW = (TLparton12.M() - wmass)*(TLparton12.M() - wmass)/14.44 + (TLparton34.M() - wmass)*(TLparton34.M() - wmass)/14.44;
                    double compareZ = (TLparton12.M() - zmass)*(TLparton12.M() - zmass)/19.36 + (TLparton34.M() - zmass)*(TLparton34.M() - zmass)/19.36;
                    
                    if(compareW < PminDif){PminDif = compareW; Pselect_combi = i; PeventType = 0; TLParton12 = TLparton12; TLParton34 = TLparton34;}
                    if(compareZ < PminDif){PminDif = compareZ; Pselect_combi = i; PeventType = 1; TLParton12 = TLparton12; TLParton34 = TLparton34;}
                    
                    cout<<"PminDif : "<<PminDif<<" Pselect_combi : "<<Pselect_combi<<" PeventType : "<<PeventType<<endl;
                    
                }
                vP_Px.push_back(TLParton12.Px()); vP_Py.push_back(TLParton12.Py()); vP_Pz.push_back(TLParton12.Pz()); vP_En.push_back(TLParton12.E());
                vP_Px.push_back(TLParton34.Px()); vP_Py.push_back(TLParton34.Py()); vP_Pz.push_back(TLParton34.Pz()); vP_En.push_back(TLParton34.E());
                cout<<"TLParton12.M() : "<<TLParton12.M()<<" TLParton34.M() : "<<TLParton34.M()<<endl;
                cout<<"vectP[Pselect_combi][0] : "<<vectP[Pselect_combi][0]<<" vectP[Pselect_combi][1] : "<<vectP[Pselect_combi][1]<<endl;
                cout<<"vectP[Pselect_combi+1][0] : "<<vectP[Pselect_combi+1][0]<<" vectP[Pselect_combi+1][1] : "<<vectP[Pselect_combi+1][1]<<endl;
                
                
                
                //the following code used to find the thrust
                double thetaMin = TMath::Pi(), phiMin = 2*TMath::Pi();
                double thetaMin2 = 0, phiMin2 = 0;
                double thetaL = 0, phiL = 0, thetaR = TMath::Pi(), phiR = 2*TMath::Pi();
                int iter = 0;
                double Told = 0;
                double Tnew = 0;
                double cut = 1;
                double thetaRange = 0, phiRange = 0;
                do{
                    iter += 1;
                    cout<<"iter : "<<iter<<endl;
                    if(iter == 1){
                        thetaRange = thetaR - thetaL, phiRange = phiR - phiL;
                    }
                    else if(iter != 1){
                        
                        thetaRange = 0.1*(thetaR - thetaL);
                        phiRange = 0.1*(phiR - phiL);
                        
                        thetaL =  thetaMin - thetaRange;
                        thetaR = thetaMin + thetaRange;
                        phiL = phiMin - phiRange;
                        phiR = phiMin + phiRange;
                        thetaRange = thetaR - thetaL, phiRange = phiR - phiL;
                        
                        cout<<"thetaL : "<<thetaL<<" thetaR : "<<thetaR<<endl;
                        cout<<"phiL : "<<phiL<<" phiR : "<<phiR<<endl;
                    }
                    
                    cout<<"thetaRange : "<<thetaRange<<" phiRange : "<<phiRange<<endl;
                    for(double theta = thetaL; theta <= thetaR; theta += 0.1*thetaRange){   //in this round, find the max T
                        for(double phi = phiL; phi <= phiR; phi += 0.1*phiRange){
                            
                            double x = sin(theta)*cos(phi);
                            double y = sin(theta)*sin(phi);
                            double z = cos(theta);
                            
                            double denominator = 0;
                            double numerator = 0;
                            for(int i = 0; i<MCPUsedForThrust.size(); i++){
                                MCParticle* temp = MCPUsedForThrust.at(i);
                                TLorentzVector TLtemp(temp->getMomentum(), temp->getEnergy());
                                TVector3 TVtemp = TLtemp.Vect();
                                denominator += TVtemp.Mag();
                                numerator += abs(x*TVtemp(0) + y*TVtemp(1) + z*TVtemp(2));
                            }
                            double Ttemp = numerator/denominator;
                            if(Ttemp > T){
                                thetaMin = theta;   phiMin = phi; T = Ttemp;
                                cout<<"*************"<<endl;
                                cout<<"T : "<<T<<"thetaMin : phiMin "<<thetaMin<<" : "<<phiMin<<endl;
                                cout<<"*************"<<endl;
                            }
                        }
                    }
                    
                    
                    if(iter == 1){Told = T; Tnew = T;}
                    else if(T >= Tnew && iter != 1){
                        Told = Tnew; Tnew = T; cut = (Tnew - Told)/Tnew;
                    }
                    cout<<"cut : "<<cut<<endl;
                }
                while(cut >= 0.2);
                
                
                TVector3 tempThrust(0,0,0);
                tempThrust.SetXYZ(sin(thetaMin)*cos(phiMin), sin(thetaMin)*sin(phiMin), cos(thetaMin));
                TVector3 backThrust(0,0,0);
                backThrust = -tempThrust;
                
                //the following code used to get Hemisphere masses
                std::vector<MCParticle* > hemisphere1;
                std::vector<MCParticle* > hemisphere2;
                double visEn = 0;
                double JetBroadeningDenominator = 0;
                for(int i = 0; i<MCPUsedForThrust.size(); i++){
                    MCParticle* a_MCP = MCPUsedForThrust.at(i);
                    TLorentzVector TLtemp(a_MCP->getMomentum(), a_MCP->getEnergy());
                    TVector3 TVtemp = TLtemp.Vect();
                    if(TVtemp.Angle(tempThrust) > 0.5*TMath::Pi()){hemisphere1.push_back(a_MCP);}
                    else {hemisphere2.push_back(a_MCP);}
                    visEn += a_MCP->getEnergy();
                    JetBroadeningDenominator += TVtemp.Mag();
                }
                
                cout<<"the number of particles in two hemispheres is "<<hemisphere1.size()+hemisphere2.size()<<endl;
                
                
                double JetBroadeningNumerator1 = 0, JetBroadeningNumerator2 =0;
                TLorentzVector TLHemi1(0,0,0,0);
                TLorentzVector TLHemi2(0,0,0,0);
                for(int i = 0; i<hemisphere1.size(); i++){
                    MCParticle* a_MCP = hemisphere1.at(i);
                    TLorentzVector TLtemp(a_MCP->getMomentum(), a_MCP->getEnergy());
                    TVector3 TVtemp = TLtemp.Vect();
                    TLHemi1 += TLtemp;
                    JetBroadeningNumerator1 += abs(TVtemp.Mag() * tempThrust.Mag() * sin(TVtemp.Angle(tempThrust)));
                    
                }
                for(int i = 0; i<hemisphere2.size(); i++){
                    MCParticle* a_MCP = hemisphere2.at(i);
                    TLorentzVector TLtemp(a_MCP->getMomentum(), a_MCP->getEnergy());
                    TVector3 TVtemp = TLtemp.Vect();
                    TLHemi2 += TLtemp;
                    JetBroadeningNumerator2 += abs(TVtemp.Mag() * tempThrust.Mag() * sin(TVtemp.Angle(tempThrust)));
                }
                hemiMass1 = (TLHemi1.M())*(TLHemi1.M())/(visEn*visEn);
                hemiMass2 = (TLHemi2.M())*(TLHemi2.M())/(visEn*visEn);
                hemiBroadening1 = JetBroadeningNumerator1/(2*JetBroadeningDenominator);
                hemiBroadening2 = JetBroadeningNumerator2/(2*JetBroadeningDenominator);
                cout<<"hemiMass1 : "<<hemiMass1<<" hemiMass2 : "<<hemiMass2<<endl;
                cout<<"hemiBroadening1 : "<<hemiBroadening1<<" hemiBroadening2 : "<<hemiBroadening2<<endl;
                
                
                //the following code used to calculate sphericity
                
                double p2Min = 1e-20;
                double denom = 0.;
                double tt[4][4];
                double eVal1 = 0, eVal2 = 0, eVal3 = 0;
                for(int j = 1; j<4; ++j){
                    for(int k = j; k<4; ++k){
                        tt[j][k] = 0.;
                    }
                }
                for(int i = 0; i<MCPUsedForThrust.size(); i++){
                    MCParticle* a_MCP = MCPUsedForThrust.at(i);
                    TLorentzVector TLtemp(a_MCP->getMomentum(), a_MCP->getEnergy());
                    TVector3 TVtemp = TLtemp.Vect();
                    double pNow[4];
                    pNow[1] = TLtemp.Px();
                    pNow[2] = TLtemp.Py();
                    pNow[3] = TLtemp.Pz();
                    double p2Now = pow(pNow[1],2) + pow(pNow[2],2) + pow(pNow[3],2);
                    double pWeight = 1./sqrt(max(p2Now, p2Min));
                    for(int j = 1; j<4; ++j){
                        for(int k=j; k<4; ++k){
                            tt[j][k] += pWeight * pNow[j] * pNow[k];
                            denom += pWeight * p2Now;
                        }
                    }
                }
                //Normalize tensor to trace = 1.
                for(int j = 1; j < 4; ++j){
                    for(int k = j; k < 4; ++k){
                        tt[j][k] /= denom;
                    }
                }
                //find eigenvalues to matrix (third degree equation)
                double qCoef = ( tt[1][1] * tt[2][2] + tt[1][1] * tt[3][3]
                                + tt[2][2] * tt[3][3] - pow(tt[1][2],2) - pow(tt[1][3],2)
                                - pow(tt[2][3],2) ) / 3. - 1./9.;
                double qCoefRt = sqrt( -qCoef);
                double rCoef = -0.5 * ( qCoef + 1./9. + tt[1][1] * pow(tt[2][3],2)
                                       + tt[2][2] * pow(tt[1][3],2) + tt[3][3] * pow(tt[1][2],2)
                                       - tt[1][1] * tt[2][2] * tt[3][3] )
                + tt[1][2] * tt[1][3] * tt[2][3] + 1./27.;
                double pTemp = max( min( rCoef / pow(qCoefRt,3), 1.), -1.);
                double pCoef = cos( acos(pTemp) / 3.);
                double pCoefRt = sqrt( 3. * (1. - pow(pCoef,2)) );
                eVal1 = 1./3. + qCoefRt * max( 2. * pCoef,  pCoefRt - pCoef);
                eVal3 = 1./3. + qCoefRt * min( 2. * pCoef, -pCoefRt - pCoef);
                eVal2 = 1. - eVal1 - eVal3;
                CParameter = 3*(eVal1*eVal2 + eVal2*eVal3 + eVal3*eVal1);
                
                cout<<"C parameter is "<<CParameter<<endl;
                
                
                
                
            }
            
            
//            cout<<"count_u : count_c "<<count_u<<" : "<<count_c<<endl;
//            cout<<"total_MCPFLParticles : "<<total_MCPFLParticles<<endl;
            
            

            
            //matching four original quarks to two groups

            
            
            cout<<"test1********"<<endl;
            
            
            // GenJet
            // select the most suited grouping method
            cout<<"GenJet *****"<<endl;
            LCCollection* MCP_Jet = evtP->getCollection( "MCPFastJets" );
            int n_GenJet = MCP_Jet->getNumberOfElements();
            for(int i = 0; i<n_GenJet; i++){
                ReconstructedParticle* jet = dynamic_cast<ReconstructedParticle*>(MCP_Jet->getElementAt(i));
                cout<<"energy of GenJet is "<<jet->getEnergy()<<endl;
            }

            
            
            LCCollection* col_Rela = evtP->getCollection( "MCPSIMURelation" );
            int NLink = col_Rela->getNumberOfElements();
            TLorentzVector GTL1(0,0,0,0);
            TLorentzVector GTL2(0,0,0,0);
            double GminDif = 9999.0;
            int GeventType = 99999;
            int Gselect_combi = 0;
            double GminiAngle = 0;
            for(int i=0; i<(vect1.size() - 1); i = i+2){
                vector<ReconstructedParticle* > group1;
                vector<ReconstructedParticle* > group2;
                group1.clear();     group2.clear();
                int num_GenJetParticles = 0;
                int charge_one = 0;
                int charge_two = 0;
                for(int j=0; j<vect1[i].size(); j++){
                    ReconstructedParticle* jet = dynamic_cast<ReconstructedParticle*>(MCP_Jet->getElementAt(vect1[i][j]));
                    group1.push_back(jet);
                    ReconstructedParticleVec components = jet->getParticles();
                    int ncomps = components.size();
                    num_GenJetParticles += ncomps;
                    for(int icomp = 0; icomp<ncomps; icomp++)
                    {
                        ReconstructedParticle* compPar = components.at(icomp);
                        for(int k=0; k<NLink; k++)
                        {
                            LCRelation* a_link = dynamic_cast<LCRelation*>(col_Rela->getElementAt(k));
                            if(a_link->getFrom() == compPar)
                            {
                                MCParticle* linkto = dynamic_cast<MCParticle*>(a_link->getTo());
                                charge_one += linkto->getCharge();
                            }
                        }
                    }
                }
                for(int j=0; j<vect1[i+1].size(); j++){
                    ReconstructedParticle* jet = dynamic_cast<ReconstructedParticle*>(MCP_Jet->getElementAt(vect1[i+1][j]));
                    group2.push_back(jet);
                    ReconstructedParticleVec components = jet->getParticles();
                    int ncomps = components.size();
                    num_GenJetParticles += ncomps;
                    for(int icomp = 0; icomp<ncomps; icomp++)
                    {
                        ReconstructedParticle* compPar = components.at(icomp);
                        for(int k=0; k<NLink;k++)
                        {
                            LCRelation* a_link = dynamic_cast<LCRelation*>(col_Rela->getElementAt(k));
                            if(a_link->getFrom() == compPar)
                            {
                                MCParticle* linkto = dynamic_cast<MCParticle*>(a_link->getTo());
                                charge_two += linkto->getCharge();
                            }
                        }
                    }
                }
//                cout<<"group1.size() : group2.size() "<<group1.size()<<" : "<<group2.size()<<endl;
//                cout<<"In MCPFastJets, total particles is "<<num_GenJetParticles<<endl;
                
                TLorentzVector TLGroup1(0,0,0,0);
                TLorentzVector TLGroup2(0,0,0,0);
                for(int k=0; k<group1.size(); k++){
                    ReconstructedParticle* recojet = group1.at(k);
                    TLorentzVector TLRecoJet(recojet->getMomentum(), recojet->getEnergy());
                    TLGroup1 += TLRecoJet;
                }
                for(int k=0; k<group2.size(); k++){
                    ReconstructedParticle* recojet = group2.at(k);
                    TLorentzVector TLRecoJet(recojet->getMomentum(), recojet->getEnergy());
                    TLGroup2 += TLRecoJet;
                }
                
                TVector3 TVGroup1 = TLGroup1.Vect();
                TVector3 TVGroup2 = TLGroup2.Vect();
//                cout<<"TVGroup1.Angle(TVGroup2) : "<<TVGroup1.Angle(TVGroup2)<<endl;
//                cout<<"TLGroup1.M() : "<<TLGroup1.M()<<" TLGroup2.M() : "<<TLGroup2.M()<<endl;
                
                double compareW = (TLGroup1.M() - wmass)*(TLGroup1.M() - wmass)/14.44 + (TLGroup2.M() - wmass)*(TLGroup2.M() - wmass)/14.44 + k1 * pow(abs(TLGroup2.M() - TLGroup1.M()),2) + k2 * pow(abs(charge_one) - 1 + abs(charge_two) - 1, 2);
                double compareZ = (TLGroup1.M() - zmass)*(TLGroup1.M() - zmass)/19.36 + (TLGroup2.M() - zmass)*(TLGroup2.M() - zmass)/19.36 + k1 * pow(abs(TLGroup2.M() - TLGroup1.M()),2) + k2 * pow(abs(charge_one) + abs(charge_two), 2);
                if(compareW < GminDif){GminDif = compareW; GeventType = 0; GTL1 = TLGroup1; GTL2 = TLGroup2; Gselect_combi = i; GminiAngle = TVGroup1.Angle(TVGroup2);}
                if(compareZ < GminDif){GminDif = compareZ; GeventType = 1; GTL1 = TLGroup1; GTL2 = TLGroup2;  Gselect_combi = i; GminiAngle = TVGroup1.Angle(TVGroup2);}
//                cout<<"the mass difference is "<<(TLGroup1.M() - TLGroup2.M())<<endl;
//                cout<<"charge_one : charge_two "<<charge_one<<" : "<<charge_two<<endl;
//                cout<<"compareW : compareZ "<<compareW<<" : "<<compareZ<<endl;
//                cout<<"chi2 : "<<GminDif<<endl;
                cout<<"the total energy of GenJets is "<<(TLGroup2 + TLGroup1).E()<<endl;
            }
            
//            cout<<"Gselect_combi : "<<Gselect_combi<<endl;
//            cout<<"GminiAngle : "<<GminiAngle<<endl;
            
            vG_Px.push_back(GTL1.Px()); vG_Py.push_back(GTL1.Py()); vG_Pz.push_back(GTL1.Pz()); vG_En.push_back(GTL1.E());
            vG_Px.push_back(GTL2.Px()); vG_Py.push_back(GTL2.Py()); vG_Pz.push_back(GTL2.Pz()); vG_En.push_back(GTL2.E());
            mass_G1 = GTL1.M(); mass_G2 = GTL2.M();
            
            
            float En_no94 = 0;
            float test_EnTotal = 0;
            En_minPart1 = 0;
            int num_group1 = vect1[Gselect_combi].size();
            int num_group2 = vect1[Gselect_combi + 1].size();
            cout<<"num_group1 : num_group2 "<<num_group1<<" : "<<num_group2<<endl;
            for(int i = 0; i < num_group1; i++){
                cout<<"the first time >>>>>>>"<<endl;
                
                float En_B1Part = 0, En_B2Part = 0;
                ReconstructedParticle* jet = dynamic_cast<ReconstructedParticle*>(MCP_Jet->getElementAt(vect1[Gselect_combi][i]));
                ReconstructedParticleVec components = jet->getParticles();
                int ncomps = components.size();
                for(int j = 0; j<ncomps; j++){
                    ReconstructedParticle* compPar = components.at(j);
                    for(int k = 0; k < NLink; k++){
                        LCRelation* a_link = dynamic_cast<LCRelation*>(col_Rela->getElementAt(k));
                        if(a_link->getFrom() == compPar){
                            MCParticle* linkto = dynamic_cast<MCParticle*>(a_link->getTo());
//                            cout<<"linkto->getEnergy() : "<<linkto->getEnergy()<<endl;
                            MCParticle* to_parent = linkto;
                            do{to_parent = to_parent->getParents()[0];}
                            while(to_parent->getParents().size() != 0 && to_parent->getPDG() != 94);
                            if(to_parent->getPDG() == 94)
                            {
//                                cout<<"the mass of 94 is "<<to_parent->getMass()<<endl;
                                if(to_parent->getMass() == MCP941->getMass()){En_B1Part += linkto->getEnergy();}
                                else {En_B2Part += linkto->getEnergy();}
                            }
                            else if(to_parent->getParents().size() == 0){En_no94 += linkto->getEnergy();}
                            
                        }
                    }
                }
                if(En_B1Part <= En_B2Part){En_minPart1 += En_B1Part;}
                else {En_minPart1 += En_B2Part;}
                test_EnTotal += En_B1Part;
                test_EnTotal += En_B2Part;
                cout<<"En_B1Part : "<<En_B1Part<<" En_B2Part : "<<En_B2Part<<endl;
            }
            
            cout<<"the second ::::::"<<endl;
            
            En_minPart2 = 0;
            for(int i = 0; i < num_group2; i++){
                cout<<"the second time >>>>>>>>"<<endl;
                float En_B1Part = 0, En_B2Part = 0;
                ReconstructedParticle* jet = dynamic_cast<ReconstructedParticle*>(MCP_Jet->getElementAt(vect1[Gselect_combi + 1][i]));
                ReconstructedParticleVec components = jet->getParticles();
                int ncomps = components.size();
                for(int j = 0; j<ncomps; j++){
                    ReconstructedParticle* compPar = components.at(j);
                    for(int k = 0; k < NLink; k++){
                        LCRelation* a_link = dynamic_cast<LCRelation*>(col_Rela->getElementAt(k));
                        if(a_link->getFrom() == compPar){
                            MCParticle* linkto = dynamic_cast<MCParticle*>(a_link->getTo());
//                            cout<<"linkto->getEnergy() : "<<linkto->getEnergy()<<endl;
                            MCParticle* to_parent = linkto;
                            do{to_parent = to_parent->getParents()[0];}
                            while(to_parent->getParents().size() != 0 && to_parent->getPDG() != 94);
                            if(to_parent->getPDG() == 94)
                            {
//                                cout<<"the mass of 94 is "<<to_parent->getMass()<<endl;
                                if(to_parent->getMass() == MCP941->getMass()){En_B1Part += linkto->getEnergy();}
                                else {En_B2Part += linkto->getEnergy();}
                            }
                            else if(to_parent->getParents().size() == 0){En_no94 += linkto->getEnergy();}
                        }
                    }
                }
                if(En_B1Part <= En_B2Part){En_minPart2 += En_B1Part;}
                else {En_minPart2 += En_B2Part;}
                test_EnTotal += En_B1Part;
                test_EnTotal += En_B2Part;
                cout<<"En_B1Part : "<<En_B1Part<<" En_B2Part : "<<En_B2Part<<endl;

            }
            
            cout<<"mass_B1 : "<<mass_B1<<" mass_B2: "<<mass_B2<<endl;
            cout<<"count_boson : En_no94 : test_Entotal : En_minPart1 : En_minPart2 "<<count_boson<<" : "<<En_no94<<" : "<<test_EnTotal<<" : "<<En_minPart1<<" : "<<En_minPart2<<endl;
            
            
            
//            cout<<"eventType : "<<eventType<<" GeventType : "<<GeventType<<endl;
//            cout<<"chi2 : "<<GminDif<<endl;
//            cout<<"mass_G1 : G2 : R1 : R2 "<<mass_G1<<" : "<<mass_G2<<" : "<<mass_R1<<" : "<<mass_R2<<endl;
            
            
            
            
            //analysis the most suited grouping method
            

            TLorentzVector TLtest(0,0,0,0);
            misMinEn = 0;
            double En11 = 0, En12 = 0, En21 = 0, En22 = 0;
            for(int j=0; j<vect1[Gselect_combi].size(); j++){
                ReconstructedParticle* jet = dynamic_cast<ReconstructedParticle*>(MCP_Jet->getElementAt(vect1[Gselect_combi][j]));
                ReconstructedParticleVec components = jet->getParticles();
                int ncomps = components.size();
                TLorentzVector TL941(0,0,0,0);   //used to store the TL beloned to each boson of each jet
                TLorentzVector TL942(0,0,0,0);
                for(int icomp = 0; icomp < ncomps; icomp++)
                {
                    ReconstructedParticle* compPar = components.at(icomp);
                    TLorentzVector TLcompPar(compPar->getMomentum(), compPar->getEnergy());
                    for(int i=0; i<NLink; i++)
                    {
                        LCRelation* a_link = dynamic_cast<LCRelation*>(col_Rela->getElementAt(i));
                        if(a_link->getFrom() == compPar)
                        {
                            MCParticle* linkto = dynamic_cast<MCParticle*>(a_link->getTo());
                            MCParticle* a_parent = linkto;
                            do{a_parent = a_parent->getParents()[0];}
                            while(a_parent->getParents().size() != 0 && a_parent->getPDG() != 94);
                            if(a_parent->getPDG() == 94)
                            {
                                double mass_test = a_parent->getMass();
                                if(mass_test == mass_941){TL941 += TLcompPar; En11 += TLcompPar.M();}
                                else if(mass_test != mass_941){TL942 += TLcompPar; En12 += TLcompPar.M();}
                            }
                        }
                    }
                }
                if(TL941.E() <= TL942.E()){misMinEn += TL941.E();}
                else if(TL941.E() > TL942.E()){misMinEn += TL942.E();}
                TLtest += TL941;
                TLtest += TL942;
            }
            for(int j=0; j<vect1[Gselect_combi + 1].size(); j++){
                ReconstructedParticle* jet = dynamic_cast<ReconstructedParticle*>(MCP_Jet->getElementAt(vect1[1+Gselect_combi][j]));
                ReconstructedParticleVec components = jet->getParticles();
                int ncomps = components.size();
                TLorentzVector TL941(0,0,0,0);   //used to store the TL beloned to each boson of each jet
                TLorentzVector TL942(0,0,0,0);
                for(int icomp = 0; icomp < ncomps; icomp++)
                {
                    ReconstructedParticle* compPar = components.at(icomp);
                    TLorentzVector TLcompPar(compPar->getMomentum(), compPar->getEnergy());
                    for(int i=0; i<NLink; i++)
                    {
                        LCRelation* a_link = dynamic_cast<LCRelation*>(col_Rela->getElementAt(i));
                        if(a_link->getFrom() == compPar)
                        {
                            MCParticle* linkto = dynamic_cast<MCParticle*>(a_link->getTo());
                            MCParticle* a_parent = linkto;
                            do{a_parent = a_parent->getParents()[0];}
                            while(a_parent->getParents().size() != 0 && a_parent->getPDG() != 94);
                            if(a_parent->getPDG() == 94)
                            {
                                double mass_test = a_parent->getMass();
                                if(mass_test == mass_941){TL941 += TLcompPar; En21 += TLcompPar.M();}
                                else if(mass_test != mass_941){TL942 += TLcompPar; En22 += TLcompPar.M();}

                            }
                        }
                    }
                }
                if(TL941.E() <= TL942.E()){misMinEn += TL941.E();}
                else if(TL941.E() > TL942.E()){misMinEn += TL942.E();}
                TLtest += TL941;
                TLtest += TL942;
            }
//            cout<<"En11 : En12 : En21 : En22 "<<En11<<" : "<<En12<<" : "<<En21<<" : "<<En22<<endl;
//            cout<<"TLtest.E() : "<<TLtest.E()<<endl;
//            cout<<"misMinEn : "<<misMinEn<<endl;
            

            
            

            }catch (lcio::DataNotAvailableException err) {  }

    }
    
    _outputTree->Fill();
    Num ++;
}



void JetClustering::end()
{
    
    if (_outputTree) {
        
        TFile *tree_file = _outputTree->GetCurrentFile(); //just in case we switched to a new file
        //tree_file->cd();
        tree_file->Write();
        delete tree_file;
        //tree_file->Close();
    }
    
}


/*
int main(){
	int n;
	cout << "Please input the number of array : ";
	cin >> n;
	vector<int> a(n);
	cout << "Please input the value of array : " << endl;
	for(int i=0; i<n; i++)
		cin >> a[i];
	vector<vector<int> > res;
	res = Combine(n, a);
	cout << "size : " << res.size() << endl;
	for(int i=0; i<res.size(); i++){
		cout << i+1 << " -- ";
		for(int j=0; j<res[0].size(); j++)
			cout << res[i][j] << " ";
		cout <<endl;
	}
	return 0;
}
*/


