#include <MCPJetClustering.hh>
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

MCPJetClustering a_MCPJetClustering_instance;

MCPJetClustering::MCPJetClustering()
	: Processor("MCPJetClustering"),
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

struct CRP {
    bool operator()(ReconstructedParticle *p1, ReconstructedParticle *p2) const
    {
        return p1->getEnergy() > p2->getEnergy();
    }
};

//std::sort(MCP_PFO_excludeISR.begin(), MCP_PFO_excludeISR.end(), CMP());

struct CMP {
    bool operator()(MCParticle *p1, MCParticle *p2) const
    {
        return p1->getEnergy() > p2->getEnergy();
    }
};

void MCPJetClustering::init() {

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
    _outputTree->Branch("select_j1j2", &select_j1j2, "select_j1j2/F");
    _outputTree->Branch("select_j3j4", &select_j3j4, "select_j3j4/F");
    _outputTree->Branch("eventType", &eventType, "eventType/I");
    _outputTree->Branch("ratio_1", &ratio_1, "ratio_1/F");
    _outputTree->Branch("ratio_2", &ratio_2, "ratio_2/F");
    _outputTree->Branch("En_WP", &En_WP, "En_WP/F");
    _outputTree->Branch("En_WM", &En_WM, "En_WM/F");
    _outputTree->Branch("angle_q1q2", &angle_q1q2, "angle_q1q2/F");
    _outputTree->Branch("angle_q3q4", &angle_q3q4, "angle_q3q4/F");
    _outputTree->Branch("angle_j1j2", &angle_j1j2, "angle_j1j2/F");
    _outputTree->Branch("angle_j3j4", &angle_j3j4, "angle_j3j4/F");
    _outputTree->Branch("qselect_q1q2", &qselect_q1q2, "qselect_q1q2/F");
    _outputTree->Branch("qselect_q3q4", &qselect_q3q4, "qselect_q3q4/F");
    _outputTree->Branch("qeventType", &qeventType, "qeventType/I");
    _outputTree->Branch("ISRPt", &ISRPt, "ISRPt/F");
    _outputTree->Branch("NuPt", &NuPt, "NuPt/F");
    _outputTree->Branch("newratio_1", &newratio_1, "newratio_1/F");
    _outputTree->Branch("newratio_2", &newratio_2, "newratio_2/F");
    _outputTree->Branch("angle_q1q3", &angle_q1q3, "angle_q1q3/F");
    _outputTree->Branch("angle_q1q4", &angle_q1q4, "angle_q1q4/F");
    _outputTree->Branch("angle_q2q3", &angle_q2q3, "angle_q2q3/F");
    _outputTree->Branch("angle_q2q4", &angle_q2q4, "angle_q2q4/F");
    _outputTree->Branch("ratio_WP_1", &ratio_WP_1, "ratio_WP_1/F");
    _outputTree->Branch("ratio_WP_2", &ratio_WP_2, "ratio_WP_2/F");
    _outputTree->Branch("ratio_WM_1", &ratio_WM_1, "ratio_WM_1/F");
    _outputTree->Branch("ratio_WM_2", &ratio_WM_2, "ratio_WM_2/F");
    _outputTree->Branch("ratio_Z1_1", &ratio_Z1_1, "ratio_Z1_1/F");
    _outputTree->Branch("ratio_Z1_2", &ratio_Z1_2, "ratio_Z1_2/F");
    _outputTree->Branch("ratio_Z2_1", &ratio_Z2_1, "ratio_Z2_1/F");
    _outputTree->Branch("ratio_Z2_2", &ratio_Z2_2, "ratio_Z2_2/F");
    _outputTree->Branch("mass_pair1", &mass_pair1, "mass_pair1/F");
    _outputTree->Branch("mass_pair2", &mass_pair2, "mass_pair2/F");
    _outputTree->Branch("costheta_1", &costheta_1, "costheta_1/F");
    _outputTree->Branch("costheta_2", &costheta_2, "costheta_2/F");
    _outputTree->Branch("costheta_3", &costheta_3, "costheta_3/F");
    _outputTree->Branch("costheta_4", &costheta_4, "costheta_4/F");
    _outputTree->Branch("num_d", &num_d, "num_d/I");
    _outputTree->Branch("num_u", &num_u, "num_u/I");
    _outputTree->Branch("num_s", &num_s, "num_s/I");
    _outputTree->Branch("num_c", &num_c, "num_c/I");
    _outputTree->Branch("num_b", &num_b, "num_b/I");
    _outputTree->Branch("num_t", &num_t, "num_t/I");
    _outputTree->Branch("ratio_P1", &ratio_P1, "ratio_P1/F");
    _outputTree->Branch("ratio_P2", &ratio_P2, "ratio_P2/F");
	Num = 0;
}

void MCPJetClustering::processEvent( LCEvent * evtP )
{		

	if (evtP) 								
	{
        try{
            
            cout<<"************************************************************************"<<endl;

            std::vector<MCParticle*> ori_quarks;
            ori_quarks.clear();
            LCCollection* col_MCP = evtP->getCollection( "MCParticle" );
            eventNr = evtP->getEventNumber();
            cout<<"eventNr "<<eventNr<<" Num "<<Num<<endl;
            int nMCP = col_MCP->getNumberOfElements();
            float want_en = 0;
            TLorentzVector want_TL(0,0,0,0);
            for(int i=0; i<nMCP; i++)
            {
                MCParticle* a_MCP = dynamic_cast<MCParticle*>(col_MCP->getElementAt(i));
                int NParents = a_MCP->getParents().size();
                int NDaughters = a_MCP->getDaughters().size();
                int PDG = a_MCP->getPDG();
                TLorentzVector temp_ori(a_MCP->getMomentum()[0],a_MCP->getMomentum()[1],a_MCP->getMomentum()[2],a_MCP->getEnergy());
                if(a_MCP->getGeneratorStatus() == 1){want_TL += temp_ori; cout<<"stable in generator, PDG is "<<PDG<<endl;}
                if(abs(PDG) == 25 && NDaughters == 2)
                {
                    MCParticle* Dau_1 = a_MCP->getDaughters()[0]; MCParticle* Dau_2 = a_MCP->getDaughters()[1]; ori_quarks.push_back(Dau_1); ori_quarks.push_back(Dau_2);
                }
                if(NParents == 0 && (abs(PDG) == 1 || abs(PDG) == 2 || abs(PDG) == 3 || abs(PDG) == 4 || abs(PDG) == 5 || abs(PDG) == 6)){ori_quarks.push_back(a_MCP);}
            }
            cout<<"stable in generator, the energy is "<<want_TL.E()<<endl;
            cout<<"test 1mmmmmmmm"<<endl;
            cout<<"ori_quarks.size() "<<ori_quarks.size()<<endl;
            if(ori_quarks.size() == 5){for(int i = 0; i< 5; i++){cout<<ori_quarks.at(i)->getPDG()<<endl;}}
            cout<<"test 2mmmmmmmmm"<<endl;
            float wmass = 80, zmass = 91, hmass = 125;

            MCParticle* quark[4];
            std::vector<TLorentzVector> TL_quark;
            cout<<"test 3mmmmmmm"<<endl;
            std::vector<TVector3> TV_quark;
            cout<<"test 4mmmmmmm"<<endl;
            for(int i=0; i< 4; i++){TLorentzVector temp(0,0,0,0); TL_quark.push_back(temp);}
            num_u = 0, num_d = 0, num_s = 0, num_c = 0, num_b = 0, num_t = 0;
            for(int i=0; i<ori_quarks.size(); i++)
            {
                quark[i] = ori_quarks.at(i);
		if(abs(quark[i]->getPDG()) == 1){num_d += 1;}
		else if(abs(quark[i]->getPDG()) == 2){num_u += 1;}
		else if(abs(quark[i]->getPDG()) == 3){num_s += 1;}
		else if(abs(quark[i]->getPDG()) == 4){num_c += 1;}
		else if(abs(quark[i]->getPDG()) == 5){num_b += 1;}
		else if(abs(quark[i]->getPDG()) == 6){num_t += 1;}
                TL_quark.at(i).SetPxPyPzE(quark[i]->getMomentum()[0], quark[i]->getMomentum()[1], quark[i]->getMomentum()[2], quark[i]->getEnergy());
//                TV_quark.at(i).SetPxPyPz(quark[i]->getMomentum()[0], quark[i]->getMomentum()[1], quark[i]->getMomentum()[2]);
            }
	    cout<<"num_d : num_u : num_s : num_c : num_b : num_t "<<num_d<<" : "<<num_u<<" : "<<num_s<<" : "<<num_c<<" : "<<num_b<<" : "<<num_t<<endl;
            cout<<"test 5mmmmmm"<<endl;
            int select_pair = 999;
            float charge_pair1 = 999, charge_pair2 = 999;
            mass_pair1 = 999, mass_pair2 = 999;
            for(int i=1; i<4; i++)
            {
                if(abs(quark[0]->getCharge() + quark[i]->getCharge()) == 1){select_pair = i;}
            }
            if(select_pair == 1)
            {
                charge_pair1 = quark[0]->getCharge() + quark[1]->getCharge();
                if(charge_pair1 == 1)
                {
                    mass_pair1 = (TL_quark.at(0) + TL_quark.at(1)).M();
                    mass_pair2 = (TL_quark.at(2) + TL_quark.at(3)).M();
                }
                else if(charge_pair1 == -1)
                {
                    mass_pair2 = (TL_quark.at(0) + TL_quark.at(1)).M();
                    mass_pair1 = (TL_quark.at(2) + TL_quark.at(3)).M();
                }
                charge_pair2 = quark[2]->getCharge() + quark[3]->getCharge();
            }
            else if(select_pair == 2)
            {
                charge_pair1 = quark[0]->getCharge() + quark[2]->getCharge();
                if(charge_pair1 == 1)
                {
                    mass_pair1 = (TL_quark.at(0) + TL_quark.at(2)).M();
                    mass_pair2 = (TL_quark.at(1) + TL_quark.at(3)).M();
                }
                else if(charge_pair1 == -1)
                {
                    mass_pair2 = (TL_quark.at(0) + TL_quark.at(2)).M();
                    mass_pair1 = (TL_quark.at(1) + TL_quark.at(3)).M();
                }
                charge_pair2 = quark[1]->getCharge() + quark[3]->getCharge();
            }
            else if(select_pair == 3)
            {
                charge_pair1 = quark[0]->getCharge() + quark[3]->getCharge();
                if(charge_pair1 == 1)
                {
                    mass_pair1 = (TL_quark.at(0) + TL_quark.at(3)).M();
                    mass_pair2 = (TL_quark.at(1) + TL_quark.at(2)).M();
                }
                else if(charge_pair1 == -1)
                {
                    mass_pair2 = (TL_quark.at(0) + TL_quark.at(3)).M();
                    mass_pair1 = (TL_quark.at(1) + TL_quark.at(2)).M();
                }
                charge_pair2 = quark[1]->getCharge() + quark[2]->getCharge();
            }
            
            cout<<"charge_pair1 : charge_pair2 : mass_pair1 : mass_pair2 "<<charge_pair1<<" : "<<charge_pair2<<" : "<<mass_pair1<<" : "<<mass_pair2<<endl;
            
            float mass_q1q2 = (TL_quark.at(0) + TL_quark.at(1)).M();
            float mass_q3q4 = (TL_quark.at(2) + TL_quark.at(3)).M();
            float mass_q1q3 = (TL_quark.at(0) + TL_quark.at(2)).M();
            float mass_q2q4 = (TL_quark.at(1) + TL_quark.at(3)).M();
            float mass_q1q4 = (TL_quark.at(0) + TL_quark.at(3)).M();
            float mass_q2q3 = (TL_quark.at(1) + TL_quark.at(2)).M();
            
            float combi_quark[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
            combi_quark[0] = (mass_q1q2 - wmass)*(mass_q1q2 - wmass) + (mass_q3q4 - wmass)*(mass_q3q4 - wmass);
            combi_quark[1] = (mass_q1q3 - wmass)*(mass_q1q3 - wmass) + (mass_q2q4 - wmass)*(mass_q2q4 - wmass);
            combi_quark[2] = (mass_q1q4 - wmass)*(mass_q1q4 - wmass) + (mass_q2q3 - wmass)*(mass_q2q3 - wmass);
            combi_quark[3] = (mass_q1q2 - zmass)*(mass_q1q2 - zmass) + (mass_q3q4 - zmass)*(mass_q3q4 - zmass);
            combi_quark[4] = (mass_q1q3 - zmass)*(mass_q1q3 - zmass) + (mass_q2q4 - zmass)*(mass_q2q4 - zmass);
            combi_quark[5] = (mass_q1q4 - zmass)*(mass_q1q4 - zmass) + (mass_q2q3 - zmass)*(mass_q2q3 - zmass);
            combi_quark[6] = (mass_q1q2 - hmass)*(mass_q1q2 - hmass) + (mass_q3q4 - zmass)*(mass_q3q4 - zmass);
            combi_quark[7] = (mass_q1q3 - hmass)*(mass_q1q3 - hmass) + (mass_q2q4 - zmass)*(mass_q2q4 - zmass);
            combi_quark[8] = (mass_q1q4 - hmass)*(mass_q1q4 - hmass) + (mass_q2q3 - zmass)*(mass_q2q3 - zmass);
            combi_quark[9] = (mass_q1q2 - zmass)*(mass_q1q2 - zmass) + (mass_q3q4 - hmass)*(mass_q3q4 - hmass);
            combi_quark[10] = (mass_q1q3 - zmass)*(mass_q1q3 - zmass) + (mass_q2q4 - hmass)*(mass_q2q4 - hmass);
            combi_quark[11] = (mass_q1q4 - zmass)*(mass_q1q4 - zmass) + (mass_q2q3 - hmass)*(mass_q2q3 - hmass);
            cout<<"mass_q1q2 "<<mass_q1q2<<endl;
            float qminimal = 1E9;
            int qMinIndex = -10;
            int qcombiFlag = 999;
            qselect_q1q2 = 0; qselect_q3q4 = 0;
            for(int i=0; i<12; i++)
            {
                if(combi_quark[i] < qminimal )
                {
                    qminimal = combi_quark[i];
                    qMinIndex = i;
                }
            }
            qeventType = 999;
            qeventType = int(qMinIndex/3); //used to store the type of event
            qcombiFlag = qMinIndex%3;    //used to record the combine method
            
            if(qcombiFlag == 0){qselect_q1q2 = mass_q1q2; qselect_q3q4 = mass_q3q4; cout<<"qcombiFlag : "<<qcombiFlag<<endl;}
            else if(qcombiFlag == 1){qselect_q1q2 = mass_q1q3; qselect_q3q4 = mass_q2q4; cout<<"qcombiFlag : "<<qcombiFlag<<endl;}
            else if(qcombiFlag == 2){qselect_q1q2 = mass_q1q4; qselect_q3q4 = mass_q2q3; cout<<"qcombiFlag : "<<qcombiFlag<<endl;}
            cout<<"qeventType : "<<qeventType<<"  qselect_q1q2 : "<<qselect_q1q2<<"  qselect_q3q4 : "<<qselect_q3q4<<" qminimal : "<<qminimal<<endl;
            
            
            LCCollection* col_MCPRECO = evtP->getCollection( "UsedGenJetParticle" );
            int nMCPRECO = col_MCPRECO->getNumberOfElements();
            cout<<"UsedGenJetParticle, the number of particles used in GenJet is "<<nMCPRECO<<endl;
            
            
    
    //for MCPFastJets
            float Mass_j1j2 = 0, Mass_j3j4 = 0, Mass_j1j3 = 0, Mass_j2j4 = 0, Mass_j1j4 = 0,  Mass_j2j3 = 0;
            float Mass_WP = 0, Mass_WM = 0;
            costheta_1 = 999, costheta_2 = 999, costheta_3 = 999, costheta_4 = 999;
            En_WP = 999; En_WM = 999;
            angle_q1q2 = 0; angle_q3q4 = 0; angle_j1j2 = 0; angle_j3j4 = 0;
            TLorentzVector combi_12(0,0,0,0);
            TLorentzVector combi_34(0,0,0,0);
            TLorentzVector combi_13(0,0,0,0);
            TLorentzVector combi_24(0,0,0,0);
            TLorentzVector combi_14(0,0,0,0);
            TLorentzVector combi_23(0,0,0,0);
            TLorentzVector TL_WP(0,0,0,0);
            TLorentzVector TL_WM(0,0,0,0);
            TVector3 TV_WP(0,0,0), TV_WM(0,0,0);
            LCCollection* col_MCJet = evtP->getCollection( "MCPFastJets" );
            int nGetJet = col_MCJet->getNumberOfElements();
            cout<<" the number of GenJet is "<<nGetJet<<endl;
            ReconstructedParticle* a_GenJet_1 = dynamic_cast<ReconstructedParticle*>(col_MCJet->getElementAt(0));
            ReconstructedParticle* a_GenJet_2 = dynamic_cast<ReconstructedParticle*>(col_MCJet->getElementAt(1));
            ReconstructedParticle* a_GenJet_3 = dynamic_cast<ReconstructedParticle*>(col_MCJet->getElementAt(2));
            ReconstructedParticle* a_GenJet_4 = dynamic_cast<ReconstructedParticle*>(col_MCJet->getElementAt(3));
            costheta_1 = (a_GenJet_1->getMomentum()[2])/(a_GenJet_1->getEnergy());
            costheta_2 = (a_GenJet_2->getMomentum()[2])/(a_GenJet_2->getEnergy());
            costheta_3 = (a_GenJet_3->getMomentum()[2])/(a_GenJet_3->getEnergy());
            costheta_4 = (a_GenJet_4->getMomentum()[2])/(a_GenJet_4->getEnergy());
            cout<<"costheta_1 : costheta_2 : costheta_3 : costheta_4"<<costheta_1<<" : "<<costheta_2<<" : "<<costheta_3<<" : "<<costheta_4<<endl;
            
            
            TLorentzVector TL_GenJet_1(a_GenJet_1->getMomentum(), a_GenJet_1->getEnergy());
            TLorentzVector TL_GenJet_2(a_GenJet_2->getMomentum(), a_GenJet_2->getEnergy());
            TLorentzVector TL_GenJet_3(a_GenJet_3->getMomentum(), a_GenJet_3->getEnergy());
            TLorentzVector TL_GenJet_4(a_GenJet_4->getMomentum(), a_GenJet_4->getEnergy());
            TVector3 TV_GenJet_1 = a_GenJet_1->getMomentum();
            TVector3 TV_GenJet_2 = a_GenJet_2->getMomentum();
            TVector3 TV_GenJet_3 = a_GenJet_3->getMomentum();
            TVector3 TV_GenJet_4 = a_GenJet_4->getMomentum();
            cout<<"the total energy in MCPFastJets is "<<(TL_GenJet_1 + TL_GenJet_2 + TL_GenJet_3 + TL_GenJet_4).E()<<endl;
            combi_12 = TL_GenJet_1 + TL_GenJet_2;   Mass_j1j2 = combi_12.M();
            combi_34 = TL_GenJet_3 + TL_GenJet_4;   Mass_j3j4 = combi_34.M();
            combi_13 = TL_GenJet_1 + TL_GenJet_3;   Mass_j1j3 = combi_13.M();
            combi_24 = TL_GenJet_2 + TL_GenJet_4;   Mass_j2j4 = combi_24.M();
            combi_14 = TL_GenJet_1 + TL_GenJet_4;   Mass_j1j4 = combi_14.M();
            combi_23 = TL_GenJet_2 + TL_GenJet_3;   Mass_j2j3 = combi_23.M();

            float combi[12] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
            combi[0] = (Mass_j1j2 - wmass)*(Mass_j1j2 - wmass) + (Mass_j3j4 - wmass)*(Mass_j3j4 - wmass);
            combi[1] = (Mass_j1j3 - wmass)*(Mass_j1j3 - wmass) + (Mass_j2j4 - wmass)*(Mass_j2j4 - wmass);
            combi[2] = (Mass_j1j4 - wmass)*(Mass_j1j4 - wmass) + (Mass_j2j3 - wmass)*(Mass_j2j3 - wmass);
            combi[3] = (Mass_j1j2 - zmass)*(Mass_j1j2 - zmass) + (Mass_j3j4 - zmass)*(Mass_j3j4 - zmass);
            combi[4] = (Mass_j1j3 - zmass)*(Mass_j1j3 - zmass) + (Mass_j2j4 - zmass)*(Mass_j2j4 - zmass);
            combi[5] = (Mass_j1j4 - zmass)*(Mass_j1j4 - zmass) + (Mass_j2j3 - zmass)*(Mass_j2j3 - zmass);
            combi[6] = (Mass_j1j2 - zmass)*(Mass_j1j2 - zmass) + (Mass_j3j4 - hmass)*(Mass_j3j4 - hmass);
            combi[7] = (Mass_j1j3 - zmass)*(Mass_j1j3 - zmass) + (Mass_j2j4 - hmass)*(Mass_j2j4 - hmass);
            combi[8] = (Mass_j1j4 - zmass)*(Mass_j1j4 - zmass) + (Mass_j2j3 - hmass)*(Mass_j2j3 - hmass);
            combi[9] = (Mass_j1j2 - hmass)*(Mass_j1j2 - hmass) + (Mass_j3j4 - zmass)*(Mass_j3j4 - zmass);
            combi[10] = (Mass_j1j3 - hmass)*(Mass_j1j3 - hmass) + (Mass_j2j4 - zmass)*(Mass_j2j4 - zmass);
            combi[11] = (Mass_j1j4 - hmass)*(Mass_j1j4 - hmass) + (Mass_j2j3 - zmass)*(Mass_j2j3 - zmass);
            
            float minimal = 1E9;
            int MinIndex = -10;
            int combiFlag = 999;
            select_j1j2 = 0, select_j3j4 = 0;
            for(int i=0; i<12; i++)
            {
                if(combi[i] < minimal )
                {
                    minimal = combi[i];
                    MinIndex = i;
                }
            }
            eventType = 999;
            eventType = int(MinIndex/3); //used to store the type of event
            combiFlag = MinIndex%3;    //used to record the combine method
            
            if(combiFlag == 0){select_j1j2 = Mass_j1j2; select_j3j4 = Mass_j3j4; cout<<"combiFlag : "<<combiFlag<<endl;}
            else if(combiFlag == 1){select_j1j2 = Mass_j1j3; select_j3j4 = Mass_j2j4; cout<<"combiFlag : "<<combiFlag<<endl;}
            else if(combiFlag == 2){select_j1j2 = Mass_j1j4; select_j3j4 = Mass_j2j3; cout<<"combiFlag : "<<combiFlag<<endl;}
            
            cout<<"minimal "<<minimal<<"   select_j1j2 : select_j3j4 "<<select_j1j2<<" : "<<select_j3j4<<endl;
            cout<<"0 : w    1 : z    2 : zh    3  : hz    "<<endl;
            cout<<" eventType :  "<<eventType<<endl;
            
            ReconstructedParticleVec components_1 = a_GenJet_1->getParticles(); int ncomp_1 = components_1.size();
            ReconstructedParticleVec components_2 = a_GenJet_2->getParticles(); int ncomp_2 = components_2.size();
            ReconstructedParticleVec components_3 = a_GenJet_3->getParticles(); int ncomp_3 = components_3.size();
            ReconstructedParticleVec components_4 = a_GenJet_4->getParticles(); int ncomp_4 = components_4.size();
            cout<<"MCPFastJets, the number used in GenJet is "<<ncomp_4+ncomp_1+ncomp_2+ncomp_3<<endl;
            TVector3 TV_quark1(0,0,0), TV_quark2(0,0,0), TV_quark3(0,0,0), TV_quark4(0,0,0);
            LCCollection* col_Rela = evtP->getCollection( "UsedRelation" ); //LCRelation used to link particles in GenJet with MCParticles
            int NLink = col_Rela->getNumberOfElements();

    //for GenJet1
            map<int, int> mapcharge_1;
            map<float, int> map_z_1;
            float mass_1 = 99999;
                for(int icomp_1=0; icomp_1<ncomp_1; icomp_1++)
                {
                    ReconstructedParticle* compPar_1 = components_1.at(icomp_1);
                    for(int i=0; i<NLink; i++)
                    {
                        LCRelation* a_link = dynamic_cast<LCRelation*>(col_Rela->getElementAt(i));
                        if(a_link->getFrom() == compPar_1)
                        {
                            MCParticle* linkto = dynamic_cast<MCParticle*>(a_link->getTo());
                            MCParticle* a_parent = linkto;
                            do{a_parent = a_parent->getParents()[0];}
                            while(a_parent->getParents().size() != 0 && a_parent->getParents()[0]->getPDG() != 94);
                         
                            MCParticle* new_parent = a_parent->getParents()[0];
                            MCParticle* Dau94_1 = new_parent->getDaughters()[0]; MCParticle* Dau94_2 = new_parent->getDaughters()[1];
                            TVector3 TV_Dau94_1 = Dau94_1->getMomentum();  TVector3 TV_Dau94_2 = Dau94_2->getMomentum();
                            TLorentzVector TL_Dau94_1(Dau94_1->getMomentum(), Dau94_1->getEnergy());
                            TLorentzVector TL_Dau94_2(Dau94_2->getMomentum(), Dau94_2->getEnergy());
                            float charge_Dau941 = Dau94_1->getCharge(); float charge_Dau942 = Dau94_2->getCharge();
                            float sum_charge = charge_Dau941 + charge_Dau942;
                            float mass_test = (TL_Dau94_1 + TL_Dau94_2).M();
                            if(icomp_1 == 0){mass_1 = mass_test;}
                            if(mass_test == mass_1){map_z_1[mass_test] += 1; angle_q1q2 = TV_Dau94_1.Angle(TV_Dau94_2); TL_WP = TL_Dau94_1 + TL_Dau94_2; TV_quark1 = TV_Dau94_1; TV_quark2 = TV_Dau94_2;}
                            else if(mass_test != mass_1){map_z_1[mass_test] += 1; angle_q3q4 = TV_Dau94_1.Angle(TV_Dau94_2); TL_WM = TL_Dau94_1 + TL_Dau94_2; TV_quark3 = TV_Dau94_1; TV_quark4 = TV_Dau94_2;}
                            if(sum_charge == 1){mapcharge_1[1] += 1; angle_q1q2 = TV_Dau94_1.Angle(TV_Dau94_2); TL_WP = TL_Dau94_1 + TL_Dau94_2; TV_quark1 = TV_Dau94_1; TV_quark2 = TV_Dau94_2;}
                            else if(sum_charge == -1){mapcharge_1[-1] += 1; angle_q3q4 = TV_Dau94_1.Angle(TV_Dau94_2); TL_WM = TL_Dau94_1 + TL_Dau94_2; TV_quark3 = TV_Dau94_1; TV_quark4 = TV_Dau94_2;}
                                int ParPDG = a_parent->getPDG();
                        }
                    }
                }
            auto it_1 = map_z_1.begin(); int num_1 = it_1->second; float mass_temp_1 = it_1->first;
            cout<<"ncomp_1 "<<ncomp_1<<endl;
            for(; it_1 != map_z_1.end(); it_1 ++){
                if(num_1 < it_1->second){mass_temp_1 = it_1->first; num_1 = it_1->second;}
                cout<<"it_1->second "<<(it_1->second)<<endl; cout<<"it_1->first "<<(it_1->first)<<endl;
                
            }
            cout<<"num_1 : mass_temp_1 "<<num_1<<" : "<<mass_temp_1<<endl;
            auto itcharge_1 = mapcharge_1.begin(); int tmpcharge_1 = itcharge_1->second; int rescharge_1 = itcharge_1->first;
            for(; itcharge_1 != mapcharge_1.end(); itcharge_1++){ if(tmpcharge_1 < itcharge_1->second){tmpcharge_1 = itcharge_1->second; rescharge_1 = itcharge_1->first;} }
            cout<<"the larger "<<mapcharge_1.size()<<" "<<rescharge_1<<" "<<tmpcharge_1<<endl;
            TLorentzVector fault_1(0,0,0,0);
            TLorentzVector yes_1(0,0,0,0);
            TLorentzVector z_fault_1(0,0,0,0);
            TLorentzVector z_yes_1(0,0,0,0);
            for(int icomp_1=0; icomp_1<ncomp_1; icomp_1++)
            {
                ReconstructedParticle* compPar_1 = components_1.at(icomp_1);
                TLorentzVector temp_1(compPar_1->getMomentum(),compPar_1->getEnergy());
                
                for(int i=0; i<NLink; i++)
                {
                    LCRelation* a_link = dynamic_cast<LCRelation*>(col_Rela->getElementAt(i));
                    if(a_link->getFrom() == compPar_1)
                    {
                        MCParticle* linkto = dynamic_cast<MCParticle*>(a_link->getTo());
                        MCParticle* a_parent = linkto;
                        do{a_parent = a_parent->getParents()[0];}
                        while(a_parent->getParents().size() != 0 && a_parent->getParents()[0]->getPDG() != 94);
                      
                        int ParPDG = a_parent->getPDG();
                        MCParticle* new_parent = a_parent->getParents()[0];
                        MCParticle* Dau94_1 = new_parent->getDaughters()[0]; MCParticle* Dau94_2 = new_parent->getDaughters()[1];
                        TLorentzVector TL_Dau94_1(Dau94_1->getMomentum(), Dau94_1->getEnergy());
                        TLorentzVector TL_Dau94_2(Dau94_2->getMomentum(), Dau94_2->getEnergy());
                        float mass_test = (TL_Dau94_1 + TL_Dau94_2).M();
                        if(mass_test == mass_temp_1){z_yes_1 += temp_1;}
                        else if(mass_test != mass_temp_1){z_fault_1 += temp_1;}
                        float charge_Dau941 = Dau94_1->getCharge(); float charge_Dau942 = Dau94_2->getCharge();
                        float sum_charge = charge_Dau941 + charge_Dau942;
                        if(sum_charge != rescharge_1){fault_1 += temp_1;}
                        else if(sum_charge == rescharge_1){yes_1 += temp_1;}
                    }
                }
            }


    //for GenJet2
            map<int, int>mapcharge_2;
            map<float, int>map_z_2;
            float mass_2 = 99999;
                for(int icomp_2=0; icomp_2<ncomp_2; icomp_2++)
                {
                    ReconstructedParticle* compPar_2 = components_2.at(icomp_2);
                    for(int i=0; i<NLink; i++)
                    {
                        LCRelation* a_link = dynamic_cast<LCRelation*>(col_Rela->getElementAt(i));
                        if(a_link->getFrom() == compPar_2)
                        {
                            MCParticle* linkto = dynamic_cast<MCParticle*>(a_link->getTo());
                            MCParticle* a_parent = linkto;
                            do{a_parent = a_parent->getParents()[0];}
                            while(a_parent->getParents().size() != 0 && a_parent->getParents()[0]->getPDG() != 94);
          
                            MCParticle* new_parent = a_parent->getParents()[0];
                            MCParticle* Dau94_1 = new_parent->getDaughters()[0]; MCParticle* Dau94_2 = new_parent->getDaughters()[1];
                            TVector3 TV_Dau94_1 = Dau94_1->getMomentum(); TVector3 TV_Dau94_2 = Dau94_2->getMomentum();
                            TLorentzVector TL_Dau94_1(Dau94_1->getMomentum(), Dau94_1->getEnergy());
                            TLorentzVector TL_Dau94_2(Dau94_2->getMomentum(), Dau94_2->getEnergy());
                            float mass_test = (TL_Dau94_1 + TL_Dau94_2).M();
                            if(icomp_2 == 0){mass_2 = mass_test;}
                            if(mass_test == mass_2){map_z_2[mass_test] += 1; angle_q1q2 = TV_Dau94_1.Angle(TV_Dau94_2); TL_WP = TL_Dau94_1 + TL_Dau94_2; TV_quark1 = TV_Dau94_1; TV_quark2 = TV_Dau94_2;}
                            else if(mass_test != mass_2){map_z_2[mass_test] += 1; angle_q3q4 = TV_Dau94_1.Angle(TV_Dau94_2); TL_WM = TL_Dau94_1 + TL_Dau94_2; TV_quark3 = TV_Dau94_1; TV_quark4 = TV_Dau94_2;}
                            float charge_Dau941 = Dau94_1->getCharge(); float charge_Dau942 = Dau94_2->getCharge();
                            float sum_charge = charge_Dau941 + charge_Dau942;
                            if(sum_charge == 1){mapcharge_2[1] += 1; angle_q1q2 = TV_Dau94_1.Angle(TV_Dau94_2); TL_WP = TL_Dau94_1 + TL_Dau94_2; TV_quark1 = TV_Dau94_1; TV_quark2 = TV_Dau94_2;}
                            else if(sum_charge == -1){mapcharge_2[-1] += 1; angle_q3q4 = TV_Dau94_1.Angle(TV_Dau94_2); TL_WM = TL_Dau94_1 + TL_Dau94_2; TV_quark3 = TV_Dau94_1; TV_quark4 = TV_Dau94_2;}
                                int ParPDG = a_parent->getPDG();
                        }
                    }
                }
            auto it_2 = map_z_2.begin(); int num_2 = it_2->second; float mass_temp_2 = it_2->first;
            cout<<"ncomp_2 "<<ncomp_2<<endl;
            for(; it_2 != map_z_2.end(); it_2 ++){
                if(num_2 < it_2->second){mass_temp_2 = it_2->first; num_2 = it_2->second;}
                cout<<"it_2->second "<<(it_2->second)<<endl; cout<<"it_2->first "<<(it_2->first)<<endl;
                
            }
            cout<<"num_2 : mass_temp_2 "<<num_2<<" : "<<mass_temp_2<<endl;
            
            auto itcharge_2 = mapcharge_2.begin(); int tmpcharge_2 = itcharge_2->second; int rescharge_2 = itcharge_2->first;
            for(; itcharge_2 != mapcharge_2.end(); itcharge_2++){ if(tmpcharge_2 < itcharge_2->second){tmpcharge_2 = itcharge_2->second; rescharge_2 = itcharge_2->first;} }
            cout<<"the larger "<<mapcharge_2.size()<<" "<<rescharge_2<<" "<<tmpcharge_2<<endl;
            TLorentzVector fault_2(0,0,0,0);
            TLorentzVector yes_2(0,0,0,0);
            TLorentzVector z_fault_2(0,0,0,0);
            TLorentzVector z_yes_2(0,0,0,0);
            for(int icomp_2=0; icomp_2<ncomp_2; icomp_2++)
            {
                ReconstructedParticle* compPar_2 = components_2.at(icomp_2);
                TLorentzVector temp_2(compPar_2->getMomentum(),compPar_2->getEnergy());
                
                for(int i=0; i<NLink; i++)
                {
                    LCRelation* a_link = dynamic_cast<LCRelation*>(col_Rela->getElementAt(i));
                    if(a_link->getFrom() == compPar_2)
                    {
                        MCParticle* linkto = dynamic_cast<MCParticle*>(a_link->getTo());
                        MCParticle* a_parent = linkto;
                        do{a_parent = a_parent->getParents()[0];}
                        while(a_parent->getParents().size() != 0 && a_parent->getParents()[0]->getPDG() != 94);
                            int ParPDG = a_parent->getPDG();
                            MCParticle* new_parent = a_parent->getParents()[0];
                            MCParticle* Dau94_1 = new_parent->getDaughters()[0]; MCParticle* Dau94_2 = new_parent->getDaughters()[1];
                        TLorentzVector TL_Dau94_1(Dau94_1->getMomentum(), Dau94_1->getEnergy());
                        TLorentzVector TL_Dau94_2(Dau94_2->getMomentum(), Dau94_2->getEnergy());
                        float mass_test = (TL_Dau94_2 + TL_Dau94_1).M();
                        if(mass_test == mass_temp_2){z_yes_2 += temp_2;}
                        else if(mass_test != mass_temp_2){z_fault_2 += temp_2;}
                            float charge_Dau941 = Dau94_1->getCharge(); float charge_Dau942 = Dau94_2->getCharge();
                            float sum_charge = charge_Dau941 + charge_Dau942;
                            if(sum_charge != rescharge_2){fault_2 += temp_2;}
                            if(sum_charge == rescharge_2){yes_2 += temp_2;}
                    }
                }
            }

            
        
    //for GenJet3
            map<int, int>mapcharge_3;
            map<float, int>map_z_3;
            float mass_3 = 99999;
            for(int icomp_3=0; icomp_3<ncomp_3; icomp_3++)
            {
                ReconstructedParticle* compPar_3 = components_3.at(icomp_3);
                for(int i=0; i<NLink; i++)
                {
                    LCRelation* a_link = dynamic_cast<LCRelation*>(col_Rela->getElementAt(i));
                    if(a_link->getFrom() == compPar_3)
                    {
                        MCParticle* linkto = dynamic_cast<MCParticle*>(a_link->getTo());
                        MCParticle* a_parent = linkto;
                        do{a_parent = a_parent->getParents()[0];}
                        while(a_parent->getParents().size() != 0 && a_parent->getParents()[0]->getPDG() != 94);
                            int ParPDG = a_parent->getPDG();
                        MCParticle* new_parent = a_parent->getParents()[0];
                        MCParticle* Dau94_1 = new_parent->getDaughters()[0]; MCParticle* Dau94_2 = new_parent->getDaughters()[1];
                        TVector3 TV_Dau94_1 = Dau94_1->getMomentum(); TVector3 TV_Dau94_2 = Dau94_2->getMomentum();
                        TLorentzVector TL_Dau94_1(Dau94_1->getMomentum(), Dau94_1->getEnergy());
                        TLorentzVector TL_Dau94_2(Dau94_2->getMomentum(), Dau94_2->getEnergy());
                        float mass_test = (TL_Dau94_1 + TL_Dau94_2).M();
                        if(icomp_3 == 0){mass_3 = mass_test;}
                        if(mass_test == mass_3){map_z_3[mass_test] += 1; angle_q1q2 = TV_Dau94_1.Angle(TV_Dau94_2); TL_WP = TL_Dau94_1 + TL_Dau94_2; TV_quark1 = TV_Dau94_1; TV_quark2 = TV_Dau94_2;}
                        else if(mass_test != mass_3){map_z_3[mass_test] += 1; angle_q3q4 = TV_Dau94_1.Angle(TV_Dau94_2); TL_WM = TL_Dau94_1 + TL_Dau94_2; TV_quark3 = TV_Dau94_1; TV_quark4 = TV_Dau94_2;}
                        float charge_Dau941 = Dau94_1->getCharge(); float charge_Dau942 = Dau94_2->getCharge();
                        float sum_charge = charge_Dau941 + charge_Dau942;
                        if(sum_charge == 1){mapcharge_3[1] += 1; angle_q1q2 = TV_Dau94_1.Angle(TV_Dau94_2); TL_WP = TL_Dau94_1 + TL_Dau94_2; TV_quark1 = TV_Dau94_1; TV_quark2 = TV_Dau94_2;}
                        else if(sum_charge == -1){mapcharge_3[-1] += 1; angle_q3q4 = TV_Dau94_1.Angle(TV_Dau94_2); TL_WM = TL_Dau94_1 + TL_Dau94_2; TV_quark3 = TV_Dau94_1; TV_quark4 = TV_Dau94_2;}
                    }
                }
            }
            auto it_3 = map_z_3.begin(); int num_3 = it_3->second; float mass_temp_3 = it_3->first;
            cout<<"ncomp_3 "<<ncomp_3<<endl;
            for(; it_3 != map_z_3.end(); it_3 ++){
                if(num_3 < it_3->second){mass_temp_3 = it_3->first; num_3 = it_3->second;}
                cout<<"it_3->second "<<(it_3->second)<<endl; cout<<"it_3->first "<<(it_3->first)<<endl;
                
            }
            cout<<"num_3 : mass_temp_3 "<<num_3<<" : "<<mass_temp_3<<endl;
            
            auto itcharge_3 = mapcharge_3.begin(); int tmpcharge_3 = itcharge_3->second; int rescharge_3 = itcharge_3->first;
            for(; itcharge_3 != mapcharge_3.end(); itcharge_3++){ if(tmpcharge_3 < itcharge_3->second){tmpcharge_3 = itcharge_3->second; rescharge_3 = itcharge_3->first;} }
            cout<<"the larger "<<mapcharge_3.size()<<" "<<rescharge_3<<" "<<tmpcharge_3<<endl;
            TLorentzVector fault_3(0,0,0,0);
            TLorentzVector yes_3(0,0,0,0);
            TLorentzVector z_fault_3(0,0,0,0);
            TLorentzVector z_yes_3(0,0,0,0);
            
            for(int icomp_3=0; icomp_3<ncomp_3; icomp_3++)
            {
                ReconstructedParticle* compPar_3 = components_3.at(icomp_3);
                TLorentzVector temp_3(compPar_3->getMomentum(),compPar_3->getEnergy());
                
                for(int i=0; i<NLink; i++)
                {
                    LCRelation* a_link = dynamic_cast<LCRelation*>(col_Rela->getElementAt(i));
                    if(a_link->getFrom() == compPar_3)
                    {
                        MCParticle* linkto = dynamic_cast<MCParticle*>(a_link->getTo());
                        MCParticle* a_parent = linkto;
                        do{a_parent = a_parent->getParents()[0];}
                        while(a_parent->getParents().size() != 0 && a_parent->getParents()[0]->getPDG() != 94);
                        int ParPDG = a_parent->getPDG();
                        MCParticle* new_parent = a_parent->getParents()[0];
                        MCParticle* Dau94_1 = new_parent->getDaughters()[0]; MCParticle* Dau94_2 = new_parent->getDaughters()[1];
                        TLorentzVector TL_Dau94_1(Dau94_1->getMomentum(), Dau94_1->getEnergy());
                        TLorentzVector TL_Dau94_2(Dau94_2->getMomentum(), Dau94_2->getEnergy());
                        float mass_test = (TL_Dau94_2 + TL_Dau94_1).M();
                        if(mass_test == mass_temp_3){z_yes_3 += temp_3;}
                        else if(mass_test != mass_temp_3){z_fault_3 += temp_3;}
                        float charge_Dau941 = Dau94_1->getCharge(); float charge_Dau942 = Dau94_2->getCharge();
                        float sum_charge = charge_Dau941 + charge_Dau942;
                        if(sum_charge != rescharge_3){fault_3 += temp_3;}
                        if(sum_charge == rescharge_3){yes_3 += temp_3;}
                    }
                }
            }

    //for GenJet4
            map<int, int>mapcharge_4;
            map<float, int>map_z_4;
            float mass_4 = 99999;
            cout<<"ncomp_4 is "<<ncomp_4<<endl;
            int count_GenJet4 = 0;
            for(int icomp_4=0; icomp_4<ncomp_4; icomp_4++)
            {
                ReconstructedParticle* compPar_4 = components_4.at(icomp_4);
                for(int i=0; i<NLink; i++)
                {
                    LCRelation* a_link = dynamic_cast<LCRelation*>(col_Rela->getElementAt(i));
                    if(a_link->getFrom() == compPar_4)
                    {
                        MCParticle* linkto = dynamic_cast<MCParticle*>(a_link->getTo());
                        MCParticle* a_parent = linkto;
                        do{a_parent = a_parent->getParents()[0];}
                        while(a_parent->getParents()[0]->getPDG() != 94 && a_parent->getParents().size() != 0);
                        MCParticle* new_parent = a_parent->getParents()[0];
                        MCParticle* Dau94_1 = new_parent->getDaughters()[0]; MCParticle* Dau94_2 = new_parent->getDaughters()[1];
                        TVector3 TV_Dau94_1 = Dau94_1->getMomentum(); TVector3 TV_Dau94_2 = Dau94_2->getMomentum();
                        TLorentzVector TL_Dau94_1(Dau94_1->getMomentum(), Dau94_1->getEnergy());
                        TLorentzVector TL_Dau94_2(Dau94_2->getMomentum(), Dau94_2->getEnergy());
                        float mass_test = (TL_Dau94_1 + TL_Dau94_2).M();
                        if(icomp_4 == 0){mass_4 = mass_test;}
                        if(mass_test == mass_4){map_z_4[mass_test] += 1; angle_q1q2 = TV_Dau94_1.Angle(TV_Dau94_2); TL_WP = TL_Dau94_1 + TL_Dau94_2; TV_quark1 = TV_Dau94_1; TV_quark2 = TV_Dau94_2;}
                        else if(mass_test != mass_4){map_z_4[mass_test] += 1; angle_q3q4 = TV_Dau94_1.Angle(TV_Dau94_2); TL_WM = TL_Dau94_1 + TL_Dau94_2; TV_quark3 = TV_Dau94_1; TV_quark4 = TV_Dau94_2;}
                        float charge_Dau941 = Dau94_1->getCharge(); float charge_Dau942 = Dau94_2->getCharge();
                        float sum_charge = charge_Dau941 + charge_Dau942;
                        if(sum_charge == 1){mapcharge_4[1] += 1; angle_q1q2 = TV_Dau94_1.Angle(TV_Dau94_2); TL_WP = TL_Dau94_2 + TL_Dau94_1; TV_quark1 = TV_Dau94_1; TV_quark2 = TV_Dau94_2;}
                        else if(sum_charge == -1){mapcharge_4[-1] += 1; angle_q3q4 = TV_Dau94_1.Angle(TV_Dau94_2); TL_WM = TL_Dau94_1 + TL_Dau94_2; TV_quark3 = TV_Dau94_1; TV_quark4 = TV_Dau94_2;}
                        int ParPDG = a_parent->getPDG();

                    }
                }
            }
            auto it_4 = map_z_4.begin(); int num_4 = it_4->second; float mass_temp_4 = it_4->first;
            cout<<"ncomp_4 "<<ncomp_4<<endl;
            for(; it_4 != map_z_4.end(); it_4 ++){
                if(num_4 < it_4->second){mass_temp_4 = it_4->first; num_4 = it_4->second;}
                cout<<"it_4->second "<<(it_4->second)<<endl; cout<<"it_4->first "<<(it_4->first)<<endl;
                
            }
            cout<<"num_4 : mass_temp_4 "<<num_4<<" : "<<mass_temp_4<<endl;
            
            auto itcharge_4 = mapcharge_4.begin(); int tmpcharge_4 = itcharge_4->second; int rescharge_4 = itcharge_4->first;
            for(; itcharge_4 != mapcharge_4.end(); itcharge_4++){ if(tmpcharge_4 < itcharge_4->second){tmpcharge_4 = itcharge_4->second; rescharge_4 = itcharge_4->first;} }
            cout<<"the larger "<<mapcharge_4.size()<<" "<<rescharge_4<<" "<<tmpcharge_4<<endl;
            TLorentzVector fault_4(0,0,0,0);
            TLorentzVector yes_4(0,0,0,0);
            TLorentzVector z_fault_4(0,0,0,0);
            TLorentzVector z_yes_4(0,0,0,0);
            
            for(int icomp_4=0; icomp_4<ncomp_4; icomp_4++)
            {
                ReconstructedParticle* compPar_4 = components_4.at(icomp_4);
                TLorentzVector temp_4(compPar_4->getMomentum(),compPar_4->getEnergy());
                
                for(int i=0; i<NLink; i++)
                {
                    LCRelation* a_link = dynamic_cast<LCRelation*>(col_Rela->getElementAt(i));
                    if(a_link->getFrom() == compPar_4)
                    {
                        MCParticle* linkto = dynamic_cast<MCParticle*>(a_link->getTo());
                        MCParticle* a_parent = linkto;
                        do{a_parent = a_parent->getParents()[0];}
                        while(a_parent->getParents().size() != 0 && a_parent->getParents()[0]->getPDG() != 94);
                        int ParPDG = a_parent->getPDG();
                        MCParticle* new_parent = a_parent->getParents()[0];
                        MCParticle* Dau94_1 = new_parent->getDaughters()[0]; MCParticle* Dau94_2 = new_parent->getDaughters()[1];
                        TLorentzVector TL_Dau94_1(Dau94_1->getMomentum(), Dau94_1->getEnergy());
                        TLorentzVector TL_Dau94_2(Dau94_2->getMomentum(), Dau94_2->getEnergy());
                        float mass_test = (TL_Dau94_2 + TL_Dau94_1).M();
                        if(mass_test == mass_temp_4){z_yes_4 += temp_4;}
                        else if(mass_test != mass_temp_4){z_fault_4 += temp_4;}
                        float charge_Dau941 = Dau94_1->getCharge(); float charge_Dau942 = Dau94_2->getCharge();
                        float sum_charge = charge_Dau941 + charge_Dau942;
                        if(sum_charge != rescharge_4){fault_4 += temp_4;}
                        if(sum_charge == rescharge_4){yes_4 += temp_4;}
                    }
                }
            }
            cout<<"(fault_1 + fault_2 + fault_3 + fault_4 + yes_1 + yes_2 + yes_3 + yes_4).E()"<<(fault_1 + fault_2 + fault_3 + fault_4 + yes_1 + yes_2 + yes_3 + yes_4).E()<<endl;
            cout<<"(z_fault_1 + z_fault_2 + z_fault_3 + z_fault_4 + z_yes_1 + z_yes_2 + z_yes_3 + z_yes_4).E()"<<(z_fault_1 + z_fault_2 + z_fault_3 + z_fault_4 + z_yes_1 + z_yes_2 + z_yes_3 + z_yes_4).E()<<endl;

            angle_q1q3 = TV_quark1.Angle(TV_quark3), angle_q1q4 = TV_quark1.Angle(TV_quark4), angle_q2q3 = TV_quark2.Angle(TV_quark3), angle_q2q4 = TV_quark2.Angle(TV_quark4);
            cout<<"angle_q1q3 : angle_q1q4  "<<angle_q1q3<<" : "<<angle_q1q4<<"           angle_q2q3 : angle_q2q4 "<<angle_q2q3<<" : "<<angle_q2q4<<endl;
            ratio_1 = 999, ratio_2 = 999;
            newratio_1 = 0, newratio_2 = 0;
            ratio_WP_1 = 0, ratio_WP_2 = 0, ratio_WM_1 = 0, ratio_WM_2 = 0;
            ratio_Z1_1 = 0, ratio_Z1_2 = 0, ratio_Z2_1 = 0, ratio_Z2_2 = 0;
            ratio_P1 = 0, ratio_P2 = 0;

            ratio_Z1_1 = z_fault_1.E()/(z_fault_1 + z_yes_1).E();
            ratio_Z1_2 = z_fault_2.E()/(z_fault_2 + z_yes_2).E();
            ratio_Z2_1 = z_fault_3.E()/(z_fault_3 + z_yes_3).E();
            ratio_Z2_2 = z_fault_4.E()/(z_fault_4 + z_yes_4).E();
            if(combiFlag == 0)
            {
//                select_j1j2 = Mass_j1j2; select_j3j4 = Mass_j3j4;  ratio_1 = ((fault_1 + fault_2).E())/((yes_1 + yes_2 + fault_1 + fault_2).E());  ratio_2 = ((fault_3 + fault_4).E())/((yes_3 + yes_4 + fault_3 + fault_4).E());  cout<<"combiFlag : "<<combiFlag<<endl;
                select_j1j2 = Mass_j1j2; select_j3j4 = Mass_j3j4;
                
                ratio_P1 = (z_fault_1 + z_fault_2).E()/(z_fault_1 + z_fault_2 + z_yes_1 + z_yes_2).E();
                ratio_P2 = (z_fault_3 + z_fault_4).E()/(z_fault_3 + z_fault_4 + z_yes_3 + z_yes_4).E();
                
                if(rescharge_1 == 1)
                {
                    angle_j1j2 = TV_GenJet_1.Angle(TV_GenJet_2); angle_j3j4 = TV_GenJet_3.Angle(TV_GenJet_4);
                    select_j1j2 = Mass_j1j2; select_j3j4 = Mass_j3j4;
                    ratio_1 = ((fault_1 + fault_2).E())/((yes_1 + yes_2 + fault_1 + fault_2).E());
                    ratio_2 = ((fault_3 + fault_4).E())/((yes_3 + yes_4 + fault_3 + fault_4).E());
                    newratio_1 = ((yes_1 + yes_2).E())/((yes_1 + yes_2 + fault_3 + fault_4).E());
                    newratio_2 = ((yes_3 + yes_4).E())/((yes_3 + yes_4 + fault_1 + fault_2).E());
                    ratio_WP_1 = fault_1.E()/(fault_1 + yes_1).E();
                    ratio_WP_2 = fault_2.E()/(fault_2 + yes_2).E();
                    ratio_WM_1 = fault_3.E()/(fault_3 + yes_3).E();
                    ratio_WM_2 = fault_4.E()/(fault_4 + yes_4).E();
                }
                else if(rescharge_1 == -1)
                {
                    angle_j3j4 = TV_GenJet_1.Angle(TV_GenJet_2); angle_j1j2 = TV_GenJet_3.Angle(TV_GenJet_4);
                    select_j3j4 = Mass_j1j2; select_j1j2 = Mass_j3j4;
                    ratio_1 = ((fault_3 + fault_4).E())/((yes_3 + yes_4 + fault_3 + fault_4).E());
                    ratio_2 = ((fault_1 + fault_2).E())/((yes_1 + yes_2 + fault_1 + fault_2).E());
                    newratio_1 = ((yes_3 + yes_4).E())/((yes_3 + yes_4 + fault_1 +fault_2).E());
                    newratio_2 = ((yes_1 + yes_2).E())/((yes_1 + yes_2 + fault_3 +fault_4).E());
                    ratio_WP_1 = fault_3.E()/(fault_3 + yes_3).E();
                    ratio_WP_2 = fault_4.E()/(fault_4 + yes_4).E();
                    ratio_WM_1 = fault_1.E()/(fault_1 + yes_1).E();
                    ratio_WM_2 = fault_2.E()/(fault_2 + yes_2).E();
                }
            }
            else if(combiFlag == 1)
            {
//                select_j1j2 = Mass_j1j3; select_j3j4 = Mass_j2j4;  ratio_1 = ((fault_1 + fault_3).E())/((yes_1 + yes_3 + fault_1 + fault_3).E()); ratio_2 = ((fault_2 + fault_4).E())/((fault_2 + fault_4 + yes_2 + yes_4).E()); cout<<"combiFlag : "<<combiFlag<<endl;
                select_j1j2 = Mass_j1j3; select_j3j4 = Mass_j2j4;

                ratio_P1 = (z_fault_1 + z_fault_3).E()/(z_fault_1 + z_fault_3 + z_yes_1 + z_yes_3).E();
                ratio_P2 = (z_fault_2 + z_fault_4).E()/(z_fault_2 + z_fault_4 + z_yes_2 + z_yes_4).E();
                if(rescharge_1 == 1)
                {
                    angle_j1j2 = TV_GenJet_1.Angle(TV_GenJet_3); angle_j3j4 = TV_GenJet_2.Angle(TV_GenJet_4);
                    select_j1j2 = Mass_j1j3; select_j3j4 = Mass_j2j4;
                    ratio_1 = ((fault_1 + fault_3).E())/((yes_1 + yes_3 + fault_1 + fault_3).E());
                    ratio_2 = ((fault_2 + fault_4).E())/((fault_2 + fault_4 + yes_2 + yes_4).E());
                    newratio_1 = ((yes_1 + yes_3).E())/((yes_1 + yes_3 + fault_2 + fault_4).E());
                    newratio_2 = ((yes_2 + yes_4).E())/((yes_2 + yes_4 + fault_1 + fault_3).E());
                    ratio_WP_1 = fault_1.E()/(fault_1 + yes_1).E();
                    ratio_WP_2 = fault_3.E()/(fault_3 + yes_3).E();
                    ratio_WM_1 = fault_2.E()/(fault_2 + yes_2).E();
                    ratio_WM_2 = fault_4.E()/(fault_4 + yes_4).E();
                }
                else if(rescharge_1 == -1)
                {
                    angle_j1j2 = TV_GenJet_2.Angle(TV_GenJet_4); angle_j3j4 = TV_GenJet_1.Angle(TV_GenJet_3);
                    select_j1j2 = Mass_j2j4; select_j3j4 = Mass_j1j3;
                    ratio_1 = ((fault_2 + fault_4).E())/((fault_2 + fault_4 + yes_2 + yes_4).E());
                    ratio_2 = ((fault_1 + fault_3).E())/((yes_1 + yes_3 + fault_1 + fault_3).E());
                    newratio_1 = ((yes_2 + yes_4).E())/((yes_2 + yes_4 + fault_1 + fault_3).E());
                    newratio_2 = ((yes_1 + yes_3).E())/((yes_1 + yes_3 + fault_2 + fault_4).E());
                    ratio_WP_1 = fault_2.E()/(fault_2 + yes_2).E();
                    ratio_WP_2 = fault_4.E()/(fault_4 + yes_4).E();
                    ratio_WM_1 = fault_1.E()/(fault_1 + yes_1).E();
                    ratio_WM_2 = fault_3.E()/(fault_3 + yes_3).E();
                    
                }
            }
            else if(combiFlag == 2)
            {
//                select_j1j2 = Mass_j1j4; select_j3j4 = Mass_j2j3; ratio_1 = ((fault_4 + fault_1).E())/((fault_1 + fault_4 + yes_4 + yes_1).E()); ratio_2 = ((fault_3 + fault_2).E())/((fault_2 + fault_3 + yes_3 + yes_2).E()); cout<<"combiFlag : "<<combiFlag<<endl;
                select_j1j2 = Mass_j1j4; select_j3j4 = Mass_j2j3;

                ratio_P1 = (z_fault_1 + z_fault_4).E()/(z_fault_1 + z_fault_4 + z_yes_1 + z_yes_4).E();
                ratio_P2 = (z_fault_3 + z_fault_2).E()/(z_fault_3 + z_fault_2 + z_yes_3 + z_yes_2).E();
               if(rescharge_1 == 1)
               {
                   angle_j1j2 = TV_GenJet_1.Angle(TV_GenJet_4); angle_j3j4 = TV_GenJet_2.Angle(TV_GenJet_3);
                   select_j1j2 = Mass_j1j4; select_j3j4 = Mass_j2j3;
                   ratio_1 = ratio_1 = ((fault_4 + fault_1).E())/((fault_1 + fault_4 + yes_4 + yes_1).E());
                   ratio_2 = ((fault_3 + fault_2).E())/((fault_2 + fault_3 + yes_3 + yes_2).E());
                   newratio_1 = ((yes_1 + yes_4).E())/((yes_1 + yes_4 + fault_2 + fault_3).E());
                   newratio_2 = ((yes_2 + yes_3).E())/((yes_2 + yes_3 + fault_1 + fault_4).E());
                   ratio_WP_1 = fault_1.E()/(fault_1 + yes_1).E();
                   ratio_WP_2 = fault_4.E()/(fault_4 + yes_4).E();
                   ratio_WM_1 = fault_2.E()/(fault_2 + yes_2).E();
                   ratio_WM_2 = fault_3.E()/(fault_3 + yes_3).E();
               }
                else if(rescharge_1 == -1)
                {
                    angle_j1j2 = TV_GenJet_2.Angle(TV_GenJet_3); angle_j3j4 = TV_GenJet_1.Angle(TV_GenJet_4);
                    select_j1j2 = Mass_j2j3; select_j3j4 = Mass_j1j4;
                    ratio_1 = ((fault_3 + fault_2).E())/((fault_2 + fault_3 + yes_3 + yes_2).E());
                    ratio_2 = ((fault_4 + fault_1).E())/((fault_1 + fault_4 + yes_4 + yes_1).E());
                    newratio_1 = ((yes_2 + yes_3).E())/((yes_2 + yes_3 + fault_1 + fault_4).E());
                    newratio_2 = ((yes_1 + yes_4).E())/((yes_1 + yes_4 + fault_2 + fault_3).E());
                    ratio_WP_1 = fault_2.E()/(fault_2 + yes_2).E();
                    ratio_WP_2 = fault_3.E()/(fault_3 + yes_3).E();
                    ratio_WM_1 = fault_1.E()/(fault_1 + yes_1).E();
                    ratio_WM_2 = fault_4.E()/(fault_4 + yes_4).E();
                }
            }
            cout<<"ratio_Z1_1 : ratio_Z1_2 : ratio_Z2_1 : ratio_Z2_2 "<<ratio_Z1_1<<" : "<<ratio_Z1_2<<" : "<<ratio_Z2_1<<" : "<<ratio_Z2_2<<endl;
            cout<<"ratio_WP_1 : ratio_WP_2 : ratio_WM_1 : ratio_WM_2 "<<ratio_WP_1<<" : "<<ratio_WP_2<<" : "<<ratio_WM_1<<" : "<<ratio_WM_2<<endl;
            En_WP = TL_WP.E(); En_WM = TL_WM.E();
            cout<<"ratio_1 "<<ratio_1<<"  ratio_2  "<<ratio_2<<endl;
            cout<<"correct classification, newratio_1 "<<newratio_1<<" newratio_2  "<<newratio_2<<endl;
            cout<<"Mass_WP : angle_q1q2 : angle_j1j2 "<<En_WP<<" : "<<angle_q1q2<<" : "<<angle_j1j2<<"   Mass_WM : angle_q3q4 : angle_j3j4 "<<En_WM<<" : "<<angle_q3q4<<" : "<<angle_j3j4<<endl;
            if(angle_q1q2 == 0 || angle_q3q4 == 0){
                cout<<"**********  ********      ***********           **********"<<endl;
                cout<<"**             **         **       **           **      **"<<endl;
                cout<<"**             **         **       **           **      **"<<endl;
                cout<<"**********     **         **       **           **********"<<endl;
                cout<<"        **     **         **       **           **        "<<endl;
                cout<<"        **     **         **       **           **"<<endl;
                cout<<"**********     **         ***********           **"<<endl;
                
            }
            
        } catch (lcio::DataNotAvailableException err) {  }
	
        _outputTree->Fill();
        Num ++;
	}  	  

}	



void MCPJetClustering::end()
{

	if (_outputTree) {

		TFile *tree_file = _outputTree->GetCurrentFile(); //just in case we switched to a new file
		//tree_file->cd();
		tree_file->Write();
		delete tree_file;
		//tree_file->Close();
	}

}



