#include <TauAna.hh>
#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/LCFloatVec.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <IMPL/MCParticleImpl.h>
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

TauAna a_TauAna_instance;

TauAna::TauAna()
	: Processor("TauAna"),
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

void TauAna::init() {

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
    _outputTree->Branch("PDG", &PDG, "PDG/I");
    _outputTree->Branch("Massj1j2", &Massj1j2, "Massj1j2/F");
    _outputTree->Branch("angle12_dif", &angle12_dif, "angle12_dif/F");
    _outputTree->Branch("angle_mcp", &angle_mcp, "angle_mcp/F");
    _outputTree->Branch("jet_num", &jet_num, "jet_num/I");
    _outputTree->Branch("angle_jet", &angle_jet, "angle_jet/F");
	Num = 0;
}

void TauAna::processEvent( LCEvent * evtP ) 
{		

	if (evtP) 								
	{
        try{
            
            
            
            TVector3 mcp_TVector3(0,0,0);
            TVector3 mcp_tauP(0,0,0);
            TVector3 mcp_tauM(0,0,0);
            TVector3 lead_TVector3(0,0,0);
            TVector3 leadTVector3_tauPD(0,0,0);
            TVector3 leadTVector3_tauMD(0,0,0);
            TVector3 sumTVector3_tauP(0,0,0);
            TVector3 sumTVector3_tauM(0,0,0);
            TVector3 ISN(0,0,0);
            TVector3 FSN(0,0,0);
            
            TLorentzVector sum_tauPD(0,0,0,0);
            TLorentzVector sum_tauMD(0,0,0,0);
            TLorentzVector ISNTL(0,0,0,0);
            TLorentzVector FSNTL(0,0,0,0);
            TLorentzVector deNeutrino(0,0,0,0);
            
            std::vector<MCParticle* > mcp_tauPDau;
            std::vector<MCParticle* > mcp_tauMDau;
            
            MCParticle* lead_mcp;
            MCParticle* lead_tauPD;
            MCParticle* lead_tauMD;
            
            float leadEn_tauPD = 0;
            float leadEn_tauMD = 0;
            LCCollection* col_MCP = evtP->getCollection( "MCParticle" );
            int nMCP = col_MCP->getNumberOfElements();
            int PDG_HDau = 999;
            eventNr = evtP->getEventNumber();
            cout<<"MCP "<<eventNr<<" Num "<<Num<<endl;
            
            //to record the information of MCParticle
            
            for(int i=0; i<nMCP; i++)
            {
                MCParticle* a_MCP = dynamic_cast<MCParticle*>(col_MCP->getElementAt(i));
                int NParents = a_MCP->getParents().size();
                int NDaughters = a_MCP->getDaughters().size();
                PDG = a_MCP->getPDG();
                if(abs(PDG) == 25 && NDaughters == 2){PDG_HDau = abs(a_MCP->getDaughters()[0]->getPDG()); cout<<a_MCP->getDaughters()[0]->getPDG()<<" : "<<a_MCP->getDaughters()[1]->getPDG()<<endl;}
            }
            
            for(int i=0; i<nMCP; i++)
            {
                MCParticle* a_MCP = dynamic_cast<MCParticle*>(col_MCP->getElementAt(i));
                MCParticle* a_parent = a_MCP;
                TLorentzVector mcp_tmp(a_MCP->getMomentum()[0], a_MCP->getMomentum()[1], a_MCP->getMomentum()[2], a_MCP->getEnergy());
                
                int NParents = a_MCP->getParents().size();
                int NDaughters = a_MCP->getDaughters().size();
                PDG = a_MCP->getPDG();
                float mcp_energy = a_MCP->getEnergy();
                mcp_TVector3 = a_MCP->getMomentum();
                
                if(NParents!= 0 && a_MCP->getParents()[0]->getPDG() == 25)
                {
                    if(PDG == PDG_HDau){mcp_tauM = mcp_TVector3;}
                    if(PDG == -PDG_HDau){mcp_tauP = mcp_TVector3;}
                }
                
                if(NDaughters == 0)
                {
                    if(abs(PDG) == 12 || abs(PDG) == 14 || abs(PDG) == 16)
                    {
                        FSNTL += mcp_tmp;
                        FSN += mcp_TVector3;
                        if(NParents == 0){ISN += mcp_TVector3; ISNTL += mcp_tmp;}
                    }

                    if(NParents == 0)continue;
//                    if(NParents == 0 || abs(PDG) == 12 || abs(PDG) == 14 || abs(PDG) == 16)continue;
                    do{
                        a_parent = a_parent->getParents()[0];
                    }while(a_parent->getParents().size() > 0 && abs(a_parent->getPDG()) != abs(PDG_HDau));
                    if((abs(PDG) == 12 || abs(PDG) == 14 || abs(PDG) == 16) && abs(a_parent->getPDG()) == PDG_HDau){deNeutrino += mcp_tmp;}
                    if(a_parent->getPDG() == PDG_HDau){
                        if(abs(PDG) != 12 && abs(PDG) != 14 && abs(PDG) != 16){mcp_tauMDau.push_back(a_MCP);}
                    }
                    else if(a_parent->getPDG() == -PDG_HDau && abs(PDG) != 12 && abs(PDG) != 14 && abs(PDG) != 16){mcp_tauPDau.push_back(a_MCP);}
                }

            }
            
            //the daughters of tau_plus
            int n_tauPD = mcp_tauPDau.size();
        
            for(int i=0; i<n_tauPD; i++)
            {
                MCParticle* tauPD = mcp_tauPDau[i];
                TLorentzVector mcp_tmp(tauPD->getMomentum()[0], tauPD->getMomentum()[1], tauPD->getMomentum()[2], tauPD->getEnergy());
                sum_tauPD += mcp_tmp;
                sumTVector3_tauP += tauPD->getMomentum();
                if(tauPD->getEnergy() > leadEn_tauPD)
                {
                    lead_tauPD = tauPD;
                    leadEn_tauPD = tauPD->getEnergy();
                    leadTVector3_tauPD = tauPD->getMomentum();
                }
            }

            
            
            //the daughters of tau_minus
            int n_tauMD = mcp_tauMDau.size();
            cout<<"the number of daughters of tau_minus "<<n_tauMD<<endl;
            for(int i=0; i<n_tauMD; i++)
            {
                MCParticle* tauMD = mcp_tauMDau[i];
                TLorentzVector mcp_tmp(tauMD->getMomentum()[0], tauMD->getMomentum()[1], tauMD->getMomentum()[2], tauMD->getEnergy());
                sum_tauMD += mcp_tmp;
 //               cout<<"one of tauMD energy "<<sum_tauMD.E()<<endl;

                sumTVector3_tauM += tauMD->getMomentum();
                if(tauMD->getEnergy() > leadEn_tauMD)
                {
                    lead_tauMD = tauMD;
                    leadEn_tauMD = tauMD->getEnergy();
                    leadTVector3_tauMD = tauMD->getMomentum();
                }
            }

            
            std::vector<ReconstructedParticle*> jet;
            ReconstructedParticle* jet_1;
            ReconstructedParticle* jet_2;
            TVector3 Tjet_1(0,0,0);
            TVector3 Tjet_2(0,0,0);
            angle_jet;
            jet_num = 0;
            
            LCCollection* col_Fastjet = evtP->getCollection( "FastJets" );
            int nJet = col_Fastjet->getNumberOfElements();
            for(int i=0; i<nJet; i++)
            {
                ReconstructedParticle * a_Jet = dynamic_cast<ReconstructedParticle*>(col_Fastjet->getElementAt(i));
                jet.push_back(a_Jet);
            }
            
            jet_num = jet.size();
            if(jet.size() == 2)
            {
                jet_1 = jet.at(0);
                jet_2 = jet.at(1);
                Tjet_1 = jet_1->getMomentum();
                Tjet_2 = jet_2->getMomentum();
                angle_jet = Tjet_1.Angle(Tjet_2);
                
            }
            
            //record the information of Reconstructed Particle
            std::vector<ReconstructedParticle* > reco;
            int type;
            int jet1_ncom = 0, jet2_ncom = 0, jet3_ncom=0, jet4_ncom=0, select=0, count=0, count1=0;
            float invM1 = 0, invM2 = 0;
            float pfo_energy = 0;
            TVector3 pfo_TVector3(0,0,0);
            TVector3 jet_first(0,0,0);
            TVector3 jet_second(0,0,0);
            TVector3 jet_third(0,0,0);
            TVector3 jet_fourth(0,0,0);
            TLorentzVector jet_1(0,0,0,0);
            TLorentzVector jet_2(0,0,0,0);
            TLorentzVector jet_3(0,0,0,0);
            TLorentzVector jet_4(0,0,0,0);
            TLorentzVector pfo_tmp(0,0,0,0);
            TLorentzVector ArborPFO(0,0,0,0);
            TLorentzVector allArborPOF(0,0,0,0);
            
            LCCollection* col_RecoP = evtP->getCollection( "ArborPFOs" );
            int nRecoP = col_RecoP->getNumberOfElements();
            for(int i=0; i<nRecoP; i++)
            {
                ReconstructedParticle * a_RecoP = dynamic_cast<ReconstructedParticle*>(col_RecoP->getElementAt(i));
                ArborPFO.SetPxPyPzE(a_RecoP->getMomentum()[0], a_RecoP->getMomentum()[1], a_RecoP->getMomentum()[2], a_RecoP->getEnergy());
                allArborPOF += ArborPFO;
                reco.push_back(a_RecoP);
                type = a_RecoP->getType();
                pfo_energy = a_RecoP->getEnergy();
            }
            
            if(reco.size() >= 2)
            {
                std::sort(reco.begin(), reco.end(), CRP());
                ReconstructedParticle *ori_reco = reco.at(0);
                TVector3 ori_recoTV = ori_reco->getMomentum();
                for(int i=1; i<reco.size(); i++)
                {
                    ReconstructedParticle *temp_reco = reco.at(i);
                    TVector3 temp_recoTV = temp_reco->getMomentum();
//                    cout<<temp_recoTV.Angle(ori_recoTV)<<endl;
                    if(temp_recoTV.Angle(ori_recoTV) > 1.9)
                    {
//                        cout<<"test*******************"<<endl;
                        count += 1;
                        if(count == 1){select = i;cout<<"select"<<select<<endl;}
                    }
                }
                
//                std::sort(reco.begin(), reco.end(), CRP());
                ReconstructedParticle *lead_first = reco.at(0);
                ReconstructedParticle *lead_second = reco.at(select);
//                ReconstructedParticle *lead_third = reco.at(2);
//                ReconstructedParticle *lead_fourth = reco.at(3);
                jet_first = lead_first->getMomentum();
                jet_second = lead_second->getMomentum();
 //               jet_third = lead_third->getMomentum();
 //               jet_fourth = lead_fourth->getMomentum();
                jet_1.SetPxPyPzE(lead_first->getMomentum()[0], lead_first->getMomentum()[1], lead_first->getMomentum()[2], lead_first->getEnergy());
                jet_2.SetPxPyPzE(lead_second->getMomentum()[0], lead_second->getMomentum()[1], lead_second->getMomentum()[2], lead_second->getEnergy());
 //               jet_3.SetPxPyPzE(lead_third->getMomentum()[0], lead_third->getMomentum()[1], lead_third->getMomentum()[2], lead_third->getEnergy());
 //               jet_4.SetPxPyPzE(lead_fourth->getMomentum()[0], lead_fourth->getMomentum()[1], lead_fourth->getMomentum()[2], lead_fourth->getEnergy());
                
                for(int i=1; i<reco.size(); i++)
                {
                    if(i == select)continue;
                    ReconstructedParticle *pfo_temp = reco.at(i);
                    pfo_tmp.SetPxPyPzE(pfo_temp->getMomentum()[0], pfo_temp->getMomentum()[1], pfo_temp->getMomentum()[2], pfo_temp->getEnergy());
                    pfo_TVector3 = pfo_temp->getMomentum();
                    float angle1 = pfo_TVector3.Angle(jet_first);
                    float angle2 = pfo_TVector3.Angle(jet_second);
 //                   float angle3 = pfo_TVector3.Angle(jet_third);
 //                   float angle4 = pfo_TVector3.Angle(jet_fourth);
                    if(angle1 < angle2){jet1_ncom += 1; jet_first += pfo_TVector3; jet_1 += pfo_tmp;}
                    else if(angle2 <= angle1){jet2_ncom += 1; jet_second += pfo_TVector3; jet_2 += pfo_tmp;}
//                    if(angle1 <= angle2 && angle1 <= angle3 && angle1 <= angle4 ){jet1_ncom +=1; jet_first += pfo_TVector3; jet_1 += pfo_tmp;}
//                    else if(angle2 <= angle1 && angle2 <= angle3 && angle2 <= angle4){jet2_ncom += 1; jet_second += pfo_TVector3; jet_2 += pfo_tmp;}
 //                   else if(angle3 <= angle1 && angle3 <= angle2 && angle3 <= angle4){jet3_ncom += 1; jet_third += pfo_TVector3; jet_3 += pfo_tmp;}
 //                   else if(angle4 <= angle1 && angle4 <= angle2 && angle4 <= angle3){jet4_ncom += 1; jet_fourth += pfo_TVector3; jet_4 += pfo_tmp;}
                }
            }

    
            Massj1j2 = 0; angle12_dif = 0; angle_mcp = 0;
            std::vector<TLorentzVector> FourJet;
            FourJet.push_back(jet_1);
            FourJet.push_back(jet_2);
   /*         FourJet.push_back(jet_3);
            FourJet.push_back(jet_4);
            int par11=999, par12=999, par21=999, par22=999;
            float cut = 999;
            for(int i=1; i<4; i++)
            {
                invM1 = (FourJet[0] + FourJet[i]).M();
                if(abs(invM1 - 80) < cut){cut=abs(invM1 - 80); par11=0; par12=i;}
            }
            if(par12 == 1){invM2 = (FourJet[2] + FourJet[3]).M();}
            else if(par12 == 2){invM2 = (FourJet[1] + FourJet[3]).M();}
            else if(par12 == 3){invM2 = (FourJet[1] + FourJet[2]).M();}
    */
            
            
            
            angle_mcp = mcp_tauM.Angle(mcp_tauP);
            Massj1j2 = (jet_1 + jet_2).M();
            angle12_dif = mcp_tauM.Angle(mcp_tauP) - jet_first.Angle(jet_second);
            cout<<(FourJet[0] + FourJet[1]).M()<<endl;
            cout<<"the energy of all ArborPFO is "<<allArborPOF.E()<<endl;
            cout<<"the energy of jet_1 : jet_2 "<<jet_1.E()<<" : "<<jet_2.E()<<" : "<<(jet_1 + jet_2).E()<<endl;
            cout<<"the energy of two tau in MCTruth "<<sum_tauPD.E()<<" : "<<sum_tauMD.E()<<" : "<<(sum_tauMD + sum_tauPD).E()<<endl;
            cout<<"the mass of two tau in MCTruth "<<sum_tauPD.M()<<" : "<<sum_tauMD.M()<<endl;
            cout<<"reco.size() : jet1_ncom : jet2_ncom "<<reco.size()<<" : "<<jet1_ncom<<" : "<<jet2_ncom<<endl;
            cout<<"the angle of two jet pfo : MCTruth "<<jet_first.Angle(jet_second)<<" : "<<mcp_tauM.Angle(mcp_tauP)<<endl;
            cout<<"the angle between tau and reco jet "<<jet_first.Angle(mcp_tauP)<<" : "<<jet_first.Angle(mcp_tauM)<<" : "<<jet_second.Angle(mcp_tauP)<<" : "<<jet_second.Angle(mcp_tauM)<<endl;
            cout<<"the transverse momentum of neutrinos decayed from tau "<<(FSN - ISN).Perp()<<endl;
            cout<<"the sum of the energy of two jet and neutrinos decayed from quarks "<<(jet_1 + jet_2 + deNeutrino).E()<<endl;
            cout<<"the angle of two jet in MCTruth "<<sumTVector3_tauM.Angle(sumTVector3_tauP)<<endl;
            cout<<"the energy of neutrinos decayed from taus "<<(FSNTL - ISNTL).E()<<" : "<<deNeutrino.E()<<endl;
            cout<<".........................................................."<<endl;
        } catch (lcio::DataNotAvailableException err) {  }
	
        _outputTree->Fill();
        Num ++;
	}  	  

}	



void TauAna::end()
{

	if (_outputTree) {

		TFile *tree_file = _outputTree->GetCurrentFile(); //just in case we switched to a new file
		//tree_file->cd();
		tree_file->Write();
		delete tree_file;
		//tree_file->Close();
	}

}



