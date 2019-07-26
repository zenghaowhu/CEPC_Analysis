#include <higgs2zz.hh>
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

higgs2zz a_higgs2zz_instance;

higgs2zz::higgs2zz()
	: Processor("higgs2zz"),
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

void higgs2zz::init() {

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
    _outputTree->Branch("mass_Hllqq", &mass_Hllqq, "mass_Hllqq/F");
    _outputTree->Branch("Z1Dau_PDG", &Z1Dau_PDG, "Z1Dau_PDG/I");
    _outputTree->Branch("Z2Dau_PDG", &Z2Dau_PDG, "Z2Dau_PDG/I");
    _outputTree->Branch("angle_muon_mcp", &angle_muon_mcp, "angle_muon_mcp/F");
    _outputTree->Branch("angle_muon_pfo", &angle_muon_pfo, "angle_muon_pfo/F");
    _outputTree->Branch("mmMiss_mass", &mmMiss_mass, "mmMiss_mass/F");
    _outputTree->Branch("mmRecoil_mass", &mmRecoil_mass, "mmRecoil_mass/F");
    _outputTree->Branch("lepton", &lepton, "lepton/I");
    _outputTree->Branch("quark", &quark, "quark/I");
    _outputTree->Branch("neutrino", &neutrino, "neutrino/I");
    
	Num = 0;

}

void higgs2zz::processEvent( LCEvent * evtP )
{		

	if (evtP) 								
	{
        try{
            
            cout<<"************************************************************************"<<endl;
            
            std::vector<MCParticle* > MCP_PFO;
            std::vector<MCParticle* > Higgs_Daus;
            std::vector<MCParticle* > Higgs_Daus_ISR;
            
            MCParticle* HDau_1; MCParticle* HDau_2;
            
            TVector3 TVHDau_1(0,0,0), TVHDau_2(0,0,0);
            
            TLorentzVector TL_ori(0,0,0,240);
            TLorentzVector TLHDau_1(0,0,0,0), TLHDau_2(0,0,0,0);
            TLorentzVector ori_generator(0,0,0,0);
            TLorentzVector TLHiggs_Daus(0,0,0,0);
            
            Z1Dau_PDG = 999, Z2Dau_PDG = 999;
            
            LCCollection* col_MCP = evtP->getCollection( "MCParticle" );
            int nMCP = col_MCP->getNumberOfElements();
            eventNr = evtP->getEventNumber();
            cout<<"MCP "<<eventNr<<" Num "<<Num<<endl;
            
            //to record the information of MCParticle
            for(int i=0; i<nMCP; i++)
            {
                MCParticle* a_MCP = dynamic_cast<MCParticle*>(col_MCP->getElementAt(i));
                int NParents = a_MCP->getParents().size();
                int NDaughters = a_MCP->getDaughters().size();
                int PDG = a_MCP->getPDG();

                //Higgs boson
                if(abs(PDG) == 25 && NDaughters == 2)
                {
                    cout<<"the PDG of the daughters of Higgs boson : "<<a_MCP->getDaughters()[0]->getPDG()<<" : "<<a_MCP->getDaughters()[1]->getPDG()<<endl;
                    HDau_1 = a_MCP->getDaughters()[0];
                    HDau_2 = a_MCP->getDaughters()[1];
                    TVHDau_1 = HDau_1->getMomentum();
                    TVHDau_2 = HDau_2->getMomentum();
                    TLorentzVector TLHDau_1(TVHDau_1,HDau_1->getEnergy());
                    TLorentzVector TLHDau_2(TVHDau_2,HDau_2->getEnergy());
                    Z1Dau_PDG = abs(HDau_1->getDaughters()[0]->getPDG());
                    Z2Dau_PDG = abs(HDau_2->getDaughters()[0]->getPDG());
                }
            }
            cout<<"the PDG of Z1 : Z2 "<<Z1Dau_PDG<<" : "<<Z2Dau_PDG<<endl;

            cout<<"test1"<<endl;
            
//            TLorentzVector test(0,0,0,0);
            for(int i=0; i<nMCP; i++)
            {
                MCParticle* a_MCP = dynamic_cast<MCParticle*>(col_MCP->getElementAt(i));
                MCParticle* a_parent = a_MCP;
                TLorentzVector mcp_tmp(a_MCP->getMomentum()[0], a_MCP->getMomentum()[1], a_MCP->getMomentum()[2], a_MCP->getEnergy());
                
                int NParents = a_MCP->getParents().size();
                int NDaughters = a_MCP->getDaughters().size();
                int PDG = abs(a_MCP->getPDG());
                
                if(NParents == 0)
                {
                    ori_generator += mcp_tmp;
                }
                
//                cout<<"test2"<<endl;
                
                if(NParents != 0 && NDaughters == 0 && a_MCP->isCreatedInSimulation() == false)
                {
                    if(PDG != 12 && PDG != 14 && PDG != 16)
                    {
                        do{
                            a_parent = a_parent->getParents()[0];
                     //       cout<<"the PDG of a_parent "<<a_parent->getPDG()<<endl;
                        }while(a_parent->getPDG() != 23 && a_parent->getParents().size() != 0);
//                        cout<<"a_parent->getPDG() "<<a_parent->getPDG()<<endl;
//                        cout<<"a_MCP->getPDG() "<<a_MCP->getPDG()<<endl;
                        if(a_parent->getPDG() == 23){Higgs_Daus.push_back(a_MCP);}
                        /*
                        if(a_parent->getPDG() == 23)
                        {
                            TLorentzVector test_(0,0,0,0);
                            test_.SetPxPyPzE(a_MCP->getMomentum()[0], a_MCP->getMomentum()[1], a_MCP->getMomentum()[2], a_MCP->getEnergy());
                            test += test_;
                        }
                        */
                    }
                }
                
            /*
                if(NDaughters == 0 && !(a_MCP->isCreatedInSimulation()) && PDG != 12 && PDG != 14 && PDG != 16)
                {
                    Higgs_Daus_ISR.push_back(a_MCP);
                }
                */
            }
//            cout<<"test.M() "<<test.M()<<endl;

//            cout<<"test3"<<endl;
            int numHDaus = Higgs_Daus.size();
//            cout<<"   "<<numHDaus<<endl;
            for(int i=0; i<numHDaus; i++)
            {
                MCParticle* HiggsDau = Higgs_Daus.at(i);
                TVector3 TVHiggsDau(HiggsDau->getMomentum());
                TLorentzVector TLHiggsDau(HiggsDau->getMomentum()[0], HiggsDau->getMomentum()[1], HiggsDau->getMomentum()[2], HiggsDau->getEnergy());
                TLHiggs_Daus += TLHiggsDau;
//                cout<<"TLHiggsDau "<<TLHiggsDau.M()<<"    "<<TLHiggs_Daus.M()<<endl;
            }
            
            cout<<"for nnhzz, the invariant mass of Higgs in MCTruth is "<<TLHiggs_Daus.M()<<endl;
            
            MCParticle* MCP_muPlus; MCParticle* MCP_muMinus;
            TVector3 TV_MCP_muMinus(0,0,0), TV_MCP_muPlus(0,0,0);
            TLorentzVector TL_MCP_muMinus(0,0,0,0), TL_MCP_muPlus(0,0,0,0);
            angle_muon_mcp = 999;
            for(int i=0; i<nMCP; i++)
            {
                MCParticle* a_MCP = dynamic_cast<MCParticle*>(col_MCP->getElementAt(i));
                int NParents = a_MCP->getParents().size();
                int NDaughters = a_MCP->getDaughters().size();
                int PDG = a_MCP->getPDG();
                
                if(NParents == 0 && abs(PDG) == 13)
                {
                    if(PDG == 13){MCP_muMinus = a_MCP; TV_MCP_muMinus = a_MCP->getMomentum(); TL_MCP_muMinus.SetPxPyPzE(a_MCP->getMomentum()[0], a_MCP->getMomentum()[1], a_MCP->getMomentum()[2], a_MCP->getEnergy());}
                    else if(PDG == -13){MCP_muPlus = a_MCP; TV_MCP_muPlus = a_MCP->getMomentum(); TL_MCP_muPlus.SetPxPyPzE(a_MCP->getMomentum()[0], a_MCP->getMomentum()[1], a_MCP->getMomentum()[2], a_MCP->getEnergy());}
                }
            }
            angle_muon_mcp = TV_MCP_muPlus.Angle(TV_MCP_muMinus);
            cout<<"the angle of two muons in MCP is "<<angle_muon_mcp<<endl;
            cout<<"the invarint mass of Z calculated with two muons in MCTruth "<<(TL_MCP_muPlus + TL_MCP_muMinus).M()<<endl;
            
            
            lepton = 0, quark = 0, neutrino = 0;
            if(Z1Dau_PDG == 11 || Z1Dau_PDG == 13){lepton = 1;}
            else if(Z1Dau_PDG == 1 || Z1Dau_PDG == 2 || Z1Dau_PDG == 3 || Z1Dau_PDG ==4 || Z1Dau_PDG == 5 || Z1Dau_PDG == 6){quark = 1;}
            else if(Z1Dau_PDG == 12 || Z1Dau_PDG == 14 || Z1Dau_PDG == 16){neutrino = 1;}
            
            
            if(Z2Dau_PDG == 11 || Z2Dau_PDG == 13){lepton = 1;}
            else if(Z2Dau_PDG == 1 || Z2Dau_PDG ==2 || Z2Dau_PDG == 3 || Z2Dau_PDG == 4 || Z2Dau_PDG ==5 || Z2Dau_PDG ==6){quark = 1;}
            else if(Z2Dau_PDG ==12 || Z2Dau_PDG == 14 || Z2Dau_PDG == 16){neutrino == 1;}
            
            LCCollection* col_Reco = evtP->getCollection( "ArborPFOs" );
            int nReco = col_Reco->getNumberOfElements();
            std::vector<ReconstructedParticle* > muPlus;
            std::vector<ReconstructedParticle* > muMinus;
            
            TLorentzVector TLReco_pfo(0,0,0,0);
            for(int i=0; i<nReco; i++)
            {
                ReconstructedParticle* a_Reco = dynamic_cast<ReconstructedParticle*>(col_Reco->getElementAt(i));
                int type = a_Reco->getType();
                TLorentzVector TLtemp_Reco(a_Reco->getMomentum()[0], a_Reco->getMomentum()[1], a_Reco->getMomentum()[2], a_Reco->getEnergy());
                TLReco_pfo += TLtemp_Reco;
                
                if(type == 13){muMinus.push_back(a_Reco);}
                else if(type == -13){muPlus.push_back(a_Reco);}
                //Higgs boson

            }
            
            mass_Hllqq = 999;
            
            mass_Hllqq = TLReco_pfo.M();
            cout<<"for nnhzz, the invariant mass of all ArborPFO is "<<mass_Hllqq<<endl;
            
            float angle_muMinus = 999, angle_muPlus = 999;
            int select_muMinus = 999, select_muPlus = 999;
            TVector3 TVmuMinus_temp(0,0,0), TVmuPlus_temp(0,0,0);
            cout<<"the number of muMinus and muPlus in ArborPFO "<<muMinus.size()<<" : "<<muPlus.size()<<endl;
            int num_muMinus_PFO = muMinus.size(), num_muPlus_PFO = muPlus.size();
            if(num_muPlus_PFO > 0 && num_muMinus_PFO > 0)
            {
                for(int i=0; i<muMinus.size(); i++)
                {
                    ReconstructedParticle* muMinus_temp = muMinus.at(i);
                    TVmuMinus_temp = muMinus_temp->getMomentum();
                    if(TVmuMinus_temp.Angle(TV_MCP_muMinus) < angle_muMinus){cout<<"judge 1"<<endl; angle_muMinus = TVmuMinus_temp.Angle(TV_MCP_muMinus); select_muMinus = i;}
                }
                for(int i=0; i<muPlus.size(); i++)
                {
                    ReconstructedParticle* muPlus_temp = muPlus.at(i);
                    TVmuPlus_temp = muPlus_temp->getMomentum();
                    if(TVmuPlus_temp.Angle(TV_MCP_muPlus) < angle_muPlus){cout<<"judge 2"<<endl; angle_muPlus = TVmuPlus_temp.Angle(TV_MCP_muPlus); select_muPlus = i;}
                }
            }
            cout<<"the angle muMinus and muPlus "<<angle_muMinus<<" : "<<angle_muPlus<<endl;
            cout<<"select the muon particle is "<<select_muMinus<<" : "<<select_muPlus<<endl;
            
            mmRecoil_mass = 999;
            angle_muon_pfo = 999;
            mmMiss_mass = 999;
           
            if(num_muMinus_PFO > 0 && num_muPlus_PFO > 0)
            {
                ReconstructedParticle* muMinus_pfo = muMinus.at(select_muMinus);
                ReconstructedParticle* muPlus_pfo = muPlus.at(select_muPlus);
                TLorentzVector TL_muMinus_pfo(0,0,0,0), TL_muPlus_pfo(0,0,0,0);
                TVector3 TV_muMinus_pfo(0,0,0), TV_muPlus_pfo(0,0,0);
                TV_muMinus_pfo = muMinus_pfo->getMomentum();
                TV_muPlus_pfo = muPlus_pfo->getMomentum();
                TL_muMinus_pfo.SetPxPyPzE(muMinus_pfo->getMomentum()[0], muMinus_pfo->getMomentum()[1], muMinus_pfo->getMomentum()[2], muMinus_pfo->getEnergy());
                TL_muPlus_pfo.SetPxPyPzE(muPlus_pfo->getMomentum()[0], muPlus_pfo->getMomentum()[1], muPlus_pfo->getMomentum()[2], muPlus_pfo->getEnergy());
                cout<<"the recoil mass(teo muons) is "<<(TL_ori - TL_muPlus_pfo - TL_muMinus_pfo).M()<<endl;
                mmRecoil_mass = (TL_ori - TL_muPlus_pfo - TL_muMinus_pfo).M();
                angle_muon_pfo = TV_muPlus_pfo.Angle(TV_muMinus_pfo);
                mmMiss_mass = (TL_ori - TLReco_pfo).M();
                cout<<"for mmhnnqq, the missing mass is "<<mmMiss_mass<<endl;
            }
            

            cout<<"the angle of two muons in ArborPFO "<<angle_muon_pfo<<endl;

            
            
            
            
            
            
            
        } catch (lcio::DataNotAvailableException err) {  }
	
        _outputTree->Fill();
        Num ++;
	}  	  

}	



void higgs2zz::end()
{

	if (_outputTree) {

		TFile *tree_file = _outputTree->GetCurrentFile(); //just in case we switched to a new file
		//tree_file->cd();
		tree_file->Write();
		delete tree_file;
		//tree_file->Close();
	}

}



