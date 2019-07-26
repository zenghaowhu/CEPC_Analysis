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

//std::sort(MCP_PFO_excludeISR.begin(), MCP_PFO_excludeISR.end(), CMP());

struct CMP {
    bool operator()(MCParticle *p1, MCParticle *p2) const
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
    _outputTree->Branch("angle_HDau", &angle_HDau, "angle_HDau/F");
    _outputTree->Branch("angle_ZDau", &angle_ZDau, "angle_ZDau/F");
    _outputTree->Branch("angle_ZpH1", &angle_ZpH1, "angle_ZpH1/F");
    _outputTree->Branch("angle_ZpH2", &angle_ZpH2, "angle_ZpH2/F");
    _outputTree->Branch("angle_ZmH1", &angle_ZmH1, "angle_ZmH1/F");
    _outputTree->Branch("angle_ZmH2", &angle_ZmH2, "angle_ZmH2/F");
    
	Num = 0;
    angle_ZDau = 999; angle_HDau = 999, angle_ZpH1 = 999, angle_ZpH2 = 999, angle_ZmH1 = 999, angle_ZmH2 = 999;
}

void TauAna::processEvent( LCEvent * evtP ) 
{		

	if (evtP) 								
	{
        try{
            
            cout<<"************************************************************************"<<endl;

            std::vector<MCParticle* > MCP_PFO;
            std::vector<MCParticle* > MCP_PFO_excludeISR;
            std::vector<MCParticle* > new_MCP_PFO_excludeISR;
            std::vector<MCParticle* > HiggsDaus;
            std::vector<MCParticle* > ZDaus;
            
            MCParticle* HDau_1;
            MCParticle* HDau_2;
            MCParticle* ZDau_plus;
            MCParticle* ZDau_minus;
            
            TVector3 TVZDau_plus(0,0,0), TVZDau_minus(0,0,0);
            TVector3 TVHDau_1(0,0,0), TVHDau_2(0,0,0);
            
            TLorentzVector ori_generator(0,0,0,0);
            
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
                    int PDG_HDau = abs(a_MCP->getDaughters()[0]->getPDG());
                    cout<<"the PDG of the daughters of Higgs boson : "<<a_MCP->getDaughters()[0]->getPDG()<<" : "<<a_MCP->getDaughters()[1]->getPDG()<<endl;
                    if(a_MCP->getDaughters()[0]->getPDG() > 0){HDau_1 = a_MCP->getDaughters()[0]; HDau_2 = a_MCP->getDaughters()[1];}
                    else if(a_MCP->getDaughters()[0]->getPDG() < 0){HDau_1 = a_MCP->getDaughters()[1]; HDau_2 = a_MCP->getDaughters()[0];}
                    TVHDau_1 = HDau_1->getMomentum();
                    TVHDau_2 = HDau_2->getMomentum();
                    TLorentzVector TLHDau_1(TVHDau_1,HDau_1->getEnergy());
                    TLorentzVector TLHDau_2(TVHDau_2,HDau_2->getEnergy());
                }
                
                //Z boson
                if(NParents == 0){
                    if(abs(PDG) == 1 || abs(PDG) == 2 || abs(PDG) == 3 || abs(PDG) == 4 || abs(PDG) == 5 || abs(PDG) == 6){
                        if(PDG > 0){ZDau_plus = a_MCP; TVZDau_plus = ZDau_plus->getMomentum();}
                        else if(PDG < 0){ZDau_minus = a_MCP; TVZDau_minus = ZDau_minus->getMomentum();}
                    }
                }
                
            }
            angle_HDau = TVHDau_1.Angle(TVHDau_2);
            angle_ZDau = TVZDau_minus.Angle(TVZDau_plus);
            angle_ZpH1 = TVZDau_plus.Angle(TVHDau_1);
            angle_ZpH2 = TVZDau_plus.Angle(TVHDau_2);
            angle_ZmH1 = TVZDau_minus.Angle(TVHDau_1);
            angle_ZmH2 = TVZDau_minus.Angle(TVHDau_2);
            
            cout<<"the angle of two Higgs daughters is "<<angle_HDau<<endl;
            cout<<"the angle of two Z boson daughters is "<<angle_ZDau<<endl;
            cout<<"the angle between ZDau_plus and HDau_1 "<<angle_ZpH1<<"  the angle between ZDau_plus and HDau_2  "<<angle_ZpH2<<endl;
            cout<<"the angle between ZDau_minus and HDau_1 "<<angle_ZmH1<<"   the angle between ZDau_minus and HDau_2   "<<angle_ZmH2<<endl;
            
            cout<<"the PDG of the daughters of Z boson : "<<ZDau_minus->getPDG()<<" : "<<ZDau_plus->getPDG()<<endl;
            
            int count_mcppfo = 0, count_ISR = 0;
            for(int i=0; i<nMCP; i++)
            {
                MCParticle* a_MCP = dynamic_cast<MCParticle*>(col_MCP->getElementAt(i));
                MCParticle* a_parent = a_MCP;
                TLorentzVector mcp_tmp(a_MCP->getMomentum()[0], a_MCP->getMomentum()[1], a_MCP->getMomentum()[2], a_MCP->getEnergy());
                
                int NParents = a_MCP->getParents().size();
                int NDaughters = a_MCP->getDaughters().size();
                int PDG = a_MCP->getPDG();
//                mcp_TVector3 = a_MCP->getMomentum();
                
                if(NParents == 0)
                {
                    ori_generator += mcp_tmp;
                }
                
                if(NDaughters == 0 && !(a_MCP->isCreatedInSimulation()) && NParents != 0 && abs(PDG) != 12 && abs(PDG) != 14 && abs(PDG) != 16)
                {
                        MCP_PFO.push_back(a_MCP);
                        count_mcppfo += 1;
                    do{
                        a_parent = a_parent->getParents()[0];
                    }while(a_parent->getPDG() != 25 && a_parent->getParents().size() != 0);
                    if(a_parent->getPDG() == 25){HiggsDaus.push_back(a_MCP); MCP_PFO_excludeISR.push_back(a_MCP);}
                    else if(a_parent->getParents().size() == 0){
                        if(a_parent->getPDG() == 22){count_ISR += 1;}
                        if(abs(a_parent->getPDG()) == 1 || abs(a_parent->getPDG()) == 2 || abs(a_parent->getPDG()) == 3 || abs(a_parent->getPDG()) == 4 || abs(a_parent->getPDG()) == 5 || abs(a_parent->getPDG()) == 6 ){
                            ZDaus.push_back(a_MCP); MCP_PFO_excludeISR.push_back(a_MCP);
                        }
                    }
                    
                }
            }
            cout<<"the total energy of generator(include neutrino and ISR) is "<<ori_generator.E()<<endl;
            cout<<"the number of MCP_PFO is "<<count_mcppfo<<endl;
            cout<<"the number of MCP_PFO_excludeISR is "<<MCP_PFO_excludeISR.size()<<endl;
            cout<<"the number of ISR stored in MCP_PFO is "<<count_ISR<<endl;
            cout<<"the number of generator particles decayed from HIggs "<<HiggsDaus.size()<<" the number of generator particles decayed from Z "<<ZDaus.size()<<endl;
            
            
            MCParticle* energetic_1; MCParticle* energetic_2; MCParticle* energetic_3; MCParticle* energetic_4;
            MCParticle* energetic_5; MCParticle* energetic_6; MCParticle* energetic_7; MCParticle* energetic_8;
            
            TVector3 TVenergetic_1(0,0,0), TVenergetic_2(0,0,0), TVenergetic_3(0,0,0), TVenergetic_4(0,0,0);
            TVector3 TVenergetic_5(0,0,0), TVenergetic_6(0,0,0), TVenergetic_7(0,0,0), TVenergetic_8(0,0,0);
            
            std::sort(MCP_PFO_excludeISR.begin(), MCP_PFO_excludeISR.end(), CMP());
            //see the eight energetic particles
            if(MCP_PFO_excludeISR.size() >= 8)
            {
                energetic_1 = MCP_PFO_excludeISR.at(0); energetic_2 = MCP_PFO_excludeISR.at(1); energetic_3 = MCP_PFO_excludeISR.at(2);
                energetic_4 = MCP_PFO_excludeISR.at(3); energetic_5 = MCP_PFO_excludeISR.at(4); energetic_6 = MCP_PFO_excludeISR.at(5);
                energetic_7 = MCP_PFO_excludeISR.at(6); energetic_8 = MCP_PFO_excludeISR.at(7);
                TVenergetic_1 = energetic_1->getMomentum(); TVenergetic_2 = energetic_2->getMomentum(); TVenergetic_3 = energetic_3->getMomentum();
                TVenergetic_4 = energetic_4->getMomentum(); TVenergetic_5 = energetic_5->getMomentum(); TVenergetic_6 = energetic_6->getMomentum();
                TVenergetic_7 = energetic_7->getMomentum(); TVenergetic_8 = energetic_8->getMomentum();
                cout<<"the energy of energetic particles in MCP_PFO_excludeISR is "<<energetic_1->getEnergy()<<" "<<energetic_2->getEnergy()<<" "<<energetic_3->getEnergy()<<" "<<energetic_4->getEnergy()<<" "<<energetic_5->getEnergy()<<" "<<energetic_6->getEnergy()<<" "<<energetic_7->getEnergy()<<" "<<energetic_8->getEnergy()<<endl;
            }
            
            //see the eight energetic particles' parent
            for(int i=0; i<8; i++){
                MCParticle* temp_MCPPFO = MCP_PFO_excludeISR.at(i);
                MCParticle* a_parent_ = temp_MCPPFO;
                //int NParents = temp_MCPPFO->getParents().size();
                //int the collection of MCP_PFO_excludeISR, all particles have have parents
                do{
                    a_parent_=a_parent_->getParents()[0];
                }while(a_parent_->getPDG() != 25 && a_parent_->getParents().size() != 0);
                if(a_parent_->getPDG() == 25){cout<<"decayed form Higgs "<<i<<endl;}
                else if(a_parent_->getParents().size() == 0){cout<<"decayed form Z boson "<<i<<endl;}
            }
            
            TVector3 TVfirst_jet(0,0,0), TVsecond_jet(0,0,0), TVthird_jet(0,0,0), TVfouth_jet(0,0,0);
            TLorentzVector TLfirst_jet(0,0,0,0), TLsecond_jet(0,0,0,0), TLthird_jet(0,0,0,0), TLfourth_jet(0,0,0,0);
            MCParticle* initial_jet = MCP_PFO_excludeISR.at(0);  TVfirst_jet = initial_jet->getMomentum();
            float angle_temp = 999; int select_MCPPFO = 999;
            for(int i=1; i<MCP_PFO_excludeISR.size(); i++)
            {
                    MCParticle* temp = MCP_PFO_excludeISR.at(i);
                    new_MCP_PFO_excludeISR.push_back(temp);
            }
            MCP_PFO_excludeISR.clear();
            for(int i=0; i<new_MCP_PFO_excludeISR.size(); i++)
            {
                MCParticle* temp = new_MCP_PFO_excludeISR.at(i);
                MCP_PFO_excludeISR.push_back(temp);
            }
            new_MCP_PFO_excludeISR.clear();
            
            do{
                for(int i=0; i<MCP_PFO_excludeISR.size(); i++)
                {
                    MCParticle* temp_MCPPFO = MCP_PFO_excludeISR.at(i);
                    TVector3 temp_TV(temp_MCPPFO->getMomentum());
                    if(TVfirst_jet.Angle(temp_TV) < angle_temp)
                    {
                        angle_temp = TVfirst_jet.Angle(temp_TV);
                        select_MCPPFO = i;
                    }
                }
                
                MCParticle* select_particle = MCP_PFO_excludeISR.at(select_MCPPFO);
                TVector3 TVselect(0,0,0); TVselect = select_particle->getMomentum();
                TLorentzVector TLselect(TVselect,select_particle->getEnergy()); // TLselect.SetPxPyPzE(TVselect, select_particle->getEnergy());
                TVector3 TVfirst_jet_new(0,0,0); TVfirst_jet_new = TVfirst_jet + TVselect;
                TLorentzVector TLfirst_jet_new(0,0,0,0); TLfirst_jet_new = TLfirst_jet + TLselect;
                float angle_iter = TVfirst_jet_new.Angle(TVfirst_jet);
                float energy_iter_difference = TLfirst_jet_new.E() - TLfirst_jet.E();
                cout<<"the iterate angle is "<<angle_iter<<endl;
                cout<<"the iterate energy difference "<<energy_iter_difference<<endl;
                TVfirst_jet += TVselect;
                TLfirst_jet += TLselect;
                
                for(int i=0; i<MCP_PFO_excludeISR.size(); i++)
                {
                    if(i != select_MCPPFO)
                    {
                        MCParticle* temp = MCP_PFO_excludeISR.at(i);
                        new_MCP_PFO_excludeISR.push_back(temp);
                    }

                }
                MCP_PFO_excludeISR.clear();
                cout<<"MCP_PFO_excludeISR.size() "<<MCP_PFO_excludeISR.size()<<endl;
                for(int i=0; i<new_MCP_PFO_excludeISR.size(); i++)
                {
                    MCParticle* temp = new_MCP_PFO_excludeISR.at(i);
                    MCP_PFO_excludeISR.push_back(temp);
                }
                new_MCP_PFO_excludeISR.clear();
                cout<<"after excluding the closet particle, MCP_PFO_excludeISR.size() "<<MCP_PFO_excludeISR.size()<<endl;

            }while(MCP_PFO_excludeISR.size() > 0);

            
                
            /*
                float angle_temp = TVfirst_jet
                TLorentzVector temp_TL(temp_MCPPFO->getMomentum()[0], temp_MCPPFO->getMomentum()[1], temp_MCPPFO->getMomentum()[2], temp_MCPPFO->getEnergy());
                TVector3 TVfirst_jet_new(0,0,0); TVfirst_jet_new = TVfirst_jet + temp_TV;
                TLorentzVector TLfirst_jet_new(0,0,0,0); TLfirst_jet_new = TLfirst_jet + temp_TL;
                float angle_iter = TVfirst_jet_new.Angle(TVfirst_jet);
                float energy_iter_difference = TLfirst_jet_new.E() - TLfirst_jet.E();
                cout<<"the iterate angle is "<<angle_iter<<endl;
                cout<<"the iterate energy difference "<<energy_iter_difference<<endl;
                TVfirst_jet += temp_TV;
                TLfirst_jet += temp_TL;
            */
            
            
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



