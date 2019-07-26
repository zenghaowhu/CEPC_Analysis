#include <softLink.hh>
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

softLink a_softLink_instance;

softLink::softLink()
	: Processor("softLink"),
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



void softLink::init() {

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

    _outputTree->Branch("v941Px", &v941Px);
    _outputTree->Branch("v941Py", &v941Py);
    _outputTree->Branch("v941Pz", &v941Pz);
    _outputTree->Branch("v941E", &v941E);
    
    _outputTree->Branch("v942Px", &v942Px);
    _outputTree->Branch("v942Py", &v942Py);
    _outputTree->Branch("v942Pz", &v942Pz);
    _outputTree->Branch("v942E", &v942E);
    
    _outputTree->Branch("vquarknear941Px", &vquarknear941Px);
    _outputTree->Branch("vquarknear941Py", &vquarknear941Py);
    _outputTree->Branch("vquarknear941Pz", &vquarknear941Pz);
    _outputTree->Branch("vquarknear941E", &vquarknear941E);

    _outputTree->Branch("vquarknear942Px", &vquarknear942Px);
    _outputTree->Branch("vquarknear942Py", &vquarknear942Py);
    _outputTree->Branch("vquarknear942Pz", &vquarknear942Pz);
    _outputTree->Branch("vquarknear942E", &vquarknear942E);
    
    _outputTree->Branch("vquark941Px", &vquark941Px);
    _outputTree->Branch("vquark941Py", &vquark941Py);
    _outputTree->Branch("vquark941Pz", &vquark941Pz);
    _outputTree->Branch("vquark941E", &vquark941E);

    _outputTree->Branch("vquark942Px", &vquark942Px);
    _outputTree->Branch("vquark942Py", &vquark942Py);
    _outputTree->Branch("vquark942Pz", &vquark942Pz);
    _outputTree->Branch("vquark942E", &vquark942E);

    
    

    
	Num = 0;
}

void softLink::processEvent( LCEvent * evtP )
{		

	if (evtP) 								
	{
        try{
            cout<<"************************************************************************"<<endl;
            eventNr = evtP->getEventNumber();
            cout<<"eventNr "<<eventNr<<" Num "<<Num<<endl;
            std::vector<MCParticle* > MCTQuark;
            std::vector<MCParticle* > col941;
            std::vector<MCParticle* > col942;
            std::vector<MCParticle* > other;
            std::vector<MCParticle* > neutrino;
            std::vector<MCParticle* > quark941;
            std::vector<MCParticle* > quark942;
            std::vector<MCParticle* > quarknear941;
            std::vector<MCParticle* > quarknear942;
            LCCollection* col_MCP = evtP->getCollection( "MCParticle" );
            int n_MCP = col_MCP->getNumberOfElements();
            int count = 0; //used to select the first 94 particle
            int counttest = 0;
            int countFinal = 0;
            TLorentzVector TLFinal(0,0,0,0);
            double mass941 = 999, mass942 = 999;

            
            for(int i = 0; i<n_MCP; i++)
            {
                MCParticle* a_MCP = dynamic_cast<MCParticle*>(col_MCP->getElementAt(i));
                int NParents = a_MCP->getParents().size();
                int NDaughters = a_MCP->getDaughters().size();
                int PDG = a_MCP->getPDG();
                TLorentzVector temp(a_MCP->getMomentum()[0], a_MCP->getMomentum()[1], a_MCP->getMomentum()[2], a_MCP->getEnergy());
                
                if(NParents == 0 && (abs(PDG) == 1 || abs(PDG) == 2 || abs(PDG) == 3 || abs(PDG) == 4 || abs(PDG) == 5 || abs(PDG) == 6))
                {
                    cout<<"the energy of original parton is "<<temp.E()<<endl;
                }
                
                
                
                if(NDaughters == 0){
                    TLFinal += temp;
                    countFinal += 1;
                }
                
                if(NDaughters == 0 && NParents != 0 && (abs(PDG) == 12 || abs(PDG) == 14 || abs(PDG) == 16)){
                    neutrino.push_back(a_MCP);
                }
                
                if(NDaughters == 0 && NParents != 0 && abs(PDG) != 12 && abs(PDG) != 14 && abs(PDG) != 16){
                    MCParticle* a_parent = a_MCP;
                    
                    do{a_parent = a_parent->getParents()[0];}
                    while(a_parent->getPDG() != 92 && a_parent->getParents().size() != 0);
                    
                    if(a_parent->getPDG() == 92){
                        if(count == 0){
                            count += 1; mass941 = a_parent->getMass();
                            int num_parent92 = a_parent->getParents().size();
                            for(int j = 0; j<num_parent92; j++){    //only used to find out the parton which decay to 92
                                MCParticle* parent92 = a_parent->getParents()[j];
                                int PDG_parent92 = parent92->getPDG();
                                if(abs(PDG_parent92) == 1 || abs(PDG_parent92) == 2 || abs(PDG_parent92) == 3 || abs(PDG_parent92) == 4 || abs(PDG_parent92) == 5 || abs(PDG_parent92) == 6){
                                    quarknear941.push_back(parent92);
                                    do{parent92 = parent92->getParents()[0];}
                                    while(parent92->getParents().size() != 0);
                                    if(parent92->getParents().size() == 0){quark941.push_back(parent92);}
                                }
                            }
                        }
                        
                        
                        if(a_parent->getMass() == mass941){
                            col941.push_back(a_MCP);
                        }
                        else if(a_parent->getMass() != mass941){
                            if(counttest == 0){
                                counttest += 1;
                                int num_parent92 = a_parent->getParents().size();
                                for(int j = 0; j<num_parent92; j++){    //only used to find out the parton which decay to 92
                                    MCParticle* parent92 = a_parent->getParents()[j];
                                    int PDG_parent92 = parent92->getPDG();
                                    if(abs(PDG_parent92) == 1 || abs(PDG_parent92) == 2 || abs(PDG_parent92) == 3 || abs(PDG_parent92) == 4 || abs(PDG_parent92) == 5 || abs(PDG_parent92) == 6){
                                        quarknear942.push_back(parent92);
                                        do{parent92 = parent92->getParents()[0];}
                                        while(parent92->getParents().size() != 0);
                                        if(parent92->getParents().size() == 0){quark942.push_back(parent92);}
                                    }
                                }
                            }
                            mass942 = a_parent->getMass();
                            col942.push_back(a_MCP);
                        }
                    }
                    else if(a_parent->getParents().size() == 0){
                        other.push_back(a_MCP);
                    }
                }
            }
            
            
            vquark941Px.clear(); vquark941Py.clear(); vquark941Pz.clear(); vquark941E.clear();
            vquark942Px.clear(); vquark942Py.clear(); vquark942Pz.clear(); vquark942E.clear();
            
            int num_quark941 = quark941.size(), num_quark942 = quark942.size();
            
            cout<<"num_quark941 : "<<num_quark941<<" num_quark942 : "<<num_quark942<<endl;
            
            for(int i = 0; i<num_quark941; i++){
                MCParticle* a_MCP = quark941.at(i);
                vquark941Px.push_back(a_MCP->getMomentum()[0]);
                vquark941Py.push_back(a_MCP->getMomentum()[1]);
                vquark941Pz.push_back(a_MCP->getMomentum()[2]);
                vquark941E.push_back(a_MCP->getEnergy());
            }
            
            for(int i = 0; i<num_quark942; i++){
                MCParticle* a_MCP = quark941.at(i);
                vquark942Px.push_back(a_MCP->getMomentum()[0]);
                vquark942Py.push_back(a_MCP->getMomentum()[1]);
                vquark942Pz.push_back(a_MCP->getMomentum()[2]);
                vquark942E.push_back(a_MCP->getEnergy());
            }
            
            

            cout<<"mass941 : mass942 "<<mass941<<" : "<<mass942<<endl;
            cout<<"quarknear941.size() : "<<quarknear941.size()<<" quarknear942.size(): "<<quarknear942.size()<<endl;
            
            vquarknear941Px.clear(); vquarknear941Py.clear(); vquarknear941Pz.clear(); vquarknear941E.clear();
            vquarknear942Px.clear(); vquarknear942Py.clear(); vquarknear942Pz.clear(); vquarknear942E.clear();
            int num_quarknear941 = quarknear941.size(), num_quarknear942 = quarknear942.size();
            
            for(int i = 0; i<num_quarknear941; i++){
                MCParticle* a_MCP = quarknear941.at(i);
                vquarknear941Px.push_back(a_MCP->getMomentum()[0]);
                vquarknear941Py.push_back(a_MCP->getMomentum()[1]);
                vquarknear941Pz.push_back(a_MCP->getMomentum()[2]);
                vquarknear941E.push_back(a_MCP->getEnergy());
            }
            
            for(int i = 0; i<num_quarknear942; i++){
                MCParticle* a_MCP = quarknear942.at(i);
                vquarknear942Px.push_back(a_MCP->getMomentum()[0]);
                vquarknear942Py.push_back(a_MCP->getMomentum()[1]);
                vquarknear942Pz.push_back(a_MCP->getMomentum()[2]);
                vquarknear942E.push_back(a_MCP->getEnergy());
            }
            
            
            
            
            v941Px.clear(); v941Py.clear(); v941Pz.clear(); v941E.clear();
            v942Px.clear(); v942Py.clear(); v942Pz.clear(); v942E.clear();
            
            int num_941 = col941.size(), num_942 = col942.size(), num_other = other.size(), num_neutrino = neutrino.size();
            TLorentzVector TL941(0,0,0,0);
            for(int i = 0; i<num_941; i++){
                MCParticle* a_MCP = col941.at(i);
                TLorentzVector temp(a_MCP->getMomentum(), a_MCP->getEnergy());
                TL941 += temp;
                v941Px.push_back(a_MCP->getMomentum()[0]);
                v941Py.push_back(a_MCP->getMomentum()[1]);
                v941Pz.push_back(a_MCP->getMomentum()[2]);
                v941E.push_back(a_MCP->getEnergy());
                
                
            }
            
            TLorentzVector TL942(0,0,0,0);
            for(int i = 0; i<num_942; i++){
                MCParticle* a_MCP = col942.at(i);
                TLorentzVector temp(a_MCP->getMomentum(), a_MCP->getEnergy());
                TL942 += temp;
                v942Px.push_back(a_MCP->getMomentum()[0]);
                v942Py.push_back(a_MCP->getMomentum()[1]);
                v942Pz.push_back(a_MCP->getMomentum()[2]);
                v942E.push_back(a_MCP->getEnergy());
            }
            
            int num_941parton = quark941.size(), num_942parton = quark942.size();
            for(int i = 0; i<num_941parton; i++){
                MCParticle* a_MCP = quark941.at(i);
                cout<<"the mass of selected parton is  "<<a_MCP->getEnergy()<<endl;
            }
            
            cout<<"......."<<endl;
            for(int i = 0; i<num_942parton; i++){
                MCParticle* a_MCP = quark942.at(i);
                cout<<"the mass of selected parton is  "<<a_MCP->getEnergy()<<endl;
            }
    
            TLorentzVector TLOther(0,0,0,0);
            for(int i = 0; i<num_other; i++){
                MCParticle* a_MCP = other.at(i);
                TLorentzVector temp(a_MCP->getMomentum(), a_MCP->getEnergy());
                TLOther += temp;
            }
            
            TLorentzVector TLNeutrino(0,0,0,0);
            for(int i = 0; i<num_neutrino; i++){
                MCParticle* a_MCP = neutrino.at(i);
                TLorentzVector temp(a_MCP->getMomentum(), a_MCP->getEnergy());
                TLNeutrino += temp;
            }

            cout<<"col941.size() + col942.size() + other.size() + neutrino.size() : "<<col941.size() + col942.size() + other.size() + neutrino.size()<<endl;
            cout<<"col941.size() : col942.size() : other.size() : neutrino.size() "<<col941.size()<<" : "<<col942.size()<<" : "<<other.size()<<" : "<<neutrino.size()<<endl;
            cout<<"countFinal : "<<countFinal<<endl;
            cout<<"(TL941 + TL942 + TLOther + TLNeutrino).M(): "<<(TL941 + TL942 + TLOther + TLNeutrino).M()<<endl;
            cout<<"TLFinal.M() : "<<TLFinal.M()<<endl;
            
            
        } catch (lcio::DataNotAvailableException err) {  }
	
        _outputTree->Fill();
        Num ++;
	}  	  

}	



void softLink::end()
{

	if (_outputTree) {

		TFile *tree_file = _outputTree->GetCurrentFile(); //just in case we switched to a new file
		//tree_file->cd();
		tree_file->Write();
		delete tree_file;
		//tree_file->Close();
	}

}

