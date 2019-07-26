#include <softLink2quark.hh>
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

softLink2quark a_softLink2quark_instance;

softLink2quark::softLink2quark()
	: Processor("softLink2quark"),
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



void softLink2quark::init() {

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

    _outputTree->Branch("vfinalPx", &vfinalPx);
    _outputTree->Branch("vfinalPy", &vfinalPy);
    _outputTree->Branch("vfinalPz", &vfinalPz);
    _outputTree->Branch("vfinalE",  &vfinalE);
    
    _outputTree->Branch("vquarkPx", &vquarkPx);
    _outputTree->Branch("vquarkPy", &vquarkPy);
    _outputTree->Branch("vquarkPz", &vquarkPz);
    _outputTree->Branch("vquarkE", &vquarkE);
    _outputTree->Branch("PDG_quark1", &PDG_quark1, "PDG_quark1/I");
    _outputTree->Branch("PDG_quark2", &PDG_quark2, "PDG_quark2/I");
    
    
	Num = 0;
}

void softLink2quark::processEvent( LCEvent * evtP )
{		

	if (evtP) 								
	{
        try{
            cout<<"************************************************************************"<<endl;
            eventNr = evtP->getEventNumber();
            cout<<"eventNr "<<eventNr<<" Num "<<Num<<endl;
            std::vector<MCParticle* > MCTQuark;
            std::vector<MCParticle* > finalFrom92;
            
            LCCollection* col_MCP = evtP->getCollection( "MCParticle" );
            int n_MCP = col_MCP->getNumberOfElements();
            
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
                    MCTQuark.push_back(a_MCP);
                }
                
                if(NParents == 0 && PDG == 22)
                {
                    cout<<"the energy of original photon is "<<temp.E()<<endl;
                }
                
                if(NParents == 0 && (abs(PDG) == 11 || abs(PDG) == 12 || abs(PDG) == 13 || abs(PDG) == 14 || abs(PDG) == 15 || abs(PDG) == 16))
                {
                    cout<<"the energy of original lepton is "<<temp.E()<<endl;
                }
                
                
                if(NDaughters == 0 && NParents != 0 && abs(PDG) != 12 && abs(PDG) != 14 && abs(PDG) != 16){
                    MCParticle* a_parent = a_MCP;
                    do{a_parent = a_parent->getParents()[0];}
                    while(a_parent->getPDG() != 92 && a_parent->getParents().size() != 0);
                    if(a_parent->getPDG() == 92){
                        finalFrom92.push_back(a_MCP);
                    }
                }
            }
            
            vquarkPx.clear(); vquarkPy.clear(); vquarkPz.clear(); vquarkE.clear();
            int num_quark = MCTQuark.size();
            cout<<"num_quark : "<<num_quark<<endl;
            for(int i = 0; i < num_quark; i++){
                MCParticle* a_MCP = MCTQuark.at(i);
                vquarkPx.push_back(a_MCP->getMomentum()[0]);
                vquarkPy.push_back(a_MCP->getMomentum()[1]);
                vquarkPz.push_back(a_MCP->getMomentum()[2]);
                vquarkE.push_back(a_MCP->getEnergy());
            }
            MCParticle* quark1 = MCTQuark.at(0);
            MCParticle* quark2 = MCTQuark.at(1);
            PDG_quark1 = quark1->getPDG();
            PDG_quark2 = quark2->getPDG();
            cout<<"PDG_quark1 : "<<PDG_quark1<<" PDG_quark2 : "<<PDG_quark2<<endl;
        
            
            
            
            vfinalPx.clear(); vfinalPy.clear(); vfinalPz.clear(); vfinalE.clear();
            int num_final = finalFrom92.size();
            for(int i = 0; i<num_final; i++){
                MCParticle* a_MCP = finalFrom92.at(i);
                vfinalPx.push_back(a_MCP->getMomentum()[0]);
                vfinalPy.push_back(a_MCP->getMomentum()[1]);
                vfinalPz.push_back(a_MCP->getMomentum()[2]);
                vfinalE.push_back(a_MCP->getEnergy());
            }
            
            
            cout<<"......."<<endl;
            
        } catch (lcio::DataNotAvailableException err) {  }
	
        _outputTree->Fill();
        Num ++;
	}  	  

}	



void softLink2quark::end()
{

	if (_outputTree) {

		TFile *tree_file = _outputTree->GetCurrentFile(); //just in case we switched to a new file
		//tree_file->cd();
		tree_file->Write();
		delete tree_file;
		//tree_file->Close();
	}

}

