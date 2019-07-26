#include <GetML.hh>
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
#include <fstream>
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

GetML a_GetML_instance;

GetML::GetML()
	: Processor("GetML"),
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

//ofstream fwccbs;
//ofstream fwccds;
//ofstream fwcuxx;
//ofstream fwuubd;
//ofstream fwuusd;
//ofstream fzccnots;
//ofstream fzdtdt;
//ofstream fzutut;
//ofstream fzuunotd;

void GetML::init() {

	printParameters();
	TFile *tree_file=new TFile(_treeFileName.c_str(),(_overwrite ? "RECREATE" : "UPDATE"));
	if (!tree_file->IsOpen()) {
		delete tree_file;
		tree_file=new TFile(_treeFileName.c_str(),"NEW");
	}

	_outputTree = new TTree(_treeName.c_str(),_treeName.c_str());
	_outputTree->SetAutoSave(32*1024*1024);  // autosave every 32MB
    _outputTree->Branch("Px_RecoJet1", &Px_RecoJet1,"Px_RecoJet1/F");
    _outputTree->Branch("Py_RecoJet1", &Py_RecoJet1,"Py_RecoJet1/F");
    _outputTree->Branch("Pz_RecoJet1", &Pz_RecoJet1,"Pz_RecoJet1/F");
    _outputTree->Branch("En_RecoJet1", &En_RecoJet1, "En_RecoJet1/F");
    
    _outputTree->Branch("Px_RecoJet2", &Px_RecoJet2,"Px_RecoJet2/F");
    _outputTree->Branch("Py_RecoJet2", &Py_RecoJet2,"Py_RecoJet2/F");
    _outputTree->Branch("Pz_RecoJet2", &Pz_RecoJet2,"Pz_RecoJet2/F");
    _outputTree->Branch("En_RecoJet2", &En_RecoJet2,"En_RecoJet2/F");
    
    _outputTree->Branch("Px_RecoJet3", &Px_RecoJet3,"Px_RecoJet3/F");
    _outputTree->Branch("Py_RecoJet3", &Py_RecoJet3,"Py_RecoJet3/F");
    _outputTree->Branch("Pz_RecoJet3", &Pz_RecoJet3,"Pz_RecoJet3/F");
    _outputTree->Branch("En_RecoJet3", &En_RecoJet3,"En_RecoJet3/F");
    
    _outputTree->Branch("Px_RecoJet4", &Px_RecoJet4,"Px_RecoJet4/F");
    _outputTree->Branch("Py_RecoJet4", &Py_RecoJet4,"Py_RecoJet4/F");
    _outputTree->Branch("Pz_RecoJet4", &Pz_RecoJet4,"Pz_RecoJet4/F");
    _outputTree->Branch("En_RecoJet4", &En_RecoJet4,"En_RecoJet4/F");
    
    _outputTree->Branch("Px_GenJet1", &Px_GenJet1,"Px_GenJet1/F");
    _outputTree->Branch("Py_GenJet1", &Py_GenJet1,"Py_GenJet1/F");
    _outputTree->Branch("Pz_GenJet1", &Pz_GenJet1,"Pz_GenJet1/F");
    _outputTree->Branch("En_GenJet1", &En_GenJet1,"En_GenJet1/F");
    
    _outputTree->Branch("Px_GenJet2", &Px_GenJet2,"Px_GenJet2/F");
    _outputTree->Branch("Py_GenJet2", &Py_GenJet2,"Py_GenJet2/F");
    _outputTree->Branch("Pz_GenJet2", &Pz_GenJet2,"Pz_GenJet2/F");
    _outputTree->Branch("En_GenJet2", &En_GenJet2,"En_GenJet2/F");
    
    _outputTree->Branch("Px_GenJet3", &Px_GenJet3,"Px_GenJet3/F");
    _outputTree->Branch("Py_GenJet3", &Py_GenJet3,"Py_GenJet3/F");
    _outputTree->Branch("Pz_GenJet3", &Pz_GenJet3,"Pz_GenJet3/F");
    _outputTree->Branch("En_GenJet3", &En_GenJet3,"En_GenJet3/F");
    
    _outputTree->Branch("Px_GenJet4", &Px_GenJet4,"Px_GenJet4/F");
    _outputTree->Branch("Py_GenJet4", &Py_GenJet4,"Py_GenJet4/F");
    _outputTree->Branch("Pz_GenJet4", &Pz_GenJet4,"Pz_GenJet4/F");
    _outputTree->Branch("En_GenJet4", &En_GenJet4,"En_GenJet4/F");
    
    _outputTree->Branch("Px_q1", &Px_q1,"Px_q1/F");
    _outputTree->Branch("Py_q1", &Py_q1,"Py_q1/F");
    _outputTree->Branch("Pz_q1", &Pz_q1,"Pz_q1/F");
    _outputTree->Branch("En_q1", &En_q1,"En_q1/F");
    
    _outputTree->Branch("Px_q2", &Px_q2,"Px_q2/F");
    _outputTree->Branch("Py_q2", &Py_q2,"Py_q2/F");
    _outputTree->Branch("Pz_q2", &Pz_q2,"Pz_q2/F");
    _outputTree->Branch("En_q2", &En_q2,"En_q2/F");
    
    _outputTree->Branch("Px_q3", &Px_q3,"Px_q3/F");
    _outputTree->Branch("Py_q3", &Py_q3,"Py_q3/F");
    _outputTree->Branch("Pz_q3", &Pz_q3,"Pz_q3/F");
    _outputTree->Branch("En_q3", &En_q3,"En_q3/F");
    
    _outputTree->Branch("Px_q4", &Px_q4,"Px_q4/F");
    _outputTree->Branch("Py_q4", &Py_q4,"Py_q4/F");
    _outputTree->Branch("Pz_q4", &Pz_q4,"Pz_q4/F");
    _outputTree->Branch("En_q4", &En_q4,"En_q4/F");

	Num = 0;
}

void GetML::processEvent( LCEvent * evtP )
{		

	if (evtP) 								
	{
        try{
            
            cout<<"************************************************************************"<<endl;


            eventNr = evtP->getEventNumber();
            cout<<"eventNr "<<eventNr<<" Num "<<Num<<endl;
            
    //for RecoJets
            cout<<".............................................RecoJet"<<endl;
            LCCollection* col_Jet = evtP->getCollection( "FastJets" );
            int nJet = col_Jet->getNumberOfElements();
            cout<<" the number of RecoJet is "<<nJet<<endl;
            if(nJet == 4)
            {
                ReconstructedParticle* a_RecoJet_1 = dynamic_cast<ReconstructedParticle*>(col_Jet->getElementAt(0));
                ReconstructedParticle* a_RecoJet_2 = dynamic_cast<ReconstructedParticle*>(col_Jet->getElementAt(1));
                ReconstructedParticle* a_RecoJet_3 = dynamic_cast<ReconstructedParticle*>(col_Jet->getElementAt(2));
                ReconstructedParticle* a_RecoJet_4 = dynamic_cast<ReconstructedParticle*>(col_Jet->getElementAt(3));
                
                Px_RecoJet1 = a_RecoJet_1->getMomentum()[0]; Py_RecoJet1 = a_RecoJet_1->getMomentum()[1];
                Pz_RecoJet1 = a_RecoJet_1->getMomentum()[2]; En_RecoJet1 = a_RecoJet_1->getEnergy();
                
                Px_RecoJet2 = a_RecoJet_2->getMomentum()[0]; Py_RecoJet2 = a_RecoJet_2->getMomentum()[1];
                Pz_RecoJet2 = a_RecoJet_2->getMomentum()[2]; En_RecoJet2 = a_RecoJet_2->getEnergy();
                
                Px_RecoJet3 = a_RecoJet_3->getMomentum()[0]; Py_RecoJet3 = a_RecoJet_3->getMomentum()[1];
                Pz_RecoJet3 = a_RecoJet_3->getMomentum()[2]; En_RecoJet3 = a_RecoJet_3->getEnergy();
                
                Px_RecoJet4 = a_RecoJet_4->getMomentum()[0]; Py_RecoJet4 = a_RecoJet_4->getMomentum()[1];
                Pz_RecoJet4 = a_RecoJet_4->getMomentum()[2]; En_RecoJet4 = a_RecoJet_4->getEnergy();
        
            }
            cout<<"En_RecoJet1+En_RecoJet4+En_RecoJet2+En_RecoJet3 = "<<En_RecoJet1+En_RecoJet4+En_RecoJet2+En_RecoJet3<<endl;
            
  //for MCParticle
            LCCollection* col_MCP = evtP->getCollection( "MCParticle" );
            int n_MCP = col_MCP->getNumberOfElements();
            std::vector<MCParticle* > MCTQuark;
            float En_ISR = 0;
            for(int i = 0; i<n_MCP; i++)
            {
                MCParticle* a_MCP = dynamic_cast<MCParticle*>(col_MCP->getElementAt(i));
                int NParents = a_MCP->getParents().size();
                int PDG = a_MCP->getPDG();
                TLorentzVector temp(a_MCP->getMomentum(), a_MCP->getEnergy());
                if(NParents == 0 && PDG != 22)
                {
                    MCTQuark.push_back(a_MCP);
                }
                if(NParents == 0 && PDG == 22)
                {
                    En_ISR += a_MCP->getEnergy();
                }
            }
            
            if(MCTQuark.size() == 4)
            {
                MCParticle* q1 = MCTQuark.at(0);
                MCParticle* q2 = MCTQuark.at(1);
                MCParticle* q3 = MCTQuark.at(2);
                MCParticle* q4 = MCTQuark.at(3);
                
                cout<<q1->getPDG()<<" "<<q2->getPDG()<<" "<<q3->getPDG()<<" "<<q4->getPDG()<<endl;
                
                Px_q1 = q1->getMomentum()[0]; Py_q1 = q1->getMomentum()[1];
                Pz_q1 = q1->getMomentum()[2]; En_q1 = q1->getEnergy();
                
                Px_q2 = q2->getMomentum()[0]; Py_q2 = q2->getMomentum()[1];
                Pz_q2 = q2->getMomentum()[2]; En_q2 = q2->getEnergy();
                
                Px_q3 = q3->getMomentum()[0]; Py_q3 = q3->getMomentum()[1];
                Pz_q3 = q3->getMomentum()[2]; En_q3 = q3->getEnergy();
                
                Px_q4 = q4->getMomentum()[0]; Py_q4 = q4->getMomentum()[1];
                Pz_q4 = q4->getMomentum()[2]; En_q4 = q4->getEnergy();

            }
            cout<<"En_q1 + En_q2 + En_q3 + En_q4 + En_ISR = "<<En_q1 + En_q2 + En_q3 + En_q4 + En_ISR<<endl;
            

            //for GenJets
            cout<<"............................................GenJet"<<endl;
            LCCollection* col_GJet = evtP->getCollection( "MCPFastJets" );
            int nGJet = col_GJet->getNumberOfElements();
            cout<<" the number of GenJet is "<<nGJet<<endl;
            if(nGJet == 4)
            {
                ReconstructedParticle* a_GenJet_1 = dynamic_cast<ReconstructedParticle*>(col_GJet->getElementAt(0));
                ReconstructedParticle* a_GenJet_2 = dynamic_cast<ReconstructedParticle*>(col_GJet->getElementAt(1));
                ReconstructedParticle* a_GenJet_3 = dynamic_cast<ReconstructedParticle*>(col_GJet->getElementAt(2));
                ReconstructedParticle* a_GenJet_4 = dynamic_cast<ReconstructedParticle*>(col_GJet->getElementAt(3));
            
                
                Px_GenJet1 = a_GenJet_1->getMomentum()[0]; Py_GenJet1 = a_GenJet_1->getMomentum()[1];
                Pz_GenJet1 = a_GenJet_1->getMomentum()[2]; En_GenJet1 = a_GenJet_1->getEnergy();
                
                Px_GenJet2 = a_GenJet_2->getMomentum()[0]; Py_GenJet2 = a_GenJet_2->getMomentum()[1];
                Pz_GenJet2 = a_GenJet_2->getMomentum()[2]; En_GenJet2 = a_GenJet_2->getEnergy();
                
                Px_GenJet3 = a_GenJet_3->getMomentum()[0]; Py_GenJet3 = a_GenJet_3->getMomentum()[1];
                Pz_GenJet3 = a_GenJet_3->getMomentum()[2]; En_GenJet3 = a_GenJet_3->getEnergy();
                
                Px_GenJet4 = a_GenJet_4->getMomentum()[0]; Py_GenJet4 = a_GenJet_4->getMomentum()[1];
                Pz_GenJet4 = a_GenJet_4->getMomentum()[2]; En_GenJet4 = a_GenJet_4->getEnergy();

            }
            cout<<"En_GenJet1 + En_GenJet2 + En_GenJet3 + En_GenJet4 = "<<En_GenJet1 + En_GenJet2 + En_GenJet3 + En_GenJet4<<endl;
            
        } catch (lcio::DataNotAvailableException err) {  }
	
        _outputTree->Fill();
        Num ++;
	}  	  

}	

void GetML::end()
{

	if (_outputTree) {

		TFile *tree_file = _outputTree->GetCurrentFile(); //just in case we switched to a new file
		//tree_file->cd();
		tree_file->Write();
		delete tree_file;
		//tree_file->Close();
	}

}



