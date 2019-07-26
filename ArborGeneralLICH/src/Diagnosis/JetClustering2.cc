#include <JetClustering2.hh>
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

JetClustering2 a_JetClustering2_instance;

JetClustering2::JetClustering2()
	: Processor("JetClustering2"),
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

float TLORENTZVECTOR(std::vector<ReconstructedParticle * > slp);
int jetCharge(std::vector<ReconstructedParticle * > slp);
float costhetaSum(std::vector<ReconstructedParticle* > slp);
TLorentzVector getTLorentzVector(std::vector<ReconstructedParticle* > slp);
int NCharge(std::vector<ReconstructedParticle* > slp);
TLorentzVector getTLOfEachVector(std::vector<MCParticle* > vector);
void fillVector(std::vector<MCParticle* > slp,  std::vector<Float_t> &Px, std::vector<Float_t> &Py, std::vector<Float_t> &Pz, std::vector<Float_t> &E);

int ISBLong(int mc);
int ISBbar(int mc);
int ISB(int mc);
int ISS(int mc);
int ISC(int mc);
int ISMuon(int mc);
int ISElectron(int mc);

void JetClustering2::init() {

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
    _outputTree->Branch("Q1CosTheta", &Q1CosTheta, "Q1CosTheta/F");
    _outputTree->Branch("Q2CosTheta", &Q2CosTheta, "Q2CosTheta/F");
    _outputTree->Branch("total_charge", &total_charge, "total_charge/I");
    _outputTree->Branch("angletwoQ", &angletwoQ, "angletwoQ/F");
    _outputTree->Branch("PDG_Q1", &PDG_Q1, "PDG_Q1/I");
    _outputTree->Branch("PDG_Q2", &PDG_Q2, "PDG_Q2/I");
    _outputTree->Branch("sum1", &sum1, "sum1/F");
    _outputTree->Branch("sum2", &sum2, "sum2/F");
    _outputTree->Branch("CosThetaSum1", &CosThetaSum1, "CosThetaSum1/F");
    _outputTree->Branch("CosThetaSum2", &CosThetaSum2, "CosThetaSum2/F");
    _outputTree->Branch("jet1charge", &jet1charge, "jet1charge/I");
    _outputTree->Branch("jet2charge", &jet2charge, "jet2charge/I");
    _outputTree->Branch("num_jet1charge", &num_jet1charge, "num_jet1charge/I");
    _outputTree->Branch("num_jet2charge", &num_jet2charge, "num_jet2charge/I");
    _outputTree->Branch("alpha", &alpha, "alpha/F");
    _outputTree->Branch("alpha1", &alpha1, "alpha1/F");
    _outputTree->Branch("alpha2", &alpha2, "alpha2/F");
    _outputTree->Branch("Benergy", &Benergy, "Benergy/F");
    _outputTree->Branch("Bbarenergy", &Bbarenergy, "Bbarenergy/F");
    _outputTree->Branch("BCenergy", &BCenergy, "BCenergy/F");
    _outputTree->Branch("BbarCenergy", &BbarCenergy, "BbarCenergy/F");
    
    _outputTree->Branch("vBElcPx", &vBElcPx);
    _outputTree->Branch("vBElcPy", &vBElcPy);
    _outputTree->Branch("vBElcPz", &vBElcPz);
    _outputTree->Branch("vBElcE", &vBElcE);
    
    _outputTree->Branch("vBPosPx", &vBPosPx);
    _outputTree->Branch("vBPosPy", &vBPosPy);
    _outputTree->Branch("vBPosPz", &vBPosPz);
    _outputTree->Branch("vBPosE", &vBPosE);
    
    _outputTree->Branch("vBMuPx", &vBMuPx);
    _outputTree->Branch("vBMuPy", &vBMuPy);
    _outputTree->Branch("vBMuPz", &vBMuPz);
    _outputTree->Branch("vBMuE", &vBMuE);
    
    _outputTree->Branch("vBMuPlusPx", &vBMuPlusPx);
    _outputTree->Branch("vBMuPlusPy", &vBMuPlusPy);
    _outputTree->Branch("vBMuPlusPz", &vBMuPlusPz);
    _outputTree->Branch("vBMuPlusE", &vBMuPlusE);
    
    _outputTree->Branch("vBKPlusPx", &vBKPlusPx);
    _outputTree->Branch("vBKPlusPy", &vBKPlusPy);
    _outputTree->Branch("vBKPlusPz", &vBKPlusPz);
    _outputTree->Branch("vBKPlusE", &vBKPlusE);
    
    _outputTree->Branch("vBKMinusPx", &vBKMinusPx);
    _outputTree->Branch("vBKMinusPy", &vBKMinusPy);
    _outputTree->Branch("vBKMinusPz", &vBKMinusPz);
    _outputTree->Branch("vBKMinusE", &vBKMinusE);
    

    
    
    _outputTree->Branch("vBCElcPx", &vBCElcPx);
    _outputTree->Branch("vBCElcPy", &vBCElcPy);
    _outputTree->Branch("vBCElcPz", &vBCElcPz);
    _outputTree->Branch("vBCElcE", &vBCElcE);
    
    _outputTree->Branch("vBCPosPx", &vBCPosPx);
    _outputTree->Branch("vBCPosPy", &vBCPosPy);
    _outputTree->Branch("vBCPosPz", &vBCPosPz);
    _outputTree->Branch("vBCPosE", &vBCPosE);
    
    _outputTree->Branch("vBCMuPx", &vBCMuPx);
    _outputTree->Branch("vBCMuPy", &vBCMuPy);
    _outputTree->Branch("vBCMuPz", &vBCMuPz);
    _outputTree->Branch("vBCMuE", &vBCMuE);
    
    _outputTree->Branch("vBCMuPlusPx", &vBCMuPlusPx);
    _outputTree->Branch("vBCMuPlusPy", &vBCMuPlusPy);
    _outputTree->Branch("vBCMuPlusPz", &vBCMuPlusPz);
    _outputTree->Branch("vBCMuPlusE", &vBCMuPlusE);
    
    _outputTree->Branch("vBCKPlusPx", &vBCKPlusPx);
    _outputTree->Branch("vBCKPlusPy", &vBCKPlusPy);
    _outputTree->Branch("vBCKPlusPz", &vBCKPlusPz);
    _outputTree->Branch("vBCKPlusE", &vBCKPlusE);
    
    _outputTree->Branch("vBCKMinusPx", &vBCKMinusPx);
    _outputTree->Branch("vBCKMinusPy", &vBCKMinusPy);
    _outputTree->Branch("vBCKMinusPz", &vBCKMinusPz);
    _outputTree->Branch("vBCKMinusE", &vBCKMinusE);
    
    

    _outputTree->Branch("vBbarElcPx", &vBbarElcPx);
    _outputTree->Branch("vBbarElcPy", &vBbarElcPy);
    _outputTree->Branch("vBbarElcPz", &vBbarElcPz);
    _outputTree->Branch("vBbarElcE", &vBbarElcE);
    
    _outputTree->Branch("vBbarPosPx", &vBbarPosPx);
    _outputTree->Branch("vBbarPosPy", &vBbarPosPy);
    _outputTree->Branch("vBbarPosPz", &vBbarPosPz);
    _outputTree->Branch("vBbarPosE", &vBbarPosE);
    
    _outputTree->Branch("vBbarMuPx", &vBbarMuPx);
    _outputTree->Branch("vBbarMuPy", &vBbarMuPy);
    _outputTree->Branch("vBbarMuPz", &vBbarMuPz);
    _outputTree->Branch("vBbarMuE", &vBbarMuE);
    
    _outputTree->Branch("vBbarMuPlusPx", &vBbarMuPlusPx);
    _outputTree->Branch("vBbarMuPlusPy", &vBbarMuPlusPy);
    _outputTree->Branch("vBbarMuPlusPz", &vBbarMuPlusPz);
    _outputTree->Branch("vBbarMuPlusE", &vBbarMuPlusE);
    
    _outputTree->Branch("vBbarKPlusPx", &vBbarKPlusPx);
    _outputTree->Branch("vBbarKPlusPy", &vBbarKPlusPy);
    _outputTree->Branch("vBbarKPlusPz", &vBbarKPlusPz);
    _outputTree->Branch("vBbarKPlusE", &vBbarKPlusE);
    
    _outputTree->Branch("vBbarKMinusPx", &vBbarKMinusPx);
    _outputTree->Branch("vBbarKMinusPy", &vBbarKMinusPy);
    _outputTree->Branch("vBbarKMinusPz", &vBbarKMinusPz);
    _outputTree->Branch("vBbarKMinusE", &vBbarKMinusE);
    

    
    
    _outputTree->Branch("vBbarCElcPx", &vBbarCElcPx);
    _outputTree->Branch("vBbarCElcPy", &vBbarCElcPy);
    _outputTree->Branch("vBbarCElcPz", &vBbarCElcPz);
    _outputTree->Branch("vBbarCElcE", &vBbarCElcE);
    
    _outputTree->Branch("vBbarCPosPx", &vBbarCPosPx);
    _outputTree->Branch("vBbarCPosPy", &vBbarCPosPy);
    _outputTree->Branch("vBbarCPosPz", &vBbarCPosPz);
    _outputTree->Branch("vBbarCPosE", &vBbarCPosE);
    
    _outputTree->Branch("vBbarCMuPx", &vBbarCMuPx);
    _outputTree->Branch("vBbarCMuPy", &vBbarCMuPy);
    _outputTree->Branch("vBbarCMuPz", &vBbarCMuPz);
    _outputTree->Branch("vBbarCMuE", &vBbarCMuE);
    
    _outputTree->Branch("vBbarCMuPlusPx", &vBbarCMuPlusPx);
    _outputTree->Branch("vBbarCMuPlusPy", &vBbarCMuPlusPy);
    _outputTree->Branch("vBbarCMuPlusPz", &vBbarCMuPlusPz);
    _outputTree->Branch("vBbarCMuPlusE", &vBbarCMuPlusE);
    
    _outputTree->Branch("vBbarCKPlusPx", &vBbarCKPlusPx);
    _outputTree->Branch("vBbarCKPlusPy", &vBbarCKPlusPy);
    _outputTree->Branch("vBbarCKPlusPz", &vBbarCKPlusPz);
    _outputTree->Branch("vBbarCKPlusE", &vBbarCKPlusE);
    
    _outputTree->Branch("vBbarCKMinusPx", &vBbarCKMinusPx);
    _outputTree->Branch("vBbarCKMinusPy", &vBbarCKMinusPy);
    _outputTree->Branch("vBbarCKMinusPz", &vBbarCKMinusPz);
    _outputTree->Branch("vBbarCKMinusE", &vBbarCKMinusE);
    

    
    
    _outputTree->Branch("vOElcPx", &vOElcPx);
    _outputTree->Branch("vOElcPy", &vOElcPy);
    _outputTree->Branch("vOElcPz", &vOElcPz);
    _outputTree->Branch("vOElcE", &vOElcE);
    
    _outputTree->Branch("vOPosPx", &vOPosPx);
    _outputTree->Branch("vOPosPy", &vOPosPy);
    _outputTree->Branch("vOPosPz", &vOPosPz);
    _outputTree->Branch("vOPosE", &vOPosE);
    
    _outputTree->Branch("vOMuPx", &vOMuPx);
    _outputTree->Branch("vOMuPy", &vOMuPy);
    _outputTree->Branch("vOMuPz", &vOMuPz);
    _outputTree->Branch("vOMuE", &vOMuE);
    
    _outputTree->Branch("vOMuPlusPx", &vOMuPlusPx);
    _outputTree->Branch("vOMuPlusPy", &vOMuPlusPy);
    _outputTree->Branch("vOMuPlusPz", &vOMuPlusPz);
    _outputTree->Branch("vOMuPlusE", &vOMuPlusE);
    
    _outputTree->Branch("vOKPlusPx", &vOKPlusPx);
    _outputTree->Branch("vOKPlusPy", &vOKPlusPy);
    _outputTree->Branch("vOKPlusPz", &vOKPlusPz);
    _outputTree->Branch("vOKPlusE", &vOKPlusE);
    
    _outputTree->Branch("vOKMinusPx", &vOKMinusPx);
    _outputTree->Branch("vOKMinusPy", &vOKMinusPy);
    _outputTree->Branch("vOKMinusPz", &vOKMinusPz);
    _outputTree->Branch("vOKMinusE", &vOKMinusE);

    
	Num = 0;
}

void JetClustering2::processEvent( LCEvent * evtP )
{		

	if (evtP) 								
	{
        try{
            cout<<"************************************************************************"<<endl;
            eventNr = evtP->getEventNumber();
            cout<<"eventNr "<<eventNr<<" Num "<<Num<<endl;
            std::vector<MCParticle* > MCTQuark;
            LCCollection* col_MCP = evtP->getCollection( "MCParticle" );
            int n_MCP = col_MCP->getNumberOfElements();
            for(int i = 0; i<n_MCP; i++)
            {
                MCParticle* a_MCP = dynamic_cast<MCParticle*>(col_MCP->getElementAt(i));
                int NParents = a_MCP->getParents().size();
                int NDaughters = a_MCP->getDaughters().size();
                int PDG = a_MCP->getPDG();
                float charge = a_MCP->getCharge();
                TLorentzVector temp(a_MCP->getMomentum()[0], a_MCP->getMomentum()[1], a_MCP->getMomentum()[2], a_MCP->getEnergy());
                float costheta = temp.CosTheta();
                TVector3 temp_TV = a_MCP->getMomentum();
                
                
                if(NParents == 0 && ( abs(PDG) == 1 || abs(PDG) == 2 || abs(PDG) == 3 || abs(PDG) == 4 || abs(PDG) == 5 || abs(PDG) == 6))
                {
                    MCTQuark.push_back(a_MCP);
                }

            }
            
//            cout<<"MCTQuark.size() : "<<MCTQuark.size()<<endl;
            MCParticle* Q1 = MCTQuark.at(0);
            MCParticle* Q2 = MCTQuark.at(1);
            PDG_Q1 = Q1->getPDG();
            PDG_Q2 = Q2->getPDG();
            TLorentzVector TLQ1(Q1->getMomentum(), Q1->getEnergy());
            TLorentzVector TLQ2(Q2->getMomentum(), Q2->getEnergy());
            TVector3 TVQ1 = Q1->getMomentum();
            TVector3 TVQ2 = Q2->getMomentum();
            angletwoQ = TVQ1.Angle(TVQ2);
            Q1CosTheta = TLQ1.CosTheta();
            Q2CosTheta = TLQ2.CosTheta();
            cout<<"PDG_Q1 : PDG_Q2 "<<PDG_Q1<<" : "<<PDG_Q2<<endl;
//            cout<<"Q1CosTheta : Q2CosTheta "<<Q1CosTheta<<" : "<<Q2CosTheta<<endl;
            float halfCostheta = 0.5*(Q1CosTheta + Q2CosTheta);
            
            
            std::vector<MCParticle* >B;
            std::vector<MCParticle* >Bbar;
            
            std::vector<MCParticle* >BKPlus;
            std::vector<MCParticle* >BKMinus;
            std::vector<MCParticle* >BElectron;
            std::vector<MCParticle* >BPositron;
            std::vector<MCParticle* >BMuon;
            std::vector<MCParticle* >BMuonPlus;

            std::vector<MCParticle* >BbarKPlus;
            std::vector<MCParticle* >BbarKMinus;
            std::vector<MCParticle* >BbarElectron;
            std::vector<MCParticle* >BbarPositron;
            std::vector<MCParticle* >BbarMuon;
            std::vector<MCParticle* >BbarMuonPlus;
            
            std::vector<MCParticle* >BCMuon;
            std::vector<MCParticle* >BCMuonPlus;
            std::vector<MCParticle* >BCElectron;
            std::vector<MCParticle* >BCPositron;
            std::vector<MCParticle* >BCKPlus;
            std::vector<MCParticle* >BCKMinus;
            
            std::vector<MCParticle* >BbarCMuon;
            std::vector<MCParticle* >BbarCMuonPlus;
            std::vector<MCParticle* >BbarCElectron;
            std::vector<MCParticle* >BbarCPositron;
            std::vector<MCParticle* >BbarCKPlus;
            std::vector<MCParticle* >BbarCKMinus;
            
            
            
            //find out the particles original from b ro bbar
            for(int i = 0; i<n_MCP; i++)
            {
                MCParticle* a_MCP = dynamic_cast<MCParticle*>(col_MCP->getElementAt(i));
                int NParents = a_MCP->getParents().size();
                int NDaughters = a_MCP->getDaughters().size();
                
                if( abs(PDG_Q1) == 5 && abs(PDG_Q2) == 5 && NDaughters == 0 && NParents != 0)
                {
                    int PDG = 0, parent_PDG = 0;
                    MCParticle* a_parent = a_MCP;
                    do{
                        PDG = a_parent->getPDG();
                        a_parent = a_parent->getParents()[0];
                        parent_PDG = a_parent->getPDG();
//                        cout<<"PDG : parent_PDG "<<PDG<<" : "<<parent_PDG<<endl;
                    }
                    while (a_parent->getParents().size() != 0 && parent_PDG != 92);
                    if(parent_PDG == 92)
                    {
 //                       cout<<PDG<<endl;
                        if(ISB(PDG) == 1){B.push_back(a_MCP);}
                        else if(ISBbar(PDG) == 1){Bbar.push_back(a_MCP);}
                    }
                }
                
            }
            cout<<"B.size() : "<<B.size()<<" Bbar.size() : "<<Bbar.size()<<endl;

            
            
            
            // store e+ e- u+ u- K+ K-, which doesn't pass through c quark or b quark
            std::vector<MCParticle* >others;
            int finalWant = 0; //the final state particles which I want
            for(int i = 0; i<n_MCP; i++)
            {
                MCParticle* a_MCP = dynamic_cast<MCParticle*>(col_MCP->getElementAt(i));
                int NParents = a_MCP->getParents().size();
                int NDaughters = a_MCP->getDaughters().size();
                int PDG = a_MCP->getPDG();
                MCParticle* a_parent = a_MCP;
                int parent_PDG = 0;
                if((abs(PDG) == 11 || abs(PDG) == 13 || abs(PDG) == 321) && NDaughters == 0 && NParents != 0)
                {
                    finalWant += 1;
                    do{
                        a_parent = a_parent->getParents()[0];
                        parent_PDG = a_parent->getPDG();
                    }
                    while(ISC(abs(parent_PDG)) != 1 && ISBLong(abs(parent_PDG)) != 1 && abs(parent_PDG) != 92 && a_parent->getParents().size() != 0);
                    if(abs(parent_PDG) == 92 || a_parent->getParents().size() == 0)
                    {
                        others.push_back(a_MCP);
                    }
                }
            }
            
        
            cout<<"for particles which doesn't pass through b or bbar quark......"<<endl;
            std::vector<MCParticle* >otherElectron;
            std::vector<MCParticle* >otherPositron;
            std::vector<MCParticle* >otherMuon;
            std::vector<MCParticle* >otherMuonPlus;
            std::vector<MCParticle* >otherKPlus;
            std::vector<MCParticle* >otherKMinus;
            
            for(int i = 0; i<others.size(); i++)
            {
                MCParticle* a_MCP = others.at(i);
                int PDG = a_MCP->getPDG();
                cout<<"in others, PDG is "<<PDG<<endl;
                if(PDG == 11){otherElectron.push_back(a_MCP);}
                else if(PDG == -11){otherPositron.push_back(a_MCP);}
                else if(PDG == 13){otherMuon.push_back(a_MCP);}
                else if(PDG == -13){otherMuonPlus.push_back(a_MCP);}
                else if(PDG == 321){otherKPlus.push_back(a_MCP);}
                else if(PDG == -321){otherKMinus.push_back(a_MCP);}
            }
            
            vOElcPx.clear(); vOElcPy.clear(); vOElcPz.clear(); vOElcE.clear();     vOMuPx.clear(); vOMuPy.clear(); vOMuPz.clear(); vOMuE.clear();
            vOPosPx.clear(); vOPosPy.clear(); vOPosPz.clear(); vOPosE.clear();     vOMuPlusPx.clear(); vOMuPlusPy.clear(); vOMuPlusPz.clear(); vOMuPlusE.clear();
            vOKPlusPx.clear(); vOKPlusPy.clear(); vOKPlusPz.clear(); vOKPlusE.clear(); vOKMinusPx.clear(); vOKMinusPy.clear(); vOKMinusPz.clear(); vOKMinusE.clear();
            
            fillVector(otherElectron,  vOElcPx,    vOElcPy,    vOElcPz,    vOElcE);
            fillVector(otherPositron,  vOPosPx,    vOPosPy,    vOPosPz,    vOPosE);
            fillVector(otherMuon,      vOMuPx,     vOMuPy,     vOMuPz,     vOMuE);
            fillVector(otherMuonPlus,  vOMuPlusPx, vOMuPlusPy, vOMuPlusPz, vOMuPlusE);
            fillVector(otherKPlus,     vOKPlusPx,  vOKPlusPy,  vOKPlusPz,  vOKPlusE);
            fillVector(otherKMinus,    vOKMinusPx, vOKMinusPy, vOKMinusPz, vOKMinusE);

            
            
            
            //store the particles which pass through b quark
            cout<<"for particles which pass through b quark......."<<endl;
            int count = 0;
            
            for(int i = 0; i<B.size(); i++)
            {
                MCParticle* a_MCP = B.at(i);
                int PDG = a_MCP->getPDG();
                cout<<"in Bottom, PDG : "<<PDG<<endl;
                MCParticle* a_parent = a_MCP;
                int parent_PDG = 0;
                do{
                    a_parent = a_parent->getParents()[0];
                    parent_PDG = a_parent->getPDG();
                    cout<<"parent_PDG : "<<parent_PDG<<endl;
                }
                while( ISC( abs(parent_PDG) ) != 1 && ISB(parent_PDG) != 1 && a_parent->getParents().size() != 0);
                cout<<"parent_PDG : "<<parent_PDG<<endl;
                     if(PDG == 11   && ISB(parent_PDG) == 1){BElectron.push_back(a_MCP);}
                else if(PDG == -11  && ISB(parent_PDG) == 1){BPositron.push_back(a_MCP);}
                else if(PDG == 13   && ISB(parent_PDG) == 1){BMuon.push_back(a_MCP);}
                else if(PDG == -13  && ISB(parent_PDG) == 1){BMuonPlus.push_back(a_MCP);}
                else if(PDG == 11   && ISC(abs(parent_PDG)) == 1){BCElectron.push_back(a_MCP);}
                else if(PDG == -11  && ISC(abs(parent_PDG)) == 1){BCPositron.push_back(a_MCP);}
                else if(PDG == 13   && ISC(abs(parent_PDG)) == 1){BCMuon.push_back(a_MCP);}
                else if(PDG == -13  && ISC(abs(parent_PDG)) == 1){BCMuonPlus.push_back(a_MCP);}
                else if(PDG == 321  && ISB(parent_PDG) == 1){BKPlus.push_back(a_MCP);}
                else if(PDG == -321 && ISB(parent_PDG) == 1){BKMinus.push_back(a_MCP);}
                else if(PDG == 321  && ISC(abs(parent_PDG)) == 1){BCKPlus.push_back(a_MCP);}
                else if(PDG == -321 && ISC(abs(parent_PDG)) == 1){BCKMinus.push_back(a_MCP);}
                else cout<<"there is a special case !"<<endl;

                if(ISMuon(abs(PDG)) == 1 || ISElectron(abs(PDG)) == 1 || abs(PDG) == 321){count += 1;}
            }
            
            vBElcPx.clear(); vBElcPy.clear(); vBElcPz.clear(); vBElcE.clear();            vBPosPx.clear(); vBPosPy.clear(); vBPosPz.clear(); vBPosE.clear();
            vBMuPx.clear(); vBMuPy.clear(); vBMuPz.clear(); vBMuE.clear();                vBMuPlusPx.clear(); vBMuPlusPy.clear(); vBMuPlusPz.clear(); vBMuPlusE.clear();
            vBKPlusPx.clear(); vBKPlusPy.clear(); vBKPlusPz.clear(); vBKPlusE.clear();    vBKMinusPx.clear(); vBKMinusPy.clear(); vBKMinusPz.clear(); vBKMinusE.clear();
            
            fillVector(BElectron, vBElcPx,    vBElcPy,    vBElcPz,    vBElcE);
            fillVector(BPositron, vBPosPx,    vBPosPy,    vBPosPz,    vBPosE);
            fillVector(BMuon,     vBMuPx,     vBMuPy,     vBMuPz,     vBMuE);
            fillVector(BMuonPlus, vBMuPlusPx, vBMuPlusPy, vBMuPlusPz, vBMuPlusE);
            fillVector(BKPlus,    vBKPlusPx,  vBKPlusPy,  vBKPlusPz,  vBKPlusE);
            fillVector(BKMinus,   vBKMinusPx, vBKMinusPy, vBKMinusPz, vBKMinusE);

            vBCElcPx.clear();   vBCElcPy.clear();   vBCElcPz.clear();   vBCElcE.clear();      vBCPosPx.clear();    vBCPosPy.clear();    vBCPosPz.clear();    vBCPosE.clear();
            vBCMuPx.clear();    vBCMuPy.clear();    vBCMuPz.clear();    vBCMuE.clear();       vBCMuPlusPx.clear(); vBCMuPlusPy.clear(); vBCMuPlusPz.clear(); vBCMuPlusE.clear();
            vBCKPlusPx.clear(); vBCKPlusPy.clear(); vBCKPlusPz.clear(); vBCKPlusE.clear();    vBCKMinusPx.clear(); vBCKMinusPy.clear(); vBCKMinusPz.clear(); vBCKMinusE.clear();
            
            fillVector(BCElectron, vBCElcPx,    vBCElcPy,    vBCElcPz,    vBCElcE);
            fillVector(BCPositron, vBCPosPx,    vBCPosPy,    vBCPosPz,    vBCPosE);
            fillVector(BCMuon,     vBCMuPx,     vBCMuPy,     vBCMuPz,     vBCMuE);
            fillVector(BCMuonPlus, vBCMuPlusPx, vBCMuPlusPy, vBCMuPlusPz, vBCMuPlusE);
            fillVector(BCKPlus,    vBCKPlusPx,  vBCKPlusPy,  vBCKPlusPz,  vBCKPlusE);
            fillVector(BCKMinus,   vBCKMinusPx, vBCKMinusPy, vBCKMinusPz, vBCKMinusE);
            
            
            //store the particles which pass through bbar quark
            cout<<"for particles which pass through bbar quark......"<<endl;
            cout<<"Bbar*******"<<endl;
            
            int countBar = 0;
            for(int i = 0; i<Bbar.size(); i ++)
            {
                MCParticle* a_MCP = Bbar.at(i);
                int PDG = a_MCP->getPDG();
                cout<<"in Bbar, PDG : "<<PDG<<endl;
                MCParticle* a_parent = a_MCP;
                int parent_PDG = 0;
                do{
                    a_parent = a_parent->getParents()[0];
                    parent_PDG = a_parent->getPDG();
                    cout<<"parent_PDG : "<<parent_PDG<<endl;
                }
                while (ISC(abs(parent_PDG)) != 1 && ISBbar(parent_PDG) != 1 && a_parent->getParents().size() != 0);
                cout<<"parent_PDG : "<<parent_PDG<<endl;

                     if(PDG == 11   && ISBbar(parent_PDG) == 1){BbarElectron.push_back(a_MCP);}
                else if(PDG == -11  && ISBbar(parent_PDG) == 1){BbarPositron.push_back(a_MCP);}
                else if(PDG == 13   && ISBbar(parent_PDG) == 1){BbarMuon.push_back(a_MCP);}
                else if(PDG == -13  && ISBbar(parent_PDG) == 1){BbarMuonPlus.push_back(a_MCP);}
                else if(PDG == 11   && ISC(abs(parent_PDG)) == 1){BbarCElectron.push_back(a_MCP); cout<<"there is a BbarC"<<endl;}
                else if(PDG == -11  && ISC(abs(parent_PDG)) == 1){BbarCPositron.push_back(a_MCP);}
                else if(PDG == 13   && ISC(abs(parent_PDG)) == 1){BbarCMuon.push_back(a_MCP);}
                else if(PDG == -13  && ISC(abs(parent_PDG)) == 1){BbarCMuonPlus.push_back(a_MCP);}
                else if(PDG == 321  && ISBbar(parent_PDG) == 1){BbarKPlus.push_back(a_MCP);}
                else if(PDG == -321 && ISBbar(parent_PDG) == 1){BbarKMinus.push_back(a_MCP);}
                else if(PDG == 321  && ISC(abs(parent_PDG)) == 1){BbarCKPlus.push_back(a_MCP);}
                else if(PDG == -321 && ISC(abs(parent_PDG)) == 1){BbarCKMinus.push_back(a_MCP);}
                else cout<<"there is a special case !"<<endl;
                if(ISMuon(abs(PDG)) == 1 || ISElectron(abs(PDG)) == 1 || abs(PDG) == 321){countBar += 1;}

            }
            cout<<"BbarCElectron.size() : "<<BbarCElectron.size()<<endl;
            
            vBbarElcPx.clear();   vBbarElcPy.clear();   vBbarElcPz.clear();   vBbarElcE.clear();      vBbarPosPx.clear();    vBbarPosPy.clear();    vBbarPosPz.clear();    vBbarPosE.clear();
            vBbarMuPx.clear();    vBbarMuPy.clear();    vBbarMuPz.clear();    vBbarMuE.clear();       vBbarMuPlusPx.clear(); vBbarMuPlusPy.clear(); vBbarMuPlusPz.clear(); vBbarMuPlusE.clear();
            vBbarKPlusPx.clear(); vBbarKPlusPy.clear(); vBbarKPlusPz.clear(); vBbarKPlusE.clear();    vBbarKMinusPx.clear(); vBbarKMinusPy.clear(); vBbarKMinusPz.clear(); vBbarKMinusE.clear();
            
            fillVector(BbarElectron, vBbarElcPx,    vBbarElcPy,    vBbarElcPz,    vBbarElcE);
            fillVector(BbarPositron, vBbarPosPx,    vBbarPosPy,    vBbarPosPz,    vBbarPosE);
            fillVector(BbarMuon,     vBbarMuPx,     vBbarMuPy,     vBbarMuPz,     vBbarMuE);
            fillVector(BbarMuonPlus, vBbarMuPlusPx, vBbarMuPlusPy, vBbarMuPlusPz, vBbarMuPlusE);
            fillVector(BbarKPlus,    vBbarKPlusPx,  vBbarKPlusPy,  vBbarKPlusPz,  vBbarKPlusE);
            fillVector(BbarKMinus,   vBbarKMinusPx, vBbarKMinusPy, vBbarKMinusPz, vBbarKMinusE);

            vBbarCElcPx.clear();   vBbarCElcPy.clear();   vBbarCElcPz.clear();   vBbarCElcE.clear();      vBbarCPosPx.clear();    vBbarCPosPy.clear();    vBbarCPosPz.clear();   vBbarCPosE.clear();
            vBbarCMuPx.clear();    vBbarCMuPy.clear();    vBbarCMuPz.clear();    vBbarCMuE.clear();       vBbarCMuPlusPx.clear(); vBbarCMuPlusPy.clear(); vBbarCMuPlusPz.clear(); vBbarCMuPlusE.clear();
            vBbarCKPlusPx.clear(); vBbarCKPlusPy.clear(); vBbarCKPlusPz.clear(); vBbarCKPlusE.clear();    vBbarCKMinusPx.clear(); vBbarCKMinusPy.clear(); vBbarCKMinusPz.clear(); vBbarCKMinusE.clear();
            
            fillVector(BbarCElectron, vBbarCElcPx,    vBbarCElcPy,    vBbarCElcPz,    vBbarCElcE);
            fillVector(BbarCPositron, vBbarCPosPx,    vBbarCPosPy,    vBbarCPosPz,    vBbarCPosE);
            fillVector(BbarCMuon,     vBbarCMuPx,     vBbarCMuPy,     vBbarCMuPz,     vBbarCMuE);
            fillVector(BbarCMuonPlus, vBbarCMuPlusPx, vBbarCMuPlusPy, vBbarCMuPlusPz, vBbarCMuPlusE);
            fillVector(BbarCKPlus,    vBbarCKPlusPx,  vBbarCKPlusPy,  vBbarCKPlusPz,  vBbarCKPlusE);
            fillVector(BbarCKMinus,   vBbarCKMinusPx, vBbarCKMinusPy, vBbarCKMinusPz, vBbarCKMinusE);
            
            if(abs(PDG_Q1) == 5 && abs(PDG_Q2) == 5){
                cout<<"finalWant : others "<<finalWant<<" : "<<others.size()<<endl;
                cout<<"Electron, B : Bbar "<<BElectron.size()<<" : "<<BbarElectron.size()<<endl;
                cout<<"Positron, B : Bbar "<<BPositron.size()<<" : "<<BbarPositron.size()<<endl;
                cout<<"Muon, B : Bbar     "<<BMuon.size()<<" : "<<BbarMuon.size()<<endl;
                cout<<"MuonPlus, B : Bbar "<<BMuonPlus.size()<<" : "<<BbarMuonPlus.size()<<endl;
                cout<<"KPlus, B : Bbar    "<<BKPlus.size()<<" : "<<BbarKPlus.size()<<endl;
                cout<<"KMinus, B : Bbar   "<<BKMinus.size()<<" : "<<BbarKMinus.size()<<endl;
                
                
                cout<<"Electron, BC : BbarC "<<BCElectron.size()<<" : "<<BbarCElectron.size()<<endl;
                cout<<"Positron, BC : BbarC "<<BCPositron.size()<<" : "<<BbarCPositron.size()<<endl;
                cout<<"Muon, BC : BbarC     "<<BCMuon.size()<<" : "<<BbarCMuon.size()<<endl;
                cout<<"MuonPlus, BC : BbarC "<<BCMuonPlus.size()<<" : "<<BbarCMuonPlus.size()<<endl;
                cout<<"KPlus, BC : BbarC    "<<BCKPlus.size()<<" : "<<BbarCKPlus.size()<<endl;
                cout<<"KMinus, BC : BbarC   "<<BCKMinus.size()<<" : "<<BbarCKMinus.size()<<endl;
            }
            
            

            
            
            /*
            
            std::vector<ReconstructedParticle* > iterate1;    //pseudoJet, iterate among reconstruction
            std::vector<ReconstructedParticle* > iterate2;
            TLorentzVector TLONE(0,0,0,0);
            TLorentzVector TLTWO(0,0,0,0);
            TVector3 TVONE = TLONE.Vect();
            TVector3 TVTWO = TLTWO.Vect();
            
            std::vector<ReconstructedParticle* > half1;     //cheating, using the costheta of MCTruth quark
            std::vector<ReconstructedParticle* > half2;
            
            std::vector<ReconstructedParticle* > halfCostheta1;     //  half ,non-cheating, without iterating
            std::vector<ReconstructedParticle* > halfCostheta2;
            TVector3 first(0,0,0);
            TVector3 second(0,0,0);
            
            map<double, int> map_energy;
            LCCollection* col_PFO = evtP->getCollection( "ArborPFOs" );
            int n_PFO = col_PFO->getNumberOfElements();
            for(int i = 0; i<n_PFO; i++)
            {
                ReconstructedParticle* a_PFO = dynamic_cast<ReconstructedParticle*>(col_PFO->getElementAt(i));
                TLorentzVector temp(a_PFO->getMomentum(), a_PFO->getEnergy());
                double pfo_energy = temp.E();
                map_energy[pfo_energy] = i;
            }
            
            total_charge = 0;
            auto it = map_energy.rbegin(); float energy = it->first; int number = it->second;
            for(; it != map_energy.rend(); it++)
//            for(int i = 0; i<n_PFO; i++)
            {
                int i = it->second;
                ReconstructedParticle* a_PFO = dynamic_cast<ReconstructedParticle*>(col_PFO->getElementAt(i));
                int pfo_charge = a_PFO->getCharge();
                total_charge += pfo_charge;
                TLorentzVector temp(a_PFO->getMomentum(), a_PFO->getEnergy());
                float tempCostheta = temp.CosTheta();
                TVector3 TVtemp = a_PFO->getMomentum();
                //iterate among reconstruction
                if(i == 0){TVONE = TVtemp; TVTWO = -TVtemp;}
                if(TVtemp.Angle(TVONE) > TVtemp.Angle(TVTWO)){
                    TLTWO += temp;
                    TVTWO = TLTWO.Vect();
                    iterate2.push_back(a_PFO);
                }
                else if(TVtemp.Angle(TVONE) <= TVtemp.Angle(TVTWO)){
                    TLONE += temp;
                    TVONE = TLONE.Vect();
                    iterate1.push_back(a_PFO);
                }
                // half ,non-cheating, without iterating
                if(i == 0){first = TVtemp; second = -TVtemp;}
                if(TVtemp.Angle(first) < TVtemp.Angle(second)){
                    halfCostheta1.push_back(a_PFO);
                }
                else {halfCostheta2.push_back(a_PFO);}
                
                // cheating
                if(tempCostheta > halfCostheta){
                    half1.push_back(a_PFO);
                }
                else {half2.push_back(a_PFO);}
            }
            
            alpha = 999;
            alpha1 = 999, alpha2 = 999;
            
            TLorentzVector TLJet1(0,0,0,0);
            TLorentzVector TLJet2(0,0,0,0);
            TLJet1 = getTLorentzVector(iterate1);
            TLJet2 = getTLorentzVector(iterate2);
            
            TVector3 TVJet1 = TLJet1.Vect();
            TVector3 TVJet2 = TLJet2.Vect();
            
            float compare1 = TVQ1.Angle(TVJet1) * TVQ2.Angle(TVJet2);
            float compare2 = TVQ1.Angle(TVJet2) * TVQ2.Angle(TVJet1);
            
            std::vector<ReconstructedParticle* >jet1;
            std::vector<ReconstructedParticle* >jet2;
            
            if(compare1 < compare2){jet1 = iterate1; jet2 = iterate2; alpha = compare1; alpha1 = TVQ1.Angle(TVJet1); alpha2 = TVQ2.Angle(TVJet2);    }
            else if(compare2 <= compare1){jet1 = iterate2; jet2 = iterate1; alpha = compare2; alpha1 = TVQ1.Angle(TVJet2); alpha2 = TVQ2.Angle(TVJet1);  }
            
            
            sum1 = 0; sum2 = 0;
            jet1charge = 0; jet2charge = 0;
            CosThetaSum1 = 0; CosThetaSum2 = 0;
            num_jet1charge = 0; num_jet2charge = 0;
            
            num_jet1charge = NCharge(jet1);
            num_jet2charge = NCharge(jet2);
            
            sum1 = TLORENTZVECTOR(jet1);
            sum2 = TLORENTZVECTOR(jet2);
            
            jet1charge = jetCharge(jet1);
            jet2charge = jetCharge(jet2);

            CosThetaSum1 = costhetaSum(jet1);
            CosThetaSum2 = costhetaSum(jet2);
            
            cout<<"sum1 : sum2 "<<sum1<<" : "<<sum2<<endl;
            cout<<"jet1charge : jet2charge "<<jet1charge<<" : "<<jet2charge<<endl;
            */

    

        } catch (lcio::DataNotAvailableException err) {  }
	
        _outputTree->Fill();
        Num ++;
	}  	  

}	



void JetClustering2::end()
{

	if (_outputTree) {

		TFile *tree_file = _outputTree->GetCurrentFile(); //just in case we switched to a new file
		//tree_file->cd();
		tree_file->Write();
		delete tree_file;
		//tree_file->Close();
	}

}

int NCharge(std::vector<ReconstructedParticle* > slp){
    int num_charge = 0;
    for(int i = 0; i < slp.size(); i++)
    {
        ReconstructedParticle* pfo = slp.at(i);
        int pfo_charge = pfo->getCharge();
        if(pfo_charge != 0){num_charge += 1;}
    }
    return num_charge;
}

int jetCharge(std::vector<ReconstructedParticle * > slp){
    int total_charge = 0;
    for(int i = 0; i< slp.size(); i++)
    {
        ReconstructedParticle* pfo = slp.at(i);
        int pfo_charge = pfo->getCharge();
        total_charge += pfo_charge;
    }
    return total_charge;
}


TLorentzVector getTLorentzVector(std::vector<ReconstructedParticle* > slp){
    TLorentzVector TLtotal(0,0,0,0);
    for(int i = 0; i < slp.size(); i++)
    {
        ReconstructedParticle* pfo = slp.at(i);
        TLorentzVector temp(pfo->getMomentum(), pfo->getEnergy());
        TLtotal += temp;
    }
    return TLtotal;
}

float costhetaSum(std::vector<ReconstructedParticle* > slp){
    float costhetasum = 0;
    TLorentzVector TLtotal(0,0,0,0);
    for(int i = 0; i< slp.size(); i++)
    {
        ReconstructedParticle* pfo = slp.at(i);
        TLorentzVector temp(pfo->getMomentum(),  pfo->getEnergy());
        TLtotal += temp;
    }
    TVector3 TVtotal = TLtotal.Vect();
    for(int i = 0; i< slp.size(); i++)
    {
        ReconstructedParticle* pfo = slp.at(i);
        TLorentzVector temp(pfo->getMomentum(),  pfo->getEnergy());
        TVector3 TVtemp = temp.Vect();
        costhetasum += fabs(cos(TVtotal.Angle(TVtemp)));
    }
    return costhetasum;
}


float TLORENTZVECTOR(std::vector<ReconstructedParticle * > slp){
    TLorentzVector TLtotal(0,0,0,0);
    for(int i = 0; i< slp.size(); i++)
    {
        ReconstructedParticle* pfo = slp.at(i);
        TLorentzVector temp(pfo->getMomentum(),  pfo->getEnergy());
        TLtotal += temp;
    }
    TVector3 TVtotal = TLtotal.Vect();
    float sum = 0;
    for(int i = 0; i< slp.size(); i++)
    {
        ReconstructedParticle* pfo = slp.at(i);
        int charge = pfo->getCharge();
        TLorentzVector temp(pfo->getMomentum(),  pfo->getEnergy());
        TVector3 TVtemp = temp.Vect();
        sum += charge * fabs(cos(TVtotal.Angle(TVtemp)));
    }
    return sum;
}


int ISBbar(int mc){
    if(mc == 511 || mc == 521 || mc == 10511 || mc == 10521 || mc == 513 || mc == 523 || mc == 10513 || mc == 10523 || mc == 20513 || mc == 20523 || mc == 515 || mc == 525 || mc == 531 || mc == 10531 || mc == 533 || mc == 10533 || mc == 20533 || mc == 535 || mc == 541 || mc == 10541 || mc == 543 || mc == 10543 || mc == 20543 || mc == 545 || mc == -5122 || mc == -5112 || mc == -5212 || mc == -5222 || mc == -5114 || mc == -5214 || mc == -5224 || mc == -5132 || mc == -5232 || mc == -5312 || mc == -5322 || mc == -5314 || mc == -5324 || mc == -5332 || mc == -5334 || mc == -5142 || mc == -5242 || mc == -5412 || mc == -5422 || mc == -5414 || mc == -5424 || mc == -5342 || mc == -5432 || mc == -5434 || mc == -5442 || mc == -5444 || mc == -5512 || mc == -5522 || mc == -5514 || mc == -5524 || mc == -5532 || mc == -5534 || mc == -5542 || mc == -5544 || mc == -5544 ) {return 1;}
    else {return 0;}
}

int ISB(int mc){
    if(mc == -511 || mc == -521 || mc == -10511 || mc == -10521 || mc == -513 || mc == -523 || mc == -10513 || mc == -10523 || mc == -20513 || mc == -20523 || mc == -515 || mc == -525 || mc == -531 || mc == -10531 || mc == -533 || mc == -10533 || mc == -20533 || mc == -535 || mc == -541 || mc == -10541 || mc == -543 || mc == -10543 || mc == -20543 || mc == -545 || mc == 5122 || mc == 5112 || mc == 5212 || mc == 5222 || mc == 5114 || mc == 5214 || mc == 5224 || mc == 5132 || mc == 5232 || mc == 5312 || mc == 5322 || mc == 5314 || mc == 5324 || mc == 5332 || mc == 5334 || mc == 5142 || mc == 5242 || mc == 5412 || mc == 5422 || mc == 5414 || mc == 5424 || mc == 5342 || mc == 5432 || mc == 5434 || mc == 5442 || mc == 5444 || mc == 5512 || mc == 5522 || mc == 5514 || mc == 5524 || mc == 5532 || mc == 5534 || mc == 5542 || mc == 5544 || mc == 5544 ) {return 1;}
    else {return 0;}

}

int ISBLong(int mc){
    if(mc == 511 || mc == 521 || mc == 531 || mc == 541 || mc == 5112 || mc == 5122 || mc == 5132 || mc == 5232 || mc == 5332)
    {return 1;}
    else {return 0;}

}

int ISC(int mc){
    if(mc == 411 || mc == 421 || mc == 10411 || mc == 10421 || mc == 413 || mc == 423 || mc == 10413 || mc == 10423 || mc == 20413 || mc == 20423 || mc == 415 || mc == 425 || mc == 431 || mc == 10431 || mc == 433 || mc == 10433 || mc == 20433 || mc == 435      || mc == 4122 || mc == 4222 || mc == 4212 || mc == 4112 || mc == 4224 || mc == 4214 || mc == 4114 || mc == 4232 || mc == 4132 || mc == 4322 || mc == 4312 || mc == 4324 || mc == 4314 || mc == 4332 || mc == 4334 || mc == 4412 || mc == 4422 || mc == 4414 || mc == 4424 || mc == 4432 || mc == 4434 || mc == 4444) {return 1;}
    else {return 0;}
}

int ISS(int mc){
    if(mc == 130 || mc == 310 || mc == 311 || mc == 321 || mc == 10311 || mc == 10321 || mc == 100311 || mc == 100321 || mc == 200311 || mc == 200321 || mc == 9000311 || mc == 9000321 || mc == 313 || mc == 323 || mc == 10313 || mc == 10323 || mc == 20313 || mc == 20323 || mc == 100313 || mc == 100323 || mc == 9000313 || mc == 9000323 || mc == 30313 || mc == 30323 || mc == 315 || mc == 325 || mc == 9000315 || mc == 9000325 || mc == 10315 || mc == 10325 || mc == 20315 || mc == 20325 || mc == 100315 || mc == 100325 || mc == 9010315 || mc == 9010325 || mc == 317 || mc == 327 || mc == 9010317 || mc == 9010327 || mc == 319 || mc == 329 || mc == 9000319 || mc == 9000329       || mc == 3122 || mc == 3222 || mc == 3212 || mc == 3112 || mc == 3224 || mc == 3214 || mc == 3114 || mc == 3322 || mc == 3312 || mc == 3324 || mc == 3314 || mc == 3334) {return 1;}
    else {return 0;}

}

int ISMuon(int mc){
    if(abs(mc) == 13) {return 1;}
    else {return 0;}
}

int ISElectron(int mc){
    if(abs(mc) == 11) {return 1;}
    else {return 0;}
}


TLorentzVector getTLOfEachVector(std::vector<MCParticle* > vector){
    TLorentzVector want(0,0,0,0);
    if(vector.size() != 0)
    {
        for(int i = 0; i<vector.size(); i++)
        {
            MCParticle* mcp = vector.at(i);
            cout<<mcp->getPDG()<<endl;
            TLorentzVector temp(mcp->getMomentum(), mcp->getEnergy());
            want += temp;
        }
    }
    return want;
}


void fillVector(std::vector<MCParticle* > slp,  std::vector<Float_t> &Px, std::vector<Float_t> &Py, std::vector<Float_t> &Pz, std::vector<Float_t> &E){
    if(slp.size() != 0)
    {
        cout<<"the number of particles in this vector is "<<slp.size()<<endl;
        for(int i = 0; i<slp.size(); i++)
        {
            MCParticle* a_MCP = slp.at(i);
            cout<<"it's not empty"<<endl;
            Px.push_back(a_MCP->getMomentum()[0]);
            Py.push_back(a_MCP->getMomentum()[1]);
            Pz.push_back(a_MCP->getMomentum()[2]);
            E.push_back(a_MCP->getEnergy());
        }
        cout<<E.size()<<endl;
    }
}
