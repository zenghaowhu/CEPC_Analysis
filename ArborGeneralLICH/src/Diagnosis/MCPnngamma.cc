#include <MCPnngamma.hh>
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

MCPnngamma a_MCPnngamma_instance;

MCPnngamma::MCPnngamma()
	: Processor("MCPnngamma"),
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

void MCPnngamma::init() {

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
    _outputTree->Branch("En_charge", &En_charge, "En_charge/D");
    _outputTree->Branch("En_neutron", &En_neutron, "En_neutron/D");
    _outputTree->Branch("En_ISR", &En_ISR, "En_ISR/D");
    _outputTree->Branch("En_neutrino", &En_neutrino, "En_neutrino/D");
    _outputTree->Branch("MCP_energy", &MCP_energy, "MCP_energy[100]/D");
    _outputTree->Branch("MCP_costheta", &MCP_costheta, "MCP_costheta[100]/D");
    _outputTree->Branch("gamma_size", &gamma_size,"gamma_size/I");
    
	Num = 0;
}

void MCPnngamma::processEvent( LCEvent * evtP )
{		

	if (evtP) 								
	{
        try{
            
       //     cout<<"************************************************************************"<<endl;


            eventNr = evtP->getEventNumber();
       //     cout<<"eventNr "<<eventNr<<" Num "<<Num<<endl;


    //for MCParticle
            TLorentzVector final_TL(0,0,0,0);
            TLorentzVector ISR_TL(0,0,0,0);
            TLorentzVector Neutrino_TL(0,0,0,0);
            TVector3 ISR_TV(0,0,0);
            TVector3 Neutrino_TV(0,0,0);
            TLorentzVector charge_TL(0,0,0,0);
            TLorentzVector gamma_TL(0,0,0,0);
            TLorentzVector neutron_TL(0,0,0,0);
            TLorentzVector barrel_TL(0,0,0,0);
            TLorentzVector overlap_TL(0,0,0,0);
            TLorentzVector endcap_TL(0,0,0,0);
            LCCollection* col_MCP = evtP->getCollection( "MCParticle" );
            int n_MCP = col_MCP->getNumberOfElements();
            int n_gamma = 0, n_neutron = 0, n_ISR = 0, n_charge = 0, n_total = 0, n_neutrino = 0;
            std::vector<MCParticle*> MCP_gamma;
            int count_low = 0;
            for(int i = 0; i<n_MCP; i++)
            {
                MCParticle* a_MCP = dynamic_cast<MCParticle*>(col_MCP->getElementAt(i));
                int NParents = a_MCP->getParents().size();
                int NDaughters = a_MCP->getDaughters().size();
                int PDG = a_MCP->getPDG();
                float charge = a_MCP->getCharge();
                TLorentzVector temp(a_MCP->getMomentum(), a_MCP->getEnergy());
                double costheta = temp.CosTheta();
                TVector3 temp_TV = a_MCP->getMomentum();

                if(a_MCP->getGeneratorStatus() == 1)
                {
                    n_total += 1;
                    final_TL += temp;
                    if(abs(PDG) == 12 || abs(PDG) == 14 || abs(PDG) == 16)
                    {
                        n_neutrino += 1;
                        Neutrino_TL += temp;
                    }
                    if(charge != 0){charge_TL += temp; n_charge +=1;}
                    if(charge == 0){neutron_TL += temp; n_neutron += 1;}
                    if(PDG == 22 && NParents != 0){
			if(fabs(temp.CosTheta()) < 0.01 && a_MCP->getEnergy() < 0.75)
			{count_low += 1;   //  printf("energy very los   %e, %e\n", a_MCP->getEnergy(),  temp.CosTheta());
				}
		    gamma_TL += temp; n_gamma += 1; MCP_gamma.push_back(a_MCP);}
                    if(PDG == 22 && NParents == 0){n_ISR += 1; ISR_TL += temp;}
                    
                }
                
            }
            if(count_low != 0)
		{
            cout<<"count_low : "<<count_low<<endl;
                }
	    gamma_size = MCP_gamma.size();
       //     cout<<"gamma_size "<<gamma_size<<endl;
            memset(MCP_costheta, 999.0, sizeof(MCP_costheta));
            memset(MCP_energy, 999.0, sizeof(MCP_energy));
            //MCP_energy[] = {999.0, 999.0, 999.0, 999.0, 999.0};
            for(int i = 0; i< gamma_size; i++)
            {
                MCParticle* gamma = MCP_gamma.at(i);
                TLorentzVector TL_gamma(gamma->getMomentum(), gamma->getEnergy());
                MCP_costheta[i] = TL_gamma.CosTheta();
                MCP_energy[i] = TL_gamma.E();
         //       cout<<"MCP_energy[i] "<<MCP_energy[i]<<endl;
        //        cout<<"MCP_costheta[i] "<<MCP_costheta[i]<<endl;
                       if(fabs(MCP_costheta[i]) < 0.01 && MCP_energy[i] < 0.75)
                        {
//printf("energy very los   %e, %e\n", MCP_energy[i],  MCP_costheta[i]);
			}

            }
            
            En_charge = charge_TL.E(); En_neutron = neutron_TL.E();
            En_ISR = ISR_TL.E();
            En_neutrino = Neutrino_TL.E();
       //     cout<<"............................................MCP"<<endl;
       //     cout<<"the energy of final state particle is "<<final_TL.E()<<endl;
      //      cout<<"the energy of ISR is "<<ISR_TL.E()<<endl;
        //    cout<<"the energy of neutrino is "<<Neutrino_TL.E()<<endl;
        //    cout<<"the energy of charge particles in MCP is "<<charge_TL.E()<<endl;
       //     cout<<"the energy of neutron particles in MCP is "<<neutron_TL.E()<<endl;
       //     cout<<"the energy of gamma in MCP is "<<gamma_TL.E()<<endl;
       //     cout<<"the total number of MCParticles is "<<n_total<<endl;
       //     cout<<"the number of neutrino particles is "<<n_neutrino<<endl;
       //     cout<<"the number of charge particles is "<<n_charge<<endl;
       //     cout<<"the number of neutron particles is "<<n_neutron<<endl;
       //     cout<<"the number of ISR is "<<n_ISR<<endl;
       //     cout<<"the number fo photon is "<<n_gamma<<endl;

        } catch (lcio::DataNotAvailableException err) {  }
	
        _outputTree->Fill();
        Num ++;
	}  	  

}	



void MCPnngamma::end()
{

	if (_outputTree) {

		TFile *tree_file = _outputTree->GetCurrentFile(); //just in case we switched to a new file
		//tree_file->cd();
		tree_file->Write();
		delete tree_file;
		//tree_file->Close();
	}

}



