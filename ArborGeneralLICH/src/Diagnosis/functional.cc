#include <functional.hh>
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
#include <TH1F.h>
#include <TVector3.h>
#include <TRandom.h>
#include <Rtypes.h> 
#include <sstream>		
#include <cmath>
#include <vector>
#include <TMath.h>
#include "TLorentzVector.h"
#include <UTIL/CellIDDecoder.h>

using namespace std;

const string ECALCellIDDecoder = "M:3,S-1:3,I:9,J:9,K-1:6";

functional a_functional_instance;

functional::functional()
	: Processor("functional"),
	_output(0)
{
	_description = "Print MC Truth" ;

	_treeFileName="MCTruth.root";
	registerProcessorParameter( "TreeOutputFile" , 
			"The name of the file to which the ROOT tree will be written" ,
			_treeFileName ,
			_treeFileName);
	
	_treeName="MCPart";
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

void functional::init() {

	printParameters();

	TFile *tree_file=new TFile(_treeFileName.c_str(),(_overwrite ? "RECREATE" : "UPDATE"));

	if (!tree_file->IsOpen()) {
		delete tree_file;
		tree_file=new TFile(_treeFileName.c_str(),"NEW");
	}

	_outputTree = new TTree(_treeName.c_str(),_treeName.c_str());
	_outputTree->SetAutoSave(32*1024*1024);  // autosave every 32MB
	_outputTree->Branch("EventNr", &_eventNr, "EventNr/I");
	_outputTree->Branch("Num", &_Num, "Num/I");
    
    //MCParticle
    _outputTree->Branch("ISREn", &ISREn, "ISREn/F");
    _outputTree->Branch("ISRPt", &ISRPt, "ISRPt/F");
    _outputTree->Branch("NPt_mcp", &NPt_mcp, "NPt_mcp/F");
    _outputTree->Branch("NEn_mcp", &NEn_mcp, "NEn_mcp/F");
    _outputTree->Branch("NEn", &NEn, "NEn/F");
    _outputTree->Branch("NPt", &NPt, "NPt/F");
    _outputTree->Branch("MaxJetCosTheta", &MaxJetCosTheta, "MaxJetCosTheta/F");
    _outputTree->Branch("HDir", &HDir, "HDir/F");
    _outputTree->Branch("OriQuarkID", &OriQuarkID, "OriQuarkID/I");
    _outputTree->Branch("type", &type, "type/I");
    _outputTree->Branch("PID1_decayedfromH", &PID1_decayedfromH, "PID1_decayedfromH/I");
    _outputTree->Branch("PID2_decayedfromH", &PID2_decayedfromH, "PID2_decayedfromH/I");
    _outputTree->Branch("PID1_decayedfromZ1", &PID1_decayedfromZ1, "PID1_decayedfromZ1/I");
    _outputTree->Branch("PID2_decayedfromZ1", &PID2_decayedfromZ1, "PID2_decayedfromZ1/I");
    _outputTree->Branch("PID1_decayedfromZ2", &PID1_decayedfromZ2, "PID1_decayedfromZ2/I");
    _outputTree->Branch("PID2_decayedfromZ2", &PID2_decayedfromZ2, "PID2_decayedfromZ2/I");
    _outputTree->Branch("mcp_Z1D1_energy", &mcp_Z1D1_energy, "mcp_Z1D1_energy/F");
    _outputTree->Branch("mcp_Z1D2_energy", &mcp_Z1D2_energy, "mcp_Z1D2_energy/F");
    _outputTree->Branch("mcp_Z2D1_energy", &mcp_Z2D1_energy, "mcp_Z2D1_energy/F");
    _outputTree->Branch("mcp_Z2D2_energy", &mcp_Z2D2_energy, "mcp_Z2D2_energy/F");
    _outputTree->Branch("mcp_muminus_energy", &mcp_muminus_energy, "mcp_muminus_energy/F");
    _outputTree->Branch("mcp_muplus_energy", &mcp_muplus_energy, "mcp_muplus_energy/F");
    _outputTree->Branch("mcp_HDplus_energy", &mcp_HDplus_energy, "mcp_HDplus_energy/F");
    _outputTree->Branch("mcp_HDminus_energy", &mcp_HDminus_energy, "mcp_HDminus_energy/F");
    _outputTree->Branch("mass_dimuon1_mcp", &mass_dimuon1_mcp, "mass_dimuon1_mcp/F");
    //ArborPFO
    _outputTree->Branch("NMuplusCount", &NMuplusCount, "NMuplusCount/I");
    _outputTree->Branch("NMuminusCount", &NMuminusCount, "NMuminusCount/I");
    _outputTree->Branch("NNeutronCount", &NNeutronCount, "NNeutronCount/I");
    _outputTree->Branch("energy_neutron_second", &energy_neutron_second, "energy_neutron_second/F");
    _outputTree->Branch("energy_neutron_first", &energy_neutron_first, "energy_neutron_first/F");
    _outputTree->Branch("Mass_neutron_value", &Mass_neutron_value, "Mass_neutron_value/F");
    _outputTree->Branch("energy_muminus", &energy_muminus, "energy_muminus/F");
    _outputTree->Branch("energy_muplus", &energy_muplus, "energy_muplus/F");
    _outputTree->Branch("energy_muminus_H", &energy_muminus_H, "energy_muminus_H/F");
    _outputTree->Branch("energy_muplus_H", &energy_muplus_H, "energy_muplus_H/F");
    _outputTree->Branch("angle_minus", &angle_minus, "angle_minus/F");
    _outputTree->Branch("angle_plus", &angle_plus, "angle_plus/F");
    _outputTree->Branch("H_angle_minus", &H_angle_minus, "H_angle_minus/F");
    _outputTree->Branch("H_angle_plus", &H_angle_plus, "H_angle_plus/F");
	_outputTree->Branch("Mass_a", &Mass_a, "Mass_a/F");
    _outputTree->Branch("mass_recoil1", &mass_recoil1, "mass_recoil1/F");
    _outputTree->Branch("mass_recoil2", &mass_recoil2, "mass_recoil2/F");
    _outputTree->Branch("mass_dilepton1_pfo", &mass_dilepton1_pfo, "mass_dilepton1_pfo/F");
    _outputTree->Branch("mass_dilepton2_pfo", &mass_dilepton2_pfo, "mass_dilepton2_pfo/F");
    _outputTree->Branch("missing_mass", &missing_mass, "missing_mass/F");
    _outputTree->Branch("visible_mass", &visible_mass, "visible_mass/F");

	_outputPFO = new TTree("PFO","PFO");
	_outputPFO->Branch("EventNr", &_eventNr, "EventNr/I");
    _outputPFO->Branch("Num", &_Num, "Num/I");


	_Num = 0;
}

struct CRP {
    
    bool operator()(ReconstructedParticle *p1, ReconstructedParticle *p2) const
    {
        return p1->getEnergy() > p2->getEnergy();
    }
    
};

struct CMP {
    
    bool operator()(MCParticle *p1, MCParticle *p2) const
    {
        return p1->getEnergy() > p2->getEnergy();
    }
    
};


void functional::processEvent( LCEvent * evtP )
{		

	if (evtP) 								
	{
        _eventNr = evtP->getEventNumber();

        TVector3 mcpT0(0,0,0);
        TVector3 mcpT1(0,0,0);
        TVector3 mcpT2(0,0,0);
        TVector3 mcpT3(0,0,0);
        TVector3 mcpT_HD1(0,0,0);
        TVector3 mcpT_HD2(0,0,0);
        TVector3 mcpT_Z1D1(0,0,0);
        TVector3 mcpT_Z1D2(0,0,0);
        TVector3 mcpT_Z2D1(0,0,0);
        TVector3 mcpT_Z2D2(0,0,0);
        TVector3 mcpT_ZDM(0,0,0);
        TVector3 mcpT_ZDP(0,0,0);
        
        TLorentzVector center(0,0,0,240);
        TLorentzVector mcp_muminus(0,0,0,0);
        TLorentzVector mcp_muplus(0,0,0,0);
        TLorentzVector ArborTotalP(0, 0, 0, 0);
        
        NEn = 0;
        NPt = 0;
        NEn_mcp = 0;
        float NEn_ori = 0;
        NPt_mcp = 0;
        float NPt_ori = 0;
        
        type = 0;
        PID1_decayedfromH = 0;
        PID2_decayedfromH = 0;
        PID1_decayedfromZ1 = 0;
        PID2_decayedfromZ2 = 0;
        PID2_decayedfromZ1 = 0;
        PID1_decayedfromZ2 = 0;
        
        mcp_Z1D1_energy = 0;
        mcp_Z1D2_energy = 0;
        mcp_Z2D1_energy = 0;
        mcp_Z2D2_energy = 0;
        mcp_HDplus_energy = 0;
        mcp_HDminus_energy = 0;
        mcp_muplus_energy = 0;
        mcp_muminus_energy = 0;
        mass_dimuon1_mcp = 0;
        
        NMuminusCount = 0;
        NMuplusCount = 0;
        energy_muplus = 0;
        energy_muminus = 0;
        energy_muplus_H = 0;
        energy_muminus_H = 0;
        
        NNeutronCount = 0;
        energy_neutron_first = 0;
        energy_neutron_second = 0;
        Mass_neutron_value = 0;
        

		OriQuarkID = 0;
		HDir = -10;
        
        MaxJetCosTheta = 10;
        float _J1CosTheta = 10;
        float _J2CosTheta = 10;
	
		try{
			LCCollection * MCP = evtP->getCollection("MCParticle");
			TVector3 tmpP; 
            TLorentzVector generator(0,0,0,0);
            
			for(int s0 = 0; s0 < MCP->getNumberOfElements(); s0++)
			{
				MCParticle *a1_MCP = dynamic_cast<EVENT::MCParticle *>(MCP->getElementAt(s0));
				int tmpPID = a1_MCP->getPDG();
				int NParent = a1_MCP->getParents().size();
				int NDaughter = a1_MCP->getDaughters().size();
				TVector3 VTX = a1_MCP->getVertex();
                TVector3 EndP = a1_MCP->getEndpoint();
                TLorentzVector mcpP( a1_MCP->getMomentum()[0], a1_MCP->getMomentum()[1], a1_MCP->getMomentum()[2], a1_MCP->getEnergy());

                if(NParent == 0 && s0<6){generator += mcpP;}
                
                if(tmpPID == 13 && NParent == 0 && s0 < 6)
                {
                    mcpT0 = a1_MCP->getMomentum();
                    mcp_muminus = mcpP;
                    mcp_muminus_energy = mcpP.E();
                }
                
                if(tmpPID == -13 && NParent == 0 && s0 < 6)
                {
                    mcpT1 = a1_MCP->getMomentum();
                    mcp_muplus = mcpP;
                    mcp_muplus_energy = mcpP.E();
                }
                
				if(tmpPID == 22 && NParent == 0 && s0 < 6)
				{
					tmpP = a1_MCP->getMomentum();
                    ISREn += tmpP.Mag();
                    ISRPt += tmpP.Perp();
				}

				if( (abs(tmpPID) == 12 || abs(tmpPID) == 14 || abs(tmpPID) == 16) && NDaughter == 0 )
				{
                    tmpP = a1_MCP->getMomentum();
                    NEn_mcp += tmpP.Mag();
                    NPt_mcp += tmpP.Perp();
                    
                    if(NParent == 0)
                    {
                        NEn_ori += tmpP.Mag();
                        NPt_ori += tmpP.Perp();
                    }

				}

				if(tmpPID == 25 && NDaughter > 1 && NParent !=0 ) //Higgs, in mctruth, there is a case Higgs decay to Higgs, so we ask NDaughter > 1
				{
					HDir = a1_MCP->getMomentum()[2]/a1_MCP->getEnergy();

                    if(NDaughter == 4)
                    {
                        cout<<"Higgs boson has two daughters........"<<endl;
                        PID1_decayedfromH = a1_MCP->getDaughters()[0]->getPDG();
                        PID2_decayedfromH = a1_MCP->getDaughters()[1]->getPDG();
                        mcp_HD1_energy = a1_MCP->getDaughters()[0]->getEnergy();
                        mcp_HD2_energy = a1_MCP->getDaughters()[1]->getEnergy();
                        mcpT_HD1 = a1_MCP->getDaughters()[0]->getMomentum();
                        mcpT_HD2 = a1_MCP->getDaughters()[1]->getMomentum();
                        MCParticle *D1 = a1_MCP->getDaughters()[0];
                        _J1CosTheta = D1->getMomentum()[2]/D1->getEnergy();
                        MCParticle *D2 = a1_MCP->getDaughters()[1];
                        _J2CosTheta = D2->getMomentum()[2]/D2->getEnergy();
                        type = PID1_decayedfromH;
                        
                        if(abs(PID1_decayedfromH) == 13)
                        {
                            if(PID1_decayedfromH == 13)
                            {mcpT2 = mcpT_HD1; mcp_HDminus_energy = mcp_HD1_energy;}
                            else if(PID1_decayedfromH == -13)
                            {mcpT3 = mcpT_HD1; mcp_HDplus_energy = mcp_HD1_energy;}
                            if(PID2_decayedfromH == 13)
                            {mcpT2 = mcpT_HD2; mcp_HDminus_energy = mcp_HD2_energy;}
                            else if(a1_MCP->getDaughters()[1]->getPDG() == -13)
                            {mcpT3 = mcpT_HD2; mcp_HDplus_energy = mcp_HD2_energy;}
                        }
                        
                        if(abs(a1_MCP->getDaughters()[0]->getPDG()) == 23)
                        {

                            PID1_decayedfromZ1 = a1_MCP->getDaughters()[0]->getDaughters()[0]->getPDG();
                            PID2_decayedfromZ1 = a1_MCP->getDaughters()[0]->getDaughters()[1]->getPDG();
                            PID1_decayedfromZ2 = a1_MCP->getDaughters()[1]->getDaughters()[0]->getPDG();
                            PID2_decayedfromZ2 = a1_MCP->getDaughters()[1]->getDaughters()[1]->getPDG();
                            
                            mcp_Z1D1_energy = a1_MCP->getDaughters()[0]->getDaughters()[0]->getEnergy();
                            mcpT_Z1D1 = a1_MCP->getDaughters()[0]->getDaughters()[0]->getMomentum();
                            
                            mcp_Z1D2_energy = a1_MCP->getDaughters()[0]->getDaughters()[1]->getEnergy();
                            mcpT_Z1D2 = a1_MCP->getDaughters()[0]->getDaughters()[1]->getMomentum();
                            
                            mcp_Z2D1_energy = a1_MCP->getDaughters()[1]->getDaughters()[0]->getEnergy();
                            mcpT_Z2D1 = a1_MCP->getDaughters()[1]->getDaughters()[0]->getMomentum();
                            
                            mcp_Z2D2_energy = a1_MCP->getDaughters()[1]->getDaughters()[1]->getEnergy();
                            mcpT_Z2D2 = a1_MCP->getDaughters()[1]->getDaughters()[1]->getMomentum();

                        }
                        
                    }
                    
				}
                
               
				if(abs(tmpPID)<7 && NParent == 0)
				{
                    OriQuarkID = abs(tmpPID);
				}

				if(abs(_J1CosTheta) > abs(_J2CosTheta))
					MaxJetCosTheta = _J1CosTheta;
				else
					MaxJetCosTheta = _J2CosTheta;

			}
            NEn = NEn_mcp - NEn_ori;
            NPt = NPt_mcp - NPt_ori;
            mass_dimuon1_mcp = (mcp_muplus +mcp_muminus).M();

            cout<<"generater energy is "<<generator.E()<<endl;
		}catch(lcio::DataNotAvailableException err) { }

        std::vector<ReconstructedParticle *> muminusList;
        std::vector<ReconstructedParticle *> muplusList;
        std::vector<ReconstructedParticle *> neutronList;
        
		try{
			LCCollection* col_RecoNeP = evtP->getCollection( "ArborPFOs" );

            
			for(int i0 = 0; i0 < col_RecoNeP->getNumberOfElements(); i0++)
			{
				ReconstructedParticle *a_RecoP = dynamic_cast<EVENT::ReconstructedParticle *>(col_RecoNeP->getElementAt(i0));	
				TLorentzVector currP( a_RecoP->getMomentum()[0], a_RecoP->getMomentum()[1], a_RecoP->getMomentum()[2], a_RecoP->getEnergy());
				ArborTotalP += currP;	
				TVector3 currMom = a_RecoP->getMomentum();


                if(a_RecoP->getCharge() == 0)
                {
                    neutronList.push_back(a_RecoP);
                }
                if(a_RecoP->getType() == 13)
                {
                    muminusList.push_back(a_RecoP);
                }
                if(a_RecoP->getType() == -13)
                {
                    muplusList.push_back(a_RecoP);
                }

                
				_outputPFO->Fill();

			}
		}catch (lcio::DataNotAvailableException err) { }


		Mass_a = 0;
        missing_mass = 0;
        visible_mass = 0;
        mass_dilepton1_pfo = 0;
        mass_dilepton2_pfo = 0;
        mass_recoil1 = 0;
        mass_recoil2 = 0;
        
        angle_minus = 999;
        angle_plus = 999;
        H_angle_plus = 999;
        H_angle_minus = 999;
        
        TLorentzVector first_neutron(0,0,0,0);
        TLorentzVector second_neutron(0,0,0,0);
        TLorentzVector muminus(0,0,0,0);
        TLorentzVector muplus(0,0,0,0);
        TLorentzVector Hmuminus(0,0,0,0);
        TLorentzVector Hmuplus(0,0,0,0);
        
        NNeutronCount = neutronList.size();
        if(neutronList.size() >= 2)
        {
            std::sort(neutronList.begin(), neutronList.end(), CRP());
            ReconstructedParticle *n1 = neutronList.at(0);
            ReconstructedParticle *n2 = neutronList.at(1);
            first_neutron.SetPxPyPzE(n1->getMomentum()[0], n1->getMomentum()[1], n1->getMomentum()[2], n1->getEnergy());
            second_neutron.SetPxPyPzE(n2->getMomentum()[0], n2->getMomentum()[1], n2->getMomentum()[2], n2->getEnergy());
            energy_neutron_first = n1->getEnergy();
            energy_neutron_second = n2->getEnergy();
        }
        Mass_neutron_value = (first_neutron + second_neutron).M();

        
        float angle1 = 999, angle2 = 999, angle11 = 999, angle22 = 999;
        int cor_minus1 = 999, cor_plus1 = 999, cor_minus2 = 999, cor_plus2 = 999;
        
        NMuminusCount = muminusList.size();
        NMuplusCount = muplusList.size();
        
        cout<<"NMuminusCount "<<NMuminusCount<<" NMuplusCount "<<NMuplusCount<<endl;
        
//        if(abs(PID1_decayedfromZ1) != 13 && abs(PID2_decayedfromZ2) != 13 && NMuplusCount >= 1 && NMuminusCount >= 1)
        if(NMuplusCount >= 1 && NMuminusCount >= 1)
        {
            for(int j = 0; j<NMuminusCount; j++)
            {
                ReconstructedParticle *n1 = muminusList.at(j);
                float recoT1 = mcpT0.Angle(n1->getMomentum());
                if(recoT1 < angle1){cor_minus1 = j; angle1 = recoT1;}
            }
            
            for(int k = 0; k<NMuplusCount; k++)
            {
                ReconstructedParticle *n2 = muplusList.at(k);
                float recoT2 = mcpT1.Angle(n2->getMomentum());
                if(recoT2 < angle2){cor_plus1 = k; angle2 = recoT2;}
            }
            
            
            ReconstructedParticle *n_minus = muminusList.at(cor_minus1);
            ReconstructedParticle *n_plus = muplusList.at(cor_plus1);
            energy_muminus = n_minus->getEnergy();
            energy_muplus = n_plus->getEnergy();
            muminus.SetPxPyPzE(n_minus->getMomentum()[0], n_minus->getMomentum()[1], n_minus->getMomentum()[2], n_minus->getEnergy());
            muplus.SetPxPyPzE(n_plus->getMomentum()[0], n_plus->getMomentum()[1], n_plus->getMomentum()[2], n_plus->getEnergy());
            Mass_a = (ArborTotalP - muminus - muplus).M();
            missing_mass = (center - ArborTotalP).M();
            visible_mass = ArborTotalP.M();
            mass_dilepton1_pfo = (muminus + muplus).M();
            mass_recoil1 = (center - muminus - muplus).M();
            cout<<"no muon decayed from Z boson "<<Mass_a<<" "<<mass_dilepton1_pfo<<" "<<mass_recoil1<<endl;
            cout<<"visible_mass "<<visible_mass<<" missing_mass "<<missing_mass<<endl;
        }
        
        
        if((abs(PID1_decayedfromZ1) == 13 || abs(PID2_decayedfromZ2) == 13) && NMuplusCount >= 2 && NMuminusCount >= 2)
        {
            if(PID1_decayedfromZ1 == 13){mcpT_ZDM = mcpT_Z1D1; mcpT_ZDP = mcpT_Z1D2;}
            else if(PID2_decayedfromZ1 == 13){mcpT_ZDM = mcpT_Z1D2; mcpT_ZDP = mcpT_Z1D1;}
            if(PID1_decayedfromZ2 == 13){mcpT_ZDM = mcpT_Z2D1; mcpT_ZDP = mcpT_Z2D2;}
            else if(PID2_decayedfromZ2 == 13){mcpT_ZDM = mcpT_Z2D2; mcpT_ZDP = mcpT_Z2D1;}
            
            for(int j=0; j<NMuminusCount; j++)
            {
                ReconstructedParticle *n1 = muminusList.at(j);
                float recoT1 = mcpT0.Angle(n1->getMomentum());
                if(recoT1 < angle1){cor_minus1 = j; angle1 = recoT1;}
            }
            
            for(int j = 0; j<NMuminusCount; j++)
            {
                if(j != cor_minus1)
                {
                ReconstructedParticle *n1 = muminusList.at(j);
                float recoT11 = mcpT_ZDM.Angle(n1->getMomentum());
                if(recoT11 < angle11){cor_minus2 = j; angle11 = recoT11;}
                }
            }
            for(int k=0; k<NMuplusCount; k++)
            {
                ReconstructedParticle *n2 = muplusList.at(k);
                float recoT2 = mcpT1.Angle(n2->getMomentum());
                if(recoT2 < angle2){cor_plus1 = k; angle2 = recoT2;}
            }
            
            for(int k=0; k<NMuplusCount; k++)
            {
                if(k != cor_plus1){
                ReconstructedParticle *n2 = muplusList.at(k);
                float recoT22 = mcpT_ZDP.Angle(n2->getMomentum());
                if(recoT22 < angle22){cor_plus2 = k; angle22 = recoT22;}
                }
            }
            
//            cout<<"............................................."<<endl;
            cout<<"there is muon decayed from Z boson"<<endl;
            cout<<cor_minus1<<" "<<cor_minus2<<" "<<cor_plus1<<" "<<cor_plus2<<endl;
            cout<<NMuplusCount<<" "<<NMuminusCount<<endl;
            cout<<angle1<<" "<<angle11<<" "<<angle2<<" "<<angle22<<endl;
            
            ReconstructedParticle *n_minus = muminusList.at(cor_minus1);
            ReconstructedParticle *n_plus = muplusList.at(cor_plus1);
            energy_muminus = n_minus->getEnergy();
            energy_muplus = n_plus->getEnergy();
            muminus.SetPxPyPzE(n_minus->getMomentum()[0], n_minus->getMomentum()[1], n_minus->getMomentum()[2], n_minus->getEnergy());
            muplus.SetPxPyPzE(n_plus->getMomentum()[0], n_plus->getMomentum()[1], n_plus->getMomentum()[2], n_plus->getEnergy());
            Mass_a = (ArborTotalP - muminus - muplus).M();
            mass_dilepton1_pfo = (muminus + muplus).M();
            mass_recoil1 = (center - muminus - muplus).M();
            
            ReconstructedParticle *nH_minus = muminusList.at(cor_minus2);
            ReconstructedParticle *nH_plus = muplusList.at(cor_plus2);
            energy_muminus_H = nH_minus->getEnergy();
            energy_muplus_H = nH_plus->getEnergy();
            Hmuminus.SetPxPyPzE(nH_minus->getMomentum()[0], nH_minus->getMomentum()[1], nH_minus->getMomentum()[2], nH_minus->getEnergy());
            Hmuplus.SetPxPyPzE(nH_plus->getMomentum()[0], nH_plus->getMomentum()[1], nH_plus->getMomentum()[2], nH_plus->getEnergy());
            mass_dilepton2_pfo = (Hmuplus + Hmuminus).M();
            mass_recoil2 = (center - muplus - muminus - Hmuminus - Hmuplus).M();
            
            TVector3 nT_minus =n_minus->getMomentum();
            TVector3 nT_plus =n_plus->getMomentum();
            TVector3 nHT_minus = nH_minus->getMomentum();
            TVector3 nHT_plus =nH_plus->getMomentum();
            cout<<nT_minus.Angle(nT_plus)<<" "<<mcpT0.Angle(mcpT1)<<endl;
            cout<<nHT_minus.Angle(nHT_plus)<<" "<<mcpT_ZDM.Angle(mcpT_ZDP)<<endl;
            
            angle_minus = nT_minus.Angle(mcpT0);
            angle_plus = nT_plus.Angle(mcpT1);
            H_angle_minus = nHT_minus.Angle(mcpT_ZDM);
            H_angle_plus = nHT_plus.Angle(mcpT_ZDP);
            
            cout<<angle_minus<<" "<<angle_plus<<" "<<H_angle_minus<<" "<<H_angle_plus<<endl;
            cout<<"......................................................"<<endl;
        }

		_outputTree->Fill();
		_Num++;
	}  	  

}	

void functional::end()
{

	if (_outputTree) {

		TFile *tree_file = _outputTree->GetCurrentFile(); //just in case we switched to a new file
		tree_file->Write();
		delete tree_file;
	}
}



