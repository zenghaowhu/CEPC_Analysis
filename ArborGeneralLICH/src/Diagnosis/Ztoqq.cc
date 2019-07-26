#include <Ztoqq.hh>
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

Ztoqq a_Ztoqq_instance;

Ztoqq::Ztoqq()
	: Processor("Ztoqq"),
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

void Ztoqq::init() {

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
    _outputTree->Branch("CParameter", &CParameter, "CParameter/D");
    _outputTree->Branch("hemiMass1", &hemiMass1, "hemiMass1/D");
    _outputTree->Branch("hemiMass2", &hemiMass2, "hemiMass2/D");
    _outputTree->Branch("hemiBroadening1", &hemiBroadening1, "hemiBroadening1/D");
    _outputTree->Branch("hemiBroadening2", &hemiBroadening2, "hemiBroadening2/D");
    _outputTree->Branch("T", &T, "T/D");
    _outputTree->Branch("angleTwoQ", &angleTwoQ, "angleTwoQ/D");
    _outputTree->Branch("En_ISR", &En_ISR, "En_ISR/D");
    _outputTree->Branch("Q1BackThrust", &Q1BackThrust, "Q1BackThrust/D");
    _outputTree->Branch("Q1Thrust", &Q1Thrust, "Q1Thrust/D");
    _outputTree->Branch("Q1PDG", &Q1PDG, "Q1PDG/I");
    _outputTree->Branch("Q2PDG", &Q2PDG, "Q2PDG/I");
    _outputTree->Branch("En_twoQ", &En_twoQ, "En_twoQ/D");
	Num = 0;
}

void Ztoqq::processEvent( LCEvent * evtP )
{		

	if (evtP) 								
	{
        try{
            
            

            eventNr = evtP->getEventNumber();
            cout<<"**************************************"<<endl;
            cout<<"eventNr "<<eventNr<<" Num "<<Num<<endl;

            
            LCCollection* col_PFO = evtP->getCollection( "ArborPFOs" );
            int nPFO = col_PFO->getNumberOfElements();

            TLorentzVector Roverlap_TL(0,0,0,0);
            
            for(int i = 0; i<nPFO; i++)
            {
                ReconstructedParticle* a_Reco = dynamic_cast<ReconstructedParticle*>(col_PFO->getElementAt(i));
                TLorentzVector temp(a_Reco->getMomentum(), a_Reco->getEnergy());
                float Rcostheta = temp.CosTheta();
                int type = a_Reco->getType();
                int Rcharge = a_Reco->getCharge();
            }


            LCCollection* col_MCP = evtP->getCollection( "MCParticle" );
            int n_MCP = col_MCP->getNumberOfElements();
            std::vector<MCParticle* >FinalStateParticle;
            std::vector<MCParticle* >MCPGamma;
            std::vector<MCParticle* >MCPGammaExcISR;
            std::vector<MCParticle* >ISR;
            std::vector<MCParticle* >MCPNeutrino;
            std::vector<MCParticle* >MCPCharge;
            std::vector<MCParticle* >MCPNeutron;
            std::vector<MCParticle* >MCPVisibleParticle;
            std::vector<MCParticle* >MCPQuark;
            std::vector<MCParticle* >MCPVisibleParticleExcISR;
            std::vector<MCParticle* >MCPUsedForThrust;
            
            FinalStateParticle.clear();
            MCPCharge.clear();
            MCPNeutron.clear();
            
            TLorentzVector TLFinalStateParticle(0,0,0,0);
            TLorentzVector TLMCPGamma(0,0,0,0);
            TLorentzVector TLMCPGammaExcISR(0,0,0,0);
            TLorentzVector TLISR(0,0,0,0);
            TLorentzVector TLMCPNeutrino(0,0,0,0);
            TLorentzVector TLMCPCharge(0,0,0,0);
            TLorentzVector TLMCPNeutron(0,0,0,0);
            TLorentzVector TLMCPVisibleParticle(0,0,0,0);
            TLorentzVector TLMCPVisibleParticleExcISR(0,0,0,0);
    
            
            int num_charge = 0, num_neutron = 0, num_neutrino = 0;
            
            for(int i = 0; i < n_MCP; i++){
                MCParticle* a_MCP = dynamic_cast<MCParticle*>(col_MCP->getElementAt(i));
                int NParents = a_MCP->getParents().size();
                int NDaughters = a_MCP->getDaughters().size();
                int PDG = a_MCP->getPDG();
                double charge = a_MCP->getCharge();
                TLorentzVector TLMCPtemp(a_MCP->getMomentum(), a_MCP->getEnergy());
                double costheta = TLMCPtemp.CosTheta();
                TVector3 TVMCPtemp = TLMCPtemp.Vect();
                
                
                if(NParents == 0){
                    if(abs(PDG) == 1 || abs(PDG) == 2 || abs(PDG) == 3 || abs(PDG) == 4 || abs(PDG) == 5 || abs(PDG) == 6){
                        MCPQuark.push_back(a_MCP);             //the original quarks
                    }
                    
                    if(PDG == 22){                    //ISR
                        TLISR += TLMCPtemp;
                        ISR.push_back(a_MCP);
                    }
                }
                
                if(a_MCP->getGeneratorStatus() == 1){
                    TLFinalStateParticle += TLMCPtemp;
                    FinalStateParticle.push_back(a_MCP);    //final state particles include ISR and neutrinos， mainly used for checking
                    
                    if(PDG != 22){MCPUsedForThrust.push_back(a_MCP);} //used for store the final state particles doesn't decay from ISR
                    else if(PDG == 22){
                        MCParticle* b_parent = a_MCP;
                        do{b_parent = b_parent->getParents()[0];}
                        while(b_parent->getParents().size() != 0);
                        if(b_parent->getParents().size() == 0 && (b_parent->getPDG()) != 22){
                            MCPUsedForThrust.push_back(a_MCP);
                        }
                    }
                    
                    
                    if(charge == 0 && abs(PDG) != 12 && abs(PDG) != 14 && abs(PDG) != 16){MCPNeutron.push_back(a_MCP);} //store neutron particles include gamma exclude neutrino
                    
                    if(charge != 0){MCPCharge.push_back(a_MCP);}
                    
                    if(abs(PDG) == 12 || abs(PDG) == 14 || abs(PDG) == 16){
                        TLMCPNeutrino += TLMCPtemp;
                        MCPNeutrino.push_back(a_MCP);
                    }
                    
                    else {MCPVisibleParticle.push_back(a_MCP);  TLMCPVisibleParticle += TLMCPtemp;} //store all visible final state particles invlude ISR
                    
                    if(PDG != 22 && abs(PDG) != 12 && abs(PDG) != 14 && abs(PDG) != 16){MCPVisibleParticleExcISR.push_back(a_MCP);}
                    MCParticle* a_parent = a_MCP;
                    do{a_parent = a_parent->getParents()[0];}
                    while(a_parent->getParents().size() != 0);
                    if(PDG == 22 && a_parent->getPDG() != 22 ){ // this photon is not a ISR
                        MCPVisibleParticleExcISR.push_back(a_MCP);
                    }
                }
            }
            
            cout<<"FinalStateParticle.size() : "<<FinalStateParticle.size()<<" MCPUsedForThrust.size() : "<<MCPUsedForThrust.size()<<endl;
            
            // stydy original quarks
            MCParticle* quark1 = MCPQuark.at(0);
            MCParticle* quark2 = MCPQuark.at(1);
            Q1PDG = quark1->getPDG();
            Q2PDG = quark2->getPDG();
            TLorentzVector TLquark1(quark1->getMomentum(), quark1->getEnergy());
            TLorentzVector TLquark2(quark2->getMomentum(), quark2->getEnergy());
            TVector3 TVquark1 = TLquark1.Vect();
            TVector3 TVquark2 = TLquark2.Vect();
            En_twoQ = (TLquark1 + TLquark2).E();
            angleTwoQ = TVquark1.Angle(TVquark2);
            cout<<"the angle between two quarks is "<<TVquark1.Angle(TVquark2)<<endl;
            cout<<"the energy of quark1 : quark2 "<<TLquark1.E()<<" : "<<TLquark2.E()<<endl;
            cout<<"the mass of quark1 : quark2 "<<TLquark1.M()<<" : "<<TLquark2.M()<<endl;
            
            
            //the following code used to find the thrust
            double thetaMin = TMath::Pi(), phiMin = 2*TMath::Pi();
            double thetaMin2 = 0, phiMin2 = 0;
            double thetaL = 0, phiL = 0, thetaR = TMath::Pi(), phiR = 2*TMath::Pi();
            int iter = 0;
            double Told = 0;
            double Tnew = 0;
            double cut = 1;
            double thetaRange = 0, phiRange = 0;
            do{
                iter += 1;
                cout<<"iter : "<<iter<<endl;
                if(iter == 1){
                    thetaRange = thetaR - thetaL, phiRange = phiR - phiL;
                }
                else if(iter != 1){
                    
                    thetaRange = 0.1*(thetaR - thetaL);
                    phiRange = 0.1*(phiR - phiL);
                    
                    thetaL =  thetaMin - thetaRange;
                    thetaR = thetaMin + thetaRange;
                    phiL = phiMin - phiRange;
                    phiR = phiMin + phiRange;
                    thetaRange = thetaR - thetaL, phiRange = phiR - phiL;

                    cout<<"thetaL : "<<thetaL<<" thetaR : "<<thetaR<<endl;
                    cout<<"phiL : "<<phiL<<" phiR : "<<phiR<<endl;
                }
                
                cout<<"thetaRange : "<<thetaRange<<" phiRange : "<<phiRange<<endl;
                for(double theta = thetaL; theta <= thetaR; theta += 0.1*thetaRange){   //in this round, find the max T
                    for(double phi = phiL; phi <= phiR; phi += 0.1*phiRange){
                        
                        double x = sin(theta)*cos(phi);
                        double y = sin(theta)*sin(phi);
                        double z = cos(theta);

                        double denominator = 0;
                        double numerator = 0;
                        for(int i = 0; i<MCPUsedForThrust.size(); i++){
                            MCParticle* temp = MCPUsedForThrust.at(i);
                            TLorentzVector TLtemp(temp->getMomentum(), temp->getEnergy());
                            TVector3 TVtemp = TLtemp.Vect();
                            denominator += TVtemp.Mag();
                            numerator += abs(x*TVtemp(0) + y*TVtemp(1) + z*TVtemp(2));
                        }
                        double Ttemp = numerator/denominator;
                        if(Ttemp > T){
                            thetaMin = theta;   phiMin = phi; T = Ttemp;
                            cout<<"*************"<<endl;
                            cout<<"T : "<<T<<"thetaMin : phiMin "<<thetaMin<<" : "<<phiMin<<endl;
                            cout<<"*************"<<endl;
                        }
                    }
                }
                
                
                if(iter == 1){Told = T; Tnew = T;}
                else if(T >= Tnew && iter != 1){
                    Told = Tnew; Tnew = T; cut = (Tnew - Told)/Tnew;
                }
                cout<<"cut : "<<cut<<endl;
            }
            while(cut >= 0.2);

            
            TVector3 tempThrust(0,0,0);
            tempThrust.SetXYZ(sin(thetaMin)*cos(phiMin), sin(thetaMin)*sin(phiMin), cos(thetaMin));
            TVector3 backThrust(0,0,0);
            backThrust = -tempThrust;
            Q1Thrust = TVquark1.Angle(tempThrust);
            Q1BackThrust = TVquark1.Angle(backThrust);
            cout<<"TVquark1.Angle(tempThrust) : "<<TVquark1.Angle(tempThrust)<<" TVquark1.Angle(backThrust) : "<<TVquark1.Angle(backThrust)<<endl;

            //the following code used to get Hemisphere masses
            std::vector<MCParticle* > hemisphere1;
            std::vector<MCParticle* > hemisphere2;
            double visEn = 0;
            double JetBroadeningDenominator = 0;
            for(int i = 0; i<MCPUsedForThrust.size(); i++){
                MCParticle* a_MCP = MCPUsedForThrust.at(i);
                TLorentzVector TLtemp(a_MCP->getMomentum(), a_MCP->getEnergy());
                TVector3 TVtemp = TLtemp.Vect();
                if(TVtemp.Angle(tempThrust) > 0.5*TMath::Pi()){hemisphere1.push_back(a_MCP);}
                else {hemisphere2.push_back(a_MCP);}
                visEn += a_MCP->getEnergy();
                JetBroadeningDenominator += TVtemp.Mag();
            }
            
            cout<<"the number of particles in two hemispheres is "<<hemisphere1.size()+hemisphere2.size()<<endl;
            

            double JetBroadeningNumerator1 = 0, JetBroadeningNumerator2 =0;
            TLorentzVector TLHemi1(0,0,0,0);
            TLorentzVector TLHemi2(0,0,0,0);
            for(int i = 0; i<hemisphere1.size(); i++){
                MCParticle* a_MCP = hemisphere1.at(i);
                TLorentzVector TLtemp(a_MCP->getMomentum(), a_MCP->getEnergy());
                TVector3 TVtemp = TLtemp.Vect();
                TLHemi1 += TLtemp;
                JetBroadeningNumerator1 += abs(TVtemp.Mag() * tempThrust.Mag() * sin(TVtemp.Angle(tempThrust)));
                
            }
            for(int i = 0; i<hemisphere2.size(); i++){
                MCParticle* a_MCP = hemisphere2.at(i);
                TLorentzVector TLtemp(a_MCP->getMomentum(), a_MCP->getEnergy());
                TVector3 TVtemp = TLtemp.Vect();
                TLHemi2 += TLtemp;
                JetBroadeningNumerator2 += abs(TVtemp.Mag() * tempThrust.Mag() * sin(TVtemp.Angle(tempThrust)));
            }
            hemiMass1 = (TLHemi1.M())*(TLHemi1.M())/(visEn*visEn);
            hemiMass2 = (TLHemi2.M())*(TLHemi2.M())/(visEn*visEn);
            hemiBroadening1 = JetBroadeningNumerator1/(2*JetBroadeningDenominator);
            hemiBroadening2 = JetBroadeningNumerator2/(2*JetBroadeningDenominator);
            cout<<"hemiMass1 : "<<hemiMass1<<" hemiMass2 : "<<hemiMass2<<endl;
            cout<<"hemiBroadening1 : "<<hemiBroadening1<<" hemiBroadening2 : "<<hemiBroadening2<<endl;
            
            
            //the following code used to calculate sphericity
            
            double p2Min = 1e-20;
            double denom = 0.;
            double tt[4][4];
            double eVal1 = 0, eVal2 = 0, eVal3 = 0;
            for(int j = 1; j<4; ++j){
                for(int k = j; k<4; ++k){
                    tt[j][k] = 0.;
                }
            }
            for(int i = 0; i<MCPUsedForThrust.size(); i++){
                MCParticle* a_MCP = MCPUsedForThrust.at(i);
                TLorentzVector TLtemp(a_MCP->getMomentum(), a_MCP->getEnergy());
                TVector3 TVtemp = TLtemp.Vect();
                double pNow[4];
                pNow[1] = TLtemp.Px();
                pNow[2] = TLtemp.Py();
                pNow[3] = TLtemp.Pz();
                double p2Now = pow(pNow[1],2) + pow(pNow[2],2) + pow(pNow[3],2);
                double pWeight = 1./sqrt(max(p2Now, p2Min));
                for(int j = 1; j<4; ++j){
                    for(int k=j; k<4; ++k){
                        tt[j][k] += pWeight * pNow[j] * pNow[k];
                        denom += pWeight * p2Now;
                    }
                }
            }
            //Normalize tensor to trace = 1.
            for(int j = 1; j < 4; ++j){
                for(int k = j; k < 4; ++k){
                    tt[j][k] /= denom;
                }
            }
            //find eigenvalues to matrix (third degree equation)
            double qCoef = ( tt[1][1] * tt[2][2] + tt[1][1] * tt[3][3]
                            + tt[2][2] * tt[3][3] - pow(tt[1][2],2) - pow(tt[1][3],2)
                            - pow(tt[2][3],2) ) / 3. - 1./9.;
            double qCoefRt = sqrt( -qCoef);
            double rCoef = -0.5 * ( qCoef + 1./9. + tt[1][1] * pow(tt[2][3],2)
                                   + tt[2][2] * pow(tt[1][3],2) + tt[3][3] * pow(tt[1][2],2)
                                   - tt[1][1] * tt[2][2] * tt[3][3] )
                                   + tt[1][2] * tt[1][3] * tt[2][3] + 1./27.;
            double pTemp = max( min( rCoef / pow(qCoefRt,3), 1.), -1.);
            double pCoef = cos( acos(pTemp) / 3.);
            double pCoefRt = sqrt( 3. * (1. - pow(pCoef,2)) );
            eVal1 = 1./3. + qCoefRt * max( 2. * pCoef,  pCoefRt - pCoef);
            eVal3 = 1./3. + qCoefRt * min( 2. * pCoef, -pCoefRt - pCoef);
            eVal2 = 1. - eVal1 - eVal3;
            CParameter = 3*(eVal1*eVal2 + eVal2*eVal3 + eVal3*eVal1);
            
            cout<<"C parameter is "<<CParameter<<endl;
            
            
            
            
            std::vector<MCParticle* >quark1Par;
            std::vector<MCParticle* >quark2Par;
            
            for(int i = 0; i<FinalStateParticle.size(); i++){
                MCParticle* temp = FinalStateParticle.at(i);
                int NParents = temp->getParents().size();
                int NDaughters = temp->getDaughters().size();
                int PDG = temp->getPDG();
                
                if(NParents != 0){
                    MCParticle* a_parent = temp;
                    do{a_parent = a_parent->getParents()[0];}
                    while(a_parent->getParents().size() != 0);
                    if(a_parent->getParents().size() == 0){
                        TLorentzVector TLtemp(a_parent->getMomentum(), a_parent->getEnergy());
                        if(TLtemp.E() == TLquark1.E()){quark1Par.push_back(temp);}
                        else if(TLtemp.E() == TLquark2.E()){quark2Par.push_back(temp);}
                        }
                    }
                
                else if(NParents == 0){
                    TLorentzVector TLtemp(temp->getMomentum(), temp->getEnergy());
                    if(TLtemp.E() == TLquark1.E()){quark1Par.push_back(temp);}
                    else if(TLtemp.E() == TLquark2.E()){quark2Par.push_back(temp);}
                }
            }
            
            //study the particles belong to each quark
            TLorentzVector TLquark1Par(0,0,0,0);
            TLorentzVector TLquark2Par(0,0,0,0);
            int quark1Neutrino = 0, quark1Neutron = 0, quark1Charge = 0;
            int quark2Neutrino = 0, quark2Neutron = 0, quark2Charge = 0;
            std::vector<MCParticle* > quark1ChargePar;
            std::vector<MCParticle* > quark1NeutronPar;
            std::vector<MCParticle* > quark1VisPar;
            std::vector<MCParticle* > quark1NeutrinoPar;
            std::vector<MCParticle* > quark2ChargePar;
            std::vector<MCParticle* > quark2NeutronPar;
            std::vector<MCParticle* > quark2VisPar;
            std::vector<MCParticle* > quark2NeutrinoPar;
            
            for(int i = 0; i<quark1Par.size(); i++){
                MCParticle* temp = quark1Par.at(i);
                double charge = temp->getCharge();
                int PDG = temp->getPDG();
                if(charge != 0){quark1ChargePar.push_back(temp);}
                if(abs(PDG) != 12 && abs(PDG) != 14 && abs(PDG) != 16){
                    quark1VisPar.push_back(temp);
                    if(charge == 0){quark1NeutronPar.push_back(temp);}
                }
                else if(abs(PDG) == 12 || abs(PDG) == 14 || abs(PDG) == 16){quark1NeutrinoPar.push_back(temp);}
                TLorentzVector TLtemp(temp->getMomentum(), temp->getEnergy());
                TLquark1Par += TLtemp;
            }
            for(int i = 0; i<quark2Par.size(); i++){
                MCParticle* temp = quark2Par.at(i);
                double charge = temp->getCharge();
                int PDG = temp->getPDG();
                if(charge != 0){quark2ChargePar.push_back(temp);}
                if(abs(PDG) != 12 && abs(PDG) != 14 && abs(PDG) != 16){
                    quark2VisPar.push_back(temp);
                    if(charge == 0){quark2NeutronPar.push_back(temp);}
                }
                else if(abs(PDG) == 12 || abs(PDG) == 14 || abs(PDG) == 16){quark2NeutrinoPar.push_back(temp);}
                TLorentzVector TLtemp(temp->getMomentum(), temp->getEnergy());
                TLquark2Par += TLtemp;
            }
            cout<<"the total energy of particles belong to each quark is "<<(TLquark1Par + TLquark2Par).E()<<endl;
            cout<<"the total number of particles belong to each quark is "<<quark1Par.size() + quark2Par.size()<<endl;
            
            
            for(int i = 0; i<MCPVisibleParticleExcISR.size(); i++){
                MCParticle* temp = MCPVisibleParticleExcISR.at(i);
                TLorentzVector TLtemp(temp->getMomentum(), temp->getEnergy());
                TLMCPVisibleParticleExcISR += TLtemp;
            }
            
            En_ISR = TLISR.E();
            cout<<"the energy of final state particles is "<<TLFinalStateParticle.E()<<endl;
            cout<<"the invariant mass of final state particles is "<<TLFinalStateParticle.M()<<endl;
            cout<<"the energy of visible final state particles is "<<TLMCPVisibleParticle.E()<<endl;
            cout<<"the energy of neutrino is "<<TLMCPNeutrino.E()<<endl;
            cout<<"the energy of ISR is "<<TLISR.E()<<endl;
            cout<<"the energy of visible final state particles exclude ISR is "<<TLMCPVisibleParticleExcISR.E()<<endl;
            cout<<"the invariant mass of visible final state particles exclude ISR is "<<TLMCPVisibleParticleExcISR.M()<<endl;
            
            
            map<double, int> map_chargeMCP;
            map<double, int> map_neutronMCP;
            for(int i = 0; i < MCPNeutron.size(); i++){
                MCParticle* temp_particle = MCPNeutron.at(i);
                TLorentzVector TLtemp(temp_particle->getMomentum(), temp_particle->getEnergy());
                double energy = TLtemp.E();
                map_neutronMCP[energy] = i;
            }
            
            
            for(int i = 0; i<MCPCharge.size(); i++){
                MCParticle* temp_particle = MCPCharge.at(i);
                TLorentzVector TLtemp(temp_particle->getMomentum(), temp_particle->getEnergy());
                double energy = TLtemp.E();
                map_chargeMCP[energy] = i;
            }
            
            cout<<"ISR : "<<ISR.size()<<endl;
            cout<<"MCPQuark : "<<MCPQuark.size()<<endl;
            cout<<"FinalStateParticle : "<<FinalStateParticle.size()<<endl;
            cout<<"num_charge : "<<MCPCharge.size()<<" num_neutron : "<<MCPNeutron.size()<<endl;
            cout<<"num_visibleEcxISR : "<<MCPVisibleParticleExcISR.size()<<endl;
            
            //we can develop the algorithm to identify ISR, so map_visibleMCP doesn't include ISR
            map<double, int> map_visibleMCP;
            map_visibleMCP.clear();
            
            for(int i = 0; i < MCPVisibleParticleExcISR.size(); i++){
                MCParticle* temp_particle = MCPVisibleParticleExcISR.at(i);
                TLorentzVector TLtemp(temp_particle->getMomentum(), temp_particle->getEnergy());
                double energy = temp_particle->getEnergy();
                map_visibleMCP[energy] = i;
            }
            
            
            std::vector<MCParticle* > jet1_particle;
            std::vector<MCParticle* > jet2_particle;
            auto it = map_visibleMCP.rbegin(); int num = it->second; double mcp_energy = it->first;
            MCParticle* jet1 = MCPVisibleParticleExcISR.at(num);
            jet1_particle.push_back(jet1);
            TLorentzVector TLjet1(jet1->getMomentum(), jet1->getEnergy());
            TLorentzVector TLjet2(0,0,0,0);
            TVector3 TVjet1(0,0,0);
            TVector3 TVjet2(0,0,0);
            TVjet1 = TLjet1.Vect();
            TVjet2 = -TVjet1;
            
            
            
//            cout<<"it->second : "<<it->second<<" it->first : "<<it->first<<endl;
            cout<<"the number of particles in map_visibleMCP is "<<map_visibleMCP.size()<<endl;
//            cout<<"after sorting according to energy "<<endl;
            for(++it; it != map_visibleMCP.rend(); it ++){
                MCParticle* iter = MCPVisibleParticleExcISR.at(it->second);
                TLorentzVector TLiter(iter->getMomentum(), iter->getEnergy());
                TVector3 TViter = iter->getMomentum();
                double compangle1 = TViter.Angle(TVjet1);
                double compangle2 = TViter.Angle(TVjet2);
                if(compangle1 < compangle2){TLjet1 += TLiter; TVjet1 = TLjet1.Vect(); jet1_particle.push_back(iter);}
                else if(compangle1 > compangle2){TLjet2 += TLiter; TVjet2 = TLjet2.Vect(); jet2_particle.push_back(iter);}
            }
            
            cout<<"TLjet1.E() : "<<TLjet1.E()<<" TLjet2.E() : "<<TLjet2.E()<<" TVjet1.Angle(TVjet2) : "<<TVjet1.Angle(TVjet2)<<endl;
            cout<<"the total energy of two jet is "<<(TLjet1 + TLjet2).E()<<endl;
            cout<<"TLjet1.M() : "<<TLjet1.M()<<" TLjet2.M() : "<<TLjet2.M()<<endl;
            
            

        
        } catch (lcio::DataNotAvailableException err) {  }
	
        _outputTree->Fill();
        Num ++;
	}  	  

}	



void Ztoqq::end()
{

	if (_outputTree) {

		TFile *tree_file = _outputTree->GetCurrentFile(); //just in case we switched to a new file
		//tree_file->cd();
		tree_file->Write();
		delete tree_file;
		//tree_file->Close();
	}

}



