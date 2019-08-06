#include <eventShape.hh>
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

eventShape a_eventShape_instance;

eventShape::eventShape()
	: Processor("eventShape"),
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

void CalcuThrustMCP(std::vector<MCParticle* > UsedForThrust, std::vector<double> &result){
    result.clear();
    double T = 0;
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
//        cout<<"iter : "<<iter<<endl;
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
            
//            cout<<"thetaL : "<<thetaL<<" thetaR : "<<thetaR<<endl;
//            cout<<"phiL : "<<phiL<<" phiR : "<<phiR<<endl;
        }
//        cout<<"thetaRange : "<<thetaRange<<" phiRange : "<<phiRange<<endl;
        for(double theta = thetaL; theta <= thetaR; theta += 0.1*thetaRange){   //in this round, find the max T
            for(double phi = phiL; phi <= phiR; phi += 0.1*phiRange){
                
                double x = sin(theta)*cos(phi);
                double y = sin(theta)*sin(phi);
                double z = cos(theta);
                
                double denominator = 0;
                double numerator = 0;
                for(int i = 0; i<UsedForThrust.size(); i++){
                    MCParticle* temp = UsedForThrust.at(i);
                    TLorentzVector TLtemp(temp->getMomentum(), temp->getEnergy());
                    TVector3 TVtemp = TLtemp.Vect();
                    denominator += TVtemp.Mag();
                    numerator += abs(x*TVtemp(0) + y*TVtemp(1) + z*TVtemp(2));
                }
                double Ttemp = numerator/denominator;
                if(Ttemp > T){
                    thetaMin = theta;   phiMin = phi; T = Ttemp;
 //                   cout<<"*************"<<endl;
 //                   cout<<"T : "<<T<<"thetaMin : phiMin "<<thetaMin<<" : "<<phiMin<<endl;
 //                   cout<<"*************"<<endl;
                }
            }
        }
        if(iter == 1){Told = T; Tnew = T;}
        else if(T >= Tnew && iter != 1){
            Told = Tnew; Tnew = T; cut = (Tnew - Told)/Tnew;
        }
//        cout<<"cut : "<<cut<<endl;
    }
    while(cut >= 0.2);
    
//    result[3] = {T, phiMin, thetaMin};
    result.push_back(T);

    TVector3 tempThrust(0,0,0);
    tempThrust.SetXYZ(sin(thetaMin)*cos(phiMin), sin(thetaMin)*sin(phiMin), cos(thetaMin));

    //the following code used to get Hemisphere masses
    std::vector<MCParticle* > hemisphere1;
    std::vector<MCParticle* > hemisphere2;
    hemisphere1.clear(); hemisphere2.clear();
    double visEn = 0;
    double JetBroadeningDenominator = 0;
    TVector3 TVtotal(0,0,0);
    for(int i = 0; i<UsedForThrust.size(); i++){
        MCParticle* a_MCP = UsedForThrust.at(i);
        TLorentzVector TLtemp(a_MCP->getMomentum(), a_MCP->getEnergy());
        TVector3 TVtemp = TLtemp.Vect();
        TVtotal += TVtemp;
        if(TVtemp.Angle(tempThrust) > 0.5*TMath::Pi()){hemisphere1.push_back(a_MCP);}
        else {hemisphere2.push_back(a_MCP);}
        visEn += a_MCP->getEnergy();
        JetBroadeningDenominator += TVtemp.Mag();
    }
    
//    cout<<"the number of particles in two hemispheres is "<<hemisphere1.size()+hemisphere2.size()<<endl;
    
    
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
    double hemiMass1 = (TLHemi1.M())*(TLHemi1.M())/(visEn*visEn);
    double hemiMass2 = (TLHemi2.M())*(TLHemi2.M())/(visEn*visEn);
    double hemiBroadening1 = JetBroadeningNumerator1/(2*JetBroadeningDenominator);
    double hemiBroadening2 = JetBroadeningNumerator2/(2*JetBroadeningDenominator);
//    cout<<"hemiMass1 : "<<hemiMass1<<" hemiMass2 : "<<hemiMass2<<endl;
//    cout<<"hemiBroadening1 : "<<hemiBroadening1<<" hemiBroadening2 : "<<hemiBroadening2<<endl;
    result.push_back(hemiMass1);
    result.push_back(hemiMass2);
    result.push_back(hemiBroadening1);
    result.push_back(hemiBroadening2);
    
    
    
    //the following code used to get rapidity
    TVector3 TVParaThrust = tempThrust.Orthogonal();
    double transverseM = TVtotal.Perp(TVParaThrust);
    double rapidity = 0.5*log((visEn + transverseM)/(visEn - transverseM));
    result.push_back(rapidity);
    
    cout<<"visEn : "<<visEn<<endl;
    
    
}

//double CalcuMagnitude(std::vector<short> &sign, )

void CalcuThrustReco(std::vector<ReconstructedParticle* > UsedForThrust, std::vector<double> &result)
{
    result.clear();
    double T = 0;
    TVector3 n_T(0,0,0);
    TVector3 n_Ttmp(0,0,0);
    //the following code used to calculate the thrust, the calculation method at https://doi.org/10.1016/0021-9991(83)90010-4
    std::vector<short> s(UsedForThrust.size(),1);      //this vector used to store the sign of P_i
    int i, j, k = 0;

    for(i = 0; i<UsedForThrust.size(); i++)
    {
        ReconstructedParticle* particle_i = UsedForThrust.at(i);
        TVector3 P_i = particle_i->getMomentum();
        for(j=i+1; j<UsedForThrust.size(); j++)
        {
            ReconstructedParticle* particle_j = UsedForThrust.at(j);
            TVector3 P_j = particle_j->getMomentum();
            n_Ttmp.SetXYZ(0,0,0);
            for(k=0; k<UsedForThrust.size(); k++)
            {
                if (k == i || k == j) continue;
                ReconstructedParticle* particle_k = UsedForThrust.at(k);
                TVector3 P_k = particle_k->getMomentum();
                //calculate P_k * (P_i X P_j)
                double stmp = 0;
                stmp = P_k(1)*(P_i(2)*P_j(3) - P_i(3)*P_j(2)) + P_k(2)*(P_i(3)*P_j(1) - P_i(1)*P_j(3)) + P_k(3)*(P_i(1)*P_j(2) - P_i(2)*P_j(1));
                if (stmp >= 0) s[k] = 1;
                if (stmp < 0) s[k] = -1;
                //s[i] = 1; s[j] = -1;
                n_Ttmp += s[k] * P_k;


            }
            //4 combinations of s_i and s_j
            s[i]=1;s[j]=1; n_Ttmp += s[i] * P_i + s[j] * P_j; if (n_Ttmp.Mag()>T) {T = n_Ttmp.Mag();n_T = n_Ttmp.Unit();}
            s[i]=1;s[j]=-1; n_Ttmp += s[i] * P_i + s[j] * P_j; if (n_Ttmp.Mag()>T) {T = n_Ttmp.Mag();n_T = n_Ttmp.Unit();}
            s[i]=-1;s[j]=1; n_Ttmp += s[i] * P_i + s[j] * P_j; if (n_Ttmp.Mag()>T) {T = n_Ttmp.Mag();n_T = n_Ttmp.Unit();}
            s[i]=-1;s[j]=-1; n_Ttmp += s[i] * P_i + s[j] * P_j; if (n_Ttmp.Mag()>T) {T = n_Ttmp.Mag();n_T = n_Ttmp.Unit();}
        }
    }
    double normlization = 0;
    for (i = 0; i<UsedForThrust.size(); i++)
    {
        ReconstructedParticle* particle_a = UsedForThrust.at(i);
        TVector3 P_a = particle_a->getMomentum();
        normlization += P_a.Mag();
    }
    T = T / normlization;




    result.push_back(T);
    
    TVector3 tempThrust(0,0,0);
    tempThrust = n_T;
    
    //the following code used to get Hemisphere masses
    std::vector<ReconstructedParticle* > hemisphere1;
    std::vector<ReconstructedParticle* > hemisphere2;
    hemisphere1.clear(); hemisphere2.clear();
    double visEn = 0;
    double JetBroadeningDenominator = 0;
    TVector3 TVtotal(0,0,0);
    for(int i = 0; i<UsedForThrust.size(); i++){
        ReconstructedParticle* a_MCP = UsedForThrust.at(i);
        TLorentzVector TLtemp(a_MCP->getMomentum(), a_MCP->getEnergy());
        TVector3 TVtemp = TLtemp.Vect();
        TVtotal += TVtemp;
        if(TVtemp.Angle(tempThrust) > 0.5*TMath::Pi()){hemisphere1.push_back(a_MCP);}
        else {hemisphere2.push_back(a_MCP);}
        visEn += a_MCP->getEnergy();
        JetBroadeningDenominator += TVtemp.Mag();
    }
    
//    cout<<"the number of particles in two hemispheres is "<<hemisphere1.size()+hemisphere2.size()<<endl;
    
    
    double JetBroadeningNumerator1 = 0, JetBroadeningNumerator2 =0;
    TLorentzVector TLHemi1(0,0,0,0);
    TLorentzVector TLHemi2(0,0,0,0);
    for(int i = 0; i<hemisphere1.size(); i++){
        ReconstructedParticle* a_MCP = hemisphere1.at(i);
        TLorentzVector TLtemp(a_MCP->getMomentum(), a_MCP->getEnergy());
        TVector3 TVtemp = TLtemp.Vect();
        TLHemi1 += TLtemp;
        JetBroadeningNumerator1 += abs(TVtemp.Mag() * tempThrust.Mag() * sin(TVtemp.Angle(tempThrust)));
        
    }
    for(int i = 0; i<hemisphere2.size(); i++){
        ReconstructedParticle* a_MCP = hemisphere2.at(i);
        TLorentzVector TLtemp(a_MCP->getMomentum(), a_MCP->getEnergy());
        TVector3 TVtemp = TLtemp.Vect();
        TLHemi2 += TLtemp;
        JetBroadeningNumerator2 += abs(TVtemp.Mag() * tempThrust.Mag() * sin(TVtemp.Angle(tempThrust)));
    }
    double hemiMass1 = (TLHemi1.M())*(TLHemi1.M())/(visEn*visEn);
    double hemiMass2 = (TLHemi2.M())*(TLHemi2.M())/(visEn*visEn);
    double hemiBroadening1 = JetBroadeningNumerator1/(2*JetBroadeningDenominator);
    double hemiBroadening2 = JetBroadeningNumerator2/(2*JetBroadeningDenominator);
    //    cout<<"hemiMass1 : "<<hemiMass1<<" hemiMass2 : "<<hemiMass2<<endl;
    //    cout<<"hemiBroadening1 : "<<hemiBroadening1<<" hemiBroadening2 : "<<hemiBroadening2<<endl;
    result.push_back(hemiMass1);
    result.push_back(hemiMass2);
    result.push_back(hemiBroadening1);
    result.push_back(hemiBroadening2);
    
    
    TVector3 TVParaThrust = tempThrust.Orthogonal();
    double transverseM = TVtotal.Perp(TVParaThrust);
    double rapidity = 0.5*log((visEn + transverseM)/(visEn - transverseM));
    result.push_back(rapidity);
//    cout<<"visEn : "<<visEn<<endl;
    
}


double CalcuSphericityMCP(std::vector<MCParticle*> UsedForSphericity, std::vector<double> &result){
    result.clear();
    //the following code used to calculate sphericity
    double CPara = 0;
    double p2Min = 1e-20;
    double denom = 0.;
    double tt[4][4];
    double eVal1 = 0, eVal2 = 0, eVal3 = 0;
    for(int j = 1; j<4; ++j){
        for(int k = j; k<4; ++k){
            tt[j][k] = 0.;
        }
    }
    for(int i = 0; i<UsedForSphericity.size(); i++){
        MCParticle* a_MCP = UsedForSphericity.at(i);
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
    CPara = 3*(eVal1*eVal2 + eVal2*eVal3 + eVal3*eVal1);
    
//    cout<<"C parameter is "<<CPara<<endl;
    result.push_back(CPara);
    double DPara = 27*eVal1*eVal2*eVal3;
    result.push_back(DPara);
//    return CPara;
}

double CalcuSphericityReco(std::vector<ReconstructedParticle*> UsedForSphericity, std::vector<double> &result){
    
    //the following code used to calculate sphericity
    double CPara = 0;
    double p2Min = 1e-20;
    double denom = 0.;
    double tt[4][4];
    double eVal1 = 0, eVal2 = 0, eVal3 = 0;
    for(int j = 1; j<4; ++j){
        for(int k = j; k<4; ++k){
            tt[j][k] = 0.;
        }
    }
    for(int i = 0; i<UsedForSphericity.size(); i++){
        ReconstructedParticle* a_MCP = UsedForSphericity.at(i);
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
    CPara = 3*(eVal1*eVal2 + eVal2*eVal3 + eVal3*eVal1);
    
//    cout<<"C parameter is "<<CPara<<endl;
    result.push_back(CPara);
    double DPara = 27*eVal1*eVal2*eVal3;
    result.push_back(DPara);
    
}






void eventShape::init() {
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
    _outputTree->Branch("Mark", &Mark, "Mark/I");
    _outputTree->Branch("CParameter", &CParameter, "CParameter/D");
    _outputTree->Branch("hemiMass1", &hemiMass1, "hemiMass1/D");
    _outputTree->Branch("hemiMass2", &hemiMass2, "hemiMass2/D");
    _outputTree->Branch("hemiBroadening1", &hemiBroadening1, "hemiBroadening1/D");
    _outputTree->Branch("hemiBroadening2", &hemiBroadening2, "hemiBroadening2/D");
    _outputTree->Branch("En_ISR", &En_ISR, "En_ISR/D");
    _outputTree->Branch("Mass_ISR", &Mass_ISR, "Mass_ISR/D");
    _outputTree->Branch("En_quarks", &En_quarks, "En_quarks/D");
    _outputTree->Branch("Mass_quarks", &Mass_quarks, "Mass_quarks/D");
    _outputTree->Branch("vQuarkPDG", &vQuarkPDG);
    _outputTree->Branch("vHiggsDausQuarks", &vHiggsDausQuarks);
    _outputTree->Branch("vHiggsDaus", &vHiggsDaus);
    _outputTree->Branch("vHiggsDausPDG", &vHiggsDausPDG);
    _outputTree->Branch("num_94", &num_94, "num_94/I");
    _outputTree->Branch("num_quark", &num_quark ,"num_quark/I");
    _outputTree->Branch("angle_twoPlan", &angle_twoPlan, "angle_twoPlan/D");
    _outputTree->Branch("lightQuark", &lightQuark, "lightQuark/I");
    _outputTree->Branch("En_MCPVis", &En_MCPVis, "En_MCPVis/D");
    _outputTree->Branch("MCPEEC", &MCPEEC);
    _outputTree->Branch("ArborEEC", &ArborEEC);
    _outputTree->Branch("gluonAngle", &gluonAngle);
    _outputTree->Branch("gluonEn", &gluonEn);
    _outputTree->Branch("num_gluon", &num_gluon);
    
    _outputTree->Branch("energyNeutronHadron", &energyNeutronHadron, "energyNeutronHadron/D");
    _outputTree->Branch("energyChargeHadron", &energyChargeHadron, "energyChargeHadron/D");
    _outputTree->Branch("energyGamma", &energyGamma, "energyGamma/D");
    _outputTree->Branch("avNeutronHadron", &avNeutronHadron, "avNeutronHadron/D");
    _outputTree->Branch("avChargeHadron", &avChargeHadron, "avChargeHadron/D");
    _outputTree->Branch("avGamma", &avGamma, "avGamma/D");
    _outputTree->Branch("num_chargeHadrom", &num_chargeHadrom, "num_chargeHadrom/I");
    _outputTree->Branch("num_neutronHadron", &num_neutronHadron, "num_neutronHadron/I");
    _outputTree->Branch("num_gamma", &num_gamma, "num_gamma/I");
    _outputTree->Branch("num_chargeLight", &num_chargeLight, "num_chargeLight/I");
    _outputTree->Branch("energyLight", &energyLight, "energyLight/D");
    _outputTree->Branch("avLight", &avLight, "avLight/D");
    
    
    
    _outputTree->Branch("arborThrust", &arborThrust);
    _outputTree->Branch("arborHemiMass1", &arborHemiMass1, "arborHemiMass1/D");
    _outputTree->Branch("arborHemiMass2", &arborHemiMass2, "arborHemiMass2/D");
    _outputTree->Branch("arborHemiBroadening1", &arborHemiBroadening1, "arborHemiBroadening1/D");
    _outputTree->Branch("arborHemiBroadening2", &arborHemiBroadening2, "arborHemiBroadening2/D");
    _outputTree->Branch("arborRapidity", &arborRapidity, "arborRapidity/D");
    _outputTree->Branch("En_arbor", &En_arbor, "En_arbor/D");
    _outputTree->Branch("Pt_arbor", &Pt_arbor, "Pt_arbor/D");
    _outputTree->Branch("Mag_arbor", &Mag_arbor, "Mag_arbor/D");
    _outputTree->Branch("Mass_arbor", &Mass_arbor, "Mass_arbor/D");
    
    _outputTree->Branch("energy_aLight", &energy_aLight, "energy_aLight/D");
    _outputTree->Branch("energy_aGama", &energy_aGama, "energy_aGama/D");
    _outputTree->Branch("energy_aNeuHad", &energy_aNeuHad, "energy_aNeuHad/D");
    _outputTree->Branch("energy_aChgHad", &energy_aChgHad, "energy_aChgHad/D");
    _outputTree->Branch("num_aLight", &num_aLight, "num_aLight/I");
    _outputTree->Branch("num_aGama", &num_aGama, "num_aGama/I");
    _outputTree->Branch("num_aNeuHad", &num_aNeuHad, "num_aNeuHad/I");
    _outputTree->Branch("num_aChgHad", &num_aChgHad, "num_aChgHad/I");
    _outputTree->Branch("av_aLight", &av_aLight, "av_aLight/D");
    _outputTree->Branch("av_aGama", &av_aGama, "av_aGama/D");
    _outputTree->Branch("av_aNeuHad", &av_aNeuHad, "av_aNeuHad/D");
    _outputTree->Branch("av_aChgHad", &av_aChgHad, "av_aChgHad/D");

    

    _outputTree->Branch("MCPVisISRThrust", &MCPVisISRThrust, "MCPVisISRThrust/D");
    _outputTree->Branch("MCPVisISRHemiMass1", &MCPVisISRHemiMass1, "MCPVisISRHemiMass1/D");
    _outputTree->Branch("MCPVisISRHemiMass2", &MCPVisISRHemiMass2, "MCPVisISRHemiMass2/D");
    _outputTree->Branch("MCPVisISRHemiBroadening1", &MCPVisISRHemiBroadening1, "MCPVisISRHemiBroadening1/D");
    _outputTree->Branch("MCPVisISRHemiBroadening2", &MCPVisISRHemiBroadening2, "MCPVisISRHemiBroadening2/D");
    _outputTree->Branch("MCPVisISRRapidity", &MCPVisISRRapidity, "MCPVisISRRapidity/D");
    
    _outputTree->Branch("MCPVisExcISRThrust", &MCPVisExcISRThrust, "MCPVisExcISRThrust/D");
    _outputTree->Branch("MCPVisExcISRHemiMass1", &MCPVisExcISRHemiMass1, "MCPVisExcISRHemiMass1/D");
    _outputTree->Branch("MCPVisExcISRHemiMass2", &MCPVisExcISRHemiMass2, "MCPVisExcISRHemiMass2/D");
    _outputTree->Branch("MCPVisExcISRHemiBroadening1", &MCPVisExcISRHemiBroadening1, "MCPVisExcISRHemiBroadening1/D");
    _outputTree->Branch("MCPVisExcISRHemiBroadening2", &MCPVisExcISRHemiBroadening2, "MCPVisExcISRHemiBroadening2/D");

    _outputTree->Branch("MCPExcISRThrust", &MCPExcISRThrust, "MCPExcISRThrust/D");
    _outputTree->Branch("MCPExcISRHemiMass1", &MCPExcISRHemiMass1, "MCPExcISRHemiMass1/D");
    _outputTree->Branch("MCPExcISRHemiMass2", &MCPExcISRHemiMass2, "MCPExcISRHemiMass2/D");
    _outputTree->Branch("MCPExcISRHemiBroadening1", &MCPExcISRHemiBroadening1, "MCPExcISRHemiBroadening1/D");
    _outputTree->Branch("MCPExcISRHemiBroadening2", &MCPExcISRHemiBroadening2, "MCPExcISRHemiBroadening2/D");
    
    _outputTree->Branch("MCPVisISRSphericity", &MCPVisISRSphericity, "MCPVisISRSphericity/D");
    _outputTree->Branch("MCPVisISRDPara", &MCPVisISRDPara, "MCPVisISRDPara/D");
    _outputTree->Branch("MCPVisExcIsrSphericity", &MCPVisExcIsrSphericity, "MCPVisExcIsrSphericity/D");
    _outputTree->Branch("MCPExcISRSphericity", &MCPExcISRSphericity, "MCPExcISRSphericity/D");
    _outputTree->Branch("arborSphericity", &arborSphericity, "arborSphericity/D");
    _outputTree->Branch("arborDPara", &arborDPara, "arborDPara/D");
	Num = 0;
}

void eventShape::processEvent( LCEvent * evtP )
{		

	if (evtP) 								
	{
        try{
            vQuarkPDG.clear();
            eventNr = evtP->getEventNumber();
            cout<<"**************************************"<<endl;
            cout<<"eventNr "<<eventNr<<" Num "<<Num<<endl;
            Mark = 3;  //ZH->4quarks: 1    ZZ->4quarks: 2     WW->4quarks: 3    ZH->WW*->6quarks: 4  ZH->ZZ*->6quarks: 5
            
            
            
            //for ArborPFO
            LCCollection* col_PFO = evtP->getCollection( "ArborPFOs" );
            int nPFO = col_PFO->getNumberOfElements();
            TLorentzVector TLArborPFO(0,0,0,0);
            TLorentzVector TLArborCharge(0,0,0,0);
            std::vector<ReconstructedParticle* > arborPFO;
            std::vector<ReconstructedParticle* > arborHadron;
            std::vector<ReconstructedParticle* > arborCharge;
            std::vector<ReconstructedParticle* > arborGamma;
            arborHadron.clear();
            arborCharge.clear();
            arborGamma.clear();
            for(int i = 0; i<nPFO; i++)
            {
                ReconstructedParticle* a_Reco = dynamic_cast<ReconstructedParticle*>(col_PFO->getElementAt(i));
                arborPFO.push_back(a_Reco);
                TLorentzVector temp(a_Reco->getMomentum(), a_Reco->getEnergy());
                TVector3 TVtemp = temp.Vect();
                TLArborPFO += temp;
                float Rcostheta = temp.CosTheta();
                int type = a_Reco->getType();
                int Rcharge = a_Reco->getCharge();
                if(Rcharge != 0){TLArborCharge += temp;}
                if(type == 22){arborGamma.push_back(a_Reco);}
                if(type == 11 || type == 13 || type == 15){arborCharge.push_back(a_Reco);}
                if(type != 11 && type != 13 && type != 15 && type != 22){arborHadron.push_back(a_Reco);}
            }
            En_arbor = TLArborPFO.E();
            Mass_arbor = TLArborPFO.M();
            TVector3 TVArborPFO = TLArborPFO.Vect();
            Pt_arbor = TVArborPFO.Perp();
            Mag_arbor = TVArborPFO.Mag();
            double En_ArborCharge = TLArborCharge.E();
            cout<<"En_ArborCharge : "<<En_ArborCharge<<endl;
            TVector3 TVArborCharge = TLArborCharge.Vect();
            cout<<"TVArborCharge.Perp() : "<<TVArborCharge.Perp()<<endl;
            
            
            ArborEEC.clear();
            
            double arborEEC[100] = {0};
            for(int i =0; i<nPFO; i++){
                for(int j = 0; j<nPFO; j++){
                    ReconstructedParticle* a_Reco = dynamic_cast<ReconstructedParticle*>(col_PFO->getElementAt(i));
                    ReconstructedParticle* b_Reco = dynamic_cast<ReconstructedParticle*>(col_PFO->getElementAt(j));
                    TLorentzVector TLaReco(a_Reco->getMomentum(), a_Reco->getEnergy());
                    TLorentzVector TLbReco(b_Reco->getMomentum(), b_Reco->getEnergy());
                    TVector3 TVaReco = TLaReco.Vect();
                    TVector3 TVbReco = TLbReco.Vect();
                    double angle_ab = TVaReco.Angle(TVbReco);
                    int cos_ab = cos(angle_ab)*50;
                    arborEEC[cos_ab + 49] += (a_Reco->getEnergy() * b_Reco->getEnergy());
                }
            }
            
            double arborsum = 0;
            for(int i = 0; i<100; i++){
                ArborEEC.push_back( arborEEC[i] / pow(TLArborPFO.E(), 2) );
                arborsum += arborEEC[i]/pow(TLArborPFO.E(), 2);
            }
            cout<<"arborsum : "<<arborsum<<endl;
            
            
            
            
            num_aNeuHad = 0;
            num_aChgHad = 0;
            num_aGama = arborGamma.size();
            num_aLight = arborCharge.size();
            energy_aNeuHad = 0;
            energy_aChgHad = 0;
            av_aNeuHad = 0;
            av_aChgHad = 0;
            energy_aGama = 0;
            av_aGama = 0;
            energy_aLight = 0;
            av_aLight = 0;
            
            for(int i = 0; i<arborHadron.size(); i++){
                ReconstructedParticle* a_Reco = arborHadron.at(i);
                double charge = a_Reco->getCharge();
                double energy = a_Reco->getEnergy();
                if(charge == 0){
                    num_aNeuHad += 1;
                    av_aNeuHad += energy;
                    if(energy > energy_aNeuHad){energy_aNeuHad = energy;}
                }
                if(charge != 0){
                    num_aChgHad += 1;
                    av_aChgHad += energy;
                    if(energy > energy_aChgHad){energy_aChgHad = energy;}
                }
            }
            cout<<"av_aNeuHad : "<<av_aNeuHad<<" av_aChgHad : "<<av_aChgHad<<endl;
            if(num_aNeuHad != 0){ av_aNeuHad = av_aNeuHad/num_aNeuHad;}
            if(num_aChgHad != 0) {av_aChgHad = av_aChgHad/num_aChgHad;}
            cout<<"num_aNeuHad : "<<num_aNeuHad<<" num_aChgHad : "<<num_aChgHad<<endl;
            cout<<"av_aNeuHad : "<<av_aNeuHad<<" av_aChgHad : "<<av_aChgHad<<endl;

            for(int i = 0; i<arborGamma.size(); i++){
                ReconstructedParticle* a_Reco = arborGamma.at(i);
                double charge = a_Reco->getCharge();
                double energy = a_Reco->getEnergy();
                av_aGama += energy;
                if(energy > energy_aGama){energy_aGama = energy;}
            }
            cout<<"energy_aGama : "<<energy_aGama<<" energy_aNeuHad : "<<energy_aNeuHad<<" energy_aChgHad : "<<energy_aChgHad<<endl;
            cout<<"av_aGama"<<av_aGama<<endl;
            if(num_aGama != 0){ av_aGama = av_aGama/num_aGama; }
            cout<<"av_aGama : "<<av_aGama<<" num_aGama : "<<num_aGama<<endl;
            cout<<""<<num_aLight<<endl;
            
            for(int i = 0; i<arborCharge.size(); i++){
                ReconstructedParticle* a_Reco = arborCharge.at(i);
                double energy = a_Reco->getEnergy();
                av_aLight += energy;
                if(energy> energy_aLight){energy_aLight = energy;}
            }
            if(num_aLight != 0){av_aLight = av_aLight/num_aLight;}
            cout<<"av_aLight : "<<av_aLight<<" energy_aLight : "<<energy_aLight<<endl;
            
            cout<<"the number of arborPFOs is "<<arborHadron.size() + arborCharge.size() + arborGamma.size()<<endl;
            cout<<"nPFO : "<<nPFO<<endl;
            
            
//            double arbor[3] = {0,0,0};
            std::vector<double> arbor; arbor.clear();
            CalcuThrustReco(arborPFO, arbor);
            arborThrust = arbor.at(0);
            arborHemiMass1 = arbor.at(1);
            arborHemiMass2 = arbor.at(2);
            arborHemiBroadening1 = arbor.at(3);
            arborHemiBroadening2 = arbor.at(4);
            arborRapidity = arbor.at(5);
//            cout<<"arborThrust : "<<arborThrust<<endl;
            std::vector<double> arborS; arborS.clear();
            CalcuSphericityReco(arborPFO, arborS);
            arborSphericity = arborS.at(0);
            arborDPara = arborS.at(1);
            cout<<"arborDPara : "<<arborDPara<<endl;
            cout<<"arborRapidity : "<<arborRapidity<<endl;
            
            
            //for MCParticle
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
            std::vector<MCParticle* >MCPBoson1;
            std::vector<MCParticle* >MCPBoson2;
            std::vector<MCParticle* >MCPHadron;
            std::vector<MCParticle* >MCPChargeLight;
            
            
            
            FinalStateParticle.clear();
            MCPCharge.clear();
            MCPNeutron.clear();
            MCPUsedForThrust.clear();
            MCPVisibleParticleExcISR.clear();
            MCPVisibleParticle.clear();
            MCPQuark.clear();
            ISR.clear();
            MCPGamma.clear();
            MCPGammaExcISR.clear();
            MCPHadron.clear();
            MCPChargeLight.clear();
            
            
            TLorentzVector TLFinalStateParticle(0,0,0,0);
            TLorentzVector TLMCPGamma(0,0,0,0);
            TLorentzVector TLMCPGammaExcISR(0,0,0,0);
            TLorentzVector TLISR(0,0,0,0);
            TLorentzVector TLMCPNeutrino(0,0,0,0);
            TLorentzVector TLMCPCharge(0,0,0,0);
            TLorentzVector TLMCPNeutron(0,0,0,0);
            TLorentzVector TLMCPVisibleParticle(0,0,0,0);
            TLorentzVector TLMCPVisibleParticleExcISR(0,0,0,0);
    
            
            int num_neutron = 0, num_neutrino = 0;
            num_94 = 0;
            MCParticle* first_94;
            MCParticle* second_94;
            MCParticle* third_94;
            vHiggsDausPDG.clear();
            vHiggsDaus.clear();
            vHiggsDausQuarks.clear();
            
            for(int i = 0; i < n_MCP; i++){
                MCParticle* a_MCP = dynamic_cast<MCParticle*>(col_MCP->getElementAt(i));
                int NParents = a_MCP->getParents().size();
                int NDaughters = a_MCP->getDaughters().size();
                int PDG = a_MCP->getPDG();
                double charge = a_MCP->getCharge();
                TLorentzVector TLMCPtemp(a_MCP->getMomentum(), a_MCP->getEnergy());
                double costheta = TLMCPtemp.CosTheta();
                TVector3 TVMCPtemp = TLMCPtemp.Vect();
                
                
                
                if(abs(PDG) == 25 && NDaughters == 2){
                    MCParticle* HiggsDau1 = a_MCP->getDaughters()[0];
                    MCParticle* HiggsDau2 = a_MCP->getDaughters()[1];
                    
                    cout<<"PDG of Higgs Daughters is "<<HiggsDau1->getPDG()<<" : "<<HiggsDau2->getPDG()<<endl;
                    
                    if(HiggsDau1->getPDG() == 23 || HiggsDau1->getPDG() == 24)
                    {
                        MCParticle* Boson1Dau1 = HiggsDau1->getDaughters()[0];
                        MCParticle* Boson1Dau2 = HiggsDau1->getDaughters()[1];
                        
                        MCParticle* Boson2Dau1 = HiggsDau2->getDaughters()[0];
                        MCParticle* Boson2Dau2 = HiggsDau2->getDaughters()[1];
                        
                        cout<<"the energy of Higgs Daughter 1 is "<<HiggsDau1->getEnergy()<<", daughters' PDG is "<<Boson1Dau1->getPDG()<<" : "<<Boson1Dau2->getPDG()<<endl;
                        cout<<"the energy of Higgs Daughter 2 is "<<HiggsDau2->getEnergy()<<", daughters' PDG is "<<Boson2Dau1->getPDG()<<" : "<<Boson2Dau2->getPDG()<<endl;
                        
                        vHiggsDausPDG.push_back(HiggsDau1->getPDG());
                        vHiggsDausPDG.push_back(HiggsDau2->getPDG());
                        vHiggsDausPDG.push_back(Boson1Dau1->getPDG());
                        vHiggsDausPDG.push_back(Boson1Dau2->getPDG());
                        vHiggsDausPDG.push_back(Boson2Dau1->getPDG());
                        vHiggsDausPDG.push_back(Boson2Dau2->getPDG());
                        
                        vHiggsDaus.push_back(HiggsDau1->getEnergy());
                        vHiggsDaus.push_back(HiggsDau2->getEnergy());
                        
                        MCPBoson1.push_back(HiggsDau1);
                        MCPBoson1.push_back(Boson1Dau1);
                        MCPBoson1.push_back(Boson1Dau2);
                        
                        MCPBoson2.push_back(HiggsDau2);
                        MCPBoson2.push_back(Boson2Dau1);
                        MCPBoson2.push_back(Boson2Dau2);
                    }
                    
                    if(HiggsDau1->getPDG() != 23 && HiggsDau1->getPDG() != 24){
                        vHiggsDausQuarks.push_back(HiggsDau1->getPDG());
                        vHiggsDausQuarks.push_back(HiggsDau2->getPDG());
                    }

                    
                }
                
                if(PDG == 94){
                    num_94 += 1;
                    if(num_94 == 1){first_94 = a_MCP;}
                    else if(num_94 == 2){second_94 = a_MCP;}
                    else {third_94 = a_MCP;}
                }
                
                if(PDG == 21 && NParents == 1){
                    MCParticle* parent = a_MCP->getParents()[0];
                    int parentPDG = parent->getPDG();
                    if(abs(parentPDG) == 1 || abs(parentPDG) == 2 || abs(parentPDG) == 3 || abs(parentPDG) == 4 || abs(parentPDG) == 5 || abs(parentPDG) == 6){
                        
                    }
                }
                
                
                if(NParents == 0){
                    //the original quarks
                    if(abs(PDG) == 1 || abs(PDG) == 2 || abs(PDG) == 3 || abs(PDG) == 4 || abs(PDG) == 5 || abs(PDG) == 6){
                        MCPQuark.push_back(a_MCP);
                    }
                    //ISR
                    if(PDG == 22){
                        TLISR += TLMCPtemp;
                        ISR.push_back(a_MCP);
                    }
                }
                
                if(a_MCP->getGeneratorStatus() == 1){
                    TLFinalStateParticle += TLMCPtemp;
                    //final state particles include ISR and neutrinos， mainly used for checking
                    FinalStateParticle.push_back(a_MCP);
                    
                    if(abs(PDG) != 11 && abs(PDG) != 12 && abs(PDG) != 13 && abs(PDG) != 14 && abs(PDG) != 15 && abs(PDG) != 16 && abs(PDG) != 22){
                        MCPHadron.push_back(a_MCP);
                    }
                    
                    if(abs(PDG) == 22){MCPGamma.push_back(a_MCP);}
                    if(abs(PDG) == 11 || abs(PDG) == 13 || abs(PDG) == 15){MCPChargeLight.push_back(a_MCP);}
                    
                    
                    
                    if(PDG != 22){MCPUsedForThrust.push_back(a_MCP);} //used for store the final state particles doesn't decay from ISR
                    else if(PDG == 22){
                        MCParticle* b_parent = a_MCP;
                        do{b_parent = b_parent->getParents()[0];}
                        while(b_parent->getParents().size() != 0);
                        if(b_parent->getParents().size() == 0 && (b_parent->getPDG()) != 22){
                            MCPUsedForThrust.push_back(a_MCP);
                        }
                    }
                    
                    //store neutron particles include gamma exclude neutrino
                    if(charge == 0 && abs(PDG) != 12 && abs(PDG) != 14 && abs(PDG) != 16){MCPNeutron.push_back(a_MCP);}
                    
                    if(charge != 0){MCPCharge.push_back(a_MCP);}
                    
                    if(abs(PDG) == 12 || abs(PDG) == 14 || abs(PDG) == 16){
                        TLMCPNeutrino += TLMCPtemp;
                        MCPNeutrino.push_back(a_MCP);
                    }
                    //store all visible final state particles include ISR
                    else {MCPVisibleParticle.push_back(a_MCP);  TLMCPVisibleParticle += TLMCPtemp;}
                    
                    if(PDG != 22 && abs(PDG) != 12 && abs(PDG) != 14 && abs(PDG) != 16){MCPVisibleParticleExcISR.push_back(a_MCP);}
                    MCParticle* a_parent = a_MCP;
                    do{a_parent = a_parent->getParents()[0];}
                    while(a_parent->getParents().size() != 0);
                    if(PDG == 22 && a_parent->getPDG() != 22 && a_parent->getParents().size() == 0){ // this photon is not a ISR
                        MCPVisibleParticleExcISR.push_back(a_MCP);
                    }
                }
            }
            
            cout<<"num_94 : "<<num_94<<endl;
//            cout<<"FinalStateParticle.size() : "<<FinalStateParticle.size()<<" MCPUsedForThrust.size() : "<<MCPUsedForThrust.size()<<endl;
            
            
            gluonEn.clear();
            gluonAngle.clear();
            num_gluon = 0;
            for(int i = 0; i < n_MCP; i++){
                MCParticle* a_MCP = dynamic_cast<MCParticle*>(col_MCP->getElementAt(i));
                int NParents = a_MCP->getParents().size();
                int NDaughters = a_MCP->getDaughters().size();
                int PDG = a_MCP->getPDG();
                TVector3 TVparticle = a_MCP->getMomentum();
                MCParticle* quark1 = MCPQuark.at(0);
                MCParticle* quark2 = MCPQuark.at(1);
                TVector3 TVquark1 = quark1->getMomentum();
                TVector3 TVquark2 = quark2->getMomentum();
                
                if(PDG == 21 && NParents == 1){
                    MCParticle* parent = a_MCP->getParents()[0];
                    int parentPDG = parent->getPDG();
                    if(abs(parentPDG) == 1 || abs(parentPDG) == 2 || abs(parentPDG) == 3 || abs(parentPDG) == 4 || abs(parentPDG) == 5 || abs(parentPDG) == 6){
                        num_gluon += 1;
                        double angleProduct = TVparticle.Angle(TVquark1) * TVparticle.Angle(TVquark2);
                        double energygluon = a_MCP->getEnergy();
                        gluonAngle.push_back(angleProduct);
                        gluonEn.push_back(energygluon);
                        cout<<"energygluon : "<<energygluon<<" angleProduct : "<<angleProduct<<endl;
                    }
                }
                
            }
            cout<<"num_gluon : "<<num_gluon<<endl;

            
            
            //study ZH, Higgs to WW or ZZ
            int num_charge = MCPCharge.size();
            int first_Charge_PDG = 0, second_Charge_PDG = 0;
            double first_Charge_energy = 0, second_Charge_energy = 0;
            for(int i = 0; i<num_charge; i++){
                MCParticle* chargePar = MCPCharge.at(i);
                int charge_PDG = chargePar->getPDG();
                double charge_energy = chargePar->getEnergy();
                if(charge_energy > first_Charge_energy){
                    second_Charge_energy = first_Charge_energy;
                    second_Charge_PDG = first_Charge_PDG;
                    first_Charge_energy = charge_energy;
                    first_Charge_PDG = charge_PDG;
                }
            }
            
            cout<<"first_Charge_energy : "<<first_Charge_energy<<" second_Charge_energy : "<<second_Charge_energy<<endl;
            cout<<"first_Charge_PDG : "<<first_Charge_PDG<<" second_Charge_PDG : "<<second_Charge_PDG<<endl;
            
            
            
            //study hadron
            num_chargeHadrom = 0;
            num_neutronHadron = 0;
            num_gamma = MCPGamma.size();
            int num_hadron = MCPHadron.size();
            energyNeutronHadron = 0;
            energyChargeHadron = 0;
            avNeutronHadron = 0;
            avChargeHadron = 0;
            energyGamma = 0;
            avGamma = 0;
            
            for(int i = 0; i<num_hadron; i++){
                MCParticle* hpar = MCPHadron.at(i);
                double hcharge = hpar->getCharge();
                double henergy = hpar->getEnergy();
                if(hcharge == 0){
                    num_neutronHadron += 1;
                    avNeutronHadron += henergy;
                    if(henergy > energyNeutronHadron){energyNeutronHadron = henergy; cout<<"energyNeutronHadron : "<<energyNeutronHadron<<endl;}
                }
                else if(hcharge != 0){
                    num_chargeHadrom += 1;
                    avChargeHadron += henergy;
                    if(henergy > energyChargeHadron){energyChargeHadron = henergy; cout<<"energyChargeHadron : "<<energyChargeHadron<<endl;}
                }
            }
            cout<<"the energy of hadron is : "<<avChargeHadron+avNeutronHadron<<endl;
            if(num_chargeHadrom != 0){avChargeHadron = avChargeHadron/num_chargeHadrom;}
            if(num_neutronHadron != 0){avNeutronHadron = avNeutronHadron/num_neutronHadron;}
            else avNeutronHadron = 0;
            cout<<"the average energy of hadron for two conditions :"<<avChargeHadron<<" : "<<avNeutronHadron<<endl;
            
            for(int i = 0; i<num_gamma; i++){
                MCParticle* Gpar = MCPGamma.at(i);
                double Genergy = Gpar->getEnergy();
                avGamma += Genergy;
                if(Genergy > energyGamma){energyGamma = Genergy; cout<<"energyGamma : "<<energyGamma<<endl;}
            }
            cout<<"avGamma"<<avGamma<<endl;
            if(num_gamma != 0){ avGamma = avGamma/num_gamma;}
            cout<<"avGamma"<<avGamma<<endl;

            num_chargeLight = MCPChargeLight.size();
            energyLight = 0;
            avLight = 0;
            for(int i = 0; i<num_chargeLight; i++){
                MCParticle* Lpar = MCPChargeLight.at(i);
                double Lenergy = Lpar->getEnergy();
                avLight += Lenergy;
                if(Lenergy > energyLight){energyLight = Lenergy; cout<<"energyLight : "<<energyLight<<endl;;}
            }
            cout<<"avLight : "<<avLight<<endl;
            if(num_chargeLight != 0){avLight = avLight/num_chargeLight;}
            else {avLight = 0;}
            cout<<"avLight : "<<avLight<<endl;

            
            // stydy original quarks
            TLorentzVector TLMCPQuark(0,0,0,0);
            num_quark = MCPQuark.size();
            lightQuark = 0;
            for(int i = 0; i<num_quark; i++){
                MCParticle* quark = MCPQuark.at(i);
                vQuarkPDG.push_back(quark->getPDG());
                TLorentzVector TLquarktemp(quark->getMomentum(), quark->getEnergy());
                TLMCPQuark += TLquarktemp;
                
                if(abs(quark->getPDG()) == 1 || abs(quark->getPDG()) == 2 || abs(quark->getPDG()) == 3){lightQuark += 1;}
            }
            cout<<"lightQuark : "<<lightQuark<<endl;
            En_quarks = TLMCPQuark.E();
            Mass_quarks = TLMCPQuark.M();
            //the following code used to calculate the angle between two boson plan
            std::vector<MCParticle* >quark_Group1;
            std::vector<MCParticle* >quark_Group2;
            if(num_94 == 2 && num_quark == 4){
                MCParticle* first1 = first_94->getParents()[0];
                MCParticle* first2 = first_94->getParents()[1];
                MCParticle* second1 = second_94->getParents()[0];
                MCParticle* second2 = second_94->getParents()[1];
                if(first1->getParents().size() != 0){
                    do{first1 = first1->getParents()[0];}
                    while(first1->getParents().size() != 0);
                    if(first1->getParents().size() == 0){quark_Group1.push_back(first1);}
                    }
                else {quark_Group1.push_back(first1);}
                if(first2->getParents().size() != 0){
                    do{first2 = first2->getParents()[0];}
                    while(first2->getParents().size() != 0);
                    if(first2->getParents().size() == 0){quark_Group1.push_back(first2);}
                }
                else {quark_Group1.push_back(first2);}
                
                if(second1->getParents().size() != 0){
                    do{second1 = second1->getParents()[0];}
                    while(second1->getParents().size() != 0);
                    if(second1->getParents().size() == 0){quark_Group2.push_back(second1);}
                }
                else {quark_Group2.push_back(second1);}
                if(second2->getParents().size() != 0){
                    do{second2 = second2->getParents()[0];}
                    while(second2->getParents().size() != 0);
                    if(second2->getParents().size() == 0){quark_Group2.push_back(second2);}
                }
                else {quark_Group2.push_back(second2);}
                
            }
            cout<<"Group1.size() : Group2.size() "<<quark_Group1.size()<<" : "<<quark_Group2.size()<<endl;
            
            angle_twoPlan = 999;
            
            if(num_94 == 2 && quark_Group1.size() == 2 && quark_Group2.size() == 2){
                MCParticle* first1B = quark_Group1.at(0);
                MCParticle* first2B = quark_Group1.at(1);
                int PDGf1 = first1B->getPDG();
                int PDGf2 = first2B->getPDG();
                MCParticle* first1;
                MCParticle* first2;
                first1 = (PDGf1 > PDGf2) ? first1B : first2B;
                first2 = (PDGf1 > PDGf2) ? first2B : first1B;
                TLorentzVector TLfirst1(first1->getMomentum(), first1->getEnergy());
                TLorentzVector TLfirst2(first2->getMomentum(), first2->getEnergy());
                TLorentzVector TL941(0,0,0,0);
                TL941 = TLfirst1 + TLfirst2;
                TVector3 TV941 = TL941.Vect();
                TVector3 TVfirst1 = first1->getMomentum();
                TVector3 TVfirst2 = first2->getMomentum();
                double firstPx = TVfirst1.Py()*TVfirst2.Pz() - TVfirst1.Pz()*TVfirst2.Py();
                double firstPy = TVfirst1.Pz()*TVfirst2.Px() - TVfirst2.Pz()*TVfirst1.Px();
                double firstPz = TVfirst1.Px()*TVfirst2.Py() - TVfirst1.Py()*TVfirst2.Px();
                TVector3 TVper1(firstPx, firstPy, firstPz);
                
                MCParticle* second1B = quark_Group2.at(0);
                MCParticle* second2B = quark_Group2.at(1);
                int PDGs1 = second1B->getPDG();
                int PDGs2 = second2B->getPDG();
                MCParticle* second1;
                MCParticle* second2;
                second1 = (PDGs1 > PDGs2) ? second1B : second2B;
                second2 = (PDGs1 > PDGs2) ? second2B : second1B;
                TLorentzVector TLsecond1(second1->getMomentum(), second1->getEnergy());
                TLorentzVector TLsecond2(second2->getMomentum(), second2->getEnergy());
                TLorentzVector TL942(0,0,0,0);
                TL942 = TLsecond1 + TLsecond2;
                TVector3 TV942 = TL942.Vect();
                TVector3 TVsecond1 = second1->getMomentum();
                TVector3 TVsecond2 = second2->getMomentum();
                double secondPx = TVsecond1.Py()*TVsecond2.Pz() - TVsecond1.Pz()*TVsecond2.Py();
                double secondPy = TVsecond1.Pz()*TVsecond2.Px() - TVsecond2.Pz()*TVsecond1.Px();
                double secondPz = TVsecond1.Px()*TVsecond2.Py() - TVsecond1.Py()*TVsecond2.Px();
                TVector3 TVper2(secondPx, secondPy, secondPz);
                
//                double angle1 = TVper1.Angle(TVper2);
//                double angle2 = TMath::Pi() - angle1;
//                double angleMax = (angle1 >= angle2) ? angle1 : angle2;
//                double angleMin = (angle1 >= angle2) ? angle2 : angle1;
//                if(TV941.Angle(TV942) > 0.5*TMath::Pi()){angle_twoPlan = angleMax;}
//                else {angle_twoPlan = angleMin;}
//                cout<<"angleMax : "<<angleMax<<" angleMin : "<<angleMin<<" TV941.Angle(TV942) : "<<TV941.Angle(TV942)<<endl;
                angle_twoPlan = TVper1.Angle(TVper2);
            }
            cout<<"angle_twoPlan : "<<angle_twoPlan<<endl;
            

            
            

            std::vector<double> MCPVisISR;
            std::vector<double> MCPVisExcISR;
            std::vector<double> MCPExcISR;
            MCPVisExcISR.clear(); MCPVisISR.clear(); MCPExcISR.clear();
            

            CalcuThrustMCP(MCPVisibleParticle, MCPVisISR);
            MCPVisISRThrust          = MCPVisISR.at(0);
            MCPVisISRHemiMass1       = MCPVisISR.at(1);
            MCPVisISRHemiMass2       = MCPVisISR.at(2);
            MCPVisISRHemiBroadening1 = MCPVisISR.at(3);
            MCPVisISRHemiBroadening2 = MCPVisISR.at(4);
            MCPVisISRRapidity        = MCPVisISR.at(5);
            std::vector<double> MCPVisISRS;
            CalcuSphericityMCP(MCPVisibleParticle, MCPVisISRS);
            MCPVisISRSphericity = MCPVisISRS.at(0);
            MCPVisISRDPara = MCPVisISRS.at(1);

            cout<<"MCPVisISRRapidity : "<<MCPVisISRRapidity<<endl;
            cout<<"MCPVisISRDPara : "<<MCPVisISRDPara<<endl;
            
            
            
            
            
            
            
            
            
//            std::vector<double> MCPEEC; MCPVisibleParticle
            MCPEEC.clear();
            
            double EEC[100] = {0};
            for(int i = 0; i<MCPVisibleParticle.size(); i++){
                for(int j = 0; j<MCPVisibleParticle.size(); j++){
                    MCParticle* MCP_i = MCPVisibleParticle.at(i);
                    MCParticle* MCP_j = MCPVisibleParticle.at(j);
                    TLorentzVector TLi(MCP_i->getMomentum(), MCP_i->getEnergy());
                    TLorentzVector TLj(MCP_j->getMomentum(), MCP_j->getEnergy());
                    TVector3 TVi = TLi.Vect();
                    TVector3 TVj = TLj.Vect();
                    double angle_ij = TVi.Angle(TVj);
                    int cos_ij = cos(angle_ij) * 50;
//                    cout<<"cos_ij : "<<cos_ij<<endl;
                    EEC[cos_ij + 49] += (MCP_i->getEnergy() * MCP_j->getEnergy());

                }
            }
            
            
            double sum = 0;
            for(int i=0; i<100; i++){
                MCPEEC.push_back( EEC[i] / pow(TLMCPVisibleParticle.E(), 2) );
                sum += EEC[i]/pow(TLMCPVisibleParticle.E(), 2);
            }
            cout<<"sum : "<<sum<<endl;
            
            
            
            
            
            
            
//            CalcuThrustMCP(MCPVisibleParticleExcISR, MCPVisExcISR);
//            MCPVisExcISRThrust          = MCPVisExcISR.at(0);
//            MCPVisExcISRHemiMass1       = MCPVisExcISR.at(1);
//            MCPVisExcISRHemiMass2       = MCPVisExcISR.at(2);
//            MCPVisExcISRHemiBroadening1 = MCPVisExcISR.at(3);
//            MCPVisExcISRHemiBroadening2 = MCPVisExcISR.at(4);
//            MCPVisExcIsrSphericity = CalcuSphericityMCP(MCPVisibleParticleExcISR);

//            CalcuThrustMCP(MCPUsedForThrust, MCPExcISR);
//            MCPExcISRThrust          = MCPExcISR.at(0);
//            MCPExcISRHemiMass1       = MCPExcISR.at(1);
//            MCPExcISRHemiMass2       = MCPExcISR.at(2);
//            MCPExcISRHemiBroadening1 = MCPExcISR.at(3);
//            MCPExcISRHemiBroadening2 = MCPExcISR.at(4);
//            MCPExcISRSphericity = CalcuSphericityMCP(MCPUsedForThrust);

            
//            cout<<"MCPVisISRThrust : MCPVisExcISRThrust : "<<MCPVisISRThrust<<" : "<<MCPVisExcISRThrust<<endl;
//            cout<<"MCPExcISRThrust : "<<MCPExcISRThrust<<endl;
            
//            cout<<"MCPExcISRHemiMass1 : "<<MCPExcISRHemiMass1<<" MCPExcISRHemiMass2 : "<<MCPExcISRHemiMass2<<" MCPExcISRHemiBroadening1 : "<<MCPExcISRHemiBroadening1<<" MCPExcISRHemiBroadening2 : "<<MCPExcISRHemiBroadening2<<endl;
            
//            cout<<"MCPVisISRSphericity : MCPVisExcIsrSphericity : "<<MCPVisISRSphericity<<" : "<<MCPVisExcIsrSphericity<<endl;
//            cout<<"MCPExcISRSphericity : "<<MCPExcISRSphericity<<endl;
            
            
            
            for(int i = 0; i<MCPVisibleParticleExcISR.size(); i++){
                MCParticle* temp = MCPVisibleParticleExcISR.at(i);
                TLorentzVector TLtemp(temp->getMomentum(), temp->getEnergy());
                TLMCPVisibleParticleExcISR += TLtemp;
            }
            
            En_ISR = TLISR.E();
            Mass_ISR = TLISR.M();
            En_MCPVis = TLMCPVisibleParticle.E();
            cout<<"the energy of final state particles is "<<TLFinalStateParticle.E()<<endl;
            cout<<"the invariant mass of final state particles is "<<TLFinalStateParticle.M()<<endl;
            cout<<"the energy of visible final state particles is "<<TLMCPVisibleParticle.E()<<endl;
            cout<<"the energy of neutrino is "<<TLMCPNeutrino.E()<<endl;
            cout<<"the energy of ISR is "<<TLISR.E()<<endl;
            cout<<"the energy of visible final state particles exclude ISR is "<<TLMCPVisibleParticleExcISR.E()<<endl;
            cout<<"the invariant mass of visible final state particles exclude ISR is "<<TLMCPVisibleParticleExcISR.M()<<endl;
            
            
            
            cout<<"ISR : "<<ISR.size()<<endl;
            cout<<"MCPQuark : "<<MCPQuark.size()<<endl;
            cout<<"FinalStateParticle : "<<FinalStateParticle.size()<<endl;
            cout<<"num_charge : "<<MCPCharge.size()<<" num_neutron : "<<MCPNeutron.size()<<endl;
            cout<<"num_visibleEcxISR : "<<MCPVisibleParticleExcISR.size()<<endl;

        
        } catch (lcio::DataNotAvailableException err) {  }
	
        _outputTree->Fill();
        Num ++;
	}  	  

}	



void eventShape::end()
{

	if (_outputTree) {

		TFile *tree_file = _outputTree->GetCurrentFile(); //just in case we switched to a new file
		//tree_file->cd();
		tree_file->Write();
		delete tree_file;
		//tree_file->Close();
	}

}



