#ifndef _mmhzzTotalInvMass_hh_
#define _mmhzzTotalInvMass_hh_

#include <string>
#include <iostream>
#include <fstream>
#include <marlin/Processor.h>
#include <TNtuple.h>
#include <TObject.h>
#include <TTree.h>
#include <TFile.h>

class TTree;

class mmhzzTotalInvMass  : public marlin::Processor
{
	public:

		Processor*  newProcessor() { return new mmhzzTotalInvMass ; }

		mmhzzTotalInvMass();

		~mmhzzTotalInvMass() {};

		void init();

		void processEvent( LCEvent * evtP );

		void end();

	protected:
		std::string _treeFileName;
		std::string _treeName;
		std::string _colName;
		std::string _colAdcVals;

		int _overwrite, _leptonID;
		float _cmsE; 
		TTree *_outputTree, *_outputPFO;
		float _ISRP[3];
		float _ISREn, _ISRPt; 
		float _N3En, _N3Pt; 
		float _NEn, _NPt;
    float NEn_t, NPt_t;
    float vtxmag, endpmag;
		int type;
		int _HDPID, _OriQuarkID; 
		float _OQDir, _HDir;
		int _NMuP, _NMuM, _NChP, _NChM;
		float _P_MuP[4], _P_MuM[4], _P_DL[4];
		int _EventType; 
		float _InvMass, _RecoilMass; 
		float _J1CosTheta, _J2CosTheta, _MaxJetCosTheta; 
		float _OriJ1CosTheta, _OriJ2CosTheta, _MaxOriJetCosTheta; 
		float _Mass_a, _Mass_a_Pdimuon, _Mass_p, H_Mass_a;
    float mcp_muminus_energy, mcp_muplus_energy, mcp_HDminus_energy, mcp_HDplus_energy;
    
		int _PID1, _PID2;
    int NMuplusCount, NMuminusCount, NNeutronCount;
    float energy_neutron_first, energy_neutron_second, energy_muminus, energy_muplus,energy_muminus_H, energy_muplus_H;
    float Mass_a_E, Mass_a_Px, Mass_a_Py, Mass_a_Pz;
    float angle_minus, angle_plus, H_angle_minus, H_angle_plus;
    float mcp_muminus_E, mcp_muminus_Px, mcp_muminus_Py,  mcp_muminus_Pz;
    float mcp_muplus_E, mcp_muplus_Px, mcp_muplus_Py, mcp_muplus_Pz;
    float Mass_neutron_value;
		float _PL1[4], _PL2[4], _RPL1[4], _RPL2[4], _SM[4], _P_allCharged[4], _P_allNeutral[4], _P_Higgs[4], _P_allReco[4];
		float _Hmass; 
		int _Num;
		int _NHDaug; 
		int _HdaughterPID; 
		int _ZdaughterPID; 
		float _Pz[4], _Ph[4], _PzD1[4], _PzD2[4], _PhD1[4], _PhD2[4], _RPzD1[4], _RPzD2[4], _RPhD1[4], _RPhD2[4];
		float _P[4], _SumP[4], _VisP[4], _MissP[4];
		int _PID, _NFMCP, _MotherFlag, _NNeutrino; 
		float _ENeutrino, _DiPhMass, _DiPhMassCorr; 
//		float _CosTheta, _Phi, _Charge;
		float _Mz, _Mrecoil, _MzReco, _MhReco, _MrecoilReco; 
		float KthEn[7][9];
		unsigned int _eventNr;

		float _Mass_p_Pisr, _Mass_a_Pisr, _Mass_a_Plcal;

		float TotalP_a[4], TotalP_p[4];
		int nCHPFO_a, nCHPFO_p, nNEPFO_a, nNEPFO_p;
		float NeCaloE_a[2], NeCaloE_p[2];
		float ElargeP[2], EequP[2], EsmallP[2];

		float ChP[4], FrP[4], PhP[4], NeP[4], UdP[4], FrPh[4], FrNe[4], KPF[4];
		float _EcalTotalE, _HcalTotalE, _EcalCluE, _HcalCluE, _EcalCluE_p, _HcalCluE_p; 
		float _HcalEn1, _HcalEn2,_HcalEn3,_HcalEn4,_HcalEn5;
		float _EcalEn1, _EcalEn2,_EcalEn3,_EcalEn4,_EcalEn5;

		int Type, Charge;
		float Energy, TrkSumEn; 
		float P[3], CluEnCom[2];
		float CluEn;

		int TrackHit;
		float StartPos[3], EndPos[3];
		float _visE;
    int PID2_decayedfromZ2, PID1_decayedfromZ1, PID2_decayedfromZ1, PID1_decayedfromZ2;
    float mcp_Z1D1_energy, mcp_Z1D2_energy, mcp_Z2D1_energy, mcp_Z2D2_energy;
    float mass_dilepton1_pfo, mass_dilepton2_pfo, mass_recoil1, mass_recoil2;
    float Z1_mass, Z2_mass, Z1_energy, Z2_energy;
    float missing_mass, visible_mass;
    float nnhInvisible_mass, nnhInmissing_mass;
    
		std::string _fileName;
		std::ostream *_output;
		std::string _histFileName;
};

#endif


