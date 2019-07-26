#ifndef _MassVeto_hh_
#define _MassVeto_hh_

#include <string>
#include <iostream>
#include <fstream>
#include <marlin/Processor.h>
#include <TNtuple.h>
#include <TObject.h>
#include <TTree.h>
#include <TFile.h>

class TTree;

class MassVeto  : public marlin::Processor
{
	public:

		Processor*  newProcessor() { return new MassVeto ; }
		MassVeto();
		~MassVeto() {};

		void init();
		void processEvent( LCEvent * evtP );
		void end();

	protected:
		std::string _treeFileName;

		int _overwrite;

		TTree *_outputTree, *_outputPFO;
		float _ISRP[3];
		float _ISREn, _ISRPt; 
		float _N3En, _N3Pt; 
		float _NEn, _NPt; 

		float _J1CosTheta, _J2CosTheta, _MaxJetCosTheta; 
		float _Mass_a; 

		int _Num;
		unsigned int _eventNr;
		float _Mass_a_Pisr;
		float _EcalTotalE, _HcalTotalE;
		float TotalP_a[4];
		int Type, Charge;
		float Energy, TrkSumEn; 
		float P[3];
		float CluEn;

		int TrackHit;
		float StartPos[3], EndPos[3];
		float _visE;

		float MinMomentumDiff,MomentumDiff,VetoMCPEn;
		float MaxLeptonEn;
		float VetoMCPP[3],VetoPFOP[3];
		float _Mass_a_Veto;
		int NeutrinoFlavor;

		std::string _fileName;
		std::ostream *_output;
};

#endif


