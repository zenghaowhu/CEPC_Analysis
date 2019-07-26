#ifndef _RecoAna_hh_
#define _RecoAna_hh_

//#include <RConfigure.h>
#include <string>
#include <iostream>
#include <fstream>
#include <marlin/Processor.h>
#include <TNtuple.h>
#include <TObject.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>

#include <TTree.h>
#include <TFile.h>

class TTree;

class RecoAna  : public marlin::Processor
{
	public:

		Processor*  newProcessor() { return new RecoAna ; }

		RecoAna();

		~RecoAna() {};

		void init();

		void processEvent( LCEvent * evtP );

		void end();

	protected:
		std::string _treeFileName;
		std::string _treeName;
		std::string _colName;
		std::string _colAdcVals;

		int _Num;
		int _overwrite;
		int _FLAGDIAG; 
		std::vector<float> _ParaOpt; 
		TTree *_outputEvt;
		TTree *_outputMCP; 
		TTree *_outputPFO; 

		unsigned int _eventNr;
		int _PID, _GenStatus, _SimStatus, _parentnum, _daughternum, _Type, _nMCP, _nSelMCP, _nPFOs;
		float _mass, _Charge, _energy, _Px, _Py, _Pz, _CluEn, _CluEnCorr; 
		float _Posx, _Posy, _Posz, _Vex_x, _Vex_y, _Vex_z; 
		int _NNePFO, _NChPFO; 
		int _PairType, _NClu, _NHit;
		float _RecoDepth, _RecoFD, _RecoT0;
		float _DisCC, _CCTrkEn, _CCCluEn, _CCEnB; 
		int _CCType; 	
		float _TotalEnRecoP, _TotalEnMCP, _Mass_MCP, _Mass_RecoP; 
		float _P_RecoP[3], _P_MCP[3];
		int _PairIndex, _MCPID, _RecoPID;
		float _MCPEn, _MCPCharge, _RecoPEn, _RecoPCharge, _DR, _DE; 
		float _RecoPP[3], _MCPP[3];	
		float _ISRPt, _ISREn, _N3Pt, _N3En;
		float _MaxCT; 
		float _Mass_a; 
		int count;
		std::string _fileName;
		std::ostream *_output;
};



#endif


