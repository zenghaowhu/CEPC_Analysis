#ifndef _HitDis_hh_
#define _HitDis_hh_

//#include <RConfigure.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <marlin/Processor.h>
#include <EVENT/CalorimeterHit.h>
#include <TNtuple.h>
#include <TObject.h>

#include <TTree.h>
#include <TFile.h>
#include <TH3.h>
class TTree;

class HitDis  : public marlin::Processor
{
	public:

		Processor*  newProcessor() { return new HitDis ; }

		HitDis();

		~HitDis() {};

		void init();

		void processEvent( LCEvent * evtP );

		void end();

	protected:
		std::string _treeFileName;
		std::string _treeName;
		std::string _colName;
		std::vector<std::string> _EcalInputHits; 
		std::vector<std::string> _HcalInputHits; 
		std::ostream _output; 

		int _overwrite;
		TTree *_outputTree;	// Used only for single particle analysis
		//TTree *_outputHits;
		
		int _Num, _eventNr; 
		int _MCPID;
		float _MCPEn;
		float _MCP[3];
		
		int _NNeu;
		int _NIsoHit;
		int _NHit_E, _NHit_H;
		int _NCluHit;
		int _NLCHit;
		
		float _En_E, _En_H; 
		float _CluEn;
		float _LCEn;

		float _IsoDis[999];
};

#endif


