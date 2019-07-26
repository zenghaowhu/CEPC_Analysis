#ifndef _ClusterAna_hh_
#define _ClusterAna_hh_

//#include <RConfigure.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <marlin/Processor.h>
#include <EVENT/CalorimeterHit.h>
#include <TNtuple.h>
#include <TObject.h>
#include "TLorentzVector.h"

#include <TTree.h>
#include <TFile.h>
#include <TH3.h>
class TTree;

class ClusterAna  : public marlin::Processor
{
	public:

		Processor*  newProcessor() { return new ClusterAna ; }

		ClusterAna();

		~ClusterAna() {};

		void init();

		void processEvent( LCEvent * evtP );

		void end();

	protected:
		std::string _treeFileName;
		std::string _treeName;
		std::string _colName;

		int _overwrite;
		TTree *_outputTree;

		unsigned int _eventNr;
		int _NumCluster;
		int _ClusterIndex; 
		int _CluSize; 
		float _T0, _AvT_10per, _AvT_half, _Te; 
		float _Pos[3];
		float _NaiveClusterEn, _ClusterFD, _ClusterDepth; 
		float _TotalCluEn,_EcalEn,_HcalEn;

		float _VisibleEn, _VisibleMass; 

		float _DisPro, _MCPEn;
                int _MCPID;
                float _MCP_P[3];

		std::string _fileName;
		std::ostream *_output;
		std::string _histFileName;

      float _shenT0,_shenClusterFD,_shenDepth;
            float _shenCluSize;
            int _shenid;

   
   float _shenEcalnhit,_shenHcalnhit,_shenEcalen,_shenHcalen;


    float _shenEoH;float _shenmaxDepth;float _shenminDepth;
            float _shenroE;float _shenRoE;float _shenpoE;float _shen10oE;
float _shenEE;
float _shengraDepth;
float _shenFD_ECAL,_shenFD_all,_shenFD_HCAL,_shenFD_ECALF10,_shenNH_ECALF10,_shenFD_ECALL20,_shenNH_ECALL20;

float _shenmaxDisHtoL,_shenminDisHtoL;
float _shencluDepth2,_shengraAbsDepth,_shenavDisHtoL,_shenavEnDisHtoL;
int _iscon;
};

#endif


