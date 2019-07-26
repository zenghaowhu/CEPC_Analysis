#ifndef _MCPJetClustering_hh_
#define _MCPJetClustering_hh_

#include <string>
#include <iostream>
#include <fstream>
#include <marlin/Processor.h>
#include <TNtuple.h>
#include <TObject.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>

class TTree;

class MCPJetClustering  : public marlin::Processor
{
	public:

		Processor*  newProcessor() { return new MCPJetClustering ; }

		MCPJetClustering();

		~MCPJetClustering() {};

		void init();

		void processEvent( LCEvent * evtP );

		void end();

	protected:
		std::string _treeFileName;
		std::string _treeName;
		std::string _colName;
		std::string _colAdcVals;
//		TFile *tree_file;

		int _overwrite;
		TTree *_outputTree;

		unsigned int eventNr;
		int Num;
        float select_j1j2, select_j3j4;
        int eventType;
        float ratio_1, ratio_2;
    float En_WP, En_WM, angle_q1q2, angle_q3q4, angle_j1j2, angle_j3j4;
    float qselect_q1q2, qselect_q3q4;
    int qeventType;
    float newratio_1, newratio_2;
    float ISRPt, NuPt;
    float angle_q1q3, angle_q1q4, angle_q2q3, angle_q2q4;
    float ratio_WP_1, ratio_WP_2, ratio_WM_1, ratio_WM_2;
    float ratio_Z1_1, ratio_Z1_2, ratio_Z2_1, ratio_Z2_2;
    float mass_pair1, mass_pair2;
    float costheta_1, costheta_2, costheta_3, costheta_4;
    float ratio_P1, ratio_P2;
    int num_d, num_b, num_c, num_s, num_t, num_u;
		std::string _fileName;
		std::ostream *_output;
		std::string _histFileName;
};

#endif


