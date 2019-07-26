#ifndef _JetClustering_hh_
#define _JetClustering_hh_

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

class JetClustering  : public marlin::Processor
{
	public:

		Processor*  newProcessor() { return new JetClustering ; }

		JetClustering();

		~JetClustering() {};

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
        int eventType;
    
    int GeventType;
    int count_u, count_c;
    float mass_B1, mass_B2, mass_G1, mass_G2, mass_R1, mass_R2;
    double misMinEn;
    double ISR_En;
    int count_boson;
    int Pselect_combi;
    float En_minPart1;
    float En_minPart2;
    
    std::vector<Double_t> vB_Px;
    std::vector<Double_t> vB_Py;
    std::vector<Double_t> vB_Pz;
    std::vector<Double_t> vB_En;
    std::vector<Double_t> vG_Px;
    std::vector<Double_t> vG_Py;
    std::vector<Double_t> vG_Pz;
    std::vector<Double_t> vG_En;
    std::vector<Double_t> vR_Px;
    std::vector<Double_t> vR_Py;
    std::vector<Double_t> vR_Pz;
    std::vector<Double_t> vR_En;
    std::vector<Double_t> vP_Px;
    std::vector<Double_t> vP_Py;
    std::vector<Double_t> vP_Pz;
    std::vector<Double_t> vP_En;
    
    
    double CParameter;
    double hemiMass1, hemiMass2;
    double hemiBroadening1, hemiBroadening2;
    double T;
    double Q1BackThrust;
    double Q1Thrust;
    
    
		std::string _fileName;
		std::ostream *_output;
		std::string _histFileName;
};

#endif


