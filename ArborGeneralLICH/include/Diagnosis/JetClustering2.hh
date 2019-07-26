#ifndef _JetClustering2_hh_
#define _JetClustering2_hh_

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

class JetClustering2  : public marlin::Processor
{
	public:

		Processor*  newProcessor() { return new JetClustering2 ; }

		JetClustering2();

		~JetClustering2() {};

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
    int total_charge;
    float Q1CosTheta;
    float Q2CosTheta;
    float angletwoQ;
    float sum1, sum2;
    int jet1charge, jet2charge;
    int PDG_Q1, PDG_Q2;
    float CosThetaSum1, CosThetaSum2;
    int num_jet1charge, num_jet2charge;
    float alpha;
    float alpha1, alpha2;
    
    float Benergy, Bbarenergy, BCenergy, BbarCenergy;
    
    std::vector<Float_t> vOElcPx;
    std::vector<Float_t> vOElcPy;
    std::vector<Float_t> vOElcPz;
    std::vector<Float_t> vOElcE;
    
    std::vector<Float_t> vOMuPx;
    std::vector<Float_t> vOMuPy;
    std::vector<Float_t> vOMuPz;
    std::vector<Float_t> vOMuE;
    
    std::vector<Float_t> vOPosPx;
    std::vector<Float_t> vOPosPy;
    std::vector<Float_t> vOPosPz;
    std::vector<Float_t> vOPosE;
    
    std::vector<Float_t> vOMuPlusPx;
    std::vector<Float_t> vOMuPlusPy;
    std::vector<Float_t> vOMuPlusPz;
    std::vector<Float_t> vOMuPlusE;
    
    std::vector<Float_t> vOKPlusPx;
    std::vector<Float_t> vOKPlusPy;
    std::vector<Float_t> vOKPlusPz;
    std::vector<Float_t> vOKPlusE;
    
    std::vector<Float_t> vOKMinusPx;
    std::vector<Float_t> vOKMinusPy;
    std::vector<Float_t> vOKMinusPz;
    std::vector<Float_t> vOKMinusE;
    
    
    
    std::vector<Float_t> vBElcPx;
    std::vector<Float_t> vBElcPy;
    std::vector<Float_t> vBElcPz;
    std::vector<Float_t> vBElcE;
    std::vector<Float_t> vBPosPx;
    std::vector<Float_t> vBPosPy;
    std::vector<Float_t> vBPosPz;
    std::vector<Float_t> vBPosE;
    
    std::vector<Float_t> vBMuPx;
    std::vector<Float_t> vBMuPy;
    std::vector<Float_t> vBMuPz;
    std::vector<Float_t> vBMuE;
    std::vector<Float_t> vBMuPlusPx;
    std::vector<Float_t> vBMuPlusPy;
    std::vector<Float_t> vBMuPlusPz;
    std::vector<Float_t> vBMuPlusE;
    
    std::vector<Float_t> vBKPlusPx;
    std::vector<Float_t> vBKPlusPy;
    std::vector<Float_t> vBKPlusPz;
    std::vector<Float_t> vBKPlusE;
    std::vector<Float_t> vBKMinusPx;
    std::vector<Float_t> vBKMinusPy;
    std::vector<Float_t> vBKMinusPz;
    std::vector<Float_t> vBKMinusE;
    

    
    std::vector<Float_t> vBCElcPx;
    std::vector<Float_t> vBCElcPy;
    std::vector<Float_t> vBCElcPz;
    std::vector<Float_t> vBCElcE;
    std::vector<Float_t> vBCPosPx;
    std::vector<Float_t> vBCPosPy;
    std::vector<Float_t> vBCPosPz;
    std::vector<Float_t> vBCPosE;
    
    std::vector<Float_t> vBCMuPx;
    std::vector<Float_t> vBCMuPy;
    std::vector<Float_t> vBCMuPz;
    std::vector<Float_t> vBCMuE;
    std::vector<Float_t> vBCMuPlusPx;
    std::vector<Float_t> vBCMuPlusPy;
    std::vector<Float_t> vBCMuPlusPz;
    std::vector<Float_t> vBCMuPlusE;
    
    std::vector<Float_t> vBCKPlusPx;
    std::vector<Float_t> vBCKPlusPy;
    std::vector<Float_t> vBCKPlusPz;
    std::vector<Float_t> vBCKPlusE;
    std::vector<Float_t> vBCKMinusPx;
    std::vector<Float_t> vBCKMinusPy;
    std::vector<Float_t> vBCKMinusPz;
    std::vector<Float_t> vBCKMinusE;
    

    
    std::vector<Float_t> vBbarElcPx;
    std::vector<Float_t> vBbarElcPy;
    std::vector<Float_t> vBbarElcPz;
    std::vector<Float_t> vBbarElcE;
    std::vector<Float_t> vBbarPosPx;
    std::vector<Float_t> vBbarPosPy;
    std::vector<Float_t> vBbarPosPz;
    std::vector<Float_t> vBbarPosE;
    
    std::vector<Float_t> vBbarMuPx;
    std::vector<Float_t> vBbarMuPy;
    std::vector<Float_t> vBbarMuPz;
    std::vector<Float_t> vBbarMuE;
    std::vector<Float_t> vBbarMuPlusPx;
    std::vector<Float_t> vBbarMuPlusPy;
    std::vector<Float_t> vBbarMuPlusPz;
    std::vector<Float_t> vBbarMuPlusE;
    
    std::vector<Float_t> vBbarKPlusPx;
    std::vector<Float_t> vBbarKPlusPy;
    std::vector<Float_t> vBbarKPlusPz;
    std::vector<Float_t> vBbarKPlusE;
    std::vector<Float_t> vBbarKMinusPx;
    std::vector<Float_t> vBbarKMinusPy;
    std::vector<Float_t> vBbarKMinusPz;
    std::vector<Float_t> vBbarKMinusE;
    
    std::vector<Float_t> vBbarCElcPx;
    std::vector<Float_t> vBbarCElcPy;
    std::vector<Float_t> vBbarCElcPz;
    std::vector<Float_t> vBbarCElcE;
    std::vector<Float_t> vBbarCPosPx;
    std::vector<Float_t> vBbarCPosPy;
    std::vector<Float_t> vBbarCPosPz;
    std::vector<Float_t> vBbarCPosE;
    
    std::vector<Float_t> vBbarCMuPx;
    std::vector<Float_t> vBbarCMuPy;
    std::vector<Float_t> vBbarCMuPz;
    std::vector<Float_t> vBbarCMuE;
    std::vector<Float_t> vBbarCMuPlusPx;
    std::vector<Float_t> vBbarCMuPlusPy;
    std::vector<Float_t> vBbarCMuPlusPz;
    std::vector<Float_t> vBbarCMuPlusE;
    
    std::vector<Float_t> vBbarCKPlusPx;
    std::vector<Float_t> vBbarCKPlusPy;
    std::vector<Float_t> vBbarCKPlusPz;
    std::vector<Float_t> vBbarCKPlusE;
    std::vector<Float_t> vBbarCKMinusPx;
    std::vector<Float_t> vBbarCKMinusPy;
    std::vector<Float_t> vBbarCKMinusPz;
    std::vector<Float_t> vBbarCKMinusE;
    
		std::string _fileName;
		std::ostream *_output;
		std::string _histFileName;
};

#endif


