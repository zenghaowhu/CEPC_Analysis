#ifndef _eventShape_hh_
#define _eventShape_hh_

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

class eventShape  : public marlin::Processor
{
	public:

		Processor*  newProcessor() { return new eventShape ; }

		eventShape();

		~eventShape() {};

		void init();

		void processEvent( LCEvent * evtP );

		void end();

	protected:
		std::string _treeFileName;
		std::string _treeName;
		std::string _colName;
		std::string _colAdcVals;

		int _overwrite;
		TTree *_outputTree;

		unsigned int eventNr;
		int Num;
    
    int Mark;
    
        //about event shape
    double CParameter;
    double hemiMass1, hemiMass2;
    double hemiBroadening1, hemiBroadening2;
    double T;
    double angleTwoQ;
    double En_quarks;
    double En_ISR;
    double Mass_ISR;
    double Mass_quarks;
    double En_MCPVis;
    double MCPVisISRThrust, MCPVisExcISRThrust, MCPExcISRThrust, arborThrust;
    double MCPVisISRSphericity, MCPVisExcIsrSphericity, MCPExcISRSphericity, arborSphericity;
    double arborHemiMass1, arborHemiMass2, arborHemiBroadening1,arborHemiBroadening2;
    double MCPVisISRRapidity, arborRapidity;
    double MCPVisISRDPara, arborDPara;
    double Pt_arbor, En_arbor;
    int num_94, num_quark;
    double angle_twoPlan;
    double Mag_arbor;
    double Mass_arbor;
    int lightQuark;
    double energyNeutronHadron, energyChargeHadron, energyGamma, avNeutronHadron, avChargeHadron, avGamma, energyLight, avLight;
    int num_chargeHadrom, num_neutronHadron, num_gamma, num_chargeLight;
    
    double energy_aNeuHad, energy_aChgHad, energy_aGama, av_aNeuHad, av_aChgHad, av_aGama, energy_aLight, av_aLight;
    int num_aNeuHad, num_aChgHad, num_aGama, num_aLight;
    
    
    double MCPVisExcISRHemiMass1, MCPVisExcISRHemiMass2, MCPVisExcISRHemiBroadening1, MCPVisExcISRHemiBroadening2;
    double MCPVisISRHemiMass1, MCPVisISRHemiMass2, MCPVisISRHemiBroadening1, MCPVisISRHemiBroadening2;
    double MCPExcISRHemiMass1, MCPExcISRHemiMass2, MCPExcISRHemiBroadening1, MCPExcISRHemiBroadening2;
    
    int num_gluon;
    std::vector<Int_t> vQuarkPDG;
    std::vector<Int_t> vHiggsDausPDG;
    std::vector<Double_t> vHiggsDaus;
    std::vector<Int_t> vHiggsDausQuarks;
    std::vector<Double_t> MCPEEC;
    std::vector<Double_t> ArborEEC;
    std::vector<Double_t> gluonEn;
    std::vector<Double_t> gluonAngle;
    
		std::string _fileName;
		std::ostream *_output;
		std::string _histFileName;
};

#endif


