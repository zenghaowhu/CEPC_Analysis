/****************************************************************
* 			Filename      : 	FilterProcessor.h
*
*			Time	      :		March, 01, 2016, Tuesday
*
* 			Author        : 	Li-bo Liao
*
*			Mail          : 	liaolb@ihep.ac.cn
*
*			Description   :     
***************************************************************/
#ifndef _FilterProcessor_
#define _FilterProcessor_

#include "marlin/Processor.h"
#include "IO/LCWriter.h"
#include "lcio.h"
#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/ReconstructedParticle.h>
#include "UTIL/LCRelationNavigator.h"
#include <string>
#include <vector>
#include <TFile.h>
#include <TTree.h>
using namespace lcio ;
using namespace marlin ;
using namespace std;

class FilterProcessor : public Processor
{
	public:

		virtual Processor* newProcessor() { return new FilterProcessor ; }
		FilterProcessor() ;

		virtual void init();
		virtual void processRunHeader( LCRunHeader* run);
		virtual void processEvent( LCEvent* evtP);
		virtual void end();

	protected:

		LCWriter* _lcWrt_u;
		LCWriter* _lcWrt_d;
		LCWriter* _lcWrt_s;

		std::string _inputMCCollection;
		std::string _lcioOutputFile_u;
		std::string _lcioOutputFile_d;
		std::string _lcioOutputFile_s;
		int _nRun, _nEvt;
		int _nMC;

};

#endif
