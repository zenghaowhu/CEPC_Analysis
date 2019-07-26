#include <ClusterAna.hh>
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/CalorimeterHit.h>
#include <EVENT/SimCalorimeterHit.h>
#include <values.h>
#include <string>
#include <iostream>
#include <cmath>
#include <TMath.h>
#include <vector>
#include <EVENT/LCFloatVec.h>
#include <EVENT/LCParameters.h>
#include <EVENT/Cluster.h>
#include <EVENT/LCRelation.h>
#include <stdexcept>
#include <TFile.h> 
#include <TTree.h>
#include <Rtypes.h> 
#include <sstream>		
#include <TH1.h>
#include <TVector3.h>
#include "UTIL/CellIDDecoder.h"



#include "gear/BField.h"
#include "gear/CalorimeterParameters.h"
#include <marlin/Global.h>



using namespace std;
/*
const float BField = 3.5;
//ILD_Default
const float DHCALBarrelRadius = 2058.0;         //Octo
const float DHCALEndCapInnerZ = 2650.0;
const float ECALEndCapInnerZ = 2450.0;
const float ECALHalfZ = 2350.0; // mm, Endcap Ente
const float ECALRadius = 1847.4; // mm... minimal part for the octagnle.
const float TPCRadius = 1808.0 ;

const float TPCOuterRadius = 1808.0; 
const float TPCInnerRadius = 325.0; 
const float LStar = 1500.0; 
const string ECALCellIDDecoder = "M:3,S-1:3,I:9,J:9,K-1:6";
const double pi = acos(-1.0);

*/
//const float zmaxBarrelEcal=2350.0;
//const float rminBarrelEcal=1843;
//const float zminEndCapEcal=2450;
const string ECALCellIDDecoder = "M:3,S-1:3,I:9,J:9,K-1:6";
const float ECALHalfZ=2350.0;
const float ECALRadius=1843;
const double pi = acos(-1.0); 

	float rminBarrelEcal;
	float rmaxBarrelEcal;
	float zminBarrelEcal;
	float zmaxBarrelEcal;
	
	float rminEndCapEcal;
	float rmaxEndCapEcal;
	float zminEndCapEcal;
	float zmaxEndCapEcal;
	
	float rminBarrelHcal;
	float rmaxBarrelHcal;
	float zminBarrelHcal;
	float zmaxBarrelHcal;
	
	float rminEndCapHcal;
	float rmaxEndCapHcal;
	float zminEndCapHcal;
	float zmaxEndCapHcal;



//ILD_v2
struct ClusterMeasure {

	int ClusterSize;
	float T_zero; 
	float AverageTime_10Percent;
	float AverageTime_Half; 
	float T_end; 
	float NaiveClusterEn; 
	float ClusterFD; 
	float ClusterDepth; 

	TVector3 SeedPos;
	
	// Define Cluster Nature As the Course of the Most Leading Component of the Most Energetic/Inner Most Hit??
} ;

ClusterAna aClusterAna ;
ClusterAna::ClusterAna()
	: Processor("ClusterAna"),
	_output(0)
{
	_description = "Print Cluster Information" ;

	_treeFileName="MCTruth.root";
	registerProcessorParameter( "TreeOutputFile" , 
			"The name of the file to which the ROOT tree will be written" ,
			_treeFileName ,
			_treeFileName);

	_treeName="BHit";
	registerProcessorParameter( "TreeName" , 
			"The name of the ROOT tree" ,
			_treeName ,
			_treeName);

	_overwrite=0;
	registerProcessorParameter( "OverwriteFile" , 
			"If zero an already existing file will not be overwritten." ,
			_overwrite ,
			_overwrite);

}

std::vector<int>SortMeasure( std::vector<float> Measure, int ControlOrderFlag)
{
	//ControlOrderFlag = 0, Ascend order; 1, descend order

	std::vector<int> objindex; 
	int Nobj = Measure.size();

	for(int k = 0; k < Nobj; k++)
	{
		objindex.push_back(k);
	}

	int FlagSwapOrder = 1;
	float SwapMeasure = 0;
	int SwapIndex = 0; 

	for(int i = 0; i < Nobj && FlagSwapOrder; i++)
	{
		FlagSwapOrder = 0;
		for(int j = 0; j < Nobj - 1; j++)
		{
			if((Measure[j] < Measure[j+1] && ControlOrderFlag) || (Measure[j] > Measure[j+1] && !ControlOrderFlag) )
			{
				FlagSwapOrder = 1; 
				SwapMeasure = Measure[j];
				Measure[j] = Measure[j+1];
				Measure[j+1] = SwapMeasure; 

				SwapIndex = objindex[j];
				objindex[j] = objindex[j+1];
				objindex[j+1] = SwapIndex; 
			}
		}
	}

	return objindex; 
}


void gearPara(){

	const gear::CalorimeterParameters &ecalBarrelParameters = marlin::Global::GEAR->getEcalBarrelParameters();
	const gear::CalorimeterParameters &ecalEndCapParameters = marlin::Global::GEAR->getEcalEndcapParameters();
	const gear::CalorimeterParameters &hcalBarrelParameters = marlin::Global::GEAR->getHcalBarrelParameters();
	const gear::CalorimeterParameters &hcalEndCapParameters = marlin::Global::GEAR->getHcalEndcapParameters();
	
	 rminBarrelEcal=(ecalBarrelParameters.getExtent()[0]);
	 rmaxBarrelEcal=(ecalBarrelParameters.getExtent()[1]);
	 zminBarrelEcal=(ecalBarrelParameters.getExtent()[2]);
	 zmaxBarrelEcal=(ecalBarrelParameters.getExtent()[3]);
	
	 rminEndCapEcal=(ecalEndCapParameters.getExtent()[0]);
	 rmaxEndCapEcal=(ecalEndCapParameters.getExtent()[1]);
	 zminEndCapEcal=(ecalEndCapParameters.getExtent()[2]);
	 zmaxEndCapEcal=(ecalEndCapParameters.getExtent()[3]);
	
	 rminBarrelHcal=(hcalBarrelParameters.getExtent()[0]);
	 rmaxBarrelHcal=(hcalBarrelParameters.getExtent()[1]);
	 zminBarrelHcal=(hcalBarrelParameters.getExtent()[2]);
	 zmaxBarrelHcal=(hcalBarrelParameters.getExtent()[3]);
	
	 rminEndCapHcal=(hcalEndCapParameters.getExtent()[0]);
	 rmaxEndCapHcal=(hcalEndCapParameters.getExtent()[1]);
	 zminEndCapHcal=(hcalEndCapParameters.getExtent()[2]);
	 zmaxEndCapHcal=(hcalEndCapParameters.getExtent()[3]);
	
}

int SubDeFlag(TVector3 inputPos){


	int FlagD(-1);
	if( fabs(inputPos[2]) > zmaxEndCapHcal || fabs(inputPos.Perp()) > rmaxBarrelHcal )
	{
		FlagD = 2;
	}
	else if( fabs(inputPos[2]) > zminEndCapHcal || fabs(inputPos.Perp()) > rminBarrelHcal)
	{
		FlagD = 1;          // Position outsider than DHCAL Region
	}
	else if( fabs(inputPos[2]) > zmaxEndCapEcal || fabs(inputPos.Perp()) > rmaxBarrelEcal)
	{
		FlagD = 10;
	}
	else if( fabs(inputPos[2]) > zminEndCapEcal || fabs(inputPos.Perp()) > rminBarrelEcal)
	{
		FlagD = 0;
	}

	else
	{
		FlagD = 11;         // Position inside Calo... Problematic for Seeds... But could be PreShower hits.
	}

	return FlagD;

}




int NHScaleV3( const std::string& encoder_str, Cluster * clu0, int RatioX, int RatioY, int RatioZ )
{

	int ReScaledNH = 0;
	int NumHit = clu0->getCalorimeterHits().size();
	int tmpI = 0;
	int tmpJ = 0;
	int tmpK = 0;
	float tmpEn = 0;
	int NewCellID0 = 0;
	int NewCellID1 = 0;

	CellIDDecoder<CalorimeterHit> idDecoder(encoder_str);      //Input Hits here refer to AllCleanHits collection

	std::map <double, float> testIDtoEnergy;
	double testlongID = 0;

	for(int i = 0; i < NumHit; i++)
	{
		CalorimeterHit *hit = dynamic_cast<CalorimeterHit*>( clu0->getCalorimeterHits()[i]);

		tmpI = idDecoder(hit)["I"]/RatioX;
		tmpJ = idDecoder(hit)["J"]/RatioY;
		tmpK = (idDecoder(hit)["K-1"]+1)/RatioZ;
		tmpEn = hit->getEnergy();

		NewCellID0 = (tmpK<<24) + (tmpJ<<12) + tmpI;

		testlongID = NewCellID1*1073741824 + NewCellID0;
		if(testIDtoEnergy.find(testlongID) == testIDtoEnergy.end() )
		{
			testIDtoEnergy[testlongID] = tmpEn;
		}
		else
		{
			testIDtoEnergy[testlongID] += tmpEn;
		}
	}

	ReScaledNH = testIDtoEnergy.size();

	return ReScaledNH;

}

int NHScaleV2( const std::string& encoder_str, std::vector<CalorimeterHit*> clu0, int RatioX, int RatioY, int RatioZ )
{

	int ReScaledNH = 0;
	int NumHit = clu0.size();
	int tmpI = 0;
	int tmpJ = 0;
	int tmpK = 0;
	float tmpEn = 0;
	int NewCellID0 = 0;

	CellIDDecoder<CalorimeterHit> idDecoder(encoder_str);      //Input Hits here refer to AllCleanHits collection

	std::map <double, float> testIDtoEnergy;

	for(int i = 0; i < NumHit; i++)
	{
		CalorimeterHit *hit = dynamic_cast<CalorimeterHit*>( clu0[i]);

		tmpI = idDecoder(hit)["I"]/RatioX;
		tmpJ = idDecoder(hit)["J"]/RatioY;
		tmpK = (idDecoder(hit)["K-1"]+1)/RatioZ;
		tmpEn = hit->getEnergy();

		NewCellID0 = (tmpK<<24) + (tmpJ<<12) + tmpI;

		if(testIDtoEnergy.find(NewCellID0) == testIDtoEnergy.end() )
		{
			testIDtoEnergy[NewCellID0] = tmpEn;
		}
		else
		{
			testIDtoEnergy[NewCellID0] += tmpEn;
		}
	}

	ReScaledNH = testIDtoEnergy.size();

	return ReScaledNH;
}



float FDV2(std::vector<CalorimeterHit*> clu, const std::string& encoder_str)
{
	float FractalDim = 0;
	int NReSizeHit[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	int Scale[10] = {20, 30, 40, 50, 60, 70, 80, 90, 100, 200};
	int OriNHit = clu.size();

	for(int j = 0; j < 10; j++)
	{
		NReSizeHit[j] = NHScaleV2(encoder_str, clu, Scale[j], Scale[j], 1);
		FractalDim += 0.1 * TMath::Log(float(OriNHit)/NReSizeHit[j])/TMath::Log(float(Scale[j]));
	}
	if(clu.size() == 0)
		FractalDim = -1;
	return FractalDim;
}



float FDV3( Cluster * clu, const std::string& encoder_str )
{
	float FractalDim = -1;
	int NReSizeHit[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	int Scale[5] = {2, 3, 4, 5, 6};
	int OriNHit = clu->getCalorimeterHits().size();

	if(clu->getCalorimeterHits().size() > 0)
	{
		FractalDim = 0.0;
		for(int j = 0; j < 5; j++)
		{
			NReSizeHit[j] = NHScaleV3(encoder_str, clu, Scale[j], Scale[j], 1);
			FractalDim += 0.2 * TMath::Log(float(OriNHit)/NReSizeHit[j])/TMath::Log(float(Scale[j]));
		}
	}

	return FractalDim;
}


float ClusterT0(Cluster * a_Clu)
{
	float T0 = 1E9; 
	float tmpTime = 0; 
	TVector3 CluHitPos; 
	for(unsigned int i = 0; i < a_Clu->getCalorimeterHits().size(); i++)
	{
		CalorimeterHit * a_hit = a_Clu->getCalorimeterHits()[i];
		CluHitPos = a_hit->getPosition();
		tmpTime = a_hit->getTime() - 1.0/300*CluHitPos.Mag();
		if(tmpTime < T0 && tmpTime > 0)
		{
			T0 = tmpTime; 
		}
	}
	return T0;
}


TVector3 ClusterCoG(Cluster * inputCluster)
{
	TVector3 CenterOfGravity; 

	int inputClusterSize = inputCluster->getCalorimeterHits().size();

	TVector3 tmphitPos; 
	float tmphitEnergy;
	float sumhitEnergy = 0; 

	for(int i = 0; i < inputClusterSize; i++)
	{
		CalorimeterHit * tmpHit = inputCluster->getCalorimeterHits()[i];
		tmphitPos = tmpHit->getPosition();
		tmphitEnergy = tmpHit->getEnergy();

		CenterOfGravity += tmphitPos*tmphitEnergy;
		sumhitEnergy += tmphitEnergy; 
	}

	CenterOfGravity = 1.0/sumhitEnergy * CenterOfGravity; 

	return CenterOfGravity; 
}


float DisSeedSurfaceHit( TVector3 SeedPos )
{

	float DisSS = 0;

	if( fabs(SeedPos[2]) > zmaxBarrelEcal )         //EcalEndcap hit start from 2350 + 100 = 2450
	{

		if( SeedPos.Perp() > rminBarrelEcal )
		{
			if( fabs(SeedPos[2])/SeedPos.Perp() > (zminEndCapEcal+3)/(rminBarrelEcal+ 100) )
			{
				DisSS = ( fabs(SeedPos[2]) - zminEndCapEcal - 3 ) * SeedPos.Mag()/fabs(SeedPos[2]);
			}
			else
			{
				DisSS = (SeedPos.Perp() - rminBarrelEcal - 100 )*SeedPos.Mag()/SeedPos.Perp();
			}
		}
		else
		{
			DisSS = fabs(SeedPos[2]) - zminEndCapEcal - 3;
		}

	}
	else if( SeedPos.Perp() > rminBarrelEcal + 400 )
	{
		DisSS = SeedPos.Perp() - rminBarrelEcal - 100;
	}
	else if( (SeedPos.Phi() > 0 && int(SeedPos.Phi() * 4/TMath::Pi() + 0.5) % 2 == 0 ) || (SeedPos.Phi() < 0 && int(SeedPos.Phi() * 4/TMath::Pi() + 8.5) % 2 == 0 ))
	{
		DisSS = min( fabs(fabs(SeedPos[0]) - rminBarrelEcal), fabs(fabs(SeedPos[1]) - rminBarrelEcal ) );
	}
	else
	{
		DisSS = min( fabs(fabs(SeedPos[0] + SeedPos[1])/TMath::Pi() -rminBarrelEcal), fabs(fabs(SeedPos[0] - SeedPos[1])/TMath::Pi() - rminBarrelEcal) );
	}

	return DisSS;
}


float DisSeedSurface( TVector3 SeedPos )	//ECAL, HCAL, EndCapRing...
{

	float DisSS = 0;

	if( fabs(SeedPos[2]) > ECALHalfZ )         //EcalEndcap hit start from 2350 + 100 = 2450
	{

		if( SeedPos.Perp() > ECALRadius )
		{
			if( fabs(SeedPos[2])/SeedPos.Perp() > (ECALHalfZ + 103)/(ECALRadius + 100) )
			{
				DisSS = ( fabs(SeedPos[2]) - ECALHalfZ - 103 ) * SeedPos.Mag()/fabs(SeedPos[2]);
			}
			else
			{
				DisSS = (SeedPos.Perp() - ECALRadius - 100 )*SeedPos.Mag()/SeedPos.Perp();
			}
		}
		else
		{
			DisSS = fabs(SeedPos[2]) - ECALHalfZ - 103;
		}
	}
	else if( SeedPos.Perp() > ECALRadius + 400 )
	{
		DisSS = SeedPos.Perp() - ECALRadius - 100;
	}
	else if( (SeedPos.Phi() > 0 && int(SeedPos.Phi() * 4/pi + 0.5) % 2 == 0 ) || (SeedPos.Phi() < 0 && int(SeedPos.Phi() * 4/pi + 8.5) % 2 == 0 ))
	{
		DisSS = min( fabs(fabs(SeedPos[0]) - ECALRadius), fabs(fabs(SeedPos[1]) - ECALRadius ) );
	}
	else
	{
		DisSS = min( fabs(fabs(SeedPos[0] + SeedPos[1])/1.414214 -ECALRadius), fabs(fabs(SeedPos[0] - SeedPos[1])/1.414214 - ECALRadius) );
	}

	return DisSS;
}


ClusterMeasure MyMeasureCluster(Cluster * a_Clu)
{
	ClusterMeasure b_Measure; 
	int CluSize = a_Clu->getCalorimeterHits().size();
    cout<<"CluSize : "<<CluSize<<endl;
	b_Measure.ClusterSize = CluSize;
	b_Measure.NaiveClusterEn = a_Clu->getEnergy();
	b_Measure.SeedPos = a_Clu->getPosition();
	b_Measure.AverageTime_10Percent = 0;
        b_Measure.AverageTime_Half = 0;

	std::vector<float> HitTime;
	TVector3 CluHitPos; 
	TVector3 CluCoG=ClusterCoG(a_Clu); 
	for(int i = 0; i < CluSize; i++)
	{
		CalorimeterHit * a_hit = a_Clu->getCalorimeterHits()[i];
		CluHitPos = a_hit->getPosition();
		float tmpTime = a_hit->getTime() - 1.0/300*CluHitPos.Mag();
		HitTime.push_back(tmpTime);
	}
	std::vector<int> HitSeq = SortMeasure(HitTime, 0);

	for(int j = 0; j < CluSize/2; j++)
	{
		if(CluSize > 9)
		{
			int Weight1 = CluSize/10; 
			if(j < Weight1)
				b_Measure.AverageTime_10Percent += 1.0/Weight1*HitTime[HitSeq[j]];
		}
		int Weight2 = CluSize/2; 
		b_Measure.AverageTime_Half += 1.0/Weight2*HitTime[HitSeq[j]];
	}

	b_Measure.T_zero = HitTime[HitSeq[0]]; 
	b_Measure.T_end = HitTime[HitSeq[CluSize - 1]];
	TVector3 CluPos = a_Clu->getPosition();
	b_Measure.ClusterFD = FDV3(a_Clu, ECALCellIDDecoder);
	b_Measure.ClusterDepth =DisSeedSurface(CluHitPos); 

	return b_Measure; 
}

TVector3 ImpactPoint( TVector3 MCP_Mom, int MCP_Charge, float Track_Half_Z, float Track_Radius, float B_Field ) //B_Field always along Z+ direction
{
        TVector3 ImpactP(0, 0, 0);
        float impact_En = MCP_Mom.Mag();
        float impact_Pt = MCP_Mom.Perp();
        float impact_Phi = MCP_Mom.Phi();
        float ScaleFactor = 1.0;
        float Tau = 0;
        float R_helix = 1000*impact_Pt/(0.3*B_Field);
        float Ratio_Z_Pz = 1000/(0.3*B_Field);

        if(MCP_Charge == 0)
        {
                if( fabs(MCP_Mom[2]/impact_En) > Track_Half_Z/sqrt(Track_Half_Z*Track_Half_Z + Track_Radius*Track_Radius) )     // Endcap
                {
                        ScaleFactor = fabs(Track_Half_Z/MCP_Mom[2]);
                }
                else
                {
                        ScaleFactor = Track_Radius/impact_Pt;
                }
                ImpactP.SetXYZ( ScaleFactor*MCP_Mom.X(), ScaleFactor*MCP_Mom.Y(), ScaleFactor*MCP_Mom.Z() );
        }
        else{
                Tau = std::min( double(fabs(1.0/Ratio_Z_Pz*Track_Half_Z/MCP_Mom.Z())), double(fabs(Track_Radius/R_helix)));
                ImpactP.SetXYZ( R_helix*(MCP_Charge*sin(MCP_Charge*Tau - impact_Phi) + MCP_Charge*sin(impact_Phi)), R_helix*(MCP_Charge*cos(MCP_Charge*Tau - impact_Phi) - MCP_Charge*cos(impact_Phi)), Ratio_Z_Pz*MCP_Mom.Z()*Tau);
        }

        return ImpactP;
}


void ClusterAna::init() {

	printParameters();

	TFile *tree_file=new TFile(_treeFileName.c_str(),(_overwrite ? "RECREATE" : "UPDATE"));

	if (!tree_file->IsOpen()) {
		delete tree_file;
		tree_file=new TFile(_treeFileName.c_str(),"NEW");
	}

	_outputTree = new TTree(_treeName.c_str(),_treeName.c_str());
	_outputTree->SetAutoSave(32*1024*1024);  

	_outputTree->Branch("EventNr", &_eventNr, "EventNr/I");
	_outputTree->Branch("NumCluster",&_NumCluster,"NumCluster/I");
	_outputTree->Branch("VisibleEn", &_VisibleEn, "VisibleEn/F");
	_outputTree->Branch("EcalEn", &_EcalEn, "EcalEn/F");
	_outputTree->Branch("HcalEn", &_HcalEn, "HcalEn/F");
	_outputTree->Branch("TotalCluEn", &_TotalCluEn, "TotalCluEn/F");
	_outputTree->Branch("VisibleMass", &_VisibleMass, "VisibleMass/F");


	_outputTree->Branch("ClusterIndex", &_ClusterIndex, "ClusterIndex/I");
	_outputTree->Branch("CluSize", &_CluSize, "CluSize/I");
	_outputTree->Branch("shenCluSize", &_shenCluSize, "shenCluSize/F");
	_outputTree->Branch("T0", &_T0, "T0/F");
	_outputTree->Branch("shenT0", &_shenT0, "shenT0/F");
	_outputTree->Branch("AvT_10per", &_AvT_10per, "AvT_10per/F");
	_outputTree->Branch("AvT_half", &_AvT_half, "AvT_half/F");
	_outputTree->Branch("Te", &_Te, "Te/F");
	_outputTree->Branch("Pos", _Pos, "Pos[3]/F");
	_outputTree->Branch("NaiveClusterEn", &_NaiveClusterEn, "NaiveClusterEn/F");
	_outputTree->Branch("ClusterFD", &_ClusterFD, "ClusterFD/F");
	_outputTree->Branch("shenClusterFD", &_shenClusterFD, "shenClusterFD/F");
	_outputTree->Branch("ClusterDepth", &_ClusterDepth, "ClusterDepth/F");
	_outputTree->Branch("shenDepth", &_shenDepth, "shenDepth/F");
	_outputTree->Branch("shenid", &_shenid, "shenid/I");

	//MCP that most probably induces this cluster
	_outputTree->Branch("MCPID", &_MCPID, "MCPID/I");
	_outputTree->Branch("MCPEn", &_MCPEn, "MCPEn/F");
	_outputTree->Branch("MCP_P", _MCP_P, "MCP_P[3]/F");
	_outputTree->Branch("DisPro", &_DisPro, "DisPro/F");


    _outputTree->Branch("shenEcalnhit", &_shenEcalnhit, "shenEcalnhit/F");
   _outputTree->Branch("shenHcalnhit", &_shenHcalnhit, "shenHcalnhit/F");
   _outputTree->Branch("shenEcalen", &_shenEcalen, "shenEcalen/F");
   _outputTree->Branch("shenHcalen", &_shenHcalen, "shenHcalen/F");
   
   _outputTree->Branch("shenEoH", &_shenEoH, "shenEoH/F");
   _outputTree->Branch("shenmaxDepth", &_shenmaxDepth, "shenmaxDepth/F");
   _outputTree->Branch("shenminDepth", &_shenminDepth, "shenminDepth/F");
   _outputTree->Branch("shenroE", &_shenroE, "shenroE/F");
   _outputTree->Branch("shenRoE", &_shenRoE, "shenRoE/F");
   _outputTree->Branch("shenpoE", &_shenpoE, "shenpoE/F");
   _outputTree->Branch("shen10oE", &_shen10oE, "shen10oE/F");
   _outputTree->Branch("shenEE", &_shenEE, "shenEE/F");
   _outputTree->Branch("shengraDepth", &_shengraDepth, "shengraDepth/F");
   _outputTree->Branch("shenFD_ECAL", &_shenFD_ECAL, "shenFD_ECAL/F");
   _outputTree->Branch("shenFD_all", &_shenFD_all, "shenFD_all/F");
   _outputTree->Branch("shenFD_HCAL", &_shenFD_HCAL, "shenFD_HCAL/F");
   _outputTree->Branch("shenFD_ECALF10", &_shenFD_ECALF10, "shenFD_ECALF10/F");
   _outputTree->Branch("shenNH_ECALF10", &_shenNH_ECALF10, "shenNH_ECALF10/F");
   _outputTree->Branch("shenFD_ECALL20", &_shenFD_ECALL20, "shenFD_ECALL20/F");
   _outputTree->Branch("shenNH_ECALL20", &_shenNH_ECALL20, "shenNH_ECALL20/F");
   
   _outputTree->Branch("shenmaxDisHtoL", &_shenmaxDisHtoL, "shenmaxDisHtoL/F");
   _outputTree->Branch("shenminDisHtoL", &_shenminDisHtoL, "shenminDisHtoL/F");
   _outputTree->Branch("shencluDepth2", &_shencluDepth2, "shencluDepth2/F");
   _outputTree->Branch("shengraAbsDepth", &_shengraAbsDepth, "shengraAbsDepth/F");
   _outputTree->Branch("shenavDisHtoL", &_shenavDisHtoL, "shenavDisHtoL/F");
   _outputTree->Branch("shenavEnDisHtoL", &_shenavEnDisHtoL, "shenavEnDisHtoL/F");
   _outputTree->Branch("iscon", &_iscon, "iscon/I");

}

void ClusterAna::processEvent( LCEvent * evtP ) 
{		

	if (evtP) 								
	{	

		try 	
		{    
			_eventNr=evtP->getEventNumber();

            std::cout<<_eventNr<<" events processed"<<std::endl;
            gearPara();

			LCCollection *MCParticleCol = evtP->getCollection("MCParticle");
            int _NMCP = MCParticleCol->getNumberOfElements();
			TVector3 VTX, EndP, MCP_P, ImpactPos, SumP;
			float charge = 0;  
			int PDG = 0; 
			std::vector<MCParticle*> ImpactMCP;	
			std::map<MCParticle*, TVector3> MCP_Impact; 		
			_VisibleEn = 0; 
			_VisibleMass = 0; 
			SumP.SetXYZ(0, 0, 0);
            _iscon=0;
            cout<<"MCPartilce*********"<<endl;
			for(int i0 = 0; i0 < _NMCP; i0++)
			{
				MCParticle * a_MCP = dynamic_cast<MCParticle*>(MCParticleCol->getElementAt(i0));
                int PDG_mcp = a_MCP->getPDG();
//                if(a_MCP->getGeneratorStatus() == 1){
//                    cout<<"mcpartilce pdg is "<<PDG_mcp<<endl;
//                }
                
                int nParents = a_MCP->getParents().size();
                int nDaughters = a_MCP->getDaughters().size();
                if(nParents == 0){
                    cout<<"nParents == 0, pdg is "<<PDG_mcp<<endl;
                    for(int j= 0; j<nDaughters; j++){
                        MCParticle* pi0Daughter = a_MCP->getDaughters()[j];
                        cout<<"pi0 daughter is "<<pi0Daughter->getPDG()<<endl;
                    }
                }
                
                
                if(a_MCP->getPDG()==22 && a_MCP->getParents().size()==0){//shen
                    VTX = a_MCP->getVertex();
                    EndP = a_MCP->getEndpoint();
                    MCP_P = a_MCP->getMomentum();
                    charge = a_MCP->getCharge();
                    PDG = a_MCP->getPDG();
                    if(a_MCP->isDecayedInTracker()){ _iscon=1;}
                    else _iscon=0;
                    continue;
                }
			}
            cout<<"********"<<endl;
            
            cout<<"ArborPFOs*******"<<endl;
            LCCollection* col_pfo = evtP->getCollection("ArborPFOs");
            int nPFO = col_pfo->getNumberOfElements();
            int clusterSize_pfo = 999;
            float clusterFD_pfo = 999;
            int PDG_pfo = -999;
            ClusterMeasure pfo_Measure;
            for(int i = 0; i<nPFO; i++){
                ReconstructedParticle* pfo = dynamic_cast<ReconstructedParticle* >(col_pfo->getElementAt(i));
                PDG_pfo = pfo->getType();
                cout<<"the pdg of pfo is "<<PDG_pfo<<endl;
                
                if(pfo->getClusters().size() > 0 && PDG_pfo != 0){
                    Cluster *a_Clu = pfo->getClusters()[0];
                    pfo_Measure = MyMeasureCluster(a_Clu);
                    clusterSize_pfo = pfo_Measure.ClusterSize;
                    clusterFD_pfo = pfo_Measure.ClusterFD;
                }

                
            }
            
            if(PDG_pfo == 22 && nPFO == 1)
            {
//                pi0toGamma += 1;
                cout<<"clusterFD_pfo : "<<clusterFD_pfo<<" clusterSize_pfo :"<<clusterSize_pfo<<endl;
            }
            cout<<"**********"<<endl;
            
            
            
			LCCollection *ClusterCol = evtP->getCollection("ArborNeutral");
			_NumCluster = ClusterCol->getNumberOfElements(); 
//            cout<<"   shenshen"<<_NumCluster<<endl;
			ClusterMeasure a_Measure;
			TVector3 PosSeed; 
			TVector3 CandiMCPImpactPos; 
			float tmpProjDis = 0; 
			float MinProjDis = 1.0E10; 
			_MCPID=0;
			std::vector<int> cluIndex;
			std::vector<float> cluEn;
			_TotalCluEn=0;
			for(int i = 0; i < _NumCluster; i++)
			{
				Cluster * a_clu = dynamic_cast<Cluster*>(ClusterCol->getElementAt(i));
 //               cout<<"________cluen:"<<a_clu->getEnergy()<<endl;
                cluEn.push_back(a_clu->getEnergy());
				_TotalCluEn+=a_clu->getEnergy();
			}
			cluIndex=SortMeasure(cluEn,1);
			_ClusterDepth=0;
			_shenDepth=-99;
			_ClusterFD=0;
			_shenClusterFD=-99;
            _shenT0=-99;
            _shenCluSize=-99;
            _shenid=-99;

            if(_NumCluster>0){
				Cluster * a_Clu = dynamic_cast<Cluster*>(ClusterCol->getElementAt(cluIndex[0]));//shen

				//Cluster * a_Clu = dynamic_cast<Cluster*>(ClusterCol->getElementAt(cluIndex[i]));
				a_Measure = MyMeasureCluster(a_Clu);
				_CluSize = a_Measure.ClusterSize; 
                _shenCluSize= float(a_Clu->getCalorimeterHits().size());
				//_EcalEn=a_Clu->getSubdetectorEnergies()[0];
				//_HcalEn=a_Clu->getSubdetectorEnergies()[1];
				_EcalEn=0;
				_HcalEn=0;
                _shenid=0;

				for(int j = 0; j < _CluSize; j++)
				{
					CalorimeterHit * a_hit = a_Clu->getCalorimeterHits()[j];
					//cout<<a_hit->getEnergy()<<" "<<abs(a_hit->getEnergy()-0.11)<<endl;
					if(abs(a_hit->getEnergy()-0.11)<1e-6)_HcalEn+=a_hit->getEnergy();
					else
                    _EcalEn+=a_hit->getEnergy();

				}
//				cout<<"hcal "<<_HcalEn<<endl;
//==============================shen lich var====================
                TVector3 HitPos;
                _shenHcalnhit=0;
                _shenEcalnhit=0;
                _shenHcalen=0;
                _shenEcalen=0;
                _shenEoH=0;_shenmaxDepth=-100;_shenminDepth=1E6;
                _shenroE=0;_shenRoE=0;_shenpoE=0;_shen10oE=0;

                float EEClu_r=0,EEClu_R=0,EEClu_p=0,EEClu_L10=0;
                float currDepth=0;

                TVector3 CluPos = a_Clu->getPosition();
                TVector3 IntDir = ClusterCoG(a_Clu)-CluPos;
                std::vector<CalorimeterHit*> Ecalhits;
                std::vector<CalorimeterHit*> Hcalhits;
                std::vector<CalorimeterHit*> allhits;

                std::vector<CalorimeterHit*> Ecalf10hits;
                std::vector<CalorimeterHit*> Ecall20hits;

      
      
                allhits.clear();
                Ecalhits.clear();
                Hcalhits.clear();

                Ecalf10hits.clear();
                Ecall20hits.clear();

                int index1 = 0, index2 = 0;



                for(int j = 0; j < _CluSize; j++){
                    CalorimeterHit * a_hit = a_Clu->getCalorimeterHits()[j];
            
                    CellIDDecoder<CalorimeterHit> idDecoder(ECALCellIDDecoder);
                    int NLayer = idDecoder(a_hit)["K-1"];
                    allhits.push_back(a_hit);

                    HitPos = a_hit->getPosition();

                    currDepth = DisSeedSurfaceHit(HitPos);
                    if(currDepth > _shenmaxDepth)
                    {
                        _shenmaxDepth = currDepth;
                        index1 = j;
                    }
                    if(currDepth < _shenminDepth)
                    {
                        _shenminDepth = currDepth;
                        index2 = j;
                    }

                    float	crdis = (CluPos-HitPos).Mag()*sin((CluPos-HitPos).Angle(IntDir));


                    if(SubDeFlag(HitPos)==1){
                        _shenHcalnhit=_shenHcalnhit+1;
                        _shenHcalen=_shenHcalen+a_hit->getEnergy();
                        Hcalhits.push_back(a_hit);
                    }
                    if(SubDeFlag(HitPos)==0){
                        _shenEcalnhit=_shenEcalnhit+1;
                        _shenEcalen=_shenEcalen+a_hit->getEnergy();
                        Ecalhits.push_back(a_hit);

                        if(crdis < 22) EEClu_R += a_hit->getEnergy();
                        if(crdis < 11) EEClu_r += a_hit->getEnergy();
                        if(crdis < 6) EEClu_p += a_hit->getEnergy();
                        if(NLayer < 10)  EEClu_L10 += a_hit->getEnergy();
                        if(NLayer< 10) Ecalf10hits.push_back(a_hit);
                        else Ecall20hits.push_back(a_hit);
                    }
            }

                if(_shenEcalen!=0){
                    _shenEoH=_shenHcalen/_shenEcalen;
                    _shenroE=EEClu_r/_shenEcalen;
                    _shenRoE=EEClu_R/_shenEcalen;
                    _shenpoE=EEClu_p/_shenEcalen;
                    _shen10oE=EEClu_L10/_shenEcalen;
                }
                _shenEE=_shenEcalen/a_Clu->getEnergy();
                if(a_Clu->getCalorimeterHits().size()>0){
                    CalorimeterHit * maxdis_hit = a_Clu->getCalorimeterHits()[index1];
                    CalorimeterHit * mindis_hit = a_Clu->getCalorimeterHits()[index2];
                    TVector3 maxpos = maxdis_hit->getPosition();
                    TVector3 minpos = mindis_hit->getPosition();
                    _shencluDepth2 = (maxpos-minpos).Mag();
                    TVector3 GraPos = ClusterCoG(a_Clu);
                    _shengraDepth = DisSeedSurfaceHit(GraPos);
                    _shengraAbsDepth = (GraPos-minpos).Mag();
                    _shenmaxDisHtoL = -100;
                    _shenminDisHtoL = 1E6;
                    float totDisHtoL = 0;

                    float totHitEn = 0;
                    float totHitEnDis = 0;
                    float HitEn;


                    for(int s2 = 0; s2 < _CluSize; s2++)
                    {
                        CalorimeterHit * a_hit2 = a_Clu->getCalorimeterHits()[s2];
                        HitPos = a_hit2->getPosition();
                        HitEn  = a_hit2->getEnergy();
                        TVector3 par1 = GraPos-minpos;
                        TVector3 par2 = minpos-HitPos;
                        TVector3 par3 = par1.Cross(par2);
                        float disHtoL = par3.Mag()/par1.Mag();
                        totDisHtoL+=disHtoL;
                        totHitEn+=HitEn;
                        totHitEnDis+=HitEn*disHtoL;
                        if (disHtoL > _shenmaxDisHtoL) _shenmaxDisHtoL = disHtoL;
                        if (disHtoL < _shenminDisHtoL) _shenminDisHtoL = disHtoL;
                    }
      
                    _shenavDisHtoL = totDisHtoL/_CluSize;
                    _shenavEnDisHtoL = totHitEnDis/totHitEn;
                    _shenFD_all = FDV2(allhits, ECALCellIDDecoder);
                    _shenFD_ECAL = FDV2(Ecalhits, ECALCellIDDecoder);
                    _shenFD_HCAL = FDV2(Hcalhits, ECALCellIDDecoder);
        	        _shenFD_ECALF10 = FDV2(Ecalf10hits, ECALCellIDDecoder);
                    _shenNH_ECALF10 = Ecalf10hits.size();
        	        _shenFD_ECALL20 = FDV2(Ecall20hits, ECALCellIDDecoder);
        	        _shenNH_ECALL20 = Ecall20hits.size();
                }

      //===============================
				_T0 = a_Measure.T_zero; 
				_shenT0 = ClusterT0(a_Clu); 
				_AvT_10per = a_Measure.AverageTime_10Percent;
				_AvT_half = a_Measure.AverageTime_Half;
				_Te = a_Measure.T_end;
				PosSeed = a_Measure.SeedPos; 
				_Pos[0] = PosSeed.X();
				_Pos[1] = PosSeed.Y();
				_Pos[2] = PosSeed.Z();
				_NaiveClusterEn = a_Measure.NaiveClusterEn;
				_ClusterFD = a_Measure.ClusterFD;
				_shenClusterFD = FDV3(a_Clu,ECALCellIDDecoder);
                TVector3 PosClu = a_Clu->getPosition();//shen
                _shenDepth = DisSeedSurface(PosClu);
                _ClusterDepth = a_Measure.ClusterDepth;
				_ClusterIndex = 0;//shen 
				//_ClusterIndex = i; 

                if(_shenClusterFD>0.18*TMath::Log(_shenCluSize)-0.53&&_shenClusterFD<0.16*TMath::Log(_shenCluSize)+0.025&&_shenClusterFD>-0.2*TMath::Log(_shenCluSize)+0.4&&((log10(_shenT0)<-2&&log10(_shenDepth)<2&&log10(_shenCluSize)>2)||(log10(_shenT0)<-1.5&&log10(_shenCluSize)<2))) {_shenid=22;}

//                cout<<"shenid---"<<_shenid<<"-----cov---"<<_iscon<<endl;

				//if(fabs(cos(PosSeed.Theta()))>0.95)continue;
				for(unsigned int j = 0; j < ImpactMCP.size(); j++)
				{	
					MCParticle * a_MCP = ImpactMCP[j];
					CandiMCPImpactPos = MCP_Impact[a_MCP];
					tmpProjDis = (PosSeed - CandiMCPImpactPos).Mag()*CandiMCPImpactPos.Angle(PosSeed - CandiMCPImpactPos);
					if(tmpProjDis < MinProjDis)
					{
						MinProjDis = tmpProjDis; 
						_DisPro = tmpProjDis;
						_MCPID = a_MCP->getPDG();
						_MCPEn = a_MCP->getEnergy();
						_MCP_P[0] = a_MCP->getMomentum()[0];
						_MCP_P[1] = a_MCP->getMomentum()[1];
						_MCP_P[2] = a_MCP->getMomentum()[2];

					}
				}

				_outputTree->Fill();
			} //shentest

		}		
		catch (lcio::DataNotAvailableException err) { }

	}  	

}	

void ClusterAna::end()
{

	if (_outputTree) {

		TFile *tree_file = _outputTree->GetCurrentFile(); //just in case we switched to a new file
		tree_file->Write();
		delete tree_file;

	}

}


