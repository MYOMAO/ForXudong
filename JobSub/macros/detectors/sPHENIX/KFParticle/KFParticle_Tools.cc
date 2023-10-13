/*
 * This file is part of KFParticle package
 * Copyright ( C ) 2007-2019 FIAS Frankfurt Institute for Advanced Studies
 *               2007-2019 Goethe University of Frankfurt
 *               2007-2019 Ivan Kisel <I.Kisel@compeng.uni-frankfurt.de>
 *               2007-2019 Maksym Zyzak
 *
 * KFParticle is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * ( at your option ) any later version.
 *
 * KFParticle is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

/*****************/
/* Cameron Dean  */
/*   LANL 2020   */
/* cdean@bnl.gov */
/*****************/

#include "KFParticle_Tools.h"

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxVertex.h>
#include <trackbase_historic/SvtxVertexMap.h>

#include <phool/getClass.h>

//KFParticle stuff
#include <KFParticle.h>
#include <KFVertex.h>

#include <Rtypes.h>
#include <TDatabasePDG.h>
#include <TMatrixD.h>
#include <TMatrixDfwd.h>  // for TMatrixD
#include <TMatrixT.h>     // for TMatrixT, operator*

#include <Eigen/Dense>

#include <algorithm>  // for max, remove, minmax_el...
#include <cmath>      // for sqrt, pow, M_PI
#include <cstdlib>    // for abs, NULL
#include <iostream>   // for operator<<, basic_ostream
#include <iterator>   // for end
#include <map>        // for _Rb_tree_iterator, map
#include <memory>     // for allocator_traits<>::va...


#include <fun4all/SubsysReco.h>

#include <trackbase/ActsTrackingGeometry.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrDefs.h>

#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxVertexMap.h>

#include <tpc/TpcClusterZCrossingCorrection.h>
#include <tpc/TpcDistortionCorrection.h>

#include <Acts/Definitions/Algebra.hpp>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#pragma GCC diagnostic ignored "-Wunused-value"
#include <Acts/Propagator/Propagator.hpp>
#pragma GCC diagnostic pop

#include <Acts/EventData/TrackParameters.hpp>
#include <Acts/Surfaces/CylinderSurface.hpp>
#include <Acts/Utilities/Result.hpp>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include <ActsExamples/EventData/Trajectories.hpp>
#pragma GCC diagnostic pop

#include <Eigen/Dense>
#include <trackbase_historic/ActsTransformations.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <Acts/Propagator/EigenStepper.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <utility>



class TFile;
class TH1D;
class TNtuple;
class TTree;


#include "TFile.h"
#include "TTree.h"
#include "TMath.h"


TFile * KFFile;
TFile * ACTSFile;

TTree * ACTSTree; 
TTree * KFTree; 


int ACTSIndex;

float ACTSPropX;
float ACTSPropY;
float ACTSPropZ;


float ACTSX;
float ACTSY;
float ACTSZ;
float ACTSPx;
float ACTSPy;
float ACTSPz;
float ACTSPt;


int KFIndex;


float KFPrePropX;
float KFPrePropY;
float KFPrePropZ;
float KFPrePropPx;
float KFPrePropPy;
float KFPrePropPz;


float KFPropX;
float KFPropY;
float KFPropZ;
float KFPropPx;
float KFPropPy;
float KFPropPz;



float KFTruthX;
float KFTruthY;
float KFTruthZ;
float KFTruthPx;
float KFTruthPy;
float KFTruthPz;





float KFX;
float KFY;
float KFZ;
float KFPx;
float KFPy;
float KFPz;
float KFPt;


float KFBx;
float KFBy;
float KFBz;

float Hamil;
float Theta;

int Trial;

using BoundTrackParam = const Acts::BoundTrackParameters;
using BoundTrackParamResult = Acts::Result<BoundTrackParam>;
using SurfacePtr = std::shared_ptr<const Acts::Surface>;
using Trajectory = ActsExamples::Trajectories;

class PHCompositeNode;
class SvtxTrack;
class SvtxTrackMap;
class SvtxVertexMap;


std::vector<float> ACTSHitsX;
std::vector<float> ACTSHitsY;
std::vector<float> ACTSHitsZ;


std::vector<float> ACTSHitsPX;
std::vector<float> ACTSHitsPY;
std::vector<float> ACTSHitsPZ;


//#include <trackingdiagnostics/KshortReconstruction.h>
//#include "../trackingdiagnostics/KshortReconstruction.cc"

Acts::BoundTrackParameters makeTrackParams(SvtxTrack* track);
Acts::Vector3 getVertex(SvtxTrack* track);
ActsGeometry* _tGeometry = nullptr;

SvtxVertexMap* m_vertexMap = nullptr;
Acts::Result<BoundTrackParam> propagateTrack(const Acts::BoundTrackParameters& params, const SurfacePtr& targetSurf);


int NTry = 500;
float StepMove = 0.01;
float PerturbationX = 0;		
float PerturbationY = 0;		

int NTrial = 200;
int NPert = 100;
float PertRange = 0.30;

int Initial = 20;

float AngleRatio = 1; 

std::ofstream PrintInfo("PrintInfo.log");

//float ComparePx[20] = {1.5082383,1.5903980,0.5210004,-1.055958,-2.288107,2.8571836,-1.883892,-0.751257,-2.375273,1.6151963,2.4277606,1.2849296,1.9204852,-0.805020,-0.965543,0.3492338,-2.570624,-0.109778,-1.249456,0.7592150};
//float KFParTestX[20] = {24.814907,0.6694057,18.505744,-7.542235,-11.63196,17.549827,-28.93202,-4.482494,-15.56186,4.1796078,24.818105,0.6692505,18.506456,-7.542699,-11.63092,17.548925,-28.92981,-4.485909,-15.56326,4.1799278};
//float KFParTestY[20] = {-20.06811,0.8657092,37.149852,-15.23564,-12.41074,18.29916,2.6157825,-20.96750,12.294669,-6.643179,-20.06574,0.8657371,37.150939,-15.23553,-12.41175,18.297948,2.6091358,-20.96584,12.299227,-6.643013};
//float KFParTestZ[20] = {16.059873,-0.931764,-22.33364,-5.245448,5.9941983,-1.548370,-23.73804,-22.18427,-23.04204,0.2520581,16.053299,-0.931898,-22.33109, -5.244925, 5.9925980, -1.548754,-23.74473,-22.18704,-23.04225,0.2467427};


float SelectTrackX =  0.8954519;

/// KFParticle constructor
KFParticle_Tools::KFParticle_Tools()
	: m_has_intermediates(false)
	, m_min_mass(0)
	, m_max_mass(0)
	, m_min_decayTime(-1 * FLT_MAX)
	, m_max_decayTime(FLT_MAX)
	, m_min_decayLength(-1 * FLT_MAX)
	, m_max_decayLength(FLT_MAX)
	, m_track_pt(0.25)
	, m_track_ptchi2(FLT_MAX)
	, m_track_ip(-1.)
	, m_track_ipchi2(10.)
	, m_track_chi2ndof(4.)
	, m_nMVTXHits(3)
	, m_nTPCHits(20)
	, m_comb_DCA(0.05)
	, m_vertex_chi2ndof(15.)
	, m_fdchi2(50.)
	, m_dira_min(0.95)
	, m_dira_max(1.01)
	, m_mother_pt(0.)
	, m_mother_ipchi2(FLT_MAX)
	, m_get_charge_conjugate(true)
	  , m_extrapolateTracksToSV(true)
	  , m_vtx_map_node_name("SvtxVertexMap")
	, m_trk_map_node_name("SvtxTrackMap")
	, m_dst_vertexmap()
	, m_dst_trackmap()
	, m_dst_vertex()
	  , m_dst_track()
{

	//std::cout << "STARTing from 0 to NTry ZZ Bro OK  - ACTS + + - Vec Prop ACTSHitsX - KFPropX = TestKF.Get_X();  - ALL out: No extct  Perturbation: "  << "   Trials - " << NTrial << "   NTry = " << NTry << "   No Radius Added" << "   Back to cm units bro" << std::endl;
	std::cout << "Check Covariance Matrices Requested by Cameron - ACTS Updated - Use Other nOw 76316" << std::endl;
}

KFParticle KFParticle_Tools::makeVertex(PHCompositeNode * /*topNode*/)
{
	float f_vertexParameters[6] = {m_dst_vertex->get_x(),
		m_dst_vertex->get_y(),
		m_dst_vertex->get_z(), 0, 0, 0};


	//	std::cout << "WHat the HECK is - makeVertex: m_dst_vertex->get_x() =  " << m_dst_vertex->get_x() << "  m_dst_vertex->get_y() =  " << m_dst_vertex->get_y() << "  m_dst_vertex->get_z() =  " << m_dst_vertex->get_z() << std::endl;


	float f_vertexCovariance[21];
	unsigned int iterate = 0;
	for (unsigned int i = 0; i < 3; ++i)
		for (unsigned int j = 0; j <= i; ++j)
		{
			f_vertexCovariance[iterate] = m_dst_vertex->get_error(i, j);
			++iterate;
		}

	KFParticle kfp_vertex;
	kfp_vertex.Create(f_vertexParameters, f_vertexCovariance, 0, -1);
	kfp_vertex.NDF() = m_dst_vertex->get_ndof();
	kfp_vertex.Chi2() = m_dst_vertex->get_chisq();

	return kfp_vertex;
}

std::vector<KFParticle> KFParticle_Tools::makeAllPrimaryVertices(PHCompositeNode *topNode, const std::string &vertexMapName)
{
	std::string vtxMN;
	if (vertexMapName.empty()){
		vtxMN = m_vtx_map_node_name;
		//std::cout << "VTX Map is Empty Why" << std::endl;
	}
	else{
		vtxMN = vertexMapName;
		//std::cout << "VTX Map is NOT Empty - " << "   It is " << vtxMN << std::endl;

	}

	std::vector<KFParticle> primaryVertices;
	m_dst_vertexmap = findNode::getClass<SvtxVertexMap>(topNode, vtxMN);
	unsigned int vertexID = 0;

	//	std::cout << "m_dst_vertexmap.size() = " <<  m_dst_vertexmap->size() << std::endl;

	for (SvtxVertexMap::ConstIter iter = m_dst_vertexmap->begin(); iter != m_dst_vertexmap->end(); ++iter)
	{
		m_dst_vertex = iter->second;

		primaryVertices.push_back(makeVertex(topNode));
		primaryVertices[vertexID].SetId(iter->first);

		//		std::cout << "primaryVertices[vertexID].X() = " << primaryVertices[vertexID].X() << "    primaryVertices[vertexID].Y() = " << primaryVertices[vertexID].Y() << "   primaryVertices[vertexID].Z() = " << primaryVertices[vertexID].Z() << std::endl; 
		++vertexID;
	}

	return primaryVertices;
}

KFParticle KFParticle_Tools::makeParticle(PHCompositeNode * /*topNode*/)  ///Return a KFPTrack from track vector and covariance matrix. No mass or vertex constraints
{
	KFParticle kfp_particle;

	float f_trackParameters[6] = {m_dst_track->get_x(),
		m_dst_track->get_y(),
		m_dst_track->get_z(),
		m_dst_track->get_px(),
		m_dst_track->get_py(),
		m_dst_track->get_pz()};

	float f_trackCovariance[21];
	unsigned int iterate = 0;
	for (unsigned int i = 0; i < 6; ++i)
		for (unsigned int j = 0; j <= i; ++j)
		{
			f_trackCovariance[iterate] = m_dst_track->get_error(i, j);
			++iterate;
		}

	kfp_particle.Create(f_trackParameters, f_trackCovariance, (Int_t) m_dst_track->get_charge(), -1);
	kfp_particle.NDF() = m_dst_track->get_ndf();
	kfp_particle.Chi2() = m_dst_track->get_chisq();
	kfp_particle.SetId(m_dst_track->get_id());


	return kfp_particle;
}

std::vector<KFParticle> KFParticle_Tools::makeAllDaughterParticles(PHCompositeNode *topNode)
{

	//std::cout << "FUCKING Make KFParticles"  << std::endl;

	std::vector<KFParticle> daughterParticles;
	m_dst_trackmap = findNode::getClass<SvtxTrackMap>(topNode, m_trk_map_node_name.c_str());

	unsigned int trackID = 0;


	for (SvtxTrackMap::Iter iter = m_dst_trackmap->begin(); iter != m_dst_trackmap->end(); ++iter)
	{
		m_dst_track = iter->second;

		//First check if we have the required number of MVTX and TPC hits
		TrackSeed *tpcseed = m_dst_track->get_tpc_seed();
		TrackSeed *silseed = m_dst_track->get_silicon_seed();
		int MVTX_hits = 0;
		int TPC_hits = 0;

		//std::cout << "FUCKING Make KFParticles: m_dst_trackmap->get_x() =  " << m_dst_track->get_x() << "  m_dst_trackmap->get_y() =  " << m_dst_track->get_y() << "  m_dst_trackmap->get_z() =  " << m_dst_track->get_z() << std::endl;

		if (silseed)
		{
			for (auto cluster_iter = silseed->begin_cluster_keys();
					cluster_iter != silseed->end_cluster_keys(); ++cluster_iter)
			{
				const auto &cluster_key = *cluster_iter;
				const auto trackerID = TrkrDefs::getTrkrId(cluster_key);

				if (trackerID == TrkrDefs::mvtxId) ++MVTX_hits;
			}
			if (MVTX_hits < m_nMVTXHits) continue;
		}
		if (tpcseed)
		{
			for (auto cluster_iter = tpcseed->begin_cluster_keys(); cluster_iter != tpcseed->end_cluster_keys(); ++cluster_iter)
			{
				const auto &cluster_key = *cluster_iter;
				const auto trackerID = TrkrDefs::getTrkrId(cluster_key);

				if (trackerID == TrkrDefs::tpcId) ++TPC_hits;
			}
			if (TPC_hits < m_nTPCHits) continue;
		}

		daughterParticles.push_back(makeParticle(topNode));  ///Turn all dst tracks in KFP tracks
		daughterParticles[trackID].SetId(iter->first);

		//std::cout << "KFParticle: X =  " << daughterParticles[trackID].GetX() << "  Y =  " << daughterParticles[trackID].GetY() << "  Z =  " << daughterParticles[trackID].GetZ() << "    Px = "  << daughterParticles[trackID].GetPx() << "  Py =  " << daughterParticles[trackID].GetPy() << "  Pz =  " << daughterParticles[trackID].GetPz()  << std::endl;

		++trackID;
	}

	//std::cout << "------------------------------ Now I Do ZZ Shit ----------------------------------" << std::endl; 


	m_dst_trackmap = findNode::getClass<SvtxTrackMap>(topNode, m_trk_map_node_name.c_str());

	m_vertexMap = findNode::getClass<SvtxVertexMap>(topNode, m_vtx_map_node_name);


	_tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");

			
	float FinalMomX = 1.76316;
	
	for (auto& [key, track] : *m_dst_trackmap)
	{


		std::cout << "Track Here - " << track->get_x() << "   px = " << track->get_px()  <<  std::endl;

		if(abs(track->get_x()-SelectTrackX) < 0.001){
	
			std::cout << "------------------------------------------------ Do ACTS First BRO ------------------------------------------------" << std::endl;

		for (auto iter = track->begin_states(); iter != track->end_states(); ++iter)
		{



			std::cout  << "x =  "  << iter->second->get_x() <<  "  y =  "  << iter->second->get_y()  <<  "  z =  "  << iter->second->get_z()  << "   px =  "  << iter->second->get_px() <<  "  py =  "  << iter->second->get_py()  <<  "  pz =  "  << iter->second->get_pz() << std::endl;

			float TrackStatX = iter->second->get_x();
			float TrackStatY = iter->second->get_y();
			float TrackStatZ = iter->second->get_z();
			float TrackStatPx = iter->second->get_px();
			float TrackStatPy = iter->second->get_py();
			float TrackStatPz = iter->second->get_pz();


			if(abs(TrackStatPx - FinalMomX) < 0.001){



				PrintInfo << "------------------------------------------------ Record Info ACTS Propagating to KF Mother - After Propagation ------------------------------------------------" << std::endl;
				PrintInfo << "ACTS: " << "  After Prop: X =  "  <<   TrackStatX << "   Y = "  << TrackStatY  << "   Z = "   << TrackStatZ     << "   Px =  "  <<   TrackStatPx << "   Py = "  << TrackStatPy  << "   Pz = "   <<  TrackStatPz  << std::endl;

				PrintInfo << "---------------------------------------------------------------------- ATCS Final Track Covariance Matrix BEGIN ----------------------------------------------------------------------" << std::endl;

				for(int r = 0; r < 6; r++){

					for(int s = 0; s < 6; s++){

						PrintInfo << "ACTS Matched: r =  " << r << "   s = " << s << "     iter->second->get_error(r,s) " <<  iter->second->get_error(r,s) << std::endl;

					}
				}

				PrintInfo << "---------------------------------------------------------------------- ATCS Final Track Covariance Matrix END ----------------------------------------------------------------------" << std::endl;

			}



		}

		PrintInfo << "------------------------------------------------ DONE ACTS Record ------------------------------------------------" << std::endl;

	}





	}



	for (SvtxTrackMap::Iter iter = m_dst_trackmap->begin(); iter != m_dst_trackmap->end(); ++iter)
	{
		m_dst_track = iter->second;


		float trackpx = m_dst_track->get_px();
		float trackpy = m_dst_track->get_py();
		float trackpz = m_dst_track->get_pz();


		float trackx = m_dst_track->get_x();
		float tracky = m_dst_track->get_y();
		float trackz = m_dst_track->get_z();



		//KFParticle Mother Vertex Position//
		float ACTSProX =  18.5062;
		float ACTSProY =  37.1506;
		float ACTSProZ =  -22.332;



		auto perigee = Acts::Surface::makeShared<Acts::PerigeeSurface>(Acts::Vector3(ACTSProX * Acts::UnitConstants::cm ,ACTSProY * Acts::UnitConstants::cm ,ACTSProZ * Acts::UnitConstants::cm));
		auto params = makeTrackParams(m_dst_track);		
		auto result = propagateTrack(params, perigee);

		if (result.ok())
		{
			auto projectionPos = result.value().position(_tGeometry->geometry().getGeoContext());
			const auto momentum = result.value().momentum();

			float pos0 = projectionPos.x() / Acts::UnitConstants::cm;
			float pos1 = projectionPos.y() / Acts::UnitConstants::cm;
			float pos2 = projectionPos.z() / Acts::UnitConstants::cm;


			float mom0 = momentum.x();
			float mom1 = momentum.y();
			float mom2 = momentum.z();

			if(abs(trackx - SelectTrackX) < 0.0001 ){	




				PrintInfo << "------------------------------------------------ Record Info ACTS Propagating to KF Mother - Before Propagation ------------------------------------------------" << std::endl;

				PrintInfo << "ACTS: " << "  Before Prop: X =  "  <<   trackx << "   Y = "  << tracky  << "   Z = "   << trackz     << "   Px =  "  <<   trackpx << "   Py = "  << trackpy  << "   Pz = "   <<  trackpz  << std::endl;
				PrintInfo << "ACTS Propagation to Vertex: X =  "  <<   ACTSProX << "   Y = "  << ACTSProY  << "   Z = "   << ACTSProZ    << std::endl;
				PrintInfo << "ACTS Actually Gives Us: " << "  After Prop: X =  "  <<   pos0 << "   Y = "  << pos1  << "   Z = "   << pos2     << "   Px =  "  <<   mom0 << "   Py = "  << mom1  << "   Pz = "   <<  mom2  << std::endl;

				PrintInfo << "---------------------------------------------------------------------- Daughter Track Covariance Matrix BEGIN ----------------------------------------------------------------------" << std::endl;

				for(int r = 0; r < 6; r++){
					for(int s = 0; s < 6; s++){

						PrintInfo << "r =  " << r << "   s = " << s  <<  "     m_dst_track->get_error(r,s) = "  <<  m_dst_track->get_error(r,s) << std::endl;		

					}
				}

				PrintInfo << "---------------------------------------------------------------------- Daughter Track Covariance Matrix END ----------------------------------------------------------------------" << std::endl;

			}

		}

	}


	return daughterParticles;
}

int KFParticle_Tools::getTracksFromVertex(PHCompositeNode *topNode, KFParticle vertex, const std::string &vertexMapName)
{
	std::string vtxMN;
	if (vertexMapName.empty())
		vtxMN = m_vtx_map_node_name;
	else
		vtxMN = vertexMapName;

	SvtxVertex *associatedVertex = NULL;
	m_dst_vertexmap = findNode::getClass<SvtxVertexMap>(topNode, vtxMN);

	associatedVertex = m_dst_vertexmap->find(vertex.Id())->second;

	return associatedVertex->size_tracks();
}

/*const*/ bool KFParticle_Tools::isGoodTrack(KFParticle particle, const std::vector<KFParticle> &primaryVertices)
{
	bool goodTrack = false;

	float min_ip = 0;
	float min_ipchi2 = 0;

	float pt = particle.GetPt();
	float pterr = particle.GetErrPt();
	float ptchi2 = pow(pterr / pt, 2);
	float trackchi2ndof = particle.GetChi2() / particle.GetNDF();
	calcMinIP(particle, primaryVertices, min_ip, min_ipchi2);

	if (pt >= m_track_pt 
			&& ptchi2 <= m_track_ptchi2 
			&& min_ip >= m_track_ip 
			&& min_ipchi2 >= m_track_ipchi2 
			&& trackchi2ndof <= m_track_chi2ndof)
		goodTrack = true;

	return goodTrack;
}

int KFParticle_Tools::calcMinIP(KFParticle track, std::vector<KFParticle> PVs,
		float &minimumIP, float &minimumIPchi2)
{
	std::vector<float> ip, ipchi2;

	for (unsigned int i_verts = 0; i_verts < PVs.size(); ++i_verts)
	{
		ip.push_back(track.GetDistanceFromVertex(PVs[i_verts]));
		ipchi2.push_back(track.GetDeviationFromVertex(PVs[i_verts]));
	}

	auto minmax_ip = minmax_element(ip.begin(), ip.end());  //Order the IP from small to large
	minimumIP = *minmax_ip.first;
	auto minmax_ipchi2 = minmax_element(ipchi2.begin(), ipchi2.end());  //Order the IP chi2 from small to large
	minimumIPchi2 = *minmax_ipchi2.first;

	return 0;
}

std::vector<int> KFParticle_Tools::findAllGoodTracks(std::vector<KFParticle> daughterParticles, const std::vector<KFParticle> primaryVertices)
{
	std::vector<int> goodTrackIndex;

	for (unsigned int i_parts = 0; i_parts < daughterParticles.size(); ++i_parts)
	{
		if (isGoodTrack(daughterParticles[i_parts], primaryVertices)) goodTrackIndex.push_back(i_parts);
	}

	removeDuplicates(goodTrackIndex);

	return goodTrackIndex;
}

std::vector<std::vector<int>> KFParticle_Tools::findTwoProngs(std::vector<KFParticle> daughterParticles, std::vector<int> goodTrackIndex, int nTracks)
{
	std::vector<std::vector<int>> goodTracksThatMeet;

	for (std::vector<int>::iterator i_it = goodTrackIndex.begin(); i_it != goodTrackIndex.end(); ++i_it)
	{
		for (std::vector<int>::iterator j_it = goodTrackIndex.begin(); j_it != goodTrackIndex.end(); ++j_it)
		{
			if (i_it < j_it)
			{
				if (daughterParticles[*i_it].GetDistanceFromParticle(daughterParticles[*j_it]) <= m_comb_DCA)
				{
					KFVertex twoParticleVertex;
					twoParticleVertex += daughterParticles[*i_it];
					twoParticleVertex += daughterParticles[*j_it];
					float vertexchi2ndof = twoParticleVertex.GetChi2() / twoParticleVertex.GetNDF();
					std::vector<int> combination = {*i_it, *j_it};

					if (nTracks == 2 && vertexchi2ndof <= m_vertex_chi2ndof)
					{
						goodTracksThatMeet.push_back(combination);
					}
					else if (nTracks == 2 && vertexchi2ndof > m_vertex_chi2ndof)
					{
						continue;
					}
					else
					{
						goodTracksThatMeet.push_back(combination);
					}
				}
			}
		}
	}

	return goodTracksThatMeet;
}

std::vector<std::vector<int>> KFParticle_Tools::findNProngs(std::vector<KFParticle> daughterParticles,
		std::vector<int> goodTrackIndex,
		std::vector<std::vector<int>> goodTracksThatMeet,
		int nRequiredTracks, unsigned int nProngs)
{
	unsigned int nGoodProngs = goodTracksThatMeet.size();

	for (std::vector<int>::iterator i_it = goodTrackIndex.begin(); i_it != goodTrackIndex.end(); ++i_it)
	{
		for (unsigned int i_prongs = 0; i_prongs < nGoodProngs; ++i_prongs)
		{
			bool trackNotUsedAlready = true;
			for (unsigned int i_trackCheck = 0; i_trackCheck < nProngs - 1; ++i_trackCheck)
				if (*i_it == goodTracksThatMeet[i_prongs][i_trackCheck]) trackNotUsedAlready = false;
			if (trackNotUsedAlready)
			{
				bool dcaMet = 1;
				for (unsigned int i = 0; i < nProngs - 1; ++i)
				{
					if (daughterParticles[*i_it].GetDistanceFromParticle(daughterParticles[goodTracksThatMeet[i_prongs][i]]) > m_comb_DCA)
					{
						dcaMet = 0;
					}
				}

				if (dcaMet)
				{
					KFVertex particleVertex;
					particleVertex += daughterParticles[*i_it];
					std::vector<int> combination;
					combination.push_back(*i_it);
					for (unsigned int i = 0; i < nProngs - 1; ++i)
					{
						particleVertex += daughterParticles[goodTracksThatMeet[i_prongs][i]];
						combination.push_back(goodTracksThatMeet[i_prongs][i]);
					}
					float vertexchi2ndof = particleVertex.GetChi2() / particleVertex.GetNDF();

					if ((unsigned int) nRequiredTracks == nProngs && vertexchi2ndof <= m_vertex_chi2ndof)
					{
						goodTracksThatMeet.push_back(combination);
					}
					else if ((unsigned int) nRequiredTracks == nProngs && vertexchi2ndof > m_vertex_chi2ndof)
					{
						continue;
					}
					else
					{
						goodTracksThatMeet.push_back(combination);
					}
				}
			}
		}
	}

	goodTracksThatMeet.erase(goodTracksThatMeet.begin(), goodTracksThatMeet.begin() + nGoodProngs);
	for (unsigned int i = 0; i < goodTracksThatMeet.size(); ++i) sort(goodTracksThatMeet[i].begin(), goodTracksThatMeet[i].end());
	removeDuplicates(goodTracksThatMeet);

	return goodTracksThatMeet;
}

std::vector<std::vector<int>> KFParticle_Tools::appendTracksToIntermediates(KFParticle intermediateResonances[], std::vector<KFParticle> daughterParticles, std::vector<int> goodTrackIndex, int num_remaining_tracks)
{
	std::vector<std::vector<int>> goodTracksThatMeet, goodTracksThatMeetIntermediates;  //, vectorOfGoodTracks;
	if (num_remaining_tracks == 1)
	{
		for (std::vector<int>::iterator i_it = goodTrackIndex.begin(); i_it != goodTrackIndex.end(); ++i_it)
		{
			std::vector<KFParticle> v_intermediateResonances(intermediateResonances, intermediateResonances + m_num_intermediate_states);
			std::vector<std::vector<int>> dummyTrackList;
			std::vector<int> dummyTrackID;  //I already have the track ids stored in goodTracksThatMeet[i]
			v_intermediateResonances.insert(end(v_intermediateResonances), daughterParticles[*i_it]);
			for (unsigned int k = 0; k < v_intermediateResonances.size(); ++k)
			{
				dummyTrackID.push_back(k);
			}
			dummyTrackList = findTwoProngs(v_intermediateResonances, dummyTrackID, (int) v_intermediateResonances.size());
			if (v_intermediateResonances.size() > 2)
			{
				for (unsigned int p = 3; p <= v_intermediateResonances.size(); ++p) dummyTrackList = findNProngs(v_intermediateResonances,
						dummyTrackID, dummyTrackList,
						(int) v_intermediateResonances.size(), (int) p);
			}

			if (dummyTrackList.size() != 0)
			{
				std::vector<int> goodTrack{*i_it};
				goodTracksThatMeetIntermediates.push_back(goodTrack);
			}
		}
	}
	else
	{
		goodTracksThatMeet = findTwoProngs(daughterParticles, goodTrackIndex, num_remaining_tracks);

		for (int p = 3; p <= num_remaining_tracks; ++p)
		{
			goodTracksThatMeet = findNProngs(daughterParticles, goodTrackIndex, goodTracksThatMeet, num_remaining_tracks, p);
		}

		for (unsigned int i = 0; i < goodTracksThatMeet.size(); ++i)
		{
			std::vector<KFParticle> v_intermediateResonances(intermediateResonances, intermediateResonances + m_num_intermediate_states);
			std::vector<std::vector<int>> dummyTrackList;
			std::vector<int> dummyTrackID;  //I already have the track ids stored in goodTracksThatMeet[i]
			for (unsigned int j = 0; j < goodTracksThatMeet[i].size(); ++j)
			{
				v_intermediateResonances.push_back(daughterParticles[goodTracksThatMeet[i][j]]);
			}
			for (unsigned int k = 0; k < v_intermediateResonances.size(); ++k)
			{
				dummyTrackID.push_back(k);
			}
			dummyTrackList = findTwoProngs(v_intermediateResonances, dummyTrackID, (int) v_intermediateResonances.size());
			for (unsigned int p = 3; p <= v_intermediateResonances.size(); ++p)
			{
				dummyTrackList = findNProngs(v_intermediateResonances, dummyTrackID, dummyTrackList, (int) v_intermediateResonances.size(), (int) p);
			}

			if (dummyTrackList.size() != 0) goodTracksThatMeetIntermediates.push_back(goodTracksThatMeet[i]);
		}
	}

	return goodTracksThatMeetIntermediates;
}

float KFParticle_Tools::eventDIRA(KFParticle particle, KFParticle vertex)
{
	TMatrixD flightVector(3, 1);
	TMatrixD momVector(3, 1);
	flightVector(0, 0) = particle.GetX() - vertex.GetX();
	flightVector(1, 0) = particle.GetY() - vertex.GetY();
	flightVector(2, 0) = particle.GetZ() - vertex.GetZ();

	momVector(0, 0) = particle.GetPx();
	momVector(1, 0) = particle.GetPy();
	momVector(2, 0) = particle.GetPz();

	TMatrixD momDotFD(1, 1);  //Calculate momentum dot flight distance
	momDotFD = TMatrixD(momVector, TMatrixD::kTransposeMult, flightVector);
	float f_momDotFD = momDotFD(0, 0);

	TMatrixD sizeOfMom(1, 1);  //Calculates the size of the momentum vector
	sizeOfMom = TMatrixD(momVector, TMatrixD::kTransposeMult, momVector);
	float f_sizeOfMom = sqrt(sizeOfMom(0, 0));

	TMatrixD sizeOfFD(1, 1);  //Calculates the size of the flight distance vector
	sizeOfFD = TMatrixD(flightVector, TMatrixD::kTransposeMult, flightVector);
	float f_sizeOfFD = sqrt(sizeOfFD(0, 0));

	return f_momDotFD / (f_sizeOfMom * f_sizeOfFD);
}

float KFParticle_Tools::flightDistanceChi2(KFParticle particle, KFParticle vertex)
{
	TMatrixD flightVector(3, 1);
	TMatrixD flightDistanceCovariance(3, 3);

	KFParticle kfp_vertex(vertex);

	flightVector(0, 0) = particle.GetX() - kfp_vertex.GetX();
	flightVector(1, 0) = particle.GetY() - kfp_vertex.GetY();
	flightVector(2, 0) = particle.GetZ() - kfp_vertex.GetZ();

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			flightDistanceCovariance(i, j) = particle.GetCovariance(i, j) + kfp_vertex.GetCovariance(i, j);
		}
	}

	TMatrixD anInverseMatrix(3, 3);
	anInverseMatrix = flightDistanceCovariance.Invert();
	TMatrixD m_chi2Value(1, 1);
	m_chi2Value = TMatrixD(flightVector, TMatrixD::kTransposeMult, anInverseMatrix * flightVector);
	return m_chi2Value(0, 0);
}

std::tuple<KFParticle, bool> KFParticle_Tools::buildMother(KFParticle vDaughters[], std::string daughterOrder[],
		bool isIntermediate, int intermediateNumber, int nTracks,
		bool constrainMass, float required_vertexID)
{
	KFParticle mother;
	KFParticle inputTracks[nTracks];

	mother.SetConstructMethod(2);

	bool daughterMassCheck = true;
	float unique_vertexID = 0;


	//Figure out if the decay has reco. tracks mixed with resonances
	int num_tracks_used_by_intermediates = 0;
	for (int i = 0; i < m_num_intermediate_states; ++i) num_tracks_used_by_intermediates += m_num_tracks_from_intermediate[i];
	int num_remaining_tracks = m_num_tracks - num_tracks_used_by_intermediates;

	//	KFParticle TestKF;




	for (int i = 0; i < nTracks; ++i)
	{
		float daughterMass = 0;
		daughterMass = constrainMass ? getParticleMass(daughterOrder[i].c_str()) : vDaughters[i].GetMass();
		if ((num_remaining_tracks > 0 && i >= m_num_intermediate_states) || isIntermediate)
		{
			daughterMass = getParticleMass(daughterOrder[i].c_str());
		}

		//	std::cout << "Initial Daughters: " << "   X = " << vDaughters[i].GetX() << "   Y = " << vDaughters[i].GetY() << "   Z = " << vDaughters[i].GetZ() << "   Px = " <<  vDaughters[i].GetPx()  << std::endl;

		inputTracks[i].Create(vDaughters[i].Parameters(),
				vDaughters[i].CovarianceMatrix(),
				(Int_t) vDaughters[i].GetQ(),
				daughterMass);
		mother.AddDaughter(inputTracks[i]);
		unique_vertexID += vDaughters[i].GetQ() * getParticleMass(daughterOrder[i].c_str());



	}

	if (isIntermediate) mother.SetPDG(getParticleID(m_intermediate_name[intermediateNumber].c_str()));
	if (!isIntermediate && !m_mother_name_Tools.empty()) mother.SetPDG(getParticleID(m_mother_name_Tools));

	bool chargeCheck;
	if (m_get_charge_conjugate)
		chargeCheck = std::abs(unique_vertexID) == std::abs(required_vertexID) ? 1 : 0;
	else
		chargeCheck = unique_vertexID == required_vertexID ? 1 : 0;


	bool printed = false;

	for(int j =0; j< nTracks;++j){



		float trackx = inputTracks[j].GetX();
		float tracky = inputTracks[j].GetY();
		float trackz = inputTracks[j].GetZ();
		float trackpx = inputTracks[j].GetPx();
		float trackpy = inputTracks[j].GetPy();
		float trackpz = inputTracks[j].GetPz();


		if(abs(trackx - SelectTrackX) < 0.0001 && printed == false ){	

			PrintInfo << "---------------------------------------------------------------------- KFParticle Information Before Propagation ----------------------------------------------------------------------" << std::endl;


			PrintInfo << "KFParticle: " << "  Before Prop: X =  "  <<   trackx << "   Y = "  << tracky  << "   Z = "   << trackz     << "   Px =  "  <<   trackpx << "   Py = "  << trackpy  << "   Pz = "   <<  trackpz  << std::endl;

			PrintInfo << "---------------------------------------------------------------------- Daughter Track Covariance Matrix BEGIN ----------------------------------------------------------------------" << std::endl;

			for(int r = 0; r < 6; r++){
				for(int s = 0; s < 6; s++){

					PrintInfo << "r =  " << r << "   s = " << s  <<  "     inputTracks[j].GetCovariance(r,s) = "  << inputTracks[j].GetCovariance(r,s) << std::endl;		

				}
			}

			PrintInfo << "---------------------------------------------------------------------- Daughter Track Covariance Matrix END ----------------------------------------------------------------------" << std::endl;



			PrintInfo << "Mother Seconary Vertex Position: " << "  X =  "  <<   mother.GetX() << "   Y = "  << mother.GetY()  << "   Z = "   << mother.GetZ()    << std::endl;		
			PrintInfo << "---------------------------------------------------------------------- Mother Seconary Vertex Covariance Matrix BEGIN ----------------------------------------------------------------------" << std::endl;

			for(int r = 0; r < 6; r++){
				for(int s = 0; s < 6; s++){

					PrintInfo << "r =  " << r << "   s = " << s  <<  "     mother.GetCovariance(r,s) = "  << mother.GetCovariance(r,s) << std::endl;		

				}
			}

			PrintInfo << "---------------------------------------------------------------------- Mother Seconary Vertex Covariance Matrix END ----------------------------------------------------------------------" << std::endl;

			printed = true;
		}





		inputTracks[j].SetProductionVertex(mother);

		if(abs(trackx - SelectTrackX) < 0.0001 ){	

			PrintInfo << "---------------------------------------------------------------------- KFParticle Information After Propagation ----------------------------------------------------------------------" << std::endl;

			PrintInfo << "KFParticle: " << "  After Prop: X =  "  <<    inputTracks[j].GetX() << "   Y = "  <<  inputTracks[j].GetY()   << "   Z = "   <<  inputTracks[j].GetZ()      << "   Px =  "  <<    inputTracks[j].GetPx()  << "   Py = "  <<  inputTracks[j].GetPy()   << "   Pz = "   <<  trackpz  << std::endl;

			PrintInfo << "---------------------------------------------------------------------- Daughter Track Covariance Matrix BEGIN ----------------------------------------------------------------------" << std::endl;

			for(int r = 0; r < 6; r++){
				for(int s = 0; s < 6; s++){

					PrintInfo << "r =  " << r << "   s = " << s  <<  "     inputTracks[j].GetCovariance(r,s) = "  << inputTracks[j].GetCovariance(r,s) << std::endl;		

				}
			}

			PrintInfo << "---------------------------------------------------------------------- Daughter Track Covariance Matrix END ----------------------------------------------------------------------" << std::endl;

		}


		if (!m_allowZeroMassTracks)
		{
			if (inputTracks[j].GetMass() == 0) daughterMassCheck = false;
		}
	}

	float calculated_mass, calculated_mass_err;
	mother.GetMass(calculated_mass, calculated_mass_err);
	float calculated_pt = mother.GetPt();

	float min_mass = isIntermediate ? m_intermediate_mass_range[intermediateNumber].first : m_min_mass;
	float max_mass = isIntermediate ? m_intermediate_mass_range[intermediateNumber].second : m_max_mass;
	float min_pt = isIntermediate ? m_intermediate_min_pt[intermediateNumber] : m_mother_pt;

	bool goodCandidate = false;
	if (calculated_mass >= min_mass && calculated_mass <= max_mass &&
			calculated_pt >= min_pt && daughterMassCheck && chargeCheck)
		goodCandidate = true;

	// Check the requirements of an intermediate states against this mother and re-do goodCandidate
	if (goodCandidate && m_has_intermediates && !isIntermediate)  //The decay has intermediate states and we are now looking at the mother
	{
		for (int k = 0; k < m_num_intermediate_states; ++k)
		{
			float intermediate_DIRA = eventDIRA(vDaughters[k], mother);
			float intermediate_FDchi2 = flightDistanceChi2(vDaughters[k], mother);
			if (intermediate_DIRA < m_intermediate_min_dira[k] ||
					intermediate_FDchi2 < m_intermediate_min_fdchi2[k])
				goodCandidate = false;
		}
	}

	/*
	   foutDebug->cd();
	   ACTSTree->Write();
	   KFTree->Write();
	   foutDebug->Close();
	   */


	return std::make_tuple(mother, goodCandidate);
}

void KFParticle_Tools::constrainToVertex(KFParticle &particle, bool &goodCandidate, KFParticle &vertex)
{
	KFParticle particleCopy = particle;
	particleCopy.SetProductionVertex(vertex);
	particleCopy.TransportToDecayVertex();

	float calculated_decayTime;
	float calculated_decayTimeErr;
	float calculated_decayLength;
	float calculated_decayLengthErr;

	particleCopy.GetLifeTime(calculated_decayTime, calculated_decayTimeErr);
	particleCopy.GetDecayLength(calculated_decayLength, calculated_decayLengthErr);

	float calculated_fdchi2 = flightDistanceChi2(particle, vertex);
	float calculated_dira = eventDIRA(particle, vertex);
	float calculated_ipchi2 = particle.GetDeviationFromVertex(vertex);

	goodCandidate = false;

	const float speed = 2.99792458e-1;
	calculated_decayTime /= speed;

	if (calculated_fdchi2 >= m_fdchi2 && calculated_ipchi2 <= m_mother_ipchi2 && isInRange(m_dira_min, calculated_dira, m_dira_max) && isInRange(m_min_decayTime, calculated_decayTime, m_max_decayTime) && isInRange(m_min_decayLength, calculated_decayLength, m_max_decayLength))
		goodCandidate = true;
}

std::tuple<KFParticle, bool> KFParticle_Tools::getCombination(KFParticle vDaughters[], std::string daughterOrder[], KFParticle vertex, bool constrain_to_vertex, bool isIntermediate, int intermediateNumber, int nTracks, bool constrainMass, float required_vertexID)
{
	KFParticle candidate;
	bool isGoodCandidate;

	std::tie(candidate, isGoodCandidate) = buildMother(vDaughters, daughterOrder, isIntermediate, intermediateNumber, nTracks, constrainMass, required_vertexID);

	if (constrain_to_vertex && isGoodCandidate && !isIntermediate) constrainToVertex(candidate, isGoodCandidate, vertex);

	return std::make_tuple(candidate, isGoodCandidate);
}

std::vector<std::vector<std::string>> KFParticle_Tools::findUniqueDaughterCombinations(int start, int end)
{
	std::vector<int> vect_permutations;
	std::vector<std::vector<std::string>> uniqueCombinations;
	std::map<int, std::string> daughterMap;
	for (int i = start; i < end; i++)
	{
		daughterMap.insert(std::pair<int, std::string>(i, m_daughter_name[i]));
		vect_permutations.push_back(i);
	}
	int *permutations = &vect_permutations[0];

	do
	{
		std::vector<std::string> combination;
		for (int i = 0; i < (end - start); i++) combination.push_back(daughterMap.find(permutations[i])->second);
		uniqueCombinations.push_back(combination);
	} while (std::next_permutation(permutations, permutations + vect_permutations.size()));

	removeDuplicates(uniqueCombinations);

	return uniqueCombinations;
}

double KFParticle_Tools::calculateEllipsoidRadius(int posOrNeg, double sigma_ii, double sigma_jj, double sigma_ij)
{  //Note - Only works for a 2D ellipsoid OR rotated nD ellipsoid to avoid projections
	if (std::abs(posOrNeg) != 1)
	{
		std::cout << "You have set posOrNeg to " << posOrNeg << ". This value must be  +/- 1! Skipping" << std::endl;
		return 0;
	}

	double r_ij = sqrt((sigma_ii + sigma_jj) / 2 + posOrNeg * (sqrt(pow((sigma_ii - sigma_jj) / 2, 2) + pow(sigma_ij, 2))));

	return r_ij;
}

float KFParticle_Tools::calculateEllipsoidVolume(KFParticle particle)
{
	TMatrixD cov_matrix(3, 3);

	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			cov_matrix(i, j) = particle.GetCovariance(i, j);

	float volume;
	if (cov_matrix(0, 0) * cov_matrix(1, 1) * cov_matrix(2, 2) == 0)
		volume = 0;
	else
		volume = (4. / 3.) * M_PI * sqrt((std::abs(cov_matrix.Determinant())));  //The covariance matrix is error-squared

	return volume;
}

float KFParticle_Tools::calculateJT(KFParticle mother, KFParticle daughter)
{
	Eigen::Vector3f motherP = Eigen::Vector3f(mother.GetPx(), mother.GetPy(), mother.GetPz());
	Eigen::Vector3f daughterP = Eigen::Vector3f(daughter.GetPx(), daughter.GetPy(), daughter.GetPz());

	Eigen::Vector3f motherP_X_daughterP = motherP.cross(daughterP);

	float jT = (motherP_X_daughterP.norm()) / motherP.norm();

	return jT;
}

bool KFParticle_Tools::isInRange(float min, float value, float max)
{
	return min <= value && value <= max;
}

void KFParticle_Tools::removeDuplicates(std::vector<double> &v)
{
	auto end = v.end();
	for (auto it = v.begin(); it != end; ++it)
	{
		end = remove(it + 1, end, *it);
	}
	v.erase(end, v.end());
}

void KFParticle_Tools::removeDuplicates(std::vector<int> &v)
{
	auto end = v.end();
	for (auto it = v.begin(); it != end; ++it)
	{
		end = remove(it + 1, end, *it);
	}
	v.erase(end, v.end());
}

void KFParticle_Tools::removeDuplicates(std::vector<std::vector<int>> &v)
{
	auto end = v.end();
	for (auto it = v.begin(); it != end; ++it)
	{
		end = remove(it + 1, end, *it);
	}
	v.erase(end, v.end());
}

void KFParticle_Tools::removeDuplicates(std::vector<std::vector<std::string>> &v)
{
	auto end = v.end();
	for (auto it = v.begin(); it != end; ++it)
	{
		end = remove(it + 1, end, *it);
	}
	v.erase(end, v.end());
}

bool KFParticle_Tools::findParticle(const std::string &particle)
{
	bool particleFound = true;
	if (!TDatabasePDG::Instance()->GetParticle(particle.c_str()))
	{
		std::cout << "The particle, " << particle << " is not recognized" << std::endl;
		std::cout << "Check TDatabasePDG for a list of available particles" << std::endl;
		particleFound = false;
	}

	return particleFound;
}

int KFParticle_Tools::getParticleID(const std::string &particle)
{
	return TDatabasePDG::Instance()->GetParticle(particle.c_str())->PdgCode();
}

float KFParticle_Tools::getParticleMass(const std::string &particle)
{
	return TDatabasePDG::Instance()->GetParticle(particle.c_str())->Mass();
}

float KFParticle_Tools::getParticleMass(const int PDGID)
{
	return TDatabasePDG::Instance()->GetParticle(PDGID)->Mass();
}

void KFParticle_Tools::identify(KFParticle particle)
{
	std::cout << "Track ID: " << particle.Id() << std::endl;
	std::cout << "PDG ID: " << particle.GetPDG() << ", charge: " << (int) particle.GetQ() << ", mass: " << particle.GetMass() << " GeV" << std::endl;
	std::cout << "(px,py,pz) = (" << particle.GetPx() << " +/- " << sqrt(particle.GetCovariance(3, 3)) << ", ";
	std::cout << particle.GetPy() << " +/- " << sqrt(particle.GetCovariance(4, 4)) << ", ";
	std::cout << particle.GetPz() << " +/- " << sqrt(particle.GetCovariance(5, 5)) << ") GeV" << std::endl;
	std::cout << "(x,y,z) = (" << particle.GetX() << " +/- " << sqrt(particle.GetCovariance(0, 0)) << ", ";
	std::cout << particle.GetY() << " +/- " << sqrt(particle.GetCovariance(1, 1)) << ", ";
	std::cout << particle.GetZ() << " +/- " << sqrt(particle.GetCovariance(2, 2)) << ") cm\n"
		<< std::endl;
}


inline Acts::BoundTrackParameters makeTrackParams(SvtxTrack* track)
{
	Acts::Vector3 momentum(track->get_px(), track->get_py(), track->get_pz());

	auto actsVertex = getVertex(track);
	auto perigee = Acts::Surface::makeShared<Acts::PerigeeSurface>(actsVertex);
	auto actsFourPos =
		Acts::Vector4(track->get_x() * Acts::UnitConstants::cm,
				track->get_y() * Acts::UnitConstants::cm,
				track->get_z() * Acts::UnitConstants::cm,
				10 * Acts::UnitConstants::ns);

	ActsTransformations transformer;

	Acts::BoundSymMatrix cov = transformer.rotateSvtxTrackCovToActs(track);

	return ActsExamples::TrackParameters::create(perigee, _tGeometry->geometry().getGeoContext(),
			actsFourPos, momentum,
			track->get_charge() / track->get_p(),
			cov)
		.value();
}

inline Acts::Vector3 getVertex(SvtxTrack* track)
{
	auto vertexId = track->get_vertex_id();
	const SvtxVertex* svtxVertex = m_vertexMap->get(vertexId);
	Acts::Vector3 vertex = Acts::Vector3::Zero();
	if (svtxVertex)
	{
		vertex(0) = svtxVertex->get_x() * Acts::UnitConstants::cm;
		vertex(1) = svtxVertex->get_y() * Acts::UnitConstants::cm;
		vertex(2) = svtxVertex->get_z() * Acts::UnitConstants::cm;
	}

	return vertex;
}


inline BoundTrackParamResult propagateTrack(
		const Acts::BoundTrackParameters& params,
		const SurfacePtr& targetSurf)
{

	//	std::cout << "Inside propagateTrack -- Pass 0 " << std::endl;
	using Stepper = Acts::EigenStepper<>;
	using Propagator = Acts::Propagator<Stepper>;
	/*
	   if(_tGeometry){
	   std::cout << "Yes - _tGeometry" << std::endl;

	   }else{
	   std::cout << "No - _tGeometry" << std::endl;
	   }
	   */
	auto field = _tGeometry->geometry().magField;
	//	std::cout << "Inside propagateTrack -- Pass 1 " << std::endl;

	Stepper stepper(field);
	Propagator propagator(stepper);

	Acts::Logging::Level logLevel = Acts::Logging::INFO;

	//	std::cout << "Inside propagateTrack -- Pass 2 " << std::endl;

	auto logger = Acts::getDefaultLogger("PHActsTrackProjection", logLevel);

	Acts::PropagatorOptions<> options(_tGeometry->geometry().getGeoContext(),
			_tGeometry->geometry().magFieldContext,
			Acts::LoggerWrapper{*logger});
	//	std::cout << "Inside propagateTrack -- Pass 3 " << std::endl;

	auto result = propagator.propagate(params, *targetSurf, options);
	if (result.ok())
	{
		return Acts::Result<BoundTrackParam>::success(std::move((*result).endParameters.value()));
	}

	//	std::cout << "Inside propagateTrack -- Pass 4 " << std::endl;

	return result.error();
}
