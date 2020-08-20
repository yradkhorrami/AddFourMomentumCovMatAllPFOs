#include "AddFourMomentumCovMatAllPFOs.h"
#include <iostream>
#include <EVENT/LCCollection.h>
#include "EVENT/LCCollection.h"
#include "IMPL/LCCollectionVec.h"
#include "EVENT/MCParticle.h"
#include "EVENT/ReconstructedParticle.h"
#include <IMPL/ReconstructedParticleImpl.h>
#include "IMPL/ParticleIDImpl.h"
#include "UTIL/PIDHandler.h"
#include "marlin/VerbosityLevels.h"
#include <GeometryUtil.h>
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TTree.h"

using namespace lcio ;
using namespace marlin ;
using namespace std ;

AddFourMomentumCovMatAllPFOs aAddFourMomentumCovMatAllPFOs;

AddFourMomentumCovMatAllPFOs::AddFourMomentumCovMatAllPFOs() :

Processor("AddFourMomentumCovMatAllPFOs"),
m_nRun(0),
m_nEvt(0),
m_nRunSum(0),
m_nEvtSum(0),
m_Bfield(0.f),
c(0.),
mm2m(0.),
eV2GeV(0.),
eB(0.),
proton_mass(0.),
proton_mass_sq(0.),
kaon_mass(0.),
kaon_mass_sq(0.),
pion_mass(0.),
pion_mass_sq(0.),
m_pTFile(NULL),
m_pTTree(NULL),
h_nTracks_PFOCharge(NULL),
h_nClusters_nTracks(NULL),
h_SigmaPx2nT(NULL),
h_SigmaPxPynT(NULL),
h_SigmaPy2nT(NULL),
h_SigmaPxPznT(NULL),
h_SigmaPyPznT(NULL),
h_SigmaPz2nT(NULL),
h_SigmaPxEnT(NULL),
h_SigmaPyEnT(NULL),
h_SigmaPzEnT(NULL),
h_SigmaE2nT(NULL),
h_SigmaPx2T(NULL),
h_SigmaPxPyT(NULL),
h_SigmaPy2T(NULL),
h_SigmaPxPzT(NULL),
h_SigmaPyPzT(NULL),
h_SigmaPz2T(NULL),
h_SigmaPxET(NULL),
h_SigmaPyET(NULL),
h_SigmaPzET(NULL),
h_SigmaE2T(NULL),
h_SigmaPx2(NULL),
h_SigmaPxPy(NULL),
h_SigmaPy2(NULL),
h_SigmaPxPz(NULL),
h_SigmaPyPz(NULL),
h_SigmaPz2(NULL),
h_SigmaPxE(NULL),
h_SigmaPyE(NULL),
h_SigmaPzE(NULL),
h_SigmaE2(NULL)

{
	_description = "Set the convariance matrix in (P,E) for all pfos (charged particles, neutral hadrons and photons)";

	registerInputCollection(	LCIO::RECONSTRUCTEDPARTICLE,
					"inputPfoCollection",
					"Name of input pfo collection",
					m_inputPfoCollection,
					std::string("PandoraPFOs")
				);

	registerInputCollection(	LCIO::TRACK,
					"MarlinTrkTracksCollection" ,
					"Name of the MarlinTrkTracks collection"  ,
					m_MarlinTrkTracks,
					std::string("MarlinTrkTracks")
				);

	registerInputCollection(	LCIO::TRACK,
					"MarlinTrkTracksCollectionKaon" ,
					"Name of the MarlinTrkTracks collection"  ,
					m_MarlinTrkTracksKAON,
					std::string("MarlinTrkTracksKaon")
				);

	registerInputCollection(	LCIO::TRACK,
					"MarlinTrkTracksCollectionProton" ,
					"Name of the MarlinTrkTracks collection"  ,
					m_MarlinTrkTracksPROTON,
					std::string("MarlinTrkTracksProton")
				);

	registerOutputCollection(	LCIO::RECONSTRUCTEDPARTICLE,
					"outputPfoCollection",
					"Name of output pfo collection",
					m_outputPfoCollection,
					std::string("CorrectedPfoCollection")
				);

	registerProcessorParameter(	"useClusterPositionError",
					"true: use cluster position error for CovMat, false: use cluster direction error for CovMat",
					m_useClusterPositionError,
					bool(true)
				);

	registerProcessorParameter(	"RootFile",
	                                "Name of the output root file",
					m_rootFile,
					std::string("FourMomentumCovMatAllPFOs.root")
				);

}

void AddFourMomentumCovMatAllPFOs::init()
{

	streamlog_out(DEBUG) << "   init called  " << std::endl;
	m_Bfield = MarlinUtil::getBzAtOrigin();
	printParameters();
	streamlog_out(DEBUG) << " BField =  "<< m_Bfield << " Tesla" << std::endl ;
	m_nRun = 0 ;
	m_nEvt = 0 ;
	m_nRunSum = 0;
	m_nEvtSum = 0;
	c = 2.99792458e8;
	mm2m = 1e-3;
	eV2GeV = 1e-9;
	eB = m_Bfield * c * mm2m * eV2GeV;

	proton_mass = 0.938272088;
	proton_mass_sq = proton_mass * proton_mass;
	kaon_mass = 0.493677;
	kaon_mass_sq = kaon_mass * kaon_mass;
	pion_mass = 0.13957018;
	pion_mass_sq = pion_mass * pion_mass;

//	this->Clear();
	m_pTFile = new TFile(m_rootFile.c_str(), "recreate");
	m_pTTree = new TTree("CovMatAllPFOsTree", "CovMatAllPFOsTree");
	m_pTTree->SetDirectory(m_pTFile);
	h_nTracks_PFOCharge = new TH2I("h_nTracks_PFOCharge", "All PFOs; charge of PFO; nTracks", 5, -2.5, 2.5, 5, -0.5, 4.5);
	h_nTracks_PFOCharge->SetDirectory(m_pTFile);
	h_nClusters_nTracks = new TH2I("h_nClusters_nTracks", "Neutral PFOs; nClusters; nTracks", 5, -0.5, 4.5, 5, -0.5, 4.5);
	h_nClusters_nTracks->SetDirectory(m_pTFile);
	h_SigmaPx2nT = new TH2F("h_SigmaPx2nT", "Neutral PFOs (nTracks = 0); #sigma_{p_{x}}^{2} (new PFO); #sigma_{p_{x}}^{2} (old PFO)", 400, 0.0, 0.1, 400, 0.0, 0.1);
	h_SigmaPx2nT->SetDirectory(m_pTFile);
	h_SigmaPxPynT = new TH2F("h_SigmaPxPynT", "Neutral PFOs (nTracks = 0); #sigma_{p_{x}p_{y}} (new PFO); #sigma_{p_{x}p_{y}} (old PFO)", 400, -0.1, 0.1, 400, -0.1, 0.1);
	h_SigmaPxPynT->SetDirectory(m_pTFile);
	h_SigmaPy2nT = new TH2F("h_SigmaPy2nT", "Neutral PFOs (nTracks = 0); #sigma_{p_{y}}^{2} (new PFO); #sigma_{p_{y}}^{2} (old PFO)", 400, 0.0, 0.1, 400, 0.0, 0.1);
	h_SigmaPy2nT->SetDirectory(m_pTFile);
	h_SigmaPxPznT = new TH2F("h_SigmaPxPznT", "Neutral PFOs (nTracks = 0); #sigma_{p_{x}p_{z}} (new PFO); #sigma_{p_{x}p_{z}} (old PFO)", 400, -0.1, 0.1, 400, -0.1, 0.1);
	h_SigmaPxPznT->SetDirectory(m_pTFile);
	h_SigmaPyPznT = new TH2F("h_SigmaPyPznT", "Neutral PFOs (nTracks = 0); #sigma_{p_{y}p_{z}} (new PFO); #sigma_{p_{y}p_{z}} (old PFO)", 400, -0.1, 0.1, 400, -0.1, 0.1);
	h_SigmaPyPznT->SetDirectory(m_pTFile);
	h_SigmaPz2nT = new TH2F("h_SigmaPz2nT", "Neutral PFOs (nTracks = 0); #sigma_{p_{z}}^{2} (new PFO); #sigma_{p_{z}}^{2} (old PFO)", 400, 0.0, 0.1, 400, 0.0, 0.1);
	h_SigmaPz2nT->SetDirectory(m_pTFile);
	h_SigmaPxEnT = new TH2F("h_SigmaPxEnT", "Neutral PFOs (nTracks = 0); #sigma_{p_{x}E} (new PFO); #sigma_{p_{x}E} (old PFO)", 400, -0.1, 0.1, 400, -0.1, 0.1);
	h_SigmaPxEnT->SetDirectory(m_pTFile);
	h_SigmaPyEnT = new TH2F("h_SigmaPyEnT", "Neutral PFOs (nTracks = 0); #sigma_{p_{y}E} (new PFO); #sigma_{p_{y}E} (old PFO)", 400, -0.1, 0.1, 400, -0.1, 0.1);
	h_SigmaPyEnT->SetDirectory(m_pTFile);
	h_SigmaPzEnT = new TH2F("h_SigmaPzEnT", "Neutral PFOs (nTracks = 0); #sigma_{p_{z}E} (new PFO); #sigma_{p_{z}E} (old PFO)", 400, -0.1, 0.1, 400, -0.1, 0.1);
	h_SigmaPzEnT->SetDirectory(m_pTFile);
	h_SigmaE2nT = new TH2F("h_SigmaE2nT", "Neutral PFOs (nTracks = 0); #sigma_{E}^{2} (new PFO); #sigma_{E}^{2} (old PFO)", 400, 0.0, 0.6, 400, 0.0, 0.6);
	h_SigmaE2nT->SetDirectory(m_pTFile);
	h_SigmaPx2T = new TH2F("h_SigmaPx2T", "Neutral PFOs (nTracks = 2); #sigma_{p_{x}}^{2} (new PFO); #sigma_{p_{x}}^{2} (old PFO)", 400, 0.0, 0.1, 400, 0.0, 0.1);
	h_SigmaPx2T->SetDirectory(m_pTFile);
	h_SigmaPxPyT = new TH2F("h_SigmaPxPyT", "Neutral PFOs (nTracks = 2); #sigma_{p_{x}p_{y}} (new PFO); #sigma_{p_{x}p_{y}} (old PFO)", 400, -0.1, 0.1, 400, -0.1, 0.1);
	h_SigmaPxPyT->SetDirectory(m_pTFile);
	h_SigmaPy2T = new TH2F("h_SigmaPy2T", "Neutral PFOs (nTracks = 2); #sigma_{p_{y}}^{2} (new PFO); #sigma_{p_{y}}^{2} (old PFO)", 400, 0.0, 0.1, 400, 0.0, 0.1);
	h_SigmaPy2T->SetDirectory(m_pTFile);
	h_SigmaPxPzT = new TH2F("h_SigmaPxPzT", "Neutral PFOs (nTracks = 2); #sigma_{p_{x}p_{z}} (new PFO); #sigma_{p_{x}p_{z}} (old PFO)", 400, -0.1, 0.1, 400, -0.1, 0.1);
	h_SigmaPxPzT->SetDirectory(m_pTFile);
	h_SigmaPyPzT = new TH2F("h_SigmaPyPzT", "Neutral PFOs (nTracks = 2); #sigma_{p_{y}p_{z}} (new PFO); #sigma_{p_{y}p_{z}} (old PFO)", 400, -0.1, 0.1, 400, -0.1, 0.1);
	h_SigmaPyPzT->SetDirectory(m_pTFile);
	h_SigmaPz2T = new TH2F("h_SigmaPz2T", "Neutral PFOs (nTracks = 2); #sigma_{p_{z}}^{2} (new PFO); #sigma_{p_{z}}^{2} (old PFO)", 400, 0.0, 0.1, 400, 0.0, 0.1);
	h_SigmaPz2T->SetDirectory(m_pTFile);
	h_SigmaPxET = new TH2F("h_SigmaPxET", "Neutral PFOs (nTracks = 2); #sigma_{p_{x}E} (new PFO); #sigma_{p_{x}E} (old PFO)", 400, -0.1, 0.1, 400, -0.1, 0.1);
	h_SigmaPxET->SetDirectory(m_pTFile);
	h_SigmaPyET = new TH2F("h_SigmaPyET", "Neutral PFOs (nTracks = 2); #sigma_{p_{y}E} (new PFO); #sigma_{p_{y}E} (old PFO)", 400, -0.1, 0.1, 400, -0.1, 0.1);
	h_SigmaPyET->SetDirectory(m_pTFile);
	h_SigmaPzET = new TH2F("h_SigmaPzET", "Neutral PFOs (nTracks = 2); #sigma_{p_{z}E} (new PFO); #sigma_{p_{z}E} (old PFO)", 400, -0.1, 0.1, 400, -0.1, 0.1);
	h_SigmaPzET->SetDirectory(m_pTFile);
	h_SigmaE2T = new TH2F("h_SigmaE2T", "Neutral PFOs (nTracks = 2); #sigma_{E}^{2} (new PFO); #sigma_{E}^{2} (old PFO)", 400, 0.0, 0.6, 400, 0.0, 0.6);
	h_SigmaE2T->SetDirectory(m_pTFile);
	h_SigmaPx2 = new TH2F("h_SigmaPx2", "Charged PFOs; #sigma_{p_{x}}^{2} (new PFO); #sigma_{p_{x}}^{2} (old PFO)", 400, 0.0, 0.1, 400, 0.0, 0.1);
	h_SigmaPx2->SetDirectory(m_pTFile);
	h_SigmaPxPy = new TH2F("h_SigmaPxPy", "Charged PFOs; #sigma_{p_{x}p_{y}} (new PFO); #sigma_{p_{x}p_{y}} (old PFO)", 400, -0.1, 0.1, 400, -0.1, 0.1);
	h_SigmaPxPy->SetDirectory(m_pTFile);
	h_SigmaPy2 = new TH2F("h_SigmaPy2", "Charged PFOs; #sigma_{p_{y}}^{2} (new PFO); #sigma_{p_{y}}^{2} (old PFO)", 400, 0.0, 0.1, 400, 0.0, 0.1);
	h_SigmaPy2->SetDirectory(m_pTFile);
	h_SigmaPxPz = new TH2F("h_SigmaPxPz", "Charged PFOs; #sigma_{p_{x}p_{z}} (new PFO); #sigma_{p_{x}p_{z}} (old PFO)", 400, -0.1, 0.1, 400, -0.1, 0.1);
	h_SigmaPxPz->SetDirectory(m_pTFile);
	h_SigmaPyPz = new TH2F("h_SigmaPyPz", "Charged PFOs; #sigma_{p_{y}p_{z}} (new PFO); #sigma_{p_{y}p_{z}} (old PFO)", 400, -0.1, 0.1, 400, -0.1, 0.1);
	h_SigmaPyPz->SetDirectory(m_pTFile);
	h_SigmaPz2 = new TH2F("h_SigmaPz2", "Charged PFOs; #sigma_{p_{z}}^{2} (new PFO); #sigma_{p_{z}}^{2} (old PFO)", 400, 0.0, 0.1, 400, 0.0, 0.1);
	h_SigmaPz2->SetDirectory(m_pTFile);
	h_SigmaPxE = new TH2F("h_SigmaPxE", "Charged PFOs; #sigma_{p_{x}E} (new PFO); #sigma_{p_{x}E} (old PFO)", 400, -0.1, 0.1, 400, -0.1, 0.1);
	h_SigmaPxE->SetDirectory(m_pTFile);
	h_SigmaPyE = new TH2F("h_SigmaPyE", "Charged PFOs; #sigma_{p_{y}E} (new PFO); #sigma_{p_{y}E} (old PFO)", 400, -0.1, 0.1, 400, -0.1, 0.1);
	h_SigmaPyE->SetDirectory(m_pTFile);
	h_SigmaPzE = new TH2F("h_SigmaPzE", "Charged PFOs; #sigma_{p_{z}E} (new PFO); #sigma_{p_{z}E} (old PFO)", 400, -0.1, 0.1, 400, -0.1, 0.1);
	h_SigmaPzE->SetDirectory(m_pTFile);
	h_SigmaE2 = new TH2F("h_SigmaE2", "Charged PFOs; #sigma_{E}^{2} (new PFO); #sigma_{E}^{2} (old PFO)", 400, 0.0, 0.6, 400, 0.0, 0.6);
	h_SigmaE2->SetDirectory(m_pTFile);

}

void AddFourMomentumCovMatAllPFOs::Clear()
{

}

void AddFourMomentumCovMatAllPFOs::processRunHeader()
{

	m_nRun = 0;
	m_nEvt = 0;
	++m_nRunSum;

}

void AddFourMomentumCovMatAllPFOs::processEvent( EVENT::LCEvent *pLCEvent )
{

	m_nRun = pLCEvent->getRunNumber();
	m_nEvt = pLCEvent->getEventNumber();
	++m_nEvtSum;

	LCCollection *inputPfoCollection{};
	this->Clear();

        try
        {
		inputPfoCollection = pLCEvent->getCollection(m_inputPfoCollection);
		int n_PFO = inputPfoCollection->getNumberOfElements();
		LCCollectionVec *m_col_outputPfo = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
		for ( int i_pfo = 0; i_pfo < n_PFO ; ++i_pfo )
		{
			ReconstructedParticle* inputPFO = dynamic_cast<ReconstructedParticleImpl*>(inputPfoCollection->getElementAt(i_pfo));
			ReconstructedParticleImpl* outputPFO = new ReconstructedParticleImpl;
			double outputPFOMomentum[3]{0., 0., 0.};
			double outputPFOEnergy = 0.;
			float pfoCharge = inputPFO->getCharge();
			
			TrackVec pfoTracks	= inputPFO->getTracks();
			int nTrackspfo		= pfoTracks.size();
			int nClusterspfo	= (inputPFO->getClusters()).size();
			streamlog_out(DEBUG) << "Number of tracks assigned to PFO : " << nTrackspfo << std::endl;
			streamlog_out(DEBUG) << "Number of clusters assigned to PFO : " << nClusterspfo << std::endl;
			h_nTracks_PFOCharge->Fill( pfoCharge , nTrackspfo );
			TLorentzVector pfoFourMomentum(0.,0.,0.,0.);
			TVector3 clusterPosition(0.,0.,0.);
			float pfoMass		= inputPFO->getMass();
			std::vector<float> outputCovMatrix( 10, 0.0 );
			std::vector<float> inputCovMatrix = inputPFO->getCovMatrix();//( 10, 0.0 );
			if ( pfoCharge == 0)
			{
				h_nClusters_nTracks->Fill( nClusterspfo , nTrackspfo );
				if ( nTrackspfo == 0 )
				{
					float pfoEnergy		= inputPFO->getEnergy();
//					float clusterE		= ( inputPFO->getClusters()[0] )->getEnergy();
					float clusterTheta	= ( inputPFO->getClusters()[0] )->getITheta();
					float clusterPhi	= ( inputPFO->getClusters()[0] )->getIPhi();
					clusterPosition		= TVector3( ( inputPFO->getClusters()[0] )->getPosition()[0] , ( inputPFO->getClusters()[0] )->getPosition()[1] , ( inputPFO->getClusters()[0] )->getPosition()[2] );
					float pfoMomentumMag	= std::sqrt( pow( pfoEnergy , 2 ) - pow( pfoMass , 2 ) );
					float pfoPx		= pfoMomentumMag * sin( clusterTheta ) * cos( clusterPhi );
					float pfoPy		= pfoMomentumMag * sin( clusterTheta ) * sin( clusterPhi );
					float pfoPz		= pfoMomentumMag * cos( clusterTheta );
					TVector3 pfoMomentum( pfoPx , pfoPy , pfoPz );
					pfoFourMomentum		= TLorentzVector( pfoMomentum , pfoEnergy );
					std::vector<float> clusterDirectionError = ( inputPFO->getClusters()[0] )->getDirectionError();
					std::vector<float> clusterPositionError = ( inputPFO->getClusters()[0] )->getPositionError();
					float clusterEnergyError= ( inputPFO->getClusters()[0] )->getEnergyError();
					if ( m_useClusterPositionError )
					{
						outputCovMatrix		= this->UpdateNeutralPFOCovMatPosError( clusterPosition , pfoEnergy , pfoMass , clusterPositionError , clusterEnergyError );
					}
					else
					{
						outputCovMatrix		= this->UpdateNeutralPFOCovMatDirError( pfoFourMomentum , clusterDirectionError , clusterEnergyError );
					}
					streamlog_out(DEBUG) << "PFO is neutral (without track), CovMatrix is set using cluster information" << std::endl;
					h_SigmaPx2nT->Fill( outputCovMatrix[0] , inputCovMatrix[0] );
					h_SigmaPxPynT->Fill( outputCovMatrix[1] , inputCovMatrix[1] );
					h_SigmaPy2nT->Fill( outputCovMatrix[2] , inputCovMatrix[2] );
					h_SigmaPxPznT->Fill( outputCovMatrix[3] , inputCovMatrix[3] );
					h_SigmaPyPznT->Fill( outputCovMatrix[4] , inputCovMatrix[4] );
					h_SigmaPz2nT->Fill( outputCovMatrix[5] , inputCovMatrix[5] );
					h_SigmaPxEnT->Fill( outputCovMatrix[6] , inputCovMatrix[6] );
					h_SigmaPyEnT->Fill( outputCovMatrix[7] , inputCovMatrix[7] );
					h_SigmaPzEnT->Fill( outputCovMatrix[8] , inputCovMatrix[8] );
					h_SigmaE2nT->Fill( outputCovMatrix[9] , inputCovMatrix[9] );
				}
				else if ( nTrackspfo == 2 )
				{
					std::vector<float> pfoTrackCovMat( 10, 0.0 );
					for ( int i_trk = 0 ; i_trk < nTrackspfo ; ++i_trk )
					{
						Track* mytrack		= pfoTracks[ i_trk ];
						double track_mass = this->getTrackMass( pLCEvent , mytrack );
						if ( track_mass == 0.0 )
						{
							track_mass = pion_mass;
							streamlog_out(MESSAGE) << "Couldn't Find track mass, default mass (pion mass) is set" << std::endl;
						}
						streamlog_out(DEBUG) << "PFO is neutral (with two tracks), CovMatrix is set using cluster information" << std::endl;
						double track_omega	= mytrack->getOmega();
						double track_tanL	= mytrack->getTanLambda();
						double track_phi	= mytrack->getPhi();
						double track_pt		= eB / fabs(track_omega);
						double track_px		= track_pt * TMath::Cos(track_phi);
						double track_py		= track_pt * TMath::Sin(track_phi);
						double track_pz		= track_pt * track_tanL;
						double track_energy	= std::sqrt( pow( track_px , 2 ) + pow( track_py , 2 ) + pow( track_pz , 2 ) + pow( track_mass , 2 ) );
						TLorentzVector trackFourMomentum( track_px , track_py , track_pz , track_energy );
						pfoTrackCovMat		= this->UpdateChargedPFOCovMat( mytrack , track_mass );
						for ( int i = 0 ; i < 10 ; ++i )
						{
							outputCovMatrix[i] += pfoTrackCovMat[i];
						}
						pfoFourMomentum		+= trackFourMomentum;
					}
					h_SigmaPx2T->Fill( outputCovMatrix[0] , inputCovMatrix[0] );
					h_SigmaPxPyT->Fill( outputCovMatrix[1] , inputCovMatrix[1] );
					h_SigmaPy2T->Fill( outputCovMatrix[2] , inputCovMatrix[2] );
					h_SigmaPxPzT->Fill( outputCovMatrix[3] , inputCovMatrix[3] );
					h_SigmaPyPzT->Fill( outputCovMatrix[4] , inputCovMatrix[4] );
					h_SigmaPz2T->Fill( outputCovMatrix[5] , inputCovMatrix[5] );
					h_SigmaPxET->Fill( outputCovMatrix[6] , inputCovMatrix[6] );
					h_SigmaPyET->Fill( outputCovMatrix[7] , inputCovMatrix[7] );
					h_SigmaPzET->Fill( outputCovMatrix[8] , inputCovMatrix[8] );
					h_SigmaE2T->Fill( outputCovMatrix[9] , inputCovMatrix[9] );
				}
			}
			else
			{
				std::vector<float> pfoTrackCovMat( 10, 0.0 );
				for ( int i_trk = 0 ; i_trk < nTrackspfo ; ++i_trk )
				{
					Track* mytrack		= pfoTracks[ i_trk ];
					double track_mass = this->getTrackMass( pLCEvent , mytrack );
					if ( track_mass == 0.0 )
					{
						track_mass = pion_mass;
						streamlog_out(MESSAGE) << "Couldn't Find track mass, default mass (pion mass) is set" << std::endl;
					}
					streamlog_out(DEBUG) << "PFO is neutral (with two tracks), CovMatrix is set using cluster information" << std::endl;
					double track_omega	= mytrack->getOmega();
					double track_tanL	= mytrack->getTanLambda();
					double track_phi	= mytrack->getPhi();
					double track_pt		= eB / fabs(track_omega);
					double track_px		= track_pt * TMath::Cos(track_phi);
					double track_py		= track_pt * TMath::Sin(track_phi);
					double track_pz		= track_pt * track_tanL;
					double track_energy	= std::sqrt( pow( track_px , 2 ) + pow( track_py , 2 ) + pow( track_pz , 2 ) + pow( track_mass , 2 ) );
					TLorentzVector trackFourMomentum( track_px , track_py , track_pz , track_energy );
					pfoTrackCovMat		= this->UpdateChargedPFOCovMat( mytrack , track_mass );
					for ( int i = 0 ; i < 10 ; ++i )
					{
						outputCovMatrix[i] += pfoTrackCovMat[i];
					}
					pfoFourMomentum		+= trackFourMomentum;
				}
				h_SigmaPx2->Fill( outputCovMatrix[0] , inputCovMatrix[0] );
				h_SigmaPxPy->Fill( outputCovMatrix[1] , inputCovMatrix[1] );
				h_SigmaPy2->Fill( outputCovMatrix[2] , inputCovMatrix[2] );
				h_SigmaPxPz->Fill( outputCovMatrix[3] , inputCovMatrix[3] );
				h_SigmaPyPz->Fill( outputCovMatrix[4] , inputCovMatrix[4] );
				h_SigmaPz2->Fill( outputCovMatrix[5] , inputCovMatrix[5] );
				h_SigmaPxE->Fill( outputCovMatrix[6] , inputCovMatrix[6] );
				h_SigmaPyE->Fill( outputCovMatrix[7] , inputCovMatrix[7] );
				h_SigmaPzE->Fill( outputCovMatrix[8] , inputCovMatrix[8] );
				h_SigmaE2->Fill( outputCovMatrix[9] , inputCovMatrix[9] );
			}
			outputPFOMomentum[0]	=	pfoFourMomentum.Px();
			outputPFOMomentum[1]	=	pfoFourMomentum.Py();
			outputPFOMomentum[2]	=	pfoFourMomentum.Pz();
			outputPFOEnergy		=	pfoFourMomentum.E();
			outputPFO->setType(inputPFO->getType());
			outputPFO->setMomentum(outputPFOMomentum);
			outputPFO->setEnergy(outputPFOEnergy);
			outputPFO->setCovMatrix(outputCovMatrix);
			outputPFO->setMass(inputPFO->getMass());
			outputPFO->setCharge(inputPFO->getCharge());
			outputPFO->setReferencePoint(inputPFO->getReferencePoint());
			for (unsigned int j=0; j<inputPFO->getParticleIDs().size(); ++j)
			{
				ParticleIDImpl* inPID = dynamic_cast<ParticleIDImpl*>(inputPFO->getParticleIDs()[j]);
			        ParticleIDImpl* outPID = new ParticleIDImpl;
			        outPID->setType(inPID->getType());
			        outPID->setPDG(inPID->getPDG());
			        outPID->setLikelihood(inPID->getLikelihood());
			        outPID->setAlgorithmType(inPID->getAlgorithmType()) ;
			        for (unsigned int k=0; k<inPID->getParameters().size()  ; ++k) outPID->addParameter(inPID->getParameters()[k]) ;
			        outputPFO->addParticleID(outPID);
			}
			outputPFO->setParticleIDUsed(inputPFO->getParticleIDUsed());
			outputPFO->setGoodnessOfPID(inputPFO->getGoodnessOfPID());
			for (unsigned int j=0; j<inputPFO->getParticles().size(); ++j)
			{
				outputPFO->addParticle(inputPFO->getParticles()[j]);
			}
			for (unsigned int j=0; j<inputPFO->getClusters().size(); ++j)
			{
				outputPFO->addCluster(inputPFO->getClusters()[j]);
			}
			for (unsigned int j=0; j<inputPFO->getTracks().size(); ++j)
			{
				outputPFO->addTrack(inputPFO->getTracks()[j]);
			}
			outputPFO->setStartVertex(inputPFO->getStartVertex());
			m_col_outputPfo->addElement( outputPFO );
		}
		pLCEvent->addCollection( m_col_outputPfo , m_outputPfoCollection );
        }
        catch(DataNotAvailableException &e)
        {
          streamlog_out(MESSAGE) << "Input collection not found in event " << m_nEvt << std::endl;
        }

}

std::vector<float> AddFourMomentumCovMatAllPFOs::UpdateNeutralPFOCovMatDirError( TLorentzVector pfoFourMomentum , std::vector<float> clusterDirectionError , float clusterEnergyError )
{

//	Obtain covariance matrix on (px,py,pz,E) from the
//	covariance matrix on cluster parameters.
//
//	define the jacobian as the 3x4 matrix:
//
//
//
//			Dpx/DTheta	Dpy/DTheta	Dpz/DTheta	DE/DTheta
//
//	 J =		Dpx/DPhi	Dpy/DPhi	Dpz/DPhi	DE/DPhi
//
//			Dpx/DE		Dpy/DE		Dpz/DE		DE/DE
//
//
//
//			Px.Pz/Pt	Py.Pz/Pt	-Pt		0
//
//	J =		-Py		Px		0		0
//
//			Px.E/P2		Py.E/P2		Pz.E/P2		1
//
//
//
//	CovMatrix elements in terms of cluster direction error and cluster energy error:
//
//			Theta.Theta	Theta.phi	Theta.E
//
//	Cov =		phi.Theta	phi.phi		phi.E
//
//			E.Theta		E.phi		E.E
//
//
//

	const int rows			= 3; // n rows jacobian
	const int columns		= 4; // n columns jacobian
	const int kspace_time_dim	= 4;

	TMatrixD covMatrixMomenta(kspace_time_dim,kspace_time_dim);
	std::vector<float> covP;

	float pfoPx		=	pfoFourMomentum.Px();
	float pfoPy		=	pfoFourMomentum.Py();
	float pfoPz		=	pfoFourMomentum.Pz();
	float pfoE		=	pfoFourMomentum.E();
	float pfoPt		=	std::sqrt( pow( pfoPx , 2 ) + pow( pfoPy , 2 ) );
	float pfoP2		=	std::sqrt( pow( pfoPx , 2 ) + pow( pfoPy , 2 ) + pow( pfoPz , 2 ) );
	float SigmaTheta2	=	clusterDirectionError[ 0 ];
	float SigmaThetaSigmaPhi=	clusterDirectionError[ 1 ];
	float SigmaPhi2		=	clusterDirectionError[ 2 ];
	float SigmaE2		=	pow( clusterEnergyError , 2 );

	streamlog_out(DEBUG) << "Cluster information obtained" << std::endl;

//	Define array with jacobian matrix elements by rows
	double jacobian_by_rows[rows*columns] =
			{
				pfoPx * pfoPz / pfoPt	,	pfoPy * pfoPz / pfoPt	,	-pfoPt			,	0	,
				-pfoPy			,	pfoPx			,	0			,	0	,
				pfoPx * pfoE / pfoP2	,	pfoPy * pfoE / pfoP2	,	pfoPz * pfoE / pfoP2	,	1	 
			};

	streamlog_out(DEBUG) << "Jacobian array formed by rows" << std::endl;

//	construct the Jacobian using previous array ("F" if filling by columns, "C", if filling by rows, $ROOTSYS/math/matrix/src/TMatrixT.cxx)
	TMatrixD jacobian(rows,columns, jacobian_by_rows, "C");
	streamlog_out(DEBUG) << "Jacobian array converted to Jacobian matrix" << std::endl;

//	cluster covariance matrix by rows
	double cluster_cov_matrix_by_rows[rows*rows] =
			{
				SigmaTheta2		,	SigmaThetaSigmaPhi	,	0	,
				SigmaThetaSigmaPhi	,	SigmaPhi2		,	0	,
				0			,	0			,	SigmaE2	
			};
	streamlog_out(DEBUG) << "cluster covariance matrix array formed by rows" << std::endl;

	TMatrixD covMatrix_cluster(rows,rows, cluster_cov_matrix_by_rows, "C");
	streamlog_out(DEBUG) << "cluster covariance matrix array converted to cluster covariance matrix matrix" << std::endl;

	covMatrixMomenta.Mult( TMatrixD( jacobian ,
					TMatrixD::kTransposeMult ,
					covMatrix_cluster) ,
					jacobian
					);
	streamlog_out(DEBUG) << "cluster covariance matrix array converted to FourMomentumCovariance matrix" << std::endl;

	covP.push_back( covMatrixMomenta(0,0) ); // x-x
	covP.push_back( covMatrixMomenta(1,0) ); // y-x
	covP.push_back( covMatrixMomenta(1,1) ); // y-y
	covP.push_back( covMatrixMomenta(2,0) ); // z-x
	covP.push_back( covMatrixMomenta(2,1) ); // z-y
	covP.push_back( covMatrixMomenta(2,2) ); // z-z
	covP.push_back( covMatrixMomenta(3,0) ); // e-x
	covP.push_back( covMatrixMomenta(3,1) ); // e-y
	covP.push_back( covMatrixMomenta(3,2) ); // e-z
	covP.push_back( covMatrixMomenta(3,3) ); // e-e
	streamlog_out(DEBUG) << "FourMomentumCovarianceMatrix Filled succesfully" << std::endl;
	
	return covP;

}


std::vector<float> AddFourMomentumCovMatAllPFOs::UpdateNeutralPFOCovMatPosError( TVector3 clusterPosition , float pfoEnergy , float pfoMass , std::vector<float> clusterPositionError , float clusterEnergyError )
{

//	Obtain covariance matrix on (px,py,pz,E) from the
//	covariance matrix on cluster parameters.
//
//	define the jacobian as the 4x4 matrix:
//
//
//
//			Dpx/Dx			Dpy/Dx			Dpz/Dx			DE/Dx
//
//			Dpx/Dy			Dpy/Dy			Dpz/Dy			DE/Dy
//	J =
//			Dpx/Dz			Dpy/Dz			Dpz/Dz			DE/Dz
//
//			Dpx/DE			Dpy/DE			Dpz/DE			DE/DE
//
//
//
//
//
//			P.(r2-x2)/r3		-P.x.y/r3		-P.x.z/r3		0
//
//			-P.y.x/r3		P.(r2-y2)/r3		-P.y.z/r3		0
//	J =
//			-P.z.x/r3		-P.z.y/r3		P.(r2-z2)/r3		0
//
//			(E/P).(x/r)		(E/P).(y/r)		(E/P).(z/r)		1
//
//
//
//
//	CovMatrix elements in terms of cluster position error and cluster energy error:
//
//			x.x		x.y		x.z		x.E
//
//			y.x		y.y		y.z		y.E
//	Cov =
//			z.x		z.y		z.z		z.E
//
//			E.x		E.y		E.z		E.E
//
//
//

	const int rows			= 4; // n rows jacobian
	const int columns		= 4; // n columns jacobian
	const int kspace_time_dim	= 4;

	TMatrixD covMatrixMomenta(kspace_time_dim,kspace_time_dim);
	std::vector<float> covP;

	float pfoX		=	clusterPosition.X();
	float pfoY		=	clusterPosition.Y();
	float pfoZ		=	clusterPosition.Z();
	float pfoR		=	std::sqrt( pow( pfoX , 2 ) + pow( pfoY , 2 ) + pow( pfoZ , 2 ) );
	float pfoX2		=	pow( pfoX , 2 );
	float pfoY2		=	pow( pfoY , 2 );
	float pfoZ2		=	pow( pfoZ , 2 );
	float pfoR2		=	pow( pfoR , 2 );
	float pfoR3		=	pow( pfoR , 3 );
	float pfoE		=	pfoEnergy;
	float pfoP		=	std::sqrt( pow( pfoE , 2 ) - pow( pfoMass , 2 ) );
	float SigmaX2		=	clusterPositionError[ 0 ];
	float SigmaXY		=	clusterPositionError[ 1 ];
	float SigmaY2		=	clusterPositionError[ 2 ];
	float SigmaXZ		=	clusterPositionError[ 3 ];
	float SigmaYZ		=	clusterPositionError[ 4 ];
	float SigmaZ2		=	clusterPositionError[ 5 ];
	float SigmaE2		=	pow( clusterEnergyError , 2 );

	streamlog_out(DEBUG) << "Cluster information obtained" << std::endl;

//	Define array with jacobian matrix elements by rows
	double jacobian_by_rows[rows*columns] =
			{
				pfoP * ( pfoR2 - pfoX2 ) / pfoR3	,	-pfoP * pfoX * pfoY / pfoR3		,	-pfoP * pfoX * pfoZ / pfoR3		,	0	,
				-pfoP * pfoY * pfoX / pfoR3		,	pfoP * ( pfoR2 - pfoY2 ) / pfoR3	,	-pfoP * pfoY * pfoZ / pfoR3		,	0	,
				-pfoP * pfoZ * pfoX / pfoR3		,	-pfoP * pfoZ * pfoY / pfoR3		,	pfoP * ( pfoR2 - pfoZ2 ) / pfoR3	,	0	,
				pfoE * pfoX / ( pfoP * pfoR )		,	pfoE * pfoY / ( pfoP * pfoR )		,	pfoE * pfoZ / ( pfoP * pfoR )		,	1
			};

	streamlog_out(DEBUG) << "Jacobian array formed by rows" << std::endl;

//	construct the Jacobian using previous array ("F" if filling by columns, "C", if filling by rows, $ROOTSYS/math/matrix/src/TMatrixT.cxx)
	TMatrixD jacobian(rows,columns, jacobian_by_rows, "C");
	streamlog_out(DEBUG) << "Jacobian array converted to Jacobian matrix" << std::endl;

//	cluster covariance matrix by rows
	double cluster_cov_matrix_by_rows[rows*rows] =
			{
				SigmaX2		,	SigmaXY		,	SigmaXZ		,	0	,
				SigmaXY		,	SigmaY2		,	SigmaYZ		,	0	,
				SigmaXZ		,	SigmaYZ		,	SigmaZ2		,	0	,
				0		,	0		,	0		,	SigmaE2	
			};
	streamlog_out(DEBUG) << "cluster covariance matrix array formed by rows" << std::endl;

	TMatrixD covMatrix_cluster(rows,rows, cluster_cov_matrix_by_rows, "C");
	streamlog_out(DEBUG) << "cluster covariance matrix array converted to cluster covariance matrix matrix" << std::endl;

	covMatrixMomenta.Mult( TMatrixD( jacobian ,
					TMatrixD::kTransposeMult ,
					covMatrix_cluster) ,
					jacobian
					);
	streamlog_out(DEBUG) << "cluster covariance matrix array converted to FourMomentumCovariance matrix" << std::endl;

	covP.push_back( covMatrixMomenta(0,0) ); // x-x
	covP.push_back( covMatrixMomenta(1,0) ); // y-x
	covP.push_back( covMatrixMomenta(1,1) ); // y-y
	covP.push_back( covMatrixMomenta(2,0) ); // z-x
	covP.push_back( covMatrixMomenta(2,1) ); // z-y
	covP.push_back( covMatrixMomenta(2,2) ); // z-z
	covP.push_back( covMatrixMomenta(3,0) ); // e-x
	covP.push_back( covMatrixMomenta(3,1) ); // e-y
	covP.push_back( covMatrixMomenta(3,2) ); // e-z
	covP.push_back( covMatrixMomenta(3,3) ); // e-e
	streamlog_out(DEBUG) << "FourMomentumCovarianceMatrix Filled succesfully" << std::endl;
	
	return covP;

}


std::vector<float> AddFourMomentumCovMatAllPFOs::UpdateChargedPFOCovMat( EVENT::Track* MyTrack , float trackMass )
{

//	Obtain covariance matrix on (px,py,pz,E) from the
//	covariance matrix on track parameters.
//
//	define the jacobian as the 3x4 matrix:
//
//
//
//			Dpx/DTanLambda		Dpy/DTanLambda		Dpz/DTanLambda		DE/DTanLambda
//
//	 J =		Dpx/DPhi		Dpy/DPhi		Dpz/DPhi		DE/DPhi
//
//			Dpx/DOmega		Dpy/DOmega		Dpz/DOmega		DE/DOmega
//
//
//
//			0			0			Pt			Pz.Pt/E
//
//	J =		-Py			Px			0			0
//
//			-Px/Omega		-Py/Omega		-Pz/Omega		-P2/E.Omega
//
//
//
//	Order in the covariance matrix on helix parameters:
//
//			TanLambda.TanLambda	TanLambda.phi		TanLambda.Omega
//
//	Cov =		phi.TanLambda		phi.phi			phi.Omega
//
//			Omega.TanLambda		Omega.phi		Omega.Omega
//
//
//

	const int rows			= 3; // n rows jacobian
	const int columns		= 4; // n columns jacobian
	const int kspace_time_dim	= 4;

	TMatrixD covMatrixMomenta(kspace_time_dim,kspace_time_dim);
	std::vector<float> covP;

	double trackTanLambda		=	MyTrack->getTanLambda();
	double trackPhi			=	MyTrack->getPhi();
	double trackOmega		=	MyTrack->getOmega();
	std::vector<float> trackCovMat	=	MyTrack->getCovMatrix();
	double trackPt			=	eB / fabs( trackOmega );
	double trackPx			= 	trackPt * TMath::Cos( trackPhi );
	double trackPy			= 	trackPt * TMath::Sin( trackPhi );
	double trackPz			= 	trackPt * trackTanLambda;
	double trackP			= 	std::sqrt( pow( trackPt , 2 ) + pow( trackPz , 2 ) );
	double trackE			= 	std::sqrt( pow( trackP , 2 ) + pow( trackMass , 2 ) );
	
	float SigmaPhi2			=	trackCovMat[  2 ];
	float SigmaPhiSigmaOmega	=	trackCovMat[  4 ];
	float SigmaOmega2		=	trackCovMat[  5 ];
	float SigmaTanLambdaSigmaPhi	=	trackCovMat[ 11 ];
	float SigmaTanLambdaSigmaOmega	=	trackCovMat[ 12 ];
	float SigmaTanLambda2		=	trackCovMat[ 14 ];

	streamlog_out(DEBUG) << "track information obtained" << std::endl;

//	Define array with jacobian matrix elements by rows
	double jacobian_by_rows[rows*columns] =
	{
		0			,	0			,	trackPt			,	trackPz * trackPt / trackE			,
		-trackPy		,	trackPx			,	0			,	0						,
		-trackPx / trackOmega	,	-trackPy / trackOmega	,	-trackPz / trackOmega	,	-trackP * trackP / ( trackE * trackOmega )	 
	};

	streamlog_out(DEBUG) << "Jacobian array formed by rows" << std::endl;

//	construct the Jacobian using previous array ("F" if filling by columns, "C", if filling by rows, $ROOTSYS/math/matrix/src/TMatrixT.cxx)
	TMatrixD jacobian(rows,columns, jacobian_by_rows, "C");
	streamlog_out(DEBUG) << "Jacobian array converted to Jacobian matrix" << std::endl;

//	track covariance matrix by rows
	double track_cov_matrix_by_rows[rows*rows] =
	{
		SigmaTanLambda2			,	SigmaTanLambdaSigmaPhi		,	SigmaTanLambdaSigmaOmega	,
		SigmaTanLambdaSigmaPhi		,	SigmaPhi2			,	SigmaPhiSigmaOmega	,
		SigmaTanLambdaSigmaOmega	,	SigmaPhiSigmaOmega		,	SigmaOmega2	
	};
	streamlog_out(DEBUG) << "track covariance matrix array formed by rows" << std::endl;

	TMatrixD covMatrix_track(rows,rows, track_cov_matrix_by_rows, "C");
	streamlog_out(DEBUG) << "track covariance matrix array converted to track covariance matrix matrix" << std::endl;

	covMatrixMomenta.Mult( TMatrixD( jacobian ,
					TMatrixD::kTransposeMult ,
					covMatrix_track) ,
					jacobian
					);
	streamlog_out(DEBUG) << "track covariance matrix array converted to FourMomentumCovariance matrix" << std::endl;

	covP.push_back( covMatrixMomenta(0,0) ); // x-x
	covP.push_back( covMatrixMomenta(1,0) ); // y-x
	covP.push_back( covMatrixMomenta(1,1) ); // y-y
	covP.push_back( covMatrixMomenta(2,0) ); // z-x
	covP.push_back( covMatrixMomenta(2,1) ); // z-y
	covP.push_back( covMatrixMomenta(2,2) ); // z-z
	covP.push_back( covMatrixMomenta(3,0) ); // e-x
	covP.push_back( covMatrixMomenta(3,1) ); // e-y
	covP.push_back( covMatrixMomenta(3,2) ); // e-z
	covP.push_back( covMatrixMomenta(3,3) ); // e-e
	streamlog_out(DEBUG) << "FourMomentumCovarianceMatrix Filled succesfully" << std::endl;
	
	return covP;

}

double AddFourMomentumCovMatAllPFOs::getTrackMass( EVENT::LCEvent *pLCEvent, EVENT::Track* inputTrk )
{
	LCCollection *MarlinTrkTracks{};
	LCCollection *MarlinTrkTracksKAON{};
	LCCollection *MarlinTrkTracksPROTON{};
 	int n_TRK = -1;
	int n_TRKp = -1;
	int n_TRKk = -1;
	double trackMass = 0.0;
       try
        {
		MarlinTrkTracks = pLCEvent->getCollection(m_MarlinTrkTracks);
		MarlinTrkTracksKAON = pLCEvent->getCollection(m_MarlinTrkTracksKAON);
		MarlinTrkTracksPROTON = pLCEvent->getCollection(m_MarlinTrkTracksPROTON);
		n_TRK = MarlinTrkTracks->getNumberOfElements();
		n_TRKk = MarlinTrkTracksKAON->getNumberOfElements();
		n_TRKp = MarlinTrkTracksPROTON->getNumberOfElements();
		
		if ( ( n_TRKk != n_TRK ) || ( n_TRKp != n_TRK ) )
		{
			streamlog_out(WARNING) << " Number of tracks in refitted MarlinTrkTrack Collections (pion, proton and kaon) mis-match" << std::endl;
			return 0.0; 
		}
		for (int i_trk = 0 ; i_trk < n_TRK ; ++i_trk )
		{
			Track *protonTrack = dynamic_cast<EVENT::Track*>(MarlinTrkTracksPROTON->getElementAt(i_trk)); 
			Track *kaonTrack = dynamic_cast<EVENT::Track*>(MarlinTrkTracksKAON->getElementAt(i_trk));
			if ( inputTrk == protonTrack )
			{
				trackMass = proton_mass;
			}
			else if ( inputTrk == kaonTrack )
			{
				trackMass = kaon_mass;
			}
			else
			{
				trackMass = pion_mass;
			}
		}
		streamlog_out(DEBUG) << "track mass found: trackMass = " << trackMass << " GeV" << std::endl;
		return trackMass;
	}
	catch(DataNotAvailableException &e)
        {
          streamlog_out(MESSAGE) << "Input Track collections not found in event " << m_nEvt << std::endl;
          return 0.0;
        }

}


void AddFourMomentumCovMatAllPFOs::check(EVENT::LCEvent *pLCEvent)
{

	LCCollection *inputPfoCollection{};
	LCCollection *outputPfoCollection{};
	try
	{
		inputPfoCollection = pLCEvent->getCollection(m_inputPfoCollection);
		outputPfoCollection = pLCEvent->getCollection(m_outputPfoCollection);
		int n_inputPFOs = inputPfoCollection->getNumberOfElements();
		int n_outputPFOs = outputPfoCollection->getNumberOfElements();
		streamlog_out(DEBUG) << " CHECK : processed events: " << m_nEvtSum << " (Number of inputPFOS: " << n_inputPFOs << " , Number of outputPFOs: " << n_outputPFOs <<")" << std::endl;
	}
	catch(DataNotAvailableException &e)
        {
          streamlog_out(MESSAGE) << "Input/Output collection not found in event " << m_nEvt << std::endl;
        }

}

void AddFourMomentumCovMatAllPFOs::end()
{

	m_pTFile->cd();
	m_pTTree->Write();
	h_nTracks_PFOCharge->Write();
	h_nClusters_nTracks->Write();
	h_SigmaPx2nT->Write();
	h_SigmaPxPynT->Write();
	h_SigmaPy2nT->Write();
	h_SigmaPxPznT->Write();
	h_SigmaPyPznT->Write();
	h_SigmaPz2nT->Write();
	h_SigmaPxEnT->Write();
	h_SigmaPyEnT->Write();
	h_SigmaPzEnT->Write();
	h_SigmaE2nT->Write();
	h_SigmaPx2T->Write();
	h_SigmaPxPyT->Write();
	h_SigmaPy2T->Write();
	h_SigmaPxPzT->Write();
	h_SigmaPyPzT->Write();
	h_SigmaPz2T->Write();
	h_SigmaPxET->Write();
	h_SigmaPyET->Write();
	h_SigmaPzET->Write();
	h_SigmaE2T->Write();
	h_SigmaPx2->Write();
	h_SigmaPxPy->Write();
	h_SigmaPy2->Write();
	h_SigmaPxPz->Write();
	h_SigmaPyPz->Write();
	h_SigmaPz2->Write();
	h_SigmaPxE->Write();
	h_SigmaPyE->Write();
	h_SigmaPzE->Write();
	h_SigmaE2->Write();
	m_pTFile->Close();
	delete m_pTFile;

//	std::cout << " END : processed events: " << m_nEvtSum << std::endl;

}
