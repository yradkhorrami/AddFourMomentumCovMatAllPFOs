#include "AddFourMomentumCovMatAllPFOs.h"
#include <iostream>
#include <EVENT/LCCollection.h>
#include "EVENT/LCCollection.h"
#include "IMPL/LCCollectionVec.h"
#include <UTIL/LCRelationNavigator.h>
#include "EVENT/MCParticle.h"
#include "EVENT/Cluster.h"
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
m_Histograms(NULL),
m_CovMatElements(NULL),
m_NeutralPFOswithoutTrak(NULL),
m_NeutralPFOswith2Trak(NULL),
m_ChargedPFOs(NULL),
m_ErrorParameterization(NULL),
m_Photon(NULL),
m_NeutralPFO(NULL),
m_ChargedPFO(NULL),
m_TotalPFOs(NULL),
m_TotalNHs(NULL),
m_TotalRest(NULL),
h_nTracks_PFOCharge(NULL),
h_nClusters_nTracks(NULL),
h_clusterE_pfoE(NULL),
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
h_SigmaE2(NULL),
h_NeutPFO_PDG(NULL),
h_NeutPFO_TYPE(NULL),
h_NeutPFO_IDasPhoton(NULL),
h_NeutPFO_IDasOther(NULL),
h_NeutPFO_Mass(NULL),
h_EP_photons(NULL),
h_EP_NeutralHadrons(NULL),
h_NeutPFO_Weight(NULL),
h_ResidualEnergy_ph(NULL),
h_ResidualTheta_ph(NULL),
h_ResidualPhi_ph(NULL),
h_ErrorEnergy_ph(NULL),
h_ErrorTheta_ph(NULL),
h_ErrorPhi_ph(NULL),
h_NormalizedResidualEnergy_ph(NULL),
h_NormalizedResidualTheta_ph(NULL),
h_NormalizedResidualPhi_ph(NULL),
h_ResidualEnergy_NH(NULL),
h_ResidualTheta_NH(NULL),
h_ResidualPhi_NH(NULL),
h_ErrorEnergy_NH(NULL),
h_ErrorTheta_NH(NULL),
h_ErrorPhi_NH(NULL),
h_NormalizedResidualEnergy_NH(NULL),
h_NormalizedResidualTheta_NH(NULL),
h_NormalizedResidualPhi_NH(NULL),
h_ResidualEnergy_CH(NULL),
h_ResidualTheta_CH(NULL),
h_ResidualPhi_CH(NULL),
h_ErrorEnergy_CH(NULL),
h_ErrorTheta_CH(NULL),
h_ErrorPhi_CH(NULL),
h_NormalizedResidualEnergy_CH(NULL),
h_NormalizedResidualTheta_CH(NULL),
h_NormalizedResidualPhi_CH(NULL),
h_ResidualTotalEnergy_NH(NULL),
h_ResidualTotalTheta_NH(NULL),
h_ResidualTotalPhi_NH(NULL),
h_ErrorTotalEnergy_NH(NULL),
h_ErrorTotalTheta_NH(NULL),
h_ErrorTotalPhi_NH(NULL),
h_NormalizedResidualTotalEnergy_NH(NULL),
h_NormalizedResidualTotalTheta_NH(NULL),
h_NormalizedResidualTotalPhi_NH(NULL),
h_ResidualTotalEnergy_R(NULL),
h_ResidualTotalTheta_R(NULL),
h_ResidualTotalPhi_R(NULL),
h_ErrorTotalEnergy_R(NULL),
h_ErrorTotalTheta_R(NULL),
h_ErrorTotalPhi_R(NULL),
h_NormalizedResidualTotalEnergy_R(NULL),
h_NormalizedResidualTotalTheta_R(NULL),
h_NormalizedResidualTotalPhi_R(NULL),
h_ResidualTotalEnergy(NULL),
h_ResidualTotalTheta(NULL),
h_ResidualTotalPhi(NULL),
h_ErrorTotalEnergy(NULL),
h_ErrorTotalTheta(NULL),
h_ErrorTotalPhi(NULL),
h_NormalizedResidualTotalEnergy(NULL),
h_NormalizedResidualTotalTheta(NULL),
h_NormalizedResidualTotalPhi(NULL),
h_NH_EclusterPlusMass_Emcp(NULL),
h_NHEnergy(NULL)
{
	_description = "Set the convariance matrix in (P,E) for all pfos (charged particles, neutral hadrons and photons)";

	registerInputCollection(	LCIO::RECONSTRUCTEDPARTICLE,
					"inputPfoCollection",
					"Name of input pfo collection",
					m_inputPfoCollection,
					std::string("PandoraPFOs")
				);

	registerInputCollection(	LCIO::LCRELATION,
					"ClusterMCTruthLinkCollection",
					"Name of input m_ClusterMCTruthLink Collection",
					m_ClusterMCTruthLinkCollection,
					std::string("ClusterMCTruthLink")
				);

	registerInputCollection(	LCIO::LCRELATION,
					"MCTruthClusterLinkCollection",
					"Name of input MCTruthClusterLink Collection",
					m_MCTruthClusterLinkCollection,
					std::string("MCTruthClusterLink")
				);

	registerInputCollection(	LCIO::LCRELATION,
					"TrackMCTruthLinkCollection",
					"Name of input TrackMCTruthLink Collection",
					m_TrackMCTruthLinkCollection,
					std::string("MarlinTrkTracksMCTruthLink")
				);

	registerInputCollection(	LCIO::LCRELATION,
					"MCTruthTrackLinkCollection",
					"Name of input MCTruthTrackLink Collection",
					m_MCTruthTrackLinkCollection,
					std::string("MCTruthMarlinTrkTracksLink")
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

	registerProcessorParameter(	"updateNormalNeutrals",
					"true: update CovMat for neutrals without track, false: do not update CovMat for neutrals without track",
					m_updateNormalNeutrals,
					bool(false)
				);

	registerProcessorParameter(	"updateNeutrals_wTrack",
					"true: update CovMat for neutrals with track, false: do not update CovMat for neutrals with track",
					m_updateNeutrals_wTrack,
					bool(false)
				);

	registerProcessorParameter(	"updateCharged",
					"true: update CovMat for charged particles, false: do not update CovMat for charged particles",
					m_updateCharged,
					bool(false)
				);

	registerProcessorParameter(	"AssumeNeutralPFOMassive",
					"true: Neutral PFOs are taken massive, false: Neutral PFOs are taken massless",
					m_AssumeNeutralPFOMassive,
					bool(false)
				);

	registerProcessorParameter(	"isClusterEnergyKinEnergy",
					"true: the cluster energy is interpreted as kinetic energy of PFO, false: the cluster energy is interpreted as momentum magnitude of PFO",
					m_isClusterEnergyKinEnergy,
					bool(false)
				);

	registerProcessorParameter(	"useClusterPositionError",
					"true: use cluster position error for CovMat, false: use cluster direction error for CovMat",
					m_useClusterPositionError,
					bool(true)
				);

	registerProcessorParameter(	"updatePFO4Momentum",
					"true: Update 4-momentum of PFOs, false: set 4-momentum for PFOs same as input PFO",
					m_updatePFO4Momentum,
					bool(false)
				);

	registerProcessorParameter(	"RootFile",
	                                "Name of the output root file",
					m_rootFile,
					std::string("FourMomentumCovMatAllPFOs.root")
				);

}

void AddFourMomentumCovMatAllPFOs::init()
{

	streamlog_out(MESSAGE) << "   init called  " << std::endl;
	m_Bfield = MarlinUtil::getBzAtOrigin();
	printParameters();
	streamlog_out(MESSAGE) << " BField =  "<< m_Bfield << " Tesla" << std::endl ;
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
	m_Histograms = m_pTFile->mkdir("Histograms");
	m_CovMatElements = m_Histograms->mkdir("CovMatElements");
	m_NeutralPFOswithoutTrak = m_CovMatElements->mkdir("NeutralPFOswithoutTrak");
	m_NeutralPFOswith2Trak = m_CovMatElements->mkdir("NeutralPFOswith2Trak");
	m_ChargedPFOs = m_CovMatElements->mkdir("ChargedPFOs");
	m_ErrorParameterization = m_Histograms->mkdir("ErrorParameterization");
	m_Photon = m_ErrorParameterization->mkdir("Photon");
	m_NeutralPFO = m_ErrorParameterization->mkdir("NeutralPFO");
	m_ChargedPFO = m_ErrorParameterization->mkdir("ChargedPFO");
	m_TotalPFOs = m_ErrorParameterization->mkdir("TotalPFOs");
	m_TotalNHs = m_ErrorParameterization->mkdir("TotalNHs");
	m_TotalRest = m_ErrorParameterization->mkdir("TotalRest");
	h_nTracks_PFOCharge = new TH2I("All PFOs", "; charge of PFO; n_{Tracks}", 5, -2.5, 2.5, 5, -0.5, 4.5);
	h_nClusters_nTracks = new TH2I("Neutral PFOs", "; n_{Clusters}; n_{Tracks}", 5, -0.5, 4.5, 5, -0.5, 4.5);
	h_clusterE_pfoE = new TH2F("Neutral PFOs (n_{Tracks} = 0)", "; E_{Cluster} [GeV]; E_{PFO} [GeV]", 10000, 0.0, 100., 10000, 0.0, 100.);
	h_SigmaPx2nT = new TH2F("Neutral PFOs (n_{Tracks} = 0)", "; #sigma_{p_{x}}^{2} (new PFO) [GeV^{2}]; #sigma_{p_{x}}^{2} (old PFO) [GeV^{2}]", 400, 0.0, 0.1, 400, 0.0, 0.1);
	h_SigmaPxPynT = new TH2F("Neutral PFOs (n_{Tracks} = 0)", "; #sigma_{p_{x}p_{y}} (new PFO) [GeV^{2}]; #sigma_{p_{x}p_{y}} (old PFO) [GeV^{2}]", 400, -0.1, 0.1, 400, -0.1, 0.1);
	h_SigmaPy2nT = new TH2F("Neutral PFOs (n_{Tracks} = 0)", "; #sigma_{p_{y}}^{2} (new PFO) [GeV^{2}]; #sigma_{p_{y}}^{2} (old PFO) [GeV^{2}]", 400, 0.0, 0.1, 400, 0.0, 0.1);
	h_SigmaPxPznT = new TH2F("Neutral PFOs (n_{Tracks} = 0)", "; #sigma_{p_{x}p_{z}} (new PFO) [GeV^{2}]; #sigma_{p_{x}p_{z}} (old PFO) [GeV^{2}]", 400, -0.1, 0.1, 400, -0.1, 0.1);
	h_SigmaPyPznT = new TH2F("Neutral PFOs (n_{Tracks} = 0)", "; #sigma_{p_{y}p_{z}} (new PFO) [GeV^{2}]; #sigma_{p_{y}p_{z}} (old PFO) [GeV^{2}]", 400, -0.1, 0.1, 400, -0.1, 0.1);
	h_SigmaPz2nT = new TH2F("Neutral PFOs (n_{Tracks} = 0)", "; #sigma_{p_{z}}^{2} (new PFO) [GeV^{2}]; #sigma_{p_{z}}^{2} (old PFO) [GeV^{2}]", 400, 0.0, 0.1, 400, 0.0, 0.1);
	h_SigmaPxEnT = new TH2F("Neutral PFOs (n_{Tracks} = 0)", "; #sigma_{p_{x}E} (new PFO) [GeV^{2}]; #sigma_{p_{x}E} (old PFO) [GeV^{2}]", 400, -0.1, 0.1, 400, -0.1, 0.1);
	h_SigmaPyEnT = new TH2F("Neutral PFOs (n_{Tracks} = 0)", "; #sigma_{p_{y}E} (new PFO) [GeV^{2}]; #sigma_{p_{y}E} (old PFO) [GeV^{2}]", 400, -0.1, 0.1, 400, -0.1, 0.1);
	h_SigmaPzEnT = new TH2F("Neutral PFOs (n_{Tracks} = 0)", "; #sigma_{p_{z}E} (new PFO) [GeV^{2}]; #sigma_{p_{z}E} (old PFO) [GeV^{2}]", 400, -0.1, 0.1, 400, -0.1, 0.1);
	h_SigmaE2nT = new TH2F("Neutral PFOs (n_{Tracks} = 0)", "; #sigma_{E}^{2} (new PFO) [GeV^{2}]; #sigma_{E}^{2} (old PFO) [GeV^{2}]", 400, 0.0, 0.6, 400, 0.0, 0.6);
	h_SigmaPx2T = new TH2F("Neutral PFOs (n_{Tracks} = 2)", "; #sigma_{p_{x}}^{2} (new PFO) [GeV^{2}]; #sigma_{p_{x}}^{2} (old PFO) [GeV^{2}]", 400, 0.0, 0.1, 400, 0.0, 0.1);
	h_SigmaPxPyT = new TH2F("Neutral PFOs (n_{Tracks} = 2)", "; #sigma_{p_{x}p_{y}} (new PFO) [GeV^{2}]; #sigma_{p_{x}p_{y}} (old PFO) [GeV^{2}]", 400, -0.1, 0.1, 400, -0.1, 0.1);
	h_SigmaPy2T = new TH2F("Neutral PFOs (n_{Tracks} = 2)", "; #sigma_{p_{y}}^{2} (new PFO) [GeV^{2}]; #sigma_{p_{y}}^{2} (old PFO) [GeV^{2}]", 400, 0.0, 0.1, 400, 0.0, 0.1);
	h_SigmaPxPzT = new TH2F("Neutral PFOs (n_{Tracks} = 2)", "; #sigma_{p_{x}p_{z}} (new PFO) [GeV^{2}]; #sigma_{p_{x}p_{z}} (old PFO) [GeV^{2}]", 400, -0.1, 0.1, 400, -0.1, 0.1);
	h_SigmaPyPzT = new TH2F("Neutral PFOs (n_{Tracks} = 2)", "; #sigma_{p_{y}p_{z}} (new PFO) [GeV^{2}]; #sigma_{p_{y}p_{z}} (old PFO) [GeV^{2}]", 400, -0.1, 0.1, 400, -0.1, 0.1);
	h_SigmaPz2T = new TH2F("Neutral PFOs (n_{Tracks} = 2)", "; #sigma_{p_{z}}^{2} (new PFO) [GeV^{2}]; #sigma_{p_{z}}^{2} (old PFO) [GeV^{2}]", 400, 0.0, 0.1, 400, 0.0, 0.1);
	h_SigmaPxET = new TH2F("Neutral PFOs (n_{Tracks} = 2)", "; #sigma_{p_{x}E} (new PFO) [GeV^{2}]; #sigma_{p_{x}E} (old PFO) [GeV^{2}]", 400, -0.1, 0.1, 400, -0.1, 0.1);
	h_SigmaPyET = new TH2F("Neutral PFOs (n_{Tracks} = 2)", "; #sigma_{p_{y}E} (new PFO) [GeV^{2}]; #sigma_{p_{y}E} (old PFO) [GeV^{2}]", 400, -0.1, 0.1, 400, -0.1, 0.1);
	h_SigmaPzET = new TH2F("Neutral PFOs (n_{Tracks} = 2)", "; #sigma_{p_{z}E} (new PFO) [GeV^{2}]; #sigma_{p_{z}E} (old PFO) [GeV^{2}]", 400, -0.1, 0.1, 400, -0.1, 0.1);
	h_SigmaE2T = new TH2F("Neutral PFOs (n_{Tracks} = 2)", "; #sigma_{E}^{2} (new PFO) [GeV^{2}]; #sigma_{E}^{2} (old PFO) [GeV^{2}]", 400, 0.0, 0.6, 400, 0.0, 0.6);
	h_SigmaPx2 = new TH2F("Charged PFOs", "; #sigma_{p_{x}}^{2} (new PFO) [GeV^{2}]; #sigma_{p_{x}}^{2} (old PFO) [GeV^{2}]", 400, 0.0, 0.1, 400, 0.0, 0.1);
	h_SigmaPxPy = new TH2F("Charged PFOs", "; #sigma_{p_{x}p_{y}} (new PFO) [GeV^{2}]; #sigma_{p_{x}p_{y}} (old PFO) [GeV^{2}]", 400, -0.1, 0.1, 400, -0.1, 0.1);
	h_SigmaPy2 = new TH2F("Charged PFOs", "; #sigma_{p_{y}}^{2} (new PFO) [GeV^{2}]; #sigma_{p_{y}}^{2} (old PFO) [GeV^{2}]", 400, 0.0, 0.1, 400, 0.0, 0.1);
	h_SigmaPxPz = new TH2F("Charged PFOs", "; #sigma_{p_{x}p_{z}} (new PFO) [GeV^{2}]; #sigma_{p_{x}p_{z}} (old PFO) [GeV^{2}]", 400, -0.1, 0.1, 400, -0.1, 0.1);
	h_SigmaPyPz = new TH2F("Charged PFOs", "; #sigma_{p_{y}p_{z}} (new PFO) [GeV^{2}]; #sigma_{p_{y}p_{z}} (old PFO) [GeV^{2}]", 400, -0.1, 0.1, 400, -0.1, 0.1);
	h_SigmaPz2 = new TH2F("Charged PFOs", "; #sigma_{p_{z}}^{2} (new PFO) [GeV^{2}]; #sigma_{p_{z}}^{2} (old PFO) [GeV^{2}]", 400, 0.0, 0.1, 400, 0.0, 0.1);
	h_SigmaPxE = new TH2F("Charged PFOs", "; #sigma_{p_{x}E} (new PFO) [GeV^{2}]; #sigma_{p_{x}E} (old PFO) [GeV^{2}]", 400, -0.1, 0.1, 400, -0.1, 0.1);
	h_SigmaPyE = new TH2F("Charged PFOs", "; #sigma_{p_{y}E} (new PFO) [GeV^{2}]; #sigma_{p_{y}E} (old PFO) [GeV^{2}]", 400, -0.1, 0.1, 400, -0.1, 0.1);
	h_SigmaPzE = new TH2F("Charged PFOs", "; #sigma_{p_{z}E} (new PFO) [GeV^{2}]; #sigma_{p_{z}E} (old PFO) [GeV^{2}]", 400, -0.1, 0.1, 400, -0.1, 0.1);
	h_SigmaE2 = new TH2F("Charged PFOs", "; #sigma_{E}^{2} (new PFO) [GeV^{2}]; #sigma_{E}^{2} (old PFO) [GeV^{2}]", 400, 0.0, 0.6, 400, 0.0, 0.6);
	h_NeutPFO_PDG = new TH1I("Neutral PFOs PDG", "; PDG Code", 200001, -100000.5, 100000.5);
	h_NeutPFO_TYPE = new TH1I("Neutral PFOs TYPE", "; True Part. Type", 15, 0, 15);
	h_NeutPFO_TYPE->GetXaxis()->SetBinLabel(1,"e^{#pm}");
	h_NeutPFO_TYPE->GetXaxis()->SetBinLabel(2,"#mu^{#pm}");
	h_NeutPFO_TYPE->GetXaxis()->SetBinLabel(3,"#gamma");
	h_NeutPFO_TYPE->GetXaxis()->SetBinLabel(4,"K^{0}_{L}");
	h_NeutPFO_TYPE->GetXaxis()->SetBinLabel(5,"#pi^{#pm}");
	h_NeutPFO_TYPE->GetXaxis()->SetBinLabel(6,"K^{0}_{S}");
	h_NeutPFO_TYPE->GetXaxis()->SetBinLabel(7,"K^{#pm}");
	h_NeutPFO_TYPE->GetXaxis()->SetBinLabel(8,"n");
	h_NeutPFO_TYPE->GetXaxis()->SetBinLabel(9,"p");
	h_NeutPFO_TYPE->GetXaxis()->SetBinLabel(10,"#Sigma^{-}");
	h_NeutPFO_TYPE->GetXaxis()->SetBinLabel(11,"#Lambda");
	h_NeutPFO_TYPE->GetXaxis()->SetBinLabel(12,"#Sigma^{+}");
	h_NeutPFO_TYPE->GetXaxis()->SetBinLabel(13,"#Xi^{-}");
	h_NeutPFO_TYPE->GetXaxis()->SetBinLabel(14,"#Xi");
	h_NeutPFO_TYPE->GetXaxis()->SetBinLabel(15,"Others");
	h_NeutPFO_IDasPhoton = new TH1I("Photons", "; True Part. Type", 15, 0, 15);
	h_NeutPFO_IDasPhoton->GetXaxis()->SetBinLabel(1,"e^{#pm}");
	h_NeutPFO_IDasPhoton->GetXaxis()->SetBinLabel(2,"#mu^{#pm}");
	h_NeutPFO_IDasPhoton->GetXaxis()->SetBinLabel(3,"#gamma");
	h_NeutPFO_IDasPhoton->GetXaxis()->SetBinLabel(4,"K^{0}_{L}");
	h_NeutPFO_IDasPhoton->GetXaxis()->SetBinLabel(5,"#pi^{#pm}");
	h_NeutPFO_IDasPhoton->GetXaxis()->SetBinLabel(6,"K^{0}_{S}");
	h_NeutPFO_IDasPhoton->GetXaxis()->SetBinLabel(7,"K^{#pm}");
	h_NeutPFO_IDasPhoton->GetXaxis()->SetBinLabel(8,"n");
	h_NeutPFO_IDasPhoton->GetXaxis()->SetBinLabel(9,"p");
	h_NeutPFO_IDasPhoton->GetXaxis()->SetBinLabel(10,"#Sigma^{-}");
	h_NeutPFO_IDasPhoton->GetXaxis()->SetBinLabel(11,"#Lambda");
	h_NeutPFO_IDasPhoton->GetXaxis()->SetBinLabel(12,"#Sigma^{+}");
	h_NeutPFO_IDasPhoton->GetXaxis()->SetBinLabel(13,"#Xi^{-}");
	h_NeutPFO_IDasPhoton->GetXaxis()->SetBinLabel(14,"#Xi");
	h_NeutPFO_IDasPhoton->GetXaxis()->SetBinLabel(15,"Others");
	h_NeutPFO_IDasOther = new TH1I("Other Neutal PFOs", "; True Part. Type", 15, 0, 15);
	h_NeutPFO_IDasOther->GetXaxis()->SetBinLabel(1,"e^{#pm}");
	h_NeutPFO_IDasOther->GetXaxis()->SetBinLabel(2,"#mu^{#pm}");
	h_NeutPFO_IDasOther->GetXaxis()->SetBinLabel(3,"#gamma");
	h_NeutPFO_IDasOther->GetXaxis()->SetBinLabel(4,"K^{0}_{L}");
	h_NeutPFO_IDasOther->GetXaxis()->SetBinLabel(5,"#pi^{#pm}");
	h_NeutPFO_IDasOther->GetXaxis()->SetBinLabel(6,"K^{0}_{S}");
	h_NeutPFO_IDasOther->GetXaxis()->SetBinLabel(7,"K^{#pm}");
	h_NeutPFO_IDasOther->GetXaxis()->SetBinLabel(8,"n");
	h_NeutPFO_IDasOther->GetXaxis()->SetBinLabel(9,"p");
	h_NeutPFO_IDasOther->GetXaxis()->SetBinLabel(10,"#Sigma^{-}");
	h_NeutPFO_IDasOther->GetXaxis()->SetBinLabel(11,"#Lambda");
	h_NeutPFO_IDasOther->GetXaxis()->SetBinLabel(12,"#Sigma^{+}");
	h_NeutPFO_IDasOther->GetXaxis()->SetBinLabel(13,"#Xi^{-}");
	h_NeutPFO_IDasOther->GetXaxis()->SetBinLabel(14,"#Xi");
	h_NeutPFO_IDasOther->GetXaxis()->SetBinLabel(15,"Others");
	h_NeutPFO_Mass = new TH1F("Neutral PFOs Mass", "; PFO Mass [GeV]", 200, 0.0, 10.0);
	h_EP_photons = new TH2F("Photons", "; |#vec{p}_{PFO}|^{2} [GeV^{2}]; E_{PFO}^{2} [GeV^{2}]", 1000, 0.0, 10.0, 1000, 0.0, 10.0);
	h_EP_NeutralHadrons = new TH2F("Neutral Hadrons", "; |#vec{p}_{PFO}|^{2} [GeV^{2}]; E_{PFO}^{2} [GeV^{2}]", 1000, 0.0, 10.0, 1000, 0.0, 10.0);
	h_NeutPFO_Weight = new TH1F("Neutral Hadrons MCP Link Weight", "; Link weight", 100, 0.0, 1.0);
	h_ResidualEnergy_ph = new TH1F("Photons", "; E_{REC} - E_{MCP} [GeV]", 200, -10.0, 10.0);
	h_ResidualTheta_ph = new TH1F("Photons", "; #theta_{REC} - #theta_{MCP} [radian]", 6800, -3.4, 3.4);
	h_ResidualPhi_ph = new TH1F("Photons", "; #phi_{REC} - #phi_{MCP} [radian]", 6800, -3.4, 3.4);
	h_ErrorEnergy_ph = new TH1F("Photons", "; #sigma_{E} [GeV]", 1000, 0.0, 10.0);
	h_ErrorTheta_ph = new TH1F("Photons", "; #sigma_{#theta} [radian]", 10000, 0.0, 1.0);
	h_ErrorPhi_ph = new TH1F("Photons", "; #sigma_{#phi} [radian]", 10000, 0.0, 1.0);
	h_NormalizedResidualEnergy_ph = new TH1F("Photons", "; (E_{REC} - E_{MCP}) / #sigma_{E}", 200, -10.0, 10.0);
	h_NormalizedResidualTheta_ph = new TH1F("Photons", "; (#theta_{REC} - #theta_{MCP}) / #sigma_{#theta}", 200, -10.0, 10.0);
	h_NormalizedResidualPhi_ph = new TH1F("Photons", "; (#phi_{REC} - #phi_{MCP}) / #sigma_{#phi}", 200, -10.0, 10.0);
	h_ResidualEnergy_NH = new TH1F("Neutral Hadrons", "; E_{REC} - E_{MCP} [GeV]", 200, -10.0, 10.0);
	h_ResidualTheta_NH = new TH1F("Neutral Hadrons", "; #theta_{REC} - #theta_{MCP} [radian]", 6800, -3.4, 3.4);
	h_ResidualPhi_NH = new TH1F("Neutral Hadrons", "; #phi_{REC} - #phi_{MCP} [radian]", 6800, -3.4, 3.4);
	h_ErrorEnergy_NH = new TH1F("Neutral Hadrons", "; #sigma_{E} [GeV]", 1000, 0.0, 10.0);
	h_ErrorTheta_NH = new TH1F("Neutral Hadrons", "; #sigma_{#theta} [radian]", 10000, 0.0, 1.0);
	h_ErrorPhi_NH = new TH1F("Neutral Hadrons", "; #sigma_{#phi} [radian]", 10000, 0.0, 1.0);
	h_NormalizedResidualEnergy_NH = new TH1F("Neutral Hadrons", "; (E_{REC} - E_{MCP}) / #sigma_{E}", 200, -10.0, 10.0);
	h_NormalizedResidualTheta_NH = new TH1F("Neutral Hadrons", "; (#theta_{REC} - #theta_{MCP}) / #sigma_{#theta}", 200, -10.0, 10.0);
	h_NormalizedResidualPhi_NH = new TH1F("Neutral Hadrons", "; (#phi_{REC} - #phi_{MCP}) / #sigma_{#phi}", 200, -10.0, 10.0);
	h_ResidualEnergy_CH = new TH1F("Charged PFOs", "; E_{REC} - E_{MCP} [GeV]", 200, -10.0, 10.0);
	h_ResidualTheta_CH = new TH1F("Charged PFOs", "; #theta_{REC} - #theta_{MCP} [radian]", 6800, -3.4, 3.4);
	h_ResidualPhi_CH = new TH1F("Charged PFOs", "; #phi_{REC} - #phi_{MCP} [radian]", 6800, -3.4, 3.4);
	h_ErrorEnergy_CH = new TH1F("Charged PFOs", "; #sigma_{E} [GeV]", 1000, 0.0, 10.0);
	h_ErrorTheta_CH = new TH1F("Charged PFOs", "; #sigma_{#theta} [radian]", 10000, 0.0, 1.0);
	h_ErrorPhi_CH = new TH1F("Charged PFOs", "; #sigma_{#phi} [radian]", 10000, 0.0, 1.0);
	h_NormalizedResidualEnergy_CH = new TH1F("Charged PFOs", "; (E_{REC} - E_{MCP}) / #sigma_{E}", 200, -10.0, 10.0);
	h_NormalizedResidualTheta_CH = new TH1F("Charged PFOs", "; (#theta_{REC} - #theta_{MCP}) / #sigma_{#theta}", 200, -10.0, 10.0);
	h_NormalizedResidualPhi_CH = new TH1F("Charged PFOs", "; (#phi_{REC} - #phi_{MCP}) / #sigma_{#phi}", 200, -10.0, 10.0);
	h_ResidualTotalEnergy_NH = new TH1F("Neutral Hadrons", "; E_{REC} - E_{MCP} [GeV]", 200, -10.0, 10.0);
	h_ResidualTotalTheta_NH = new TH1F("Neutral Hadrons", "; #theta_{REC} - #theta_{MCP} [radian]", 6800, -3.4, 3.4);
	h_ResidualTotalPhi_NH = new TH1F("Neutral Hadrons", "; #phi_{REC} - #phi_{MCP} [radian]", 6800, -3.4, 3.4);
	h_ErrorTotalEnergy_NH = new TH1F("Neutral Hadrons", "; #sigma_{E} [GeV]", 1000, 0.0, 10.0);
	h_ErrorTotalTheta_NH = new TH1F("Neutral Hadrons", "; #sigma_{#theta} [radian]", 10000, 0.0, 1.0);
	h_ErrorTotalPhi_NH = new TH1F("Neutral Hadrons", "; #sigma_{#phi} [radian]", 10000, 0.0, 1.0);
	h_NormalizedResidualTotalEnergy_NH = new TH1F("Neutral Hadrons", "; (E_{REC} - E_{MCP}) / #sigma_{E}", 200, -10.0, 10.0);
	h_NormalizedResidualTotalTheta_NH = new TH1F("Neutral Hadrons", "; (#theta_{REC} - #theta_{MCP}) / #sigma_{#theta}", 200, -10.0, 10.0);
	h_NormalizedResidualTotalPhi_NH = new TH1F("Neutral Hadrons", "; (#phi_{REC} - #phi_{MCP}) / #sigma_{#phi}", 200, -10.0, 10.0);
	h_ResidualTotalEnergy_R = new TH1F("photon + charged PFOs", "; E_{REC} - E_{MCP} [GeV]", 200, -10.0, 10.0);
	h_ResidualTotalTheta_R = new TH1F("photon + charged PFOs", "; #theta_{REC} - #theta_{MCP} [radian]", 6800, -3.4, 3.4);
	h_ResidualTotalPhi_R = new TH1F("photon + charged PFOs", "; #phi_{REC} - #phi_{MCP} [radian]", 6800, -3.4, 3.4);
	h_ErrorTotalEnergy_R = new TH1F("photon + charged PFOs", "; #sigma_{E} [GeV]", 1000, 0.0, 10.0);
	h_ErrorTotalTheta_R = new TH1F("photon + charged PFOs", "; #sigma_{#theta} [radian]", 10000, 0.0, 1.0);
	h_ErrorTotalPhi_R = new TH1F("photon + charged PFOs", "; #sigma_{#phi} [radian]", 10000, 0.0, 1.0);
	h_NormalizedResidualTotalEnergy_R = new TH1F("photon + charged PFOs", "; (E_{REC} - E_{MCP}) / #sigma_{E}", 200, -10.0, 10.0);
	h_NormalizedResidualTotalTheta_R = new TH1F("photon + charged PFOs", "; (#theta_{REC} - #theta_{MCP}) / #sigma_{#theta}", 200, -10.0, 10.0);
	h_NormalizedResidualTotalPhi_R = new TH1F("photon + charged PFOs", "; (#phi_{REC} - #phi_{MCP}) / #sigma_{#phi}", 200, -10.0, 10.0);
	h_ResidualTotalEnergy = new TH1F("all PFOs", "; E_{REC} - E_{MCP} [GeV]", 200, -10.0, 10.0);
	h_ResidualTotalTheta = new TH1F("all PFOs", "; #theta_{REC} - #theta_{MCP} [radian]", 6800, -3.4, 3.4);
	h_ResidualTotalPhi = new TH1F("all PFOs", "; #phi_{REC} - #phi_{MCP} [radian]", 6800, -3.4, 3.4);
	h_ErrorTotalEnergy = new TH1F("all PFOs", "; #sigma_{E} [GeV]", 1000, 0.0, 10.0);
	h_ErrorTotalTheta = new TH1F("all PFOs", "; #sigma_{#theta} [radian]", 10000, 0.0, 1.0);
	h_ErrorTotalPhi = new TH1F("all PFOs", "; #sigma_{#phi} [radian]", 10000, 0.0, 1.0);
	h_NormalizedResidualTotalEnergy = new TH1F("all PFOs", "; (E_{REC} - E_{MCP}) / #sigma_{E}", 200, -10.0, 10.0);
	h_NormalizedResidualTotalTheta = new TH1F("all PFOs", "; (#theta_{REC} - #theta_{MCP}) / #sigma_{#theta}", 200, -10.0, 10.0);
	h_NormalizedResidualTotalPhi = new TH1F("all PFOs", "; (#phi_{REC} - #phi_{MCP}) / #sigma_{#phi}", 200, -10.0, 10.0);
	h_NH_EclusterPlusMass_Emcp = new TH2F("Neutral Hadrons", "; E_{MCP} [GeV]; E_{cluster} + m [GeV]", 1000, 0.0, 10.0, 1000, 0.0, 10.0);
	h_NHEnergy = new TH2F("Neutral Hadrons", "; E_{MCP} [GeV]; E_{PFO} [GeV]", 1000, 0.0, 10.0, 1000, 0.0, 10.0);

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
	streamlog_out(MESSAGE) << "processed event 	" << m_nEvtSum << std::endl;

	LCCollection *inputPfoCollection{};
	this->Clear();

	try
	{
		inputPfoCollection = pLCEvent->getCollection(m_inputPfoCollection);
		int n_PFO = inputPfoCollection->getNumberOfElements();
		LCCollectionVec *m_col_outputPfo = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
		TLorentzVector TotalPFOFourMomentum(0.,0.,0.,0.);
		TLorentzVector TotalMCPFourMomentum(0.,0.,0.,0.);
		std::vector<float> TotaloutputCovMatrix(10,0.);
		TLorentzVector TotalPFOFourMomentum_NH(0.,0.,0.,0.);
		TLorentzVector TotalMCPFourMomentum_NH(0.,0.,0.,0.);
		std::vector<float> TotaloutputCovMatrix_NH(10,0.);
		TLorentzVector TotalPFOFourMomentum_R(0.,0.,0.,0.);
		TLorentzVector TotalMCPFourMomentum_R(0.,0.,0.,0.);
		std::vector<float> TotaloutputCovMatrix_R(10,0.);
		std::vector<float> TotalPFOResidual( 3 , 0.0 );
		std::vector<float> TotalPFOCovMatPolar( 10 , 0.0 );
		std::vector<float> TotalPFOResidual_NH( 3 , 0.0 );
		std::vector<float> TotalPFOCovMatPolar_NH( 10 , 0.0 );
		std::vector<float> TotalPFOResidual_R( 3 , 0.0 );
		std::vector<float> TotalPFOCovMatPolar_R( 10 , 0.0 );
		
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
			std::vector<float> inputCovMatrix( 10, 0.0 );
			inputCovMatrix = inputPFO->getCovMatrix();
			outputCovMatrix = inputCovMatrix;
			TLorentzVector mcpFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
			std::vector<float> PFOResidual( 3 , 0.0 );
			std::vector<float> PFOCovMatPolar( 10 , 0.0 );
			std::vector<float> PFOCoordinateError( 6 , 0.0 );
			mcpFourMomentum = this->getLinkedMCP( pLCEvent , inputPFO, nTrackspfo , nClusterspfo );
			if ( mcpFourMomentum.E() != 0 ) TotalMCPFourMomentum += mcpFourMomentum;
			if ( pfoCharge == 0)
			{
				h_nClusters_nTracks->Fill( nClusterspfo , nTrackspfo );
				if ( nTrackspfo == 0 )
				{

					
					streamlog_out(DEBUG) << "PFO is neutral without track, CovMatrix is set using cluster information" << std::endl;
					float pfoP2 = pow( inputPFO->getMomentum()[0] , 2 ) + pow( inputPFO->getMomentum()[1] , 2 ) + pow( inputPFO->getMomentum()[2] , 2 );
					float pfoE2 = pow ( inputPFO->getEnergy() , 2 );
					h_NeutPFO_Mass->Fill( inputPFO->getMass() );

					if ( !m_AssumeNeutralPFOMassive ) pfoMass = 0.0;
					if ( inputPFO->getType() == 22 && mcpFourMomentum.E() != 0 )
					{
						h_EP_photons->Fill( pfoP2 , pfoE2 );
						TotalMCPFourMomentum_R += mcpFourMomentum;
					}
					else if ( mcpFourMomentum.E() != 0 )
					{
						h_EP_NeutralHadrons->Fill( pfoP2 , pfoE2 );
						TotalMCPFourMomentum_NH += mcpFourMomentum;
					}
					float clusterEnergy	= ( inputPFO->getClusters()[0] )->getEnergy();
					float clusterX		= ( inputPFO->getClusters()[0] )->getPosition()[0];
					float clusterY		= ( inputPFO->getClusters()[0] )->getPosition()[1];
					float clusterZ		= ( inputPFO->getClusters()[0] )->getPosition()[2];
					clusterPosition		= TVector3( clusterX , clusterY , clusterZ );
					float clusterDistance	= sqrt( pow( clusterX , 2 ) + pow( clusterY , 2 ) + pow( clusterZ , 2 ) );
					float pfoMomentumMag	= 0;
					float pfoEnergy	= inputPFO->getEnergy();
					float pfoE;
					if ( m_isClusterEnergyKinEnergy )
					{
						pfoMomentumMag	= sqrt( pow( pfoEnergy , 2 ) + 2 * pfoMass * pfoEnergy );
					}
					else
					{
						pfoMomentumMag	= pfoEnergy;
					}
					float pfoPx;
					float pfoPy;
					float pfoPz;
					if ( m_updatePFO4Momentum )
					{
						pfoPx	= pfoMomentumMag * clusterX / clusterDistance;
						pfoPy	= pfoMomentumMag * clusterY / clusterDistance;
						pfoPz	= pfoMomentumMag * clusterZ / clusterDistance;
						if ( m_isClusterEnergyKinEnergy )
						{
							pfoE	= pfoEnergy + pfoMass;
						}
						else
						{
							pfoE	= sqrt( pow( pfoMomentumMag , 2 ) + pow( pfoMass , 2 ) );
						}
					}
					else
					{
						pfoPx	= inputPFO->getMomentum()[ 0 ];
						pfoPy	= inputPFO->getMomentum()[ 1 ];
						pfoPz	= inputPFO->getMomentum()[ 2 ];
						pfoE	= inputPFO->getEnergy();
					}

					std::vector<float> clusterDirectionError = ( inputPFO->getClusters()[0] )->getDirectionError();
					std::vector<float> clusterPositionError = ( inputPFO->getClusters()[0] )->getPositionError();
					float clusterEnergyError= ( inputPFO->getClusters()[0] )->getEnergyError();
					streamlog_out(DEBUG) << "cluster / PFO Energy : " << clusterEnergy << " / " << pfoE << std::endl;
					h_clusterE_pfoE->Fill( clusterEnergy , pfoE );
					TVector3 pfoMomentum( pfoPx , pfoPy , pfoPz );
					pfoFourMomentum	= TLorentzVector( pfoMomentum , pfoE );
					if ( mcpFourMomentum.E() != 0 ) TotalPFOFourMomentum += pfoFourMomentum;
					if ( inputPFO->getType() == 22 && mcpFourMomentum.E() != 0 )
					{
						TotalPFOFourMomentum_R += pfoFourMomentum;
					}
					else if ( mcpFourMomentum.E() != 0 )
					{
						TotalPFOFourMomentum_NH += pfoFourMomentum;
					}

					
//					outputCovMatrix	= this->UpdateNeutralPFOCovMat( clusterPosition , clusterEnergy , pfoMass , clusterPositionError , clusterEnergyError );
					if ( m_updateNormalNeutrals ) outputCovMatrix	= this->UpdateNeutralPFOCovMat( clusterPosition , pfoEnergy , pfoMass , clusterPositionError , clusterEnergyError );
					if ( mcpFourMomentum.E() != 0 )
					{
						for ( int i = 0 ; i < 10 ; ++i )
						{
							if (inputPFO->getType() == 22)
							{
								TotaloutputCovMatrix_R[ i ] += outputCovMatrix[ i ];
							}
							else
							{
								TotaloutputCovMatrix_NH[ i ] += outputCovMatrix[ i ];
							}
						}
						PFOResidual = this->getPFOResidual( pfoFourMomentum , mcpFourMomentum );
						PFOCovMatPolar = this->getPFOCovMatPolarCoordinate( pfoFourMomentum , outputCovMatrix );
						PFOCoordinateError = this->getClusterDirectionError( clusterPosition , clusterPositionError );
						float EnergyError;
						float ThetaError;
						float PhiError;
						if ( m_useClusterPositionError )
						{
							ThetaError = sqrt ( PFOCovMatPolar[ 0 ] );
							PhiError = sqrt ( PFOCovMatPolar[ 2 ] );
							EnergyError = sqrt ( PFOCovMatPolar[ 9 ] );
						}
						else
						{
							ThetaError = sqrt ( PFOCoordinateError[ 2 ] );
							PhiError = sqrt ( PFOCoordinateError[ 5 ] );
							EnergyError = clusterEnergyError;
						}
						if ( inputPFO->getType() == 22 )
						{
							h_ResidualEnergy_ph->Fill( PFOResidual[ 0 ] );
							h_ResidualTheta_ph->Fill( PFOResidual[ 1 ] );
							h_ResidualPhi_ph->Fill( PFOResidual[ 2 ] );
							h_ErrorEnergy_ph->Fill( EnergyError );
							h_ErrorTheta_ph->Fill( ThetaError );
							h_ErrorPhi_ph->Fill( PhiError );
							h_NormalizedResidualEnergy_ph->Fill( PFOResidual[ 0 ] / EnergyError );
							h_NormalizedResidualTheta_ph->Fill( PFOResidual[ 1 ] / ThetaError );
							h_NormalizedResidualPhi_ph->Fill( PFOResidual[ 2 ] / PhiError );
						}
						else
						{
							h_ResidualEnergy_NH->Fill( PFOResidual[ 0 ] );
							h_ResidualTheta_NH->Fill( PFOResidual[ 1 ] );
							h_ResidualPhi_NH->Fill( PFOResidual[ 2 ] );
							h_ErrorEnergy_NH->Fill( EnergyError );
							h_ErrorTheta_NH->Fill( ThetaError );
							h_ErrorPhi_NH->Fill( PhiError );
							h_NormalizedResidualEnergy_NH->Fill( PFOResidual[ 0 ] / EnergyError );
							h_NormalizedResidualTheta_NH->Fill( PFOResidual[ 1 ] / ThetaError );
							h_NormalizedResidualPhi_NH->Fill( PFOResidual[ 2 ] / PhiError );
						}
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
					if ( mcpFourMomentum.E() != 0 ) TotalMCPFourMomentum_R += mcpFourMomentum;
					streamlog_out(DEBUG) << "PFO is neutral (with two tracks), CovMatrix is set using track information" << std::endl;
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
							if ( m_updateNeutrals_wTrack ) outputCovMatrix[i] += pfoTrackCovMat[i];
						}
						if ( m_updatePFO4Momentum )
						{
							pfoFourMomentum		+= trackFourMomentum;
						}
						else
						{
							pfoFourMomentum	= TLorentzVector( inputPFO->getMomentum()[0] , inputPFO->getMomentum()[1] , inputPFO->getMomentum()[2] , inputPFO->getEnergy() );
						}
					}
					if ( mcpFourMomentum.E() != 0 )
					{
						TotalPFOFourMomentum += pfoFourMomentum;
						TotalPFOFourMomentum_R += pfoFourMomentum;
						for ( int i = 0 ; i < 10 ; ++i )
						{
							TotaloutputCovMatrix_R[ i ] += outputCovMatrix[ i ];
						}
					}
					if ( mcpFourMomentum.E() != 0 && m_updateNeutrals_wTrack )
					{
						PFOResidual = this->getPFOResidual( pfoFourMomentum , mcpFourMomentum );
						PFOCovMatPolar = this->getPFOCovMatPolarCoordinate( pfoFourMomentum , outputCovMatrix );
						float ThetaError = sqrt( PFOCovMatPolar[ 0 ] );
						float PhiError = sqrt( PFOCovMatPolar[ 2 ] );
						float EnergyError = sqrt( PFOCovMatPolar[ 9 ] );
						h_ResidualEnergy_CH->Fill( PFOResidual[ 0 ] );
						h_ResidualTheta_CH->Fill( PFOResidual[ 1 ] );
						h_ResidualPhi_CH->Fill( PFOResidual[ 2 ] );
						h_ErrorEnergy_CH->Fill( EnergyError );
						h_ErrorTheta_CH->Fill( ThetaError );
						h_ErrorPhi_CH->Fill( PhiError );
						h_NormalizedResidualEnergy_CH->Fill( PFOResidual[ 0 ] / EnergyError );
						h_NormalizedResidualTheta_CH->Fill( PFOResidual[ 1 ] / ThetaError );
						h_NormalizedResidualPhi_CH->Fill( PFOResidual[ 2 ] / PhiError );
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
				if ( mcpFourMomentum.E() != 0 ) TotalMCPFourMomentum_R += mcpFourMomentum;
				std::vector<float> pfoTrackCovMat( 10, 0.0 );
				streamlog_out(DEBUG) << "PFO is charged (with " << nTrackspfo << " tracks), CovMatrix is set using track information" << std::endl;
				for ( int i_trk = 0 ; i_trk < nTrackspfo ; ++i_trk )
				{
					Track* mytrack		= pfoTracks[ i_trk ];
					double track_mass = this->getTrackMass( pLCEvent , mytrack );
					if ( track_mass == 0.0 )
					{
						track_mass = pion_mass;
						streamlog_out(MESSAGE) << "Couldn't Find track mass, default mass (pion mass) is set" << std::endl;
					}
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
						if ( m_updateCharged ) outputCovMatrix[i] += pfoTrackCovMat[i];
					}
					if ( m_updatePFO4Momentum )
					{
						pfoFourMomentum		+= trackFourMomentum;
					}
					else
					{
						pfoFourMomentum	= TLorentzVector( inputPFO->getMomentum()[0] , inputPFO->getMomentum()[1] , inputPFO->getMomentum()[2] , inputPFO->getEnergy() );
					}
				}
				if ( mcpFourMomentum.E() != 0 )
				{
					TotalPFOFourMomentum += pfoFourMomentum;
					TotalPFOFourMomentum_R += pfoFourMomentum;
					for ( int i = 0 ; i < 10 ; ++i )
					{
						TotaloutputCovMatrix_R[ i ] += outputCovMatrix[ i ];
					}
				}
				if ( mcpFourMomentum.E() != 0 )
				{
					PFOResidual = this->getPFOResidual( pfoFourMomentum , mcpFourMomentum );
					PFOCovMatPolar = this->getPFOCovMatPolarCoordinate( pfoFourMomentum , outputCovMatrix );
					float ThetaError = sqrt( PFOCovMatPolar[ 0 ] );
					float PhiError = sqrt( PFOCovMatPolar[ 2 ] );
					float EnergyError = sqrt( PFOCovMatPolar[ 9 ] );
					h_ResidualEnergy_CH->Fill( PFOResidual[ 0 ] );
					h_ResidualTheta_CH->Fill( PFOResidual[ 1 ] );
					h_ResidualPhi_CH->Fill( PFOResidual[ 2 ] );
					h_ErrorEnergy_CH->Fill( EnergyError );
					h_ErrorTheta_CH->Fill( ThetaError );
					h_ErrorPhi_CH->Fill( PhiError );
					h_NormalizedResidualEnergy_CH->Fill( PFOResidual[ 0 ] / EnergyError );
					h_NormalizedResidualTheta_CH->Fill( PFOResidual[ 1 ] / ThetaError );
					h_NormalizedResidualPhi_CH->Fill( PFOResidual[ 2 ] / PhiError );
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
			if ( mcpFourMomentum.E() != 0 )
			for ( int i = 0 ; i < 10 ; ++i )
			{
				TotaloutputCovMatrix[ i ] += outputCovMatrix[ i ];
			}
			outputPFOMomentum[0]	=	pfoFourMomentum.Px();
			outputPFOMomentum[1]	=	pfoFourMomentum.Py();
			outputPFOMomentum[2]	=	pfoFourMomentum.Pz();
			outputPFOEnergy	=	pfoFourMomentum.E();
			outputPFO->setType(inputPFO->getType());
			outputPFO->setMomentum( outputPFOMomentum );
			outputPFO->setEnergy( outputPFOEnergy );
			outputPFO->setMass( pfoMass );
			outputPFO->setCovMatrix(outputCovMatrix);
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
		TotalPFOResidual_NH = this->getPFOResidual( TotalPFOFourMomentum_NH , TotalMCPFourMomentum_NH );
		TotalPFOCovMatPolar_NH = this->getPFOCovMatPolarCoordinate( TotalPFOFourMomentum_NH , TotaloutputCovMatrix_NH );
		h_ResidualTotalEnergy_NH->Fill( TotalPFOResidual_NH[ 0 ] );
		h_ResidualTotalTheta_NH->Fill( TotalPFOResidual_NH[ 1 ] );
		h_ResidualTotalPhi_NH->Fill( TotalPFOResidual_NH[ 2 ] );
		h_ErrorTotalEnergy_NH->Fill( sqrt( TotaloutputCovMatrix_NH[ 9 ] ) );
		h_ErrorTotalTheta_NH->Fill( sqrt( TotalPFOCovMatPolar_NH[ 0 ] ) );
		h_ErrorTotalPhi_NH->Fill( sqrt( TotalPFOCovMatPolar_NH[ 2 ] ) );
		h_NormalizedResidualTotalEnergy_NH->Fill( TotalPFOResidual_NH[ 0 ] / sqrt( TotaloutputCovMatrix_NH[ 9 ] ) );
		h_NormalizedResidualTotalTheta_NH->Fill( TotalPFOResidual_NH[ 1 ] / sqrt( TotalPFOCovMatPolar_NH[ 0 ] ) );
		h_NormalizedResidualTotalPhi_NH->Fill( TotalPFOResidual_NH[ 2 ] / sqrt( TotalPFOCovMatPolar_NH[ 2 ] ) );

		TotalPFOResidual = this->getPFOResidual( TotalPFOFourMomentum , TotalMCPFourMomentum );
		TotalPFOCovMatPolar = this->getPFOCovMatPolarCoordinate( TotalPFOFourMomentum , TotaloutputCovMatrix );
		h_ResidualTotalEnergy->Fill( TotalPFOResidual[ 0 ] );
		h_ResidualTotalTheta->Fill( TotalPFOResidual[ 1 ] );
		h_ResidualTotalPhi->Fill( TotalPFOResidual[ 2 ] );
		h_ErrorTotalEnergy->Fill( sqrt( TotaloutputCovMatrix[ 9 ] ) );
		h_ErrorTotalTheta->Fill( sqrt( TotalPFOCovMatPolar[ 0 ] ) );
		h_ErrorTotalPhi->Fill( sqrt( TotalPFOCovMatPolar[ 2 ] ) );
		h_NormalizedResidualTotalEnergy->Fill( TotalPFOResidual[ 0 ] / sqrt( TotaloutputCovMatrix[ 9 ] ) );
		h_NormalizedResidualTotalTheta->Fill( TotalPFOResidual[ 1 ] / sqrt( TotalPFOCovMatPolar[ 0 ] ) );
		h_NormalizedResidualTotalPhi->Fill( TotalPFOResidual[ 2 ] / sqrt( TotalPFOCovMatPolar[ 2 ] ) );

		TotalPFOResidual_R = this->getPFOResidual( TotalPFOFourMomentum_R , TotalMCPFourMomentum_R );
		TotalPFOCovMatPolar_R = this->getPFOCovMatPolarCoordinate( TotalPFOFourMomentum_R , TotaloutputCovMatrix_R );
		h_ResidualTotalEnergy_R->Fill( TotalPFOResidual_R[ 0 ] );
		h_ResidualTotalTheta_R->Fill( TotalPFOResidual_R[ 1 ] );
		h_ResidualTotalPhi_R->Fill( TotalPFOResidual_R[ 2 ] );
		h_ErrorTotalEnergy_R->Fill( sqrt( TotaloutputCovMatrix_R[ 9 ] ) );
		h_ErrorTotalTheta_R->Fill( sqrt( TotalPFOCovMatPolar_R[ 0 ] ) );
		h_ErrorTotalPhi_R->Fill( sqrt( TotalPFOCovMatPolar_R[ 2 ] ) );
		h_NormalizedResidualTotalEnergy_R->Fill( TotalPFOResidual_R[ 0 ] / sqrt( TotaloutputCovMatrix_R[ 9 ] ) );
		h_NormalizedResidualTotalTheta_R->Fill( TotalPFOResidual_R[ 1 ] / sqrt( TotalPFOCovMatPolar_R[ 0 ] ) );
		h_NormalizedResidualTotalPhi_R->Fill( TotalPFOResidual_R[ 2 ] / sqrt( TotalPFOCovMatPolar_R[ 2 ] ) );

		pLCEvent->addCollection( m_col_outputPfo , m_outputPfoCollection );
        }
        catch(DataNotAvailableException &e)
        {
          streamlog_out(MESSAGE) << "Input collection not found in event " << m_nEvt << std::endl;
        }

}

TLorentzVector AddFourMomentumCovMatAllPFOs::getLinkedMCP( EVENT::LCEvent *pLCEvent, EVENT::ReconstructedParticle* inputPFO , int nTrackspfo , int nClusterspfo )
{
	LCRelationNavigator navClusterMCTruth(pLCEvent->getCollection(m_ClusterMCTruthLinkCollection));
	LCRelationNavigator navMCTruthCluster(pLCEvent->getCollection(m_MCTruthClusterLinkCollection));
	LCRelationNavigator navTrackMCTruth(pLCEvent->getCollection(m_TrackMCTruthLinkCollection));
	LCRelationNavigator navMCTruthTrack(pLCEvent->getCollection(m_MCTruthTrackLinkCollection));
	streamlog_out(DEBUG) << "PFO has " << nTrackspfo << " tracks and " << nClusterspfo << " clusters" << std::endl;

	int NeutralsPDGCode[14]{11,13,22,130,211,310,321,2112,2212,3112,3122,3222,3312,3322};
//	int ChargedPDGCode[14]{11,13,22,211,321,2212,3112,3122,3222,3312,3322};
	bool PFOlinkedtoMCP = false;
	TLorentzVector mcpFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
	if ( nTrackspfo == 0 )
	{
		const EVENT::LCObjectVec& mcpvec = navClusterMCTruth.getRelatedToObjects(inputPFO->getClusters()[0]);
		const EVENT::FloatVec&  mcpweightvec = navClusterMCTruth.getRelatedToWeights(inputPFO->getClusters()[0]);
		MCParticle *linkedMCP;
		double maxweightPFOtoMCP = 0.0;
		int iPFOtoMCPmax = -1;
		int iMCPtoPFOmax = -1;
		streamlog_out(DEBUG) << "PFO is neutral (without track), pfoType: " << inputPFO->getType() << " , looking for linked " << mcpvec.size() << " MCPs" << std::endl;
		for ( unsigned int i_mcp = 0; i_mcp < mcpvec.size(); i_mcp++ )
		{
			double mcp_weight = mcpweightvec.at(i_mcp);
			MCParticle *testMCP = (MCParticle *) mcpvec.at(i_mcp);
			streamlog_out(DEBUG) << "checking linked MCP at " << i_mcp << " , MCP PDG = " << testMCP->getPDG() << " , link weight = " << mcp_weight << std::endl;
			if ( mcp_weight > maxweightPFOtoMCP && mcp_weight >= 0.9 )
			{
				maxweightPFOtoMCP = mcp_weight;
				iPFOtoMCPmax = i_mcp;
				streamlog_out(DEBUG) << "linkedMCP: " << i_mcp << " has PDG: " << testMCP->getPDG() << " and PFO to MCP Link has weight = " << mcp_weight << std::endl;
			}
		}
		if ( iPFOtoMCPmax != -1 )
		{
			h_NeutPFO_Weight->Fill( maxweightPFOtoMCP );
			linkedMCP = (MCParticle *) mcpvec.at(iPFOtoMCPmax);
			streamlog_out(DEBUG) << "Found linked MCP, MCP PDG: " << linkedMCP->getPDG() << " , link weight = " << maxweightPFOtoMCP << std::endl;
			Cluster *testCluster;
			const EVENT::LCObjectVec& clustervec = navMCTruthCluster.getRelatedToObjects(linkedMCP);
			const EVENT::FloatVec&  clusterweightvec = navMCTruthCluster.getRelatedToWeights(linkedMCP);
			double maxweightMCPtoPFO = 0.;
			for ( unsigned int i_cluster = 0; i_cluster < clustervec.size(); i_cluster++ )
			{
				double cluster_weight = clusterweightvec.at(i_cluster);
				testCluster = (Cluster *) clustervec.at(i_cluster);
				if ( cluster_weight > maxweightMCPtoPFO && cluster_weight >= 0.9 )
				{
					maxweightMCPtoPFO = cluster_weight;
					iMCPtoPFOmax = i_cluster;
				}
			}
			if ( iMCPtoPFOmax != -1 && testCluster == inputPFO->getClusters()[0] )
			{
				PFOlinkedtoMCP = true;
				h_NeutPFO_PDG->Fill( linkedMCP->getPDG() );
				if ( inputPFO->getType() != 22 ) streamlog_out(DEBUG) << "Initial PFO type: " << inputPFO->getType() << "	, linked MCP PDG(weight): " << linkedMCP->getPDG() << " (" << maxweightPFOtoMCP << ")	, linked-back PFO type(weight): " << inputPFO->getType() << " (" << maxweightMCPtoPFO << ")" << std::endl;
				bool KnownPFO = false;
				for ( int l = 0 ; l < 14 ; ++l)
				{
					if ( abs( linkedMCP->getPDG() ) == NeutralsPDGCode[ l ] )
					{
						h_NeutPFO_TYPE->Fill( l );
						KnownPFO = true;
						if ( inputPFO->getType() == 22 )
						{
							h_NeutPFO_IDasPhoton->Fill( l );
						}
						else
						{
							h_NeutPFO_IDasOther->Fill( l );
						}
					}
				}
				if ( !KnownPFO )
				{
					h_NeutPFO_TYPE->Fill( 14 );
					if ( inputPFO->getType() == 22 )
					{
						h_NeutPFO_IDasPhoton->Fill( 14 );
					}
					else
					{
						h_NeutPFO_IDasOther->Fill( 14 );
					}
				}
				mcpFourMomentum = TLorentzVector( linkedMCP->getMomentum()[0] , linkedMCP->getMomentum()[1] , linkedMCP->getMomentum()[2] , linkedMCP->getEnergy() );
				h_NHEnergy->Fill( mcpFourMomentum.E() , inputPFO->getEnergy() );
				h_NH_EclusterPlusMass_Emcp->Fill( mcpFourMomentum.E() , ( inputPFO->getClusters()[0] )->getEnergy() + inputPFO->getMass() );
			}
		}
	}
	else
	{
		for ( int i_trk = 0 ; i_trk < nTrackspfo ; ++i_trk )
		{
//			Track* mytrack		= ( inputPFO->getTracks() )[ i_trk ];
			const EVENT::LCObjectVec& mcpvec = navTrackMCTruth.getRelatedToObjects(inputPFO->getTracks()[ i_trk ]);
			const EVENT::FloatVec&  mcpweightvec = navTrackMCTruth.getRelatedToWeights(inputPFO->getTracks()[ i_trk ]);
			MCParticle *linkedMCP;
			double maxweightPFOtoMCP = 0.0;
			int iPFOtoMCPmax = -1;
			int iMCPtoPFOmax = -1;
//			streamlog_out(DEBUG) << "PFO is neutral (without track), pfoType: " << inputPFO->getType() << " , looking for linked " << mcpvec.size() << " MCPs" << std::endl;
			int n_mcp = mcpvec.size();
			for ( int i_mcp = 0; i_mcp < n_mcp; i_mcp++ )
			{
				double mcp_weight = mcpweightvec.at(i_mcp);
				MCParticle *testMCP = (MCParticle *) mcpvec.at(i_mcp);
				streamlog_out(DEBUG) << "checking linked MCP at " << i_mcp << " , MCP PDG = " << testMCP->getPDG() << " , link weight = " << mcp_weight << std::endl;
				if ( mcp_weight > maxweightPFOtoMCP && mcp_weight >= 0.9 )
				{
					maxweightPFOtoMCP = mcp_weight;
					iPFOtoMCPmax = i_mcp;
					streamlog_out(DEBUG) << "linkedMCP: " << i_mcp << " has PDG: " << testMCP->getPDG() << " and PFO to MCP Link has weight = " << mcp_weight << std::endl;
				}
			}

			if ( iPFOtoMCPmax != -1 )
			{
//				h_NeutPFO_Weight->Fill( maxweightPFOtoMCP );
				linkedMCP = (MCParticle *) mcpvec.at(iPFOtoMCPmax);
				streamlog_out(DEBUG) << "Found linked MCP, MCP PDG: " << linkedMCP->getPDG() << " , link weight = " << maxweightPFOtoMCP << std::endl;
				Track *testTrack;
				const EVENT::LCObjectVec& trackvec = navMCTruthTrack.getRelatedToObjects(linkedMCP);
				const EVENT::FloatVec&  trackweightvec = navMCTruthTrack.getRelatedToWeights(linkedMCP);
				double maxweightMCPtoPFO = 0.;
				for ( unsigned int i_track = 0; i_track < trackvec.size(); i_track++ )
				{
					double Track_weight = trackweightvec.at(i_track);
					testTrack = (Track *) trackvec.at(i_track);
					if ( Track_weight > maxweightMCPtoPFO && Track_weight >= 0.9 )
					{
						maxweightMCPtoPFO = Track_weight;
						iMCPtoPFOmax = i_track;
					}
				}
				if ( iMCPtoPFOmax != -1 && testTrack == inputPFO->getTracks()[ i_trk ] )
				{
					PFOlinkedtoMCP = true;
					h_NeutPFO_PDG->Fill( linkedMCP->getPDG() );
					if ( inputPFO->getType() != 22 ) streamlog_out(DEBUG) << "Initial PFO type: " << inputPFO->getType() << "	, linked MCP PDG(weight): " << linkedMCP->getPDG() << " (" << maxweightPFOtoMCP << ")	, linked-back PFO type(weight): " << inputPFO->getType() << " (" << maxweightMCPtoPFO << ")" << std::endl;
					mcpFourMomentum += TLorentzVector( linkedMCP->getMomentum()[0] , linkedMCP->getMomentum()[1] , linkedMCP->getMomentum()[2] , linkedMCP->getEnergy() );
/*
					bool KnownPFO = false;
					for ( int l = 0 ; l < 14 ; ++l)
					{
						if ( abs( linkedMCP->getPDG() ) == NeutralsPDGCode[ l ] )
						{
							h_NeutPFO_TYPE->Fill( l );
							KnownPFO = true;
							if ( inputPFO->getType() == 22 )
							{
								h_NeutPFO_IDasPhoton->Fill( l );
							}
							else
							{
								h_NeutPFO_IDasOther->Fill( l );
							}
						}
					}
					if ( !KnownPFO )
					{
						h_NeutPFO_TYPE->Fill( 14 );
						if ( inputPFO->getType() == 22 )
						{
							h_NeutPFO_IDasPhoton->Fill( 14 );
						}
						else
						{
							h_NeutPFO_IDasOther->Fill( 14 );
						}
					}
*/				}
			}
		}
	}

	if ( PFOlinkedtoMCP )
	{
		return mcpFourMomentum;
	}
	else
	{
		return TLorentzVector( 0.0 , 0.0 , 0.0 , 0.0 );
	}

}

std::vector<float> AddFourMomentumCovMatAllPFOs::UpdateNeutralPFOCovMat( TVector3 clusterPosition , float pfoEc , float pfoMass , std::vector<float> clusterPositionError , float clusterEnergyError )
{

//	Obtain covariance matrix on (px,py,pz,E) from the
//	covariance matrix on cluster parameters (px,py,pz,Ek=Ec=E-E0).
//	=> E = Ek + m	;	|p| = sqrt( Ek^2 + 2mEk )
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
//			Dpx/DEk		Dpy/DEk		Dpz/DEk		DE/DEk
//
//
//
//
//
//			|P|.(r2-x2)/r3		-|P|.x.y/r3		-|P|.x.z/r3		0
//
//			-|P|.y.x/r3		|P|.(r2-y2)/r3		-|P|.y.z/r3		0
//	J =
//			-|P|.z.x/r3		-|P|.z.y/r3		|P|.(r2-z2)/r3		0
//
//			(E/|p|).(x/r)		(E/|p|).(y/r)		(E/|p|).(z/r)		1
//
//
//
//
//	CovMatrix elements in terms of cluster position error and cluster energy error:
//
//			x.x			x.y			x.z			x.Ec
//
//			y.x			y.y			y.z			y.Ec
//	Cov =
//			z.x			z.y			z.z			z.Ec
//
//			Ec.x			Ec.y			Ec.z			Ec.Ec
//
//
//

	const int rows			= 4; // n rows jacobian
	const int columns		= 4; // n columns jacobian
	const int kspace_time_dim	= 4;

	TMatrixD covMatrixMomenta(kspace_time_dim,kspace_time_dim);
	std::vector<float> covP;

//	pfoMass			= 0.0;

	float pfoX		=	clusterPosition.X();
	float pfoY		=	clusterPosition.Y();
	float pfoZ		=	clusterPosition.Z();
	float pfoR		=	std::sqrt( pow( pfoX , 2 ) + pow( pfoY , 2 ) + pow( pfoZ , 2 ) );
	float pfoX2		=	pow( pfoX , 2 );
	float pfoY2		=	pow( pfoY , 2 );
	float pfoZ2		=	pow( pfoZ , 2 );
	float pfoR2		=	pow( pfoR , 2 );
	float pfoR3		=	pow( pfoR , 3 );
	float SigmaX2		=	clusterPositionError[ 0 ];
	float SigmaXY		=	clusterPositionError[ 1 ];
	float SigmaY2		=	clusterPositionError[ 2 ];
	float SigmaXZ		=	clusterPositionError[ 3 ];
	float SigmaYZ		=	clusterPositionError[ 4 ];
	float SigmaZ2		=	clusterPositionError[ 5 ];
	float SigmaE2		=	pow( clusterEnergyError , 2 );

	float pfoE;
	float pfoP;
	float derivative_coeff	= 1.0;
	if ( m_isClusterEnergyKinEnergy )
	{
		pfoE			= pfoEc + pfoMass;
		pfoP			= sqrt( pow( pfoEc , 2 ) + 2 * pfoMass * pfoEc );
//		derivative_coeff	= 1.;
	}
	else
	{
		pfoP			= pfoEc;
		pfoE			= sqrt( pow( pfoP , 2 ) + pow( pfoMass , 2 ) );
//		derivative_coeff	= pfoP / pfoE;
	}

	streamlog_out(DEBUG) << "Cluster information obtained" << std::endl;

//	Define array with jacobian matrix elements by rows
	double jacobian_by_rows[rows*columns] =
	{
		pfoP * ( pfoR2 - pfoX2 ) / pfoR3			,	-pfoP * pfoX * pfoY / pfoR3				,	-pfoP * pfoX * pfoZ / pfoR3				,	0			,
		-pfoP * pfoY * pfoX / pfoR3				,	pfoP * ( pfoR2 - pfoY2 ) / pfoR3			,	-pfoP * pfoY * pfoZ / pfoR3				,	0			,
		-pfoP * pfoZ * pfoX / pfoR3				,	-pfoP * pfoZ * pfoY / pfoR3				,	pfoP * ( pfoR2 - pfoZ2 ) / pfoR3			,	0			,
		derivative_coeff * pfoE * pfoX / ( pfoP * pfoR )	,	derivative_coeff * pfoE * pfoY / ( pfoP * pfoR )	,	derivative_coeff * pfoE * pfoZ / ( pfoP * pfoR )	,	derivative_coeff * 1
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
				0			,	0			,	0			,	SigmaE2
			};
	streamlog_out(DEBUG) << "cluster covariance matrix array formed by rows" << std::endl;

	TMatrixD covMatrix_cluster(rows,rows, cluster_cov_matrix_by_rows, "C");
	streamlog_out(DEBUG) << "cluster covariance matrix array converted to cluster covariance matrix" << std::endl;

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


std::vector<float> AddFourMomentumCovMatAllPFOs::getClusterDirectionError( TVector3 clusterPosition , std::vector<float> clusterPositionError )
{

//	Obtain covariance matrix on (R,Theta,Phi) from the
//	covariance matrix on cluster parameters (x,y,z).
//	=> R^2 = x^2 + y^2 + z^2	;	tan(Theta) = sqrt( x^2 + y^2 ) / z	;	tan(Phi) = y / x
//	define the jacobian as the 4x4 matrix:
//
//
//
//			DR/Dx			DTheta/Dx		DPhi/Dx
//
//	J =		DR/Dy			DTheta/Dy		DPhi/Dy
//
//			DR/Dz			DTheta/Dz		DPhi/Dz
//
//
//
//
//			x/R			x.z/(R3.sqrt(1-(z2/R2)))		-y/(x2+y2)
//
//	J =		y/R			y.z/(R3.sqrt(1-(z2/R2)))		x/(x2+y2)
//
//			z/R			sqrt(1-(z2/R2))/R			0
//
//
//
//
//
//	CovMatrix elements in terms of cluster position error:
//
//			x.x			x.y			x.z
//
//	Cov =		y.x			y.y			y.z
//
//			z.x			z.y			z.z
//
//
//
//

	const int rows			= 3; // n rows jacobian
	const int columns		= 3; // n columns jacobian
	const int kspace_time_dim	= 3;

	TMatrixD covMatrixMomenta(kspace_time_dim,kspace_time_dim);
	std::vector<float> covR;

//	pfoMass			= 0.0;

	float pfoX		=	clusterPosition.X();
	float pfoY		=	clusterPosition.Y();
	float pfoZ		=	clusterPosition.Z();
	float pfoR		=	std::sqrt( pow( pfoX , 2 ) + pow( pfoY , 2 ) + pow( pfoZ , 2 ) );
	float pfoR2		=	pow( pfoR , 2 );
	float pfoR3		=	pow( pfoR , 3 );
	float pfoX2		=	pow( pfoX , 2 );
	float pfoY2		=	pow( pfoY , 2 );
	float pfoZ2		=	pow( pfoZ , 2 );
	float SigmaX2		=	clusterPositionError[ 0 ];
	float SigmaXY		=	clusterPositionError[ 1 ];
	float SigmaY2		=	clusterPositionError[ 2 ];
	float SigmaXZ		=	clusterPositionError[ 3 ];
	float SigmaYZ		=	clusterPositionError[ 4 ];
	float SigmaZ2		=	clusterPositionError[ 5 ];

	streamlog_out(DEBUG) << "Cluster information obtained" << std::endl;

//	Define array with jacobian matrix elements by rows
	double jacobian_by_rows[rows*columns] =
			{
				pfoX / pfoR		,		pfoX * pfoZ / ( pfoR3 * sqrt( 1 - pfoZ2 / pfoR2 ) )		,	-pfoY / ( pfoY2 + pfoX2 )		,
				pfoY / pfoR		,		pfoY * pfoZ / ( pfoR3 * sqrt( 1 - pfoZ2 / pfoR2 ) )		,	pfoX / ( pfoY2 + pfoX2 )		,
				pfoZ / pfoR		,		sqrt( 1 - pfoZ2 / pfoR2 ) / pfoR				,		0
			};

	streamlog_out(DEBUG) << "Jacobian array formed by rows" << std::endl;

//	construct the Jacobian using previous array ("F" if filling by columns, "C", if filling by rows, $ROOTSYS/math/matrix/src/TMatrixT.cxx)
	TMatrixD jacobian(rows,columns, jacobian_by_rows, "C");
	streamlog_out(DEBUG) << "Jacobian array converted to Jacobian matrix" << std::endl;

//	cluster covariance matrix by rows
	double cluster_cov_matrix_by_rows[rows*rows] =
			{
				SigmaX2		,	SigmaXY		,	SigmaXZ		,
				SigmaXY		,	SigmaY2		,	SigmaYZ		,
				SigmaXZ		,	SigmaYZ		,	SigmaZ2
			};
	streamlog_out(DEBUG) << "cluster covariance matrix array formed by rows" << std::endl;

	TMatrixD covMatrix_cluster(rows,rows, cluster_cov_matrix_by_rows, "C");
	streamlog_out(DEBUG) << "cluster covariance matrix array converted to cluster covariance matrix" << std::endl;

	covMatrixMomenta.Mult( TMatrixD( jacobian ,
					TMatrixD::kTransposeMult ,
					covMatrix_cluster) ,
					jacobian
					);
	streamlog_out(DEBUG) << "cluster covariance matrix array in cartesian coordinate system converted to cluster covariance matrix array in spherical (polar) coordinate system" << std::endl;

	covR.push_back( covMatrixMomenta(0,0) ); // R-R
	covR.push_back( covMatrixMomenta(1,0) ); // Theta-R
	covR.push_back( covMatrixMomenta(1,1) ); // Theta-Theta
	covR.push_back( covMatrixMomenta(2,0) ); // Phi-R
	covR.push_back( covMatrixMomenta(2,1) ); // Phi-Theta
	covR.push_back( covMatrixMomenta(2,2) ); // Phi-Phi
	streamlog_out(DEBUG) << "FourMomentumCovarianceMatrix Filled succesfully" << std::endl;

	return covR;

}

std::vector<float> AddFourMomentumCovMatAllPFOs::getPFOCovMatPolarCoordinate( TLorentzVector pfoFourMomentum , std::vector<float> pfoCovMat )
{

//	Obtain covariance matrix on (Theta,Phi,P,E) from the
//	covariance matrix on (Px,Py,Pz,E).
//	=> P^2 = Px^2 + Py^2 + Pz^2	;	tan(Theta) = sqrt( Px^2 + Py^2 ) / Pz	;	tan(Phi) = Py / Px
//	define the jacobian as the 4x4 matrix:
//
//
//
//			DTheta/DPx			DPhi/DPx		DP/DPx			DE/DPx
//
//			DTheta/DPy			DPhi/DPy		DP/DPy			DE/DPy
//	J =
//			DTheta/DPz			DPhi/DPz		DP/DPz			DE/DPz
//
//			DTheta/DE			DPhi/DE		DP/DE			DE/DE
//
//
//
//
//			(Px.Pz)/(Pt.P2)		-Py/Pt2		Px/P			0
//
//			(Py.Pz)/(Pt.P2)		Px/Pt2			Py/P			0
//	J =
//			-Pt/P2				0			Pz/P			0
//
//			0				0			0			1
//
//
//
//
//
//	CovMatrix elements in terms of cluster position error:
//
//			Px.Px			Px.Py			Px.Pz			Px.E
//
//			Py.Px			Py.Py			Py.Pz			Py.E
//	Cov =	
//			Pz.Px			Pz.Py			Pz.Pz			Pz.E
//
//			E.Px			E.Py			E.Pz			E.E
//
//
//
//

	const int rows			= 4; // n rows jacobian
	const int columns		= 4; // n columns jacobian
	const int kspace_time_dim	= 4;

	TMatrixD covMatrixMomenta(kspace_time_dim,kspace_time_dim);
	std::vector<float> covP;

//	pfoMass			= 0.0;

	float totPx		=	pfoFourMomentum.Px();
	float totPy		=	pfoFourMomentum.Py();
	float totPz		=	pfoFourMomentum.Pz();
	float totP		=	std::sqrt( pow( totPx , 2 ) + pow( totPy , 2 ) + pow( totPz , 2 ) );
	float totP2		=	pow( totP , 2 );
	float totPt		=	std::sqrt( pow( totPx , 2 ) + pow( totPy , 2 ) );
	float totPt2		=	pow( totPt , 2 );
	float SigmaPx2		=	pfoCovMat[ 0 ];
	float SigmaPxPy	=	pfoCovMat[ 1 ];
	float SigmaPy2		=	pfoCovMat[ 2 ];
	float SigmaPxPz	=	pfoCovMat[ 3 ];
	float SigmaPyPz	=	pfoCovMat[ 4 ];
	float SigmaPz2		=	pfoCovMat[ 5 ];
	float SigmaPxE		=	pfoCovMat[ 6 ];
	float SigmaPyE		=	pfoCovMat[ 7 ];
	float SigmaPzE		=	pfoCovMat[ 8 ];
	float SigmaE2		=	pfoCovMat[ 9 ];

	streamlog_out(DEBUG) << "Cluster information obtained" << std::endl;

//	Define array with jacobian matrix elements by rows
	double jacobian_by_rows[rows*columns] =
	{
		( totPx * totPz ) / ( totPt * totP2 )		,	-totPy / totPt2	,	totPx / totP		,	0	,
		( totPy * totPz ) / ( totPt * totP2 )		,	-totPx / totPt2	,	totPy / totP		,	0	,
		-totPt / totP2					,	0			,	totPz / totP		,	0	,
		0						,	0			,	0			,	1
	};

	streamlog_out(DEBUG) << "Jacobian array formed by rows" << std::endl;

//	construct the Jacobian using previous array ("F" if filling by columns, "C", if filling by rows, $ROOTSYS/math/matrix/src/TMatrixT.cxx)
	TMatrixD jacobian(rows,columns, jacobian_by_rows, "C");
	streamlog_out(DEBUG) << "Jacobian array converted to Jacobian matrix" << std::endl;

//	cluster covariance matrix by rows
	double totPFO_cov_matrix_by_rows[rows*rows] =
			{
				SigmaPx2		,	SigmaPxPy		,	SigmaPxPz		,	SigmaPxE	,
				SigmaPxPy		,	SigmaPy2		,	SigmaPyPz		,	SigmaPyE	,
				SigmaPxPz		,	SigmaPyPz		,	SigmaPz2		,	SigmaPzE	,
				SigmaPxE		,	SigmaPyE		,	SigmaPzE		,	SigmaE2
			};
	streamlog_out(DEBUG) << "cluster covariance matrix array formed by rows" << std::endl;

	TMatrixD covMatrix(rows,rows, totPFO_cov_matrix_by_rows, "C");
	streamlog_out(DEBUG) << "cluster covariance matrix array converted to cluster covariance matrix" << std::endl;

	covMatrixMomenta.Mult( TMatrixD( jacobian ,
					TMatrixD::kTransposeMult ,
					covMatrix) ,
					jacobian
					);
	streamlog_out(DEBUG) << "cluster covariance matrix array in cartesian coordinate system converted to cluster covariance matrix array in spherical (polar) coordinate system" << std::endl;

	covP.push_back( covMatrixMomenta(0,0) ); // Theta-Theta
	covP.push_back( covMatrixMomenta(1,0) ); // Theta-Phi
	covP.push_back( covMatrixMomenta(1,1) ); // Phi-Phi
	covP.push_back( covMatrixMomenta(2,0) ); // Theta-P
	covP.push_back( covMatrixMomenta(2,1) ); // Phi-P
	covP.push_back( covMatrixMomenta(2,2) ); // P-P
	covP.push_back( covMatrixMomenta(3,0) ); // Theta-E
	covP.push_back( covMatrixMomenta(3,1) ); // Phi-E
	covP.push_back( covMatrixMomenta(3,2) ); // P-E
	covP.push_back( covMatrixMomenta(3,3) ); // E-E
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
	double trackPhi		=	MyTrack->getPhi();
	double trackOmega		=	MyTrack->getOmega();
	std::vector<float> trackCovMat =	MyTrack->getCovMatrix();
	double trackPt			=	eB / fabs( trackOmega );
	double trackPx			= 	trackPt * TMath::Cos( trackPhi );
	double trackPy			= 	trackPt * TMath::Sin( trackPhi );
	double trackPz			= 	trackPt * trackTanLambda;
	double trackP			= 	std::sqrt( pow( trackPt , 2 ) + pow( trackPz , 2 ) );
	double trackE			= 	std::sqrt( pow( trackP , 2 ) + pow( trackMass , 2 ) );

	float SigmaPhi2		=	trackCovMat[  2 ];
	float SigmaPhiSigmaOmega	=	trackCovMat[  4 ];
	float SigmaOmega2		=	trackCovMat[  5 ];
	float SigmaTanLambdaSigmaPhi	=	trackCovMat[ 11 ];
	float SigmaTanLambdaSigmaOmega =	trackCovMat[ 12 ];
	float SigmaTanLambda2		=	trackCovMat[ 14 ];

	streamlog_out(DEBUG) << "track information obtained" << std::endl;

//	Define array with jacobian matrix elements by rows
	double jacobian_by_rows[rows*columns] =
	{
		0			,	0			,	trackPt			,	trackPz * trackPt / trackE			,
		-trackPy		,	trackPx		,	0				,	0						,
		-trackPx / trackOmega	,	-trackPy / trackOmega	,	-trackPz / trackOmega		,	-trackP * trackP / ( trackE * trackOmega )
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

std::vector<float> AddFourMomentumCovMatAllPFOs::getPFOResidual( TLorentzVector pfoFourMomentum , TLorentzVector mcpFourMomentum )
{
	std::vector<float> pfoResidual;

	float pfoPx		= pfoFourMomentum.Px();
	float pfoPy		= pfoFourMomentum.Py();
	float pfoPz		= pfoFourMomentum.Pz();
	float pfoE		= pfoFourMomentum.E();
	TVector3 pfoPvec( pfoPx , pfoPy , pfoPz ); pfoPvec.SetMag(1.0);
	TVector3 pfoPTvec( pfoPx , pfoPy , 0.0 ); pfoPTvec.SetMag(1.0);
	float pfoTheta		= pfoPvec.Theta();
	float pfoPhi		= pfoPvec.Phi();
	streamlog_out(DEBUG) << "PFO 4-Momentum ( px , py , pz , E , Theta , Phi ) = ( " << pfoPx << " , " << pfoPy << " , " << pfoPz << " , " << pfoE << " , " << pfoTheta << " , " << pfoPhi << " )" << std::endl;

	float mcpPx		= mcpFourMomentum.Px();
	float mcpPy		= mcpFourMomentum.Py();
	float mcpPz		= mcpFourMomentum.Pz();
	float mcpE		= mcpFourMomentum.E();
	TVector3 mcpPvec( mcpPx , mcpPy , mcpPz ); mcpPvec.SetMag(1.0);
	TVector3 mcpPTvec( mcpPx , mcpPy , 0.0 ); mcpPTvec.SetMag(1.0);
	float mcpTheta		= mcpPvec.Theta();
	float mcpPhi		= mcpPvec.Phi();
	streamlog_out(DEBUG) << "MCP 4-Momentum ( px , py , pz , E , Theta , Phi ) = ( " << mcpPx << " , " << mcpPy << " , " << mcpPz << " , " << mcpE << " , " << mcpTheta << " , " << mcpPhi << " )" << std::endl;

	float ResidualEnergy	= pfoE - mcpE;
	float ResidualTheta	= pfoTheta - mcpTheta;
	float ResidualPhi	= 0.0;
	if ( pfoPhi > mcpPhi )
	{
		ResidualPhi	= acos( pfoPTvec.Dot( mcpPTvec ) );
	}
	else
	{
		ResidualPhi	= -acos( pfoPTvec.Dot( mcpPTvec ) );
	}
	streamlog_out(DEBUG) << "	Residuals	( deltaE , deltaTheta , deltaPhi ) = ( " << ResidualEnergy << "	,	" << ResidualTheta << "	, " << ResidualPhi << "	)" << std::endl;

	pfoResidual.push_back( ResidualEnergy );
	pfoResidual.push_back( ResidualTheta );
	pfoResidual.push_back( ResidualPhi );

	return pfoResidual;

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
	m_Histograms->cd();
	h_nTracks_PFOCharge->Write();
	h_nClusters_nTracks->Write();
	h_clusterE_pfoE->Write();
	h_NeutPFO_PDG->Write();
	h_NeutPFO_TYPE->Write();
	h_NeutPFO_IDasPhoton->Write();
	h_NeutPFO_IDasOther->Write();
	h_NeutPFO_Mass->Write();
	h_EP_photons->Write();
	h_EP_NeutralHadrons->Write();
	h_NeutPFO_Weight->Write();
	h_NH_EclusterPlusMass_Emcp->Write();
	h_NHEnergy->Write();
	m_CovMatElements->cd();
	m_NeutralPFOswithoutTrak->cd();
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
	m_CovMatElements->cd();
	m_NeutralPFOswith2Trak->cd();
	h_SigmaPx2T->Write();
	h_SigmaPxPyT->Write();
	h_SigmaPy2T->Write();
	h_SigmaPxPzT->Write();
	h_SigmaPyPzT->Write();
	h_SigmaPz2T->Write();
	h_SigmaPxET->Write();
	h_SigmaPyET->Write();
	h_SigmaPzET->Write();
	m_CovMatElements->cd();
	m_ChargedPFOs->cd();
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
	m_CovMatElements->cd();
	m_Histograms->cd();
	m_ErrorParameterization->cd();
	m_Photon->cd();
	h_ResidualEnergy_ph->Write();
	h_ResidualTheta_ph->Write();
	h_ResidualPhi_ph->Write();
	h_ErrorEnergy_ph->Write();
	h_ErrorTheta_ph->Write();
	h_ErrorPhi_ph->Write();
	h_NormalizedResidualEnergy_ph->Write();
	h_NormalizedResidualTheta_ph->Write();
	h_NormalizedResidualPhi_ph->Write();
	m_ErrorParameterization->cd();
	m_NeutralPFO->cd();
	h_ResidualEnergy_NH->Write();
	h_ResidualTheta_NH->Write();
	h_ResidualPhi_NH->Write();
	h_ErrorEnergy_NH->Write();
	h_ErrorTheta_NH->Write();
	h_ErrorPhi_NH->Write();
	h_NormalizedResidualEnergy_NH->Write();
	h_NormalizedResidualTheta_NH->Write();
	h_NormalizedResidualPhi_NH->Write();
	m_ErrorParameterization->cd();
	m_ChargedPFO->cd();
	h_ResidualEnergy_CH->Write();
	h_ResidualTheta_CH->Write();
	h_ResidualPhi_CH->Write();
	h_ErrorEnergy_CH->Write();
	h_ErrorTheta_CH->Write();
	h_ErrorPhi_CH->Write();
	h_NormalizedResidualEnergy_CH->Write();
	h_NormalizedResidualTheta_CH->Write();
	h_NormalizedResidualPhi_CH->Write();
	m_TotalPFOs->cd();
	h_ResidualTotalEnergy->Write();
	h_ResidualTotalTheta->Write();
	h_ResidualTotalPhi->Write();
	h_ErrorTotalEnergy->Write();
	h_ErrorTotalTheta->Write();
	h_ErrorTotalPhi->Write();
	h_NormalizedResidualTotalEnergy->Write();
	h_NormalizedResidualTotalTheta->Write();
	h_NormalizedResidualTotalPhi->Write();
	m_ErrorParameterization->cd();
	m_TotalNHs->cd();
	h_ResidualTotalEnergy_NH->Write();
	h_ResidualTotalTheta_NH->Write();
	h_ResidualTotalPhi_NH->Write();
	h_ErrorTotalEnergy_NH->Write();
	h_ErrorTotalTheta_NH->Write();
	h_ErrorTotalPhi_NH->Write();
	h_NormalizedResidualTotalEnergy_NH->Write();
	h_NormalizedResidualTotalTheta_NH->Write();
	h_NormalizedResidualTotalPhi_NH->Write();
	m_ErrorParameterization->cd();
	m_TotalRest->cd();
	h_ResidualTotalEnergy_R->Write();
	h_ResidualTotalTheta_R->Write();
	h_ResidualTotalPhi_R->Write();
	h_ErrorTotalEnergy_R->Write();
	h_ErrorTotalTheta_R->Write();
	h_ErrorTotalPhi_R->Write();
	h_NormalizedResidualTotalEnergy_R->Write();
	h_NormalizedResidualTotalTheta_R->Write();
	h_NormalizedResidualTotalPhi_R->Write();
	m_ErrorParameterization->cd();
	m_pTFile->Close();
	delete m_pTFile;

//	std::cout << " END : processed events: " << m_nEvtSum << std::endl;

}
