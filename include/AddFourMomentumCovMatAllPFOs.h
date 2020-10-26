#ifndef AddFourMomentumCovMatAllPFOs_h
#define AddFourMomentumCovMatAllPFOs_h 1
#include <marlin/Processor.h>
#include <marlin/Global.h>
#include "UTIL/LCRelationNavigator.h"
#include "UTIL/LCIterator.h"
#include "UTIL/Operators.h"
#include <UTIL/BitField64.h>
#include "UTIL/LCTrackerConf.h"
#include <UTIL/ILDConf.h>
#include <UTIL/BitSet32.h>
#include "EVENT/LCStrVec.h"
#include <EVENT/MCParticle.h>
#include <EVENT/Track.h>
#include "IMPL/LCCollectionVec.h"
#include <IMPL/ReconstructedParticleImpl.h>
#include "TLorentzVector.h"
#include "TVector3.h"
#include "lcio.h"
#include <string>
#include <vector>
#include <math.h>
#include <set>
#include <vector>
#include "TMatrixD.h"
#include <TFile.h>
#include <TTree.h>

class TFile;
class TDirectory;
class TH1F;
class TH1I;
class TH2I;
class TH2F;
class TTree;

using namespace lcio ;
using namespace marlin ;

class AddFourMomentumCovMatAllPFOs : public Processor
{
	public:

		virtual Processor*  newProcessor()
		{
			return new AddFourMomentumCovMatAllPFOs;
		}
		AddFourMomentumCovMatAllPFOs();
		virtual ~AddFourMomentumCovMatAllPFOs() = default;
		AddFourMomentumCovMatAllPFOs(const AddFourMomentumCovMatAllPFOs&) = delete;
		AddFourMomentumCovMatAllPFOs& operator=(const AddFourMomentumCovMatAllPFOs&) = delete;
		virtual void init();
		virtual void Clear();
		virtual void processRunHeader();
		virtual void processEvent( EVENT::LCEvent *pLCEvent );
		std::vector<float> UpdateNeutralPFOCovMat( TVector3 clusterPosition , float pfoEc , float pfoMass , std::vector<float> clusterPositionError , float clusterEnergyError );
		std::vector<float> getClusterDirectionError( TVector3 clusterPosition , std::vector<float> clusterPositionError );
		std::vector<float> UpdateChargedPFOCovMat( EVENT::Track* MyTrack , float trackMass );
		std::vector<float> getPFOResidual( TLorentzVector pfoFourMomentum , TLorentzVector mcpFourMomentum );
		std::vector<float> getPFOCovMatPolarCoordinate( TLorentzVector pfoFourMomentum , std::vector<float> pfoCovMat );
		double getTrackMass( EVENT::LCEvent *pLCEvent, EVENT::Track* inputTrk );
		TLorentzVector getLinkedMCP( EVENT::LCEvent *pLCEvent, EVENT::ReconstructedParticle* inputPFO , int nTrackspfo , int nClusterspfo );
		virtual void check( EVENT::LCEvent *pLCEvent );
		virtual void end();

	private:
		std::string				m_inputPfoCollection{};
		std::string				m_ClusterMCTruthLinkCollection{};
		std::string				m_MCTruthClusterLinkCollection{};
		std::string				m_TrackMCTruthLinkCollection{};
		std::string				m_MCTruthTrackLinkCollection{};
		std::string				m_outputPfoCollection{};
		std::string				m_rootFile{};
		std::string				m_MarlinTrkTracks{};
		std::string				m_MarlinTrkTracksKAON{};
		std::string				m_MarlinTrkTracksPROTON{};

		bool					m_updateNormalNeutrals = false;
		bool					m_updateNeutrals_wTrack = false;
		bool					m_updateCharged = false;
		bool					m_AssumeNeutralPFOMassive = false;
		bool					m_isClusterEnergyKinEnergy = false;
		bool					m_useClusterPositionError = true;
		bool					m_updatePFO4Momentum = false;
		int					m_nRun;
		int					m_nEvt;
		int					m_nRunSum;
		int					m_nEvtSum;
		float					m_Bfield;
		double					c;
		double					mm2m;
		double					eV2GeV;
		double					eB;
		double					proton_mass;
		double					proton_mass_sq;
		double					kaon_mass;
		double					kaon_mass_sq;
		double					pion_mass;
		double					pion_mass_sq;
		TFile					*m_pTFile;
	        TTree					*m_pTTree;
	        TDirectory				*m_Histograms;
	        TDirectory				*m_CovMatElements;
	        TDirectory				*m_NeutralPFOswithoutTrak;
	        TDirectory				*m_NeutralPFOswith2Trak;
	        TDirectory				*m_ChargedPFOs;
	        TDirectory				*m_ErrorParameterization;
	        TDirectory				*m_Photon;
	        TDirectory				*m_NeutralPFO;
	        TDirectory				*m_ChargedPFO;
	        TDirectory				*m_TotalPFOs;
	        TDirectory				*m_TotalNHs;
	        TDirectory				*m_TotalRest;
		TH2I					*h_nTracks_PFOCharge{};
		TH2I					*h_nClusters_nTracks{};
		TH2F					*h_clusterE_pfoE{};
		TH2F					*h_SigmaPx2nT{};
		TH2F					*h_SigmaPxPynT{};
		TH2F					*h_SigmaPy2nT{};
		TH2F					*h_SigmaPxPznT{};
		TH2F					*h_SigmaPyPznT{};
		TH2F					*h_SigmaPz2nT{};
		TH2F					*h_SigmaPxEnT{};
		TH2F					*h_SigmaPyEnT{};
		TH2F					*h_SigmaPzEnT{};
		TH2F					*h_SigmaE2nT{};
		TH2F					*h_SigmaPx2T{};
		TH2F					*h_SigmaPxPyT{};
		TH2F					*h_SigmaPy2T{};
		TH2F					*h_SigmaPxPzT{};
		TH2F					*h_SigmaPyPzT{};
		TH2F					*h_SigmaPz2T{};
		TH2F					*h_SigmaPxET{};
		TH2F					*h_SigmaPyET{};
		TH2F					*h_SigmaPzET{};
		TH2F					*h_SigmaE2T{};
		TH2F					*h_SigmaPx2{};
		TH2F					*h_SigmaPxPy{};
		TH2F					*h_SigmaPy2{};
		TH2F					*h_SigmaPxPz{};
		TH2F					*h_SigmaPyPz{};
		TH2F					*h_SigmaPz2{};
		TH2F					*h_SigmaPxE{};
		TH2F					*h_SigmaPyE{};
		TH2F					*h_SigmaPzE{};
		TH2F					*h_SigmaE2{};
		TH1I					*h_NeutPFO_PDG{};
		TH1I					*h_NeutPFO_TYPE{};
		TH1I					*h_NeutPFO_IDasPhoton{};
		TH1I					*h_NeutPFO_IDasOther{};
		TH1F					*h_NeutPFO_Mass{};
		TH2F					*h_EP_photons{};
		TH2F					*h_EP_NeutralHadrons{};
		TH1F					*h_NeutPFO_Weight{};
		TH1F					*h_ResidualEnergy_ph{};
		TH1F					*h_ResidualTheta_ph{};
		TH1F					*h_ResidualPhi_ph{};
		TH1F					*h_ErrorEnergy_ph{};
		TH1F					*h_ErrorTheta_ph{};
		TH1F					*h_ErrorPhi_ph{};
		TH1F					*h_NormalizedResidualEnergy_ph{};
		TH1F					*h_NormalizedResidualTheta_ph{};
		TH1F					*h_NormalizedResidualPhi_ph{};
		TH1F					*h_ResidualEnergy_NH{};
		TH1F					*h_ResidualTheta_NH{};
		TH1F					*h_ResidualPhi_NH{};
		TH1F					*h_ErrorEnergy_NH{};
		TH1F					*h_ErrorTheta_NH{};
		TH1F					*h_ErrorPhi_NH{};
		TH1F					*h_NormalizedResidualEnergy_NH{};
		TH1F					*h_NormalizedResidualTheta_NH{};
		TH1F					*h_NormalizedResidualPhi_NH{};
		TH1F					*h_ResidualEnergy_CH{};
		TH1F					*h_ResidualTheta_CH{};
		TH1F					*h_ResidualPhi_CH{};
		TH1F					*h_ErrorEnergy_CH{};
		TH1F					*h_ErrorTheta_CH{};
		TH1F					*h_ErrorPhi_CH{};
		TH1F					*h_NormalizedResidualEnergy_CH{};
		TH1F					*h_NormalizedResidualTheta_CH{};
		TH1F					*h_NormalizedResidualPhi_CH{};
		TH1F					*h_ResidualTotalEnergy_NH{};
		TH1F					*h_ResidualTotalTheta_NH{};
		TH1F					*h_ResidualTotalPhi_NH{};
		TH1F					*h_ErrorTotalEnergy_NH{};
		TH1F					*h_ErrorTotalTheta_NH{};
		TH1F					*h_ErrorTotalPhi_NH{};
		TH1F					*h_NormalizedResidualTotalEnergy_NH{};
		TH1F					*h_NormalizedResidualTotalTheta_NH{};
		TH1F					*h_NormalizedResidualTotalPhi_NH{};
		TH1F					*h_ResidualTotalEnergy_R{};
		TH1F					*h_ResidualTotalTheta_R{};
		TH1F					*h_ResidualTotalPhi_R{};
		TH1F					*h_ErrorTotalEnergy_R{};
		TH1F					*h_ErrorTotalTheta_R{};
		TH1F					*h_ErrorTotalPhi_R{};
		TH1F					*h_NormalizedResidualTotalEnergy_R{};
		TH1F					*h_NormalizedResidualTotalTheta_R{};
		TH1F					*h_NormalizedResidualTotalPhi_R{};
		TH1F					*h_ResidualTotalEnergy{};
		TH1F					*h_ResidualTotalTheta{};
		TH1F					*h_ResidualTotalPhi{};
		TH1F					*h_ErrorTotalEnergy{};
		TH1F					*h_ErrorTotalTheta{};
		TH1F					*h_ErrorTotalPhi{};
		TH1F					*h_NormalizedResidualTotalEnergy{};
		TH1F					*h_NormalizedResidualTotalTheta{};
		TH1F					*h_NormalizedResidualTotalPhi{};
		TH2F					*h_NH_EclusterPlusMass_Emcp{};
		TH2F					*h_NHEnergy{};

};
#endif
