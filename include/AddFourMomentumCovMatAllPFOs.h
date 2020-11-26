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
		std::vector<float> UpdateNeutralPFOCovMatCartesian( TVector3 clusterPosition , float pfoEc , float pfoMass , std::vector<float> clusterPositionError , float clusterEnergyError );
		std::vector<float> UpdateNeutralPFOCovMatPolar( float pfoTheta , float pfoPhi , float pfoEnergy , float pfoMass , std::vector<float> PFOCoordinateError , float clusterEnergyError );
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

		bool					m_updateNormalNeutrals = true;
		bool					m_updateNeutrals_wTrack = true;
		bool					m_updateCharged = false;
		bool					m_AssumeNeutralPFOMassive = true;
		bool					m_isClusterEnergyKinEnergy = false;
		bool					m_useClusterPositionError = true;
		bool					m_updatePFO4Momentum = false;
		bool					m_useTrueJacobian = false;
		bool					m_scaleAngularUncertainty = false;
		float					m_AngularUncertaintyScaleFactor_NH = 1.0;
		float					m_AngularUncertaintyScaleFactor_Ph = 1.0;
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
	        TTree					*m_pTTree1;
	        TTree					*m_pTTree2;
	        TTree					*m_pTTree3;
	        TTree					*m_pTTree4;
	        TTree					*m_pTTree5;
		std::vector<int>			m_foundLinkedMCP{};
		std::vector<float>			m_mcEnergy{};
		std::vector<float>			m_mcTheta{};
		std::vector<float>			m_mcPhi{};
		std::vector<float>			m_RecoEnergy{};
		std::vector<float>			m_RecoTheta{};
		std::vector<float>			m_RecoPhi{};
		std::vector<float>			m_ResidualEnergy{};
		std::vector<float>			m_ResidualTheta{};
		std::vector<float>			m_ResidualPhi{};
		std::vector<float>			m_ErrorEnergy{};
		std::vector<float>			m_ErrorTheta{};
		std::vector<float>			m_ErrorPhi{};
		std::vector<float>			m_NormalizedResidualEnergy{};
		std::vector<float>			m_NormalizedResidualTheta{};
		std::vector<float>			m_NormalizedResidualPhi{};
		std::vector<int>			m_foundLinkedMCP_Ph{};
		std::vector<float>			m_mcEnergy_Ph{};
		std::vector<float>			m_mcTheta_Ph{};
		std::vector<float>			m_mcPhi_Ph{};
		std::vector<float>			m_RecoEnergy_Ph{};
		std::vector<float>			m_RecoTheta_Ph{};
		std::vector<float>			m_RecoPhi_Ph{};
		std::vector<float>			m_ResidualEnergy_Ph{};
		std::vector<float>			m_ResidualTheta_Ph{};
		std::vector<float>			m_ResidualPhi_Ph{};
		std::vector<float>			m_ErrorEnergy_Ph{};
		std::vector<float>			m_ErrorTheta_Ph{};
		std::vector<float>			m_ErrorPhi_Ph{};
		std::vector<float>			m_NormalizedResidualEnergy_Ph{};
		std::vector<float>			m_NormalizedResidualTheta_Ph{};
		std::vector<float>			m_NormalizedResidualPhi_Ph{};
		std::vector<int>			m_foundLinkedMCP_NH{};
		std::vector<float>			m_mcEnergy_NH{};
		std::vector<float>			m_mcTheta_NH{};
		std::vector<float>			m_mcPhi_NH{};
		std::vector<float>			m_RecoEnergy_NH{};
		std::vector<float>			m_RecoTheta_NH{};
		std::vector<float>			m_RecoPhi_NH{};
		std::vector<float>			m_ResidualEnergy_NH{};
		std::vector<float>			m_ResidualTheta_NH{};
		std::vector<float>			m_ResidualPhi_NH{};
		std::vector<float>			m_ErrorEnergy_NH{};
		std::vector<float>			m_ErrorTheta_NH{};
		std::vector<float>			m_ErrorPhi_NH{};
		std::vector<float>			m_NormalizedResidualEnergy_NH{};
		std::vector<float>			m_NormalizedResidualTheta_NH{};
		std::vector<float>			m_NormalizedResidualPhi_NH{};
		std::vector<int>			m_foundLinkedMCP_NH2T{};
		std::vector<float>			m_mcEnergy_NH2T{};
		std::vector<float>			m_mcTheta_NH2T{};
		std::vector<float>			m_mcPhi_NH2T{};
		std::vector<float>			m_RecoEnergy_NH2T{};
		std::vector<float>			m_RecoTheta_NH2T{};
		std::vector<float>			m_RecoPhi_NH2T{};
		std::vector<float>			m_ResidualEnergy_NH2T{};
		std::vector<float>			m_ResidualTheta_NH2T{};
		std::vector<float>			m_ResidualPhi_NH2T{};
		std::vector<float>			m_ErrorEnergy_NH2T{};
		std::vector<float>			m_ErrorTheta_NH2T{};
		std::vector<float>			m_ErrorPhi_NH2T{};
		std::vector<float>			m_NormalizedResidualEnergy_NH2T{};
		std::vector<float>			m_NormalizedResidualTheta_NH2T{};
		std::vector<float>			m_NormalizedResidualPhi_NH2T{};
		std::vector<int>			m_foundLinkedMCP_Ch{};
		std::vector<float>			m_mcEnergy_Ch{};
		std::vector<float>			m_mcTheta_Ch{};
		std::vector<float>			m_mcPhi_Ch{};
		std::vector<float>			m_RecoEnergy_Ch{};
		std::vector<float>			m_RecoTheta_Ch{};
		std::vector<float>			m_RecoPhi_Ch{};
		std::vector<float>			m_ResidualEnergy_Ch{};
		std::vector<float>			m_ResidualTheta_Ch{};
		std::vector<float>			m_ResidualPhi_Ch{};
		std::vector<float>			m_ErrorEnergy_Ch{};
		std::vector<float>			m_ErrorTheta_Ch{};
		std::vector<float>			m_ErrorPhi_Ch{};
		std::vector<float>			m_NormalizedResidualEnergy_Ch{};
		std::vector<float>			m_NormalizedResidualTheta_Ch{};
		std::vector<float>			m_NormalizedResidualPhi_Ch{};
	        TDirectory				*m_Histograms;
	        TDirectory				*m_CovMatElements;
	        TDirectory				*m_NeutralPFOswithoutTrak;
	        TDirectory				*m_NeutralPFOswith2Trak;
	        TDirectory				*m_ChargedPFOs;
	        TDirectory				*m_ErrorParameterization;
	        TDirectory				*m_Photon;
	        TDirectory				*m_NeutralPFO;
	        TDirectory				*m_ChargedPFO;
	        TDirectory				*m_CovMatinPolar;
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
		TH1F					*h_CovMatPolar_ThetaTheta{};
		TH1F					*h_CovMatPolar_ThetaPhi{};
		TH1F					*h_CovMatPolar_PhiPhi{};
		TH1F					*h_CovMatPolar_ThetaMomentum{};
		TH1F					*h_CovMatPolar_PhiMomentum{};
		TH1F					*h_CovMatPolar_MomentumMomentum{};
		TH1F					*h_CovMatPolar_ThetaEnergy{};
		TH1F					*h_CovMatPolar_PhiEnergy{};
		TH1F					*h_CovMatPolar_MomentumEnergy{};
		TH1F					*h_CovMatPolar_EnergyEnergy{};
		TH2F					*h_NH_EclusterPlusMass_Emcp{};
		TH2F					*h_NHEnergy{};

};
#endif
