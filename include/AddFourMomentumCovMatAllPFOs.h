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
		std::vector<float> UpdateNeutralPFOCovMatDirError( TLorentzVector pfoFourMomentum , std::vector<float> clusterDirectionError , float clusterEnergyError );
		std::vector<float> UpdateNeutralPFOCovMatPosError( TVector3 clusterPosition , float pfoEnergy , float pfoMass , std::vector<float> clusterPositionError , float clusterEnergyError );
		std::vector<float> UpdateChargedPFOCovMat( EVENT::Track* MyTrack , float trackMass );
		double getTrackMass( EVENT::LCEvent *pLCEvent, EVENT::Track* inputTrk );
		virtual void check( EVENT::LCEvent *pLCEvent );
		virtual void end();

	private:
		std::string				m_inputPfoCollection{};
		std::string				m_outputPfoCollection{};
		std::string				m_rootFile{};
		std::string				m_MarlinTrkTracks{};
		std::string				m_MarlinTrkTracksKAON{};
		std::string				m_MarlinTrkTracksPROTON{};

		bool					m_useClusterPositionError = true;
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
		TH2I					*h_nTracks_PFOCharge{};
		TH2I					*h_nClusters_nTracks{};
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

};
#endif
