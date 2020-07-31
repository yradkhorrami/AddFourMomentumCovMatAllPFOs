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
#include "lcio.h"
#include "TMatrixD.h"
#include <string>
#include <vector>
#include <math.h>
#include <set>
#include <vector>

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
		virtual void processRunHeader();
		virtual void processEvent( EVENT::LCEvent *pLCEvent );
		virtual void check( EVENT::LCEvent *pLCEvent );
		virtual void end();

	private:
		std::string				m_inputPfoCollection{};
		std::string				m_outputPfoCollection{};

		int					m_nRun;
		int					m_nEvt;
		int					m_nRunSum;
		int					m_nEvtSum;


#endif
