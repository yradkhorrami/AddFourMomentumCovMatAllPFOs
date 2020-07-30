#ifndef AddFourMomentumCovMatAllPFOs_h
#define AddFourMomentumCovMatAllPFOs_h 1


#include "marlin/Processor.h"
#include "lcio.h"
#include <EVENT/ReconstructedParticle.h>
#include "TMatrixD.h"

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
		PFOCorrection(const AddFourMomentumCovMatAllPFOs&) = delete;
		AddFourMomentumCovMatAllPFOs& operator=(const AddFourMomentumCovMatAllPFOs&) = delete;
		virtual void init();
		virtual void processRunHeader();
		virtual void processEvent( EVENT::LCEvent *pLCEvent );
		virtual void check( EVENT::LCEvent *pLCEvent );
		virtual void end();


#endif
