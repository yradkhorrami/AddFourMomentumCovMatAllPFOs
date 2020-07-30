#include "TMath.h"
#include "streamlog/streamlog.h"
#include "marlin/Global.h"
#include <GeometryUtil.h>

#include "AddFourMomentumCovMatAllPFOs.h"

using namespace lcio;
using namespace marlin;

AddFourMomentumCovMatAllPFOs aAddFourMomentumCovMatAllPFOs;

AddFourMomentumCovMatAllPFOs::AddFourMomentumCovMatAllPFOs() :

Processor("AddFourMomentumCovMatAllPFOs"),
m_nRun(0),
m_nEvt(0),
m_nRunSum(0),
m_nEvtSum(0)

{
	_description = "AddFourMomentumCovMatAllPFOs calculates Four-Momentum CovMatrix for all PFOs (charged particles, neutral hadrons and photons)";

	registerInputCollection(	LCIO::RECONSTRUCTEDPARTICLE,
					"PfoCollection",
					"Name of input pfo collection",
					m_inputPfoCollection,
					std::string("PandoraPFOs")
				);

	registerOutputCollection(	LCIO::RECONSTRUCTEDPARTICLE,
					"CorrectedPfoCollection",
					"Name of output pfo collection",
					m_outputPfoCollection,
					std::string("CorrectedPfoCollection")
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
	this->Clear();
	c = 2.99792458e8;
	mm2m = 1e-3;
	eV2GeV = 1e-9;
	eB = m_Bfield * c * mm2m * eV2GeV;
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
		n_PFO = inputPfoCollection->getNumberOfElements();
		LCCollectionVec *m_col_outputPfo = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
        }
        catch(DataNotAvailableException &e)
        {
          streamlog_out(MESSAGE) << "Input collection not found in event " << m_nEvt << std::endl;
        }
}

void AddFourMomentumCovMatAllPFOs::Clear()
{

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

	std::cout << " END : processed events: " << m_nEvtSum << std::endl;

}





