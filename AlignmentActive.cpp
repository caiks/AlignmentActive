#include "AlignmentActive.h"

#include "AlignmentApprox.h"
#include "AlignmentAeson.h"
#include "AlignmentAesonRepa.h"
#include "AlignmentRandomRepa.h"
#include "AlignmentPracticableRepa.h"
#include "AlignmentPracticableIORepa.h"

#include <stdlib.h>
#include <sstream>
#include <fstream>
#include <chrono>
#include <ctime>

using namespace Alignment;


void log_default(const std::string& str)
{
	std::cout << str << std::endl;
	return;
};

Active::Active() : terminate(false), log_f(log_default), historyOverflow(false), historyEvent(0), applicationFudId(0), applicationFudIdPersistent(0)
{
}

#define UNLOG ; log_str.flush(); this->log_f(log_str.str());}
#define LOG { std::ostringstream log_str; log_str <<

bool Alignment::Active::log()
{
	bool ok = true;
	if (this->terminate)
		return true;
	std::lock_guard<std::recursive_mutex> guard(this->mutex);
	
	try {
		
		LOG "active terminate: " << (this->terminate ? "true" : "false") UNLOG
		LOG "active underlying system size: " << (this->systemUnder ? this->systemUnder->listVarSizePair.size() : 0) UNLOG
		LOG "active system size: " << (this->system ? this->system->listVarSizePair.size() : 0) UNLOG		
		LOG "active history dimension: " << (this->history ? this->history->dimension : 0) UNLOG		
		LOG "active history size: " << (this->history ? this->history->size : 0) UNLOG	
		LOG "active history overflow: " << (this->history && this->historyOverflow ? "true" : "false") UNLOG	
		LOG "active history event: " << (this->history ? this->historyEvent : 0) UNLOG	
		LOG "active underlying model size: " << (this->applicationUnder ? fudRepasSize(*this->applicationUnder->fud) : 0)  UNLOG	
		LOG "active underlying model underlying size: " << (this->applicationUnder ? fudRepasUnderlying(*this->applicationUnder->fud)->size() : 0)  UNLOG	
		LOG "active underlying model slices size: " << (this->applicationUnder ? treesSize(*this->applicationUnder->slices) : 0)  UNLOG	
		LOG "active underlying model leaf slices size: " << (this->applicationUnder ? treesLeafElements(*this->applicationUnder->slices)->size() : 0)  UNLOG
		LOG "active model size: " << (this->application ? fudRepasSize(*this->application->fud) : 0)  UNLOG	
		LOG "active model underlying size: " << (this->application ? fudRepasUnderlying(*this->application->fud)->size() : 0)  UNLOG	
		LOG "active model slices size: " << (this->application ? treesSize(*this->application->slices) : 0)  UNLOG	
		LOG "active model leaf slices size: " << (this->application ? treesLeafElements(*this->application->slices)->size() : 0)  UNLOG
		LOG "active persistent fud id: " << this->applicationFudIdPersistent  UNLOG
		LOG "active fud id: " << this->applicationFudId  UNLOG	
	} 
	catch (const std::exception& e) 
	{
		LOG "active log error: " << e.what()  UNLOG
		ok = false;
	}
	
	return ok;
}

bool Alignment::Active::slicesSync()
{
	bool ok = true;
	if (this->terminate)
		return true;
	std::lock_guard<std::recursive_mutex> guard(this->mutex);
	
	try {

	} 
	catch (const std::exception& e) 
	{
		LOG "active slicesSync error: " << e.what()  UNLOG
		ok = false;
	}
	
	return ok;
}
