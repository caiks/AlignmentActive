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

#define UNLOG ; log_str.flush(); log(log_str.str());}
#define LOG { std::ostringstream log_str; log_str <<

Active::Active() : terminate(false), historyOverflow(false), historyEvent(0), applicationFudId(0), applicationFudIdPersistent(0)
{
}

void Alignment::activesLog(Active& act, void (*log)(const std::string&))
{
	std::lock_guard<std::mutex> guard(act.mutex);
	
	LOG "active terminate: " << (act.terminate ? "true" : "false") UNLOG
	LOG "active underlying system size: " << (act.systemUnder ? act.systemUnder->listVarSizePair.size() : 0) UNLOG
	LOG "active system size: " << (act.system ? act.system->listVarSizePair.size() : 0) UNLOG		
	LOG "active history dimension: " << (act.history ? act.history->dimension : 0) UNLOG		
	LOG "active history size: " << (act.history ? act.history->size : 0) UNLOG	
	LOG "active history overflow: " << (act.history && act.historyOverflow ? "true" : "false") UNLOG	
	LOG "active history event: " << (act.history ? act.historyEvent : 0) UNLOG	
	LOG "active underlying model size: " << (act.applicationUnder ? fudRepasSize(*act.applicationUnder->fud) : 0)  UNLOG	
	LOG "active underlying model underlying size: " << (act.applicationUnder ? fudRepasUnderlying(*act.applicationUnder->fud)->size() : 0)  UNLOG	
	LOG "active underlying model slices size: " << (act.applicationUnder ? treesSize(*act.applicationUnder->slices) : 0)  UNLOG	
	LOG "active underlying model leaf slices size: " << (act.applicationUnder ? treesLeafElements(*act.applicationUnder->slices)->size() : 0)  UNLOG
	LOG "active model size: " << (act.application ? fudRepasSize(*act.application->fud) : 0)  UNLOG	
	LOG "active model underlying size: " << (act.application ? fudRepasUnderlying(*act.application->fud)->size() : 0)  UNLOG	
	LOG "active model slices size: " << (act.application ? treesSize(*act.application->slices) : 0)  UNLOG	
	LOG "active model leaf slices size: " << (act.application ? treesLeafElements(*act.application->slices)->size() : 0)  UNLOG
	LOG "active persistent fud id: " << act.applicationFudIdPersistent  UNLOG
	LOG "active fud id: " << act.applicationFudId  UNLOG
}

