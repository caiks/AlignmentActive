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

typedef std::chrono::duration<double> sec;
typedef std::chrono::high_resolution_clock clk;

void log_default(const std::string& str)
{
	std::cout << str << std::endl;
	return;
};

Active::Active() : terminate(false), log(log_default), historyOverflow(false), historyEvent(0), applicationFudId(0), applicationFudIdPersistent(0)
{
}

#define UNLOG ; log_str.flush(); this->log(log_str.str());}
#define LOG { std::ostringstream log_str; log_str <<

bool Alignment::Active::report()
{
	bool ok = true;
	if (this->terminate)
		return true;
	
	try {
		std::lock_guard<std::recursive_mutex> guard(this->mutex);
		
		auto t0 = clk::now();
		LOG "report\tterminate: " << (this->terminate ? "true" : "false") UNLOG
		LOG "report\tunderlying system size: " << (this->systemUnder ? this->systemUnder->listVarSizePair.size() : 0) UNLOG
		LOG "report\tsystem size: " << (this->system ? this->system->listVarSizePair.size() : 0) UNLOG		
		LOG "report\thistory dimension: " << (this->history ? this->history->dimension : 0) UNLOG		
		LOG "report\thistory size: " << (this->history ? this->history->size : 0) UNLOG	
		LOG "report\thistory overflow: " << (this->history && this->historyOverflow ? "true" : "false") UNLOG	
		LOG "report\thistory event: " << (this->history ? this->historyEvent : 0) UNLOG	
		LOG "report\tunderlying model size: " << (this->applicationUnder ? fudRepasSize(*this->applicationUnder->fud) : 0)  UNLOG	
		LOG "report\tunderlying model underlying size: " << (this->applicationUnder ? fudRepasUnderlying(*this->applicationUnder->fud)->size() : 0)  UNLOG	
		LOG "report\tunderlying model slices size: " << (this->applicationUnder ? treesSize(*this->applicationUnder->slices) : 0)  UNLOG	
		LOG "report\tunderlying model leaf slices size: " << (this->applicationUnder ? treesLeafElements(*this->applicationUnder->slices)->size() : 0)  UNLOG
		LOG "report\tmodel size: " << (this->application ? fudRepasSize(*this->application->fud) : 0)  UNLOG	
		LOG "report\tmodel underlying size: " << (this->application ? fudRepasUnderlying(*this->application->fud)->size() : 0)  UNLOG	
		LOG "report\tmodel slices size: " << (this->application ? treesSize(*this->application->slices) : 0)  UNLOG	
		LOG "report\tmodel leaf slices size: " << (this->application ? treesLeafElements(*this->application->slices)->size() : 0)  UNLOG
		LOG "report\tpersistent fud id: " << this->applicationFudIdPersistent  UNLOG
		LOG "report\tfud id: " << this->applicationFudId  UNLOG	
		auto t1 = clk::now();
		LOG "report\ttime: " << ((sec)(t1 - t0)).count() << "s" UNLOG	
	} 
	catch (const std::exception& e) 
	{
		LOG "report\terror: " << e.what()  UNLOG
		ok = false;
	}
	
	return ok;
}

bool Alignment::Active::slicesSync(std::size_t tint)
{
	auto hrsel = eventsHistoryRepasHistoryRepaSelection_u;
	auto frmul = historyRepasFudRepasMultiply_up;

	bool ok = true;
	if (this->terminate)
		return true;
	
	try {
		std::lock_guard<std::recursive_mutex> guard(this->mutex);
		if (!this->history || !this->application)
			return ok;
		auto t0 = clk::now();
		if (this->eventsSlice.size() != this->history->size)
			this->eventsSlice.resize(this->history->size);
		auto listSliceLeaf = treesLeafElements(*this->application->slices);
		SizeSet setSliceLeaf(listSliceLeaf->begin(),listSliceLeaf->end());
		SizeList listEvent;
		SizeSet setSlice;
		auto activeSize = this->historyOverflow ? this->history->size : this->historyEvent + 1;
		for (std::size_t i = 0; i < activeSize; i++)
		{
			auto sl = this->eventsSlice[i];
			if (!sl || setSliceLeaf.find(sl) == setSliceLeaf.end())
			{
				listEvent.push_back(i);
				setSlice.insert(sl);		
			}
		}
		auto hr = hrsel(listEvent.size(), listEvent.data(), *this->history);
		auto t1 = clk::now();
		LOG "slicesSync\tselection size: " << hr->size << "\ttime: " << ((sec)(t1 - t0)).count() << "s" UNLOG	
		
		auto m = listSliceLeaf->size();
		
	} 
	catch (const std::exception& e) 
	{
		LOG "slicesSync error: " << e.what()  UNLOG
		ok = false;
	}
	
	return ok;
}
