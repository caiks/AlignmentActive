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

ActiveEventsRepa::ActiveEventsRepa() : references(0)
{
}

ActiveEventsArray::ActiveEventsArray() : references(0)
{
}

// Active::Active() : terminate(false), log(log_default), historyOverflow(false), historyEvent(0), applicationFudId(0), applicationFudIdPersistent(0)
Active::Active() : terminate(false), log(log_default), historyOverflow(false), historyEvent(0), historySize(0)
{
}

// std::size_t Alignment::Active::size() const
// {
	// return this->historyOverflow ? this->history->size : this->historyEvent + 1;
// }

#define UNLOG ; log_str.flush(); this->log(log_str.str());}
#define LOG { std::ostringstream log_str; log_str <<

// bool Alignment::Active::report()
// {
	// bool ok = true;
	// if (this->terminate)
		// return true;
	
	// try {
		// std::lock_guard<std::recursive_mutex> guard(this->mutex);
		
		// auto t0 = clk::now();
		// LOG "report\tterminate: " << (this->terminate ? "true" : "false") UNLOG
		// LOG "report\tunderlying system size: " << (this->systemUnder ? this->systemUnder->listVarSizePair.size() : 0) UNLOG
		// LOG "report\tsystem size: " << (this->system ? this->system->listVarSizePair.size() : 0) UNLOG		
		// LOG "report\thistory dimension: " << (this->history ? this->history->dimension : 0) UNLOG		
		// LOG "report\thistory size: " << (this->history ? this->history->size : 0) UNLOG	
		// LOG "report\thistory overflow: " << (this->history && this->historyOverflow ? "true" : "false") UNLOG	
		// LOG "report\thistory event: " << (this->history ? this->historyEvent : 0) UNLOG	
		// LOG "report\tunderlying model size: " << (this->applicationUnder ? fudRepasSize(*this->applicationUnder->fud) : 0)  UNLOG	
		// LOG "report\tunderlying model underlying size: " << (this->applicationUnder ? fudRepasUnderlying(*this->applicationUnder->fud)->size() : 0)  UNLOG	
		// LOG "report\tunderlying model slices size: " << (this->applicationUnder ? treesSize(*this->applicationUnder->slices) : 0)  UNLOG	
		// LOG "report\tunderlying model leaf slices size: " << (this->applicationUnder ? treesLeafElements(*this->applicationUnder->slices)->size() : 0)  UNLOG
		// LOG "report\tmodel size: " << (this->application ? fudRepasSize(*this->application->fud) : 0)  UNLOG	
		// LOG "report\tmodel underlying size: " << (this->application ? fudRepasUnderlying(*this->application->fud)->size() : 0)  UNLOG	
		// LOG "report\tmodel slices size: " << (this->application ? treesSize(*this->application->slices) : 0)  UNLOG	
		// LOG "report\tmodel leaf slices size: " << (this->application ? treesLeafElements(*this->application->slices)->size() : 0)  UNLOG
		// LOG "report\tpersistent fud id: " << this->applicationFudIdPersistent  UNLOG
		// LOG "report\tfud id: " << this->applicationFudId  UNLOG	
		// LOG "report\teventsSlice size: " << this->eventsSlice.size()  UNLOG	
		// LOG "report\tslicesSetEvent size: " << this->slicesSetEvent.size()  UNLOG	
		// LOG "report\tsetSliceLeaf size: " << this->setSliceLeaf.size()  UNLOG	

		// auto t1 = clk::now();
		// LOG "report\ttime: " << ((sec)(t1 - t0)).count() << "s" UNLOG	
	// } 
	// catch (const std::exception& e) 
	// {
		// LOG "report\terror: " << e.what()  UNLOG
		// ok = false;
	// }
	
	// return ok;
// }

// bool Alignment::Active::slicesSync(std::size_t tint)
// {
	// auto hrsel = eventsHistoryRepasHistoryRepaSelection_u;
	// auto frmul = historyRepasFudRepasMultiply_up;

	// bool ok = true;
	// if (this->terminate)
		// return true;
	
	// try {
		// std::lock_guard<std::recursive_mutex> guard(this->mutex);
		// if (ok && (!this->history || !this->application))
		// {
			// LOG "slicesSync\terror: active is not set" UNLOG
			// ok = false;
		// }
		// if (ok && this->history->evient)
		// {
			// LOG "slicesSync\terror: active history is evient" UNLOG
			// ok = false;
		// }
		// auto t0 = clk::now();
		// if (ok)
		// {
			// if (this->eventsSlice.size() != this->history->size)
				// this->eventsSlice.resize(this->history->size);
		// }
		// SizeList listEvent;
		// SizeSet setSlice;	
		// if (ok)
		// {
			// auto mark = clk::now();
			// auto listSliceLeaf = treesLeafElements(*this->application->slices);
			// this->setSliceLeaf.clear();
			// this->setSliceLeaf.insert(listSliceLeaf->begin(),listSliceLeaf->end());
			// auto activeSize = this->size();
			// for (std::size_t i = 0; i < activeSize; i++)
			// {
				// auto sl = this->eventsSlice[i];
				// if (!sl || this->setSliceLeaf.find(sl) == this->setSliceLeaf.end())
				// {
					// listEvent.push_back(i);
					// setSlice.insert(sl);		
				// }
			// }
			// LOG "slicesSync\tslice leaf size: " << this->setSliceLeaf.size() 
				// << "\ttime: " << ((sec)(clk::now() - mark)).count() << "s" UNLOG	
		// }
		// std::unique_ptr<HistoryRepa> historySelection;
		// if (ok)
		// {
			// auto mark = clk::now();
			// historySelection = hrsel(listEvent.size(), listEvent.data(), *this->history);
			// LOG "slicesSync\tselection size: " << historySelection->size 
				// << "\ttime: " << ((sec)(clk::now() - mark)).count() << "s" UNLOG	
		// }
		// if (ok)
		// {
			// auto mark = clk::now();
			// if (this->applicationUnder)
				// historySelection = frmul(tint, *historySelection, *this->applicationUnder->fud);
			// historySelection = frmul(tint, *historySelection, *this->application->fud);
			// LOG "slicesSync\tapplication time: " << ((sec)(clk::now() - mark)).count() << "s" UNLOG	
		// }
		// if (ok)
		// {			
			// auto mark = clk::now();
			// auto& hr = historySelection;
			// auto z = hr->size;
			// auto m = this->setSliceLeaf.size();
			// auto& mvv = hr->mapVarInt();
			// auto rr = hr->arr;
			// SizeList pvv;
			// pvv.reserve(m);
			// for (auto v : this->setSliceLeaf)
				// pvv.push_back(mvv[v]);
			// std::size_t eventCount = 0;
			// for (std::size_t j = 0; j < z; j++)
			// {
				// std::size_t i = 0;
				// for (auto sl : this->setSliceLeaf)
				// {
					// std::size_t u = rr[pvv[i]*z + j];
					// if (u)
					// {
						// auto ev = listEvent[j];
						// this->eventsSlice[ev] = sl;
						// this->slicesSetEvent[sl].insert(ev);
						// eventCount++;
						// break;
					// }
					// i++;
				// }
			// }
			// if (eventCount < z)
			// {
				// LOG "slicesSync\twarning: only " << eventCount << " events updated" UNLOG				
			// }
			// LOG "slicesSync\tupdate time: " << ((sec)(clk::now() - mark)).count() << "s" UNLOG	
		// }
		// if (ok)
		// {			
			// auto mark = clk::now();
			// for (auto sl : setSlice)
				// this->slicesSetEvent.erase(sl);
			// LOG "slicesSync\ttidy time: " << ((sec)(clk::now() - mark)).count() << "s" UNLOG		
		// }
		// if (ok)
		// {	
			// LOG "slicesSync\ttime: " << ((sec)(clk::now() - t0)).count() << "s" UNLOG	
		// }
	// } 
	// catch (const std::exception& e) 
	// {
		// LOG "slicesSync error: " << e.what()  UNLOG
		// ok = false;
	// }
	
	// return ok;
// }

// bool Alignment::Active::slicesUpdate(std::size_t tint)
// {
	// auto hrsel = eventsHistoryRepasHistoryRepaSelection_u;
	// auto frmul = historyRepasFudRepasMultiply_up;

	// bool ok = true;
	// if (this->terminate)
		// return true;
	
	// try {
		// std::lock_guard<std::recursive_mutex> guard(this->mutex);
		// if (ok && (!this->history || !this->application))
		// {
			// LOG "slicesUpdate\terror: active is not set" UNLOG
			// ok = false;
		// }
		// auto t0 = clk::now();
		// if (ok && this->historyEvent >= this->history->size)
		// {
			// LOG "slicesUpdate\terror: invalid historyEvent " << this->historyEvent UNLOG
			// ok = false;
		// }
		// std::unique_ptr<HistoryRepa> historySelection;
		// if (ok)
		// {
			// SizeList listEvent;
			// listEvent.push_back(this->historyEvent);
			// historySelection = hrsel(listEvent.size(), listEvent.data(), *this->history);
			// LOG "slicesUpdate\tselection size: " << historySelection->size 
				// << "\ttime: " << ((sec)(clk::now() - t0)).count() << "s" UNLOG	
		// }
		// if (ok)
		// {
			// auto mark = clk::now();
			// if (this->applicationUnder)
				// historySelection = frmul(tint, *historySelection, *this->applicationUnder->fud);
			// historySelection = frmul(tint, *historySelection, *this->application->fud);
			// LOG "slicesUpdate\tapplication time: " << ((sec)(clk::now() - mark)).count() << "s" UNLOG	
		// }
		// if (ok)
		// {			
			// auto mark = clk::now();
			// auto ev = this->historyEvent;
			// this->slicesSetEvent[this->eventsSlice[ev]].erase(ev);
			// this->eventsSlice[ev] = 0;
			// auto& hr = historySelection;
			// auto m = this->setSliceLeaf.size();
			// auto& mvv = hr->mapVarInt();
			// auto rr = hr->arr;
			// SizeList pvv;
			// pvv.reserve(m);
			// for (auto v : this->setSliceLeaf)
				// pvv.push_back(mvv[v]);
			// std::size_t eventCount = 0;
			// std::size_t i = 0;
			// for (auto sl : this->setSliceLeaf)
			// {
				// std::size_t u = rr[pvv[i]];
				// if (u)
				// {
					// this->eventsSlice[ev] = sl;
					// this->slicesSetEvent[sl].insert(ev);
					// eventCount++;
					// break;
				// }
				// i++;
			// }
			// if (!eventCount)
			// {
				// LOG "slicesUpdate\twarning: no events updated" UNLOG				
			// }
			// else
			// {
				// LOG "slicesUpdate\tevent: " << ev << " slice: " << this->eventsSlice[ev] UNLOG				
			// }
			// LOG "slicesUpdate\tupdate time: " << ((sec)(clk::now() - mark)).count() << "s" UNLOG	
		// }
		// if (ok)
		// {	
			// LOG "slicesUpdate\ttime: " << ((sec)(clk::now() - t0)).count() << "s" UNLOG	
		// }
	// } 
	// catch (const std::exception& e) 
	// {
		// LOG "slicesUpdate error: " << e.what()  UNLOG
		// ok = false;
	// }
	
	// return ok;
// }
