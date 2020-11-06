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

ActiveSystem::ActiveSystem() : blockBits(16), blockCurrent(1)
{
}

ActiveEventsRepa::ActiveEventsRepa(std::size_t referencesA) : references(referencesA)
{
}

std::ostream& operator<<(std::ostream& out, const ActiveEventsRepa& ev)
{
	out << "(" << ev.references << ",[";
	bool first = true;
	for (auto& pp : ev.mapIdEvent)
	{
		if (!first)
			out << ",";
		first = false;
		out << "(" << pp.first << "," << *pp.second.first << "," << pp.second.second << ")";
	}
	out << "])";
	return out;
}

ActiveEventsArray::ActiveEventsArray(std::size_t referencesA) : references(referencesA)
{
}

std::ostream& operator<<(std::ostream& out, const ActiveEventsArray& ev)
{
	out << "(" << ev.references << ",[";
	bool first = true;
	for (auto& pp : ev.mapIdEvent)
	{
		if (!first)
			out << ",";
		first = false;
		out << "(" << pp.first << "," << *pp.second.first << "," << pp.second.second << ")";
	}
	out << "])";
	return out;
}

Active::Active() : terminate(false), log(log_default), historyOverflow(false), historyEvent(0), historySize(0), blockBits(16), blockCurrent(1), induceThreshold(0)
{
}

#define UNLOG ; log_str.flush(); this->log(log_str.str());}
#define LOG { std::ostringstream log_str; log_str <<

// event ids should be monotonic and updated no more than once
bool Alignment::Active::update(std::size_t mapCapacity)
{
	auto drmul = historyRepaPtrListsHistorySparseArrayPtrListsDecompFudSlicedRepasEventsPathSlice_u;

	bool ok = true;
	if (this->terminate)
		return true;
	try 
	{
		while (ok)
		{
			std::lock_guard<std::mutex> guard(this->mutex);
			// check consistent underlying
			if (ok)
			{
				for (auto& ev : this->underlyingEventsRepa)
					ok = ok && ev;
				for (auto& ev : this->underlyingEventsSparse)
					ok = ok && ev;
				for (auto& hr : this->underlyingHistoryRepa)
					ok = ok && hr && hr->size == historySize && hr->dimension > 0;
				for (auto& hr : this->underlyingHistorySparse)
					ok = ok && hr && hr->size == historySize && hr->capacity == 1;
				ok = ok && this->underlyingEventsRepa.size() == this->underlyingHistoryRepa.size();
				ok = ok && this->underlyingEventsSparse.size() == this->underlyingHistorySparse.size();
				if (!ok)
				{
					LOG "Active::update\terror: inconsistent underlying" UNLOG
					break;
				}	
			}
			// find set of intersecting event ids
			SizeSet eventsA;
			SizeSet eventsUpdatedA(this->eventsUpdated);
			if (ok)
			{
				bool eventsFirst = true;
				for (auto& ev : this->underlyingEventsRepa)
				{
					if (!eventsFirst && !eventsA.size())
						break;				
					std::lock_guard<std::mutex> guard(ev->mutex);
					SizeSet eventsB;
					for (auto& qq : ev->mapIdEvent)	
					{
						if (eventsUpdatedA.find(qq.first) == eventsUpdatedA.end() 
							&& (eventsFirst || eventsA.find(qq.first) != eventsA.end()))
							eventsB.insert(qq.first);						
					}	
					eventsFirst = false;		
					eventsA = eventsB;
				}
				for (auto& ev : this->underlyingEventsSparse)
				{
					if (!eventsFirst && !eventsA.size())
						break;				
					std::lock_guard<std::mutex> guard(ev->mutex);
					SizeSet eventsB;
					for (auto& qq : ev->mapIdEvent)	
					{
						if (eventsUpdatedA.find(qq.first) == eventsUpdatedA.end() 
							&& (eventsFirst || eventsA.find(qq.first) != eventsA.end()))
							eventsB.insert(qq.first);						
					}	
					eventsFirst = false;		
					eventsA = eventsB;
				}
			}		
			// if there is an event then process it otherwise return
			if (!(ok && eventsA.size() && (!eventsUpdatedA.size() || *eventsA.end() > *eventsUpdatedA.end()))) 
				break;
			// get first new event id
			auto eventA = *eventsA.begin();
			if (ok && eventsUpdatedA.size())
			{
				for (auto eventB : eventsA)
					if (eventB > *eventsUpdatedA.end())
					{
						eventA = eventB;
						break;
					}
			}
			// copy events to active history
			if (ok)
			{		
				HistoryRepaPtrList hrs;
				hrs.reserve(this->underlyingEventsRepa.size());
				for (auto& ev : this->underlyingEventsRepa)
					if (ok)
					{
						std::lock_guard<std::mutex> guard(ev->mutex);
						auto it = ev->mapIdEvent.find(eventA);
						ok = ok && it != ev->mapIdEvent.end() && it->second.first;
						if (ok)
							hrs.push_back(it->second.first);
						else
						{
							LOG "Active::update\terror: lost event repa " << eventA UNLOG
						}									
					}
				HistorySparseArrayPtrList has;
				has.reserve(this->underlyingEventsSparse.size());		
				for (auto& ev : this->underlyingEventsSparse)
					if (ok)
					{
						std::lock_guard<std::mutex> guard(ev->mutex);
						auto it = ev->mapIdEvent.find(eventA);
						ok = ok && it != ev->mapIdEvent.end() && it->second.first;
						if (ok)
							has.push_back(it->second.first);
						else
						{
							LOG "Active::update\terror: lost event sparse" << eventA UNLOG
						}									
					}
				if (!ok)
					break;					
				ok = ok && this->historyEvent < this->historySize;
				if (!ok) 
				{
					LOG "Active::update\terror: inconsistent historyEvent " << this->historyEvent << " compared to historySize " << this->historySize UNLOG		
					break;
				}
				ok = ok && hrs.size() == this->underlyingHistoryRepa.size();
				for (std::size_t h = 0; ok && h < hrs.size(); h++)
					ok = ok && hrs[h]->size == 1;
				ok = ok && has.size() == this->underlyingHistorySparse.size();
				for (std::size_t h = 0; ok && h < has.size(); h++)
					ok = ok && has[h]->size == 1;
				if (!ok) 
				{
					LOG "Active::update\terror: inconsistent underlying " UNLOG		
					break;
				}
				for (std::size_t h = 0; ok && h < hrs.size(); h++)
				{
					auto& hr = *this->underlyingHistoryRepa[h];
					auto z = hr.size;
					auto n = hr.dimension;
					auto vv = hr.vectorVar;
					auto rr = hr.arr;
					auto& hr1 = *hrs[h];
					auto n1 = hr1.dimension;
					auto vv1 = hr1.vectorVar;
					auto rr1 = hr1.arr;
					auto j = this->historyEvent;
					bool equiv = n == n1;
					for (std::size_t i = 0; equiv && i < n; i++)
						equiv = equiv && vv[i] == vv1[i];
					if (equiv)
					{
						if (hr.evient)
						{
							std::size_t jn = j*n;
							for (std::size_t i = 0; i < n; i++)
								rr[jn + i] = rr1[i];
						}
						else
						{
							for (std::size_t i = 0; i < n; i++)
								rr[i*z + j] = rr1[i];
						}							
					}
					else
					{
						auto& mvv = hr.mapVarInt();
						if (hr.evient)
						{
							std::size_t jn = j*n;
							for (std::size_t i = 0; i < n; i++)
								rr[jn + i] = 0;
							for (std::size_t i = 0; i < n1; i++)
							{
								auto v = vv1[i];
								auto it = mvv.find(v);
								if (it != mvv.end())
									rr[jn + it->second] = rr1[i];
							}
						}
						else
						{
							for (std::size_t i = 0; i < n; i++)
								rr[i*z + j] = 0;
							for (std::size_t i = 0; i < n1; i++)
							{
								auto v = vv1[i];
								auto it = mvv.find(v);
								if (it != mvv.end())
									rr[it->second * z + j] = rr1[i];
							}
						}								
					}
				}
				for (std::size_t h = 0; ok && h < has.size(); h++)
				{
					auto& hr = *this->underlyingHistorySparse[h];
					auto rr = hr.arr;
					auto hr1 = has[h];
					auto n = hr1->capacity;
					auto rr1 = hr1->arr;
					auto j = this->historyEvent;
					std::size_t v = 0;
					for (std::size_t i = 0; i < n; i++)
						if (rr1[i])
							v = rr1[i];
					rr[j] = v;
					if (v && this->slicesPath.find(v) == this->slicesPath.end())
						this->slicesPath[v] = hr1;
				}
			}
			if (ok)
			{
				ok = ok && this->decomp;
				if (!ok)
				{
					LOG "Active::update\terror: no decomp set" UNLOG
				}				
			}
			// apply the model
			if (ok)
			{
				auto ll = drmul(this->underlyingHistoryRepa,this->underlyingHistorySparse,*this->decomp,this->historyEvent,mapCapacity);	
				ok = ok && ll;
				if (!ok)
				{
					LOG "Active::update\terror: drmul failed to return a list" UNLOG
				}
				// sync active slices
				if (ok)
				{
					std::size_t sliceA = 0;
					if (ll->size())
						sliceA = ll->back();
					if (this->eventsSlice.size() != this->historySize)
						this->eventsSlice.resize(this->historySize);
					std::size_t sliceB = this->eventsSlice[this->historyEvent];
					if (!sliceA || sliceA != sliceB)
					{
						this->eventsSlice[this->historyEvent] = sliceA;
						auto& setA = this->slicesSetEvent[sliceA];
						setA.insert(this->historyEvent);
						if (this->induceThreshold && setA.size() == this->induceThreshold)
							this->slicesInduce.insert(sliceA);
						if (sliceA)
						{
							auto& setB = this->slicesSetEvent[sliceB];
							setB.erase(this->historyEvent);
							if (this->induceThreshold && setB.size() == this->induceThreshold-1)
								this->slicesInduce.erase(sliceB);
						}
					}					
				}
				// create overlying event
				if (ok && this->eventsSparse)
				{
					auto& ev = this->eventsSparse;
					std::lock_guard<std::mutex> guard(ev->mutex);
					if (ev->references)
					{
						std::size_t n = ll->size() ? ll->size() : 1;
						auto hr = std::make_shared<HistorySparseArray>(1,n);
						if (ll->size())
						{
							auto rr = hr->arr;	
							for (std::size_t i = 0; i < n; i++)						
								rr[i] = (*ll)[i];
						}
						ev->mapIdEvent.insert_or_assign(eventA,std::pair<HistorySparseArrayPtr,std::size_t>(hr,ev->references));
					}			
				}
			}
			// increment historyEvent
			if (ok)
			{
				this->historyEvent++;
				if (this->historyEvent >= this->historySize)
				{
					this->historyEvent = 0;
					this->historyOverflow = true;
				}
				for (auto eventB : eventsA)
					if (eventB <= eventA)
						this->eventsUpdated.insert(eventB);
			}
			// unreference underlying up to update event
			std::size_t eventLeast = eventA;
			if (ok)
			{
				for (auto& ev : this->underlyingEventsRepa)
				{			
					std::lock_guard<std::mutex> guard(ev->mutex);
					if (ev->mapIdEvent.size() && ev->mapIdEvent.begin()->first < eventLeast)
						eventLeast = ev->mapIdEvent.begin()->first;
					SizeSet eventsB;
					for (auto& qq : ev->mapIdEvent)	
					{
						if (qq.first <= eventA 
							&& eventsA.find(qq.first) != eventsA.end()
							&& eventsUpdatedA.find(qq.first) == eventsUpdatedA.end())
						{
							auto& refs = qq.second.second;
							if (refs)
								refs--;									
							if (!refs)
								eventsB.insert(qq.first);									
						}
					}	
					for (auto eventB : eventsB)	
						ev->mapIdEvent.erase(eventB);
				}
				for (auto& ev : this->underlyingEventsSparse)
				{			
					std::lock_guard<std::mutex> guard(ev->mutex);
					if (ev->mapIdEvent.size() && ev->mapIdEvent.begin()->first < eventLeast)
						eventLeast = ev->mapIdEvent.begin()->first;
					SizeSet eventsB;
					for (auto& qq : ev->mapIdEvent)	
					{
						if (qq.first <= eventA 
							&& eventsA.find(qq.first) != eventsA.end()
							&& eventsUpdatedA.find(qq.first) == eventsUpdatedA.end())
						{
							auto& refs = qq.second.second;
							if (refs)
								refs--;									
							if (!refs)
								eventsB.insert(qq.first);									
						}
					}	
					for (auto eventB : eventsB)	
						ev->mapIdEvent.erase(eventB);
				}
			}
			// remove all but last eventsUpdated before the first event of any underlying
			if (ok)
			{
				if (this->eventsUpdated.size())
				{
					SizeSet eventsB;
					for (auto eventB : this->eventsUpdated)
						if (eventB < eventLeast)
							eventsB.insert(eventB);	
					eventsB.erase(*this->eventsUpdated.end());
					for (auto eventB : eventsB)
						this->eventsUpdated.erase(eventB);							
				}
			}			
		}
	} 
	catch (const std::exception& e) 
	{
		LOG "Active::update error: " << e.what()  UNLOG
		ok = false;
	}
	
	return ok;
}