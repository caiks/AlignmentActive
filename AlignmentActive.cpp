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

#define ECHO(x) std::cout << #x << std::endl; x
#define EVAL(x) std::cout << #x << ": " << (x) << std::endl
#define EVALL(x) std::cout << #x << ": " << std::endl << (x) << std::endl
#define TRUTH(x) std::cout << #x << ": " << ((x) ? "true" : "false") << std::endl

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

Active::Active() : terminate(false), log(log_default), historyOverflow(false), historyEvent(0), historySize(0), bits(16), var(0), varSlice(0), induceThreshold(100), logging(false), slicesPathLenMax(1)
{
}

#define UNLOG ; log_str.flush(); this->log(log_str.str());}
#define LOG { std::ostringstream log_str; log_str <<

// event ids should be monotonic and updated no more than once
bool Alignment::Active::update(const ActiveUpdateParameters& pp)
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
			if (!(ok && eventsA.size() && (!eventsUpdatedA.size() || *eventsA.rbegin() > *eventsUpdatedA.rbegin()))) 
				break;
			// get first new event id
			auto eventA = *eventsA.begin();
			if (ok && eventsUpdatedA.size())
			{
				for (auto eventB : eventsA)
					if (eventB > *eventsUpdatedA.rbegin())
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
					for (int i = n-1; i >= 0; i--)
						if (rr1[i])
						{
							v = rr1[i];
							break;
						}
					rr[j] = v;
					if (v && this->slicesPath.find(v) == this->slicesPath.end())
					{
						this->slicesPath[v] = hr1;
						if (n > this->slicesPathLenMax)
							this->slicesPathLenMax = n;
					}
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
				auto mark = (ok && this->logging) ? clk::now() : std::chrono::time_point<clk>();
				auto ll = drmul(this->underlyingHistoryRepa,this->underlyingHistorySparse,*this->decomp,this->historyEvent,pp.mapCapacity);	
				ok = ok && ll;
				if (!ok)
				{
					LOG "Active::update\terror: drmul failed to return a list" UNLOG
				}
				// sync active slices
				std::size_t sliceA = 0;
				if (ok)
				{
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
				if (ok && this->logging)
				{
					LOG "Active::update apply\tevent id: " << eventA << "\thistory id: " << this->historyEvent << "\tslice: " << sliceA << "\tslice size: " << this->slicesSetEvent[sliceA].size() << "\ttime " << ((sec)(clk::now() - mark)).count() << "s" UNLOG
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
					eventsB.erase(*this->eventsUpdated.rbegin());
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

// assume that only one inducer so that can modify the var refs without locking active
bool Alignment::Active::induce(const ActiveInduceParameters& pp)
{
	bool ok = true;
	if (this->terminate)
		return true;
	try 
	{
		while (ok)
		{
			std::size_t sliceA = 0;
			std::size_t sliceSizeA = 0;			
			auto hrr = std::make_unique<HistoryRepa>();
			auto haa = std::make_unique<HistorySparseArray>();
			std::unordered_map<std::size_t,HistorySparseArrayPtr> slppa;
			if (ok)
			{			
				std::lock_guard<std::mutex> guard(this->mutex);		
				auto& llr = this->underlyingHistoryRepa;
				auto& lla = this->underlyingHistorySparse;
				// check consistent underlying
				if (ok)
				{
					ok = ok && (llr.size() || lla.size());
					for (auto& hr : llr)
						ok = ok && hr && hr->size == historySize && hr->dimension > 0 && hr->evient;
					for (auto& hr : lla)
						ok = ok && hr && hr->size == historySize && hr->capacity == 1;
					if (!ok)
					{
						LOG "Active::induce\terror: inconsistent underlying" UNLOG
						break;
					}	
				}			
				// get largest slice
				if (ok)
				{
					for (auto sliceB : this->slicesInduce)
					{
						auto sliceSizeB = this->slicesSetEvent[sliceB].size();
						if (sliceSizeB > sliceSizeA)
						{
							sliceA = sliceB;
							sliceSizeA = sliceSizeB;
						}
					}				
				}
				// if there is a slice then process it otherwise return
				if (!(ok && sliceSizeA)) 
					break;
				// copy events from evient active history to varient selection
				if (ok)
				{
					auto& setEventsA = this->slicesSetEvent[sliceA];
					SizeList eventsA(setEventsA.begin(),setEventsA.end());			
					SizeSet qqr;
					if (ok && llr.size())
					{
						auto& qqx = this->induceVarExlusions;					
						for (auto& hr : llr)
						{
							auto n = hr->dimension;
							auto vv = hr->vectorVar;
							for (std::size_t i = 0; i < n; i++)
							{
								auto v = vv[i];
								if (qqx.find(v) == qqx.end())
									qqr.insert(v);
							}
						}
					}
					if (ok && qqr.size())
					{
						hrr->dimension = qqr.size();
						auto nr = hrr->dimension;
						hrr->vectorVar = new std::size_t[nr];
						auto vvr = hrr->vectorVar;
						hrr->shape = new std::size_t[nr];
						auto shr = hrr->shape;
						hrr->size = eventsA.size();
						auto zr = hrr->size;
						hrr->evient = false;
						hrr->arr = new unsigned char[zr*nr];
						auto rrr = hrr->arr;		
						auto ev = eventsA.data();
						{
							std::size_t i = 0;
							for (auto v : qqr)
							{
								vvr[i] = v;
								for (auto& hr : llr)
								{
									auto& mvv = hr->mapVarInt();
									auto it = mvv.find(v);
									if (it != mvv.end())
									{
										auto n = hr->dimension;
										auto rr = hr->arr;
										auto k = it->second;
										shr[i] = hr->shape[k];
										auto izr = i*zr;
										for (std::size_t j = 0; j < zr; j++)
											rrr[izr + j] = rr[ev[j]*n + k];
										break;
									}
								}
								i++;
							}
						}
					}
					SizeSizeUMap qqa;
					if (ok && lla.size())
					{
						auto za = eventsA.size(); 
						auto na = lla.size();
						auto ev = eventsA.data();
						auto& slpp = this->slicesPath;
						haa->size = za;
						haa->capacity = na;
						haa->arr = new std::size_t[za*na];
						auto raa = haa->arr;
						slppa.reserve(za*na);
						for (std::size_t i = 0; i < na; i++)
						{
							auto& hr = lla[i];
							auto rr = hr->arr;
							for (std::size_t j = 0; j < za; j++)
							{
								auto v = rr[ev[j]];
								raa[j*na + i] = v;
								if (v)
								{
									auto it = slpp.find(v);
									if (it != slpp.end())
										slppa.insert_or_assign(it->first, it->second);
								}
							}
						}
					}
				}
			}
			EVAL(sliceA);
			EVAL(sliceSizeA);
			EVAL(*hrr);
			EVAL(*haa);
			EVAL(slppa.size());
			ok = false;
			
				// if (ok && (qqr.size() || qqa.size()))
				// {
					// auto za = eventsA.size(); 
					// auto nmax = (std::size_t)std::sqrt(pp.znnmax / (double)(2*za));
					// nmax = std::max(nmax, pp.bmax);
					// if (nmax < qqr.size() + qqa.size())
					// {
						// if (qqr.size())
						// {
							// // auto ee = prents(*hrpr(vv.size(), vv.data(), *hr));
							
						// }
						// if (qqa.size())
						// {
							
						// }
						// // std::sort(ee->begin(), ee->end());
					// }				
				// }					

		}
	} 
	catch (const std::exception& e) 
	{
		LOG "Active::induce error: " << e.what()  UNLOG
		ok = false;
	}
	
	return ok;
}