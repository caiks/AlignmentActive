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
#include <cstring>

#define ECHO(x) std::cout << #x << std::endl; x
#define EVAL(x) std::cout << #x << ": " << (x) << std::endl
#define EVALL(x) std::cout << #x << ": " << std::endl << (x) << std::endl
#define TRUTH(x) std::cout << #x << ": " << ((x) ? "true" : "false") << std::endl

using namespace Alignment;

const double repaRounding = 1e-6;

typedef std::chrono::duration<double> sec;
typedef std::chrono::high_resolution_clock clk;

void log_default(const std::string& str)
{
	std::cout << str << std::endl;
	return;
};

ActiveSystem::ActiveSystem() : bits(16), block(0)
{
}

std::size_t Alignment::ActiveSystem::next(int bitsA)
{
	std::lock_guard<std::mutex> guard(this->mutex);
	if (bitsA <= this->bits)
	{
		this->block++;
		if (this->block > (std::size_t(-1) >> bits))
			throw std::out_of_range("ActiveSystem::next");
		return this->block << bits;
	}
	auto bitsDiff = bitsA - this->bits;
	auto blockA = ((this->block >> bitsDiff) + 1) << bitsDiff;
	this->block  = blockA + (1 << bitsDiff) - 1;
	if (this->block > (std::size_t(-1) >> bits))
		throw std::out_of_range("ActiveSystem::next");
	return blockA << bits;
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

Active::Active() : terminate(false), log(log_default), historyOverflow(false), historyEvent(0), historySize(0), bits(16), var(0), varSlice(0), induceThreshold(100), logging(false), pathLenMax(1), updateCallback(0)
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
			SizeSet eventsA;
			std::size_t eventA = 0;
			std::size_t historyEventA = 0;
			std::size_t sliceA = 0;
			std::lock_guard<std::mutex> guard(this->mutex);
			if (ok)
			{				
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
						LOG "update\terror: inconsistent underlying" UNLOG
						break;
					}	
				}
				// find set of intersecting event ids
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
				eventA = *eventsA.begin();
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
								LOG "update\terror: lost event repa " << eventA UNLOG
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
								LOG "update\terror: lost event sparse" << eventA UNLOG
							}									
						}
					if (!ok)
						break;					
					ok = ok && this->historyEvent < this->historySize;
					if (!ok) 
					{
						LOG "update\terror: inconsistent historyEvent " << this->historyEvent << " compared to historySize " << this->historySize UNLOG		
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
						LOG "update\terror: inconsistent underlying " UNLOG		
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
							if (n > this->pathLenMax)
								this->pathLenMax = n;
						}
					}
				}
				if (ok)
				{
					ok = ok && this->decomp;
					if (!ok)
					{
						LOG "update\terror: no decomp set" UNLOG
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
						LOG "update\terror: drmul failed to return a list" UNLOG
					}
					// sync active slices
					if (ok)
					{
						if (ll->size())
							sliceA = ll->back();
						if (this->historySparse.size != this->historySize)
						{
							auto z = this->historySize;
							auto& hr = this->historySparse;
							delete[] hr.arr;
							hr.size = z;
							hr.capacity = 1;
							hr.arr = new std::size_t[z];
							memset(hr.arr, 0, z*sizeof(std::size_t));
						}
						std::size_t sliceB = this->historySparse.arr[this->historyEvent];
						if (!sliceA || sliceA != sliceB)
						{
							this->historySparse.arr[this->historyEvent] = sliceA;
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
						LOG "update apply\tevent id: " << eventA << "\thistory id: " << this->historyEvent << "\tslice: " << sliceA << "\tslice size: " << this->slicesSetEvent[sliceA].size() << "\ttime " << ((sec)(clk::now() - mark)).count() << "s" UNLOG
					}
				}
				// increment historyEvent
				if (ok)
				{
					historyEventA = this->historyEvent;
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
			if (ok && eventsA.size() && updateCallback)
			{
				ok = ok && updateCallback(eventsA,eventA,historyEventA,sliceA);
			}
		}
	} 
	catch (const std::exception& e) 
	{
		LOG "update error: " << e.what()  UNLOG
		ok = false;
	}
	
	return ok;
}

bool Alignment::Active::induce(const ActiveInduceParameters& pp)
{
	auto hrhrred = setVarsHistoryRepasHistoryRepaReduced_u;
	auto hrjoin = vectorHistoryRepasJoin_u;	
	auto hrpr = setVarsHistoryRepasRed_u;
	auto prents = histogramRepaRedsListEntropy;
	auto hrshuffle = historyRepasShuffle_us;
	auto hashuffle = historySparseArrayShuffle_us;
	auto frvars = fudRepasSetVar;
	auto frder = fudRepasDerived;
	auto frund = fudRepasUnderlying;
	auto llfr = setVariablesListTransformRepasFudRepa_u;
	auto frmul = historyRepasFudRepasMultiply_up;
	auto frdep = fudRepasSetVarsDepends;
	auto layerer = parametersLayererMaxRollByMExcludedSelfHighestLogIORepa_up;
		
	bool ok = true;
	if (this->terminate)
		return true;
	try 
	{
		while (ok)
		{
			std::size_t varA = 0;
			std::size_t sliceA = 0;
			std::size_t sliceSizeA = 0;			
			std::unique_ptr<HistoryRepa> hrr;
			std::unique_ptr<HistorySparseArray> haa;
			SizeSet qqr;
			std::unordered_map<std::size_t,HistorySparseArrayPtr> slppa;
			std::size_t slppalen = 0;		
			// copy repa and sparse from locked active
			if (ok)
			{			
				auto mark = (ok && this->logging) ? clk::now() : std::chrono::time_point<clk>();
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
						LOG "induce\terror: inconsistent underlying" UNLOG
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
							auto it = this->sliceFailsSize.find(sliceB);
							if (it == this->sliceFailsSize.end() || it->second < sliceSizeB)
							{
								sliceA = sliceB;
								sliceSizeA = sliceSizeB;							
							}
						}
					}				
				}
				// if there is a slice then process it otherwise return
				if (!(ok && sliceSizeA)) 
					break;
				// copy events from evient active history to varient selection
				if (ok)
				{
					varA = this->var;
					auto& setEventsA = this->slicesSetEvent[sliceA];
					SizeList eventsA(setEventsA.begin(),setEventsA.end());			
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
						hrr = std::make_unique<HistoryRepa>();
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
					if (ok && lla.size())
					{
						auto za = eventsA.size(); 
						auto na = lla.size();
						auto ev = eventsA.data();
						auto& slpp = this->slicesPath;
						haa = std::make_unique<HistorySparseArray>();
						haa->size = za;
						haa->capacity = na;
						haa->arr = new std::size_t[za*na];
						auto raa = haa->arr;
						slppa.reserve(za*na);
						slppalen = this->pathLenMax;
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
				if (ok && this->logging)
				{
					LOG "induce copy\tslice: " << sliceA << "\tslice size: " << sliceSizeA << "\trepa dimension: " << (hrr ? hrr->dimension : 0) << "\tsparse capacity: " << (haa ? haa->capacity : 0) << "\tsparse paths: " << slppa.size() << "\ttime " << ((sec)(clk::now() - mark)).count() << "s" UNLOG
				}	
			}
			// check consistent copy
			if (ok)
			{
				ok = ok && (hrr || haa);
				ok = ok && (!hrr || (hrr->dimension > 0 && hrr->size == sliceSizeA && hrr->arr));
				ok = ok && (!haa || (haa->capacity > 0 && haa->size == sliceSizeA && haa->arr));
				if (!ok)
				{
					LOG "induce\terror: inconsistent copy" UNLOG
					break;
				}	
			}
			bool fail = false;
			std::unique_ptr<HistoryRepa> hr;
			std::unique_ptr<FudRepa> fr;
			SizeList kk;
			// induce while unlocked
			if (ok)
			{
				auto mark = (ok && this->logging) ? clk::now() : std::chrono::time_point<clk>();
				SizeSizeUMap qqa;
				std::unordered_map<std::size_t, SizeSet> mma;
				// prepare for the sparse entropy calculations
				if (ok && haa && haa->size && haa->capacity)
				{
					auto za = haa->size; 
					auto na = haa->capacity; 
					auto raa = haa->arr;
					qqa.reserve((slppalen > 0 ? slppalen : 1) * za * na * 3);
					mma.reserve((slppalen > 0 ? slppalen : 1) * za * na * 3);
					for (std::size_t k = 0; k < na; k++)
					{
						for (std::size_t j = 0; j < za; j++)
						{
							auto v = raa[j*na + k];
							auto it = slppa.find(v);
							if (it != slppa.end())
							{
								auto& hr1 = it->second;
								auto n1 = hr1->capacity;
								auto rr1 = hr1->arr;
								for (std::size_t i = 0; i < n1; i++)
								{
									auto w = rr1[i];
									if (w)
									{
										qqa[w]++;
										for (std::size_t m = i+1; m < n1; m++)
										{
											auto x = rr1[m];
											if (x)
												mma[w].insert(x);
										}
									}									
								}
							}
							else
							{
								qqa[v]++;
							}
						}
					}
				}
				// get top nmax vars by entropy
				// remove any sparse parents with same entropy as children
				if (ok && (qqr.size() || qqa.size()))
				{
					// EVAL(qqr.size());
					// EVAL(qqr);
					// EVAL(qqa.size());
					// EVAL(qqa);
					auto nmax = (std::size_t)std::sqrt(pp.znnmax / (double)(2*sliceSizeA));
					// EVAL(nmax);
					nmax = std::max(nmax, pp.bmax);
					// EVAL(nmax);
					DoubleSizePairList ee;
					ee.reserve(qqr.size() + qqa.size());
					if (qqr.size())
					{
						SizeList vv(qqr.begin(),qqr.end());
						auto eer = prents(*hrpr(vv.size(), vv.data(), *hrr));
						// EVAL(*hrr);
						// EVAL(*hrpr(vv.size(), vv.data(), *hrr));
						// EVAL(*eer);
						for (auto p : *eer)
							if (p.first > repaRounding)
								ee.push_back(DoubleSizePair(-p.first,p.second));
					}
					// EVAL(ee);
					if (qqa.size())
					{
						std::map<std::size_t, SizeSet> eem;
						for (auto p : qqa)		
						{
							if (p.second > 0 && p.second < sliceSizeA)
							{
								auto v = p.first; 
								auto& xx = eem[p.second];
								if (xx.size())
								{
									auto iv = mma.find(v);
									if (iv != mma.end())
									{
										bool found = false;
										for (auto w : xx)
										{
											if (iv->second.find(w) != iv->second.end())
											{
												found = true;
												break;
											}
										}
										if (!found)								
											xx.insert(v);
									}
									else
										xx.insert(v);
									SizeSet yy(xx);
									for (auto w : yy)
									{
										auto iw = mma.find(w);
										if (iw != mma.end() && iw->second.find(v) != iw->second.end())
											xx.erase(w);
									}
								}
								else
									xx.insert(v);
							}
						}	
						// EVAL(eem);
						double f = 1.0/(double)sliceSizeA;
						for (auto& p : eem)	
							for (auto& q : p.second)	
							{
								double a = (double)p.first * f;
								double e = -(a * std::log(a) + (1.0-a) * std::log(1.0-a));
								if (e > repaRounding)
									ee.push_back(DoubleSizePair(-e,q));
							}
					}
					std::sort(ee.begin(), ee.end());
					// EVAL(ee.size());
					// EVAL(ee);
					SizeUSet qq;
					qq.reserve(ee.size());
					for (std::size_t i = 0; i < nmax && i < ee.size(); i++)
						qq.insert(ee[i].second);
					SizeSet qqr1(qqr);
					for (auto v : qqr1)
						if (qq.find(v) == qq.end())
							qqr.erase(v);
					SizeUSet qqa1;
					qqa1.reserve(qqa.size());
					for (auto p : qqa)
						qqa1.insert(p.first);
					for (auto v : qqa1)
						if (qq.find(v) == qq.end())
							qqa.erase(v);
					for (auto v : qqr)
						qqa.erase(v);
					fail = ok && !qqr.size() && !qqa.size();
					if (ok && this->logging)
					{
						if (!fail)
						{
							LOG "induce model\trepa dimension: " << qqr.size() << "\tsparse dimension: " << qqa.size() UNLOG							
						}
						else
						{
							LOG "induce model\tno entropy"  UNLOG
						}
					}	
					// EVAL(qqr.size());
					// EVAL(qqr);
					// EVAL(qqa.size());
					// EVAL(qqa);
				}	
				if (ok && !fail)
				{
					std::unique_ptr<HistoryRepa> hrs;
					std::ranlux48_base gen((unsigned int)pp.seed);
					if (ok && qqr.size())
					{
						// EVAL(*hrr);
						if (qqr.size() < hrr->dimension)
						{
							SizeList vv(qqr.begin(),qqr.end());
							hr = hrhrred(vv.size(), vv.data(), *hrr);
						}
						else
							hr = std::move(hrr);
						// EVAL(*hr);
						hrs = hrshuffle(*hr,gen);
						// EVAL(*hrs);
					}
					if (ok && qqa.size())
					{
						auto has = hashuffle(*haa,gen);				
						auto hra = std::make_unique<HistoryRepa>();
						auto hras = std::make_unique<HistoryRepa>();	
						{					
							auto n = qqa.size();
							hra->dimension = n;
							hras->dimension = n;
							hra->vectorVar = new std::size_t[n];
							hras->vectorVar = new std::size_t[n];
							hra->shape = new std::size_t[n];
							hras->shape = new std::size_t[n];
							auto za = sliceSizeA;
							hra->size = za;
							hras->size = za;
							hra->evient = false;
							hras->evient = false;
							hra->arr = new unsigned char[za*n];
							hras->arr = new unsigned char[za*n];
							auto rra = hra->arr;
							auto rras = hras->arr;
							std::memset(rra, 0, za*n);
							std::memset(rras, 0, za*n);
							auto na = haa->capacity;
							auto raa = haa->arr;						
							auto ras = has->arr;
							{
								std::size_t i = 0;
								for (auto p : qqa)
								{
									auto v = p.first;
									hra->vectorVar[i] = v;	
									hras->vectorVar[i] = v;	
									hra->shape[i] = 2;
									hras->shape[i] = 2;
									i++;						
								}
							}
							auto& mvv = hra->mapVarInt();						
							for (std::size_t i = 0; i < na; i++)
							{
								for (std::size_t j = 0; j < za; j++)
								{
									{
										auto v = raa[j*na + i];
										if (v)
										{
											auto iv = slppa.find(v);
											if (iv != slppa.end())
											{
												auto& hr1 = iv->second;
												auto nr1 = hr1->capacity;
												auto rr1 = hr1->arr;
												for (std::size_t k = 0; k < nr1; k++)
												{
													auto w = rr1[k];
													if (w)
													{
														auto iw = mvv.find(w);
														if (iw != mvv.end())
															rra[iw->second * za + j] = 1;
													}
												}
											}
											else
											{
												auto iw = mvv.find(v);
												if (iw != mvv.end())
													rra[iw->second * za + j] = 1;
											}
										}								
									}
									{
										auto v = ras[j*na + i];
										if (v)
										{
											auto iv = slppa.find(v);
											if (iv != slppa.end())
											{
												auto& hr1 = iv->second;
												auto nr1 = hr1->capacity;
												auto rr1 = hr1->arr;
												for (std::size_t k = 0; k < nr1; k++)
												{
													auto w = rr1[k];
													if (w)
													{
														auto iw = mvv.find(w);
														if (iw != mvv.end())
															rras[iw->second * za + j] = 1;
													}
												}
											}
											else
											{
												auto iw = mvv.find(v);
												if (iw != mvv.end())
													rras[iw->second * za + j] = 1;
											}
										}								
									}
								}						
							}													
						}					
						// EVAL(*hra);					
						// EVAL(*hras);					
						if (hr)
						{
							hr = hrjoin(HistoryRepaPtrList{std::move(hr),std::move(hra)});
							hrs = hrjoin(HistoryRepaPtrList{std::move(hrs),std::move(hras)});
						}
						else
						{
							hr = std::move(hra);
							hrs = std::move(hras);
						}	
						// EVAL(*hr);		
						// EVAL(*hrs);					
					}
					// check consistent reduction
					if (ok)
					{
						ok = ok && hr && hrs 
							&& hr->dimension == (qqr.size() + qqa.size()) 
							&& hr->dimension == hrs->dimension 
							&& hr->size == sliceSizeA
							&& hr->size == hrs->size;
						if (!ok)
						{
							LOG "induce\terror: inconsistent reduction" UNLOG
							break;
						}	
					}
					// check active system
					if (ok)
					{
						ok = ok && this->system;
						if (!ok)
						{
							LOG "induce\terror: no system" UNLOG
							break;
						}	
					}
					if (ok && this->logging)
					{
						LOG "induce model\tdimension: " << hr->dimension << "\tsize: " << hr->size UNLOG
					}						
					// layerer
					if (ok)
					{
						try
						{
							SizeList vv;
							{
								auto n = hr->dimension;
								auto vv1 = hr->vectorVar;
								vv.reserve(n);
								for (std::size_t i = 0; i < n; i++)
									vv.push_back(vv1[i]);
							}
							auto t = layerer(pp.wmax, pp.lmax, pp.xmax, pp.omax, pp.bmax, pp.mmax, pp.umax, pp.pmax, pp.tint, vv, *hr, *hrs, this->log, this->logging, varA);
							fr = std::move(std::get<0>(t));
							auto mm = std::move(std::get<1>(t));
							fail = !fr || (!mm || !mm->size());
							TRUTH(fail);
							if (ok && !fail)
							{
								kk = mm->back().second;
								SizeUSet kk1(kk.begin(), kk.end());
								SizeUSet vv1(vv.begin(), vv.end());
								fr = llfr(vv1, *frdep(*fr, kk1));
								EVAL(mm->size());			
								EVAL(mm->back());
								EVAL(*mm);
								auto& a = mm->back().first;
								auto m = kk.size();
								auto z = hr->size;
								EVAL(m);
								EVAL(a);			
								EVAL(z);	
								EVAL(100.0*(exp(a/z/(m-1))-1.0));
								EVAL(fudRepasSize(*fr));
								EVAL(frvars(*fr)->size());
								EVAL(frder(*fr)->size());
								EVAL(frund(*fr)->size());
								EVAL(sorted(*frund(*fr)));
							}
						}
						catch (const std::out_of_range& e)
						{
							ok = false;
							LOG "induce\tout of range exception: " << e.what() UNLOG
							break;
						}
					}
				}
				if (ok && this->logging)
				{
					LOG "induce model\tslice: " << sliceA << "\tslice size: " << sliceSizeA << "\ttime " << ((sec)(clk::now() - mark)).count() << "s" UNLOG
				}				
			}
			// add new fud to locked active
			if (ok && !fail)	
			{
				auto mark = (ok && this->logging) ? clk::now() : std::chrono::time_point<clk>();
				std::lock_guard<std::mutex> guard(this->mutex);		
				
				// remap kk and fr with block ids
				if (ok)
				{
					auto frSize = fudRepasSize(*fr);
					if (frSize > (1 << this->bits))
					{
						ok = false;
						LOG "induce\terror: block too small" << "\tfud size: " << frSize << "\tblock: " << (1 << this->bits) UNLOG
						break;								
					}
					if (((this->var + frSize) >> this->bits) > (this->var >> this->bits))
						this->var = this->system->next(this->bits);
					SizeSizeUMap nn;
					nn.reserve(frSize);
					for (auto& ll : fr->layers)
						for (auto& tr : ll)
						{
							nn[tr->derived] = this->var;					
							this->var++;
						}
					fr->reframe_u(nn);
					for (std::size_t i = 0; i < kk.size(); i++)	
						kk[i] = nn[kk[i]];
				}
					
				EVAL(kk);
				EVAL(sorted(*frvars(*fr)));			
				
				// auto dr = std::make_unique<ApplicationRepa>();
				// {
					// dr->substrate = vv;
					// dr->fud = std::make_shared<FudRepa>();
					// dr->slices = std::make_shared<SizeTree>();	
					// auto vd = std::make_shared<Variable>(d);
					// auto vl = std::make_shared<Variable>("s");
					// auto vf = std::make_shared<Variable>((int)f);
					// auto vdf = std::make_shared<Variable>(vd, vf);
					// auto vfl = std::make_shared<Variable>(vdf, vl);
					// SizeUSet kk1(kk.begin(), kk.end());
					// SizeUSet vv1(vv.begin(), vv.end());
					// auto gr = llfr(vv1, *frdep(*fr, kk1));
					// auto ar = hrred(1.0, m, kk.data(), *frmul(tint, *hr, *gr));
					// SizeList sl;
					// TransformRepaPtrList ll;
					// std::size_t sz = 1;
					// auto skk = ar->shape;
					// auto rr0 = ar->arr;
					// for (std::size_t i = 0; i < m; i++)
						// sz *= skk[i];
					// sl.reserve(sz);
					// ll.reserve(sz);
					// bool remainder = false;
					// std::size_t b = 1;
					// auto& llu = ur->listVarSizePair;
					// for (std::size_t i = 0; i < sz; i++)
					// {
						// if (rr0[i] <= 0.0)
						// {
							// remainder = true;
							// continue;
						// }
						// auto tr = std::make_shared<TransformRepa>();
						// tr->dimension = m;
						// tr->vectorVar = new std::size_t[m];
						// auto ww = tr->vectorVar;
						// tr->shape = new std::size_t[m];
						// auto sh = tr->shape;
						// for (std::size_t j = 0; j < m; j++)
						// {
							// ww[j] = kk[j];
							// sh[j] = skk[j];
						// }
						// tr->arr = new unsigned char[sz];
						// auto rr = tr->arr;
						// for (std::size_t j = 0; j < sz; j++)
							// rr[j] = 0;
						// rr[i] = 1;
						// tr->valency = 2;
						// auto vb = std::make_shared<Variable>((int)b++);
						// auto vflb = std::make_shared<Variable>(vfl, vb);
						// llu.push_back(VarSizePair(vflb, 2));
						// auto w = llu.size() - 1;
						// tr->derived = w;
						// sl.push_back(w);
						// ll.push_back(tr);
					// }
					// if (remainder)
					// {
						// auto tr = std::make_shared<TransformRepa>();
						// tr->dimension = m;
						// tr->vectorVar = new std::size_t[m];
						// auto ww = tr->vectorVar;
						// tr->shape = new std::size_t[m];
						// auto sh = tr->shape;
						// for (std::size_t j = 0; j < m; j++)
						// {
							// ww[j] = kk[j];
							// sh[j] = skk[j];
						// }
						// tr->arr = new unsigned char[sz];
						// auto rr = tr->arr;
						// for (std::size_t j = 0; j < sz; j++)
							// rr[j] = rr0[j] <= 0.0 ? 1 : 0;
						// tr->valency = 2;
						// auto vb = std::make_shared<Variable>((int)b++);
						// auto vflb = std::make_shared<Variable>(vfl, vb);
						// llu.push_back(VarSizePair(vflb, 2));
						// auto w = llu.size() - 1;
						// tr->derived = w;
						// sl.push_back(w);
						// ll.push_back(tr);
					// }
					// dr->fud->layers.insert(dr->fud->layers.end(), gr->layers.begin(), gr->layers.end());
					// dr->fud->layers.push_back(ll);
					// dr->slices->_list.reserve(sz);
					// for (auto& s : sl)
						// dr->slices->_list.push_back(SizeSizeTreePair(s, std::make_shared<SizeTree>()));			
				// }
							
				// update this decomp mapVarParent and mapVarInt
				// update historySparse and slicesSetEvent
				// tidy slicesInduce and sliceFailsSize
				// tidy new events
				
				if (ok && this->logging)
				{
					LOG "induce copy\tslice: " << sliceA << "\tslice size: " << sliceSizeA << "\trepa dimension: " << (hrr ? hrr->dimension : 0) << "\tsparse capacity: " << (haa ? haa->capacity : 0) << "\tsparse paths: " << slppa.size() << "\ttime " << ((sec)(clk::now() - mark)).count() << "s" UNLOG
				}	
			}

			EVAL(sliceA);
			EVAL(sliceSizeA);
			// EVAL(*hrr);
			// EVAL(*haa);
			EVAL(slppa.size());
			ok = false;
			

		}
	} 
	catch (const std::exception& e) 
	{
		LOG "induce error: " << e.what()  UNLOG
		ok = false;
	}
	
	return ok;
}