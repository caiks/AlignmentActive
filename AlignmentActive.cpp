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

Active::Active() : terminate(false), log(log_default), historyOverflow(false), historyEvent(0), historySize(0), bits(16), var(0), varSlice(0), induceThreshold(100), logging(false),  updateCallback(0)
{
}

#define UNLOG ; log_str.flush(); this->log(log_str.str());}
#define LOG { std::ostringstream log_str; log_str <<

// event ids should be monotonic and updated no more than once
bool Alignment::Active::update(ActiveUpdateParameters pp)
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
				SizeSet eventsUpdatedA(this->underlyingEventUpdateds);
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
						if (v && this->underlyingSlicesParent.find(v) == this->underlyingSlicesParent.end())
							for (int i = n-1; i > 0; i--)
								if (rr1[i])
									this->underlyingSlicesParent[rr1[i]] = rr1[i-1];
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
							// memset(hr.arr, 0, z*sizeof(std::size_t));
						}
						std::size_t sliceB = this->historySparse.arr[this->historyEvent];
						if (!sliceA || sliceA != sliceB)
						{
							this->historySparse.arr[this->historyEvent] = sliceA;
							auto& setA = this->historySlicesSetEvent[sliceA];
							setA.insert(this->historyEvent);
							if (this->induceThreshold && setA.size() == this->induceThreshold)
								this->induceSlices.insert(sliceA);
							if (sliceA)
							{
								auto& setB = this->historySlicesSetEvent[sliceB];
								setB.erase(this->historyEvent);
								if (this->induceThreshold && setB.size() == this->induceThreshold-1)
									this->induceSlices.erase(sliceB);
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
						LOG "update apply\tevent id: " << eventA << "\thistory id: " << this->historyEvent << "\tslice: " << sliceA << "\tslice size: " << this->historySlicesSetEvent[sliceA].size() << "\ttime " << ((sec)(clk::now() - mark)).count() << "s" UNLOG
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
							this->underlyingEventUpdateds.insert(eventB);
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
				// remove all but last underlyingEventUpdateds before the first event of any underlying
				if (ok)
				{
					if (this->underlyingEventUpdateds.size())
					{
						SizeSet eventsB;
						for (auto eventB : this->underlyingEventUpdateds)
							if (eventB < eventLeast)
								eventsB.insert(eventB);	
						eventsB.erase(*this->underlyingEventUpdateds.rbegin());
						for (auto eventB : eventsB)
							this->underlyingEventUpdateds.erase(eventB);							
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

bool Alignment::Active::induce(ActiveInduceParameters pp, ActiveUpdateParameters ppu)
{
	auto hrred = setVarsHistoryRepasReduce_u;
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
	auto drmul = historyRepaPtrListsHistorySparseArrayPtrListsDecompFudSlicedRepasEventsPathSlice_u;
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
			SizeList eventsA;			
			std::unique_ptr<HistoryRepa> hrr;
			std::unique_ptr<HistorySparseArray> haa;
			SizeSet qqr;
			SizeSizeUMap slppa;
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
					for (auto sliceB : this->induceSlices)
					{
						auto sliceSizeB = this->historySlicesSetEvent[sliceB].size();
						if (sliceSizeB > sliceSizeA)
						{
							auto it = this->induceSliceFailsSize.find(sliceB);
							if (it == this->induceSliceFailsSize.end() || it->second < sliceSizeB)
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
					auto& setEventsA = this->historySlicesSetEvent[sliceA];
					eventsA.insert(eventsA.end(),setEventsA.begin(),setEventsA.end());
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
						auto& slpp = this->underlyingSlicesParent;
						haa = std::make_unique<HistorySparseArray>();
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
									while (it != slpp.end())
									{
										slppa.insert_or_assign(it->first, it->second);
										it = slpp.find(it->second);
									}
								}
							}
						}
						// EVAL(sorted(slppa));
					}
				}
				if (ok && this->logging)
				{
					LOG "induce copy\tslice: " << sliceA << "\tslice size: " << sliceSizeA << "\trepa dimension: " << (hrr ? hrr->dimension : 0) << "\tsparse capacity: " << (haa ? haa->capacity : 0) << "\tsparse paths: " << slppa.size() << "\tvariable: " << varA << "\ttime " << ((sec)(clk::now() - mark)).count() << "s" UNLOG
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
			std::size_t frSize = 0;
			SizeList kk;
			// induce model while unlocked
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
					qqa.reserve(slppa.size());
					mma.reserve(slppa.size());
					for (std::size_t k = 0; k < na; k++)
					{
						for (std::size_t j = 0; j < za; j++)
						{
							auto v = raa[j*na + k];
							SizeList ll {v};
							qqa[v]++;
							auto it = slppa.find(v);
							while (it != slppa.end())
							{
								ll.push_back(it->second);
								qqa[it->second]++;
								it = slppa.find(it->second);
							}								
							for (int i = ll.size() - 1; i > 0; i--)
								for (int m = i-1; m >= 0; m--)
									mma[ll[i]].insert(ll[m]);
						}
					}
					// EVAL(sorted(qqa));
					// EVAL(sorted(mma));
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
											{
												auto iw = mvv.find(v);
												if (iw != mvv.end())
													rra[iw->second * za + j] = 1;
											}
											auto iv = slppa.find(v);
											while (iv != slppa.end())
											{
												auto iw = mvv.find(iv->second);
												if (iw != mvv.end())
													rra[iw->second * za + j] = 1;
												iv = slppa.find(iv->second);
											}
										}								
									}
									{
										auto v = ras[j*na + i];
										if (v)
										{
											{
												auto iw = mvv.find(v);
												if (iw != mvv.end())
													rras[iw->second * za + j] = 1;
											}
											auto iv = slppa.find(v);
											while (iv != slppa.end())
											{
												auto iw = mvv.find(iv->second);
												if (iw != mvv.end())
													rras[iw->second * za + j] = 1;
												iv = slppa.find(iv->second);
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
					if (ok && this->logging)
					{
						LOG "induce model\tdimension: " << hr->dimension << "\tsize: " << hr->size UNLOG
					}						
					// layerer
					if (ok)
					{
						double algn = 0.0;
						double diagonal = 0.0;
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
							auto t = layerer(pp.wmax, pp.lmax, pp.xmax, pp.omax, pp.bmax, pp.mmax, pp.umax, pp.pmax, pp.tint, vv, *hr, *hrs, this->log, this->logging && pp.logging, varA);
							fr = std::move(std::get<0>(t));
							auto mm = std::move(std::get<1>(t));
							fail = !fr || (!mm || !mm->size());
							if (ok && !fail)
							{
								kk = mm->back().second;
								SizeUSet kk1(kk.begin(), kk.end());
								SizeUSet vv1(vv.begin(), vv.end());
								fr = llfr(vv1, *frdep(*fr, kk1));
								frSize = fudRepasSize(*fr);
								algn = mm->back().first;
								auto m = kk.size();
								auto z = hr->size;
								diagonal = 100.0*(exp(algn/z/(m-1))-1.0);
								// EVAL(frvars(*fr)->size());
								// EVAL(frder(*fr)->size());
								// EVAL(frund(*fr)->size());
								// EVAL(sorted(*frund(*fr)));
							}
						}
						catch (const std::out_of_range& e)
						{
							ok = false;
							LOG "induce\tout of range exception: " << e.what() UNLOG
							break;
						}
						if (ok && this->logging)
						{
							if (!fail)
							{
								LOG "induce model\tder vars algn density: " << algn << "\timpl bi-valency percent: " << diagonal << "\tder vars cardinality: " << kk.size() << "\tfud cardinality: " << frSize UNLOG							
							}
							else
							{
								LOG "induce model\tno alignment"  UNLOG
							}
						}	
					}
				}
				if (ok && this->logging)
				{
					LOG "induce model\ttime " << ((sec)(clk::now() - mark)).count() << "s" UNLOG
				}				
			}
			// add new fud to locked active and update
			if (ok && !fail)	
			{
				auto mark = (ok && this->logging) ? clk::now() : std::chrono::time_point<clk>();
				std::lock_guard<std::mutex> guard(this->mutex);		
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
				// remap kk and fr with block ids
				if (ok)
				{
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
					// EVAL(kk);
					// EVAL(sorted(*frvars(*fr)));			
				}				
				std::size_t v = 0;
				SizeList sl;
				// create the slices
				if (ok)
				{		
					if (!this->decomp)
						this->decomp = std::make_shared<DecompFudSlicedRepa>();
					if (sliceA)
					{
						auto& mm = this->decomp->mapVarParent();
						auto it = mm.find(sliceA);
						if (it != mm.end())
							v = it->second;
						else
						{
							ok = false;
							LOG "induce update\terror: cannot find parent slice of slice " << sliceA UNLOG
							break;							
						}
					}
					auto m = kk.size();
					auto ar = hrred(1.0, m, kk.data(), *frmul(pp.tint, *hr, *fr));
					std::size_t sz = 1;
					auto skk = ar->shape;
					auto rr0 = ar->arr;
					for (std::size_t i = 0; i < m; i++)
						sz *= skk[i];
					sl.reserve(sz);
					fr->layers.push_back(TransformRepaPtrList());
					auto& ll = fr->layers.back();
					ll.reserve(sz);					
					if (sz > (1 << this->bits))
					{
						ok = false;
						LOG "induce\terror: block too small" << "\tslice size: " << sz << "\tblock: " << (1 << this->bits) UNLOG
						break;								
					}
					if (((this->varSlice + sz) >> this->bits) > (this->varSlice >> this->bits))
						this->varSlice = this->system->next(this->bits);					
					bool remainder = false;
					for (std::size_t i = 0; i < sz; i++)
					{
						if (rr0[i] <= 0.0)
						{
							remainder = true;
							continue;
						}
						auto tr = std::make_shared<TransformRepa>();
						if (v)
						{
							tr->dimension = m + 1;
							tr->vectorVar = new std::size_t[m + 1];
							auto ww = tr->vectorVar;
							tr->shape = new std::size_t[m + 1];
							auto sh = tr->shape;
							ww[0] = v;
							sh[0] = 2;
							for (std::size_t j = 0; j < m; j++)
							{
								ww[j + 1] = kk[j];
								sh[j + 1] = skk[j];
							}
							tr->arr = new unsigned char[2 * sz];
							auto rr = tr->arr;
							for (std::size_t j = 0; j < 2 * sz; j++)
								rr[j] = 0;
							rr[sz + i] = 1;
						}
						else
						{
							tr->dimension = m;
							tr->vectorVar = new std::size_t[m];
							auto ww = tr->vectorVar;
							tr->shape = new std::size_t[m];
							auto sh = tr->shape;
							for (std::size_t j = 0; j < m; j++)
							{
								ww[j] = kk[j];
								sh[j] = skk[j];
							}
							tr->arr = new unsigned char[sz];
							auto rr = tr->arr;
							for (std::size_t j = 0; j < sz; j++)
								rr[j] = 0;
							rr[i] = 1;
						}
						tr->valency = 2;						
						auto w = this->varSlice;
						this->varSlice++;
						tr->derived = w;
						sl.push_back(w);
						ll.push_back(tr);
					}
					if (remainder)
					{
						auto tr = std::make_shared<TransformRepa>();
						if (v)
						{
							tr->dimension = m + 1;
							tr->vectorVar = new std::size_t[m + 1];
							auto ww = tr->vectorVar;
							tr->shape = new std::size_t[m + 1];
							auto sh = tr->shape;
							ww[0] = v;
							sh[0] = 2;
							for (std::size_t j = 0; j < m; j++)
							{
								ww[j + 1] = kk[j];
								sh[j + 1] = skk[j];
							}
							tr->arr = new unsigned char[2 * sz];
							auto rr = tr->arr;
							for (std::size_t j = 0; j < 2 * sz; j++)
								rr[j] = j >= sz && rr0[j - sz] <= 0.0 ? 1 : 0;
						}
						else
						{
							tr->dimension = m;
							tr->vectorVar = new std::size_t[m];
							auto ww = tr->vectorVar;
							tr->shape = new std::size_t[m];
							auto sh = tr->shape;
							for (std::size_t j = 0; j < m; j++)
							{
								ww[j] = kk[j];
								sh[j] = skk[j];
							}
							tr->arr = new unsigned char[sz];
							auto rr = tr->arr;
							for (std::size_t j = 0; j < sz; j++)
								rr[j] = rr0[j] <= 0.0 ? 1 : 0;
						}
						tr->valency = 2;
						auto w = this->varSlice;
						this->varSlice++;
						tr->derived = w;
						sl.push_back(w);
						ll.push_back(tr);
					}
				}
				// update this decomp mapVarParent and mapVarInt
				if (ok)
				{
					auto& dr = *this->decomp;
					dr.fuds.push_back(FudSlicedStruct());
					auto& fs = dr.fuds.back();
					fs.parent = sliceA;
					fs.children = sl;
					fs.fud.reserve(frSize + sl.size());
					for (auto& ii : fr->layers)
						for (auto& tr : ii)
							fs.fud.push_back(tr);
					auto& vi = dr.mapVarInt();
					vi[sliceA] = dr.fuds.size() - 1;
					auto& cv = dr.mapVarParent();
					for (auto s : sl)
						cv[s] = sliceA;
				}
				// update historySparse and historySlicesSetEvent
				if (ok)
				{
					if (v)
					{
						auto z = hr->size;
						auto hrr = std::make_unique<HistoryRepa>();
						hrr->dimension = 1;
						hrr->vectorVar = new std::size_t[1];
						hrr->vectorVar[0] = v;
						hrr->shape = new std::size_t[1];
						hrr->shape[0] = 2;
						hrr->size = z;
						hrr->evient = hr->evient;
						hrr->arr = new unsigned char[z];
						auto rrr = hrr->arr;		
						for (std::size_t j = 0; j < z; j++)
							rrr[j] = 1;
						hr = hrjoin(HistoryRepaPtrList{std::move(hr),std::move(hrr)});
					}
					hr = hrhrred(sl.size(), sl.data(), *frmul(pp.tint, *hr, *fr));
					auto n = hr->dimension;
					auto vv = hr->vectorVar;
					auto z = hr->size;
					auto rr = hr->arr;	
					auto ev = eventsA.data();
					SizeSet slices;
					for (std::size_t j = 0; j < z; j++)
						for (std::size_t i = 0; i < n; i++)
						{
							auto u = rr[i*z + j];
							if (u)
							{
								auto sliceB = vv[i];
								auto eventA = ev[j];
								this->historySparse.arr[eventA] = sliceB;
								this->historySlicesSetEvent[sliceA].erase(eventA);
								this->historySlicesSetEvent[sliceB].insert(eventA);
								slices.insert(sliceB);
								break;
							}
						}
					for (auto sliceB : slices)
						if (this->historySlicesSetEvent[sliceB].size() >= induceThreshold)
							this->induceSlices.insert(sliceB);
					this->induceSlices.erase(sliceA);
					this->induceSliceFailsSize.erase(sliceA);
				}
				// tidy new events
				if (ok)
				{
					auto eventsB = this->historySlicesSetEvent[sliceA];
					if (eventsB.size())
					{
						SizeSet slices;
						for (auto eventB : eventsB)
						{
							auto ll = drmul(this->underlyingHistoryRepa,this->underlyingHistorySparse,*this->decomp,eventB,ppu.mapCapacity);	
							ok = ok && ll;
							if (!ok)
							{
								LOG "induce update\terror: drmul failed to return a list" UNLOG
								break;
							}	
							std::size_t	sliceB = 0;						
							if (ll->size())
								sliceB = ll->back();
							ok = ok && this->historySparse.arr;
							if (!ok)
							{
								LOG "induce update\terror: historySparse not initialised" UNLOG
								break;
							}
							this->historySparse.arr[eventB] = sliceB;
							this->historySlicesSetEvent[sliceB].insert(eventB);	
							slices.insert(sliceB);							
						}
						if (!ok)
							break;
						for (auto sliceB : slices)
							if (this->historySlicesSetEvent[sliceB].size() >= induceThreshold)
								this->induceSlices.insert(sliceB);
					}
					this->historySlicesSetEvent.erase(sliceA);
				}
				if (ok && this->logging)
				{
					LOG "induce update\tslice: " << sliceA << "\tparent slice: " << v << "\tchildren cardinality: " << sl.size() << "\tchildren slices: " << sl<< "\tmodel cardinality: " << this->decomp->fuds.size() << "\ttime " << ((sec)(clk::now() - mark)).count() << "s" UNLOG
				}	
			}
			if (ok && fail)
			{
				auto mark = (ok && this->logging) ? clk::now() : std::chrono::time_point<clk>();
				std::lock_guard<std::mutex> guard(this->mutex);		
				this->induceSliceFailsSize.insert_or_assign(sliceA, sliceSizeA);
				if (ok && this->logging)
				{
					LOG "induce update fail\tslice: " << sliceA << "\tslice size: " << sliceSizeA  << "\ttime " << ((sec)(clk::now() - mark)).count() << "s" UNLOG
				}					
			}
			// EVAL(sliceA);
			// EVAL(sliceSizeA);
			// EVAL(*hrr);
			// EVAL(*haa);
			// EVAL(slppa.size());
			// if (this->decomp) { EVAL(*this->decomp);}
			// ok = false;
		}
	} 
	catch (const std::exception& e) 
	{
		LOG "induce error: " << e.what()  UNLOG
		ok = false;
	}
	
	return ok;
}