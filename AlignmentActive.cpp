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
#define EVALH(x) std::cout << #x << ": " << std::hex << (x) << std::dec << std::endl
#define EVALL(x) std::cout << #x << ": " << std::endl << (x) << std::endl
#define TRUTH(x) std::cout << #x << ": " << ((x) ? "true" : "false") << std::endl

using namespace Alignment;

const double repaRounding = 1e-6;

typedef std::chrono::duration<double> sec;
typedef std::chrono::system_clock clk;

void log_default(Active& active, const std::string& str)
{
	std::cout << str << std::endl;
	return;
};

void layerer_log_default(const std::string& str)
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
		if (this->block > (std::size_t(-1) >> this->bits))
			throw std::out_of_range("ActiveSystem::next");
		return this->block << this->bits;
	}
	std::size_t blockA = this->block + 1;
	this->block  = blockA + ((std::size_t)1 << (bitsA - this->bits)) - 1;
	if (this->block > (std::size_t(-1) >> this->bits))
		throw std::out_of_range("ActiveSystem::next");
	return blockA << this->bits;
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

Active::Active(std::string nameA) : name(nameA), terminate(false), log(log_default), layerer_log(layerer_log_default), historyOverflow(false), historyEvent(0), historySize(0), continousIs(false), bits(16), var(0), varSlice(0), induceThreshold(100), updateProhibit(false), logging(false), summary(false), updateCallback(0),  induceCallback(0), client(0), historySliceCachingIs(false), historySliceCumulativeIs(false), frameUnderlyingDynamicIs(false), frameHistoryDynamicIs(false), underlyingOffsetIs(false)
{
}

inline void Alignment::Active::varPromote(SizeSizeUMap& mm, std::size_t& v)
{
	auto x = v >> this->bits << this->bits;
	auto it = mm.find(x);
	if (it != mm.end())
		v += it->second;
	else
	{
		auto y = this->system->next(this->bits);
		mm[x] = y-x;
		v += y-x;
	}
}

std::size_t Alignment::Active::varDemote(const SizeSizeUMap& mm, std::size_t v) const
{
	auto x = v >> this->bits << this->bits;
	for (auto& p : mm)
	{
		if (x == p.first + p.second)
		{
			return v - p.second;
		}
	}
	return v;
}

std::size_t Alignment::Active::varMax() const
{
	std::size_t v = this->var;
	v = std::max(v,this->varComputedMax());	
	v = std::max(v,this->varSlice);	
	if (this->decomp)
		v = std::max(v,this->decomp->varMax());
	for (auto& p : this->underlyingsVarsOffset)
		for (auto& q : p.second)	
			v = std::max(v,q.first+q.second);
	for (auto& p : this->framesVarsOffset)
		for (auto& q : p.second)	
			v = std::max(v,q.first+q.second);
	return (((v >> this->bits) + 1) << this->bits) - 1;
}

std::size_t Alignment::Active::varComputedMax() const
{
	if (!induceVarComputeds.size())
		return 0;	
	std::size_t v = *induceVarComputeds.rbegin();
    return (v << 12) + (8ull << 8) + 255 + (1ull << this->bits);
}

#define UNLOG ; log_str.flush(); this->log(*this, log_str.str());}
#define LOG { std::ostringstream log_str; log_str << (this->name.size() ? this->name + "\t" : "") << 

// event ids should be monotonic and updated no more than once
bool Alignment::Active::update(ActiveUpdateParameters pp)
{
	auto drmul = listVarValuesDecompFudSlicedRepasPathSlice_u;

	bool ok = true;
	try 
	{
		while (ok && !this->terminate)
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
					ok = ok && this->historySize > 0;
					for (auto& ev : this->underlyingEventsRepa)
						ok = ok && ev;
					for (auto& ev : this->underlyingEventsSparse)
						ok = ok && ev;
					for (auto& hr : this->underlyingHistoryRepa)
						ok = ok && hr && hr->size == this->historySize && hr->dimension > 0  && hr->evient;
					for (auto& hr : this->underlyingHistorySparse)
						ok = ok && hr && hr->size == this->historySize && hr->capacity == 1;
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
				bool continuousA = false;
				if (ok && eventsUpdatedA.size())
				{
					for (auto eventB : eventsA)
						if (eventB > *eventsUpdatedA.rbegin())
						{
							eventA = eventB;
							break;
						}
					// check for discontinuity
					continuousA = eventA == *eventsUpdatedA.rbegin() + 1;
				}
				// copy events to active history
				if (ok)
				{		
					auto& comp = this->induceVarComputeds;
					auto& slpp = this->underlyingSlicesParent;
					std::size_t block1 = (std::size_t)1 << this->bits;
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
						auto sh1 = hr1.shape;
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
						// if computed add to underlyingSlicesParent
						for (std::size_t i = 0; i < n1; i++)
							if (comp.count(vv1[i]))
							{
								std::size_t s = sh1[i];
								std::size_t b = 0; 
								if (s)
								{
									s--;
									while (s >> b)
										b++;
								}
								if (b > 1)
								{
									std::size_t v = block1 + (vv1[i] << 12) + (b << 8) + rr1[i];
                                    for (int k = (int)(b-1); k > 0 && !slpp.count(v); k--)
									{
                                        std::size_t v1 = block1 + (vv1[i] << 12) + (k << 8) + (rr1[i] >> (b-k));
										slpp[v] = v1;
										v = v1;
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
						for (int i = (int)n-1; i >= 0; i--)
							if (rr1[i])
							{
								v = rr1[i];
								break;
							}
						rr[j] = v;
						if (v && slpp.find(v) == slpp.end())
							for (int i = (int)n-1; i > 0; i--)
								if (rr1[i] && rr1[i-1])
									slpp[rr1[i]] = rr1[i-1];
					}
				}
				// check decomp exists
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
					if (ok && this->historySparse)
					{
						ok = ok && this->historySparse->size == this->historySize && this->historySparse->capacity == 1;
						if (!ok)
						{
							LOG "update\terror: inconsistent history" UNLOG
							break;
						}
					}
					SizeUCharStructList jj;
					if (ok)
					{
						auto& comp = this->induceVarComputeds;
						auto& slpp = this->underlyingSlicesParent;						
						auto promote = this->underlyingOffsetIs;
						auto& proms = this->underlyingsVarsOffset;
						std::size_t block1 = (std::size_t)1 << this->bits;
						SizeList frameUnderlyingsA(this->frameUnderlyings);
						if (!frameUnderlyingsA.size())
							frameUnderlyingsA.push_back(0);	
						std::size_t m = 0;
						for (auto& hr : this->underlyingHistoryRepa)
							m += 8*hr->dimension*frameUnderlyingsA.size();
						m += 50*this->underlyingHistorySparse.size()*frameUnderlyingsA.size();
						m += 50*this->frameHistorys.size();
						jj.reserve(m);
						auto z = this->historySize;
						auto over = this->historyOverflow;
						auto j = this->historyEvent;
						for (std::size_t g = 0; g < frameUnderlyingsA.size(); g++)
						{
							auto f = frameUnderlyingsA[g];
							if (g && !f)
								continue;
							auto& mm = this->framesVarsOffset[g];
							for (auto& hr : this->underlyingHistoryRepa)
							{
								auto n = hr->dimension;
								auto vv = hr->vectorVar;
								auto sh = hr->shape;
								auto rr = hr->arr;	
								for (std::size_t i = 0; i < n; i++)
								{
									SizeUCharStruct qq;
									qq.size = vv[i];
									if (f <= j)
										qq.uchar = rr[(j-f)*n + i];	
									else if (f && over && z > f)
										qq.uchar = rr[((j+z-f)%z)*n + i];	
									else
										qq.uchar = 0;
									if (comp.count(qq.size)) // computed
									{
										std::size_t s = sh[i];
										std::size_t b = 0; 
										if (s)
										{
											s--;
											while (s >> b)
												b++;
										}
										qq.size = block1 + (qq.size << 12) + (b << 8) + qq.uchar;
										qq.uchar = 1;
										auto it = slpp.find(qq.size);
										if (f)
											this->varPromote(mm, qq.size);
										jj.push_back(qq);
										while (it != slpp.end() && it->second)
										{
											qq.size = it->second;
											if (f)
												this->varPromote(mm, qq.size);
											jj.push_back(qq);
											it = slpp.find(it->second);
										}
									}
									else if (qq.uchar)
									{
										if (f)
											this->varPromote(mm, qq.size);
										jj.push_back(qq);
									}
								}
							}
							std::size_t h = 0;
							for (auto& hr : this->underlyingHistorySparse)
							{
								std::size_t v = 0;
								if (f <= j)
									v = hr->arr[j-f];
								else if (f && over && z > f)
									v = hr->arr[(j+z-f)%z]; 
								if (v)
								{
									{
										SizeUCharStruct qq;
										qq.uchar = 1;			
										qq.size = v;
										if (promote)
											this->varPromote(proms[h], qq.size);
										if (f)
											this->varPromote(mm, qq.size);
										jj.push_back(qq);
									}								
									auto it = slpp.find(v);
									while (it != slpp.end())
									{
										SizeUCharStruct qq;
										qq.uchar = 1;
										qq.size = it->second;
										if (!qq.size)
											break;
										if (promote)
											this->varPromote(proms[h], qq.size);
										if (f)
											this->varPromote(mm, qq.size);
										jj.push_back(qq);
										it = slpp.find(it->second);
									}										
								}
								h++;
							}										
						}
						if (ok && this->decomp && this->historySparse && this->frameHistorys.size())
						{
							auto& hr = this->historySparse;
							auto& slpp = this->decomp->mapVarParent();
							for (std::size_t g = 0; g < this->frameHistorys.size(); g++)
							{
								auto f = this->frameHistorys[g];
								if (!f)
									continue;
								auto& mm = this->framesVarsOffset[g];
								std::size_t v = 0;
								if (f <= j)
									v = hr->arr[j-f];
								else if (f && over && z > f)
									v = hr->arr[(j+z-f)%z]; 
								if (v)
								{
									{
										SizeUCharStruct qq;
										qq.uchar = 1;			
										qq.size = v;
										if (f)
											this->varPromote(mm, qq.size);			
										jj.push_back(qq);
									}								
									auto it = slpp.find(v);
									while (it != slpp.end())
									{
										SizeUCharStruct qq;
										qq.uchar = 1;
										qq.size = it->second;
										if (!qq.size)
											break;										
										if (f)
											this->varPromote(mm, qq.size);			
										jj.push_back(qq);
										it = slpp.find(it->second);
									}										
								}
							}
						}
					}
					std::unique_ptr<SizeList> ll;
					if (ok)
					{
						ll = drmul(jj,*this->decomp,(unsigned char)(pp.mapCapacity));	
						ok = ok && ll;
						if (!ok)
						{
							LOG "update\terror: drmul failed to return a list" UNLOG
						}						
					}
					// sync active slices
					if (ok && this->historySparse)
					{
						if (ll->size())
							sliceA = ll->back();
						ok = ok && this->historySparse->size == this->historySize && this->historySparse->capacity == 1;
						if (!ok)
						{
							LOG "update\terror: inconsistent history" UNLOG
							break;
						}	
						std::size_t sliceB = this->historyOverflow ? this->historySparse->arr[this->historyEvent] : 0;
						// update history
						if (!this->historyOverflow || sliceA != sliceB)
						{
							this->historySparse->arr[this->historyEvent] = sliceA;
							auto& setA = this->historySlicesSetEvent[sliceA];
							setA.insert(this->historyEvent);
							if (this->induceThreshold && setA.size() == this->induceThreshold)
								this->induceSlices.insert(sliceA);
							if (this->historyOverflow)
							{
								auto& setB = this->historySlicesSetEvent[sliceB];
								setB.erase(this->historyEvent);
								if (this->induceThreshold && setB.size() == this->induceThreshold-1)
								{
									this->induceSlices.erase(sliceB);
									this->induceSliceFailsSize.erase(sliceB);
								}
								if (!setB.size())
									this->historySlicesSetEvent.erase(sliceB);
							}
						}	
						// handle next transition
						if (this->historySliceCachingIs && !this->historySliceCumulativeIs 
							&& this->historyOverflow && this->continousIs)
						{
							auto& discont = this->continousHistoryEventsEvent;
							auto z = this->historySize;
							auto y = this->historyEvent;
							auto rs = this->historySparse->arr;
							auto& nexts = this->historySlicesSlicesSizeNext;
							auto& prevs = this->historySlicesSliceSetPrev;
							if (!discont.count((y+1)%z))
							{
								auto sliceC = rs[(y+1)%z];	
								if (sliceC != sliceB)
								{
									auto& c = nexts[sliceB][sliceC];
									if (c > 1)
										c--;
									else
									{		
										nexts[sliceB].erase(sliceC);
										if (!nexts[sliceB].size())
											nexts.erase(sliceB);
										prevs[sliceC].erase(sliceB);
										if (!prevs[sliceC].size())
											prevs.erase(sliceC);
									}
								}
							}
						}
						// handle discontinuities
						if (this->continousIs)
						{
							auto over = this->historyOverflow;
							auto& discont = this->continousHistoryEventsEvent;
							auto z = this->historySize;
							auto y = this->historyEvent;
							auto it = discont.find(y);
							if (it != discont.end() && over && y+1 < z && !discont.count(y+1))
								discont.insert_or_assign(y+1,it->second+1);
							else if (it != discont.end() && y+1 == z && !discont.count(0))
								discont.insert_or_assign(0,it->second+1);
							if (continuousA)
								discont.erase(y);
							else
								discont.insert_or_assign(y,eventA);
						}	
						// handle cached sizes and prev transition
						if (this->historySliceCachingIs)
						{
							auto cumulative = this->historySliceCumulativeIs;
							auto over = this->historyOverflow;
							auto cont = this->continousIs;
							auto& discont = this->continousHistoryEventsEvent;
							auto z = this->historySize;
							auto y = this->historyEvent;
							auto rs = this->historySparse->arr;
							auto& sizes = this->historySlicesSize;
							auto& nexts = this->historySlicesSlicesSizeNext;
							auto& prevs = this->historySlicesSliceSetPrev;
							auto& cv = this->decomp->mapVarParent();
							if (cumulative || !over || sliceA != sliceB)
							{
								auto sliceC = sliceA;
								while (true)
								{
									sizes[sliceC]++;
									if (!sliceC)
										break;
									sliceC = cv[sliceC];
								}								
							}
							if (!cumulative && over && sliceA != sliceB)
							{
								auto sliceC = sliceB;
								while (true)
								{
									auto& c = sizes[sliceC];
									if (c > 1)
										c--;
									else
										sizes.erase(sliceC);
									if (!sliceC)
										break;
									sliceC = cv[sliceC];
								}								
							}
							if ((over || y) && cont && !discont.count(y))
							{
								auto sliceC = rs[(y+z-1)%z];	
								if (sliceC != sliceA)
								{
									nexts[sliceC][sliceA]++;
									prevs[sliceA].insert(sliceC);
								}
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
						LOG "update apply\tevent id: " << eventA << "\thistory id: " << this->historyEvent << "\tslice: " << std::hex << sliceA << std::dec << "\tslice size: " << this->historySlicesSetEvent[sliceA].size() << "\ttime: " << ((sec)(clk::now() - mark)).count() << "s" UNLOG
					}
				}
				// increment historyEvent
				if (ok)
				{
					if (this->frameUnderlyingDynamicIs)
					{
						this->historyFrameUnderlying.reserve(this->historySize);
						if (this->historyEvent < this->historyFrameUnderlying.size())
							this->historyFrameUnderlying[this->historyEvent] = this->frameUnderlyings;
						else 
							this->historyFrameUnderlying.push_back(this->frameUnderlyings);
						this->historyFrameUnderlying[this->historyEvent].shrink_to_fit();
					}
					if (this->frameHistoryDynamicIs)
					{
						this->historyFrameHistory.reserve(this->historySize);
						if (this->historyEvent < this->historyFrameHistory.size())
							this->historyFrameHistory[this->historyEvent] = this->frameHistorys;
						else 
							this->historyFrameHistory.push_back(this->frameHistorys);
						this->historyFrameHistory[this->historyEvent].shrink_to_fit();
					}
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
				ok = ok && updateCallback(*this,eventsA,eventA,historyEventA,sliceA);
			}
		}
	} 
	catch (const std::exception& e) 
	{
		LOG "update error: " << e.what()  UNLOG
		ok = false;
	}
	if (!ok)
		this->terminate = true;

	return ok;
}

void run_induce(Active& active, std::size_t sliceA, ActiveInduceParameters pp, ActiveUpdateParameters ppu)
{
	active.induce(sliceA, pp, ppu);
	return;
};

bool Alignment::Active::induce(ActiveInduceParameters pp, ActiveUpdateParameters ppu)
{		
	bool ok = true;
	
	try 
	{
		if (!pp.asyncThreadMax) // run synchronously
		{
			while (ok && !this->terminate)
			{
				std::size_t sliceA = 0;
				std::size_t sliceSizeA = 0;	
				for (auto sliceB : this->induceSlices)
				{
					auto sliceSizeB = this->historySlicesSetEvent[sliceB].size();
					if (sliceSizeB > sliceSizeA)
					{
						auto it = this->induceSliceFailsSize.find(sliceB);
						if (it == this->induceSliceFailsSize.end() 
							|| (it->second < sliceSizeB 
								&& pp.induceThresholdExceeded(it->second, sliceSizeB)))
						{
							sliceA = sliceB;
							sliceSizeA = sliceSizeB;							
						}
					}
				}		
				if (ok && sliceSizeA) 
					this->induce(sliceA, pp, ppu);
				else
					break;		
			}			
		}
		else // run asynchronously
		{
			std::map<std::size_t,std::thread> threads;
			while (ok && !this->terminate)
			{
				// join any threads that have finished
				if (ok && threads.size())
				{
					SizeSet inducingSlicesA;
					{
						std::lock_guard<std::mutex> guard(this->mutex);		
						inducingSlicesA = this->inducingSlices;
					}
					SizeSet threadSlicesA;		
					for (auto& pp : threads)
						if (!inducingSlicesA.count(pp.first))
						{
							pp.second.join();
							threadSlicesA.insert(pp.first);
						}
					for (auto sliceA : threadSlicesA)		
						threads.erase(sliceA);
				}
				// get largest slice
				std::size_t sliceSizeMax = 0;	
				if (ok && threads.size() < pp.asyncThreadMax)
				{
					std::size_t sliceA = 0;
					std::size_t sliceSizeA = 0;	
					for (auto sliceB : this->induceSlices)
					{
						if (pp.asyncUpdateLimit || !threads.count(sliceB))
						{
							auto sliceSizeB = this->historySlicesSetEvent[sliceB].size();
							if (sliceSizeB > sliceSizeA)
							{
								auto it = this->induceSliceFailsSize.find(sliceB);
								if (it == this->induceSliceFailsSize.end() 
									|| (it->second < sliceSizeB 
										&& pp.induceThresholdExceeded(it->second, sliceSizeB)))
								{
									if (!threads.count(sliceB))
									{
										sliceA = sliceB;
										sliceSizeA = sliceSizeB;				
									}
									if (sliceSizeB > sliceSizeMax)
										sliceSizeMax = sliceSizeB;
								}
							}
						}
					}	
					if (ok && sliceSizeA) 
					{
						{
							std::lock_guard<std::mutex> guard(this->mutex);		
							this->inducingSlices.insert(sliceA);
						}
						threads.insert_or_assign(sliceA,std::thread(run_induce, std::ref(*this), sliceA, pp, ppu));
					}
				}
				if (ok)
					this->updateProhibit = 
						(pp.asyncUpdateLimit && sliceSizeMax > pp.asyncUpdateLimit) 
						|| threads.size() == pp.asyncThreadMax;
				if (ok && pp.asyncInterval)
					std::this_thread::sleep_for(std::chrono::milliseconds(pp.asyncInterval));
				else if (ok)
					std::this_thread::yield();
				else
					break;
			}
			if (ok)
			{
				for (auto& pp : threads)
					pp.second.join();
			}
		}
	} 
	catch (const std::exception& e) 
	{
		LOG "induce error: " << e.what()  UNLOG
		ok = false;
	}
	if (!ok)
		this->terminate = true;
	
	return ok;
}

bool Alignment::Active::induce(std::size_t sliceA, ActiveInduceParameters pp, ActiveUpdateParameters ppu)
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
	auto drmul = listVarValuesDecompFudSlicedRepasPathSlice_u;
	auto layerer = parametersLayererMaxRollByMExcludedSelfHighestLogIORepa_up;
		
	bool ok = true;
	try 
	{
		if (ok && !this->terminate)
		{
			std::size_t varA = 0;
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
					ok = ok && this->historySize > 0;
					ok = ok && (llr.size() || lla.size());
					for (auto& hr : llr)
						ok = ok && hr && hr->size == this->historySize && hr->dimension > 0 && hr->evient;
					for (auto& hr : lla)
						ok = ok && hr && hr->size == this->historySize && hr->capacity == 1;
					if (!ok)
					{
						LOG "induce\tslice: " << std::hex << sliceA << std::dec << "\terror: inconsistent underlying" UNLOG
					}	
				}			
				// get slice size
				if (ok)
				{
					sliceSizeA = this->historySlicesSetEvent[sliceA].size();
					ok = ok && sliceSizeA;
				}
				// copy events from evient active history to varient selection
				if (ok)
				{
					varA = this->var;
					auto& setEventsA = this->historySlicesSetEvent[sliceA];
					eventsA.insert(eventsA.end(),setEventsA.begin(),setEventsA.end());
					SizeList frameUnderlyingsA(this->frameUnderlyings);
					if (!frameUnderlyingsA.size())
						frameUnderlyingsA.push_back(0);					
					SizeSet qqc;
					if (ok && llr.size())
					{
						auto& comp = this->induceVarComputeds;
						auto& excl = this->induceVarExclusions;					
						for (auto& hr : llr)
						{
							auto n = hr->dimension;
							auto vv = hr->vectorVar;
							for (std::size_t i = 0; i < n; i++)
							{
								auto v = vv[i];
								if (!excl.count(v))
								{
									if (comp.count(v))
										qqc.insert(v);
									else
										qqr.insert(v);
								}
							}
						}
					}
					if (ok && qqr.size())
					{
						hrr = std::make_unique<HistoryRepa>();
						hrr->dimension = qqr.size()*frameUnderlyingsA.size();
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
						auto z = this->historySize;
						auto over = this->historyOverflow;
						std::size_t i = 0;
						for (std::size_t g = 0; g < frameUnderlyingsA.size(); g++)
						{
							auto f = frameUnderlyingsA[g];
							auto& mm = this->framesVarsOffset[g];
							for (auto v : qqr)
							{
								vvr[i] = v;
								if (g || f)
									this->varPromote(mm, vvr[i]);
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
										{
											if (this->frameUnderlyingDynamicIs)
											{
												f = 0;
												auto& frameUnderlyingsB = this->historyFrameUnderlying[ev[j]];
												if (g < frameUnderlyingsB.size())
													f = frameUnderlyingsB[g];
											}
											if (g && !f)
												rrr[izr + j] = 0;					
											else if (f <= ev[j])
												rrr[izr + j] = rr[(ev[j]-f)*n + k];
											else if (f && over && z > f)
												rrr[izr + j] = rr[((ev[j]+z-f)%z)*n + k];
											else
												rrr[izr + j] = 0;					
										}
										break;
									}
								}
								i++;
							}
						}
						qqr.clear();
						for (i = 0; i < nr; i++)
							qqr.insert(vvr[i]);
					}
					if (ok && (lla.size() || qqc.size() || (this->decomp && this->historySparse && this->frameHistorys.size())))
					{
						auto za = eventsA.size(); 
						auto na = (lla.size()+qqc.size())*frameUnderlyingsA.size();
						if (this->decomp && this->historySparse)
							na += this->frameHistorys.size();
						auto ev = eventsA.data();
						auto z = this->historySize;
						auto over = this->historyOverflow;
						auto promote = this->underlyingOffsetIs;
						auto& proms = this->underlyingsVarsOffset;
						haa = std::make_unique<HistorySparseArray>();
						haa->size = za;
						haa->capacity = na;
						haa->arr = new std::size_t[za*na];
						auto raa = haa->arr;
						slppa.reserve(za*na*4);
						std::size_t i = 0;
						if (ok && lla.size())
						{
							auto& slpp = this->underlyingSlicesParent;
							for (std::size_t g = 0; g < frameUnderlyingsA.size(); g++)
							{
								auto f = frameUnderlyingsA[g];
								auto& mm = this->framesVarsOffset[g];
								std::size_t h = 0;
								for (auto& hr : lla)
								{
									auto rr = hr->arr;
									for (std::size_t j = 0; j < za; j++)
									{
										if (this->frameUnderlyingDynamicIs)
										{
											f = 0;
											auto& frameUnderlyingsB = this->historyFrameUnderlying[ev[j]];
											if (g < frameUnderlyingsB.size())
												f = frameUnderlyingsB[g];
										}
										std::size_t v = 0;
										if (g && !f)
											v = 0;					
										else if (f <= ev[j])
											v = rr[ev[j]-f];
										else if (f && over && z > f)
											v = rr[(ev[j]+z-f)%z];
										raa[j*na + i] = v;
										if (v)
										{
											if (promote)
												this->varPromote(proms[h], raa[j*na + i]);
											if (f)
												this->varPromote(mm, raa[j*na + i]);
											auto it = slpp.find(v);
											while (it != slpp.end())
											{
												auto w1 = it->first;
												if (promote)
													this->varPromote(proms[h], w1);
												if (f)
													this->varPromote(mm, w1);
												auto w2 = it->second;
												if (!w2)
													break;
												if (promote)
													this->varPromote(proms[h], w2);
												if (f)
													this->varPromote(mm, w2);
												slppa.insert_or_assign(w1, w2);
												it = slpp.find(it->second);
											}
										}
									}	
									h++;
									i++;								
								}
							}							
						}
						if (ok && this->decomp && this->historySparse && this->frameHistorys.size())
						{
							auto& hr = this->historySparse;
							auto& slpp = this->decomp->mapVarParent();
							for (std::size_t g = 0; g < this->frameHistorys.size(); g++)
							{
								auto f = this->frameHistorys[g];
								auto& mm = this->framesVarsOffset[g];
								auto rr = hr->arr;
								for (std::size_t j = 0; j < za; j++)
								{
									if (this->frameHistoryDynamicIs)
									{
										f = 0;
										auto& frameHistorysB = this->historyFrameHistory[ev[j]];
										if (g < frameHistorysB.size())
											f = frameHistorysB[g];
									}
									std::size_t v = 0;
									if (!f)
										v = 0;									
									else if (f <= ev[j])
										v = rr[ev[j]-f];
									else if (f && over && z > f)
										v = rr[(ev[j]+z-f)%z];
									raa[j*na + i] = v;
									if (v)
									{
										if (f)
											this->varPromote(mm, raa[j*na + i]);
										auto it = slpp.find(v);
										while (it != slpp.end())
										{
											auto w1 = it->first;
											if (f)
												this->varPromote(mm, w1);
											auto w2 = it->second;
											if (!w2)
												break;
											if (f)
												this->varPromote(mm, w2);
											slppa.insert_or_assign(w1, w2);
											it = slpp.find(it->second);
										}
									}
								}
								i++;
							}									
						}							
						if (ok && qqc.size())
						{
							auto& slpp = this->underlyingSlicesParent;
							std::size_t block1 = (std::size_t)1 << this->bits;
							for (std::size_t g = 0; g < frameUnderlyingsA.size(); g++)
							{
								auto f = frameUnderlyingsA[g];
								auto& mm = this->framesVarsOffset[g];
								for (auto v : qqc)
								{
									for (auto& hr : llr)
									{
										auto& mvv = hr->mapVarInt();
										auto it = mvv.find(v);
										if (it != mvv.end())
										{
											auto n = hr->dimension;
											auto rr = hr->arr;
											auto k = it->second;
											auto s = hr->shape[k];
											std::size_t b = 0; 
											if (s)
											{
												s--;
												while (s >> b)
													b++;
											}
											for (std::size_t j = 0; j < za; j++)
											{
												if (this->frameUnderlyingDynamicIs)
												{
													f = 0;
													auto& frameUnderlyingsB = this->historyFrameUnderlying[ev[j]];
													if (g < frameUnderlyingsB.size())
														f = frameUnderlyingsB[g];
												}
												unsigned char u = 0;
												if (g && !f)
													u = 0;					
												else if (f <= ev[j])
													u = rr[(ev[j]-f)*n + k];
												else if (f && over && z > f)
													u = rr[((ev[j]+z-f)%z)*n + k];	
												std::size_t w = block1 + (v << 12) + (b << 8) + u;
												auto it = slpp.find(w);
												if (f)
													this->varPromote(mm, w);
												raa[j*na + i] = w;
												while (it != slpp.end() && it->second)
												{
													auto w1 = it->first;
													if (f)
														this->varPromote(mm, w1);
													auto w2 = it->second;
													if (f)
														this->varPromote(mm, w2);
													slppa.insert_or_assign(w1, w2);
													it = slpp.find(it->second);
												}
											}
											break;
										}
									}
									i++;
								}
							}
						}
					}
				}
				if (ok && this->logging)
				{
					LOG "induce copy\tslice: " << std::hex << sliceA << std::dec << "\tslice size: " << sliceSizeA << "\trepa dimension: " << (hrr ? hrr->dimension : 0) << "\tsparse capacity: " << (haa ? haa->capacity : 0) << "\tsparse paths: " << slppa.size() << "\tvariable: " << varA << "\ttime: " << ((sec)(clk::now() - mark)).count() << "s" UNLOG
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
					LOG "induce\tslice: " << std::hex << sliceA << std::dec << "\terror: inconsistent copy" UNLOG
				}	
			}
			bool fail = false;
			std::unique_ptr<HistoryRepa> hr;
			std::unique_ptr<FudRepa> fr;
			std::size_t frSize = 0;
			SizeList kk;
			double algn = 0.0;
			double diagonal = 0.0;
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
							for (int i = (int)(ll.size()) - 1; i > 0; i--)
								for (int m = i-1; m >= 0; m--)
									mma[ll[i]].insert(ll[m]);
						}
					}
				}
				// get top nmax vars by entropy
				// remove any sparse parents with same entropy as children
				if (ok && (qqr.size() || qqa.size()))
				{
					auto nmax = (std::size_t)std::sqrt(pp.znnmax / (double)(2*sliceSizeA));
					nmax = std::max(nmax, pp.bmax);
					DoubleSizePairList ee;
					ee.reserve(qqr.size() + qqa.size());
					if (qqr.size())
					{
						SizeList vv(qqr.begin(),qqr.end());
						auto eer = prents(*hrpr(vv.size(), vv.data(), *hrr));
						for (auto p : *eer)
							if (p.first > repaRounding)
								ee.push_back(DoubleSizePair(-p.first,p.second));
					}
					if (qqa.size())
					{
						std::map<std::size_t, SizeSet> eem0;
						for (auto p : qqa)	
							if (p.second > 0 && p.second < sliceSizeA)		
								eem0[p.second].insert(p.first);
						std::map<std::size_t, SizeSet> eem;
						for (auto p : eem0)		
						{
							auto e = p.first;
							auto& xx = p.second;
							for (auto v : xx)
							{
								auto iv = mma.find(v);
								if (iv != mma.end())
								{
									bool found = false;
									for (auto w : xx)
										if (iv->second.count(w))
										{
											found = true;
											break;
										}
									if (!found)								
										eem[e].insert(v);				
								}			
								else
									eem[e].insert(v);
							}
						}	
						double f = 1.0/(double)sliceSizeA;
						for (auto& p : eem)	
						{
							double a = (double)p.first * f;
							double e = -(a * std::log(a) + (1.0-a) * std::log(1.0-a));
							if (e > repaRounding)
								for (auto& q : p.second)	
									ee.push_back(DoubleSizePair(-e,q));
						}
					}
					if (ee.size()) std::sort(ee.begin(), ee.end());
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
						std::lock_guard<std::mutex> guard(this->mutex);	
						if (!fail)
						{
							LOG "induce model\tslice: " << std::hex << sliceA << std::dec << "\trepa dimension: " << qqr.size() << "\tsparse dimension: " << qqa.size() UNLOG
						}
						else
						{
							LOG "induce model\tslice: " << std::hex << sliceA << std::dec << "\tno entropy"  UNLOG
						}
					}	
				}	
				if (ok && !fail)
				{
					std::unique_ptr<HistoryRepa> hrs;
					std::ranlux48_base gen((unsigned int)(pp.seed+sliceA));
					if (ok && qqr.size())
					{
						if (qqr.size() < hrr->dimension)
						{
							SizeList vv(qqr.begin(),qqr.end());
							hr = hrhrred(vv.size(), vv.data(), *hrr);
						}
						else
							hr = std::move(hrr);
						hrs = hrshuffle(*hr,gen);
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
							LOG "induce\tslice: " << std::hex << sliceA << std::dec << "\terror: inconsistent reduction" UNLOG
						}	
					}
					if (ok && this->logging)
					{
						std::lock_guard<std::mutex> guard(this->mutex);	
						LOG "induce model\tslice: " << std::hex << sliceA << std::dec << "\tdimension: " << hr->dimension << "\tsize: " << hr->size UNLOG
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
							auto t = layerer(pp.wmax, pp.lmax, pp.xmax, pp.omax, pp.bmax, pp.mmax, pp.umax, pp.pmax, pp.tint, vv, *hr, *hrs, layerer_log, this->logging && pp.logging, varA);
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
								diagonal = 100.0*(algn/z/(m-1)/exp(1.0));
								fail = pp.diagonalMin > 0.0 && diagonal < pp.diagonalMin;
							}
						}
						catch (const std::out_of_range& e)
						{
							ok = false;
							LOG "induce\tslice: " << std::hex << sliceA << std::dec << "\tout of range exception: " << e.what() UNLOG
						}
						if (ok && this->logging)
						{
							std::lock_guard<std::mutex> guard(this->mutex);	
							if (!fail)
							{
								LOG "induce model\tslice: " << std::hex << sliceA << std::dec << "\tder vars algn density: " << algn << "\timpl bi-valency percent: " << diagonal << "\tder vars cardinality: " << kk.size() << "\tfud cardinality: " << frSize UNLOG							
							}
							else
							{
								LOG "induce model\tslice: " << std::hex << sliceA << std::dec << "\tno alignment\tder vars algn density: " << algn << "\timpl bi-valency percent: " << diagonal UNLOG
							}
						}	
					}
				}
				if (ok && this->logging)
				{
					std::lock_guard<std::mutex> guard(this->mutex);	
					LOG "induce model\tslice: " << std::hex << sliceA << std::dec << "\ttime: " << ((sec)(clk::now() - mark)).count() << "s" UNLOG
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
						LOG "induce\tslice: " << std::hex << sliceA << std::dec << "\terror: no system" UNLOG
					}	
				}
				// remap kk and fr with block ids
				if (ok)
				{
					if (frSize > ((std::size_t)1 << this->bits))
					{
						ok = false;
						LOG "induce\tslice: " << std::hex << sliceA << std::dec << "\terror: block too small" << "\tfud size: " << frSize << "\tblock: " << (1 << this->bits) UNLOG
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
				std::size_t v = 0;
				SizeList sl;
				// create the slices
				if (ok)
				{		
					if (!this->decomp)
						this->decomp = std::make_unique<DecompFudSlicedRepa>();
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
					if (sz > ((std::size_t)1 << this->bits))
					{
						ok = false;
						LOG "induce\tslice: " << std::hex << sliceA << std::dec << "\terror: block too small" << "\tslice size: " << sz << "\tblock: " << (1 << this->bits) UNLOG
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
						if (sliceA)
						{
							tr->dimension = m + 1;
							tr->vectorVar = new std::size_t[m + 1];
							auto ww = tr->vectorVar;
							tr->shape = new std::size_t[m + 1];
							auto sh = tr->shape;
							ww[0] = sliceA;
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
						if (sliceA)
						{
							tr->dimension = m + 1;
							tr->vectorVar = new std::size_t[m + 1];
							auto ww = tr->vectorVar;
							tr->shape = new std::size_t[m + 1];
							auto sh = tr->shape;
							ww[0] = sliceA;
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
					dr.fudRepasSize += fs.fud.size();
					auto& vi = dr.mapVarInt();
					vi[sliceA] = dr.fuds.size() - 1;
					auto& cv = dr.mapVarParent();
					for (auto s : sl)
						cv[s] = sliceA;
				}
				// check historySparse
				if (ok)
				{
					ok = ok && this->historySparse && this->historySparse->arr;
					if (!ok)
					{
						LOG "induce update\tslice: " << std::hex << sliceA << std::dec << "\terror: historySparse not initialised" UNLOG
					}
				}
				// update historySparse and historySlicesSetEvent
				if (ok)
				{
					if (sliceA)
					{
						auto z = hr->size;
						auto hrr = std::make_unique<HistoryRepa>();
						hrr->dimension = 1;
						hrr->vectorVar = new std::size_t[1];
						hrr->vectorVar[0] = sliceA;
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
								this->historySparse->arr[eventA] = sliceB;
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
							SizeUCharStructList jj;
							if (ok)
							{
								auto& comp = this->induceVarComputeds;
								auto& slpp = this->underlyingSlicesParent;
								auto promote = this->underlyingOffsetIs;
								auto& proms = this->underlyingsVarsOffset;		
								std::size_t block1 = (std::size_t)1 << this->bits;
								SizeList frameUnderlyingsA(this->frameUnderlyings);
								if (!frameUnderlyingsA.size())
									frameUnderlyingsA.push_back(0);	
								std::size_t m = 0;
								for (auto& hr : this->underlyingHistoryRepa)
									m += hr->dimension*frameUnderlyingsA.size();
								m += 50*this->underlyingHistorySparse.size()*frameUnderlyingsA.size();
								m += 50*this->frameHistorys.size();
								jj.reserve(m);
								auto z = this->historySize;
								auto over = this->historyOverflow;
								auto j = eventB;
								for (std::size_t g = 0; g < frameUnderlyingsA.size(); g++)
								{
									auto f = frameUnderlyingsA[g];
									if (this->frameUnderlyingDynamicIs)
									{
										auto& frameUnderlyingsB = this->historyFrameUnderlying[j];
										if (g < frameUnderlyingsB.size())
											f = frameUnderlyingsB[g];
										else
											f = 0;
									}
									if (g && !f)
										continue;
									auto& mm = this->framesVarsOffset[g];
									for (auto& hr : this->underlyingHistoryRepa)
									{
										auto n = hr->dimension;
										auto vv = hr->vectorVar;
										auto sh = hr->shape;
										auto rr = hr->arr;	
										for (std::size_t i = 0; i < n; i++)
										{
											SizeUCharStruct qq;
											qq.size = vv[i];
											if (f <= j)
												qq.uchar = rr[(j-f)*n + i];	
											else if (f && over && z > f)
												qq.uchar = rr[((j+z-f)%z)*n + i];	
											else
												qq.uchar = 0;
											if (comp.count(qq.size)) // computed
											{
												std::size_t s = sh[i];
												std::size_t b = 0; 
												if (s)
												{
													s--;
													while (s >> b)
														b++;
												}
												qq.size = block1 + (qq.size << 12) + (b << 8) + qq.uchar;
												qq.uchar = 1;
												auto it = slpp.find(qq.size);
												if (f)
													this->varPromote(mm, qq.size);
												jj.push_back(qq);
												while (it != slpp.end() && it->second)
												{
													qq.size = it->second;
													if (f)
														this->varPromote(mm, qq.size);
													jj.push_back(qq);
													it = slpp.find(it->second);
												}
											}
											else if (qq.uchar)
											{
												if (f)
													this->varPromote(mm, qq.size);
												jj.push_back(qq);
											}
										}							
									}
									auto& slpp = this->underlyingSlicesParent;
									std::size_t h = 0;									
									for (auto& hr : this->underlyingHistorySparse)
									{
										std::size_t v = 0;
										if (f <= j)
											v = hr->arr[j-f];
										else if (f && over && z > f)
											v = hr->arr[(j+z-f)%z]; 
										if (v)
										{
											{
												SizeUCharStruct qq;
												qq.uchar = 1;			
												qq.size = v;
												if (promote)
													this->varPromote(proms[h], qq.size);
												if (f)
													this->varPromote(mm, qq.size);
												jj.push_back(qq);
											}								
											auto it = slpp.find(v);
											while (it != slpp.end())
											{
												SizeUCharStruct qq;
												qq.uchar = 1;
												qq.size = it->second;
												if (!qq.size)
													break;
												if (promote)
													this->varPromote(proms[h], qq.size);
												if (f)
													this->varPromote(mm, qq.size);
												jj.push_back(qq);
												it = slpp.find(it->second);
											}										
										}
										h++;
									}										
								}
								if (ok && this->decomp && this->historySparse && this->frameHistorys.size())
								{
									auto& hr = this->historySparse;
									auto& slpp = this->decomp->mapVarParent();
									for (std::size_t g = 0; g < this->frameHistorys.size(); g++)
									{
										std::size_t f = this->frameHistorys[g];
										if (this->frameHistoryDynamicIs)
										{
											auto& frameHistorysB = this->historyFrameHistory[j];
											if (g < frameHistorysB.size())
												f = frameHistorysB[g];
											else
												f = 0;
										}
										if (!f)
											continue;
										auto& mm = this->framesVarsOffset[g];
										std::size_t v = 0;
										if (f <= j)
											v = hr->arr[j-f];
										else if (f && over && z > f)
											v = hr->arr[(j+z-f)%z]; 
										if (v)
										{
											{
												SizeUCharStruct qq;
												qq.uchar = 1;			
												qq.size = v;
												if (f)
													this->varPromote(mm, qq.size);
												jj.push_back(qq);
											}								
											auto it = slpp.find(v);
											while (it != slpp.end())
											{
												SizeUCharStruct qq;
												qq.uchar = 1;
												qq.size = it->second;
												if (!qq.size)
													break;
												if (f)
													this->varPromote(mm, qq.size);
												jj.push_back(qq);
												it = slpp.find(it->second);
											}										
										}
									}
								}
							}
							std::unique_ptr<SizeList> ll;
							if (ok)
							{
								ll = drmul(jj,*this->decomp,(unsigned char)(ppu.mapCapacity));	
								ok = ok && ll && ll->size() && ll->back();
								if (!ok)
								{
									LOG "induce update\tslice: " << std::hex << sliceA << std::dec << "\terror: drmul failed to return a list" UNLOG
									break;
								}						
							}							
							if (ok)
							{
								std::size_t	sliceB = ll->back();
								this->historySparse->arr[eventB] = sliceB;
								this->historySlicesSetEvent[sliceB].insert(eventB);	
								slices.insert(sliceB);									
							}
						}
						if (ok)
						{
							for (auto sliceB : slices)
								if (this->induceThreshold && this->historySlicesSetEvent[sliceB].size() >= induceThreshold)
									this->induceSlices.insert(sliceB);
						}
					}
					this->historySlicesSetEvent.erase(sliceA);
				}
				// handle cached sizes and transitions
				if (ok && this->historySliceCachingIs)
				{
					auto over = this->historyOverflow;
					auto cont = this->continousIs;
					auto& discont = this->continousHistoryEventsEvent;
					auto z = this->historySize;
					auto y = this->historyEvent;
					auto rs = this->historySparse->arr;
					auto& slices = this->historySlicesSetEvent;
					auto& sizes = this->historySlicesSize;
					auto& lengths = this->historySlicesLength;
					auto& nexts = this->historySlicesSlicesSizeNext;
					auto& prevs = this->historySlicesSliceSetPrev;					
					SizeSet events;
					for (auto sliceB : sl)
					{
						lengths[sliceB] = lengths[sliceA] + 1;
						auto slicesIt = slices.find(sliceB);
						if (slicesIt != slices.end())
						{
							sizes[sliceB] = slicesIt->second.size();
							if (cont)
								events.insert(slicesIt->second.begin(),slicesIt->second.end());
						}
					}
					if (cont)
					{
						for (auto sliceC : prevs[sliceA])
						{
							nexts[sliceC].erase(sliceA);
							if (!nexts[sliceC].size())
								nexts.erase(sliceC);
						}
						for (auto pp : nexts[sliceA])
						{
							auto sliceC = pp.first;
							prevs[sliceC].erase(sliceA);
							if (!prevs[sliceC].size())
								prevs.erase(sliceC);
						}
						prevs.erase(sliceA);
						nexts.erase(sliceA);
						for (auto ev : events)
						{
							auto sliceB = rs[ev];	
							if ((over || ev) && ev != y && !discont.count(ev))
							{
								auto sliceC = rs[(ev+z-1)%z];	
								if (sliceC != sliceB)
								{
									nexts[sliceC][sliceB]++;
									prevs[sliceB].insert(sliceC);
								}
							}
							if (!events.count((ev+1)%z) 
								&& ev != y-1 && !discont.count((ev+1)%z))
							{
								auto sliceC = rs[(ev+1)%z];		
								if (sliceC != sliceB)
								{
									nexts[sliceB][sliceC]++;
									prevs[sliceC].insert(sliceB);
								}
							}								
						}						
					}
				}				
				// remove from inducingSlices if running async
				if (ok && pp.asyncThreadMax)
				{
					this->inducingSlices.erase(sliceA);
				}
				if (ok && this->logging)
				{
					LOG "induce update\tslice: " << std::hex << sliceA << std::dec << "\tparent slice: " << v << "\tchildren cardinality: " << sl.size() << "\tfud size: " << this->decomp->fuds.back().fud.size() << "\tfud cardinality: " << this->decomp->fuds.size() << "\tmodel cardinality: " << this->decomp->fudRepasSize << "\ttime: " << ((sec)(clk::now() - mark)).count() << "s" UNLOG
				}	
				if (ok && this->summary)
				{
					std::size_t sizeA = this->historyOverflow ? this->historySize : this->historyEvent;
					if (sizeA)
					{
						const std::time_t t_c = std::chrono::system_clock::to_time_t(clk::now());
						std::string ts(std::ctime(&t_c));
						ts.pop_back();
						LOG "induce summary\tslice: " << std::hex << sliceA << std::dec << "\tdiagonal: " << diagonal << "\tfud cardinality: " << this->decomp->fuds.size() << "\tmodel cardinality: " << this->decomp->fudRepasSize<< "\tfuds per threshold: " << (double)this->decomp->fuds.size() * this->induceThreshold / sizeA << "\tat: " << ts.c_str() UNLOG
					}
				}	
				if (ok && induceCallback)
				{
					ok = ok && induceCallback(*this,sliceA,sliceSizeA);
				}
			}
			if (ok && fail)
			{
				auto mark = (ok && this->logging) ? clk::now() : std::chrono::time_point<clk>();
				std::lock_guard<std::mutex> guard(this->mutex);		
				this->induceSliceFailsSize.insert_or_assign(sliceA, sliceSizeA);
				// remove from inducingSlices if running async
				if (ok && pp.asyncThreadMax)
				{
					this->inducingSlices.erase(sliceA);
				}
				if (ok && this->logging)
				{
					LOG "induce update fail\tslice: " << std::hex << sliceA << std::dec << "\tslice size: " << sliceSizeA << "\tfails: " << this->induceSliceFailsSize  << "\ttime: " << ((sec)(clk::now() - mark)).count() << "s" UNLOG
				}					
			}
		}
	} 
	catch (const std::exception& e) 
	{
		LOG "induce error\tslice: " << std::hex << sliceA << std::dec << " : " << e.what()  UNLOG
		ok = false;
	}
	if (!ok)
		this->terminate = true;
	
	return ok;
}

bool Alignment::Active::dump(const ActiveIOParameters& pp)
{
	bool ok = true;
	
	std::ofstream out;
	try 
	{
		auto mark = (ok && this->logging) ? clk::now() : std::chrono::time_point<clk>();
		out.exceptions(out.failbit | out.badbit);
		std::lock_guard<std::mutex> guard(this->mutex);
		out.open(pp.filename, std::ios::binary);
		if (ok)
		{		
			std::size_t h = this->name.size();
			out.write(reinterpret_cast<char*>(&h), sizeof(std::size_t));
			out.write(reinterpret_cast<char*>((char*)this->name.data()), h);
		}
		if (ok)
		{		
			std::size_t h = this->underlyingEventUpdateds.size();
			out.write(reinterpret_cast<char*>(&h), sizeof(std::size_t));
			if (ok && h) 
			{
				std::size_t ev = *this->underlyingEventUpdateds.rbegin();
				out.write(reinterpret_cast<char*>(&ev), sizeof(std::size_t));
			}
		}
		if (ok)
		{		
			out.write(reinterpret_cast<char*>(&this->historySize), sizeof(std::size_t));
			out.write(reinterpret_cast<char*>(&this->historyOverflow), 1);
			out.write(reinterpret_cast<char*>(&this->historyEvent), sizeof(std::size_t));
		}
		if (ok)
		{		
			std::size_t hsize = this->underlyingHistoryRepa.size();
			out.write(reinterpret_cast<char*>(&hsize), sizeof(std::size_t));
			if (ok && hsize) 
			{
				for (std::size_t h = 0; ok && h < hsize; h++)	
				{
					auto& hr = this->underlyingHistoryRepa[h];
					ok = ok && hr;
					if (!ok)
					{
						LOG "dump error:\tfailed to write undefined underlying history repa to file: " << pp.filename  UNLOG
						break;
					}
					if (this->historyOverflow)
						historyRepasPersistent(*hr, out);
					else
						historyRepasPersistentInitial(*hr, this->historyEvent, out);
				}
			}		
		}
		if (ok)
		{		
			std::size_t hsize = this->underlyingHistorySparse.size();
			out.write(reinterpret_cast<char*>(&hsize), sizeof(std::size_t));
			if (ok && hsize) 
			{
				for (std::size_t h = 0; ok && h < hsize; h++)	
				{
					auto& hr = this->underlyingHistorySparse[h];
					ok = ok && hr;
					if (!ok)
					{
						LOG "dump error:\tfailed to write undefined underlying history sparse to file: " << pp.filename  UNLOG
						break;
					}
					if (this->historyOverflow)
						historySparseArraysPersistent(*hr, out);
					else
						historySparseArraysPersistentInitial(*hr, this->historyEvent, out);
				}
			}					
		}
		if (ok)
		{		
			std::size_t hsize = this->underlyingSlicesParent.size();
			out.write(reinterpret_cast<char*>(&hsize), sizeof(std::size_t));
			if (ok && hsize) 
			{
				for (auto& p : this->underlyingSlicesParent)	
				{
					out.write(reinterpret_cast<char*>((std::size_t*)&p.first), sizeof(std::size_t));
					out.write(reinterpret_cast<char*>((std::size_t*)&p.second), sizeof(std::size_t));
				}
			}					
		}
		if (ok)
		{		
			bool has = this->decomp ? true : false;
			out.write(reinterpret_cast<char*>(&has), 1);
			if (ok && has) 
				decompFudSlicedRepasPersistent(*this->decomp, out);
		}
		if (ok)
		{		
			out.write(reinterpret_cast<char*>(&this->bits), sizeof(int));
			out.write(reinterpret_cast<char*>(&this->var), sizeof(std::size_t));
			out.write(reinterpret_cast<char*>(&this->varSlice), sizeof(std::size_t));
			out.write(reinterpret_cast<char*>(&this->induceThreshold), sizeof(std::size_t));
			std::size_t hsize = this->induceVarExclusions.size();
			out.write(reinterpret_cast<char*>(&hsize), sizeof(std::size_t));
			if (ok && hsize) 
			{
				for (auto v : this->induceVarExclusions)	
					out.write(reinterpret_cast<char*>((std::size_t*)&v), sizeof(std::size_t));
			}
		}
		if (ok)
		{		
			bool has = this->historySparse ? true : false;
			out.write(reinterpret_cast<char*>(&has), 1);
			if (ok && has) 
			{
				auto& hr = this->historySparse;
				if (this->historyOverflow)
					historySparseArraysPersistent(*hr, out);
				else
					historySparseArraysPersistentInitial(*hr, this->historyEvent, out);
			}
		}
		if (ok)
		{		
			std::size_t hsize = this->frameUnderlyings.size();
			out.write(reinterpret_cast<char*>(&hsize), sizeof(std::size_t));
			for (auto v : this->frameUnderlyings)	
				out.write(reinterpret_cast<char*>((std::size_t*)&v), sizeof(std::size_t));
		}
		if (ok)
		{		
			std::size_t hsize = this->frameHistorys.size();
			out.write(reinterpret_cast<char*>(&hsize), sizeof(std::size_t));
			for (auto v : this->frameHistorys)	
				out.write(reinterpret_cast<char*>((std::size_t*)&v), sizeof(std::size_t));
		}
		if (ok)
		{		
			std::size_t hsize = this->framesVarsOffset.size();
			out.write(reinterpret_cast<char*>(&hsize), sizeof(std::size_t));
			for (auto mm : this->framesVarsOffset)	
			{
				std::size_t i = mm.first;
				out.write(reinterpret_cast<char*>(&i), sizeof(std::size_t));
				std::size_t msize = mm.second.size();
				out.write(reinterpret_cast<char*>(&msize), sizeof(std::size_t));
				for (auto p : mm.second)	
				{
					out.write(reinterpret_cast<char*>((std::size_t*)&p.first), sizeof(std::size_t));
					out.write(reinterpret_cast<char*>((std::size_t*)&p.second), sizeof(std::size_t));						
				}
			}
		}
		if (ok)
		{		
			out.write(reinterpret_cast<char*>(&this->continousIs), 1);
			if (this->continousIs)
			{
				std::size_t hsize = this->continousHistoryEventsEvent.size();
				out.write(reinterpret_cast<char*>(&hsize), sizeof(std::size_t));
				for (auto& p : this->continousHistoryEventsEvent)	
				{
					out.write(reinterpret_cast<char*>((std::size_t*)&p.first), sizeof(std::size_t));
					out.write(reinterpret_cast<char*>((std::size_t*)&p.second), sizeof(std::size_t));
				}
			}
		}
		if (ok)
		{
			out.write(reinterpret_cast<char*>(&this->historySliceCumulativeIs), 1);
			if (this->historySliceCumulativeIs)
			{
				auto& sizes = this->historySlicesSize;
				auto& nexts = this->historySlicesSlicesSizeNext;
				auto& prevs = this->historySlicesSliceSetPrev;
				{
					std::size_t hsize = sizes.size();
					out.write(reinterpret_cast<char*>(&hsize), sizeof(std::size_t));
					for (auto& p : sizes)	
					{
						out.write(reinterpret_cast<char*>((std::size_t*)&p.first), sizeof(std::size_t));
						out.write(reinterpret_cast<char*>((std::size_t*)&p.second), sizeof(std::size_t));
					}
				}
				{
					std::size_t hsize = nexts.size();
					out.write(reinterpret_cast<char*>(&hsize), sizeof(std::size_t));
					for (auto& p : nexts)	
					{
						std::size_t i = p.first;
						out.write(reinterpret_cast<char*>(&i), sizeof(std::size_t));
						std::size_t msize = p.second.size();
						out.write(reinterpret_cast<char*>(&msize), sizeof(std::size_t));
						for (auto q : p.second)	
						{
							out.write(reinterpret_cast<char*>((std::size_t*)&q.first), sizeof(std::size_t));
							out.write(reinterpret_cast<char*>((std::size_t*)&q.second), sizeof(std::size_t));
						}
					}
				}
				{
					std::size_t hsize = prevs.size();
					out.write(reinterpret_cast<char*>(&hsize), sizeof(std::size_t));
					for (auto& p : prevs)	
					{
						std::size_t i = p.first;
						out.write(reinterpret_cast<char*>(&i), sizeof(std::size_t));
						std::size_t msize = p.second.size();
						out.write(reinterpret_cast<char*>(&msize), sizeof(std::size_t));
						for (auto q : p.second)	
						{
							out.write(reinterpret_cast<char*>((std::size_t*)&q), sizeof(std::size_t));
						}
					}
				}				
			}
		}
		if (ok)
		{		
			out.write(reinterpret_cast<char*>(&this->frameUnderlyingDynamicIs), 1);
			if (this->frameUnderlyingDynamicIs)
			{
				std::size_t hsize = this->historyFrameUnderlying.size();
				out.write(reinterpret_cast<char*>(&hsize), sizeof(std::size_t));
				for (auto& frameUnderlyingA : this->historyFrameUnderlying)
				{
					std::size_t isize = frameUnderlyingA.size();
					out.write(reinterpret_cast<char*>(&isize), sizeof(std::size_t));
					for (auto v : frameUnderlyingA)	
						out.write(reinterpret_cast<char*>((std::size_t*)&v), sizeof(std::size_t));					
				}
			}	
			out.write(reinterpret_cast<char*>(&this->frameHistoryDynamicIs), 1);
			if (this->frameHistoryDynamicIs)
			{
				std::size_t hsize = this->historyFrameHistory.size();
				out.write(reinterpret_cast<char*>(&hsize), sizeof(std::size_t));
				for (auto& frameHistoryA : this->historyFrameHistory)
				{
					std::size_t isize = frameHistoryA.size();
					out.write(reinterpret_cast<char*>(&isize), sizeof(std::size_t));
					for (auto v : frameHistoryA)	
						out.write(reinterpret_cast<char*>((std::size_t*)&v), sizeof(std::size_t));					
				}
			}
		}
		if (ok)
		{		
			out.write(reinterpret_cast<char*>(&this->underlyingOffsetIs), 1);
			if (this->underlyingOffsetIs)
			{
				std::size_t hsize = this->underlyingsVarsOffset.size();
				out.write(reinterpret_cast<char*>(&hsize), sizeof(std::size_t));
				for (auto mm : this->underlyingsVarsOffset)	
				{
					std::size_t i = mm.first;
					out.write(reinterpret_cast<char*>(&i), sizeof(std::size_t));
					std::size_t msize = mm.second.size();
					out.write(reinterpret_cast<char*>(&msize), sizeof(std::size_t));
					for (auto p : mm.second)	
					{
						out.write(reinterpret_cast<char*>((std::size_t*)&p.first), sizeof(std::size_t));
						out.write(reinterpret_cast<char*>((std::size_t*)&p.second), sizeof(std::size_t));						
					}
				}
			}	
		}
		if (ok)
		{		
			std::size_t hsize = this->induceVarComputeds.size();
			out.write(reinterpret_cast<char*>(&hsize), sizeof(std::size_t));
			if (ok && hsize) 
			{
				for (auto v : this->induceVarComputeds)	
					out.write(reinterpret_cast<char*>((std::size_t*)&v), sizeof(std::size_t));
			}	
		}
		out.close();
		{
		// // trace sizes and transitions
		// if (ok && this->historySliceCachingIs)
		// {
			// auto& discont = this->continousHistoryEventsEvent;
			// EVAL(discont);
			// SizeSizeMap sizes;
			// SizeSizeMap nexts;
			// SizeSizeMap prevs;
			// for (auto& pp : this->historySlicesSize)
				// sizes[pp.first] = pp.second;
			// for (auto& pp : this->historySlicesSlicesSizeNext)
				// for (auto& qq : pp.second)
					// nexts[pp.first] += qq.second;
			// for (auto& pp : this->historySlicesSliceSetPrev)
				// prevs[pp.first] = pp.second.size();
			// EVAL(sizes);
			// std::size_t size = 0;
			// for (auto& pp : sizes)			
				// size += pp.second;
			// EVAL(size);
			// EVAL(nexts);
			// std::size_t next = 0;
			// for (auto& pp : nexts)			
				// next += pp.second;
			// EVAL(next);
			// EVAL(prevs);			
			// std::size_t prev = 0;
			// for (auto& pp : prevs)			
				// prev += pp.second;		
			// EVAL(prev);
		// }			
		}
		if (ok && this->logging)
		{
			LOG "dump\tfile name: " << pp.filename << "\ttime: " << ((sec)(clk::now() - mark)).count() << "s" UNLOG
		}	
	} 
	catch (const std::exception& e) 
	{
		LOG "dump error:\tfailed to dump to file: " << pp.filename << "\terror message: " << e.what()  UNLOG
		ok = false;
		try 
		{
			out.close();			
		}
		catch (const std::exception&) 
		{
		}
	}
	
	return ok;
}

bool Alignment::Active::load(const ActiveIOParameters& pp)
{
	bool ok = true;
	
	std::ifstream in;
	try 
	{
		auto mark = (ok && this->logging) ? clk::now() : std::chrono::time_point<clk>();
		in.exceptions(in.failbit | in.badbit | in.eofbit);
		std::lock_guard<std::mutex> guard(this->mutex);	
		in.open(pp.filename, std::ios::binary);
		if (ok)
		{		
			std::size_t h;
			in.read(reinterpret_cast<char*>(&h), sizeof(std::size_t));
			if (ok && h) 
			{
				std::string s(h,' ');
				in.read(reinterpret_cast<char*>((char*)s.data()), h);
				this->name = s;
			}		
		}
		if (ok)
		{		
			std::size_t h;
			in.read(reinterpret_cast<char*>(&h), sizeof(std::size_t));
			if (ok && h) 
			{
				std::size_t ev;
				in.read(reinterpret_cast<char*>(&ev), sizeof(std::size_t));
				this->underlyingEventUpdateds.clear();
				this->underlyingEventUpdateds.insert(ev);
			}		
		}
		if (ok)
		{		
			in.read(reinterpret_cast<char*>(&this->historySize), sizeof(std::size_t));
			in.read(reinterpret_cast<char*>(&this->historyOverflow), 1);
			in.read(reinterpret_cast<char*>(&this->historyEvent), sizeof(std::size_t));		
		}
		if (ok)
		{		
			std::size_t hsize = 0;
			in.read(reinterpret_cast<char*>(&hsize), sizeof(std::size_t));
			this->underlyingHistoryRepa.clear();
			if (ok && hsize)
			{
				this->underlyingHistoryRepa.reserve(hsize);
				for (std::size_t h = 0; ok && h < hsize; h++)	
				{
					std::unique_ptr<HistoryRepa> hr;
					if (this->historyOverflow)
						hr = persistentsHistoryRepa(in);
					else
						hr = persistentInitialsHistoryRepa(in);
					if (ok)
						this->underlyingHistoryRepa.push_back(std::move(hr));
				}				
			}		
		}
		if (ok)
		{		
			std::size_t hsize = 0;
			in.read(reinterpret_cast<char*>(&hsize), sizeof(std::size_t));
			this->underlyingHistorySparse.clear();
			if (ok && hsize)
			{
				this->underlyingHistorySparse.reserve(hsize);
				for (std::size_t h = 0; ok && h < hsize; h++)	
				{
					std::unique_ptr<HistorySparseArray> hr;
					if (this->historyOverflow)
						hr = persistentsHistorySparseArray(in);
					else
						hr = persistentInitialsHistorySparseArray(in);
					this->underlyingHistorySparse.push_back(std::move(hr));
				}				
			}		
		}
		if (ok)
		{		
			std::size_t hsize = 0;
			in.read(reinterpret_cast<char*>(&hsize), sizeof(std::size_t));
			this->underlyingSlicesParent.clear();
			if (ok && hsize)
			{
				this->underlyingSlicesParent.reserve(hsize);
				for (std::size_t h = 0; ok && h < hsize; h++)	
				{
					std::size_t first;
					in.read(reinterpret_cast<char*>(&first), sizeof(std::size_t));
					std::size_t second;
					in.read(reinterpret_cast<char*>(&second), sizeof(std::size_t));
					this->underlyingSlicesParent.insert_or_assign(first,second);
				}				
			}		
		}		
		if (ok)
		{		
			bool has = false;
			in.read(reinterpret_cast<char*>(&has), 1);
			this->decomp.reset();
			if (ok && has)
			{
				this->decomp = persistentsDecompFudSlicedRepa(in);	
				this->decomp->mapVarInt();
				this->decomp->mapVarParent();
			}
		}
		if (ok)
		{		
			in.read(reinterpret_cast<char*>(&this->bits), sizeof(int));
			in.read(reinterpret_cast<char*>(&this->var), sizeof(std::size_t));
			in.read(reinterpret_cast<char*>(&this->varSlice), sizeof(std::size_t));		
			in.read(reinterpret_cast<char*>(&this->induceThreshold), sizeof(std::size_t));
			std::size_t hsize = 0;
			in.read(reinterpret_cast<char*>(&hsize), sizeof(std::size_t));
			this->induceVarExclusions.clear();		
			for (std::size_t h = 0; ok && h < hsize; h++)	
			{
				std::size_t v;
				in.read(reinterpret_cast<char*>(&v), sizeof(std::size_t));
				this->induceVarExclusions.insert(v);
			}	
		}
		if (ok)
		{		
			bool has = false;
			in.read(reinterpret_cast<char*>(&has), 1);
			this->historySparse.reset();
			this->historySlicesSetEvent.clear();
			this->induceSlices.clear();
			this->induceSliceFailsSize.clear();
			if (ok && has)
			{
				if (this->historyOverflow)
					this->historySparse = persistentsHistorySparseArray(in);
				else
					this->historySparse = persistentInitialsHistorySparseArray(in);		
				auto& hr = this->historySparse;	
				auto& slev = this->historySlicesSetEvent;		
				auto z = this->historyOverflow ? this->historySize : this->historyEvent;
				auto rr = hr->arr;	
				for (std::size_t j = 0; j < z; j++)
					slev[rr[j]].insert(j);
				if (ok && this->induceThreshold)
					for (auto pp : slev)
						if (pp.second.size() >= induceThreshold)
							this->induceSlices.insert(pp.first);					
			}
		}
		if (ok)
		{		
			std::size_t hsize = 0;
			in.read(reinterpret_cast<char*>(&hsize), sizeof(std::size_t));
			this->frameUnderlyings.clear();		
			for (std::size_t h = 0; ok && h < hsize; h++)	
			{
				std::size_t v;
				in.read(reinterpret_cast<char*>(&v), sizeof(std::size_t));
				this->frameUnderlyings.push_back(v);
			}	
		}
		if (ok)
		{		
			std::size_t hsize = 0;
			in.read(reinterpret_cast<char*>(&hsize), sizeof(std::size_t));
			this->frameHistorys.clear();		
			for (std::size_t h = 0; ok && h < hsize; h++)	
			{
				std::size_t v;
				in.read(reinterpret_cast<char*>(&v), sizeof(std::size_t));
				this->frameHistorys.push_back(v);
			}	
		}
		if (ok)
		{		
			std::size_t hsize = 0;
			in.read(reinterpret_cast<char*>(&hsize), sizeof(std::size_t));
			this->framesVarsOffset.clear();		
			for (std::size_t h = 0; ok && h < hsize; h++)	
			{
				std::size_t i;
				in.read(reinterpret_cast<char*>(&i), sizeof(std::size_t));
				auto& mm = this->framesVarsOffset[i];				
				std::size_t msize;
				in.read(reinterpret_cast<char*>(&msize), sizeof(std::size_t));
				mm.reserve(msize);
				for (std::size_t m = 0; ok && m < msize; m++)	
				{
					std::size_t p;
					in.read(reinterpret_cast<char*>(&p), sizeof(std::size_t));
					std::size_t q;
					in.read(reinterpret_cast<char*>(&q), sizeof(std::size_t));
					mm[p] = q;
				}
			}	
		}
		if (ok)
		{		
			this->continousIs = false;
			this->continousHistoryEventsEvent.clear();
			in.exceptions(in.badbit);
			in.read(reinterpret_cast<char*>(&this->continousIs), 1);
			if (!in.eof())
			{
				in.clear();
				in.exceptions(in.failbit | in.badbit | in.eofbit);
				if (this->continousIs)
				{
					std::size_t hsize = 0;
					in.read(reinterpret_cast<char*>(&hsize), sizeof(std::size_t));
					for (std::size_t h = 0; ok && h < hsize; h++)	
					{
						std::size_t first;
						in.read(reinterpret_cast<char*>(&first), sizeof(std::size_t));
						std::size_t second;
						in.read(reinterpret_cast<char*>(&second), sizeof(std::size_t));
						this->continousHistoryEventsEvent.insert_or_assign(first,second);
					}					
				}				
			}
			in.clear();
			in.exceptions(in.failbit | in.badbit | in.eofbit);
		}	
		if (ok)
		{		
			auto& sizes = this->historySlicesSize;
			auto& nexts = this->historySlicesSlicesSizeNext;
			auto& prevs = this->historySlicesSliceSetPrev;
			sizes.clear();
			nexts.clear();
			prevs.clear();
			this->historySliceCumulativeIs = false;
			in.exceptions(in.badbit);
			in.read(reinterpret_cast<char*>(&this->historySliceCumulativeIs), 1);
			if (!in.eof())
			{
				in.clear();
				in.exceptions(in.failbit | in.badbit | in.eofbit);
				if (this->historySliceCumulativeIs)
				{
					this->historySliceCachingIs = true;
					{
						std::size_t hsize = 0;
						in.read(reinterpret_cast<char*>(&hsize), sizeof(std::size_t));
						for (std::size_t h = 0; ok && h < hsize; h++)	
						{
							std::size_t sliceA;
							in.read(reinterpret_cast<char*>(&sliceA), sizeof(std::size_t));
							std::size_t count;
							in.read(reinterpret_cast<char*>(&count), sizeof(std::size_t));
							sizes[sliceA] = count;
						}								
					}
					{
						std::size_t hsize = 0;
						in.read(reinterpret_cast<char*>(&hsize), sizeof(std::size_t));
						for (std::size_t h = 0; ok && h < hsize; h++)	
						{
							std::size_t sliceA = 0;
							in.read(reinterpret_cast<char*>(&sliceA), sizeof(std::size_t));
							std::size_t msize = 0;
							in.read(reinterpret_cast<char*>(&msize), sizeof(std::size_t));
							for (std::size_t m = 0; ok && m < msize; m++)	
							{							
								std::size_t sliceB;
								in.read(reinterpret_cast<char*>(&sliceB), sizeof(std::size_t));
								std::size_t count;
								in.read(reinterpret_cast<char*>(&count), sizeof(std::size_t));
								nexts[sliceA][sliceB] = count;
							}	
						}
					}
					{
						std::size_t hsize = 0;
						in.read(reinterpret_cast<char*>(&hsize), sizeof(std::size_t));
						for (std::size_t h = 0; ok && h < hsize; h++)	
						{
							std::size_t sliceA = 0;
							in.read(reinterpret_cast<char*>(&sliceA), sizeof(std::size_t));
							std::size_t msize = 0;
							in.read(reinterpret_cast<char*>(&msize), sizeof(std::size_t));
							for (std::size_t m = 0; ok && m < msize; m++)	
							{							
								std::size_t sliceB;
								in.read(reinterpret_cast<char*>(&sliceB), sizeof(std::size_t));
								prevs[sliceA].insert(sliceB);
							}	
						}
					}
				}				
			}
			in.clear();
			in.exceptions(in.failbit | in.badbit | in.eofbit);
		}			
		if (ok)
		{		
			this->frameUnderlyingDynamicIs = false;
			this->historyFrameUnderlying.clear();
			this->frameHistoryDynamicIs = false;
			this->historyFrameHistory.clear();
			in.exceptions(in.badbit);
			in.read(reinterpret_cast<char*>(&this->frameUnderlyingDynamicIs), 1);
			if (!in.eof())
			{
				in.clear();
				in.exceptions(in.failbit | in.badbit | in.eofbit);
				if (this->frameUnderlyingDynamicIs)
				{
					this->historyFrameUnderlying.reserve(this->historySize);
					std::size_t hsize = 0;
					in.read(reinterpret_cast<char*>(&hsize), sizeof(std::size_t));
					for (std::size_t h = 0; ok && h < hsize; h++)	
					{
						std::size_t isize = 0;
						in.read(reinterpret_cast<char*>(&isize), sizeof(std::size_t));
						SizeList frameUnderlyingsA;
						frameUnderlyingsA.reserve(isize);
						for (std::size_t i = 0; ok && i < isize; i++)	
						{
							std::size_t v;
							in.read(reinterpret_cast<char*>(&v), sizeof(std::size_t));
							frameUnderlyingsA.push_back(v);
						}	
						this->historyFrameUnderlying.push_back(frameUnderlyingsA);
					}					
				}	
				in.read(reinterpret_cast<char*>(&this->frameHistoryDynamicIs), 1);
				if (this->frameHistoryDynamicIs)
				{
					this->historyFrameHistory.reserve(this->historySize);
					std::size_t hsize = 0;
					in.read(reinterpret_cast<char*>(&hsize), sizeof(std::size_t));
					for (std::size_t h = 0; ok && h < hsize; h++)	
					{
						std::size_t isize = 0;
						in.read(reinterpret_cast<char*>(&isize), sizeof(std::size_t));
						SizeList frameHistorysA;
						frameHistorysA.reserve(isize);
						for (std::size_t i = 0; ok && i < isize; i++)	
						{
							std::size_t v;
							in.read(reinterpret_cast<char*>(&v), sizeof(std::size_t));
							frameHistorysA.push_back(v);
						}	
						this->historyFrameHistory.push_back(frameHistorysA);
					}					
				}					
			}
			in.clear();
			in.exceptions(in.failbit | in.badbit | in.eofbit);
		}	
		if (ok)
		{		
			this->underlyingOffsetIs = false;
			this->underlyingsVarsOffset.clear();
			in.exceptions(in.badbit);
			in.read(reinterpret_cast<char*>(&this->underlyingOffsetIs), 1);
			if (!in.eof())
			{
				in.clear();
				in.exceptions(in.failbit | in.badbit | in.eofbit);
				if (this->underlyingOffsetIs)
				{
					std::size_t hsize = 0;
					in.read(reinterpret_cast<char*>(&hsize), sizeof(std::size_t));	
					for (std::size_t h = 0; ok && h < hsize; h++)	
					{
						std::size_t i;
						in.read(reinterpret_cast<char*>(&i), sizeof(std::size_t));
						auto& mm = this->underlyingsVarsOffset[i];				
						std::size_t msize;
						in.read(reinterpret_cast<char*>(&msize), sizeof(std::size_t));
						mm.reserve(msize);
						for (std::size_t m = 0; ok && m < msize; m++)	
						{
							std::size_t p;
							in.read(reinterpret_cast<char*>(&p), sizeof(std::size_t));
							std::size_t q;
							in.read(reinterpret_cast<char*>(&q), sizeof(std::size_t));
							mm[p] = q;
						}
					}	
				}					
			}
			in.clear();
			in.exceptions(in.failbit | in.badbit | in.eofbit);
		}	
		if (ok)
		{		
			this->induceVarComputeds.clear();		
			in.exceptions(in.badbit);
			std::size_t hsize = 0;
			in.read(reinterpret_cast<char*>(&hsize), sizeof(std::size_t));
			if (!in.eof())
			{
				in.clear();
				in.exceptions(in.failbit | in.badbit | in.eofbit);
				for (std::size_t h = 0; ok && h < hsize; h++)	
				{
					std::size_t v;
					in.read(reinterpret_cast<char*>(&v), sizeof(std::size_t));
					this->induceVarComputeds.insert(v);
				}					
			}
			in.clear();
			in.exceptions(in.failbit | in.badbit | in.eofbit);
		}	
		in.close();
		// cache lengths
		if (ok && historySliceCachingIs && this->decomp && this->historySparse)
		{
			auto& lengths = this->historySlicesLength;
			auto& dr = *this->decomp;
			for (auto& fs : dr.fuds)
			{
				auto sliceA = fs.parent;
				for (auto sliceB : fs.children)
					lengths[sliceB] = lengths[sliceA] + 1;
			}
		}
		// cache sizes and transitions
		if (ok && historySliceCachingIs && !this->historySliceCumulativeIs 
			&& this->decomp && this->historySparse)
		{
			auto over = this->historyOverflow;
			auto cont = this->continousIs;
			auto& discont = this->continousHistoryEventsEvent;
			auto z = this->historySize;
			auto y = this->historyEvent;
			auto rs = this->historySparse->arr;
			auto& slices = this->historySlicesSetEvent;
			auto& sizes = this->historySlicesSize;
			auto& nexts = this->historySlicesSlicesSizeNext;
			auto& prevs = this->historySlicesSliceSetPrev;
			auto& cv = this->decomp->mapVarParent();
			sizes.clear();
			nexts.clear();
			prevs.clear();
			sizes.reserve(slices.size()*3);
			nexts.reserve(slices.size());
			prevs.reserve(slices.size());
			for (auto pp : slices)
			{
				auto sliceC = pp.first;
				auto a = pp.second.size();
				while (true)
				{
					sizes[sliceC] += a;
					if (!sliceC)
						break;
					sliceC = cv[sliceC];
				}								
			}	
			if (cont)
			{
				auto j = over ? y : z;	
				auto sliceB = rs[j%z];
				j++;
				while (j < y+z)
				{
					auto sliceC = rs[j%z];
					if (sliceC != sliceB)
					{
						if (!discont.count(j%z))
						{
							nexts[sliceB][sliceC]++;
							prevs[sliceC].insert(sliceB);
						}
						sliceB = sliceC;
					}
					j++;
				}					
			}
		}	
		{
		// // trace sizes and transitions
		// if (ok && historySliceCachingIs)
		// {
			// auto& discont = this->continousHistoryEventsEvent;
			// EVAL(discont);
			// SizeSizeMap sizes;
			// SizeSizeMap nexts;
			// SizeSizeMap prevs;
			// for (auto& pp : this->historySlicesSize)
				// sizes[pp.first] = pp.second;
			// for (auto& pp : this->historySlicesSlicesSizeNext)
				// for (auto& qq : pp.second)
					// nexts[pp.first] += qq.second;
			// for (auto& pp : this->historySlicesSliceSetPrev)
				// prevs[pp.first] = pp.second.size();
			// EVAL(sizes);
			// std::size_t size = 0;
			// for (auto& pp : sizes)			
				// size += pp.second;
			// EVAL(size);
			// EVAL(nexts);
			// std::size_t next = 0;
			// for (auto& pp : nexts)			
				// next += pp.second;
			// EVAL(next);
			// EVAL(prevs);			
			// std::size_t prev = 0;
			// for (auto& pp : prevs)			
				// prev += pp.second;		
			// EVAL(prev);		
		// }		
		// // frames
		// {
			// EVAL(this->frameUnderlyingDynamicIs);
			// EVAL(this->frameUnderlyings);
			// EVAL(this->historyFrameUnderlying);
			// EVAL(this->frameHistoryDynamicIs);
			// EVAL(this->frameHistorys);
			// EVAL(this->historyFrameHistory);
		// }	
		}
		if (ok && this->logging)
		{
			LOG "load\tfile name: " << pp.filename << "\ttime: " << ((sec)(clk::now() - mark)).count() << "s" UNLOG
		}	
	} 
	catch (const std::exception& e) 
	{
		LOG "load error:\tfailed to load from file: " << pp.filename << "\terror message: " << e.what()  UNLOG
		ok = false;
		try 
		{
			in.close();			
		}
		catch (const std::exception&) 
		{
		}	
	}
	
	return ok;
}

