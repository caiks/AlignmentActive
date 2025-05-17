#ifndef ALIGNMENTACTIVE_H
#define ALIGNMENTACTIVE_H

#include "AlignmentUtil.h"
#include "Alignment.h"
#include "AlignmentRepa.h"

#include <thread>
#include <mutex>

namespace Alignment
{
	struct ActiveSystem
	{
		ActiveSystem();
		std::mutex mutex;
		std::shared_ptr<SystemRepa> system;
		int bits;
		std::size_t block;
		std::size_t next(int bitsA);
	};
	
	typedef std::pair<HistoryRepaPtr,std::size_t> HistoryRepaPtrSizePair;
		
	struct ActiveEventRepa
	{
		std::size_t id = 0;
		HistoryRepaPtr state;
	};
	
	typedef std::shared_ptr<ActiveEventRepa> ActiveEventRepaPtr;
		
	struct ActiveEventSparse
	{
		std::size_t id = 0;
		HistorySparseArrayPtr state;
	};
	
	typedef std::shared_ptr<ActiveEventSparse> ActiveEventSparsePtr;
		
	struct ActiveUpdateParameters
	{
		std::size_t mapCapacity = 3;
	};
	
	struct ActiveInduceParameters
	{
		std::size_t tint = 1;		
		std::size_t wmax = 18;
		std::size_t lmax = 8;
		std::size_t xmax = 128;
		double znnmax = 60000.0 * 2.0 * 100.0 * 100.0;
		std::size_t omax = 10;
		std::size_t bmax = 10 * 3;
		std::size_t mmax = 3;
		std::size_t umax = 128;
		std::size_t pmax = 1;
		std::size_t mult = 1;
		std::size_t seed = 5;
		double diagonalMin = 0.0;
		std::set<std::size_t> induceThresholds;
		inline bool induceThresholdExceeded(std::size_t prev, std::size_t curr) const
		{
			bool exceeded = false;
			for (auto t : induceThresholds)
				if (t <= curr)
				{
					if (prev < t)
					{
						exceeded = true;
						break;
					}
				}
				else
					break;
			return exceeded;
		}
		std::size_t asyncThreadMax = 0;		
		std::size_t asyncInterval = 10;
		std::size_t asyncUpdateLimit = 0;
		bool logging = false;
	};
	
	struct ActiveIOParameters
	{
		std::string filename;
	};
	
	struct Active
	{
		Active(std::string nameA = "");
		
		std::string name;
		
		volatile bool terminate;
		void (*log)(Active& active, const std::string&);
		void (*layerer_log)(const std::string&);
		bool logging;
		bool summary;
		
		void* client;
		
		std::mutex mutex;
		
		std::vector<ActiveEventRepaPtr> underlyingEventsRepa;
		std::vector<ActiveEventSparsePtr> underlyingEventsSparse;
		std::size_t underlyingEventUpdated;
		
		std::size_t historySize;
		bool historyOverflow;
		std::size_t historyEvent;
		
		bool continousIs;
		SizeSizeMap continousHistoryEventsEvent;
		
		HistoryRepaPtrList underlyingHistoryRepa;
		HistorySparseArrayPtrList underlyingHistorySparse;
		SizeSizeUMap underlyingSlicesParent;

		std::shared_ptr<DecompFudSlicedRepa> decomp;
		
		std::unique_ptr<HistorySparseArray> historySparse;
		SizeSizeSetMap historySlicesSetEvent;
		
		bool historySliceCachingIs;
		bool historySliceCumulativeIs;
		SizeSizeUMap historySlicesSize;
		SizeSizeUMap historySlicesLength;
		std::unordered_map<std::size_t, SizeSizeMap> historySlicesSlicesSizeNext;
		std::unordered_map<std::size_t, SizeSet> historySlicesSliceSetPrev;
			
		ActiveEventSparsePtr eventSparse;
		
		std::shared_ptr<ActiveSystem> system;
		
		int bits;
		std::size_t var;
		std::size_t varSlice;
		
		std::size_t induceThreshold;
		SizeSet induceSlices;
		SizeSet induceVarExclusions;
		SizeSet induceVarComputeds;
		SizeSizeMap induceSliceFailsSize;
		SizeSet inducingSlices;
		volatile bool updateProhibit;
		
		// if dynamic pad out empty frames with 0 so that frame vector length is constant
		// if underlying current frame 0 is present in any it must be the first of all underlying frames
		bool frameUnderlyingDynamicIs;
		SizeList frameUnderlyings;
		SizeListList historyFrameUnderlying;
		bool frameHistoryDynamicIs;
		SizeList frameHistorys;
		SizeListList historyFrameHistory;
		std::map<std::size_t, SizeSizeUMap> framesVarsOffset;	
		
		// if underlyingOffsetIs then all and only sparse underlying are promoted otherwise none are
		bool underlyingOffsetIs;
		std::map<std::size_t, SizeSizeUMap> underlyingsVarsOffset;	

		void varPromote(SizeSizeUMap&, std::size_t&);		
		std::size_t varDemote(const SizeSizeUMap&, std::size_t) const;		
		
		std::size_t varMax() const;
		std::size_t varComputedMax() const;
		
		bool update(ActiveUpdateParameters pp = ActiveUpdateParameters());
		bool (*updateCallback)(Active& active, std::size_t eventA, std::size_t historyEventA, std::size_t sliceA);

		bool induce(ActiveInduceParameters pp = ActiveInduceParameters(),
					ActiveUpdateParameters ppu = ActiveUpdateParameters());
		bool induce(std::size_t sliceA, ActiveInduceParameters pp = ActiveInduceParameters(),
					ActiveUpdateParameters ppu = ActiveUpdateParameters());
		bool (*induceCallback)(Active& active, std::size_t sliceA, std::size_t sliceSizeA);	

		bool dump(const ActiveIOParameters&);
		bool load(const ActiveIOParameters&);	
	};
}

std::ostream& operator<<(std::ostream& out, const Alignment::ActiveEventRepa&);

std::ostream& operator<<(std::ostream& out, const Alignment::ActiveEventSparse&);

#endif
