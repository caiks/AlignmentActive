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
		
	struct ActiveEventsRepa
	{
		ActiveEventsRepa(std::size_t referencesA = 0);
		std::mutex mutex;
		std::size_t references;
		std::map<std::size_t,HistoryRepaPtrSizePair> mapIdEvent;
	};
	
	typedef std::shared_ptr<ActiveEventsRepa> ActiveEventsRepaPtr;
		
	typedef std::pair<HistorySparseArrayPtr,std::size_t> HistorySparseArrayPtrSizePair;

	struct ActiveEventsArray
	{
		ActiveEventsArray(std::size_t referencesA = 0);
		std::mutex mutex;
		std::size_t references;
		std::map<std::size_t,HistorySparseArrayPtrSizePair> mapIdEvent;
	};
	
	typedef std::shared_ptr<ActiveEventsArray> ActiveEventsArrayPtr;

	struct ActiveUpdateParameters
	{
		std::size_t mapCapacity = 3;
	};
	
	struct ActiveInduceParameters
	{
		size_t tint = 1;		
		size_t wmax = 18;
		size_t lmax = 8;
		size_t xmax = 128;
		double znnmax = 60000.0 * 2.0 * 100.0 * 100.0;
		size_t omax = 10;
		size_t bmax = 10 * 3;
		size_t mmax = 3;
		size_t umax = 128;
		size_t pmax = 1;
		size_t mult = 1;
		size_t seed = 5;
	};
	
	struct Active
	{
		Active();
		
		bool terminate;
		void (*log)(const std::string&);
		bool logging;
		
		std::mutex mutex;
		
		std::shared_ptr<ActiveSystem> system;

		int bits;
		std::size_t var;
		std::size_t varSlice;
		
		std::vector<ActiveEventsRepaPtr> underlyingEventsRepa;
		HistoryRepaPtrList underlyingHistoryRepa;
		std::vector<ActiveEventsArrayPtr> underlyingEventsSparse;
		HistorySparseArrayPtrList underlyingHistorySparse;
		SizeSet eventsUpdated;
		bool historyOverflow;
		std::size_t historyEvent;
		std::size_t historySize;
		std::map<std::size_t,HistorySparseArrayPtr> slicesPath;
		std::size_t pathLenMax;

		std::shared_ptr<DecompFudSlicedRepa> decomp;
		
		HistorySparseArray historySparse;
		SizeSizeSetMap historySlicesSetEvent;
		
		ActiveEventsArrayPtr eventsSparse;
		
		std::size_t induceThreshold;
		SizeSet induceSlices;
		SizeSet induceVarExlusions;
		SizeSizeMap induceSliceFailsSize;
		
		bool update(ActiveUpdateParameters pp = ActiveUpdateParameters());
		bool (*updateCallback)(const SizeSet& eventsA, std::size_t eventA, std::size_t historyEventA, std::size_t sliceA);

		bool induce(ActiveInduceParameters pp = ActiveInduceParameters(),
					ActiveUpdateParameters ppu = ActiveUpdateParameters());
		
	};
}

std::ostream& operator<<(std::ostream& out, const Alignment::ActiveEventsRepa&);

std::ostream& operator<<(std::ostream& out, const Alignment::ActiveEventsArray&);


#endif
