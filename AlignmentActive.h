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
		std::shared_ptr<SystemRepa> system;
		int blockBits;
		std::size_t blockCurrent;
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

	struct Active
	{
		Active();
		
		bool terminate;
		void (*log)(const std::string&);
		
		std::mutex mutex;
		
		std::shared_ptr<ActiveSystem> system;

		int blockBits;
		std::size_t blockCurrent;
		
		std::vector<ActiveEventsRepaPtr> underlyingEventsRepa;
		HistoryRepaPtrList underlyingHistoryRepa;
		std::vector<ActiveEventsArrayPtr> underlyingEventsSparse;
		HistorySparseArrayPtrList underlyingHistorySparse;
		SizeSet eventsUpdated;
		bool historyOverflow;
		std::size_t historyEvent;
		std::size_t historySize;
		std::map<std::size_t,HistorySparseArrayPtr> slicesPath;

		std::shared_ptr<DecompFudSlicedRepa> decomp;
		
		SizeList eventsSlice;
		SizeSizeSetMap slicesSetEvent;
		
		ActiveEventsArrayPtr eventsSparse;
		
		std::size_t induceThreshold;
		SizeSet slicesInduce;
		
		bool update(std::size_t mapCapacity = 5);
		
	};
}

std::ostream& operator<<(std::ostream& out, const Alignment::ActiveEventsRepa&);

std::ostream& operator<<(std::ostream& out, const Alignment::ActiveEventsArray&);


#endif
