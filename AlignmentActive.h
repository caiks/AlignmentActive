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
	
	struct ActiveEventsRepa
	{
		ActiveEventsRepa();
		std::mutex mutex;
		std::size_t references;
		std::map<std::size_t,std::pair<HistoryRepaPtr,std::size_t>> mapIdEvent;
	};
	
	typedef std::shared_ptr<ActiveEventsRepa> ActiveEventsRepaPtr;
		
	struct ActiveEventsArray
	{
		ActiveEventsArray();
		std::mutex mutex;
		std::size_t references;
		std::map<std::size_t,std::pair<HistorySparseArrayPtr,std::size_t>> mapIdEvent;
	};
	
	typedef std::shared_ptr<ActiveEventsArray> ActiveEventsArrayPtr;

	struct Active
	{
		Active();
		
		bool terminate;
		void (*log)(const std::string&);
		
		std::recursive_mutex mutex;
		
		// std::shared_ptr<SystemRepa> systemUnder;
		// std::shared_ptr<SystemRepa> system;	
		std::shared_ptr<ActiveSystem> system;

		int blockBits;
		std::size_t blockCurrent;
		
		// std::shared_ptr<HistoryRepa> history;
		std::vector<ActiveEventsRepaPtr> underlyingEventsRepa;
		std::vector<HistoryRepaPtr> underlyingHistoryRepa;
		std::vector<ActiveEventsArrayPtr> underlyingEventsSparse;
		std::vector<HistorySparseArrayPtr> underlyingHistorySparse;
		SizeSet eventsUpated;
		bool historyOverflow;
		std::size_t historyEvent;
		std::size_t historySize;
		std::map<std::size_t,HistorySparseArrayPtr> slicesPath;
		// std::size_t size() const;
		
		// std::shared_ptr<ApplicationRepa> applicationUnder;
		// std::shared_ptr<ApplicationRepa> application;
		// std::size_t applicationFudIdPersistent;
		// std::size_t applicationFudId;
		std::shared_ptr<DecompFudSlicedRepa> decomp;
		
		SizeList eventsSlice;
		SizeSizeSetMap slicesSetEvent;
		// SizeSet setSliceLeaf;
		
		ActiveEventsArrayPtr eventsSparse;
		
		// bool report();
		// bool slicesSync(std::size_t tint = 1);
		// bool slicesUpdate(std::size_t tint = 1);
		bool update();
		
	};

	
}

#endif
