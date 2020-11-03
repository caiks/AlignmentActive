#ifndef ALIGNMENTACTIVE_H
#define ALIGNMENTACTIVE_H

#include "AlignmentUtil.h"
#include "Alignment.h"
#include "AlignmentRepa.h"

#include <thread>
#include <mutex>

namespace Alignment
{
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
		
		// std::shared_ptr<ApplicationRepa> applicationUnder;
		// std::shared_ptr<ApplicationRepa> application;
		// std::size_t applicationFudIdPersistent;
		// std::size_t applicationFudId;
		

		// std::shared_ptr<SystemRepa> systemUnder;
		// std::shared_ptr<SystemRepa> system;
		
		std::size_t blockBits;
		SizeSet blocks;
		
		// std::shared_ptr<HistoryRepa> history;
		std::vector<std::pair<ActiveEventsRepaPtr,HistoryRepaPtr>> listEventsHistoryRepa;
		std::vector<std::pair<ActiveEventsArrayPtr,HistorySparseArrayPtr>> listEventsHistorySparse;
		bool historyOverflow;
		std::size_t historyEvent;
		std::size_t historySize;
		// std::size_t size() const;
		
		SizeList eventsSlice;
		SizeSizeSetMap slicesSetEvent;
		// SizeSet setSliceLeaf;
		
		// bool report();
		// bool slicesSync(std::size_t tint = 1);
		// bool slicesUpdate(std::size_t tint = 1);
		
		ActiveEventsArrayPtr events;
	};

	
}

#endif
