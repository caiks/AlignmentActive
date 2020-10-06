#ifndef ALIGNMENTACTIVE_H
#define ALIGNMENTACTIVE_H

#include "AlignmentUtil.h"
#include "Alignment.h"
#include "AlignmentRepa.h"

#include <thread>
#include <mutex>

namespace Alignment
{
	struct Active
	{
		Active();
		
		bool terminate;
		void (*log)(const std::string&);
		
		std::recursive_mutex mutex;
		
		std::shared_ptr<SystemRepa> systemUnder;
		std::shared_ptr<SystemRepa> system;
		
		std::shared_ptr<HistoryRepa> history;
		bool historyOverflow;
		std::size_t historyEvent;
		std::size_t size() const;
		
		std::shared_ptr<ApplicationRepa> applicationUnder;
		std::shared_ptr<ApplicationRepa> application;
		std::size_t applicationFudIdPersistent;
		std::size_t applicationFudId;
		
		SizeList eventsSlice;
		SizeSizeSetMap slicesSetEvent;
		SizeSet setSliceLeaf;
		
		bool report();
		bool slicesSync(std::size_t tint = 1);
		bool slicesUpdate(std::size_t tint = 1);
	};

	
}

#endif
