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
		
		std::shared_ptr<SystemRepa> systemUnder;
		std::shared_ptr<ApplicationRepa> applicationUnder;

		std::recursive_mutex mutex;
		
		std::shared_ptr<SystemRepa> system;
		
		std::shared_ptr<HistoryRepa> history;
		bool historyOverflow;
		std::size_t historyEvent;
		
		std::shared_ptr<ApplicationRepa> application;
		std::size_t applicationFudIdPersistent;
		std::size_t applicationFudId;
		
		SizeList eventsSlice;
		SizeSizeSetMap slicesSetEvent;
		
		bool report();
		bool slicesSync(std::size_t tint = 1);
	};

	
}

#endif
