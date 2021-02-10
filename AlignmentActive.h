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
		
		bool terminate;
		void (*log)(Active& active, const std::string&);
		void (*layerer_log)(const std::string&);
		bool logging;
		bool summary;
		
		void* client;
		
		std::mutex mutex;
		
		std::vector<ActiveEventsRepaPtr> underlyingEventsRepa;
		std::vector<ActiveEventsArrayPtr> underlyingEventsSparse;
		SizeSet underlyingEventUpdateds;
		
		std::size_t historySize;
		bool historyOverflow;
		std::size_t historyEvent;
		
		HistoryRepaPtrList underlyingHistoryRepa;
		HistorySparseArrayPtrList underlyingHistorySparse;
		SizeSizeUMap underlyingSlicesParent;

		std::unique_ptr<DecompFudSlicedRepa> decomp;
		
		std::unique_ptr<HistorySparseArray> historySparse;
		SizeSizeSetMap historySlicesSetEvent;
		
		ActiveEventsArrayPtr eventsSparse;
		
		std::shared_ptr<ActiveSystem> system;
		
		int bits;
		std::size_t var;
		std::size_t varSlice;
		
		std::size_t induceThreshold;
		SizeSet induceSlices;
		SizeSet induceVarExclusions;
		SizeSizeMap induceSliceFailsSize;
		
		SizeSet frameUnderlyings;
		SizeSet frameHistorys;
		std::map<std::size_t, SizeSizeUMap> framesVarsOffset;		
		
		std::size_t varMax() const;
		
		bool update(ActiveUpdateParameters pp = ActiveUpdateParameters());
		bool (*updateCallback)(Active& active, const SizeSet& eventsA, std::size_t eventA, std::size_t historyEventA, std::size_t sliceA);

		bool induce(ActiveInduceParameters pp = ActiveInduceParameters(),
					ActiveUpdateParameters ppu = ActiveUpdateParameters());
		bool (*induceCallback)(Active& active, std::size_t sliceA, std::size_t sliceSizeA);	

		bool dump(const ActiveIOParameters&);
		bool load(const ActiveIOParameters&);	
	};
}

std::ostream& operator<<(std::ostream& out, const Alignment::ActiveEventsRepa&);

std::ostream& operator<<(std::ostream& out, const Alignment::ActiveEventsArray&);

#endif
