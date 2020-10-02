#include "AlignmentUtil.h"
#include "Alignment.h"
#include "AlignmentApprox.h"
#include "AlignmentAeson.h"
#include "AlignmentRepa.h"
#include "AlignmentAesonRepa.h"
#include "AlignmentRandomRepa.h"
#include "AlignmentPracticableRepa.h"
#include "AlignmentPracticableIORepa.h"
#include "AlignmentActive.h"
#include <iomanip>
#include <set>
#include <unordered_set>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdio>
#include <thread>
#include <chrono>
#include <ctime>

using namespace Alignment;
using namespace std;

#define ECHO(x) cout << #x << endl; x
#define EVAL(x) cout << #x << ": " << (x) << endl
#define EVALL(x) cout << #x << ": " << endl << (x) << endl
#define TRUTH(x) cout << #x << ": " << ((x) ? "true" : "false") << endl
	
int main(int argc, char **argv)
{
	if (false)
	{
		auto log = [](const std::string& str)
		{
			std::cout << str << std::endl;
			return;
		};
		
		ECHO(Active active);
		TRUTH(active.log());
	}
	
	if (true)
	{
		auto log = [](const std::string& str)
		{
			std::cout << ">>> " << str << std::endl;
			return;
		};
		
		ECHO(Active active);
		active.log_f = log;
		TRUTH(active.log());
	}
	
	if (false)
	{
		auto reporter = [](Active& act, std::size_t sleep)
		{
			while (!act.terminate)
			{
				TRUTH(act.log());			
				std::this_thread::sleep_for(std::chrono::milliseconds(sleep));
			}
			return;
		};
		
		ECHO(Active active);
		std::thread t1(reporter, std::ref(active), 2000);
		std::this_thread::sleep_for(std::chrono::milliseconds(9000));
		active.terminate = true;
		ECHO(t1.join());
	}
	
	return 0;
}
