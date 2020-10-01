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
		
		Active active;
		activesLog(active,log);
	}
	
	if (true)
	{
		auto log = [](const std::string& str)
		{
			std::cout << str << std::endl;
			return;
		};
		
		auto reporter = [](Active& act, void (*log)(const std::string&), std::size_t sleep)
		{
			while (!act.terminate)
			{
				activesLog(act,log);				
				std::this_thread::sleep_for(std::chrono::milliseconds(sleep));
			}
			return;
		};
		
		Active active;
		cout << "starting" << endl;
		std::thread t1(reporter, std::ref(active), log, 2000);
		std::this_thread::sleep_for(std::chrono::milliseconds(9000));
		active.terminate = true;
		t1.join();
		cout << "finished" << endl;
	}
	
	return 0;
}
