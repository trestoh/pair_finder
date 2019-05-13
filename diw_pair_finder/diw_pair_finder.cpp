// full_run_subarray_sum.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <chrono>
#include <ctime>
#include <map>
#include <vector>
#include <ratio>
#include <windows.h>
#include <stdlib.h>
#include <fcntl.h>
#include <queue>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <limits>
#include <functional>
#include <numeric>
#include <tuple>

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>
#include <OpenMS/CHEMISTRY/ElementDB.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/KERNEL/Peak1D.h>

#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>


#include <boost/tokenizer.hpp>

enum ScoreType { ALL_ADDED, RATIO, TARGET_INT };

class WindowSpec {

public:
	WindowSpec(int s, std::string w, std::string c) : sister(s), width(w), center(c) {};

	int sister;
	std::string width;
	std::string center;




	friend bool operator>(const WindowSpec &w1, const WindowSpec &w2);
	friend bool operator<(const WindowSpec &w1, const WindowSpec &w2);
	friend bool operator==(const WindowSpec &w1, const WindowSpec &w2);


};

class Offset {

public:
	Offset(int s, double o) : sister(s), offset(o) {};

	int sister;
	double offset;

};

class peak
{
public:
	double mz;
	double intensity;

	peak(double m = 0.0, double i = 0.0) : mz(m), intensity(i) {};
};

bool operator> (const WindowSpec &w1, const WindowSpec &w2)
{
	if (w1.center > w2.center)
		return true;
	else if (w1.center == w2.center)
		return w1.width > w2.width;
	else
		return false;
}

bool operator< (const WindowSpec &w1, const WindowSpec &w2)
{
	if (w1.center < w2.center)
		return true;
	else if (w1.center == w2.center)
		return w1.width < w2.width;
	else
		return false;
}

bool operator== (const WindowSpec &w1, const WindowSpec &w2)
{
	//bool center_close = (abs(atof(w1.center.c_str()) - atof(w2.center.c_str())) <= 0.01);
	//bool off_close = (abs(atof(w1.width.c_str()) - atof(w2.width.c_str())) <= 0.01);
	return w1.center == w2.center && w1.width == w2.width;
}

class targetScore
{
public:
	double target;
	double score;

	targetScore(double t = 0.0, double s = 0.0) : target(t), score(s) {};
};

class family
{
public:
	int parent;
	int sister;

	family(int p = 0, int s = 0) : parent(p), sister(s) {};
	friend bool operator<(const family &f1, const family &f2);
};

bool operator< (const family &f1, const family &f2)
{
	if (f1.parent < f2.parent)
		return true;
	else if (f1.parent == f2.parent)
		return f1.sister < f2.sister;
	else
		return false;
}

int main(int argc, char * argv[])
{

	std::ifstream dyn_file(argv[1]);

	std::ifstream window_file(argv[2]);

	std::ofstream log("log.txt");

	std::map<int, bool> dyn_scans;

	std::map<int, double> inj_times;

	std::map<WindowSpec, std::vector<Offset>> offsets;

	if (dyn_file)
	{
		std::string line;
		int scan;

		while (std::getline(dyn_file, line))
		{
			std::istringstream iss(line);
			iss >> scan;

			dyn_scans.insert({ scan, true });

		}

	}

	if (window_file)
	{
		std::string line;
		int scan;
		std::string width;
		std::string center;
		double offset;

		while (std::getline(window_file, line))
		{
			std::istringstream iss(line);
			iss >> scan >> center >> offset >> width;

			center = center.substr(0, center.find('.') + 3);
			width = width.substr(0, width.find('.') + 2);

			log << "Adding " << scan << " " << width << " " << center << std::endl;


			auto offset_finder = offsets.find(WindowSpec(scan, width, center));
			if (offset_finder == offsets.end())
			{
				std::vector<Offset> temp;
				temp.push_back(Offset(scan, offset));
				offsets.insert({ WindowSpec(scan, width, center), temp });
			}

			else
			{
				offset_finder->second.push_back(Offset(scan, offset));
			}

		}

	}

	std::map<WindowSpec, std::string> dyn_parents;

	std::map<std::string, std::vector<family>> offset_family_lists;

	std::ofstream pair_f;
	pair_f.open("pairs.txt");

	OpenMS::MzMLFile mzMLDataFileProfile;
	OpenMS::MSExperiment msExperimentProfile;

	try {
		mzMLDataFileProfile.load(argv[3], msExperimentProfile);
	}
	catch (std::exception& e) {
		std::cout << e.what() << std::endl;
		return 1;
	}

	OpenMS::MSSpectrum<> s;
	int ms1_count = 0;
	int reg_ms2_count = 0;
	int dyn_ms2_count = 0;
	int found_dynamic = 0;
	int missing_dynamic = 0;

	int list_scan = 0;
	int diff = 0;

	std::string scan;

	int previous_ms1 = 1;

	for (int i = 0; i < msExperimentProfile.getNrSpectra(); i++)
	{
		s = msExperimentProfile.getSpectrum(i);

		std::string junk;
		std::istringstream natID = std::istringstream(s.getNativeID());
		natID >> junk >> junk >> scan;
		scan = scan.substr(scan.find('=') + 1);

		int scan_num = atoi(scan.c_str());

		if (s.getMSLevel() == 1)
		{
			std::string junk;
			std::istringstream natID = std::istringstream(s.getNativeID());
			natID >> junk >> junk >> scan;
			scan = scan.substr(scan.find('=') + 1);

			int scan_num = atoi(scan.c_str());

			previous_ms1 = scan_num;
		}

		else
		{
			
			auto scan_finder = dyn_scans.find(scan_num);
			if (scan_finder != dyn_scans.end())
			{
				dyn_ms2_count++;

				auto huh = s.getPrecursors().at(0).getPos();

				double window_width = 2 * s.getPrecursors().at(0).getIsolationWindowUpperOffset();
				std::string width_s = std::to_string(window_width);
				width_s = width_s.substr(0, width_s.find('.') + 2);

				std::string center_s = std::to_string(huh);
				center_s = center_s.substr(0, center_s.find('.') + 3);

				bool look_lower = false;
				bool look_higher = false;
				bool found = false;

				double offset;
				int sister;

				auto window_finder = offsets.find(WindowSpec(scan_num, width_s, center_s));
				if (window_finder != offsets.end())
				{
					int distance = 9999;
					for (int i = 0; i < window_finder->second.size() && window_finder->second.at(i).sister < scan_num; i++)
					{
						distance = scan_num - window_finder->second.at(i).sister;
						offset = window_finder->second.at(i).offset;

						sister = window_finder->second.at(i).sister;
					}
					if (distance > 50)
						look_lower = true;

					else
						found = true;
				}

				else
				{
					look_lower = true;
				}

				if (look_lower)
				{
					double temp = huh - .01;
					center_s = std::to_string(temp);
					center_s = center_s.substr(0, center_s.find('.') + 3);

					auto window_finder = offsets.find(WindowSpec(scan_num, width_s, center_s));
					if (window_finder != offsets.end())
					{
						std::cout << "Things as expected" << std::endl;
						int distance = 9999;
						for (int i = 0; i < window_finder->second.size() && window_finder->second.at(i).sister < scan_num; i++)
						{
							distance = scan_num - window_finder->second.at(i).sister;
							offset = window_finder->second.at(i).offset;

							sister = window_finder->second.at(i).sister;

						}
						if (distance > 50)
							look_higher = true;

						else
							found = true;
					}

					else
					{
						look_higher = true;
					}

				}

				if (look_higher)
				{
					double temp = huh + .01;
					center_s = std::to_string(temp);
					center_s = center_s.substr(0, center_s.find('.') + 3);

					auto window_finder = offsets.find(WindowSpec(scan_num, width_s, center_s));
					if (window_finder != offsets.end())
					{
						std::cout << "Things as expected" << std::endl;
						int distance = 9999;
						for (int i = 0; i < window_finder->second.size() && window_finder->second.at(i).sister < scan_num; i++)
						{
							distance = scan_num - window_finder->second.at(i).sister;
							offset = window_finder->second.at(i).offset;

							sister = window_finder->second.at(i).sister;
						}
						if (distance > 50)
							look_higher = true;

						else
							found = true;
					}

					else
					{
						look_higher = true;
					}
				}

				if (!found)
				{
					log << "Issues with scan: " << scan_num << " " << width_s << " " << center_s << std::endl;
				}

				pair_f << scan_num << "\t" << sister << std::endl;

			}

		}

	}

	//close files and report numbers
	pair_f.close();

	return 0;
}

