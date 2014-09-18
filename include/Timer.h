/*----------------------------------------------------------------
 *
 * 	Written:       10/10/2013
 *  Last updated:  10/10/2013
 *
 *
 *  This is a custom stopwatch. Use as follows:
 *
 *  Timer myTimer = new Timer("TimerName");
 *  // stuff you want to time
 *  myTimer.stop();
 *
 *----------------------------------------------------------------*/

#pragma once

#include <ctime>
#include <iostream>
#include <memory>
#include <string>
#include <sstream>


class Timer {

private:
	std::clock_t start;
	std::string name;

public:
	Timer(const char* timername) {
		start = std::clock();
		name = timername;
		std::cout << "Starting the process \"" << name << "\"." << std::endl;
	}

	//  resets the timer to the current time
	void reset() {
	    start = std::clock();
	}

	//  records time taken since instantiation and
	//  prints it out to the console
	void stop() {
		std::cout.flush();
		double t = (std::clock() - this->start) / (double) CLOCKS_PER_SEC;
		std::cout << "The process \"" << name << "\" took " << t << " seconds." << std::endl;
	};

    std::string toString() {
        std::ostringstream ss;
        ss << (std::clock() - this->start) / (double) CLOCKS_PER_SEC;
        return ss.str();
    };

};
