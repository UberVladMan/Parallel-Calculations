#ifndef TIMER_H
#define TIMER_H
#pragma once

#include <chrono>

class Timer
{
private:
	std::chrono::steady_clock::duration time;
	std::chrono::steady_clock::time_point clock_begin;
	std::chrono::steady_clock::time_point clock_end;
public:
	Timer() {
		clock_begin = std::chrono::steady_clock::now();
	}
	void start() {
		clock_begin = std::chrono::steady_clock::now();
	}
	void stop() {
		clock_end = std::chrono::steady_clock::now();
		time = clock_end - clock_begin;
	}
	double getTime() {
		return double(time.count()) * std::chrono::steady_clock::period::num / std::chrono::steady_clock::period::den;
	}

	double getCurrentTime() {
		auto current_time = std::chrono::system_clock::now();
		auto duration_in_seconds = std::chrono::duration<double>(current_time.time_since_epoch());

		return duration_in_seconds.count();
	}

};

#endif;