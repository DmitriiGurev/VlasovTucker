#pragma once

#include <iostream>
#include <chrono>
#include <string>

class Timer
{
    typedef std::chrono::high_resolution_clock Clock;
    typedef std::chrono::high_resolution_clock::time_point TimePoint;
    typedef std::chrono::microseconds mus;
public:
    Timer()
    {
        _startTime = Clock::now();
        _prevTime = _startTime;
    }

    ~Timer()
    {
        TimePoint stopTime = Clock::now();
        auto duration = std::chrono::duration_cast<mus>(stopTime - _startTime);

        std::cerr << "Total execution time: " <<  duration.count() / pow(10, 6) << " s\n";
    }

    void StartSection()
    {
        _prevTime = Clock::now();
    }

    void PrintSectionTime(std::string sectionName)
    {
        auto duration = std::chrono::duration_cast<mus>(Clock::now() - _prevTime);
        _prevTime = Clock::now();

        std::cerr << sectionName << ": " <<  duration.count() / pow(10, 6) << " s\n";
    }

private:
    TimePoint _startTime;
    TimePoint _prevTime;
};