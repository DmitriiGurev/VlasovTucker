#pragma once

#include <string>
#include <iostream>
#include <fstream>
#include <cassert>

namespace VlasovTucker
{
enum class LogLevel
{
    None,
    Console,
    TextFile,
    AllText
};

// TODO: Remove definitions from the class 
class Log
{
public:
    Log() : _level(LogLevel::None), _fPrefix("") {}

    Log(LogLevel level, std::string prefix = "") : 
        _level(level), _fPrefix(prefix)
    {
        if (_PrintToFile())
            _fOut.open(_fPrefix + "log.txt");
    }

    Log& operator=(const Log& other)
    {
        _fOut.close();
        
        _level = other._level;
        _fPrefix = other._fPrefix;
        
        if (_PrintToFile())
            _fOut.open(_fPrefix + "log.txt");

        return *this;
    }

    ~Log()
    {
        if (_PrintToFile())
            _fOut.close();
    }

    template<typename T>
    Log& operator<<(const T& data)
    {
        if (_PrintToConsole())
            std::cout << data;

        if (_PrintToFile())
        {
            assert(_fOut.is_open());
            _fOut << data << std::flush;
        }
        
        return *this;
    }

private:
    bool _PrintToFile() const
    {
        return _level == LogLevel::AllText || _level == LogLevel::TextFile;
    }

    bool _PrintToConsole() const
    {
        return _level == LogLevel::AllText || _level == LogLevel::Console;
    }

private:
    LogLevel _level;

    std::ofstream _fOut; 
    std::string _fPrefix;
};

static std::string Indent(int level) {
    int basicIndent = 4;
    return std::string(basicIndent * level, ' ');
}
}