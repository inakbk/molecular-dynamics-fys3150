#ifndef IO_H
#define IO_H
#include <fstream>
class System;
class StatisticsSampler;
using std::ofstream;

class IO
{
private:
    ofstream file;
    bool first_time = true;
public:
    IO();
    ~IO();

    void saveState(System *system);
    void saveStatistics(System *system, StatisticsSampler *statisticsSampler, double totalExecutionTime);
    void open(const char *filename);
    void close();

};
#endif
