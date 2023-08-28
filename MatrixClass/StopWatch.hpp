//Gabriele Bernardino - Stopwatch class

#ifndef STOPWATCH_H
#define STOPWATCH_H

#include <chrono>
using namespace std;

class StopWatch
{
public:
    StopWatch() {};
    void Start() {
        start = chrono::system_clock::now();
    }
    void Stop() {
        end = chrono::system_clock::now();
    }
    void Reset() {
        //set both start and end equal to now
        start = chrono::system_clock::now();
        end = start;
    }
    double GetTime() const {
        chrono::duration<double> elapsed = end - start;
        return elapsed.count();
    }
private:
    chrono::time_point<chrono::system_clock> start, end;
    StopWatch(const StopWatch& other) {
        this->start = other.start;
        this->end = other.end;
    }
    StopWatch& operator=(const StopWatch& other) {
        StopWatch watch(other);
        return watch;
    }
};

#endif // !STOPWATCH_H