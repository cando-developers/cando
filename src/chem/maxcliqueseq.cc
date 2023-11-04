/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 */

#include <iostream>
#include <bitset>
#include <cassert>
#include <cando/chem/maxcliqueseq.h>
/*
#define TRACE_LEVEL 1
#define TRACE_MASK TRACE_MASK_CLIQUE
*/
#define TRACE_LEVEL (0*10)
#define TRACE_MASK (0*TRACE_MASK_CLIQUE)



// #include "MaximumCliqueBase.h"
/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 */

#ifndef MAXIMUMCLIQUEBASE_H
#define MAXIMUMCLIQUEBASE_H


// #include "KillTimer.h"


/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 */

#ifndef KILLTIMER_H_INCLUDED
#define KILLTIMER_H_INCLUDED


// the following are only required for limiting execution time
#include <thread>
#include <condition_variable>
#include <chrono>


/** \brief Timer used to limit the execution time of a part of the code
    The kill happens passively, the code that is a subject to kill should periodically check
    variable @timedOut. Alternatively the onKill function could be overloaded with the desired
    way of killing
 */
class KillTimer {
    bool stop;          // signal for the timer to stop (or cancel its kill)
    double seconds;     // number of seconds to wait before the kill
    std::chrono::system_clock::time_point clockZero; // should use steady_clock instead of the monotonic_clock (when properly implemented in GCC)
    std::mutex m;       // closely bonded to cond
    std::condition_variable cond;   // cond.var. used to wait for the timout
    std::unique_ptr<std::thread> thread;  // thread that goes into the wait state for the designated time
    
public:
    /*
        todo: 
            time is for now in microsecond resolution -> should consider nano
            kill is too soft
    */
    bool timedOut;      // a function to be "killed" should periodically check this variable, upon true the function should stop
    
    KillTimer() : stop(false), seconds(0.0), timedOut(false) {}
    ~KillTimer() {
        if (thread)
            thread->join();
    }
    
    void start(double sec) {
        seconds = sec;
        timedOut = false;
        stop = false;
        clockZero = std::chrono::system_clock::now();
        thread = std::unique_ptr<std::thread>(new std::thread([this](){this->threadFunc();}));
    }
    
    void cancel() {
        std::unique_lock<std::mutex> lock(m);
        stop = true;
        cond.notify_all();
    }
    
    void forceTimeout() {
        std::unique_lock<std::mutex> lock(m);
        stop = true;
        timedOut = true;
        cond.notify_all();
        // immitate ordinary (time-induced) timeout 
    }
    
    // overload this function to change the behaviour of kill timer
    void onKill() {}
    
protected:
    void threadFunc() {
        //std::chrono::milliseconds dura((int)(seconds*1000));
        //std::this_thread::sleep_for(std::chrono::milliseconds(dura));
        std::unique_lock<std::mutex> lock(m);
        
        // to filter sporious awakenings, a loop is used; to exit the loop either stop must be true, canceling the kill, or the timeout must have passed
        using namespace std::chrono;
        // the following either waits until the set timeout or until stop is true and the condition variable has been notified
        bool condTimedout = !cond.wait_until(lock, clockZero + std::chrono::microseconds((long long int)(seconds*1000000)), [this](){return stop;});
        timedOut |= condTimedout;
        if (timedOut) onKill();
        stop = false; // reset stop
    }
};


// for windows users that have a function called KillTimer
typedef KillTimer KillTimer1;


#endif // KILLTIMER_H_INCLUDED




#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <string>
#include <utility>
#include <bitset>

// #include "Timer.h"

/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 */

#ifndef TIMER_H_INCLUDED
#define TIMER_H_INCLUDED


/**
    Timers using nanosecond resolution linux clock_gettime
    
    class PrecisionTimer
    class ScopeTimer
    class ExcludeScopeFromTimer
    
**/


#include <vector>
#include <sstream>

#define __linux
#ifndef __linux
	#include <windows.h>
	typedef LARGE_INTEGER time_type ;
	
	#define addDifference(dest, a, b) dest.QuadPart += a.QuadPart - b.QuadPart;
#else
	#include <sys/time.h>
	#include <time.h>
	typedef timespec time_type;

	#define addDifference(dest, a, b) dest = timespec_subadd(a, b, dest);
#endif //WIN32


class TimeType {
    time_type time;
    
public:
    TimeType() {}
    TimeType(const time_type& t) : time(t) {}
    
	double toSeconds() const {
#if defined WIN32 || WIN64
		LARGE_INTEGER freq; 
		QueryPerformanceFrequency(&freq);
		return (double)time.QuadPart / (double)freq.QuadPart;
#else
		return (double)time.tv_sec + time.tv_nsec*0.000000001;
#endif //WIN32
	}
	
#ifdef __linux
	void getSysTime() {clock_gettime(CLOCK_REALTIME, &time);} // when linking, must include: -lrt
	void init() {time.tv_sec = time.tv_nsec = 0l;}
#else
	void getSysTime() {QueryPerformanceCounter(&time);}
	void init() {time.QuadPart = 0;}
#endif //linux
	
	TimeType operator- (const TimeType& t) const {return TimeType(time_sub(time, t.time));}
	TimeType operator+ (const TimeType& t) const {return TimeType(time_add(time, t.time));}
	const TimeType& operator+= (const TimeType& t) {time=time_add(time, t.time); return *this;}
	const TimeType& operator-= (const TimeType& t) {time=time_sub(time, t.time); return *this;}
	
protected:
#if defined WIN32 || WIN64
    static LARGE_INTEGER time_sub(const LARGE_INTEGER& x, const LARGE_INTEGER& y) {
        LARGE_INTEGER r = x;
        r.QuadPart-=y.QuadPart;
        return r;
    }
    
    static LARGE_INTEGER time_add(const LARGE_INTEGER& x, const LARGE_INTEGER& y) {
        LARGE_INTEGER r=x;
        r.QuadPart+=y.QuadPart;
        return r;
    }
#else
    /// linux function for subtracting microsecond time values (that are composed of integer seconds and second fractions part)
    /// result = x - y
    static timeval time_sub(const timeval& x, const timeval& y) {
        timeval result = y;
        // Perform the carry for the later subtraction by updating y.
        if (x.tv_usec < y.tv_usec) {
            int nsec = (y.tv_usec - x.tv_usec) / 1000000 + 1;
            result.tv_usec -= 1000000 * nsec;
            result.tv_sec += nsec;
        }
        if (x.tv_usec - result.tv_usec > 1000000) {
            int nsec = (x.tv_usec - result.tv_usec) / 1000000;
            result.tv_usec += 1000000 * nsec;
            result.tv_sec -= nsec;
        }

        // tv_usec is now certainly positive.
        result.tv_sec = x.tv_sec - result.tv_sec;
        result.tv_usec = x.tv_usec - result.tv_usec;
        return result;
    }

    /// result = (x - y) + z
    timeval time_subadd(const timeval& x, const timeval& y, const timeval& z) {
        timeval result = y;
        // Perform the carry for the later subtraction by updating y.
        if (x.tv_usec < y.tv_usec) {
            int nsec = (y.tv_usec - x.tv_usec) / 1000000 + 1;
            result.tv_usec -= 1000000 * nsec;
            result.tv_sec += nsec;
        }
        if (x.tv_usec - result.tv_usec > 1000000) {
            int nsec = (x.tv_usec - result.tv_usec) / 1000000;
            result.tv_usec += 1000000 * nsec;
            result.tv_sec -= nsec;
        }

        // tv_usec is now certainly positive.
        result.tv_sec = z.tv_sec + x.tv_sec - result.tv_sec;
        result.tv_usec = z.tv_usec + x.tv_usec - result.tv_usec;
        if (result.tv_usec > 1000000) {
            result.tv_usec -= 1000000;
            result.tv_sec++;
        }
            
        return result;
    }

    /// result = x + y
    static timeval time_add(const timeval& x, const timeval& y) {
        timeval result;
        result.tv_sec = x.tv_sec + y.tv_sec;
        result.tv_usec = x.tv_usec + y.tv_usec;
        while (result.tv_usec > 1000000) {
            result.tv_usec -= 1000000;
            result.tv_sec++;
        }
            
        return result;
    }

#endif
    /// linux function for subtracting nanosecond timer values (that are composed of integer seconds and second fractions part)
    static timespec time_sub(const timespec& x, const timespec& y) {
        timespec result = y;
        // Perform the carry for the later subtraction by updating y.
        if (x.tv_nsec < y.tv_nsec) {
            int nsec = (y.tv_nsec - x.tv_nsec) / 1000000000l + 1;
            result.tv_nsec -= 1000000000l * nsec;
            result.tv_sec += nsec;
        }
        if (x.tv_nsec - result.tv_nsec > 1000000000l) {
            int nsec = (x.tv_nsec - result.tv_nsec) / 1000000000l;
            result.tv_nsec += 1000000000l * nsec;
            result.tv_sec -= nsec;
        }

        // tv_usec is now certainly positive.
        result.tv_sec = x.tv_sec - result.tv_sec;
        result.tv_nsec = x.tv_nsec - result.tv_nsec;
        return result;
    }

    /// result = x + y
    static timespec time_add(const timespec& x, const timespec& y) {
        timespec result;
        result.tv_sec = x.tv_sec + y.tv_sec;
        result.tv_nsec = x.tv_nsec + y.tv_nsec;
        while (result.tv_nsec > 1000000000l) {
            result.tv_nsec -= 1000000000l;
            result.tv_sec++;
        }
            
        return result;
    }


    /// result = (x - y) + z
    timespec time_subadd(const timespec& x, const timespec& y, const timespec& z) {
        timespec result;
        result.tv_sec = z.tv_sec + x.tv_sec - y.tv_sec;
        result.tv_nsec = z.tv_nsec + x.tv_nsec - y.tv_nsec;
        while (result.tv_nsec > 1000000000l) {
            result.tv_nsec -= 1000000000l;
            result.tv_sec++;
        }
        while (result.tv_nsec < 0) {
            result.tv_nsec += 1000000000l;
            result.tv_sec--;
        }
            
        return result;
    }
};


/**
    PrecisionTimer is a class for precision timing (microseconds in Windows, nanoseconds in Linux)
    It can be used by itself (not recommended) or in conjunction with ScopeTimer
**/
class PrecisionTimer {
	TimeType startTime, totalTime, endTime;
    
public:
	PrecisionTimer() {
		totalTime.init();
	}
	
	/// reset total time recorded by this timer
	inline void reset() {
	    totalTime.init();
	}
	
	/// start new timing interval that will be added to total timer value
	inline void start() {
		startTime.getSysTime();
	}
	
	/// read absolute time
	inline TimeType read() {
		TimeType temp;
		temp.getSysTime();
		return temp;
	}
	
	/// pause timer (be careful to call start and pause in pairs!)
	inline void pause() {	
		endTime.getSysTime();
		totalTime += endTime - startTime;
	}
	
	/// returns last time interval recorded by this timer (not total time recorded!)
	inline TimeType lastTime() const {
	    TimeType temp;
	    temp.init();
	    temp += endTime - startTime;
	    return temp;
	}
	
	const PrecisionTimer& operator += (const PrecisionTimer& other) {
	    totalTime += other.totalTime;
	    return *this;
	}
	
	/// returns last time interval recorded by this timer (not total time recorded!)
	inline double lastSeconds() const {
	    return lastTime().toSeconds();
	}
	
	/// returns total recorded time (if timer is started and not paused when this is called, the last - unfinished interval is not counted)
	inline operator TimeType() {
		return totalTime;
	}

    /// returns total time recorded in seconds (if timer is started and not paused when this is called, the last - unfinished interval is not counted)
	inline double totalSeconds() const {
		return totalTime.toSeconds();
	}
};


/**
    ScopeTimer creates a scope timer that adds the duration of current scope to selected PrecisionTimer
**/
class ScopeTimer {
	PrecisionTimer& pt;
	friend class ExcludeScopeFromTimer;
	
public:
	inline ScopeTimer(PrecisionTimer& t) : pt(t) {
		pt.start();
	}
	
	inline ~ScopeTimer() {
		pt.pause();
	}
};


/**
    ExcludeScopeTimer may be used after a ScopeTimer has been declared, to exclude a certain part (another scope) from timing
**/
class ExcludeScopeFromTimer {
	PrecisionTimer& pt;
	
public:
	inline ExcludeScopeFromTimer(PrecisionTimer& t) : pt(t) {
		pt.pause();
	}
	
	inline ~ExcludeScopeFromTimer() {
		pt.start();
	}
};


template<class C>
std::basic_string<C> timeToString(double seconds, size_t maxW = 7) {
    typedef std::basic_string<C> String;
    typedef std::basic_ostringstream<C> OStream;
    
    String unit("s");
    String empty("");
    size_t maxNr = maxW-1;
    
    if (seconds < 0.000001) {
        seconds *= 1000000000;
        unit = "ns";
        maxNr--;
    } else if (seconds < 0.0001) {
        seconds *= 1000000;
        unit = "µs";
        maxNr--;
    } else if (seconds < 0.1) {
        seconds *= 1000;
        unit = "ms";
        maxNr--;
    } else if (seconds > 359996400) {
        seconds /= (3600*24*365);
        unit = "y";
    } else if (seconds > 999999) {
        seconds /= 3600;
        unit = "h";
    } 
    
    OStream temp;
    temp.precision(maxNr);
    temp << seconds;
    if (temp.str().size() > maxNr) {
        temp.rdbuf()->str(empty);
        temp.precision(maxNr-1);
        temp << seconds;
        if (temp.str().size() > maxNr) {
            temp.rdbuf()->str(empty);
            temp.precision(maxNr-2);
            temp << seconds;
            if (temp.str().size() > maxNr) {
                temp.rdbuf()->str(empty);
                temp.precision(0);
                temp << std::scientific << seconds;
                if (temp.str().size() > maxNr) {
                    temp.rdbuf()->str(empty);
                    temp << std::string("---------", maxNr);
                }
            }
        }
    }
    temp << unit;
    return temp.str();
}


#undef getTime
#undef time_type_init
#undef addDifference


#endif // TIMER_H_INCLUDED


template<class CharPtr>
CharPtr skipPath(CharPtr text) {
    CharPtr lastSlash = text, t = text;
    while (*t != 0) {
        if (*t == '\\' || *t == '/') lastSlash = t+1;
        ++t; 
    }
    return lastSlash;
}

#ifdef NDEBUG1
    #define TRACE(text, mask, level) {}
    #define TRACEVAR(var, mask, level) {}
#else
    #define IF_TRACABLE(mask, level) if ((level <= TRACE_LEVEL) && ((TRACE_MASK & mask) > 0))
    #define TRACE(text, mask, level) {IF_TRACABLE(mask, level) std::cout  << "(" << level << ") " << skipPath(__FILE__) << " line " << __LINE__ << ": " << text << "\n";}
    #define TRACEVAR(var, mask, level) {IF_TRACABLE(mask, level) std::cout  << "(" << level << ") " << skipPath(__FILE__) << " line " << __LINE__ << ": " << #var << "=" << var << "\n";}
#endif

    
#define TRACE_MASK_THREAD 0x10
#define TRACE_MASK_CLIQUE 0x20
#define TRACE_MASK_INITIAL 0x40



template<class Ostream, class T>
Ostream& operator<< (Ostream& out, const std::vector<T>& vec) {
    if (vec.size() > 0) {
        out << "[l577:" << vec[0];
        for (size_t i = 1; i < vec.size(); ++i)
            out << "," << vec[i]; 
        out << "]";
    } else 
        out << "[/]";
    return out;
}


template<class Ostream, class T, class U>
Ostream& operator<< (Ostream& out, const std::vector<std::pair<T, U>>& vec) {
    if (vec.size() > 0) {
        out << "[<" << vec[0].first << "," << vec[0].second << ">";
        for (size_t i = 1; i < vec.size(); ++i)
            out << "; <" << vec[i].first << "," << vec[i].second << ">"; 
        out << "]";
    } else 
        out << "[/]";
    return out;
}

template<class T>
struct ClassName {
    static const char* getValue() {
        return "";
    }
};

#define REGISTER_CLASS_NAME(CLASS, NAME) template<> struct ClassName<CLASS> {static const char* getValue() {return NAME;} };
#define REGISTER_TEMPLATE1_CLASS_NAME(CLASS, NAME) template<class T> struct ClassName<CLASS<T>> {static const char* getValue() {return NAME;} };
#define REGISTER_TEMPLATE2_CLASS_NAME(CLASS, NAME) template<class T1, classT2> struct ClassName<CLASS<T1,T2>> {static const char* getValue() {return NAME;} };
#define REGISTER_TEMPLATE_EXT_CLASS_NAME(CLASS, NAME) template<class T> struct ClassName<CLASS<T>> {static const char* getValue() { \
    static std::unique_ptr<std::string> rs(new std::string); \
    if (rs->empty()) { std::stringstream s; s << NAME << "<" << ClassName<T>::getValue() << ">"; *rs = s.str(); } return rs->c_str(); } };


template<
    class VertexRepresentation,
    class SetRepresentation, 
    class GraphRepresentation, 
    class Sorter
>
class MaximumCliqueProblem : Sorter {
public:
    typedef VertexRepresentation VertexId;
    typedef SetRepresentation VertexSet;
    
protected:
    using Sorter::notEmpty;
    using Sorter::topNumber;
    using Sorter::topVertex;
    using Sorter::popTop;
    using Sorter::numberSort;
    using Sorter::initialSort;
    typedef typename Sorter::NumberedSet NumberedSet;
        
    std::string algorithmName;
    GraphRepresentation* graph;
    VertexId n;                         // number of vertices
    unsigned int maxSize;               // size of max clique
    VertexSet maxClique;
    PrecisionTimer timer;
    unsigned long long int steps;
    KillTimer1 killTimer;
    
public:
    VertexSet knownC;

    MaximumCliqueProblem(GraphRepresentation& graph) : graph(&graph), n(graph.getNumVertices()), maxSize(0), steps(0) {}
    
    // get the result of the search - maximal clique of the provided graph
    const VertexSet& getClique() const {return maxClique;}
    
    // output statistics of the last search (mostly timer readings); colored=true produces colored text on terminals supporting ANSI escape sequences
    void outputStatistics(core::Vector_sp clique,bool colored = true, bool verbose = false) {
        if (verbose) {
            std::ostringstream algorithmName;
            auto basefmt = std::cout.flags();
            auto baseFill = std::cout.fill();
            algorithmName << "MC(" << ClassName<VertexId>::getValue() << ","
                          << ClassName<VertexSet>::getValue() << "," 
                          << ClassName<GraphRepresentation>::getValue() << "," << ClassName<Sorter>::getValue() << ") ";
            std::cout << "-- " << std::setw(80-3) << std::setfill('-') << std::left << algorithmName.str();
            std::cout.flags(basefmt);
            std::cout << std::setfill(baseFill) << "\n";
            std::cout << "search took " << timer.totalSeconds() << "s; " << steps << " steps\n";
        }
        VertexSet c = getClique();
        graph->remap(c);
        if (verbose) {
            if (wasSearchInterrupted()) std::cout << "Warning, search has been interrupted, the results might not be correct\n";
            std::cout << "Clique (" << getClique().size() << "): " << c;
            std::cout << "getClique().size() -> " << getClique().size() << "\n";
        }
        copy_to_vector(clique,c);
    }
    
    // run the search for max clique
    void search() {
        killTimer.start(10000);
        ScopeTimer t(timer);
        VertexSet c; // clique
        VertexSet p; // working set of vertices
        NumberedSet numbers;       // ordered set of colors
        
        c.reserve(n);
        p.reserve(n);
        for (VertexId i = 0; i < n; i++) 
            p.add(i);
            
        Sorter::init(graph);
        initialSort(c, p, numbers);
        
        TRACEVAR(c, TRACE_MASK_INITIAL, 2); 

        // some initial sorts (e.g. MCR) also find a clique
        maxSize = c.size();
        if (maxSize > 0)
            saveSolution(c);
                            
        if (numbers.size() == 0)
            numberSort(c, p, numbers, maxSize);
        TRACEVAR(numbers, TRACE_MASK_INITIAL, 2); 
        c.clear();
        
        expand(c, p, numbers);
        
        killTimer.cancel();
    }
    
    bool wasSearchInterrupted() const {return killTimer.timedOut;}
    
protected:   
    // main recursive function 
    //  c ... clique candidate, 
    //  p ... set of candidate vertices that can be added to c,
    //  numbered ... an ordered and numbered set of vertices
    void expand(VertexSet& c, VertexSet& p, NumberedSet& np) {
        ++steps;
        while (notEmpty(np)) {
            if (c.size() + topNumber(np) <= maxSize || killTimer.timedOut) {return;}
            auto v = topVertex(np, p);
            //std::cout << "v" << v << " ";
            c.add(v);
            VertexSet p1;
            graph->intersectWithNeighbours(v, p, p1);
            
            if (p1.size() == 0) {
                if (c.size() > maxSize) 
                    saveSolution(c);
            } else {
                NumberedSet np1;
                numberSort(c, p1, np1, maxSize);
                expand(c, p1, np1);
            }
            
            c.remove(v);
            popTop(np, p);
            //std::cout << "p" << np.size() << " ";
        }
    }
    
    // when a clique, larger than its predecessor is found, call this function to store it
    void saveSolution(const VertexSet& c) {
        maxSize = c.size();
        maxClique = c;
        //std::cout << "new clique " << maxSize << "\n";
    }
};

template<class T>
struct isTypeASet {
    enum {value = false};
};

template<class VectorSetRepresentation>
class Graph {
public:
    std::vector<VectorSetRepresentation> adjacencyMatrix, invAdjacencyMatrix;
  gctools::Vec0<int> degrees;
    std::vector<int> mapping;
    
    typedef typename VectorSetRepresentation::VertexId VertexId;
    typedef VectorSetRepresentation VertexSet;
    
    size_t getNumVertices() const {return adjacencyMatrix.size();}
  void init(const gctools::Vec0<gctools::Vec0<char> >& adjacency, const gctools::Vec0<int>& d) {
        size_t n = adjacency.size();
        adjacencyMatrix.resize(n); 
        invAdjacencyMatrix.resize(n); 
        for (size_t i = 0; i < n; ++i) {
            adjacencyMatrix[i].resize(n, false);
            invAdjacencyMatrix[i].resize(n, false);
            invAdjacencyMatrix[i][i] = false;
            for (size_t j = i+1; j < adjacency[i].size(); ++j) {
                adjacencyMatrix[i][j] = adjacency[i][j];
                invAdjacencyMatrix[i][j] = (adjacency[i][j] == false);
            }
        }
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < i; ++j) {
                adjacencyMatrix[i][j] = adjacency[j][i];
                invAdjacencyMatrix[i][j] = (adjacency[j][i] == false);
            }
        }
        degrees = d;
        mapping.clear();
    }
    
    // perform intersection, "result" keeps the ordering of the set "vertices"
    void intersectWithNeighbours(VertexId p, const VertexSet& vertices, VertexSet& result) const {
        // global function intersectWithAdjecency(VectorSetRepresentation, VectorSetRepresentation, VectorSetRepresentation) must be specified
        intersectWithAdjecency(vertices, adjacencyMatrix[p], result);
    }
    
    bool intersectionExists(VertexId p, const VertexSet& vertices) const {
        size_t n = vertices.size();
        for (size_t i = 0; i < n; ++i) {
            if (adjacencyMatrix[p][vertices[i]])
                return true;
        }
        return false;
    }
    
    // change the order of vertices in the adjacency matrix (renumber them)
    template<class Vec>
    void orderVertices(const Vec& order) {
        // check order vector
        size_t n = getNumVertices();
        if (order.size() != n) 
            throw "Invalid size vector in orderVertices";
            
        // create a vertex mapping table that will be used to renumber vertices back to original
        if (mapping.size() == 0) {
            // create default mapping that maps i → i
            mapping.resize(n);
            for (size_t i = 0; i < n; ++i)
                mapping[i] = i;
        }
        decltype(mapping) mapping2(n);
        
        // remap to temporary adjacencyMatrix
        std::vector<VectorSetRepresentation> adjacencyMatrix2;
        adjacencyMatrix2.resize(n);
        for (size_t i = 0; i < n; ++i) {
            adjacencyMatrix2[i].resize(n);
            invAdjacencyMatrix[i].clear();
            invAdjacencyMatrix[i].resize(n);
            mapping2[i] = mapping[order[i]];
            auto& adjRowI = adjacencyMatrix[order[i]];
            for (size_t j = 0; j < n; ++j) {
                adjacencyMatrix2[i][j] = adjRowI[order[j]] == true;
                invAdjacencyMatrix[i][j] = (i != j) & !adjacencyMatrix2[i][j];
            }
            // the line above includes the condition (i != j) because:
            // adjacency inverse is used to filter out vertices (operator &) and it is useful if 
            // given a vertex, it filters out its neighbours as well as the vertex itself
            // therefore reset the edge linking vertices to themselves
        }
        std::swap(adjacencyMatrix2, adjacencyMatrix);
        std::swap(mapping2, mapping);
        auto oldDeg = degrees;
        for (size_t i = 0; i < degrees.size(); ++i) degrees[i] = oldDeg[order[i]];
    }
    
    // in-place remapping function
    void remap(VectorSetRepresentation& v) {
        if (mapping.size() == 0 || v.size() == 0) return;
        
        if (isTypeASet<VectorSetRepresentation>::value) {
            VectorSetRepresentation rv;
            rv.reserve(mapping.size());
            
            if (mapping.size() < v.size())
                throw "Mapping failed, mapping is not known for all vertices";
            
            for (size_t i = 0; rv.size() < v.size(); ++i) {
                if (v[i]) rv.add(mapping[i]);
            }
            std::swap(rv, v);
        } else
            for (size_t i = 0; i < v.size(); ++i) v[i] = mapping[v[i]];
    }
};
REGISTER_TEMPLATE1_CLASS_NAME(Graph, "Graph");

template<class G, class S>
bool isClique(const G& graph, const S& v) {
    size_t n = v.size();
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i+1; j < n; ++j) {
            if (graph.adjacencyMatrix[v[i]][v[j]] == false) {
                //std::cout << v[i] << " and " << v[j] << " are not connected\n" << graph.adjacencyMatrix[v[i]][v[j]] << "\n";
                return false;
            }
        }
    }
    return true;
}

// Set of vertices is based on std::vector
// used as VertexSetRepresentation in MaximumCliqueProblem
template<class T>
class VectorSet : std::vector<T> {
public:
    typedef T VertexId;
    
    using std::vector<T>::resize;
    using std::vector<T>::reserve;
    void add(const T& value) {push_back(value);}
    T pop() {T temp=this->back(); this->pop_back(); return temp;}
    void remove(const T& value) {
        if (this->back() == value) this->pop_back();
        else {
            auto temp = std::find(this->begin(), this->end(), value);
            if (temp != this->end())
                this->erase(temp);
        }
    }
    using std::vector<T>::size;
    using std::vector<T>::operator[];
    using std::vector<T>::clear;
    using std::vector<T>::back;
    template<class AdjSet>
    friend void intersectWithAdjecency (const VectorSet& v, const AdjSet& adj, VectorSet& result) {
        auto n = v.size();
        result.reserve(n);
        for (size_t i = 0; i < n; ++i) {
            if (adj[v[i]])
                result.add(v[i]);
        }
    }

    // only required for debugging
    bool isIntersectionOf(const VectorSet& bigSet) {
        size_t n = bigSet.size();
        if (n < 1) return false;
        --n;
        for (size_t i = 0; i < size(); ++i) {
            for (size_t j = 0; bigSet[j] != (*this)[i]; ++j) {
                if (j == n) return false;
            }
        }
        return true;
    }
};
REGISTER_TEMPLATE1_CLASS_NAME(VectorSet, "Vector based set");

template<class Ostream, class T>
Ostream& operator<< (Ostream& out, const VectorSet<T>& vec) {
    if (vec.size() > 0) {
        out << "[l929:" << vec[0];
        for (size_t i = 1; i < vec.size(); ++i)
            out << "," << vec[i]; 
        out << "]";
    } else 
        out << "[/]";
    return out;
}

REGISTER_CLASS_NAME(int, "int");
REGISTER_CLASS_NAME(short int, "short int");
REGISTER_CLASS_NAME(char, "char");

#endif // MAXIMUMCLIQUEBASE_H



// #include "ParallelMaximumClique.h"
/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 */

#ifndef HEADER_DEFAAFC421A53A4F
#define HEADER_DEFAAFC421A53A4F

#include <thread>
#include <memory>
#include <deque> // TODO: test deque vs. vector for storing a queue of jobs


using std::swap;


/**
    VertexRepresentation should be
        of the most efficient scalar type on the given hardware
    VertexSetRepresentation must support:
        resize(size)
        reserve(max_size)
        add(value) // vector.push_back or set.operator[]
        pop(value) // takes last added element (that equals to v) from the set
        pop() // returns the value of the last added element, and also removes it from the set (only in vectors)
        uint size()
        operator[] (index)
        clear() // erase all elements
        TODO: list others
    Vector should be a template class (like std::vector) for dynamicly-sized arrays 
    
    function intersect(VertexSet, VertexId, VertexSet& result)
**/
template<
    class Vertex_t,
    class Set_t, 
    class Graph_t, 
    class Sort_t,
    template<class S> class InitialSort_t
>
class ParallelMaximumCliqueProblem : public InitialSort_t<Sort_t> {
public:
    typedef Vertex_t VertexId;
    typedef Set_t VertexSet;
    typedef Graph_t Graph;
    typedef Sort_t Sorter;
    typedef InitialSort_t<Sort_t> InitialSorter;
    typedef typename Sorter::NumberedSet NumberedSet;
    
    struct Job {
        VertexSet c;            // current clique
        VertexSet vertices;     // available vertices
        NumberedSet numbers;    // available vertices
        unsigned int estimatedMax;       // used to queue jobs from the most perspective down to the least perspective
        
        Job()  {}
        Job(const VertexSet& cin, const VertexSet& vin, const NumberedSet& nin, int estimate) : c(cin), vertices(vin), numbers(nin), estimatedMax(estimate)  {}
        void set(const VertexSet& cin, const VertexSet& vin, const NumberedSet& nin, int estimate) {c = cin; vertices = vin; numbers = nin; estimatedMax = estimate;}
        friend void swap(Job& a, Job& b) {
            std::swap(a.c, b.c); 
            std::swap(a.vertices, b.vertices);
            std::swap(a.numbers, b.numbers);
            std::swap(a.estimatedMax, b.estimatedMax);
        }
    };

    struct Worker : Sorter {
        // parent pMCP - the problem for which this worker is working
        ParallelMaximumCliqueProblem* parent;
        // Timer and counters
        PrecisionTimer timer;
        unsigned long long int steps;
        // the same graph as the parent pMCP
        Graph* graph;
        // id of the worker (helps with the end-time statistics)
        int id;
        // copy of parent maxSize (size of the active maximumClique)
        unsigned int localMaxSize;
        // break the work on a job when going below the specified level (size of a clique); search will not continue on a level lower than jobLevelBreak
        unsigned int jobLevelBreak;
        
        using Sorter::notEmpty;
        using Sorter::topNumber;
        using Sorter::topVertex;
        using Sorter::popTop;
        using Sorter::numberSort;
        
        Worker() : graph(nullptr) {}
        
        void setup(ParallelMaximumCliqueProblem* pmcp, int newId) {
            parent = pmcp;
            graph = pmcp->graph;
            id = newId;
            Sorter::init(graph);
            steps = 0;
        }
        
        void threadFunc() {
            localMaxSize = parent->maxSize;
            TRACE("threadFunc start", TRACE_MASK_THREAD, 1);
            int jobsDone = 0;
            try {
                { // this scope is for the scope timer in the next line: the scope must end before the function sends the results to the parent thread
                    ScopeTimer t(timer);
                    // work while there are jobs available
                    TRACE("threadFunc: getting a job", TRACE_MASK_THREAD, 1);
                    while (true) {
                        Job job;
                        // grab a job from the queue
                        {  
                            std::lock_guard<std::mutex> lk(parent->mutexJobs); 
                            TRACEVAR(parent->activeJobs, TRACE_MASK_THREAD, 2);
                            // if there is no queued jobs then wait; if there is also no active jobs (that would be adding to the queue) then quit
                            if (parent->jobs.size() == 0) {
                                /*
                                // when dynamic tree splitting is used, try the following:
                                if (parent->activeJobs == 0)
                                    break;
                                else
                                    continue;
                                    */
                                break;
                            }
                            
                            // take a job from the top of the queue, check if a new job can be added to the queue immediately 
                            if ((parent->jobs.size() > parent->maxJobs) || (parent->jobs.back().estimatedMax <= localMaxSize + 5)) {
                                // job on the queue is a small one; do not split it further (5 - the number to check against - is arbitrary)
                                // (swap is the fastest operation, job is empty now)
                                swap(job, parent->jobs.back());
                                parent->jobs.pop_back();
                                parent->requireMoreJobs = true;
                                // mark where to end execution
                                jobLevelBreak = job.c.size();
                            } else {
                                // copy the job out of the queue
                                job = parent->jobs.back();
                                // change the job in the queue to make it available for another worker
                                this->popTop(parent->jobs.back().numbers, parent->jobs.back().vertices);
                                parent->jobs.back().estimatedMax = topNumber(parent->jobs.back().numbers) + parent->jobs.back().c.size();
                                // mark where to end execution
                                jobLevelBreak = job.c.size()+1;
                            }
                            // mark that another job has been activated
                            ++parent->activeJobs;
                        }
                        TRACE("threadFunc: calling expand", TRACE_MASK_THREAD, 1);
                        // work on the job
                        {
                            TRACEVAR(jobLevelBreak, TRACE_MASK_THREAD, 1);
                            // dig in
                            expand(job);
                            ++jobsDone;
                        }
                        TRACE("threadFunc: job completed", TRACE_MASK_THREAD, 1);
                        // deactivate the completed job
                        {
                            std::lock_guard<std::mutex> lk(parent->mutexJobs);
                            --parent->activeJobs;
                        }
                    }
                } // timer scope end           
                // send statistics to the parent
                {
                    TRACE("threadFunc: reposting stats", TRACE_MASK_THREAD, 1);
                    std::lock_guard<std::mutex> lk(parent->mutexJobs); 
                    parent->workerActiveTimes.push_back(timer.totalSeconds());
                    parent->workerSteps.push_back(steps);
                }
            } catch (const char* e) {
                std::cout << "Exception in thread: " << e << "\n";
            }
        }
        
        // main recursive function (parallel)
        void expand(Job& job) {
            ++steps;
            TRACE("expand start", TRACE_MASK_CLIQUE, 2);
            
            localMaxSize = parent->maxSize;
            TRACEVAR(localMaxSize, TRACE_MASK_CLIQUE, 2);
            while (notEmpty(job.numbers)) {
                TRACEVAR(job.numbers.size(), TRACE_MASK_CLIQUE, 2);
                if (job.estimatedMax <= localMaxSize || parent->killTimer.timedOut) {return;}
                Job newJob;
                auto v = topVertex(job.numbers, job.vertices);
                popTop(job.numbers, job.vertices);
                job.estimatedMax = job.c.size() + topNumber(job.numbers);
                graph->intersectWithNeighbours(v, job.vertices, newJob.vertices);
                
                newJob.c = job.c;
                newJob.c.add(v);
                if (newJob.c.size() > localMaxSize) { // condition (newJob.vertices.size() == 0) is left out to make maxSize up to date at all times
                    localMaxSize = parent->saveSolution(newJob.c);
                }
                
                TRACEVAR(newJob.vertices.size(), TRACE_MASK_CLIQUE, 2);
                if (newJob.vertices.size() > 0) {
                    // number vertices
//                    std::cout << "DBG " << std::flush;
                    newJob.numbers.resize(newJob.vertices.size());
//                    std::cout << "DBG " << newJob.c.size() << " " << newJob.vertices.size() << " " << newJob.numbers.size() << " " << localMaxSize << "\n";
                    numberSort(newJob.c, newJob.vertices, newJob.numbers, localMaxSize);
                    newJob.estimatedMax = newJob.c.size() + topNumber(newJob.numbers);
                    TRACEVAR(newJob.estimatedMax, TRACE_MASK_CLIQUE, 2);
                    
                    /*
                    // quick and dirty (not protected by mutex but should be atomic) check weather new job should be pushed to the queue
                    //bool jobQueueFull = (jobToFork == nullptr) || (!parent->requireMoreJobs);
                    if (parent->requireMoreJobs) { 
                        ScopeTimer tm(timers[timer_job_management]);
                        std::lock_guard<std::mutex> lk(parent->mutexJobs); 
                        // check the number of jobs in queue again (now protected by mutex)
                        if (parent->jobs.size() < parent->maxJobs) {
                            TRACE("expand adding a job to the queue");
                            TRACEVAR(parent->jobs.size());
                            // fork the pre-prepared job (far more appropriate than forking current job)
                            parent->addJob(*jobToFork);
                            // mark the forked job as off-limits to this worker
                            jobLevelBreak = jobToFork->c.size()+1;
                            TRACEVAR(jobLevelBreak);
                            // disable further forks; better alternative would be to allow further forks on a higher level
                            // it would be harder to implement and would introduce additional overhead (unless the main recursion was transformed to iteration)
                            jobToFork = nullptr;
                        }
                    }
                    */
                    
                    // continue exploration of the search tree within this job
                    TRACE("expand calling expand", TRACE_MASK_CLIQUE, 2);
                    expand(newJob);
                    if (job.c.size() < jobLevelBreak) {
                        TRACE("breaking worker search tree", TRACE_MASK_CLIQUE, 2);
                        TRACEVAR(jobLevelBreak, TRACE_MASK_CLIQUE, 2); 
                        return;
                    }
                }
            }
            TRACE("expand end", TRACE_MASK_CLIQUE, 2);
        }
    };

    
protected:
    using InitialSorter::initialSort;
        
    std::string algorithmName;
    Graph* graph;
    VertexId n;                         // number of vertices
    unsigned int maxSize;               // size of max clique
    unsigned int numThreads;            // stores the number of threads used in the last search (where this number a parameter to the function)
    VertexSet maxClique;
    std::deque<Job> jobs;
    bool requireMoreJobs;
    size_t maxJobs, activeJobs;
    PrecisionTimer timer;
    std::vector<double> workerActiveTimes;
    std::vector<unsigned long long> workerSteps;
    std::mutex mutexJobs, mutexQ;
    KillTimer1 killTimer;
    
public:
    VertexSet knownC;

    ParallelMaximumCliqueProblem(Graph& graph) : graph(&graph), n(graph.getNumVertices()), maxSize(0) {}
    
    // get the result of the search - maximal clique of the provided graph
    const VertexSet& getClique() const {return maxClique;}
    
    // output statistics of the last search (mostly timer readings); colored=true produces colored text on terminals supporting ANSI escape sequences
    void outputStatistics(core::Vector_sp clique, bool colored = true, bool verbose=false) {
        if (verbose) {
            std::ostringstream algorithmName;
            auto basefmt = std::cout.flags();
            auto baseFill = std::cout.fill();
            algorithmName << "pMC[" << numThreads << " threads, " << maxJobs << " jobs](" << ClassName<VertexId>::getValue() << ","
                          << ClassName<VertexSet>::getValue() << "," << ClassName<Graph>::getValue() << "," << ClassName<InitialSorter>::getValue() << ") ";
            std::cout << "-- " << std::setw(80-3) << std::setfill('-') << std::left << algorithmName.str();
            std::cout.flags(basefmt);
            std::cout << std::setfill(baseFill) << "\n";
            unsigned long long steps = 0;
            for (size_t i = 0; i < workerSteps.size(); ++i) {
                steps += workerSteps[i];
            }
            std::cout << "search took " << timer.totalSeconds() << "s; " << steps << " steps\n";
        }
        VertexSet c = getClique();
        graph->remap(c);
        if (verbose) {
            auto basefmt = std::cout.flags();
            auto baseFill = std::cout.fill();
            if (wasSearchInterrupted()) std::cout << "Warning, search has been interrupted, the results might not be correct\n";
            std::cout << "Clique (" << getClique().size() << "): " << c;
            std::cout.flags(basefmt);
            std::cout << std::setfill(baseFill);
            std::cout << "getClique().size() -> " << getClique().size() << "\n";
        }
        copy_to_vector(clique,c);
    }
    
    double workerEfficiency() const {
        double maxTime = timer.totalSeconds();
        double eff = 0;
        for (size_t i = 0; i < workerActiveTimes.size(); ++i) {
            eff += workerActiveTimes[i];
        }
        return eff / (maxTime*workerActiveTimes.size());
    }
    
    void sortLastJob() {
        if (jobs.size() > 1) {
            for (auto j = jobs.end()-1; (j != jobs.begin()) && ((j-1)->vertices.size() > j->vertices.size()); --j)
                swap(*j, *(j-1));
        }
    }
    
    void addJob(Job& job) {
        // LIFO queue, less copying
        jobs.push_front(Job());
        swap(jobs.front(), job);
        requireMoreJobs = jobs.size() == maxJobs;
        
        /*
        // Jobs queue is of LIFO type
        jobs.push_front(job);
        requireMoreJobs = jobs.size() == maxJobs;
        */
        
        // Job sorting, it works better without one :(
        //jobs.push_back(job);
        //sortLastJob();
    }
    
    void swapJob(Job& job) {
        if (jobs.size() > 0 && job.estimatedMax < jobs.back().estimatedMax) {
            //swap(job, jobs.back());
            //sortLastJob();
        }
        /* // DEBUG
        static int cnt = 1;
        if (cnt > 0) {
            --cnt;
            for (auto& j: jobs)
                std::cout << j.estimatedMax << " ";
            std::cout << "\n";
        }*/
    }
    
    // run the search for max clique
    void search(unsigned int numThreads, unsigned int numJobs, std::vector<int>& affinities) {
        killTimer.start(10000);
        ScopeTimer t(timer);
        VertexSet c; // clique
        VertexSet p; // working set of vertices
        NumberedSet numbers;       // ordered set of colors
        this->numThreads = numThreads;
        
        {
            // setting the order of vertices (initial sort that renumbers vertices in the input graph)
            c.reserve(n);
            p.reserve(n);
            for (VertexId i = 0; i < n; i++) 
                p.add(i);
            
            InitialSorter::init(graph);
            TRACE(typeid(InitialSorter).name(), TRACE_MASK_CLIQUE, 1);
            initialSort(c, p, numbers);
            
            // some initial sorts (e.g. MCR) also find a clique
            if (c.size() > 0) {
                saveSolution(c);
                c.clear();
            }
        }
        
        if (numbers.size() == 0) {
            // if initial sort did not setup "numbers", numberSort must be called
            this->numberSort(c, p, numbers, maxSize);
        }
        
        c.clear();
        
        // Generate root job
        maxJobs = std::min(numJobs, numThreads*1000);
        //jobs.reserve(maxJobs+1);
        jobs.push_back(Job(c, p, numbers, c.size() + this->topNumber(numbers)));
        activeJobs = 0;
        requireMoreJobs = jobs.size() < maxJobs;
        
        // Create threads and workers
        std::vector<std::unique_ptr<std::thread> > threads;
        std::vector<Worker> workers;
        threads.resize(numThreads);
        workers.resize(numThreads);
        
        // associate thread & worker pairs
        for (unsigned int i = 0; i < numThreads; ++i) {
            TRACE("setting up thread", TRACE_MASK_THREAD, 1);
            TRACEVAR(i, TRACE_MASK_THREAD, 1);
            workers[i].setup(this, i);
            threads[i] = std::unique_ptr<std::thread>(new std::thread([&workers, i](){workers[i].threadFunc();}));
        }
        TRACE("Done building threads, waiting for join", TRACE_MASK_THREAD, 1);
        // wait for all the workers to finish
        for (unsigned int i = 0; i < numThreads; ++i) {
            threads[i]->join();
            TRACE("Thread joined to main thread", TRACE_MASK_THREAD, 1);
        }
        killTimer.cancel();
        jobs.clear();
        TRACE("search: end", TRACE_MASK_THREAD, 1)
    }
    
    bool wasSearchInterrupted() const {return killTimer.timedOut;}
    
    void debug() {
        std::cout << "Parallel Maximum Clique problem DEBUG:\n";
        std::cout << "   " << algorithmName << "\n";
        std::cout << "   graph @ " << graph << "\n";
        std::cout << "   number of vertices " << n << "\n";                         // number of vertices
        std::cout << "   maxSize " << maxSize << "\n";               // size of max clique
        std::cout << "   numThreads " << numThreads << "\n";            // stores the number of threads used in the last search (where this number a parameter to the function)
        std::cout << "   maxClique " << maxClique.size() << "\n";
        std::cout << "   numJobs " << jobs.size() << "\n";
        std::cout << "   requireMoreJobs " << requireMoreJobs << "\n";
        std::cout << "   maxJobs/activeJobs " << maxJobs << "/" << activeJobs << "\n";
        std::cout << "   numWorkerStatsTimes " << workerActiveTimes.size() << "\n";
        std::cout << "   numWorkerStatsSteps " << workerSteps.size() << "\n";
        std::cout << "DEBUG END\n";
    }
    
protected:
    // when a clique, larger than its predecessor is found, call this function to store it
    unsigned int saveSolution(const VertexSet& c) {
        unsigned int ret;
        // make a copy of clique
        TRACE("Saving solution", TRACE_MASK_CLIQUE, 2);
        TRACEVAR(c.size(), TRACE_MASK_CLIQUE, 2);
        {
            std::lock_guard<std::mutex> lk(mutexQ); 
            if (maxSize < c.size()) {
                maxSize = c.size();
                maxClique = c;
            }
            ret = maxSize;
        }
        
        return ret;
    }
};


#endif // header guard 



// #include "BB_GreedyColorSort.h"
/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 */

#ifndef BB_GREEDYCOLORSORT_H_INCLUDED
#define BB_GREEDYCOLORSORT_H_INCLUDED


// #include "BitSet.h"

/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 */

#ifndef BITSET_H
#define BITSET_H


// some algorithms are from: http://graphics.stanford.edu/~seander/bithacks.html#IntegerMinOrMax


#include <sstream>
#include <vector>
#include <bitset>
#include <algorithm>
#include <numeric>
#include <iostream>


template<unsigned int I>
struct static_log2 {
    enum {
        value = static_log2<(I >> 1)>::value+1
    };
};

template<>
struct static_log2<1> {
    enum {
        value = 0
    };
};


struct LogTableSetter {
    LogTableSetter(char* logTable256) {
        logTable256[0] = logTable256[1] = 0;
        for (int i = 2; i < 256; i++) 
            logTable256[i] = 1 + logTable256[i >> 1];
        logTable256[0] = -1; // if you want log(0) to return -1
    }
};

// log2 does not detect negative numbers and returns max(char) if zero is passed to it
template<class Int>
unsigned int log2(Int v) {
    static char logTable256[256];
    static LogTableSetter dummy(logTable256);

    unsigned int shr = sizeof(Int) << 2;
    unsigned int ofs = 0;
    Int v2;

    while (shr >= 8) {
        v2 = v >> shr;
        if (v2) {
            ofs += shr;
            v = v2;
        }
        shr >>= 1;
    }
    return ofs + logTable256[v];
}

// reverses the order of bits
template<class Int>
Int reverseBits(Int v) {
    Int s = sizeof(v) << 3; // bit size
    Int mask = ~0;         
    while ((s >>= 1) > 0) {
        mask ^= (mask << s);
        v = ((v >> s) & mask) | ((v << s) & ~mask);
    }
    return v;
}

// find first bit that equals to 1 (only designed for 32 and 64 bit integers) that works on all machines and compilers
template<class Int>
unsigned int countTrailing0M(Int v) {
    if (sizeof(Int) < 8) {
        // smaller and possibly faster version for 32 and les bit numbers
        static const char multiplyDeBruijnBitPosition[32] = {
            0, 1, 28, 2, 29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4, 8, 
            31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6, 11, 5, 10, 9
        };
        return multiplyDeBruijnBitPosition[((uint32_t)((v & -v) * 0x077CB531U)) >> 27];
    } else {
        static const char multiplyDeBruijnBitPosition[64] = {
            0, 1, 2, 56, 3, 32, 57, 46, 29, 4, 20, 33, 7, 58, 11, 47, 
            62, 30, 18, 5, 16, 21, 34, 23, 53, 8, 59, 36, 25, 12, 48, 39, 
            63, 55, 31, 45, 28, 19, 6, 10, 61, 17, 15, 22, 52, 35, 24, 38, 
            54, 44, 27, 9, 60, 14, 51, 37, 43, 26, 13, 50, 42, 49, 41, 40
        };
        return multiplyDeBruijnBitPosition[((uint64_t)((v & -v) * 0x26752B916FC7B0DULL)) >> 58];
    };
}

// find first bit that equals to 1 (GCC only, undefined return if the operand equals 0 - [often 0 or the maximum number of bits])
template<class Int>
unsigned int countTrailing0(Int v) {
    if (sizeof(Int) < sizeof(int)) {return countTrailing0<int>(v);}
    if (sizeof(Int) == sizeof(int)) return __builtin_ctz(v);
    if (sizeof(Int) == sizeof(long int)) return __builtin_ctzl(v);
    if (sizeof(Int) >= sizeof(long long int)) return __builtin_ctzll(v);
}

// count number of bits that equal 1 (GCC only)
template<class Int>
unsigned int countOnes(Int v) {
    if (sizeof(Int) < sizeof(int)) {return countOnes<int>((int)v);}
    if (sizeof(Int) == sizeof(int)) return __builtin_popcount(v);
    if (sizeof(Int) == sizeof(long int)) return __builtin_popcountl(v);
    if (sizeof(Int) >= sizeof(long long int)) return __builtin_popcountll(v);
}

class BitSet {
protected:
    // how large are bit blocks [number of bits]
    static const unsigned int res = sizeof(unsigned long) * 8;
    // factor for shr when converting offset in bits to offset in data
    static const unsigned int f_shr = static_log2<res>::value;
    // mask for converting absolute offset in bits to relative offset inside a data cell
    static const unsigned int shift_mask = res-1;
    
    typedef std::bitset<res> b64;
    std::vector<b64> data;
    size_t numUsed, numAllocated;
    
    class BoolProxy {
        BitSet* parent;
        size_t index;
    
    public:
        BoolProxy(BitSet* p, size_t i) : parent(p), index(i) {}
        bool operator= (bool b) {parent->setValue(index, b); return b;}
        operator bool() const {return parent->getValue(index);}
    };
        
public:
    //~BitSet() {std::cout << data.size() << " : " << numAllocated << " : " << numUsed << "\n";}
    BitSet() : numUsed(0), numAllocated(0) {}
    BitSet(size_t n) : numUsed(0), numAllocated(0) {resize(n);}
    BitSet(const BitSet& other) : numUsed(other.numUsed), numAllocated(other.numAllocated) {copy(other);}
    const BitSet& operator= (const BitSet& other) {
        numUsed = other.numUsed; 
        numAllocated = other.numAllocated;
        copy(other); 
        return *this;
    }
    
    void resize(size_t newSize) {data.resize((newSize+shift_mask) >> f_shr); numAllocated = data.size() << f_shr; numUsed = newSize;}
    void resize(size_t newSize, bool value) {data.resize((newSize+shift_mask) >> f_shr); numAllocated = data.size() << f_shr; numUsed = newSize;}
    void reserve(size_t newSize) {numAllocated = ((newSize+shift_mask) >> f_shr) << f_shr; data.reserve(numAllocated >> f_shr);}
    size_t size() const {return numUsed;}
    
    // size of the bit blocks
    static unsigned int resolution() {return res;}
    
    // individual bit set/get
    bool operator[](size_t index) const {return getValue(index);}
    BoolProxy operator[](size_t index)  {return BoolProxy(this, index);}
    bool operator= (bool b) {
        set(0, numUsed, b);
        return b;
    }
    
    BitSet operator~ () const {
        BitSet notB;
        notB.resize(size());
        for (size_t i = 0; i < data.size(); ++i)
            notB.data[i] = ~data[i];
        notB.data.back() &= ((1 << (numUsed & shift_mask)) - 1); // fix the last data item, which might not have all the bits used -> the unused bits must remain 0
        return notB;
    }
    
    bool isZero() const {
        unsigned long long int zero = 0;
        for (auto b : data)
            zero |= b.to_ullong();
        return zero == 0;
    }
    
    void clear() {reset(0, numUsed);}
    
    // multiple bit set
    void set(size_t index_from, size_t index_to, bool val) {
        val ? set(index_from, index_to) : reset(index_from, index_to);
    }
    
    // multiple bit set
    void set(size_t index_from, size_t index_to) {
        if (index_from == index_to) return;
        size_t data_i1 = index_from >> f_shr;
        size_t data_i2 = (index_to-1) >> f_shr;
        unsigned long int data_mask1 = -1l << (index_from & shift_mask);
        unsigned long int data_mask2 = (1l << ((index_to-1) & shift_mask)) - 1;
        if (data_i1 == data_i2) {
            data[data_i1] |= (data_mask1 & data_mask2);
        } else {
            data[data_i1] |= data_mask1;
            for (auto i = data_i1 + 1; i < data_i2; ++i)
                data[i] |= -1l;
            data[data_i2] |= data_mask2;
        }
    }
    
    // multiple bit reset
    void reset(size_t index_from, size_t index_to) {
        if (index_from == index_to) return;
        size_t data_i1 = index_from >> f_shr;
        size_t data_i2 = (index_to-1) >> f_shr;
        unsigned long int data_mask1 = -1l << (index_from & shift_mask);
        unsigned long int data_mask2 = (1l << ((index_to-1) & shift_mask)) - 1;
        if (data_i1 == data_i2) {
            data[data_i1] &= ~(data_mask1 & data_mask2);       
        } else {
            data[data_i1] &= ~data_mask1;
            for (auto i = data_i1 + 1; i < data_i2; ++i)
                data[i] = 0;
            data[data_i2] &= ~data_mask2;
        }
    }
    
    const BitSet& operator&= (const BitSet& other) {
        for (size_t i = 0; i < data.size(); ++i)
            data[i] &= other.data[i];
        return *this;
    }
    
    const BitSet& operator^= (const BitSet& other) {
        //std::cout << data.size() << ", " << other.data.size() << ";\n";
        for (size_t i = 0, im = data.size(); i < im; ++i)
            data[i] ^= other.data[i];
        return *this;
    }
    
    // number of bits set
    size_t count() const {
        return std::accumulate(data.begin(), data.end(), 0, [](size_t a, const b64& b){return a + countOnes(b.to_ulong());});
    }
    
    // position of the next set bit (=1)
    int nextSetBit() const {
        for (size_t i = 0; i < data.size(); ++i) {
            if (data[i] != 0) {
                auto ul = data[i].to_ulong();
                return i*res + countTrailing0(ul);
            }
        }
        return size();
    }
    
    std::string to_string() const {
        std::ostringstream s;
        for (size_t i = 0; i < size(); ++i)
            s << (getValue(i) ? '1' : '0');
        return s.str();
    }

protected:
    bool getValue(size_t index) const {
        if (index >= numUsed) {
        std::cout << "Error in BitSet.getValue:\n" <<
                    "requested " << index << ", holding only " << numUsed << std::endl;
        }
        return data[index >> f_shr][index & shift_mask];}
    void setValue(size_t index, bool b = true) {data[index >> f_shr][index & shift_mask] = b;}
    void copy(const BitSet& other) {
        if (&other != this) {
            resize(other.size());
            std::copy(other.data.cbegin(), other.data.cend(), data.begin());
        }
    }
};

template<class Out>
Out& operator<< (Out& out, const BitSet& b) {
    if (b.size() == 0) {
        out << "[/]";
    } else {
        out << "[l1731:";
        bool comma = false;
        for (size_t i = 0; i < b.size(); ++ i) {
            if (b[i]) {
                if (comma)
                    out << ",";
                comma = true;
                out << i;  
            }
        }
        out << "]";
    }
    return out;
}

#endif // BITSET_H



// #include "MaximumCliqueBase.h"
// #include "SteadyGreedyColorSort.h" // SteadyVectorSet is reused here

/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 */

#ifndef STEADYGREEDYCOLORSORT_H
#define STEADYGREEDYCOLORSORT_H


// #include "MaximumCliqueBase.h"


// Set of vertices is based on std::vector
// every set consists of two sets with different ordering
// used as VertexSetRepresentation in MaximumCliqueProblem
template<class T>
class SteadyVectorSet : public std::vector<std::pair<T,T>> {
public:
    typedef T VertexId;
    typedef std::vector<std::pair<T,T>> VecPairT;
    
    using VecPairT::resize;
    using VecPairT::reserve;
    void add(const T& value) {push_back(std::make_pair(value));}
    T pop() {T temp=this->back(); this->pop_back(); return temp;}
    void remove(T value) {
        if (this->back() == value) this->pop_back();
        else {
            auto temp = std::find(this->begin(), this->end(), value);
            if (temp != this->end())
                this->erase(temp);
        }
    }
    using VecPairT::size;
    //VertexId operator[] (size_t i) const {return VecPairT::operator[](i).second();}
    //VertexId& operator[] (size_t i) {return VecPairT::operator[](i).second();}
    using VecPairT::operator[];
    using VecPairT::clear;
    using VecPairT::pop_back;
    VertexId back() const {return VecPairT::back().second;}
    const std::pair<T,T>& backPair() const {return VecPairT::back();}
};
REGISTER_TEMPLATE1_CLASS_NAME(SteadyVectorSet, "Vector based set of pairs");


template<class Graph>
class SteadyGreedyColorSort {
public:
    typedef typename Graph::VertexSet VertexSet;
    typedef typename VertexSet::VertexId VertexId;
    typedef SteadyVectorSet<VertexId> NumberedSet;
    typedef Graph GraphType;
    
protected:
    std::vector<VertexSet> colorSet;
    Graph* graph;
    
public:
    SteadyGreedyColorSort() : graph(nullptr) {}
    
    void init(GraphType* g) {
        graph = g;
        size_t n = graph->getNumVertices();
        colorSet.resize(n);
        for (auto cc : colorSet)
            cc.reserve(n);
    }
    
    void assignVertexNumber(VertexSet& vs, NumberedSet& ns, size_t i, VertexId vert, VertexId num) {ns[i] = std::make_pair(vert, num);}
    
    bool notEmpty(const NumberedSet& ns) const {return ns.size() > 0;}
    
    VertexId topNumber(const NumberedSet& ns) const {return ns.back();}
    
    VertexId topVertex(const NumberedSet& ns, const VertexSet& vs) const {return ns.backPair().first;}
    
    void popTop(NumberedSet& ns, VertexSet& vs) {
        vs.remove(ns.backPair().first);
        ns.pop_back(); 
    }
    
    // take clique c and candidate vertex set p as input
    // return numbered vertex set np as output
    // difference compared to regular GreedyColorSort: p is here only used as input
    // p is input, np is output (p does not change)
    void numberSort(const VertexSet& c, VertexSet& p, NumberedSet& np, unsigned int maxSize = 0) {
        size_t m = p.size();
        size_t numbers = 0;
            
        for (size_t i = 0; i < m; i++) {
            colorSet[i].clear();
            VertexId v = p[i];
            size_t k = 0;
            
            while (intersectionExists(v, colorSet[k])) 
                k++;
                
            colorSet[k].add(v);
            numbers = std::max(numbers, k+1);
        }
        
        np.resize(m);
        for (size_t k = 0, i = 0; k < numbers; k++) {
            for (size_t j = 0; j < colorSet[k].size(); j++) {
                VertexId v = colorSet[k][j];
                np[i].first = v; 
                //p[i] = v;
                np[i++].second = k+1;
            }
        }
    }

protected:    
    bool intersectionExists(VertexId p, const VertexSet& vertices) const {
        const VertexSet& neighbourhood = graph->adjacencyMatrix[p];
        size_t n = vertices.size();
        for (size_t i = 0; i < n; ++i) {
            if (neighbourhood[vertices[i]] == true)
                return true; 
        }
        return false;
    }
};
REGISTER_TEMPLATE1_CLASS_NAME(SteadyGreedyColorSort, "Steady greedy color sort");

#endif // STEADYGREEDYCOLORSORT_H



// Set of vertices is based on std::bitset
// every set consists of two sets with different ordering
// used as VertexSetRepresentation in MaximumCliqueProblem
class BitstringSet : protected BitSet {
    mutable size_t countCache;
    
public:
    typedef unsigned int VertexId;
    
    ~BitstringSet() {}
    BitstringSet() : countCache(0) {}
    using BitSet::resize;
    void reserve(size_t s) {BitSet::resize(s);};
    void add(VertexId value) {setValue(value, true); ++countCache;}
    void remove(VertexId value) {setValue(value, false); --countCache;}
    size_t recount() const {countCache = count(); return countCache;}
    size_t size() const {return countCache;}
    size_t alsize() const {return BitSet::size();}
    bool operator[](size_t index) const {return getValue(index);}
    BoolProxy operator[](size_t index)  {return BoolProxy(this, index);}
    //using BitSet::operator[];
    void clear() {BitSet::clear(); countCache = 0;};
    using BitSet::nextSetBit;
    using BitSet::operator~;
    
    void operator&=(const BitSet& other) {BitSet::operator&= (other);}
    void operator&=(const BitstringSet& other) {BitSet::operator&= (other);}
    void operator^=(const BitstringSet& other) {BitSet::operator^= (other);}
    
    friend void intersectWithAdjecency (const BitstringSet& v, const BitstringSet& adj, BitstringSet& result) {
        result = v;
        result &= adj;
        result.recount();
    }
    
    void remove(const BitstringSet& values) {
        // this function only works correctly if all the specified values are set in the target bitstring
        *this ^= values;
        countCache -= values.countCache;
    }
    
    bool isIntersectionOf(const BitstringSet& bigSet) {
        return  false;
    }
};
REGISTER_CLASS_NAME(BitstringSet, "bitstring based set");

template<> struct isTypeASet<BitstringSet> { enum {value = true}; };

template<class Out>
Out& operator<< (Out& out, const BitstringSet& b) {
    auto bCopy = b;
    if (b.size() == 0) {
        std::cout << "[/]";
    } else {
        std::cout << "[l1950:";
        bool comma = false;
        while (bCopy.size() > 0) {
            if (comma)
                std::cout << "," << std::flush;
            comma = true;
            auto v = bCopy.nextSetBit();
            std::cout << v;
            bCopy.remove(v);
        }
        std::cout << "]";
    }
    return out;
}

void copy_to_vector(core::Vector_sp result, const BitstringSet& b)
{
    auto bCopy = b;
    if (b.size() != 0) {
        while (bCopy.size() > 0) {
            auto v = bCopy.nextSetBit();
            result->vectorPushExtend(core::make_fixnum(v));
            bCopy.remove(v);
        }
    }
}


template<class Graph>
class BBGreedyColorSort {
public:
    typedef BitstringSet VertexSet;
    typedef typename VertexSet::VertexId VertexId;
    typedef SteadyVectorSet<VertexId> NumberedSet;
    typedef Graph GraphType;
    
protected:
    std::vector<VertexSet> colorSet;
    Graph* graph;
    VertexSet Ubb, Qbb;
    
public:
    BBGreedyColorSort() : graph(nullptr) {}
    
    void init(GraphType* g) {
        graph = g;
        
        size_t n = graph->getNumVertices();
        colorSet.resize(n);
        for (auto& cc : colorSet) {
            cc.reserve(n);
        }
    }
    
    void assignVertexNumber(VertexSet& vs, NumberedSet& ns, size_t i, VertexId vert, VertexId num) {ns[i] = std::make_pair(vert, num);}
    
    bool notEmpty(const NumberedSet& ns) const {return ns.size() > 0;}
    
    VertexId topNumber(const NumberedSet& ns) const {return ns.back();}
    
    VertexId topVertex(const NumberedSet& ns, const VertexSet& vs) const {return ns.backPair().first;}
    
    void popTop(NumberedSet& ns, VertexSet& vs) {auto v = ns.backPair().first; ns.pop_back(); vs.remove(v);}
    
    // take clique c and candidate vertex set p as input
    // return numbered vertex set np as output
    // difference compared to regular GreedyColorSort: p is here only used as input
    void numberSort(const VertexSet& c, const VertexSet& p, NumberedSet& np, unsigned int maxSize = 0) {
        Ubb = p;
        Qbb = Ubb;
        int k = 0;
        size_t i = 0;
        int kMin = (int)maxSize - (int)c.size();
        np.resize(p.size());
        
        while (Ubb.size() > 0) {
            while (Qbb.size() > 0) {
                auto v = Qbb.nextSetBit();
                Qbb &= graph->invAdjacencyMatrix[v];
                Qbb.recount();
                Ubb.remove(v);
                if (k >= kMin) {
                    np[i].first = v;
                    np[i++].second = k+1;
                }
            }
            Qbb = Ubb;
            ++k;
        }
        np.resize(i);
    }

protected:    
    bool intersectionExists(VertexId p, const VertexSet& vertices) const {
        const VertexSet& neighbourhood = graph->adjacencyMatrix[p];
        auto intersection = neighbourhood;
        intersection &= vertices;
        return intersection.size() > 0;
    }
};
REGISTER_TEMPLATE1_CLASS_NAME(BBGreedyColorSort, "greedy color sort on bitstrings");

#endif // BB_GREEDYCOLORSORT_H_INCLUDED



// #include "BB_ColorRSort.h"

/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 */

#ifndef HEADER_38B6D205806EECF7
#define HEADER_38B6D205806EECF7


// #include "BB_GreedyColorSort.h"

template<class Graph>
class BBColorRSort : public BBGreedyColorSort<Graph> {
public:
    typedef BBGreedyColorSort<Graph> ParentT;
    typedef typename ParentT::VertexSet VertexSet;
    typedef typename ParentT::NumberedSet NumberedSet;
    typedef typename ParentT::VertexId VertexId;
    
protected:
    // tell compiler about parent class member variables
    using ParentT::Ubb;
    using ParentT::Qbb;
    VertexSet intersection; // used in resort function
    
public:
    using ParentT::colorSet;
    using ParentT::graph;
    
    void numberSort(const VertexSet& c, const VertexSet& p, NumberedSet& np, unsigned int maxSize = 0) {
        Ubb = p;
        Qbb = Ubb;
        int k = 0;
        size_t i = 0;
        int kMin = (int)maxSize - (int)c.size();
        np.resize(p.size());
        
        while (Ubb.size() > 0) {
            bool first = false;
            colorSet[k].clear();
            Qbb = Ubb;
            while (Qbb.size() > 0) {
                auto v = Qbb.nextSetBit();
                if ((Qbb.size() == 1) && (k >= kMin) && first) {
                    //std::cout << "Recolor " << v << " of Qbb " << Qbb.size() << " " << k << " " << kMin << "\n";
                    if (recolor(v, Qbb, k, kMin)) {
                        Ubb.remove(v);
                        // because of the condition (Qbb.cachedSize() == 1), Qbb only contains v at this point
                        // so jump right out of the loop to save several cycles
                        break;
                        //Qbb.remove(v);
                    } else {
                        first = false;
                    }
                }
                colorSet[k].add(v);
                Qbb &= graph->invAdjacencyMatrix[v];
                Qbb.recount(); 
            }
            Ubb.remove(colorSet[k]);
            ++k;
        }
        
        for (int k1 = std::max(0,kMin); k1 < k; ++k1) {
            for (int j = colorSet[k1].size(); j > 0; --j) {
                auto v = colorSet[k1].nextSetBit();
                colorSet[k1].remove(v);
                np[i].first = v;
                np[i++].second = k1+1;
            }
        }
        np.resize(i);
    }
    
    bool recolor(VertexId v, VertexSet& p, int kv, int kMin) {
        for (int k1 = 0; k1 < kMin - 1; ++k1) {
            intersection = colorSet[k1];
            intersection &= graph->adjacencyMatrix[v];
            intersection.recount();
            if (intersection.size() == 1) {
                auto w = intersection.nextSetBit();
                for (int k2 = k1+1; k2 < kMin; ++k2) {
                    intersection = colorSet[k2];
                    intersection &= graph->adjacencyMatrix[w];
                    intersection.recount();
                    if (intersection.size() == 0) {
                        colorSet[k1].add(v);
                        colorSet[k1].remove(w);
                        colorSet[k2].add(w);
                        return true;
                    }
                }
            } else if (intersection.size() == 0) {
                colorSet[k1].add(v);
                return true;
            }
        }
        return false;
    }
};
REGISTER_TEMPLATE1_CLASS_NAME(BBColorRSort, "greedy color sort with resort on bitstrings");


#endif // header guard 



// #include "McrBB.h"

/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 */

#ifndef MCRBB_H_INCLUDED
#define MCRBB_H_INCLUDED

/**
    This file is a copy of Mcr.h but fixed for BB. TODO: include support for BBs to Mcr
**/

// #include "MaximumCliqueBase.h"


template<class ColorSort>
struct BBMcrSort : public ColorSort {
    using ColorSort::numberSort;
    typedef typename ColorSort::GraphType GraphType;
    typedef typename GraphType::VertexSet VertexSet;
    typedef typename VertexSet::VertexId VertexId;
    typedef typename ColorSort::NumberedSet NumberedSet;
    
    struct Vertex {
        int degree;
        int exDegree;
        int index;
        template<class Ostream>
        friend Ostream& operator<< (Ostream& out, const Vertex& v) {
            out << "<" << v.index << ":" << v.degree << ";" << v.exDegree << ">";
            return out;
        }
    };
    
    using ColorSort::init;
    using ColorSort::assignVertexNumber;
    
    BBMcrSort() {}
        
    void initialSort(VertexSet& c, VertexSet& vertices, NumberedSet& color) {
        TRACE("MCR style initial sort", TRACE_MASK_CLIQUE | TRACE_MASK_INITIAL, 1); 
        size_t n = vertices.size();
        if (n == 0) {
            TRACE("initial sort foud n == 0 and exited", TRACE_MASK_CLIQUE | TRACE_MASK_INITIAL, 1); 
            return;
        }
        auto maxDegree = this->graph->degrees[0];
        std::vector<Vertex> r(n);
        std::vector<VertexId> order(n); // reordering vector
        
        for (size_t i = 0; i < n; ++i) {
            r[i].index = i;
            r[i].degree = this->graph->degrees[r[i].index];
            r[i].exDegree = 0;
            maxDegree = std::max(maxDegree, r[i].degree);
        }
        TRACEVAR(r, TRACE_MASK_INITIAL, 2);
        
        // not sure if the following calculation is correct for ex-deg (not clearly specified in Tomita 2006)
        //  it is possible this should be done on every step of the following while loop, taking only
        //  the neighbourhood of the observed vertex into an account (but probably not)
        for (size_t i = 0; i < n; ++i) 
            for (size_t j = 0; j < n; ++j) 
                if (this->graph->adjacencyMatrix[i][j] == true)
                    r[i].exDegree += r[j].degree;
                    
        // sort by degree (descending), stable mode (respect relative order of the vertices with the same degree)
        std::stable_sort(r.begin(), r.end(), [](const Vertex& a, const Vertex& b) {return (a.degree > b.degree);});
        TRACEVAR(r, TRACE_MASK_INITIAL, 2); 
                
        // index in vertices
        size_t vi = n-1;
        
        size_t rMinIndex = r.size()-1;
        while (rMinIndex > 0) {
            // locate vertices with min degree
            // set of vertices "Rmin" is implemented as a subarray from index rMinIndex to the end of the set r
            int minDeg = r.back().degree;
            rMinIndex = r.size()-1;
            while ((rMinIndex > 0) && (r[rMinIndex-1].degree == minDeg)) 
                --rMinIndex;
            
            if (rMinIndex == 0)
                break;
                
            // if "Rmin" contains more than 1 element
            if (rMinIndex < r.size()-1) {
                // sort by ex-deg (descending - max is first)
                std::stable_sort(r.begin()+rMinIndex, r.end(), [](const Vertex& a, const Vertex& b){return (a.exDegree > b.exDegree);});
            }
            
            // vertex with min ex-deg in rMin goes into the ordered set of vertices (filled from the back towards the front)
            Vertex& p = r.back();
            order[vi] = p.index;
            --vi;
            
            // decrease the degree of remaining vertices that are adjacent to p
            auto& adjacent = this->graph->adjacencyMatrix[p.index];
            r.pop_back();
            size_t rs = r.size();
            for (size_t i = rs; i > 0; --i) {
                if (adjacent[r[i-1].index] == true) {
                    auto rid = --r[i-1].degree;
                    // sort the modified vertex immediately
                    for (size_t j = i; (j < rs) && (rid  < r[j].degree); ++j) 
                        swap(r[j-1], r[j]);
                }
            }
        }
        TRACE("degree of leftover vortices MCR initial sort:", TRACE_MASK_INITIAL, 2); 
        TRACEVAR(r[0].degree, TRACE_MASK_INITIAL, 2);
        TRACE("after calculation of ex-degree:", TRACE_MASK_INITIAL, 2); 
        TRACEVAR(r, TRACE_MASK_INITIAL, 2); 
        TRACEVAR(c, TRACE_MASK_INITIAL, 2); 
//      std::cout << "DEBUG: iterated to regular subgraph of degree " << r.front().degree << " and size " << vi << "\n";
        
        // all the vertices in r have the same degree (regular subgraph) → perform ordinary number sort (color sort) 
        c.clear();
        c.reserve(n);
        VertexSet dummySet;
        dummySet.reserve(n);
        
        for (size_t i = 0; i < r.size(); ++i) {
            c.add(r[i].index);
            dummySet.add(r[i].index);
        }
        
        numberSort(dummySet, c, color, 0);
        
        // calculate maximum number (color), then fill in the rest of the order vector (clearing colorC and c in the process)
        size_t m = color.back();
        size_t mmax = r.size() + maxDegree - m;
        std::vector<VertexId> storedColor(r.size());
        for (size_t i = r.size(); i > 0; --i) {
            order[i-1] = this->topVertex(color, c); 
            storedColor[i-1] = this->topNumber(color); 
            this->popTop(color, c);
        }

        // if the degree of vertices in r (which all have the same degree) equals r.size() - 1, then r is a clique
        if (r.size() == r[0].degree+1) {
        	// fill the clique with the new vertex numbers (these are actually numbers [0..r.size()-1]), the ones that will be set by the orderVertices function also for the rest of the graph
        	//c.resize(r.size());
        	c.reserve(n);
        	for (size_t i = 0; i < r.size(); ++i) {
        		c.add(i);
        	}
        	TRACE("initial clique found:", TRACE_MASK_INITIAL, 2); 
        	TRACEVAR(c, TRACE_MASK_INITIAL, 2); 
        }
        
        // reorder vertices in graph
        this->graph->orderVertices(order);
        
        // number
        color.resize(n); 
        vertices.clear();
        // first few colors remain as they are, the rest are filled in with up to mmax
        for (size_t i = 0; i < n; ++i) {
            if (i < r.size()) {
                assignVertexNumber(vertices, color, i, i, storedColor[i]);
            } else if (i < mmax) {
                ++m;
                assignVertexNumber(vertices, color, i, i, m);
            } else {
                assignVertexNumber(vertices, color, i, i, maxDegree + 1);
            }
            vertices.add(i);
            //std::cout << i << ": " << this->graph->mapping[color[i].first] << "," << color[i].second << "  ";
        }
        
        TRACE("initial coloring of vertices:", TRACE_MASK_INITIAL, 2); 
        TRACEVAR(color, TRACE_MASK_INITIAL, 2); 
        TRACEVAR(c, TRACE_MASK_INITIAL, 2); 
    }
};
REGISTER_TEMPLATE_EXT_CLASS_NAME(BBMcrSort, "MCR sort for Bitstrings");


#endif // MCRBB_H_INCLUDED


// #include "BitSet.h"



//////////////////////////////////////////////////////////////////////


#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <algorithm>

namespace chem {

    CL_LISPIFY_NAME(make-dimacs);
    CL_DEF_CLASS_METHOD
    Dimacs_sp Dimacs_O::make(size_t numVertices, size_t numEdges ) {
      auto dimacs = gctools::GC<Dimacs_O>::allocate(numVertices,numEdges);
        return dimacs;
    }

    Dimacs_O::Dimacs_O() : numVertices(0), adjacencyMatrixSizeLimit(1000000000 /* 1GB */) {
    }

    Dimacs_O::Dimacs_O(size_t num_vertices, size_t num_edges ) : adjacencyMatrixSizeLimit(1000000000 /* 1GB */) {
        this->numVertices = num_vertices;
        this->edges.reserve(num_edges);
        this->degrees.resize(num_vertices,0);
    }
    
    Dimacs_O::~Dimacs_O() {
    }

    bool Dimacs_O::load(const char* fname) {
        std::ifstream file(fname);
        std::string line;

        if (!file)
            return false;

        while (file) {
            std::getline(file, line);
            if (line.size() > 0) {
                std::istringstream ss(line);
                ss.ignore(2);
                if (line[0] == 'p') { // problem line
                    std::string word;
                    ss >> word;
                    if (word != "edge")
                        return false;
                    size_t nEdges = 0;
                    ss >> numVertices >> nEdges;
                    edges.reserve(nEdges);
                    degrees.resize(numVertices, 0);
                } else if (line[0] == 'e') { // edge line
                    int v1, v2;
                    ss >> v1 >> v2;
                    edges.push_back(std::make_pair(v1-1, v2-1));
                    if (numVertices < (unsigned int)std::max(v1, v2)) {
                        numVertices = std::max(v1, v2);
                        degrees.resize(numVertices, 0);
                    }
                    degrees[v1-1]++;
                    degrees[v2-1]++;
                }
            }
        }
        return true;
    }

    CL_LISPIFY_NAME(dimacs-add-edge);
    CL_DEFMETHOD
    void Dimacs_O::add_edge(size_t v1, size_t v2) {
        edges.push_back(std::make_pair(v1, v2));
        if (numVertices < (unsigned int)std::max(v1+1, v2+1)) {
            numVertices = std::max(v1+1, v2+1);
            degrees.resize(numVertices, 0);
        }
        degrees[v1]++;
        degrees[v2]++;
    }

gctools::Vec0<gctools::Vec0<char> > Dimacs_O::getAdjacencyMatrix() const {
        if (adjacencyMatrixSizeLimit < getNumVertices() * getNumVertices())
            throw "Cannot create adjacency matrix because the number of vertices is to large";
        gctools::Vec0<gctools::Vec0<char> > matrix;
        matrix.resize(getNumVertices());
        for (auto& v : matrix) {
            v.resize(getNumVertices(), 0);
        }
        
        for (unsigned int ei = 0; ei < edges.size(); ++ei) {
            matrix[edges[ei].first][edges[ei].second] = 1;
            matrix[edges[ei].second][edges[ei].first] = 1;
        }
        return matrix;
    }


    void Dimacs_O::calculateGraphStats(int& maxDegree, int& minDegree, std::vector<float>& degreeHistogram) {
        auto degCopy = degrees;
        std::sort(&*(degCopy.begin()), &*(degCopy.end()));
        size_t n = degCopy.size();
        size_t m = degreeHistogram.size();
        int cnt = 0;
        size_t di = n-1;
        for (size_t i = 0; i < m; ++i) {
            float bound = (m-i-1)*(float)n / m;
            while (degCopy[di] > bound && di > 0) {
                --di;
                ++cnt;
            }
            degreeHistogram[i] = cnt / (float)n;
            cnt = 0;
        }
        maxDegree = degCopy.front();
        minDegree = degCopy.back();
    }


    ////////////////////////////////////////////////////////////

DOCGROUP(cando);
CL_DEFUN
    core::T_sp chem__find_maximum_clique_search(Dimacs_sp dimacsGraph, int numThreads, int numJobs) {
        if (numThreads>0) {
            printf("%s:%d Parallel maximum-clique-search is broken - it calls back() with empty vectors\n", __FILE__, __LINE__ );
            numThreads = 0;
        }
        bool coloredOut = false;
#ifdef __linux
        coloredOut = true;
#endif
        Graph<BitstringSet> graphB;
        graphB.init(dimacsGraph->getAdjacencyMatrix(), dimacsGraph->getDegrees());
#if 0
        std::cout << dimacsGraph->getNumVertices() << " vertices " << dimacsGraph->getNumEdges() << " edges " 
                  << dimacsGraph.getNumEdges()*2.0/(dimacsGraph.getNumVertices()*(dimacsGraph.getNumVertices()-1)) << "\n";
#endif
        core::Vector_sp clique = core::core__make_vector(_lisp->_true(), 256, true, core::make_fixnum(0));
        if (numThreads == 0) {
            MaximumCliqueProblem<
                int,                        // vertex ID
                BitstringSet,               // vertex set
                Graph<BitstringSet>,        // graph
                BBMcrSort<BBGreedyColorSort<Graph<BitstringSet>>>        // initial & color sort
                > problem(graphB);
            problem.search();
            problem.outputStatistics(clique,coloredOut);
        } else {
            ParallelMaximumCliqueProblem<
                int,                        // vertex ID
                BitstringSet,               // vertex set
                Graph<BitstringSet>,        // graph      
                BBGreedyColorSort<Graph<BitstringSet>>,        // color sort
                BBMcrSort        // initial sort
                > problemP1(graphB);
            std::vector<int> affinities;
            problemP1.search(numThreads, numJobs, affinities);
            problemP1.outputStatistics(clique,coloredOut);
            // std::cout << "Thread efficiency = " << std::setprecision(3) << problemP1.workerEfficiency() << "\n\n";
        }
        return clique;
    }

#if 0
    int main(int argc, char** argv) {
        try {
            std::cout << "Trace set to: " << (TRACE_MASK) << " (" << (TRACE_LEVEL) << ")\n";
            std::vector<int> bindProcessors;
            int numThreads = 0, numJobs = 1;
            if (argc > 2) {
                numThreads = atoi(argv[2]);
                if (numThreads < 0)
                    numThreads = 0;
                if (argc > 3) {
                    numJobs = atoi(argv[3]);
                    if (numJobs < 0)
                        numJobs = 0;
                    if (argc > 4) {
                        for (int argi = 4; argi < argc; ++argi)
                            bindProcessors.push_back(atoi(argv[argi]));
                    }
                }
            }
            if (numJobs < 1)
                std::cout << "Warning: Since the number of jobs is set to " << numJobs << ", the algorithm will run no faster than the sequential algorithm\n";
            pBbmcTest(argc > 1 ? argv[1] : nullptr, numThreads, numJobs, bindProcessors);
        } catch (const char* e) {
            std::cout << "Terminated due to exception: " << e << std::endl;
        } catch (...) {
            std::cout << "Terminated due to unknown exception: " << std::endl;
        }
        return 0;
    }
#endif
};
