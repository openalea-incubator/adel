/*----------------------------------------------------------------------------
    Chrono.C	simple Chrono class (derived from Rudi's Chrono class)

    by Andreas Hohmann, ZIB, Hohmann@sc.ZIB-Berlin.DE
    
    last change:	Sep 10, 1994
    Compiler: 		CC/g++/Bcc
-----------------------------------------------------------------------------*/



#include <fstream>
using namespace std ;

#include <chrono.h>

#ifndef CLOCKS_PER_SEC              // should be defined in ANSI C but ...
#define CLOCKS_PER_SEC 1000000
#endif

//----------------------------------------------------------------------------
//    ChronoData: helper class for Chrono
//----------------------------------------------------------------------------

/*  UNIX version of the ChronoData class, not used any more
 

#include <sys/times.h>
#include <time.h>

class ChronoData {
    friend class Chrono;
private:
    ChronoData() {}
    void Start();
    void Stop();
    double Seconds() const;
    
    struct  tms buffer;
    clock_t tO, t;
    time_t  abstO, abst;
};

void ChronoData::Start() {
    times(&buffer);
    tO = buffer.tms_utime;
    time(&abstO);
}

void ChronoData::Stop() {
    times(&buffer);
    t = buffer.tms_utime;
    time(&abst);
    abst = abst - abstO;
}

double ChronoData::Seconds() const {
    return  floatTime = (t-tO)/60. ;
}

*/

//  ANSI C version of the ChronoData class

#include <time.h>

#ifndef _CLOCK_T
#define _CLOCK_T
typedef long int             clock_t;
#endif


class ChronoData {
    friend class Chrono;
private:
    ChronoData() {}
    void Start();
    void Stop();
    double Seconds();
    clock_t tO;
    clock_t t;
};

void ChronoData::Start() {
    tO = clock();
}

void ChronoData::Stop() {
    t = clock();
}

double ChronoData::Seconds()  {
    return (double)(t-tO) / CLOCKS_PER_SEC;
}

//----------------------------------------------------------------------------
//    Chrono
//----------------------------------------------------------------------------

Chrono::Chrono() {
    data = new ChronoData;
    Start();
}

Chrono::~Chrono() {
    delete data;
}

void Chrono::Start() {
    data->Start();
}

void Chrono::Stop() {
    data->Stop();
}
double Chrono::Seconds(){
     return data->Seconds();
}
  
ostream& Chrono::PrintOn(ostream &s) {
    double tmp = data->Seconds();
    s << tmp << " cpu seconds (";
    s << int(tmp)/60 << "' " << int(tmp)%60 << "'') ";
    return s;
}

ferrlog& Chrono::PrintOn(ferrlog &s) {
    double tmp = data->Seconds();
    s << tmp << " cpu seconds (";
    s << int(tmp)/60 << "' " << int(tmp)%60 << "'') ";
    return s;
}

ostream & operator <<(ostream& out ,Chrono& uhr) {
  return uhr.PrintOn(out);
}

// HA 01 2001
ferrlog & operator <<(ferrlog& out ,Chrono& uhr) {
  return uhr.PrintOn(out);
}
