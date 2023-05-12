
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <math.h>

using std::FILE;
using std::string;
using std::cout;
using std::endl;
using std::cerr;
using std::ifstream;
using std::stringstream;
using std::vector;
using std::pair;

enum Alloc {
    NO_WRITE_ALLOCATE, WRITE_ALLOCATE
};

struct Stat {
    double MissRate;
    double accTime;
    double numOfAcc;
};

class Cache {
private:
    unsigned MemCyc;
    unsigned BSize;
    unsigned CacheSize;
    unsigned Ways;
    unsigned Cycles;
    unsigned Allocate;

    vector <unsigned> validBit;
    vector <unsigned> Tag;
    vector <pair<unsigned, unsigned>> BlockRange;

    Stat stats;

    int searchInCache(unsigned tag);

public:
    Cache(unsigned MemCyc, unsigned BSize, unsigned CacheSize, unsigned Associative, unsigned Cycles, unsigned Allocate):
    MemCyc(MemCyc), BSize(BSize), CacheSize(CacheSize), Ways(Associative), Cycles(Cycles), Allocate(Allocate){
        //initiate vectors
        for (int i = 0; i < pow(BSize, 2); ++i) {
            validBit.push_back(0);
            Tag.push_back(0);
            pair<unsigned, unsigned> p = {0, 0};
            BlockRange.push_back(p);
        }
        stats.MissRate =0;
        stats.accTime=0;
        stats.numOfAcc=0;
    };
    bool writeToCache(unsigned tag, unsigned data);
    bool readFromCache(unsigned tag, unsigned data);

};



