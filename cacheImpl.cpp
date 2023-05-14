
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

class CacheTable{
    unsigned numEntries;
    unsigned Ways;
    unsigned Cycles;

    vector <unsigned> validBit;
    vector <unsigned> DirtyBit;
    vector <unsigned> Tag;
    vector <pair<unsigned, unsigned>> BlockRange;

    Stat stats;

public:
    CacheTable() = default;
    CacheTable(unsigned BSize, unsigned numEntries, unsigned Ways, unsigned Cycles): numEntries(numEntries), Ways(Ways), Cycles(Cycles){
        //initiate vectors
        for (int i = 0; i < numEntries; ++i) {
            validBit.push_back(0);
            Tag.push_back(0);
            pair<unsigned, unsigned> p = {0, 0};
            BlockRange.push_back(p);
        }

        //todo: calc tag

        //initiate stats
        stats.MissRate =0;
        stats.accTime=0;
        stats.numOfAcc=0;
    }
};

class Cache {
private:
    unsigned MemCyc;
    unsigned BSize;
    unsigned Allocate;


    CacheTable L1;
    CacheTable L2;

    int searchInCache(unsigned tag);

public:
    Cache() = default;

    void Cache_initiate(unsigned MemCyc, unsigned BSize, unsigned L1Size, unsigned L2Size, unsigned L1Cyc,unsigned L2Cyc,
                               unsigned L1Assoc, unsigned L2Assoc, unsigned WrAlloc);

    Cache(Cache const &) = delete; // disable copy ctor
    void operator=(Cache const &) = delete; // disable = operator

    ~Cache() = default;

    static Cache &getInstance() // make singleton
    {
        static Cache instance; // Guaranteed to be destroyed.
        // Instantiated on first use.
        return instance;
    }


    bool writeToCache(unsigned address);
    bool readFromCache(unsigned address);

    void CacheGetStats(Stat * statsL1, Stat *statsL2);

};



void Cache::Cache_initiate(unsigned MemCyc, unsigned BSize, unsigned L1Size, unsigned L2Size, unsigned L1Cyc,unsigned L2Cyc,
                           unsigned L1Assoc, unsigned L2Assoc, unsigned WrAlloc){
    this->MemCyc = MemCyc;
    this->BSize = BSize;
    this->Allocate = WrAlloc;

    //calculate the number of entries in each cache
    unsigned num_blocks1 = L1Size/BSize;
    unsigned num_blocks2 = L2Size/BSize;



    this->L1 = CacheTable(BSize, num_blocks1, L1Cyc, L1Assoc);
    this->L2 = CacheTable(BSize, num_blocks2, L2Cyc, L2Assoc);



};

/**********************************************
 * Main Functions
 **********************************************/

int  Cache_init(unsigned MemCyc, unsigned BSize, unsigned L1Size, unsigned L2Size, unsigned L1Cyc,unsigned L2Cyc,
                unsigned L1Assoc, unsigned L2Assoc, unsigned WrAlloc){
    Cache &cache = Cache::getInstance();
    if(&cache == nullptr) return -1;
    cache.Cache_initiate(MemCyc, BSize, L1Size, L2Size, L1Cyc, L2Cyc, L1Assoc, L2Assoc, WrAlloc);

}


void Cache_Access(char operation, unsigned address){
    Cache &cache = Cache::getInstance();
    switch (operation) {
        case 'r':
            cache.readFromCache(address);
            break;
        case 'w':
            cache.writeToCache(address);
            break;
    }
}


void Cache_GetStats(double * L1MissRate, double * L2MissRate, double * avgAccTime){
    Cache &cache = Cache::getInstance();
    Stat statsL1;
    Stat statsL2;

    cache.CacheGetStats(&statsL1, &statsL2);

    *L1MissRate = statsL1.MissRate;
    *L2MissRate = statsL2.MissRate;
    *avgAccTime = (statsL1.accTime + statsL2.accTime)/ (statsL1.numOfAcc+statsL2.numOfAcc);

}
