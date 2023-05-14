
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
    double numOfMiss;
    double numOfAcc;
    double accTime;
};


/**********************************************
 * CacheTable Class
 **********************************************/


class CacheTable{
    unsigned numEntries;
    unsigned Ways;
    unsigned BSize;

    vector <unsigned> validBit;
    vector <unsigned> DirtyBit;
    vector <unsigned> Tag;
    vector <unsigned> Counter;
    //vector <pair<unsigned, unsigned>> BlockRange;

    unsigned calcSet(unsigned address);
    unsigned calcTag(unsigned address);
    int searchInCache(unsigned address);
    int findLRU(unsigned address);


public:
    Stat stats;
    CacheTable() = default;
    CacheTable(unsigned BSize, unsigned numEntries, unsigned Ways, unsigned Cycles): numEntries(numEntries), Ways(Ways),
                BSize(BSize){
        //initiate vectors
        for (int i = 0; i < numEntries; ++i) {
            validBit.push_back(0);
            Tag.push_back(0);
            Counter.push_back(0);
            pair<unsigned, unsigned> p = {0, 0};
            //BlockRange.push_back(p);
        }

        //initiate stats
        stats.numOfMiss=0;
        stats.accTime=Cycles;
        stats.numOfAcc=0;
    }

    bool read(unsigned address);
    bool write(unsigned address);
};


/**********************************************
 * CacheTable Private Methods
 **********************************************/


unsigned CacheTable::calcSet(unsigned address)
{
    int mask = (1 << (int)log2(Ways)) - 1;
    mask = mask << (int)log2(BSize);
    return (mask&address) >> (int)log2(BSize);
}


unsigned CacheTable::calcTag(unsigned address)
{
    return address >> ((int)log2(Ways) + (int)log2(BSize));
}


int CacheTable::searchInCache(unsigned address)
{
    int tag = calcTag(address);
    int set = calcSet(address);
    for (int i=0; i<Ways; i++)
    {
        if(Tag[set] == tag && validBit[set] == 1)
        {
            Counter[set]++;
            return set;
        }
        set += numEntries/Ways;
    }
    return -1;
}


int CacheTable::findLRU(unsigned address)
{
    int set = calcSet(address);
    int min = Counter[set];
    int minSet = set;
    for (int i=0; i<Ways; i++)
    {
        if(Counter[set] < min)
        {
           minSet = set;
        }
        set += numEntries/Ways;
    }
    return minSet;
}


/**********************************************
 * CacheTable Public Methods
 **********************************************/


bool CacheTable::read(unsigned address)
{
    stats.numOfAcc++;
    if (searchInCache(address) != -1)
    {
        return true;
    }
    stats.numOfMiss++;
    return false;
}


bool CacheTable::write(unsigned address)
{
    stats.numOfAcc++;
    int set = searchInCache(address);
    if (set != -1)
    {
        DirtyBit[set] = 1;
        return true;
    }
    stats.numOfMiss++;
}


/**********************************************
 * Cache Class
 **********************************************/


class Cache {
private:
    unsigned MemCyc;
    unsigned BSize;
    unsigned Allocate;
    CacheTable L1;
    CacheTable L2;

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


/**********************************************
 * Cache Methods
 **********************************************/


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


bool Cache::readFromCache(unsigned address)
{
    if (!L1.read(address))
    {
        L2.read(address);
    }
}
//TODO:


bool Cache::writeToCache(unsigned address)
{
    if (!L1.write(address))
    {
        L2.write(address);
    }
}


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

    *L1MissRate = statsL1.numOfMiss/statsL1.numOfAcc;
    *L2MissRate = statsL2.numOfMiss/statsL2.numOfAcc;
    *avgAccTime = (statsL1.accTime*statsL1.numOfAcc + statsL2.accTime*statsL2.numOfAcc) / (statsL1.numOfAcc+statsL2.numOfAcc);
}
