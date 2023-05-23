
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
    Alloc Allocate;
    unsigned numSets;

    vector <unsigned> validBit;
    vector <unsigned> DirtyBit;
    vector <unsigned> Tag;
    vector <unsigned> Counter;

    unsigned calcSet(unsigned address);
    unsigned calcTag(unsigned address);



public:
    Stat stats;
    CacheTable() = default;
    CacheTable(unsigned BSize, unsigned numEntries, unsigned Ways, unsigned Cycles, Alloc WrAlloc): numEntries(numEntries), Ways(Ways),
                BSize(BSize), Allocate(WrAlloc){
        //initiate vectors
        for (unsigned i = 0; i < numEntries; ++i) {
            validBit.push_back(0);
            Tag.push_back(0);
        }

        numSets = numEntries / Ways;

        for (unsigned i = 0; i < numSets; ++i) {
            for (unsigned j = 0; j < Ways; ++j) {
                Counter[i+j] = j;
            }
        }


        //initiate stats
        stats.numOfMiss=0;
        stats.accTime=Cycles;
        stats.numOfAcc=0;
    }

    unsigned searchInCache(unsigned address);
    void LRUupdate(unsigned way, unsigned set);
    unsigned LRUfindMin(unsigned address, unsigned * victimSet, unsigned * victimTag);

    bool read(unsigned address);
    bool write(unsigned address);
    void getBlockFromMemory(unsigned address);


};


/**********************************************
 * CacheTable Private Methods
 **********************************************/


unsigned CacheTable::calcSet(unsigned address)
{
    unsigned mask = (1 << (unsigned)log2(Ways)) - 1;
    mask = mask << (unsigned)log2(BSize);
    return (mask&address) >> (unsigned)log2(BSize);
}


unsigned CacheTable::calcTag(unsigned address)
{
    return address >> ((unsigned)log2(Ways) + (unsigned)log2(BSize));
}

/***
 * search for the block in cache.
 * return value - the index of the block.
*/
unsigned CacheTable::searchInCache(unsigned address)
{
    unsigned tag = calcTag(address);
    unsigned set = calcSet(address);
    unsigned index = set;

    for (unsigned i=0; i<Ways; i++)
    {
        if((Tag[index] == tag) && validBit[index] == 1)
        {
            LRUupdate(i, set);
            return index;
        }
        index += numSets;
    }
    return -1;
}

/***
 * update the Count vector when access to cache.
 * the algorithm is based on Lecture's notes.
 */
void CacheTable::LRUupdate(unsigned way, unsigned set){
    unsigned x = Counter[way+set];
    Counter[way+set] = Ways-1;
    for (unsigned i = 0; i < Ways; ++i) {
        if (i != way && Counter[i+set] > x) Counter[i+set]--;
    }
}


unsigned CacheTable::LRUfindMin(unsigned address, unsigned * victimSet, unsigned * victimTag)
{
    unsigned set = calcSet(address);
    *victimSet = set;
    *victimTag = calcTag(address);
    unsigned minVal = Counter[set];
    unsigned minEntry = set;

    for (unsigned i=0; i<Ways; i++)
    {
        if(validBit[set] == 1 && Counter[set] < minVal)
        {
            minEntry = set;
        }
        set += numSets;
    }
    return minEntry;
}


/**********************************************
 * CacheTable Public Methods
 **********************************************/


bool CacheTable::read(unsigned address)
{
    stats.numOfAcc++;
    unsigned set = searchInCache(address);
    if (set != -1)
    {
        return true;
    }
    stats.numOfMiss++;

    return false;
}


bool CacheTable::write(unsigned address)
{
    stats.numOfAcc++;
    unsigned set = searchInCache(address);
    if (set != -1)
    {
        DirtyBit[set] = 1;
        return true;
    }
    stats.numOfMiss++;

    return false;
}



/**********************************************
 * Cache Class
 **********************************************/


class Cache {
private:
    unsigned MemCyc;
    unsigned BSize;
    Alloc Allocate;
    CacheTable L1;
    CacheTable L2;
    unsigned numMemAcc = 0;

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

    void writeToCache(unsigned address);
    void readFromCache(unsigned address);
    void bringFromDisk(unsigned address);

    void CacheGetStats(Stat * statsL1, Stat *statsL2, Stat *statsMem);
};


/**********************************************
 * Cache Methods
 **********************************************/


void Cache::Cache_initiate(unsigned MemCyc, unsigned BSize, unsigned L1Size, unsigned L2Size, unsigned L1Cyc,unsigned L2Cyc,
                           unsigned L1Assoc, unsigned L2Assoc, unsigned WrAlloc){
    this->MemCyc = MemCyc;
    this->BSize = BSize;
    this->Allocate = static_cast<Alloc>(WrAlloc);

    //calculate the number of entries in each cache
    unsigned num_blocks1 = L1Size/BSize;
    unsigned num_blocks2 = L2Size/BSize;

    this->L1 = CacheTable(BSize, num_blocks1, L1Cyc, L1Assoc, Allocate);
    this->L2 = CacheTable(BSize, num_blocks2, L2Cyc, L2Assoc, Allocate);
};

/*** tries to read from caches L1 and L2.
 * if block not in cache, bring to cache from disk.
***/
void Cache::readFromCache(unsigned address)
{
    if (!L1.read(address))
    {
        //L1 Miss, L2 Miss
        if(!L2.read(address)){
            numMemAcc++;
                bringFromDisk(address);
        }

        //L1 Miss, L2 Hit
        else{
            bringFromLower(address);
        }
    }

}


void Cache::writeToCache(unsigned address)
{
    if (!L1.write(address))
    {
        if (!L2.write(address)){
            numMemAcc++;
            if(Allocate == WRITE_ALLOCATE)
                bringFromDisk(address);
        }
    }
}

/*** bring from disk -
 * if block not in cache, bring to cache from disk.
***/
void Cache::bringFromDisk(unsigned address) {
    //find victim in L2
    unsigned victimSet, victimTag;
    unsigned victimIndex = L2.LRUfindMin(address, &victimSet, &victimTag);




    //Snoop
    unsigned snoopVictimSet, snoopVictimTag;
    bool is_found = false;
    unsigned snoopVictimIndex = L1.Snoop(victimSet, victimTag, &snoopVictimSet, &snoopVictimTag, &is_found);

    if (is_found){
        // if dirty: send write request to L2
        // invalidate block in L1
        if(L1.getDirty(snoopVictimIndex)){
            L2.setDirty(victimIndex, true);
            L1.Invalidate(snoopVictimIndex);
        }
    }

    if (L2.getDirty(victimIndex)){
        //update memory. No need to update numMemoryAcc
        L2.Invalidate(victimIndex);
    }


    //Bring Block to cache
    L2.

}


void Cache::CacheGetStats(Stat * statsL1, Stat *statsL2, Stat *statsMem){
    //todo: add L1 and L2

    //memory
    statsMem->numOfAcc = numMemAcc;
    statsMem->accTime = MemCyc;
}


/**********************************************
 * Main Functions
 **********************************************/


int Cache_init(unsigned MemCyc, unsigned BSize, unsigned L1Size, unsigned L2Size, unsigned L1Cyc,unsigned L2Cyc,
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
    Stat memory;
    cache.CacheGetStats(&statsL1, &statsL2, &memory);

    *L1MissRate = statsL1.numOfMiss/statsL1.numOfAcc;
    *L2MissRate = statsL2.numOfMiss/statsL2.numOfAcc;
    *avgAccTime =  (statsL1.accTime*statsL1.numOfAcc +
                    statsL2.accTime*statsL2.numOfAcc +
                    memory.accTime*memory.numOfAcc) / (statsL1.numOfAcc);
}
