
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


class Entry {
public:
    unsigned address;
    unsigned validBit;
    unsigned DirtyBit;

    Entry() {
        address = 0;
        validBit = 0;
        DirtyBit = 0;
    }

    ~Entry() = default;
};


/**********************************************
 * CacheTable Class
 **********************************************/


class CacheTable {
    unsigned numEntries;
    unsigned Ways;
    unsigned BSize;
    Alloc Allocate;
    unsigned numSets;

    vector<Entry *> Entries;
    vector<unsigned> Tag;
    vector<unsigned> Counter;




public:
    Stat stats;

    CacheTable() = default;

    CacheTable(unsigned BSize, unsigned numEntries, unsigned Ways, unsigned Cycles, Alloc WrAlloc) : numEntries(
            numEntries), Ways(Ways), BSize(BSize), Allocate(WrAlloc) {
        //initiate vectors
        for (unsigned i = 0; i < numEntries; ++i) {
            Entry *e = new Entry();
            Entries.push_back(e);
            Tag.push_back(0);
            Counter.push_back(0);
        }
        this-> Ways = pow(2, Ways) ;
        numSets = numEntries / this->Ways;

        for (unsigned i = 0; i < numSets; ++i) {
            for (unsigned j = 0; j < this->Ways; ++j) {
                Counter[i + j] = j;
            }
        }

        //initiate stats
        stats.numOfMiss = 0;
        stats.accTime = Cycles;
        stats.numOfAcc = 0;
    }

    ~CacheTable() = default;

    void erase();

    unsigned calcSet(unsigned address);

    unsigned calcTag(unsigned address);

    int searchInCache(unsigned address, bool is_update);

    void CountUpdate(unsigned way, unsigned set);

    unsigned LRUfindMin(unsigned address, unsigned *victimTag);

    bool read(unsigned address, bool is_update);

    bool write(unsigned address, bool is_update);

    void bringFromLower(unsigned address, Entry *oldBlock);

    void Evict(unsigned address);
};


/**********************************************
 * CacheTable Methods
 **********************************************/


void CacheTable::erase(){
    if (!Entries.empty()){
        for (Entry * entry : Entries) {
            delete entry;
        }
    }
}



bool CacheTable::read(unsigned address, bool is_update) {
    int index = searchInCache(address, is_update);
    if (index != -1) {
        return true;
    }
    return false;
}


bool CacheTable::write(unsigned address, bool is_update) {
    int index = searchInCache(address, is_update);
    if (index != -1) {
        Entries[index]->DirtyBit = 1;
        return true;
    }
    return false;
}


unsigned CacheTable::calcSet(unsigned address) {
    if (numSets == 1) return 0;
    unsigned mask = (1 << (unsigned) log2(numSets)) - 1;
    mask = mask << (unsigned) log2(BSize);
    return ((mask & address) >> (unsigned) log2(BSize));
}


unsigned CacheTable::calcTag(unsigned address) {
    return address >> ((unsigned) log2(numSets) + (unsigned) log2(BSize));
}


int CacheTable::searchInCache(unsigned address, bool is_update) {
    unsigned tag = calcTag(address);
    unsigned set = calcSet(address);
    unsigned index = set;
    for (unsigned i = 0; i < Ways; i++) {
        if ((Tag[index] == tag) && Entries[index]->validBit == 1) {
            if(is_update) {
                CountUpdate(i, set);
            }
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
void CacheTable::CountUpdate(unsigned way, unsigned set) {
    unsigned x = Counter[way*numSets + set]; //TODO: added *numSets
    Counter[way*numSets + set] = Ways - 1;
    for (unsigned i = 0; i < Ways; ++i) {
        if (i != way && Counter[i*numSets + set] > x) Counter[i*numSets + set]--;
    }
}


unsigned CacheTable::LRUfindMin(unsigned address, unsigned *victimWay) {
    unsigned set = calcSet(address);
    unsigned minVal = Counter[set];
    unsigned minEntry = set;

    for (unsigned i = 0; i < Ways; i++) {
        if (Entries[set]->validBit == 1 && Counter[set] < minVal) {
            minEntry = set;
            minVal = Counter[set];
            *victimWay = i;
        }
        set += numSets;
    }
    return minEntry;
}

// TODO: go over with Shai
void CacheTable::bringFromLower(unsigned address, Entry *oldBlock) {
    unsigned tag = calcTag(address);
    unsigned set = calcSet(address);
    unsigned index = set;
    bool is_empty = false;

    unsigned entryWay = 0;

    //try to find empty entry
    for (unsigned i = 0; i < Ways; i++) {
        if (Entries[index]->validBit == 0) { //TODO: Tag[index] == tag && (was part of if)
            is_empty = true;
            entryWay = i;
            break;
        }
        index += numSets;
    }

    //no empty block - find victim
    if (!is_empty) {
        index = LRUfindMin(address, &entryWay);
    }

    //copy old block
    oldBlock->DirtyBit = Entries[index]->DirtyBit;
    oldBlock->validBit = Entries[index]->validBit;
    oldBlock->address = Entries[index]->address;

    //update to new entry
    CountUpdate(entryWay, set);

    Entries[index]->validBit = 1;
    Entries[index]->address = address;
    Entries[index]->DirtyBit = 0;
    Tag[index] = tag;
}


void CacheTable::Evict(unsigned address) {
    int index = searchInCache(address, false);
    if (index != -1) {
        Entries[index]->validBit = 0;
    }
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

public:
    Cache() = default;

    void
    Cache_initiate(unsigned MemCyc, unsigned BSize, unsigned L1Size, unsigned L2Size, unsigned L1Cyc, unsigned L2Cyc,
                   unsigned L1Assoc, unsigned L2Assoc, unsigned WrAlloc);

    Cache(Cache const &) = delete; // disable copy ctor
    void operator=(Cache const &) = delete; // disable = operator

    ~Cache() = default;

    void erase();

    static Cache &getInstance() // make singleton
    {
        static Cache instance; // Guaranteed to be destroyed.
        // Instantiated on first use.
        return instance;
    }

    void writeToCache(unsigned address);

    void readFromCache(unsigned address);

    void Snoop(Entry *oldBlock);

    void CacheGetStats(Stat *statsL1, Stat *statsL2, Stat *statsMem);
};


/**********************************************
 * Cache Methods
 **********************************************/


void
Cache::Cache_initiate(unsigned MemCyc, unsigned BSize, unsigned L1Size, unsigned L2Size, unsigned L1Cyc, unsigned L2Cyc,
                      unsigned L1Assoc, unsigned L2Assoc, unsigned WrAlloc) {
    this->MemCyc = MemCyc;
    this->BSize = pow(2,BSize);
    this->Allocate = static_cast<Alloc>(WrAlloc);


    //calculate the number of entries in each cache
    unsigned num_blocks1 = pow(2,L1Size) / this->BSize;
    unsigned num_blocks2 = pow(2,L2Size) / this->BSize;

    this->L1 = CacheTable(this->BSize, num_blocks1, L1Assoc, L1Cyc, Allocate);
    this->L2 = CacheTable(this->BSize, num_blocks2, L2Assoc, L2Cyc, Allocate);
};

void Cache::erase(){
    L1.erase();
    L2.erase();
}


void Cache::readFromCache(unsigned address) {
    Entry *oldBlock = new Entry();

    //try to read from L1
    L1.stats.numOfAcc++;
    if (!L1.read(address, true)) {
        L1.stats.numOfMiss++;
        L2.stats.numOfAcc++;

        //try to read from L2
        if (!L2.read(address, true)) {
            L2.stats.numOfMiss++;
            L2.bringFromLower(address, oldBlock);
            if (oldBlock->validBit)
                Snoop(oldBlock);
        }
        oldBlock->address=0;
        oldBlock->DirtyBit=0;
        oldBlock->validBit=0;

        L1.bringFromLower(address, oldBlock); // fetch
        if (oldBlock->DirtyBit && oldBlock->validBit) //TODO: added validBit check (might not matter)
        {
            L2.write(oldBlock->address, true); // TODO: changed address to oldBlock->address
        }

    }
    delete oldBlock;

}


void Cache::writeToCache(unsigned int address) {
    Entry *oldBlock = new Entry();
    if (Allocate == WRITE_ALLOCATE) {
        //search in L1
        L1.stats.numOfAcc++;
        if (L1.read(address, false)) {
            L1.write(address, true);
        }
        else
        {
            L1.stats.numOfMiss++;
            L2.stats.numOfAcc++;

            //search L2
            if (!L2.read(address, true))
            {
                L2.stats.numOfMiss++;
                L2.bringFromLower(address, oldBlock);
                if (oldBlock->validBit) //TODO: added valid bit check (might not matter)
                {
                    Snoop(oldBlock);
                }
            }

            oldBlock->address=0;
            oldBlock->DirtyBit=0;
            oldBlock->validBit=0;
            //free block X from L1 and write X to L2 if dirty
            L1.bringFromLower(address, oldBlock);
            if (oldBlock->DirtyBit && oldBlock->validBit)
            {
                L2.write(oldBlock->address, true);
            }

            L1.write(address, true);

        }
    } else {
        L1.stats.numOfAcc++;
        if (L1.read(address, false)) {
            L1.write(address, true);
        } else {
            L1.stats.numOfMiss++;
            L2.stats.numOfAcc++;
            if (L2.read(address, false))
                L2.write(address, true);
            else
                L2.stats.numOfMiss++;
        }
    }

    delete oldBlock;
}


void Cache::Snoop(Entry *oldBlock) {
    //block is evicted from L2 but still in L1
    if (L1.searchInCache(oldBlock->address,false) !=-1) {
        L1.Evict(oldBlock->address);
       // L2.write(oldBlock->address);
    }
}
void Cache::CacheGetStats(Stat * statsL1, Stat *statsL2, Stat *statsMem){
    //L1
    statsL1->numOfAcc = L1.stats.numOfAcc;
    statsL1->numOfMiss = L1.stats.numOfMiss;
    statsL1->accTime = L1.stats.accTime;

    //L2
    statsL2->numOfAcc = L2.stats.numOfAcc;
    statsL2->numOfMiss = L2.stats.numOfMiss;
    statsL2->accTime = L2.stats.accTime;

    //memory
    statsMem->numOfAcc = L2.stats.numOfMiss;
    statsMem->accTime = MemCyc;
}

/**********************************************
 * Main Functions
 **********************************************/


int Cache_init(unsigned MemCyc, unsigned BSize, unsigned L1Size, unsigned L2Size, unsigned L1Cyc, unsigned L2Cyc,
               unsigned L1Assoc, unsigned L2Assoc, unsigned WrAlloc) {
    Cache &cache = Cache::getInstance();
    if (&cache == 0) return -1;
    cache.Cache_initiate(MemCyc, BSize, L1Size, L2Size, L1Cyc, L2Cyc, L1Assoc, L2Assoc, WrAlloc);
    return 1;
}



void Cache_Access(char operation, unsigned address) {
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


void Cache_GetStats(double *L1MissRate, double *L2MissRate, double *avgAccTime) {
    Cache &cache = Cache::getInstance();
    Stat statsL1;
    Stat statsL2;
    Stat memory;
    cache.CacheGetStats(&statsL1, &statsL2, &memory);

    *L1MissRate = statsL1.numOfMiss / statsL1.numOfAcc;
    *L2MissRate = statsL2.numOfMiss / statsL2.numOfAcc;
    *avgAccTime = (statsL1.accTime * statsL1.numOfAcc +
                   statsL2.accTime * statsL2.numOfAcc +
                   memory.accTime * memory.numOfAcc) / (statsL1.numOfAcc);

    cache.erase();

}
