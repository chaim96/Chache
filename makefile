# 046267 Computer Architecture - Winter 20/21 - HW #2

cacheSim: cacheSim.cpp
	g++ -std=c++11 -Wall -o cacheSim cacheSim.cpp cacheImplement.cpp

.PHONY: clean
clean:
	rm -f *.o
	rm -f cacheSim