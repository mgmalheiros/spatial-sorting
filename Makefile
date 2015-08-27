OPT = -O2 -Wall -Wextra
INC = 
LIB = -lrt

CGAL    = -lCGAL
CRYPTO  = -lcrypto
GSL     = -lgsl -lgslcblas

all: test-1-full-sort test-2-time-complexity test-3-gather-knn test-4-gather-range test-5-gather-knn-3d

test-1-full-sort: common.hpp algorithms.cpp distributions.cpp test-1-full-sort.cpp
	g++ $(OPT) $(INC) -o test-1-full-sort algorithms.cpp distributions.cpp test-1-full-sort.cpp $(LIB) $(CRYPTO) $(GSL)

test-2-time-complexity: common.hpp algorithms.cpp distributions.cpp test-2-time-complexity.cpp
	g++ $(OPT) $(INC) -o test-2-time-complexity algorithms.cpp distributions.cpp test-2-time-complexity.cpp $(LIB) $(GSL)

test-3-gather-knn: common.hpp algorithms.cpp distributions.cpp test-3-gather-knn.cpp
	g++ $(OPT) $(INC) -o test-3-gather-knn algorithms.cpp distributions.cpp test-3-gather-knn.cpp $(LIB) $(CGAL) $(GSL)

test-4-gather-range: common.hpp algorithms.cpp distributions.cpp test-4-gather-range.cpp
	g++ $(OPT) $(INC) -o test-4-gather-range algorithms.cpp distributions.cpp test-4-gather-range.cpp $(LIB) $(CGAL) $(GSL)

test-5-gather-knn-3d: common.hpp algorithms.cpp distributions.cpp test-5-gather-knn-3d.cpp
	g++ $(OPT) $(INC) -o test-5-gather-knn-3d algorithms.cpp distributions.cpp test-5-gather-knn-3d.cpp $(LIB) $(CGAL) $(GSL)

clean:
	rm -f test-1-full-sort test-2-time-complexity test-3-gather-knn test-4-gather-range test-5-gather-knn-3d

