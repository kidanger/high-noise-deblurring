# Makefile created with Create Makefile.sh 13/12/2011

COPT      = -O99
CFLAGS    = $(COPT) -fopenmp
CSTRICT   = -Wall -Wextra -ansi
CXXOPT    = -O99 -I.
CXXFLAGS  = $(CXXOPT) -fopenmp
CXXSTRICT = -Wall -Wextra -ansi
LDFLAGS   = -lpng -ljpeg -ltiff -lm -lfftw3f -lfftw3f_threads -lgomp
EXEC      = ponomarenko

default: $(EXEC)
all: $(EXEC)

# ------- C files -------
# ------- C++ files -------
./framework/CFramework.o: ./framework/CFramework.cpp ./framework/CFramework.h
	$(CXX) $(CXXFLAGS) -c ./framework/CFramework.cpp -o ./framework/CFramework.o

./framework/libparser.o: ./framework/libparser.cpp ./framework/libparser.h
	$(CXX) $(CXXFLAGS) -c ./framework/libparser.cpp -o ./framework/libparser.o

./framework/CHistogram.o: ./framework/CHistogram.cpp
	$(CXX) $(CXXFLAGS) -c ./framework/CHistogram.cpp -o ./framework/CHistogram.o

./framework/CImage.o: ./framework/CImage.cpp ./framework/CImage.h
	$(CXX) $(CXXFLAGS) -c ./framework/CImage.cpp -o ./framework/CImage.o

./framework/CInspector.o: ./framework/CInspector.cpp ./framework/CInspector.h
	$(CXX) $(CXXFLAGS) -c ./framework/CInspector.cpp -o ./framework/CInspector.o

./framework/framework.o: ./framework/framework.cpp
	$(CXX) $(CXXFLAGS) -c ./framework/framework.cpp -o ./framework/framework.o

./framework/operations.o: ./framework/operations.cpp
	$(CXX) $(CXXFLAGS) -c ./framework/operations.cpp -o ./framework/operations.o

./curve_filter.o: ./curve_filter.cpp ./curve_filter.h
	$(CXX) $(CXXFLAGS) -c ./curve_filter.cpp -o ./curve_filter.o

./algo.o: ./algo.cpp ./algo.h
	$(CXX) $(CXXFLAGS) -c ./algo.cpp -o ./algo.o

# ------- Main -------
$(EXEC):  ./framework/CFramework.o ./framework/libparser.o ./framework/CHistogram.o ./framework/CImage.o ./framework/CInspector.o ./framework/framework.o ./framework/operations.o ./curve_filter.o ./algo.o iio.o
	$(CXX)  $^ $(LDFLAGS) -o $(EXEC)

lint:
	$(MAKE) CFLAGS="$(CFLAGS) $(CSTRICT)" CXXFLAGS="$(CXXFLAGS) $(CXXSTRICT)"

clean:
	rm -f *.o

distclean: clean
	rm -f $(EXEC)

