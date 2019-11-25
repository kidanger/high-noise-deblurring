
CFLAGS+=-O3 -march=native -g
CXXFLAGS+=-O3 -march=native -g -std=c++11
CXXFLAGS+=-Wall -Wextra

LIBS+=-ljpeg -lpng -ltiff -lfftw3f

all: estimate-kernel

OBJS=iio.o downscale.o upsa.o

estimate-kernel: src/estimate_kernel.cpp src/image.cpp ${OBJS}
	${CXX} ${CXXFLAGS} $^ -o $@ ${LIBS}

%: %.cpp src/image.cpp ${OBJS}
	${CXX} ${CXXFLAGS} $^ -o $@ ${LIBS}

%.o: imscript/%.c
	${CC} ${CFLAGS} -c $^ -o $@ -DHIDE_ALL_MAINS -DOMIT_MAIN

clean:
	-rm ${OBJS} estimate-kernel

