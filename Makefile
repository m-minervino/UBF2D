#Build UBF2D by Mauro Minervino

EXE = UBF2D
DBG = UBF2D_dbg
CC = icc
CFLAGS = -std=c++11
TECIOLIBS = /usr/local/apps/tecplot2022r1/360ex_2022r1/bin/libtecio.so
MTL4LIBS = -lpthread -I/u1/mima619/MTL4_GITHUB/trunk

default_target: $(EXE)

.PHONY : clean

clean:
	-rm $(EXE)

$(EXE): UBF2D.cpp
	$(CC) $(CFLAGS) -o $(EXE) UBF2D.cpp $(TECIOLIBS) $(MTL4LIBS)

$(DBG): UBF2D.cpp
	$(CC) $(CFLAGS) -check-pointers=rw -o $(EXE) UBF2D.cpp $(TECIOLIBS) $(MTL4LIBS)

#icc -std=c++11 -check-pointers=rw -o UBF2D UBF2D.cpp /usr/local/apps/tecplot2022r1/360ex_2022r1/bin/libtecio.so -lpthread -I/u1/mima619/MTL4_GITHUB/trunk

#icc -check-pointers=rw -o UBF2D UBF2D.cpp /u1/mima619/TECIO_NEW/teciosrc/libtecio.a -lpthread -I/u1/mima619/MTL4_GITHUB/trunk
#icc -o UBF2D UBF2D.cpp /u1/mima619/TECIO_NEW/teciosrc/libtecio.a /usr/local/apps/tecplot2021r1/360ex_2021r1/bin/sys/libstdc++.so.6 -lpthread
#icc -o UBF2D UBF2D.cpp /usr/local/apps/tecplot2021r1/360ex_2021r1/bin/libtecio.so /usr/local/apps/tecplot2021r1/360ex_2021r1/bin/sys/libstdc++.so.6
#icc -check-pointers=rw -o UBF2D UBF2D.cpp /usr/local/apps/tecplot2021r1/360ex_2021r1/bin/libtecio.so /usr/local/apps/tecplot2021r1/360ex_2021r1/bin/sys/libstdc++.so.6

