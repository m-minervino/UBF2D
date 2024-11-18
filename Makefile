#Build UBF2D by Mauro Minervino

EXE = UBF2D
DBG = UBF2D_dbg
CC = icc
CFLAGS = -std=c++11
HDR = UBF2D.h
SRC = cellMetrics.cpp updateFirstTerm.cpp updateCompCorr.cpp updateFourthTerm.cpp updateParasite.cpp updateParasiteNew.cpp updateF_grad_rho.cpp compute_AF.cpp UBF2D.cpp
TECIOLIBS = /usr/local/apps/tecplot2022r1/360ex_2022r1/bin/libtecio.so
MTL4LIBS = -lpthread -I/u1/mima619/MTL4_GITHUB/trunk

default_target: $(EXE)

.PHONY : clean

clean:
	-rm $(EXE)

$(EXE): $(SRC) $(HDR) 
	$(CC) $(CFLAGS) -o $(EXE) $(SRC) $(TECIOLIBS) $(MTL4LIBS)

$(DBG): $(SRC) $(HDR)
	$(CC) $(CFLAGS) -check-pointers=rw -o $(EXE) $(SRC) $(TECIOLIBS) $(MTL4LIBS)

#icc -std=c++11 -check-pointers=rw -o UBF2D UBF2D.cpp /usr/local/apps/tecplot2022r1/360ex_2022r1/bin/libtecio.so -lpthread -I/u1/mima619/MTL4_GITHUB/trunk

#icc -check-pointers=rw -o UBF2D UBF2D.cpp /u1/mima619/TECIO_NEW/teciosrc/libtecio.a -lpthread -I/u1/mima619/MTL4_GITHUB/trunk
#icc -o UBF2D UBF2D.cpp /u1/mima619/TECIO_NEW/teciosrc/libtecio.a /usr/local/apps/tecplot2021r1/360ex_2021r1/bin/sys/libstdc++.so.6 -lpthread
#icc -o UBF2D UBF2D.cpp /usr/local/apps/tecplot2021r1/360ex_2021r1/bin/libtecio.so /usr/local/apps/tecplot2021r1/360ex_2021r1/bin/sys/libstdc++.so.6
#icc -check-pointers=rw -o UBF2D UBF2D.cpp /usr/local/apps/tecplot2021r1/360ex_2021r1/bin/libtecio.so /usr/local/apps/tecplot2021r1/360ex_2021r1/bin/sys/libstdc++.so.6

