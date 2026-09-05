VPATH = src include 
CPPFLAGS = -I src/include
CC = pgc++ 
COMPILE.c := $(CC) $(CFLAGS) $(CPPFLAGS) $(TRAGET_ARCH) -c
OUTPUT_OPTION = -o $@
#OPT = -acc -L/hpc2ssd/softwares/cuda/cuda-12.2/lib64 -lcudart -lcurand -lcufft

OPT = -acc -L//home/soft/cuda-12.8/lib64 -lcudart -lcurand -lcufft

#OPT += -lgsl -lgslcblas 
#OPT +=  -lfftw3
OPT += -D_8BLOCKLOOP_
OPT += -D_NOISE_

#--EXECUTABLE--
main: main.cpp io.o interaction.o cellList.o 

io.o: io.cpp
poisson.o: poisson.cpp
poisson_cfft_c2c_omp.o: poisson_cfft_c2c_omp.cpp
cellList.o: cellList.cpp
interaction.o: interaction.cpp
clean:
	rm -f *.o
#--SUFFIX RULES--
%.o: %.c
	$(COMPILE.c) $(OUTPUT_OPTION) $< -g -Wall
%.o: %.cpp
	$(COMPILE.c) $(OUTPUT_OPTION) $< -g -O3 -Wall $(OPT)
%: %.cpp
	$(CC) $(CPPFLAGS) $^ -o $@ -g -O3 -Wall $(OPT)