# Makefile for tspec








#--------------------------------------- Compile Time Options
# Run-time options
OPT += -DDEBUGGING
OPT += -DPROFILING

#--------------------------------------- Select Target Computer

SYSTYPE="CRC-Opteron-long"
#SYSTYPE="phillips"
#SYSTYPE="jared_home"

#--------------------------------------- System Specifics

ifeq ($(SYSTYPE),"jared_home")
CC = mpicc
OPTIMIZE = -Wall -g3
GSL_INCL = -I/usr/local/include/gsl
GSL_LIBS = -L/usr/local/lib
endif



ifeq ($(SYSTYPE),"phillips")
CC = gcc
OPTIMIZE = -Wall -g3
GSL_INCL = -I/opt/local/include
GSL_LIBS = -L/opt/local/lib -Wl
endif



ifeq ($(SYSTYPE),"CRC-Opteron-long")
CC       =  mpicc
#OPTIMIZE =  -O3 -Wall
OPTIMIZE = -Wall -O3
GSL_INCL = -I/opt/crc/g/gsl/1.16/intel/15.0/include
GSL_LIBS = -L/opt/crc/g/gsl/1.16/intel/15.0/lib -Wl,"-R /opt/crc/g/gsl/1.16/intel/15.0/lib"
#GSL_INCL =  -I/opt/crc/scilib/gsl/1.16/intel-14.0/include
#GSL_LIBS =  -L/opt/crc/scilib/gsl/1.16/intel-14.0/lib  -Wl,"-R /opt/crc/scilib/gsl/1.16/intel-14.0/lib"
FFTW_INCL=  -I/opt/crc/scilib/fftw/2.1.5_ompi/intel/1.6.5/include
FFTW_LIBS=  -L/opt/crc/scilib/fftw/2.1.5_ompi/intel/1.6.5/lib
#FFTW_INCL=  -I/opt/crc/scilib/fftw/2.1.5-intel/include
#FFTW_LIBS=  -L/opt/crc/scilib/fftw/2.1.5-intel/lib
MPICHLIB =	-L/opt/crc/openmpi/1.10.2/intel/15.0/lib
HDF5INCL =
HDF5LIB  =
endif



#-------------------------------------- Bookkeeping
PREFIX = ./src
OBJ_DIR = $(PREFIX)/obj

OPTIONS =  $(OPTIMIZE) $(OPT)

EXEC   = tspec

OBJS   = $(OBJ_DIR)/main.o $(OBJ_DIR)/allvars.o \
         $(OBJ_DIR)/de.o $(OBJ_DIR)/flag.o $(OBJ_DIR)/halos.o \
         $(OBJ_DIR)/load.o $(OBJ_DIR)/temperature.o \
         $(OBJ_DIR)/init.o $(OBJ_DIR)/debugging.o $(OBJ_DIR)/write_particle_data.o
   
INCL   = $(PREFIX)/allvars.h  $(PREFIX)/proto.h Makefile

INCLUDE = $(GSL_INCL)

CFLAGS = $(OPTIONS) $(GSL_INCL)

all: $(EXEC)  

LIBS = $(GSL_LIBS) -lm -lgsl -lgslcblas

#ALI Added 2/6/13  see http://www.apl.jhu.edu/Misc/Unix-info/make/make_10.html#SEC90
#                  for logic 
$(OBJ_DIR)/%.o : $(PREFIX)/%.c $(INCL)
	$(CC) $(OPTIONS) $(INCLUDE) -c $< -o $@

$(EXEC): $(OBJS) 
	$(CC) $(OBJS) $(LIBS) -o  $(EXEC)

clean:
	rm -f $(OBJS) *.gch 
