# File skb/makefile; TRANSIENT V2.51; 14 Aug. 2001
# (C) 2000-01 Miroslav Kolar


# This is an example of the UNIX makefile to maintain a TRANSIENT-based
# project - demonstration problem 1. It automatically updates TRANSIENT
# object modules (using 'make objects' in the main TRANSIENT directory) to
# correspond to the current state of the defs.h file before linking them
# with the demo1 PDM.
# You can use this makefile as a template for your own projects.

# First choose suitable options for your compiler/platform:

IGN = 0

# ----- gcc:
CC = gcc  -march=`uname -m` -O3 -funroll-all-loops -finline-functions -fexpensive-optimizations -ffast-math
# (Pick the CPU type, if -march with  uname  do not work for you):
#CC = gcc -m`uname -m`
#CC = gcc -mcpu=`uname -m`
#CC = gcc -mpentiumpro
#CC = gcc -march=K6
#CC = gcc -march=pentium
CFLAGS = -O3 
LDFLAGS = -s
IGN = 1  # To remind that some warnings can be ignored (for gcc only)

## ----- DEC (OSF) Unix cc:
#CC = cc -fast -arch host -tune host
#CFLAGS = -O4
#LDFLAGS =

## ----- HP Unix:
#CC = cc -Ae
#CFLAGS = -O
#LDFLAGS =

## ----- MS Visual C++ 6.0
#CC = cl
#CFLAGS =
#LDFLAGS =

LIBS = -lm

# ------- end of the section to be edited -------

MAKE = make


T_DIR = /home/martinet/transientV2.51/

OBJECTS = $(T_DIR)transien.o $(T_DIR)InOut.o $(T_DIR)dynarr.o\
 $(T_DIR)lufacv.o $(T_DIR)matinvv.o $(T_DIR)matrixv.o\
 $(T_DIR)nl3bandv.o $(T_DIR)DateTime.o $(T_DIR)signals.o

.c.o:
	$(CC) $(CFLAGS) -c $< -o $*.o

info:
	@echo "make skbbio, P. Martinet Oct 2001"
	@echo "If you type 'make ...' right now, the following settings will be used:"
	@echo ""; echo " CC		" $(CC)
	@echo " CFLAGS		" $(CFLAGS)
	@echo " LDFLAGS	" $(LDFLAGS)
	@echo " LIBS		" $(LIBS)
	@echo ""

skbbio:		skbbio.o objects
		$(CC) $(LDFLAGS) -o skbbio $(OBJECTS) skbbio.o $(LIBS)

skbbio.o:	skbbio.c skbbio.h $(T_DIR)transien.h $(T_DIR)dynarr.h $(T_DIR)incl_pri.h

skbio:		skbio.o objects
		$(CC) $(LDFLAGS) -o skbio $(OBJECTS) skbio.o $(LIBS)

skbio.o:	skbio.c skbio.h $(T_DIR)transien.h $(T_DIR)dynarr.h $(T_DIR)incl_pri.h

objects:
	( cd $(T_DIR); $(MAKE) 'CC=$(CC)' 'CFLAGS=$(CFLAGS)' 'IGN_remind=$(IGN)' objects )





