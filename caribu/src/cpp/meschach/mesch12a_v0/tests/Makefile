
LIBS = -L$(MESCHACH)/lib -lmeschach -lm

INCFLAG = -I$(MESCHACH)/include 
CFLAGS = -g -DANSI_C

C++ = g++
CC = gcc

.SUFFIXES :
.SUFFIXES : .C .c .o

.C.o : 
	$(C++) -c $(CFLAGS) $(INCFLAG) $< 


.c.o :
	$(CC) -c $(CFLAGS) $(INCFLAG) $<


SRCS = 
OBJS = ${SRCS:.c=.o}

CSRCS = iotort.c mfuntort.c torture.c  ztorture.c itertort.c \
        memtort.c sptort.c tutorial.c tutadv.c

COBJS = ${CSRCS:.c=.o}
EXEC = ${CSRCS:.c=}

$(EXEC)  : $$@.o $(OBJS)
	$(CC) $(CFLAGS) -o $@ $@.o $(OBJS) $(LIBS)

all : $(EXEC)

clean :
	-rm $(EXEC) $(COBJS) $(OBJS) 
