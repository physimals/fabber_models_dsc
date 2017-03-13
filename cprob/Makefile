PROJNAME = cprob

OBJS = bdtr.o btdtr.o chdtr.o drand.o expx2.o fdtr.o gamma.o gdtr.o \
igam.o igami.o incbet.o incbi.o mtherr.o nbdtr.o ndtr.o ndtri.o pdtr.o \
stdtr.o unity.o polevl.o const.o xmath.o

USRCFLAGS = -std=gnu99 -fPIC

all: libdscprob.a

libdscprob.a: $(OBJS) $(INCS)
	ar rv libdscprob.a $(OBJS)
	ranlib libdscprob.a

