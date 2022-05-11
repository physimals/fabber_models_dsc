include ${FSLCONFDIR}/default.mk

PROJNAME = fabber_dsc
XFILES   = fabber_dsc
SOFILES  = libfsl-fabber_models_dsc.so
AFILES   = libfabber_models_dsc.a

# The FSL build system changed
# substantially in FSL 6.0.6
# FSL >= 6.0.6
# In >=6.0.6 we dynamically link
# against the fsl-cprob library
ifeq (${FSL_GE_606}, true)
  LIBS         = -lfsl-fabberexec -lfsl-fabbercore \
                 -lfsl-newimage -lfsl-miscmaths -lfsl-utils \
                 -lfsl-NewNifti -lfsl-cprob -lfsl-znz -ldl
  USRCPPFLAGS += -DFSL_GE_606

# FSL <= 6.0.5
# In <=6.0.5 we statically link against
# our own copy of the cprob library (dscprob)
else
  ifeq ($(shell uname -s), Linux)
	MATLIB := -lopenblas
  endif
  USRINCFLAGS = -I${INC_NEWMAT} -I${INC_BOOST} -I${INC_CPROB} \
                -I${FSLDIR}/extras/include/armawrap
  USRLDFLAGS  = -L${LIB_NEWMAT} -I${LIB_CPROB} -L../fabber_core \
                -lfabbercore -lfabberexec -lnewimage \
                -lmiscmaths -lutils ${MATLIB} -lNewNifti \
                -lcprob -lznz -lm -lz -ldl
endif

# Forward models
OBJS =  fwdmodel_dsc.o fwdmodel_dsc_cpi.o spline_interpolator.o

# For debugging:
#OPTFLAGS = -ggdb

# Pass Git revision details
GIT_SHA1 := $(shell git describe --match=NeVeRmAtCh --always --abbrev=40 --dirty)
GIT_DATE := $(shell git log -1 --format=%ad --date=local)
CXXFLAGS += -DGIT_SHA1=\"${GIT_SHA1}\" -DGIT_DATE="\"${GIT_DATE}\""

#
# Build
#

clean:
	${RM} -f *.o *.so *.a dscprob/*.o dscprob/*.a depend.mk fabber_dsc

# FSL >=606 uses dynamic linking,
# and the standard fsl-cprob library.
ifeq (${FSL_GE_606}, true)

all: ${XFILES} ${SOFILES}

# models in a library
libfsl-fabber_models_dsc.so : ${OBJS}
	$(CXX) $(CXXFLAGS) -shared -o $@ $^ ${LDFLAGS}

# fabber built from the FSL fabbercore library including the models specifieid in this project
fabber_dsc : fabber_client.o | libfsl-fabber_models_dsc.so
	${CXX} ${CXXFLAGS} -o $@ $^ -lfsl-fabber_models_dsc ${LDFLAGS}

# FSL <=605 uses static linking, and
# a copy of the cprob library
else

all: ${XFILES} ${AFILES}

dscprob/libdscprob.a:
	$(MAKE) -C dscprob CC="${CC}" CFLAGS="${CFLAGS} -std=c99"

libfabber_models_dsc.a : ${OBJS}
	${AR} -r $@ $^

fabber_dsc : fabber_client.o ${OBJS} dscprob/libdscprob.a
	${CXX} ${CXXFLAGS} -o $@ $^ ${LDFLAGS}
endif
