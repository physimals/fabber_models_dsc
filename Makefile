include ${FSLCONFDIR}/default.mk

PROJNAME = fabber_dsc

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_PROB} -I${INC_BOOST} 
USRLDFLAGS = -Ldscprob -L${LIB_NEWMAT} -L${LIB_PROB} -L../fabber_core

FSLVERSION= $(shell cat ${FSLDIR}/etc/fslversion | head -c 1)
ifeq ($(FSLVERSION), 5) 
  NIFTILIB = -lfslio -lniftiio 
  LIB_NEWMAT = -lnewmat
else 
  NIFTILIB = -lNewNifti
endif

LIBS = -lutils -lnewimage -lmiscmaths -ldscprob ${LIB_NEWMAT} ${NIFTILIB} -lznz -lz -ldl

XFILES = fabber_dsc

# Forward models
OBJS =  fwdmodel_dsc.o fwdmodel_dsc_cpi.o spline_interpolator.o

# For debugging:
#OPTFLAGS = -ggdb

# Pass Git revision details
GIT_SHA1:=$(shell git describe --match=NeVeRmAtCh --always --abbrev=40 --dirty)
GIT_DATE:=$(shell git log -1 --format=%ad --date=local)
CXXFLAGS += -DGIT_SHA1=\"${GIT_SHA1}\" -DGIT_DATE="\"${GIT_DATE}\""

#
# Build
#
dscprob/libdscprob.a:
	cd dscprob && $(MAKE)

all:	${XFILES} libfabber_models_dsc.a

clean:
	${RM} -f *.o *.a dscprob/*.o dscprob/*.a depend.mk fabber_dsc

# models in a library
libfabber_models_dsc.a : ${OBJS} 
	${AR} -r $@ ${OBJS}

# fabber built from the FSL fabbercore library including the models specifieid in this project
fabber_dsc : fabber_client.o ${OBJS} dscprob/libdscprob.a
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ $< ${OBJS} -lfabbercore -lfabberexec ${LIBS}

# DO NOT DELETE
