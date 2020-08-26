include ${FSLCONFDIR}/default.mk

PROJNAME = fabber_dsc

USRLDFLAGS = -Ldscprob
LIBS =  -lfsl-fabber_models_dsc -lfsl-fabberexec -lfsl-fabbercore \
        -lfsl-dscprob -lfsl-newimage -lfsl-miscmaths -lfsl-utils \
        -lfsl-NewNifti -lfsl-dscprob -lfsl-znz -ldl

XFILES = fabber_dsc
SOFILES = libfsl-fabber_models_dsc.so dscprob/libfsl-dscprob.so

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
dscprob/libfsl-dscprob.so:
	$(MAKE) -C dscprob CC="${CC}" CFLAGS="${CFLAGS}"

all: ${XFILES} ${SOFILES}

clean:
	${RM} -f *.o *.so dscprob/*.o dscprob/*.so depend.mk fabber_dsc

# models in a library
libfsl-fabber_models_dsc.so : ${OBJS}
	$(CXX) $(CXXFLAGS) -shared -o $@ $^

# fabber built from the FSL fabbercore library including the models specifieid in this project
fabber_dsc : fabber_client.o libfsl-fabber_models_dsc.so dscprob/libfsl-dscprob.so
	${CXX} ${CXXFLAGS} -o $@ $< ${LDFLAGS}

# DO NOT DELETE
