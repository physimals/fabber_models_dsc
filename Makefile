include ${FSLCONFDIR}/default.mk

PROJNAME = fabber_dsc

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_PROB} -I${INC_BOOST} -I..
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_PROB} -L../fabber_core

LIBS = -lutils -lnewimage -lmiscmaths -lprob -lnewmat -lfslio -lniftiio -lznz -lz -ldl

XFILES = fabber_dsc

# Forward models
OBJS =  fwdmodel_dsc.o

# For debugging:
OPTFLAGS = -ggdb
#OPTFLAGS =

#
# Build
#

all:	${XFILES} libfabbermodels_dsc.a

# models in a library
libfabber_models_dsc.a : ${OBJS}
	${AR} -r $@ ${OBJS}

# fabber built from the FSL fabbercore library including the models specifieid in this project
fabber_dsc : fabber_client.o ${OBJS}
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ $< ${OBJS} -lfabbercore -lfabberexec ${LIBS}

# DO NOT DELETE
