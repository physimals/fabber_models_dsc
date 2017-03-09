SET(LIBPROB_SOURCES
  bdtr.c btdtr.c chdtr.c drand.c expx2.c fdtr.c gamma.c gdtr.c igam.c
  igami.c incbet.c incbi.c mtherr.c nbdtr.c ndtr.c ndtri.c pdtr.c
  stdtr.c unity.c polevl.c const.c xmath.c)

ADD_LIBRARY(dscprob ${LIBPROB_SOURCES})
set_property(TARGET dscprob PROPERTY POSITION_INDEPENDENT_CODE TRUE)
