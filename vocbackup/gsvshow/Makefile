#
# Application name and list of source files.
#

NAME=gsvshow
CCSRCS=$(NAME).cpp poseoptimization.cpp R3CatmullRomSpline.cpp


#
# Dependency libraries
#

USER_CFLAGS=-I/usr/local/include -I/usr/include/suitesparse -I/usr/include/eigen3 -DRN_USE_CSPARSE
PKG_LIBS=-lGSV -lR3Graphics -lR3Shapes -lR2Shapes -lRNMath -lRNBasics -ljpeg -lpng -lCSparse

#USER_CFLAGS=-DRN_USE_MINPACK -DRN_USE_SPLM -DRN_USE_CSPARSE
#PKG_LIBS=-lGSV -lR3Graphics -lR3Shapes -lR2Shapes -lRNMath -lRNBasics -ljpeg -lpng -lsplm -lCSparse -lminpack -lklt



#
# Ceres
#

#USER_CFLAGS=-Ic:/Funk/downloads/ceres/ceres-solver-1.3.0/include -I/usr/local/include/eigen3 -DRN_USE_CERES -DRN_USE_MINPACK -DRN_USE_SPLM -DRN_USE_CSPARSE
#PKG_LIBS=-lGSV -lR3Graphics -lR3Shapes -lR2Shapes -lRNMath -lRNBasics -ljpeg -lpng -lceres -lsplm -lCSparse -lminpack -lccolamd -lcamd -lcolamd -lamd -Wl,-Bdynamic /usr/local/lib/libsuitesparseconfig.a -lcholmod -lccolamd -lcamd -lcolamd -lamd -llapack -lblas -Wl,-Bstatic -lcxsparse -Wl,-Bdynamic -lgomp /usr/local/lib/libgflags.a /usr/local/lib/libglog.a /usr/local/lib/libprotobuf.dll.a -lklt



#
# R3 application makefile
#

include ../../makefiles/Makefile.apps



