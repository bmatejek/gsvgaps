#
# Application name and list of source files
#

NAME=houghvote
CCSRCS=$(NAME).cpp helper.cpp R3CatmullRomSpline.cpp poseoptimization.cpp

#
# Dependency libraries
#

USER_CFLAGS=-DRN_USE_CSPARSE
PKG_LIBS=-lGSV -lR3Graphics -lR3Shapes -lR2Shapes -lRNMath -lRNBasics -lCSparse -ljpeg -lpng 

#
# R3 application makefile
#

include ../../makefiles/Makefile.apps
