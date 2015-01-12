
PROJECT_ROOT   = ..
LIBDIR         = ${PROJECT_ROOT}/lib

# Individual Program Modifications Should Begin Here

# Program to make (binary)
EXECUTABLE     = slclust

# Source files for Executable
SRCS            = slcluster.cpp graph.cpp graphnode.cpp cmd_line_opts.cpp

# Object modules for Executable
OBJS           = slcluster.o graph.o graphnode.o cmd_line_opts.o

# LIBS needed to compile the Executable (eg: /usr/local/lib/libspecial.a)
# use ${LIBDIR} defined above if appropriate
LIBS           = 

# How they are invoked on the compile line (eg: -lspecial)
LLIBS          =

# Local additions for the CFLAG options
LOCAL_CFLAGS   = -Wall 
#-static



# ************************************************************************
# THE REST OF THE MAKEFILE SHOULD NOT NEED TO BE MODIFIED
# (EXCEPT BY THE PROJECT SOFTWARE DEVELOPMENT MANAGER)

BIN           = ${PROJECT_ROOT}/bin
INCLUDE       = ${PROJECT_ROOT}/include

DEBUG         = DEBUG
CFLAGS        = -I${INCLUDE} ${LOCAL_CFLAGS}
CC            = g++ ${CFLAGS}

MAKEFILE      = Makefile


# Suffix rules

.cc.o:
	${CC} -c $<

.cpp.o:
	${CC} -c $<

# Target dependencies


all : ${EXECUTABLE}


install : ${EXECUTABLE}
	mv ${EXECUTABLE} ${BIN}/


clean :
	rm -f ${OBJS} core a.out *~ \#* ${EXECUTABLE} ${MAKEFILE}.bak


${OBJ} : ${MAKEFILE}

${EXECUTABLE} : ${OBJS}
	${CC} ${OBJS} ${LIBS} -o ${EXECUTABLE}	
	chmod 755 ${EXECUTABLE}

