
__PROJECT_ROOT   = ..
LIBDIR         = ${__PROJECT_ROOT}/lib

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

# Local additions for the CFLAG options
LOCAL_CXXFLAGS   = -Wall 
#-static



# ************************************************************************
# THE REST OF THE MAKEFILE SHOULD NOT NEED TO BE MODIFIED
# (EXCEPT BY THE PROJECT SOFTWARE DEVELOPMENT MANAGER)

BIN           = ${__PROJECT_ROOT}/bin
INCLUDE       = ${__PROJECT_ROOT}/include

DEBUG         = DEBUG
CXXFLAGS        += -I${INCLUDE} ${LOCAL_CXXFLAGS}

MAKEFILE      = Makefile


# Target dependencies


all : ${EXECUTABLE}


install : ${EXECUTABLE}
	mv ${EXECUTABLE} ${BIN}/


clean :
	rm -f ${OBJS} core a.out *~ \#* ${EXECUTABLE} ${MAKEFILE}.bak \
		${BIN}/${EXECUTABLE}


${OBJ} : ${MAKEFILE}

${EXECUTABLE} : ${OBJS}
	${CXX} ${LDFLAGS} ${OBJS} ${LIBS} -o ${EXECUTABLE}
	chmod 755 ${EXECUTABLE}

