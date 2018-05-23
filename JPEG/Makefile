JPG_TARGET=libmyjpg.a
JPG_SOURCE=ckconfig.cpp jccolor.cpp jcmainct.cpp jcparam.cpp jctrans.cpp jdatasrc.cpp jdhuff.cpp jdmaster.cpp jdsample.cpp jfdctfst.cpp jidctint.cpp jquant1.cpp jcapimin.cpp jcdctmgr.cpp jcmarker.cpp jcphuff.cpp jdapimin.cpp jdcoefct.cpp jdinput.cpp jdmerge.cpp jdtrans.cpp jfdctint.cpp jidctred.cpp jquant2.cpp jcapistd.cpp jchuff.cpp jcmaster.cpp jcprepct.cpp jdapistd.cpp jdcolor.cpp jdmainct.cpp jdphuff.cpp jerror.cpp jidctflt.cpp jmemmgr.cpp jutils.cpp jccoefct.cpp jcinit.cpp jcomapi.cpp jcsample.cpp jdatadst.cpp jddctmgr.cpp jdmarker.cpp jdpostct.cpp jfdctflt.cpp jidctfst.cpp jmemnobs.cpp 

COMPILER = gcc
#COMPILER = clang

CFLAGS += -Wno-deprecated -Wno-write-strings

CFLAGS_DEBUG = -DDEBUG -g3
LFLAGS_DEBUG =

CFLAGS_RELEASE = -O3 -DRELEASE -funroll-loops -ffast-math -g
LFLAGS_RELEASE = -O3 -g

SRC = ./
BIN = ../Bin/Linux/
INCLUDE = ../

ifeq ($(COMPILER),gcc)
	CC=gcc
	CXX=g++
else
	CC=clang-3.8
	CXX=clang++-3.8
#	CC=clang-3.5
#	CXX=clang++-3.5
endif

MD=mkdir

JPG_OBJECTS=$(addprefix $(BIN), $(addsuffix .o, $(basename $(JPG_SOURCE))))

all: CFLAGS += $(CFLAGS_RELEASE)
all: LFLAGS += $(LFLAGS_RELEASE)
all: make_dir
all: $(BIN)
all: $(BIN)$(JPG_TARGET)

debug: CFLAGS += $(CFLAGS_DEBUG)
debug: LFLAGS += $(LFLAGS_DEBUG)
debug: $(BIN)
debug: $(BIN)$(JPG_TARGET)

jpeg: CFLAGS += $(CFLAGS_RELEASE)
jpeg: LFLAGS += $(LFLAGS_RELEASE)
jpeg: $(BIN)
jpeg: $(BIN)$(JPG_TARGET)

clean:
	rm -rf $(BIN)$(JPG_TARGET)
	rm -rf $(JPG_OBJECTS)

make_dir:
	$(MD) -p $(BIN)

$(BIN):
	$(MD) -p $(BIN)

$(BIN)$(JPG_TARGET): $(JPG_OBJECTS)
	ar rcs $(BIN)$(JPG_TARGET) $(JPG_OBJECTS)

$(BIN)%.o: $(SRC)%.c
	$(CC) -c -o $@ -I$(INCLUDE) $<

$(BIN)%.o: $(SRC)%.cpp
	$(CC) -c -o $@ $(CFLAGS) -I$(INCLUDE) $<

