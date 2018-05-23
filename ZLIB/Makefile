ZLIB_TARGET=libmyz.a
ZLIB_SOURCE=adler32.c compress.c crc32.c deflate.c gzio.c infblock.c infcodes.c inffast.c inflate.c inftrees.c infutil.c trees.c uncompr.c zutil.c

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

ZLIB_OBJECTS=$(addprefix $(BIN), $(addsuffix .o, $(basename $(ZLIB_SOURCE))))

all: CFLAGS += $(CFLAGS_RELEASE)
all: LFLAGS += $(LFLAGS_RELEASE)
all: make_dir
all: $(BIN)$(ZLIB_TARGET)

debug: CFLAGS += $(CFLAGS_DEBUG)
debug: LFLAGS += $(LFLAGS_DEBUG)
debug: make_dir
debug: $(BIN)$(ZLIB_TARGET)

zlib: CFLAGS += $(CFLAGS_RELEASE)
zlib: LFLAGS += $(LFLAGS_RELEASE)
zlib: make_dir
zlib: $(BIN)$(ZLIB_TARGET)

clean:
	rm -rf $(BIN)$(ZLIB_TARGET)
	rm -rf $(ZLIB_OBJECTS)

make_dir:
	$(MD) -p $(BIN)

$(BIN)$(ZLIB_TARGET): $(ZLIB_OBJECTS)
	ar rcs $(BIN)$(ZLIB_TARGET) $(ZLIB_OBJECTS)

$(BIN)%.o: $(SRC)%.c
	$(CC) -c -o $@ -I$(INCLUDE) $<

$(BIN)%.o: $(SRC)%.cpp
	$(CXX) -c -o $@ $(CFLAGS) -I$(INCLUDE) $<

