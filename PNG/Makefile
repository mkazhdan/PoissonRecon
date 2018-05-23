PNG_TARGET=libmypng.a
PNG_SOURCE=png.c pngerror.c pnggccrd.c pngget.c pngmem.c pngpread.c pngread.c pngrio.c pngrtran.c pngrutil.c pngset.c pngtest.c pngtrans.c pngvcrd.c pngwio.c pngwrite.c pngwtran.c pngwutil.c

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

PNG_OBJECTS=$(addprefix $(BIN), $(addsuffix .o, $(basename $(PNG_SOURCE))))


all: CFLAGS += $(CFLAGS_RELEASE)
all: LFLAGS += $(LFLAGS_RELEASE)
all: make_dir
all: $(BIN)$(PNG_TARGET)

debug: CFLAGS += $(CFLAGS_DEBUG)
debug: LFLAGS += $(LFLAGS_DEBUG)
debug: make_dir
debug: $(BIN)$(PNG_TARGET)

png: CFLAGS += $(CFLAGS_RELEASE)
png: LFLAGS += $(LFLAGS_RELEASE)
png: make_dir
png: $(BIN)$(PNG_TARGET)

clean:
	rm -rf $(BIN)$(PNG_TARGET)
	rm -rf $(PNG_OBJECTS)

make_dir:
	$(MD) -p $(BIN)

$(BIN)$(PNG_TARGET): $(PNG_OBJECTS)
	ar rcs $(BIN)$(PNG_TARGET) $(PNG_OBJECTS)

$(BIN)%.o: $(SRC)%.c
	$(CC) -c -o $@ -I$(INCLUDE) $<

$(BIN)%.o: $(SRC)%.cpp
	$(CXX) -c -o $@ $(CFLAGS) -I$(INCLUDE) $<

