PR_TARGET=PoissonRecon
SR_TARGET=SSDRecon
PI_TARGET=PointInterpolant
ST_TARGET=SurfaceTrimmer
EH_TARGET=EDTInHeat
IS_TARGET=ImageStitching
AV_TARGET=AdaptiveTreeVisualization
CP_TARGET=ChunkPLY
PR_SOURCE=PoissonRecon.cpp
SR_SOURCE=SSDRecon.cpp
PI_SOURCE=PointInterpolant.cpp
ST_SOURCE=SurfaceTrimmer.cpp
EH_SOURCE=EDTInHeat.cpp
IS_SOURCE=ImageStitching.cpp
AV_SOURCE=AdaptiveTreeVisualization.cpp
CP_SOURCE=ChunkPLY.cpp

COMPILER = gcc
#COMPILER = clang

ifeq ($(COMPILER),gcc)
	CFLAGS += -fopenmp -Wno-deprecated -std=c++11 -pthread -Wno-invalid-offsetof
	LFLAGS += -lgomp -lstdc++ -lpthread
else
# 	CFLAGS += -fopenmp=libiomp5 -Wno-deprecated -Wno-write-strings -std=c++11 -Wno-invalid-offsetof
# 	LFLAGS += -liomp5 -lstdc++
	CFLAGS += -Wno-deprecated -std=c++11 -pthread -Wno-invalid-offsetof
	LFLAGS += -lstdc++
endif
#LFLAGS += -lz -lpng -ljpeg

CFLAGS_DEBUG = -DDEBUG -g3
LFLAGS_DEBUG =

#CFLAGS_RELEASE = -O3 -DRELEASE -funroll-loops -ffast-math -g
#LFLAGS_RELEASE = -O3 -g
CFLAGS_RELEASE = -O3 -DRELEASE -funroll-loops -ffast-math -g
LFLAGS_RELEASE = -O3 -g

SRC = Src/
BIN = Bin/Linux/
#INCLUDE = /usr/include/
INCLUDE = .

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

PR_OBJECTS=$(addprefix $(BIN), $(addsuffix .o, $(basename $(PR_SOURCE))))
SR_OBJECTS=$(addprefix $(BIN), $(addsuffix .o, $(basename $(SR_SOURCE))))
PI_OBJECTS=$(addprefix $(BIN), $(addsuffix .o, $(basename $(PI_SOURCE))))
ST_OBJECTS=$(addprefix $(BIN), $(addsuffix .o, $(basename $(ST_SOURCE))))
EH_OBJECTS=$(addprefix $(BIN), $(addsuffix .o, $(basename $(EH_SOURCE))))
IS_OBJECTS=$(addprefix $(BIN), $(addsuffix .o, $(basename $(IS_SOURCE))))
AV_OBJECTS=$(addprefix $(BIN), $(addsuffix .o, $(basename $(AV_SOURCE))))
CP_OBJECTS=$(addprefix $(BIN), $(addsuffix .o, $(basename $(CP_SOURCE))))


all: CFLAGS += $(CFLAGS_RELEASE)
all: LFLAGS += $(LFLAGS_RELEASE)
all: make_dir
all: $(BIN)$(PR_TARGET)
all: $(BIN)$(SR_TARGET)
all: $(BIN)$(PI_TARGET)
all: $(BIN)$(ST_TARGET)
all: $(BIN)$(EH_TARGET)
all: $(BIN)$(IS_TARGET)
all: $(BIN)$(AV_TARGET)
all: $(BIN)$(CP_TARGET)

debug: CFLAGS += $(CFLAGS_DEBUG)
debug: LFLAGS += $(LFLAGS_DEBUG)
debug: make_dir
debug: $(BIN)$(PR_TARGET)
debug: $(BIN)$(SR_TARGET)
debug: $(BIN)$(PI_TARGET)
debug: $(BIN)$(ST_TARGET)
debug: $(BIN)$(EH_TARGET)
debug: $(BIN)$(IS_TARGET)
debug: $(BIN)$(AV_TARGET)
debug: $(BIN)$(CP_TARGET)

poissonrecon: CFLAGS += $(CFLAGS_RELEASE)
poissonrecon: LFLAGS += $(LFLAGS_RELEASE)
poissonrecon: make_dir
poissonrecon: $(BIN)$(PR_TARGET)

ssdrecon: CFLAGS += $(CFLAGS_RELEASE)
ssdrecon: LFLAGS += $(LFLAGS_RELEASE)
ssdrecon: make_dir
ssdrecon: $(BIN)$(SR_TARGET)

pointinterpolant: CFLAGS += $(CFLAGS_RELEASE)
pointinterpolant: LFLAGS += $(LFLAGS_RELEASE)
pointinterpolant: make_dir
pointinterpolant: $(BIN)$(PI_TARGET)

surfacetrimmer: CFLAGS += $(CFLAGS_RELEASE)
surfacetrimmer: LFLAGS += $(LFLAGS_RELEASE)
surfacetrimmer: make_dir
surfacetrimmer: $(BIN)$(ST_TARGET)

edtinheat: CFLAGS += $(CFLAGS_RELEASE)
edtinheat: LFLAGS += $(LFLAGS_RELEASE)
edtinheat: make_dir
edtinheat: $(BIN)$(EH_TARGET)

imagestitching: CFLAGS += $(CFLAGS_RELEASE)
imagestitching: LFLAGS += $(LFLAGS_RELEASE)
imagestitching: make_dir
imagestitching: $(BIN)$(IS_TARGET)

octreevisualization: CFLAGS += $(CFLAGS_RELEASE)
octreevisualization: LFLAGS += $(LFLAGS_RELEASE)
octreevisualization: make_dir
octreevisualization: $(BIN)$(AV_TARGET)

chunkply: CFLAGS += $(CFLAGS_RELEASE)
chunkply: LFLAGS += $(LFLAGS_RELEASE)
chunkply: make_dir
chunkply: $(BIN)$(CP_TARGET)

clean:
	rm -rf $(BIN)$(PR_TARGET)
	rm -rf $(BIN)$(SR_TARGET)
	rm -rf $(BIN)$(PI_TARGET)
	rm -rf $(BIN)$(ST_TARGET)
	rm -rf $(BIN)$(EH_TARGET)
	rm -rf $(BIN)$(IS_TARGET)
	rm -rf $(BIN)$(AV_TARGET)
	rm -rf $(BIN)$(CP_TARGET)
	rm -rf $(PR_OBJECTS)
	rm -rf $(SR_OBJECTS)
	rm -rf $(PI_OBJECTS)
	rm -rf $(ST_OBJECTS)
	rm -rf $(EH_OBJECTS)
	rm -rf $(IS_OBJECTS)
	rm -rf $(AV_OBJECTS)
	rm -rf $(CP_OBJECTS)
	cd PNG  && make clean


make_dir:
	$(MD) -p $(BIN)

$(BIN)$(PR_TARGET): $(PR_OBJECTS)
	cd PNG  && make
	$(CXX) -o $@ $(PR_OBJECTS) -L$(BIN) $(LFLAGS) -ljpeg -lmypng -lz

$(BIN)$(SR_TARGET): $(SR_OBJECTS)
	cd PNG  && make
	$(CXX) -o $@ $(SR_OBJECTS) -L$(BIN) $(LFLAGS) -ljpeg -lmypng -lz

$(BIN)$(PI_TARGET): $(PI_OBJECTS)
	cd PNG  && make
	$(CXX) -o $@ $(PI_OBJECTS) -L$(BIN) $(LFLAGS) -ljpeg -lmypng -lz

$(BIN)$(ST_TARGET): $(ST_OBJECTS)
	$(CXX) -o $@ $(ST_OBJECTS) $(LFLAGS)

$(BIN)$(EH_TARGET): $(EH_OBJECTS)
	$(CXX) -o $@ $(EH_OBJECTS) $(LFLAGS)

$(BIN)$(IS_TARGET): $(IS_OBJECTS)
	cd PNG  && make
	$(CXX) -o $@ $(IS_OBJECTS) -L$(BIN) $(LFLAGS) -ljpeg -lmypng -lz

$(BIN)$(AV_TARGET): $(AV_OBJECTS)
	cd PNG  && make
	$(CXX) -o $@ $(AV_OBJECTS) -L$(BIN) $(LFLAGS) -ljpeg -lmypng -lz

$(BIN)$(CP_TARGET): $(CP_OBJECTS)
	cd PNG  && make
	$(CXX) -o $@ $(CP_OBJECTS) -L$(BIN) $(LFLAGS) -ljpeg -lmypng -lz

$(BIN)%.o: $(SRC)%.c
	$(CC) -c -o $@ -I$(INCLUDE) $<

$(BIN)%.o: $(SRC)%.cpp
	$(CXX) -c -o $@ $(CFLAGS) -I$(INCLUDE) $<

