PR_TARGET=PoissonRecon
PRC_TARGET=PoissonReconClient
PRS_TARGET=PoissonReconServer
SR_TARGET=SSDRecon
PI_TARGET=PointInterpolant
ST_TARGET=SurfaceTrimmer
EH_TARGET=EDTInHeat
IS_TARGET=ImageStitching
AV_TARGET=AdaptiveTreeVisualization
CP_TARGET=ChunkPLY
RE_TARGET=ReconExample
PTD_TARGET=PointsToDisks
SN_TARGET=ScaleNormals

PR_SOURCE=PoissonRecon.cpp
PRC_SOURCE=PoissonReconClient.cpp
PRS_SOURCE=PoissonReconServer.cpp
SR_SOURCE=SSDRecon.cpp
PI_SOURCE=PointInterpolant.cpp
ST_SOURCE=SurfaceTrimmer.cpp
EH_SOURCE=EDTInHeat.cpp
IS_SOURCE=ImageStitching.cpp
AV_SOURCE=AdaptiveTreeVisualization.cpp
CP_SOURCE=ChunkPLY.cpp
RE_SOURCE=Reconstruction.example.cpp
PTD_SOURCE=PointsToDisks.cpp
SN_SOURCE=ScaleNormals.cpp

COMPILER ?= gcc
#COMPILER ?= clang

ifeq ($(COMPILER),gcc)
	CFLAGS += -fopenmp -Wno-deprecated -std=c++17 -pthread -Wno-invalid-offsetof
	LFLAGS += -lgomp -lstdc++ -lpthread
else ifeq ($(COMPILER),gcc-11)
	CFLAGS += -fopenmp -Wno-deprecated -std=c++17 -pthread -Wno-invalid-offsetof -Werror=strict-aliasing -Wno-nonnull
	LFLAGS += -lgomp -lstdc++ -lpthread
else
# 	CFLAGS += -fopenmp=libiomp5 -Wno-deprecated -Wno-write-strings -std=c++17 -Wno-invalid-offsetof
# 	LFLAGS += -liomp5 -lstdc++
	CFLAGS += -Wno-deprecated -std=c++17 -pthread -Wno-invalid-offsetof -Wno-dangling-else
	CFLAGS += -Wno-nan-infinity-disabled
	LFLAGS += -lstdc++
endif
LFLAGS_IMG += -lz -lpng -ljpeg
#LFLAGS += -ljpeg -lmypng -lz

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
else ifeq ($(COMPILER),gcc-11)
	CC=gcc-11
	CXX=g++-11
else
	CC=clang
	CXX=clang++
endif

MD=mkdir

PR_OBJECTS=$(addprefix $(BIN), $(addsuffix .o, $(basename $(PR_SOURCE))))
PRC_OBJECTS=$(addprefix $(BIN), $(addsuffix .o, $(basename $(PRC_SOURCE))))
PRS_OBJECTS=$(addprefix $(BIN), $(addsuffix .o, $(basename $(PRS_SOURCE))))
SR_OBJECTS=$(addprefix $(BIN), $(addsuffix .o, $(basename $(SR_SOURCE))))
PI_OBJECTS=$(addprefix $(BIN), $(addsuffix .o, $(basename $(PI_SOURCE))))
ST_OBJECTS=$(addprefix $(BIN), $(addsuffix .o, $(basename $(ST_SOURCE))))
EH_OBJECTS=$(addprefix $(BIN), $(addsuffix .o, $(basename $(EH_SOURCE))))
IS_OBJECTS=$(addprefix $(BIN), $(addsuffix .o, $(basename $(IS_SOURCE))))
AV_OBJECTS=$(addprefix $(BIN), $(addsuffix .o, $(basename $(AV_SOURCE))))
CP_OBJECTS=$(addprefix $(BIN), $(addsuffix .o, $(basename $(CP_SOURCE))))
RE_OBJECTS=$(addprefix $(BIN), $(addsuffix .o, $(basename $(RE_SOURCE))))
PTD_OBJECTS=$(addprefix $(BIN), $(addsuffix .o, $(basename $(PTD_SOURCE))))
SN_OBJECTS=$(addprefix $(BIN), $(addsuffix .o, $(basename $(SN_SOURCE))))


all: CFLAGS += $(CFLAGS_RELEASE)
all: LFLAGS += $(LFLAGS_RELEASE)
all: make_dir
all: $(BIN)$(PR_TARGET)
all: $(BIN)$(PRC_TARGET)
all: $(BIN)$(PRS_TARGET)
all: $(BIN)$(SR_TARGET)
all: $(BIN)$(PI_TARGET)
all: $(BIN)$(ST_TARGET)
all: $(BIN)$(EH_TARGET)
all: $(BIN)$(IS_TARGET)
all: $(BIN)$(AV_TARGET)
all: $(BIN)$(CP_TARGET)
all: $(BIN)$(RE_TARGET)
all: $(BIN)$(PTD_TARGET)
all: $(BIN)$(SN_TARGET)

debug: CFLAGS += $(CFLAGS_DEBUG)
debug: LFLAGS += $(LFLAGS_DEBUG)
debug: make_dir
debug: $(BIN)$(PR_TARGET)
debug: $(BIN)$(PRC_TARGET)
debug: $(BIN)$(PRS_TARGET)
debug: $(BIN)$(SR_TARGET)
debug: $(BIN)$(PI_TARGET)
debug: $(BIN)$(ST_TARGET)
debug: $(BIN)$(EH_TARGET)
debug: $(BIN)$(IS_TARGET)
debug: $(BIN)$(AV_TARGET)
debug: $(BIN)$(CP_TARGET)
debug: $(BIN)$(RE_TARGET)
debug: $(BIN)$(PTD_TARGET)
debug: $(BIN)$(SN_TARGET)

poissonrecon: CFLAGS += $(CFLAGS_RELEASE)
poissonrecon: LFLAGS += $(LFLAGS_RELEASE)
poissonrecon: make_dir
poissonrecon: $(BIN)$(PR_TARGET)

poissonreconclient: CFLAGS += $(CFLAGS_RELEASE)
poissonreconclient: LFLAGS += $(LFLAGS_RELEASE)
poissonreconclient: make_dir
poissonreconclient: $(BIN)$(PRC_TARGET)

poissonreconserver: CFLAGS += $(CFLAGS_RELEASE)
poissonreconserver: LFLAGS += $(LFLAGS_RELEASE)
poissonreconserver: make_dir
poissonreconserver: $(BIN)$(PRS_TARGET)

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

reconexample: CFLAGS += $(CFLAGS_RELEASE)
reconexample: LFLAGS += $(LFLAGS_RELEASE)
reconexample: make_dir
reconexample: $(BIN)$(RE_TARGET)

pointstodisks: CFLAGS += $(CFLAGS_RELEASE)
pointstodisks: LFLAGS += $(LFLAGS_RELEASE)
pointstodisks: make_dir
pointstodisks: $(BIN)$(PTD_TARGET)

scalenormals: CFLAGS += $(CFLAGS_RELEASE)
scalenormals: LFLAGS += $(LFLAGS_RELEASE)
scalenormals: make_dir
scalenormals: $(BIN)$(SN_TARGET)

clean:
	rm -rf $(BIN)$(PR_TARGET)
	rm -rf $(BIN)$(PRC_TARGET)
	rm -rf $(BIN)$(PRS_TARGET)
	rm -rf $(BIN)$(SR_TARGET)
	rm -rf $(BIN)$(PI_TARGET)
	rm -rf $(BIN)$(ST_TARGET)
	rm -rf $(BIN)$(EH_TARGET)
	rm -rf $(BIN)$(IS_TARGET)
	rm -rf $(BIN)$(AV_TARGET)
	rm -rf $(BIN)$(CP_TARGET)
	rm -rf $(BIN)$(RE_TARGET)
	rm -rf $(PR_OBJECTS)
	rm -rf $(PRC_OBJECTS)
	rm -rf $(PRS_OBJECTS)
	rm -rf $(SR_OBJECTS)
	rm -rf $(PI_OBJECTS)
	rm -rf $(ST_OBJECTS)
	rm -rf $(EH_OBJECTS)
	rm -rf $(IS_OBJECTS)
	rm -rf $(AV_OBJECTS)
	rm -rf $(CP_OBJECTS)
	rm -rf $(RE_OBJECTS)
	rm -rf $(PTD_OBJECTS)
	rm -rf $(SN_OBJECTS)
	cd PNG  && make clean

make_dir:
	$(MD) -p $(BIN)

$(BIN)$(PR_TARGET): $(PR_OBJECTS)
	cd PNG  && make COMPILER=$(COMPILER)
	$(CXX) -pthread -o $@ $(PR_OBJECTS) -L$(BIN) $(LFLAGS) $(LFLAGS_IMG)

$(BIN)$(PRC_TARGET): $(PRC_OBJECTS)
	$(CXX) -pthread -o $@ $(PRC_OBJECTS) -L$(BIN) -lboost_system $(LFLAGS)

$(BIN)$(PRS_TARGET): $(PRS_OBJECTS)
	$(CXX) -pthread -o $@ $(PRS_OBJECTS) -L$(BIN) -lboost_system $(LFLAGS)

$(BIN)$(SR_TARGET): $(SR_OBJECTS)
	cd PNG  && make COMPILER=$(COMPILER)
	$(CXX) -pthread -o $@ $(SR_OBJECTS) -L$(BIN) $(LFLAGS) $(LFLAGS_IMG)

$(BIN)$(PI_TARGET): $(PI_OBJECTS)
	cd PNG  && make COMPILER=$(COMPILER)
	$(CXX) -pthread -o $@ $(PI_OBJECTS) -L$(BIN) $(LFLAGS) $(LFLAGS_IMG)

$(BIN)$(ST_TARGET): $(ST_OBJECTS)
	$(CXX) -pthread -o $@ $(ST_OBJECTS) $(LFLAGS)

$(BIN)$(EH_TARGET): $(EH_OBJECTS)
	$(CXX) -pthread -o $@ $(EH_OBJECTS) $(LFLAGS)

$(BIN)$(IS_TARGET): $(IS_OBJECTS)
	cd PNG  && make COMPILER=$(COMPILER)
	$(CXX) -pthread -o $@ $(IS_OBJECTS) -L$(BIN) $(LFLAGS) $(LFLAGS_IMG)

$(BIN)$(AV_TARGET): $(AV_OBJECTS)
	cd PNG  && make COMPILER=$(COMPILER)
	$(CXX) -pthread -o $@ $(AV_OBJECTS) -L$(BIN) $(LFLAGS) $(LFLAGS_IMG)

$(BIN)$(CP_TARGET): $(CP_OBJECTS)
	cd PNG  && make COMPILER=$(COMPILER)
	$(CXX) -pthread -o $@ $(CP_OBJECTS) -L$(BIN) $(LFLAGS) $(LFLAGS_IMG)

$(BIN)$(RE_TARGET): $(RE_OBJECTS)
	$(CXX) -pthread -o $@ $(RE_OBJECTS) -L$(BIN) $(LFLAGS) $(LFLAGS_IMG)

$(BIN)$(PTD_TARGET): $(PTD_OBJECTS)
	$(CXX) -pthread -o $@ $(PTD_OBJECTS) -L$(BIN) $(LFLAGS)

$(BIN)$(SN_TARGET): $(SN_OBJECTS)
	$(CXX) -pthread -o $@ $(SN_OBJECTS) -L$(BIN) $(LFLAGS)

$(BIN)%.o: $(SRC)%.c
	$(CC) -c -o $@ -I$(INCLUDE) $<

$(BIN)%.o: $(SRC)%.cpp
	$(CXX) -c -o $@ $(CFLAGS) -I$(INCLUDE) $<

