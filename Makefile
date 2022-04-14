OBJ_DIR := objects
SRC_DIR := src
INC_DIR := include

SRC := $(wildcard $(SRC_DIR)/*.cpp)
EXE := $(patsubst $(SRC_DIR)/%.cpp, %, $(SRC))
RAT := -L$(RATROOT)/lib -lRATEvent -I$(RATROOT)/include -I$(INC_DIR) -lboost_system -lboost_filesystem
ROOT := $(shell root-config --cflags --libs)
FLAGS := -std=c++11
GCC := g++

_OBJS := $(patsubst $(SRC_DIR)/%.cc, %.o, $(wildcard $(SRC_DIR)/*.cc))
OBJS := $(patsubst %,$(OBJ_DIR)/%,$(_OBJS))

all: $(SRC) $(EXE) $(OBJS)

%: $(SRC_DIR)/%.cpp $(OBJS)
	$(GCC) $(FLAGS) $(ROOT) $(RAT) $(OBJS) $(SRC_DIR)/$@.cpp -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cc
	$(GCC) $(FLAGS) $(ROOT) $(RAT) -shared -c -fPIC $(SRC_DIR)/$*.cc -o $(OBJ_DIR)/$*.o

clean:
	rm -f $(EXE)
	rm -f $(OBJ_DIR)/*.o
