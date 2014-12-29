CXX = g++
DEBUG_FLAGS = -g
NDEBUG_FLAGS = -DNDEBUG -Ofast
CXXFLAGS = -std=c++11 -fopenmp -Wall -Wextra -I../ $(NDEBUG_FLAGS)
LDFLAGS = -lstdc++ -lassimp -lboost_system -lfreeimage
GEN_PROFILE = -fprofile-generate="build/profile"
USE_PROFILE = -fprofile-use="build/profile"
OBJS = $(patsubst %.cc, build/%.o, $(wildcard *.cc))
HEADERS = $(wildcard *.hh)

all: build build/radiate

build:
	mkdir -p build

build/radiate: $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS)

build/%.o: %.cc $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $@ -c $<

clean:
	rm -rf build

test: build/radiate
	./test/render_*.sh

all-pgo:
	make clean
	make NDEBUG_FLAGS="$(NDEBUG_FLAGS) $(GEN_PROFILE)"
	make test
	make clean
	make NDEBUG_FLAGS="$(NDEBUG_FLAGS) $(USE_PROFILE)"
