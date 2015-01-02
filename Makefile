CXX = g++
DEBUG = -g
NDEBUG = -DNDEBUG -Ofast
CXXFLAGS = -std=c++11 -fopenmp -Wall -Wextra -I../ $(NDEBUG)
LDFLAGS = -lstdc++ -lassimp -lboost_system -lfreeimage
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
	make NDEBUG="$(NDEBUG) -fprofile-generate -fprofile-correction"
	make test
	touch *.hh
	make NDEBUG="$(NDEBUG) -fprofile-use -fprofile-correction"
