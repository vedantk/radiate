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

.PHONY: test
test: build/radiate
	for file in test/render_*.sh; do \
		./$$file; \
	done

all-dbg:
	make clean
	make NDEBUG="$(DEBUG)"

all-pgo:
	make clean
	make NDEBUG="$(NDEBUG) -fprofile-generate -fprofile-correction"
	make test
	touch $(HEADERS)
	make NDEBUG="$(NDEBUG) -fprofile-use -fprofile-correction"
