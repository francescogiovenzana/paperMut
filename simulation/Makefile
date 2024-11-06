#seed for random number generator
seed:=/seed

#path for output data in simulation
sim_dir:=../

CXX:=g++ -std=c++20
LD:=g++ -std=c++20
CXXFLAGS:=-Wall -O3 -lgsl
CPPFLAGS:=
LDFLAGS:=-lgsl

.PHONY: all freqs

all: freqs

clean: cfreqs

#Random number generator
random/obj/random.o: random/random.cpp random/random.h
	mkdir -p random/obj
	${CXX} -I random -D 'SEED="${seed}"' ${CXXFLAGS} -c $< -o $@

#List frequencies with different input parameters
freqs: exec/list_freqs

exec/list_freqs: obj/list_freqs.o random/obj/random.o
	${LD} ${LDFLAGS} $^ -o $@

obj/list_freqs.o: src/list_freqs.cpp random/random.h
	mkdir -p obj
	${CXX} -I simulation -D 'S_DIR="${sim_dir}"' -I random -D 'SEED="${seed}"' ${CXXFLAGS} -c $< -o $@

cfreqs:
	rm -f obj/list_freqs.o
	rm -f exec/list_freqs
