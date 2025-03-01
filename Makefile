CXX:=g++ -std=c++20
LD:=g++ -std=c++20
CXXFLAGS:=-Wall -O3 -lz -lboost_iostreams
CPPFLAGS:=
LDFLAGS:=-lz -lboost_iostreams


.PHONY: all parsing parsing_yeast common mean_coverage modal_coverage list_mut mut_rate_clonal mut_rate_clonal_yeast mut_rate_clonal_yeast_hap list_coord sign_coord

all: parsing parsing_yeast common mean_coverage modal_coverage list_mut mut_rate_clonal mut_rate_clonal_yeast mut_rate_clonal_yeast_hap list_coord sign_coord

clean: cparsing cparsing_yeast ccommon cmean_coverage cmodal_coverage clist_mut cmut_rate_clonal cmut_rate_clonal_yeast cmut_rate_clonal_yeast_hap clist_coord csign_coord

#Parsing
parsing: parsing/exec/parsing

parsing/exec/parsing: parsing/obj/parsing.o
	${LD} $^ -o $@ ${LDFLAGS}

parsing/obj/parsing.o: parsing/src/parsing.cpp
	mkdir -p parsing/obj
	${CXX} -I pileup -D 'M_DIR="${mut_dir}"' ${CXXFLAGS} -c $< -o $@ ${CPPFLAGS}

cparsing:
	rm -f parsing/obj/parsing.o
	rm -f parsing/exec/parsing

#Parsing yeast
parsing_yeast: parsing/exec/parsing_yeast

parsing/exec/parsing_yeast: parsing/obj/parsing_yeast.o
	${LD} $^ -o $@ ${LDFLAGS}

parsing/obj/parsing_yeast.o: parsing/src/parsing_yeast.cpp
	mkdir -p parsing/obj
	${CXX} -I pileup -D 'M_DIR="${mut_dir}"' ${CXXFLAGS} -c $< -o $@ ${CPPFLAGS}

cparsing_yeast:
	rm -f parsing/obj/parsing_yeast.o
	rm -f parsing/exec/parsing_yeast

#Common lines in ancestor and endpoint
common: common/exec/common

common/exec/common: common/obj/common.o
	${LD} $^ -o $@

common/obj/common.o: common/src/common.cpp
	mkdir -p common/obj
	${CXX} -I pileup -D 'M_DIR="${mut_dir}"' -Wall -O3 -c $< -o $@

ccommon:
	rm -f common/obj/common.o
	rm -f common/exec/common

#Mean coverge and standard deviation for pseudo-pileup common files
mean_coverage: mean_coverage/exec/mean_coverage

mean_coverage/exec/mean_coverage: mean_coverage/obj/mean_coverage.o
	${LD} $^ -o $@

mean_coverage/obj/mean_coverage.o: mean_coverage/src/mean_coverage.cpp
	mkdir -p mean_coverage/obj
	${CXX} -I pileup -D 'M_DIR="${mut_dir}"' -Wall -O3 -c $< -o $@

cmean_coverage:
	rm -f mean_coverage/obj/mean_coverage.o
	rm -f mean_coverage/exec/mean_coverage

#Modal coverage and counts
modal_coverage: modal_coverage/exec/modal_coverage

modal_coverage/exec/modal_coverage: modal_coverage/obj/modal_coverage.o
	${LD} $^ -o $@

modal_coverage/obj/modal_coverage.o: modal_coverage/src/modal_coverage.cpp
	mkdir -p modal_coverage/obj
	${CXX} -I pileup -D 'M_DIR="${mut_dir}"' -Wall -O3 -c $< -o $@

cmodal_coverage:
	rm -f modal_coverage/obj/modal_coverage.o
	rm -f modal_coverage/exec/modal_coverage

#List mutations
list_mut: list_mutations/exec/list_mutations

list_mutations/exec/list_mutations: list_mutations/obj/list_mutations.o
	${LD} $^ -o $@

list_mutations/obj/list_mutations.o: list_mutations/src/list_mutations.cpp
	mkdir -p list_mutations/obj
	${CXX} -I pileup -D 'M_DIR="${mut_dir}"' -Wall -O3 -c $< -o $@

clist_mut:
	rm -f list_mutations/obj/list_mutations.o
	rm -f list_mutations/exec/list_mutations

#mutation rate clonal mutations CRC
mut_rate_clonal: CRC_clonal_mr/exec/clonal

CRC_clonal_mr/exec/clonal: CRC_clonal_mr/obj/clonal.o
	${LD} $^ -o $@

CRC_clonal_mr/obj/clonal.o: CRC_clonal_mr/src/clonal.cpp
	mkdir -p CRC_clonal_mr/obj
	${CXX} -I pileup -D 'M_DIR="${mut_dir}"' -Wall -O3 -c $< -o $@

cmut_rate_clonal:
	rm -f CRC_clonal_mr/obj/clonal.o
	rm -f CRC_clonal_mr/exec/clonal

#mutation rate clonal mutations yeast
mut_rate_clonal_yeast: YEAST_clonal_mr/exec/clonal

YEAST_clonal_mr/exec/clonal: YEAST_clonal_mr/obj/clonal.o
	${LD} $^ -o $@

YEAST_clonal_mr/obj/clonal.o: YEAST_clonal_mr/src/clonal.cpp
	mkdir -p YEAST_clonal_mr/obj
	${CXX} -I pileup -D 'M_DIR="${mut_dir}"' -Wall -O3 -c $< -o $@

cmut_rate_clonal_yeast:
	rm -f YEAST_clonal_mr/obj/clonal.o
	rm -f YEAST_clonal_mr/exec/clonal

#mutation rate clonal mutations yeast haploid
mut_rate_clonal_yeast_hap: YEAST_clonal_mr/exec/clonal_hap

YEAST_clonal_mr/exec/clonal_hap: YEAST_clonal_mr/obj/clonal_hap.o
	${LD} $^ -o $@

YEAST_clonal_mr/obj/clonal_hap.o: YEAST_clonal_mr/src/clonal_hap.cpp
	mkdir -p YEAST_clonal_mr/obj
	${CXX} -I pileup -D 'M_DIR="${mut_dir}"' -Wall -O3 -c $< -o $@

cmut_rate_clonal_yeast_hap:
	rm -f YEAST_clonal_mr/obj/clonal_hap.o
	rm -f YEAST_clonal_mr/exec/clonal_hap

#Coordinates of the mutations
list_coord: signatures/exec/list_coord_mut

signatures/exec/list_coord_mut: signatures/obj/list_coord_mut.o
	${LD} $^ -o $@

signatures/obj/list_coord_mut.o: signatures/src/list_coord_mut.cpp
	mkdir -p signatures/obj
	${CXX} -I pileup -D 'M_DIR="${mut_dir}"' -Wall -O3 -c $< -o $@

clist_coord:
	rm -f signatures/obj/list_coord_mut.o
	rm -f signatures/exec/list_coord_mut

#Signatures
sign_coord: signatures/exec/signatures_mut

signatures/exec/signatures_mut: signatures/obj/signatures_mut.o
	${LD} $^ -o $@

signatures/obj/signatures_mut.o: signatures/src/signatures_mut.cpp
	mkdir -p signatures/obj
	${CXX} -I pileup -D 'M_DIR="${mut_dir}"' -Wall -O3 -c $< -o $@

csign_coord:
	rm -f signatures/obj/signatures_mut.o
	rm -f signatures/exec/signatures_mut
