.PHONY: clean
all: wc uniq tr

CXX = g++10
CFLAGS = -Wall -std=c++17 -O3 -DNDEBUG -Icommon
CFLAGS += common/main.cpp
DEPS = common/*

wc: $(DEPS) aggregators/wc/agg.h
	${CXX} ${CFLAGS} -Iaggregators/wc -o ./bin/wc
uniq: $(DEPS) aggregators/uniq/agg.h
	${CXX} ${CFLAGS} -Iaggregators/uniq -o ./bin/uniq
tr: $(DEPS) aggregators/tr/agg.h
	${CXX} ${CFLAGS} -Iaggregators/tr -o ./bin/tr

clean:
	rm -rf ./bin/*
