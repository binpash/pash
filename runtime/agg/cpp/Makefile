.PHONY: clean
all: wc uniq tr

#CFLAGS = -Wall -std=c++17 -O3 -DNDEBUG -Icommon
CFLAGS = -Wall -std=c++17 -g -DDEBUG -Icommon
CFLAGS += common/main.cpp
DEPS = common/*

wc: $(DEPS) aggregators/wc/agg.h
	g++ ${CFLAGS} -Iaggregators/wc -o ./bin/wc
uniq: $(DEPS) aggregators/uniq/agg.h
	g++ ${CFLAGS} -Iaggregators/uniq -o ./bin/uniq
tr: $(DEPS) aggregators/tr/agg.h
	g++ ${CFLAGS} -Iaggregators/tr -o ./bin/tr

clean:
	rm -rf ./bin/*
