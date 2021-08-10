# To add a new aggregator, follow the steps below
# Let's assume the aggregator is for a command named cmd
.PHONY: clean
# Add cmd at the end of the line below
all: wc uniq tr grep

CXX = g++
CFLAGS = -Wall -std=c++17 -O3 -DNDEBUG -Icommon
CFLAGS += common/main.cpp
DEPS = common/*

ifeq ($(shell uname -s), Linux)
	AGG_PATH = aggregators/linux
endif
ifeq ($(shell uname -s), FreeBSD)
	AGG_PATH = aggregators/bsd
endif

# Copy one of these pairs of two lines and paste it at the end of the list
# And change each occurence of the command name (like wc) for cmd
wc: $(DEPS) $(AGG_PATH)/wc/agg.h
	${CXX} ${CFLAGS} -I$(AGG_PATH)/wc -o ./bin/wc
uniq: $(DEPS) $(AGG_PATH)/uniq/agg.h
	${CXX} ${CFLAGS} -I$(AGG_PATH)/uniq -o ./bin/uniq
tr: $(DEPS) $(AGG_PATH)/tr/agg.h
	${CXX} ${CFLAGS} -I$(AGG_PATH)/tr -o ./bin/tr
grep: $(DEPS) $(AGG_PATH)/grep/agg.h
	${CXX} ${CFLAGS} -I$(AGG_PATH)/grep -o ./bin/grep

# That's it for this file

clean:
	rm -rf ./bin/*
