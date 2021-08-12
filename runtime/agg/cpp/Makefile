.PHONY: all clean
all:

CXX = g++
CFLAGS = -Wall -std=c++17 -O3 -DNDEBUG -Icommon
CFLAGS += common/main.cpp
DEPS = common/*
TARGETS =

ifeq ($(shell uname -s), Linux)
	SYSTEM_NAME = linux
endif
ifeq ($(shell uname -s), FreeBSD)
	SYSTEM_NAME = bsd
endif

define add_aggregator

$1: $(DEPS) aggregators/$1/agg-${SYSTEM_NAME}.h
	mkdir -p intermediate/$1
	cp aggregators/$1/agg-${SYSTEM_NAME}.h intermediate/$1/agg.h
	${CXX} ${CFLAGS} -I intermediate/$1 -I aggregators/$1/ -o ./bin/$1
TARGETS += $1
endef

$(eval $(foreach agg,$(shell ls aggregators),$(call add_aggregator,$(agg))))

all: $(TARGETS)

clean:
	rm -rf ./bin/*
	rm -rf ./intermediate/*
