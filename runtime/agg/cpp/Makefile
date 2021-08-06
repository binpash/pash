# To add a new aggregator, follow the steps below
# Let's assume the aggregator is for a command named cmd
.PHONY: clean
# Add cmd at the end of the line below
all: wc uniq tr grep

CXX = g++-10
CFLAGS = -Wall -std=c++17 -O3 -DNDEBUG -Icommon
CFLAGS += common/main.cpp
DEPS = common/*

# Copy one of these pairs of two lines and paste it at the end of the list
# And change each occurence of the command name (like wc) for cmd
wc: $(DEPS) aggregators/wc/agg.h
	${CXX} ${CFLAGS} -Iaggregators/wc -o ./bin/wc
uniq: $(DEPS) aggregators/uniq/agg.h
	${CXX} ${CFLAGS} -Iaggregators/uniq -o ./bin/uniq
tr: $(DEPS) aggregators/tr/agg.h
	${CXX} ${CFLAGS} -Iaggregators/tr -o ./bin/tr
grep: $(DEPS) aggregators/grep/agg.h
	${CXX} ${CFLAGS} -Iaggregators/grep -o ./bin/grep

# That's it for this file

clean:
	rm -rf ./bin/*
