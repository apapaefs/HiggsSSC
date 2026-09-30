# Herwig activation scripts may export a build-time CXX (including stale
# compiler paths and embedded flags).  Keep this standalone analyzer on the
# system C++ driver by default; a command-line `make CXX=...` still overrides
# this assignment when a specific toolchain is required.
CXX := c++

ROOT_CFLAGS := $(shell root-config --cflags)
ROOT_LIBS := $(shell root-config --libs)

CXXFLAGS ?= -O2 -g
CXXFLAGS += -std=c++17 -Wall -Wextra -Wpedantic

TARGET := HwSimPostAnalysis_fourlepton
SOURCE := $(TARGET).cc
FIXTURE := make_synthetic_hw_sim_fourlepton
CORE_TEST := test_fourlepton_core

.PHONY: all fixture check clean

all: $(TARGET)

$(TARGET): $(SOURCE)
	$(CXX) $(CXXFLAGS) $(ROOT_CFLAGS) -o $@ $< $(ROOT_LIBS)

fixture: $(FIXTURE)

$(FIXTURE): $(FIXTURE).cc
	$(CXX) $(CXXFLAGS) $(ROOT_CFLAGS) -o $@ $< $(ROOT_LIBS)

check: $(TARGET) $(CORE_TEST)
	./$(CORE_TEST)

$(CORE_TEST): $(CORE_TEST).cc $(SOURCE)
	$(CXX) $(CXXFLAGS) $(ROOT_CFLAGS) -o $@ $< $(ROOT_LIBS)

clean:
	rm -f $(TARGET) $(TARGET).o $(FIXTURE) $(FIXTURE).o \
		$(CORE_TEST) $(CORE_TEST).o
