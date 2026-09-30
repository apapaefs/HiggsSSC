HERWIG_CONFIG ?= herwig-config
CXX := $(if $(HERWIG_CXX),$(HERWIG_CXX),c++)

CPPFLAGS := $(shell $(HERWIG_CONFIG) --cppflags)
CXXFLAGS ?= -O2 -g
CXXFLAGS += -std=c++17 -fPIC -Wall -Wextra -Wpedantic
LDFLAGS := $(shell $(HERWIG_CONFIG) --ldflags)
LDLIBS := $(shell $(HERWIG_CONFIG) --ldlibs)

ifeq ($(shell uname -s),Darwin)
  PLUGIN_FLAGS := -bundle -undefined dynamic_lookup
else
  PLUGIN_FLAGS := -shared
endif

TARGET := LHEBranchingRatioReweighter.so
SOURCE := LHEBranchingRatioReweighter.cc

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(SOURCE)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PLUGIN_FLAGS) -o $@ $< \
		$(LDFLAGS) $(LDLIBS)

clean:
	rm -f $(TARGET) *.o
