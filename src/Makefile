BIN = ../bin
GEN = $(BIN)/generated-files

MAIN := main

CLASSES := DiamondSquare

OBJS := $(foreach CLASS,$(CLASSES),$(GEN)/$(CLASS).o)
OBJS += $(GEN)/$(MAIN).o

SRCS := $(foreach CLASS,$(CLASSES),$(CLASS).cpp)
SRCS += $(MAIN).cpp

parallel-world : $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(BIN)/$@

$(OBJS) : $(GEN)/%.o : %.cpp
	$(CXX) $< $(CXXFLAGS) -c -o $@

.PHONY : doc

doc : Doxyfile
	doxygen $<