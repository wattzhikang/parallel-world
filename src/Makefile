BIN = ../bin
GEN = $(BIN)/generated-files

MAIN := main

CLASSES :=

OBJS := $(foreach CLASS,$(CLASSES),$(GEN)/$(CLASS).o)
OBJS += $(MAIN).o

SRCS := $(foreach CLASS,$(CLASSES),$(CLASS).cpp)
SRCS += $(MAIN).cpp

parallel-world : $(OBJS)
	$(CXX) $(OBJS) -o $(BIN)/$@

$(OBJS) : $(SRCS)
	$(CXX) -c $(SRCS) -o $@

.PHONY : doc

doc : Doxyfile
	doxygen $<