# Feel free to add additional target programs in the next line.
Target  = analyse

# No need to change anything after this line, unless you are a
# Makefile guru

ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS  = $(shell root-config --libs)

all: $(Target)

%: %.cpp Analysis.hpp
	$(CXX) -o $@ $< $(ROOTFLAGS) $(ROOTLIBS)

clean:
	$(RM) $(Target)
