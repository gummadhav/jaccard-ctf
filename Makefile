include config.mk

FCXX=$(CXX) $(CXXFLAGS)

all: jaccard

jaccard: jaccard.cxx Makefile config.mk
	$(FCXX) $< -o $@ $(INCLUDE_PATH) $(LIB_PATH) $(LIBS)

clean:
	rm -f jaccard
