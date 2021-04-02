IDIR=include
SDIR=src
PDIR=plots
DEPS=$(wildcard $(IDIR)/*.h)
SOURCESC=$(wildcard $(SDIR)/*.cpp)
CC=g++
EXECUTABLE=SpeckleSim
LIBS= -lgsl -lgslcblas

OBJECTSC=$(SOURCESC:.cpp=.o)

CFLAGS=-I$(IDIR) `root-config --cflags --glibs`

all: $(EXECUTABLE) $(SOURCESC) $(DEPS)

$(EXECUTABLE): $(OBJECTSC)
	$(CC) -o $@ $(OBJECTSC) $(CFLAGS) $(LIBS)

.cpp.o:
	$(CC) -c $(CFLAGS) $< -o $@


.PHONY: clean

clean:
	rm -f $(PDIR)/* $(SDIR)/*.o $(IDIR)/*~ $(EXECUTABLE) $(SDIR)/*~ /*~ *~ *.csv *.pdf
