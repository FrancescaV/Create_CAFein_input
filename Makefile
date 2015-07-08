OBJS = dma.o SteffenInterp.o
CC = g++
DEBUG = -g
CFLAGS = -Wall -Wunused -O6 -c $(DEBUG)
LFLAGS = -Wall -Wunused -O6 $(DEBUG)
EXECUTABLE = convertModel
LIBS=-lgsl -lgslcblas -lm

$(EXECUTABLE): $(EXECUTABLE).o $(OBJS)
	$(CC) $(LFLAGS) $(EXECUTABLE).o $(OBJS) $(LIBS) -o $(EXECUTABLE)
	
dma.o: dma.h dma.c
	$(CC) $(CFLAGS) dma.c

SteffenInterp.o: SteffenInterp.h SteffenInterp.c
	$(CC) $(CFLAGS) SteffenInterp.c

.PHONY = clean

clean:
	\rm *.o $(EXECUTABLE)
