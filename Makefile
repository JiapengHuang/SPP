# Makefile

TARGET = mkspeckles

CC = gcc

CFLAGS = -g -O0

INCLUDE_PATH = 
LIBRARY_PATH = 

LIBS = -lgsl -lgslcblas -lhdf5 -lz -lm

OBJS = mkspeckles.o

all : $(TARGET)

$(TARGET) : $(OBJS)
	$(CC) $(LDFLAGS) -o $@ $(OBJS) $(LIBRARY_PATH) $(LIBS)

.c.o:
	$(CC) $(CFLAGS) -c $*.c $(INCLUDE_PATH)

clean:
	rm -f $(OBJS) $(TARGET)

