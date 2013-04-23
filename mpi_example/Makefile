# Makefile

TARGET = mpi_example

CC = mpicc

CFLAGS = -finline

INCLUDE_PATH = -I$(GSL_INCDIR)
LIBRARY_PATH = -L$(GSL_LIBDIR)

LIBS = $(GSL_LIB) $(GSL_BLAS_LIB) -lm

OBJS = mpi_example.o

all : $(TARGET)

$(TARGET) : $(OBJS)
	$(CC) $(LDFLAGS) -o $@ $(OBJS) $(LIBRARY_PATH) $(LIBS)

.c.o:
	$(CC) $(CFLAGS) -c $*.c $(INCLUDE_PATH)

clean:
	rm -f $(OBJS) $(TARGET)

