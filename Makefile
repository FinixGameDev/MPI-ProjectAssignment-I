CC=mpicc
CFLAGS=-O2 -std=c99
TARGET=fox

all: $(TARGET)

$(TARGET): main.c
	$(CC) $(CFLAGS) -o $(TARGET) main.c

clean:
	rm -f $(TARGET) *.o
