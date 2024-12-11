CC = gcc
CFLAGS = -Wall -g
TARGET = bin/main

SRCS = code/main.c code/read.c code/hash.c code/blast.c
OBJS = $(SRCS:.c=.o)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS)

clean:
	rm -f $(TARGET) $(OBJS)