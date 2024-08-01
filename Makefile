# Variables
CC = gcc-14
GSL_PREFIX = /opt/homebrew/opt/gsl
CFLAGS = -I$(GSL_PREFIX)/include -O3
LDFLAGS = -L$(GSL_PREFIX)/lib -lgsl -lgslcblas -lm
TARGET = KerrAngleConversion
SRCS = $(wildcard *.c)
OBJS = $(SRCS:.c=.o)

# Default target
all: $(TARGET)

# Link
$(TARGET): $(OBJS)
	$(CC) -o $@ $^ $(LDFLAGS)

# Compile
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Clean up
clean:
	rm -f $(OBJS) $(TARGET)

# Phony targets
.PHONY: all clean