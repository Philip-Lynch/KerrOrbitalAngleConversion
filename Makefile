# Variables
CC = gcc-14
GSL_PREFIX = /opt/homebrew/opt/gsl
CFLAGS = -I$(GSL_PREFIX)/include -O3 -Iinclude
LDFLAGS = -L$(GSL_PREFIX)/lib -lgsl -lgslcblas -lm

SRC_DIR = src
BUILD_DIR = build
BIN_DIR = bin
INCLUDE_DIR = include


SRCS = $(wildcard $(SRC_DIR)/*.c)
OBJS = $(patsubst $(SRC_DIR)/%.c, $(BUILD_DIR)/%.o, $(SRCS))

TARGET = $(BIN_DIR)/KerrAngleConversion

# Default target
all: $(TARGET)

# Link
$(TARGET): $(OBJS)
	$(CC) -o $@ $^ $(LDFLAGS)

# Compile
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(BUILD_DIR)
	$(CC) $(CFLAGS) -c -o $@ $<

# Clean up
clean:
	rm -rf $(BUILD_DIR)/*.o $(TARGET)

# Phony targets
.PHONY: all clean