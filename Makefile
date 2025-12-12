# Compiler and flags
CC = gcc
CXX = g++
CFLAGS = -Wall -Wextra -g
CXXFLAGS = -Wall -Wextra -g -std=c++17

# Directories
SRC_DIR = src
INC_DIR = include
OBJ_DIR = obj
BIN_DIR = bin

# Define source files, object files, and executable
SOURCES = $(wildcard $(SRC_DIR)/*.c) $(wildcard $(SRC_DIR)/*.cc)
OBJECTS = $(patsubst $(SRC_DIR)/%.c,$(OBJ_DIR)/%.o,$(filter %.c,$(SOURCES))) \
          $(patsubst $(SRC_DIR)/%.cc,$(OBJ_DIR)/%.o,$(filter %.cc,$(SOURCES)))
EXECUTABLE = rareEventSim

# Default target
all: $(BIN_DIR)/$(EXECUTABLE)

# Rule to create object files from C source files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c | $(OBJ_DIR)
	$(CC) $(CFLAGS) -I$(INC_DIR) -c $< -o $@

# Rule to create object files from C++ source files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cc | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -I$(INC_DIR) -c $< -o $@

# Rule to link object files into an executable
$(BIN_DIR)/$(EXECUTABLE): $(OBJECTS) | $(BIN_DIR)
	$(CXX) $(OBJECTS) -o $@

# Create object directory if it doesn't exist
$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

# Create binary directory if it doesn't exist
$(BIN_DIR):
	mkdir -p $(BIN_DIR)

# Clean target
clean:
	rm -rf $(OBJ_DIR)/*.o $(BIN_DIR)/$(EXECUTABLE)

# Phony targets
.PHONY: all clean
