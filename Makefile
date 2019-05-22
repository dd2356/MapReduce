CC = mpic++
SRC = $(wildcard src/*.cpp)
OBJ = $(SRC:.cpp=.o)
BIN = ./bin
TARGET = mapreduce.out

LDFLAGS = -lm -fopenmp
CFLAGS = -I./include -g -Wall -O3 --std=c++11 -fopenmp

all: dir $(BIN)/$(TARGET)

dir: ${BIN}

${BIN}:
	mkdir -p $(BIN)

%.o: %.cpp
	$(CC) $(CFLAGS) -c -o $@ $<

$(BIN)/$(TARGET): $(OBJ)
	$(CC) $(LDFLAGS) -o $@ $^

.PHONY: clean
clean:
	rm -f $(OBJ) $(BIN)/$(TARGET)

