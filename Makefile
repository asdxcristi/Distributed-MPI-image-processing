MPI=mpic++
FLAGS=-std=c++11 -Wall -Wextra
BIN_NAME=filtru

all: build

build: $(BIN_NAME)

$(BIN_NAME):$(BIN_NAME).cpp
	$(MPI) -o $(BIN_NAME) $(FLAGS) $(BIN_NAME).cpp

clean:
	rm -f $(BIN_NAME)

