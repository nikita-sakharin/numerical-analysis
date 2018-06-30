CC=g++
CFLAGS=-std=c++17 -Werror -pedantic -Wall -Wextra -O3 -Os
LDFLAGS=
LIBS=-lm
SOURCES=LW_4_1.cpp LW_4_2.cpp
EXECUTABLE=$(SOURCES:%.cpp=%)

all: $(EXECUTABLE)

$(EXECUTABLE): %: %.cpp
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -frd $(EXECUTABLE)
