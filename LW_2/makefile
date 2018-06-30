CC=g++
CFLAGS=-std=c++14 -Werror -pedantic -Wall -Wextra -O3
LDFLAGS=
LIBS=-lm
SOURCES=LW_2_1.cpp LW_2_2.cpp
EXECUTABLE=$(SOURCES:%.cpp=%)

all: $(EXECUTABLE)

$(EXECUTABLE): %: %.cpp
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -frd $(EXECUTABLE)


