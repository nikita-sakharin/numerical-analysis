CC=g++
CFLAGS=-std=c++17 -Werror -pedantic -Wall -Wextra -O3 -Os
LDFLAGS=
LIBS=-lm
SOURCES=LW_3_1.cpp LW_3_2.cpp LW_3_3.cpp LW_3_4.cpp LW_3_5.cpp
EXECUTABLE=$(SOURCES:%.cpp=%)

all: $(EXECUTABLE)

$(EXECUTABLE): %: %.cpp
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -frd $(EXECUTABLE)
