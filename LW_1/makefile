CC=g++
CFLAGS=-c -std=c++14 -Werror -pedantic -Wall -Wextra -g
LDFLAGS=
LIBS=-lm
SOURCES=LW_1.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=LW_1

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) $(LIBS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -frd $(OBJECTS) $(EXECUTABLE)
