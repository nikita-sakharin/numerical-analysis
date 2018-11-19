CC=g++
CFLAGS=-std=c++17 -fexceptions -m64 -Wall -Werror -Weffc++ -Wextra -pedantic -Wfatal-errors -pedantic-errors -Winit-self -Wzero-as-null-pointer-constant -Wnon-virtual-dtor -Winline -Wmissing-declarations -Wmissing-include-dirs -Wundef -Wredundant-decls -Wfloat-equal -Wmain -Wunreachable-code -Wshadow -Wcast-align -Wswitch-default -Wswitch-enum
LDFLAGS=-m64
LIBS=-lm
SOURCES=LW_7.cpp
EXECUTABLE=$(SOURCES:%.cpp=%)

TARGET=Debug
ifeq ($(TARGET), Debug)
	CFLAGS+=-Og -g
else ifeq ($(TARGET), Release)
	CFLAGS+=-fomit-frame-pointer -fexpensive-optimizations -flto -O3
	LDFLAGS+=-O3 -flto -s
endif

all: $(EXECUTABLE)

$(EXECUTABLE): %: %.cpp
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -frd $(EXECUTABLE)
