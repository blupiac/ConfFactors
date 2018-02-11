TARGET = confFact

SRCDIR   = src
OBJDIR   = obj

SOURCES  := $(wildcard $(SRCDIR)/*.cpp)
OBJECTS  := $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)
rm = rm -f

LFLAGS = -lglut -lGLU -lGL -lGLEW -lm -lpthread

CC = g++
CPP = g++
LINKER = g++

FLAGS = -Wall -pthread -g -std=c++11 -O3

CFLAGS = $(FLAGS)
CXXFLAGS = $(FLAGS)



$(TARGET): $(OBJECTS)
	$(LINKER) $(OBJECTS) $(LFLAGS) -o $@

$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

.PHONY: clean
clean:
	$(rm) $(OBJECTS)

.PHONY: remove
remove: clean
	$(rm) $(BINDIR)/$(TARGET)