CIBLE = main
SRCS =  Main.cpp Camera.cpp Mesh.cpp GLProgram.cpp GLShader.cpp GLError.cpp
OPENGL_PATH = /usr/lib/nvidia-367 # change this for your own environment
LIBS = -L$(OPENGL_PATH) -lglut -lGLU -lGL -lGLEW -lm -lpthread

CC = g++
CPP = g++

FLAGS = -Wall -pthread -g -std=c++11 -O3

CFLAGS = $(FLAGS)
CXXFLAGS = $(FLAGS)

OBJS = $(SRCS:.cpp=.o)   

$(CIBLE): $(OBJS)
	g++ $(LDFLAGS) -o $(CIBLE) $(OBJS) $(LIBS)
clean:
	rm -f  *~  $(CIBLE) $(OBJS)

Camera.o: Camera.cpp Camera.h Vec3.h
Mesh.o: Mesh.cpp Mesh.h Vec3.h
GLError.o: GLError.cpp GLError.h Exception.h 
GLShader.o: GLShader.cpp GLShader.h GLError.h
GLProgram.o: GLProgram.cpp GLProgram.h GLShader.h GLError.h Exception.h
Main.o: Main.cpp Vec3.h Camera.h Mesh.h GLProgram.h Exception.h
