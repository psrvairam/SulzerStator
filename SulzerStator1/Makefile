# Use Prefix to define the home folder of the source code.
# It can be different from the folder in which you want to compile and run the mesh generator. 
# In the current directory ./ you only need to have the main.cpp and this Makefile

#PREFIX = /home/renep/Projects/hybridMeshGen/src/
#PREFIX = /home/gjoterorodrigu/Programs/MeshGen/src/
PREFIX = ../src/

# Definitions needed to compire Triangle
TRILIBDEFS = -DTRILIBRARY -DANSI_DECLARATORS
#CSWITCHES = -O3 -DLINUX -I/usr/X11R6/include -L/usr/X11R6/lib
CSWITCHES = -DLINUX -I/usr/X11R6/include -L/usr/X11R6/lib

# Define compiler and optimizer's flags
CC = g++
CFLAGS = -DWITH_OPENGL
LIBSOPENGL = -lGLU -lGL -lglut



# List of objects
OBJ_SRC = Param.o\
			point.o\
			lines.o\
			functions.o\
			struct.o\
			unstruct.o\
			transform.o\
			meshTools.o\
			triangle.o
			

OBJ = $(addprefix $(PREFIX), $(OBJ_SRC)) 

OBJ_SRC2 = Param.o
			

OBJ2 = $(addprefix $(PREFIX), $(OBJ_SRC2))
all:	bladeGen

bladeGen: main.cpp $(OBJ) $(PREFIX)OpenGlutVisu.h $(PREFIX)ArcBall.h
		$(CC) $(CFLAGS) -o bladeGen $(OBJ) main.cpp -lm $(LIBSOPENGL) -I$(PREFIX)

$(PREFIX)%.o: $(PREFIX)%.cpp	$(PREFIX)%.h
	$(CC) $(CFLAGS) -c $< -o $@
	
$(PREFIX)triangle.o: $(PREFIX)triangle.c $(PREFIX)triangle.h
	$(CC) $(CSWITCHES) $(TRILIBDEFS) -c $< -o $@



clean:
		rm bladeGen $(OBJ)
