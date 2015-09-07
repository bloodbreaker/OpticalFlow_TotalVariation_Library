OPT = -O0 -g -Wall
LIBS = -lm -lGL -lGLU -lglut
ARGS1 = basename ../images/yos
ARGS = basename ../../Pictures/images/yos
EXECUTABLE = of_frontend 

all: $(EXECUTABLE) #run

frontend: main_frontend.c Makefile
	gcc $(OPT) $(LIBS) -o frontend main_frontend.c

$(EXECUTABLE): *.c Makefile
	gcc $(OPT) $(LIBS) -o of_frontend of_frontend.c

# normal program run
run: of_frontend
	./of_frontend $(ARGS)

# debug in gdb
debug: 
	gdb --eval-command=run   --args ./$(EXECUTABLE) $(ARGS)

# let valgrind see...
debug_val:
	valgrind --db-attach=yes --track-origins=yes ./$(EXECUTABLE) $(ARGS)
	
clean:
	rm $(EXECUTABLE)	
	
	
