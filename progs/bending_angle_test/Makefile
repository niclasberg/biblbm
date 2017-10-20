CFLAGS := -I../../src/ 
SRC := bending_angle_test.cpp
OBJ := $(SRC:.cpp=.o)
PROG := bending_angle_test  

all: $(PROG) $(OBJ)

$(PROG) : $(OBJ) 
	g++ $(CFLAGS) $< -o $@
	
%.o: %.cpp
	 g++ $(CFLAGS) -c $< -o $@
	 
clean: 
	rm $(PROG) $(OBJ)