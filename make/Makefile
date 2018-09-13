FC = gfortran
FFLAGS = -O3 -fbounds-check -fbacktrace -ffpe-trap=zero,overflow -Wconversion -Wall -fcheck=all -g -ggdb

OBJ_DIR = ../obj
SRC_DIR = ../src
BIN_DIR = ../bin

_OBJ = types.o math.o io.o core.o
OBJ = $(patsubst %, $(OBJ_DIR)/%, $(_OBJ))

all: AADG3 test_AADG3

AADG3: $(OBJ) $(OBJ_DIR)/run.o | $(BIN_DIR)
	$(FC) $(FFLAGS) $^ -o $(BIN_DIR)/$@

test_AADG3: $(OBJ) $(OBJ_DIR)/test.o
	$(FC) $(FFLAGS) $^ -o $(BIN_DIR)/$@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90 | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(BIN_DIR):
	mkdir $(BIN_DIR)

$(OBJ_DIR):
	mkdir $(OBJ_DIR)

.PHONY : clean

clean :
	rm -f $(OBJ_DIR)/*.o *.mod $(BIN_DIR)/AADG3 $(BIN_DIR)/test_AADG3
