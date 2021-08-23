IDIR =include
CC=g++
CFLAGS=-I$(IDIR)

OUT_DIR=out
IN_DIR=src

_DEPS = rand_norm_matrix.h lanczos_algo.h lag_mult_estimator.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = rand_norm_matrix.o lanczos_algo.o lag_mult_estimator.o constrained_lso.o
OBJ = $(patsubst %,$(OUT_DIR)/%,$(_OBJ))


$(OUT_DIR)/%.o: $(IN_DIR)/%.cc $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(OUT_DIR)/constrained_lso: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

.PHONY: clean

clean:
	rm -f $(OUT_DIR)/** 

#*~ core $(INCDIR)/*.o