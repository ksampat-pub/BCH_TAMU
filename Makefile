

LIBDIR  =
INCDIR  =
EXEC    = BCHQ8S2T15R2.out
SRC     = main.c \
	preamble.c \
	gf_arith.c \
	decoder.c \
	gauss_solv.c \
	determinant.c 
                
OBJ     = $(SRC:%.c=%.o)
CFLAGS  = -I$(INCDIR) -L$(LIBDIR)
CC      = cc -g 
# $(CFLAGS)
 
KEEP_STATE=
# 
# 
$(EXEC):$(OBJ) global.h
#       $(CC) $(OBJ) $(LLIBS) -lm -o $(EXEC)
	$(CC)  $(OBJ) -o $(EXEC) -lm

