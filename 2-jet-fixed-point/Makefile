CC=gcc
CFLAGS=-c -O3
LFLAGS=-O3 -lm


SRC=ode.c poinca.c main.c 
HDR=ode.h poinca.h
OBJ=$(SRC:.c=.o)
OBJS=$(OBJ:.cpp=.o)
TGT=main.exe

all: $(TGT)

$(TGT): $(OBJS)
	$(CC) $^ -o $@ $(LFLAGS)

.c.o:
	$(CC) $< -o $@ $(CFLAGS)
%.o: %.c %.h
	$(CC) $< -o $@ $(CFLAGS)


clean:
	-rm -f $(OBJS) $(TGT)
