CFLAGS=-g #-Wall  #-Werror
CC=mpicc
OBJS=main.o
LIBS=-lm
output: ${OBJS}
	${CC} ${CFLAGS} -o output ${LDFLAGS} $^ ${LIBS}
clean:
	rm *.o output
