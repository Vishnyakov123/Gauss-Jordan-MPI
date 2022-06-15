NAME = prog

CC = mpicc
FLAGS = -Wall -Werror -Wextra -pedantic
LIBS = -lm
HEADER = header.h
RM = rm -f
LIST =	main.c jordan.c input.c

OBJ = $(LIST:.c=.o)

all: $(NAME)

%.o : %.c $(HEADER) Makefile
	$(CC) $(FLAGS) -c $< -o $@

$(NAME): $(OBJ) $(HEADER)
	$(CC) $(FLAGS) -o $(NAME) $(OBJ) $(LIBS)

clean:
	$(RM) $(OBJ) $(NAME)

#shamil is good