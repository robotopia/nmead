CC      = gcc

LDLIBS  =
OPTIM   = -O2

## Uncomment these to turn on fsanitize option
#LDLIBS = -lasan
#OPTIM   = -O2 -fsanitize=address 

CFLAGS  = -Wall -Wextra $(OPTIM)

LIBRARY = nmead
LIBFILE = lib$(LIBRARY).a
HDRFILE = $(LIBRARY).h
OBJS = nmead.o

all: $(LIBFILE)

$(LIBFILE): $(OBJS)
	ar rcs $@ $^

install:
	cp $(LIBFILE) /usr/local/lib
	cp $(HDRFILE) /usr/local/include

clean:
	$(RM) *.o $(LIBFILE)
