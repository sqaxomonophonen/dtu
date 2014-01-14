CC = clang

PKGS = sdl2 gl glew

CFLAGS = -m64 -O3 -Wall $(shell pkg-config $(PKGS) --cflags)
LINK_SDL = $(shell pkg-config $(PKGS) --libs)
LINK = $(LINK_SDL) -lm
OBJS = main.o psys.o

all: main

main.o: main.c psys.h
	$(CC) $(CFLAGS) -c main.c

psys.o: psys.c psys.h
	$(CC) $(CFLAGS) -c psys.c

main: $(OBJS) Makefile
	$(CC) $(CFLAGS) -o main \
		$(OBJS) \
		$(LINK)

clean:
	rm -f $(OBJS) main
