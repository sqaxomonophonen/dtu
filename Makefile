CC = clang

PKGS = sdl2 gl glew libpng

CFLAGS = -m64 -O3 -Wall $(shell pkg-config $(PKGS) --cflags)
LINK_SDL = $(shell pkg-config $(PKGS) --libs)
LINK = $(LINK_SDL) -lm
OBJS = main.o psys.o mud.o log.o

all: main

main.o: main.c psys.h
	$(CC) $(CFLAGS) -c main.c

psys.o: psys.c psys.h mud.h
	$(CC) $(CFLAGS) -c psys.c

mud.o: mud.c mud.h log.h
	$(CC) $(CFLAGS) -c mud.c

log.o: log.c log.h
	$(CC) $(CFLAGS) -c log.c

main: $(OBJS) Makefile
	$(CC) $(CFLAGS) -o main \
		$(OBJS) \
		$(LINK)

clean:
	rm -f *.o main
