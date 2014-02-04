CC = clang

PKGS = sdl2 gl glew libpng

CFLAGS = -m64 -O3 -Wall $(shell pkg-config $(PKGS) --cflags)
LINK_SDL = $(shell pkg-config $(PKGS) --libs)
LINK = $(LINK_SDL) -lm
OBJS = main.o psys.o rquad.o mud.o log.o a.o

all: main

main.o: main.c psys.h
	$(CC) $(CFLAGS) -c main.c

a.o: a.c a.h
	$(CC) $(CFLAGS) -c a.c

psys.o: psys.c psys.h mud.h rquad.h scratch.h
	$(CC) $(CFLAGS) -c psys.c

rquad.o: rquad.c rquad.h a.h
	$(CC) $(CFLAGS) -c rquad.c

mud.o: mud.c mud.h log.h
	$(CC) $(CFLAGS) -c mud.c

log.o: log.c log.h
	$(CC) $(CFLAGS) -c log.c

main: $(OBJS) Makefile
	$(CC) $(CFLAGS) -o main \
		$(OBJS) \
		$(LINK)

clean:
	rm -f *.o main rquad_unit_test rquad_vis_test


# tests:

tests: rquad_unit_test rquad_vis_test

rquad_unit_test: rquad.c rquad.h a.o
	$(CC) -m64 -O3 -Wall -DRQUAD_UNIT_TEST rquad.c a.o -lm -o rquad_unit_test

rquad_vis_test: rquad.c rquad.h a.o
	$(CC) -m64 -O3 -Wall -DRQUAD_VIS_TEST rquad.c a.o -lm -o rquad_vis_test
