GRAPH = gnuplot

CC        = clang
LDFLAGS   = -O2
ifeq ($(CC),clang)
  CFLAGS  = -Weverything -Wextra -pedantic $(LDFLAGS)
else
  CFLAGS  = -Wall -Wextra -pedantic $(LDFLAGS)
endif
LDLIBS    = $(shell gsl-config --libs)

.SUFFIXES:
.SUFFIXES:  .o .c .res .gp .png

.PHONY: clean veryclean all

main: main.o two-cubes-integrand2.o
	${CC} $(LFLAGS) -o $@ $^ $(LDLIBS)

res: mc-pi
	./mc-pi > mc_pi.res

mc-pi1.png: mc-pi1.gp res
	$(GRAPH) mc-pi1.gp

mc-pi2.png: mc-pi2.gp res
	$(GRAPH) mc-pi2.gp

clean:
	rm -f *.o *~ mc-pi

veryclean: clean
	rm -f res mc-pi1.png mc-pi2.png
