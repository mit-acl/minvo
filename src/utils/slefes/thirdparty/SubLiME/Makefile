# -- makefile
#
#  Use "make"         to create the library
#      "make example" to create the examples
#      "make clean"   to clear all the binary output
#
#    if you  encounter the message: "Missing separator..."
#    make sure that all indented lines are indented by using a TAB
#    rather than blanks
#
#  Tested under Linux/Solaris (June 2003)
#  Author: Xiaobin Wu (xwu@cise.ufl.edu)

CC = gcc
CFLAGS = -g -c -Wall
LDLIB = -lm

# Library linking 
LDLIB = -lpthread -L/usr/X11R6/lib -lglut  -lGL -lGLU \
		         -lXi -lXmu -lXext -lX11 -lm

# ---- Files --- 
output	: SubLiME

SUBLIME	= uniSlefe.o tpSlefe.o bspSlefe.o

SubLiME:    $(SUBLIME)
	ld -r -o libSubLiME.a $(SUBLIME)

example:   uniexample.c SubLiME bspexample.c
	$(CC) -o uniexample -g $(SUBLIME) uniexample.c  $(LDLIB) -LSubLiME
	$(CC) -o bspexample -g $(SUBLIME) bspexample.c  $(LDLIB) -LSubLiME

.c.o:
	$(CC) $(CFLAGS) $*.c

clean:
	rm -f *.o libSubLiME.a uniexample bspexample
