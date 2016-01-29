CFLAGS += -O2 -std=c99 -Werror

all:kseq.h fastool.c
	$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) fastool.c -o fastool

clean:
	rm -f *.o fastool
