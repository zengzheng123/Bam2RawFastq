OBJS = Bam2RawFastq.o
CC = g++
CFLAGS = -Wall -c
LFLAGS = -Wall -o3
GZFLAGS = -I./gzstream/include/ -L./gzstream/lib/ -lgzstream -lz
BAMFLAGS = -I./bamtools-master/include/ -L./bamtools-master/lib/ -lbamtools

Bam2RawFastq: $(OBJS)
	$(CC) $(OBJS) $(LFLAGS) -o Bam2RawFastq $(GZFLAGS) $(BAMFLAGS)

Bam2RawFastq.o: Bam2RawFastq.cpp
	$(CC) $(CFLAGS) $(GZFLAGS) $(BAMFLAGS) Bam2RawFastq.cpp

clean:
	rm -f *.o Bam2RawFastq