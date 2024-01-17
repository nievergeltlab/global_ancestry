OBJS = StringTokenizer.o Split.o IOHandler.o InferAnc.o

inferanc : $(OBJS)
	g++ -g $(OBJS) main.cpp -o inferanc -lm

StringTokenizer.o : StringTokenizer.cpp StringTokenizer.h
	g++ -c -g StringTokenizer.cpp -o StringTokenizer.o

Split.o : Split.cpp Split.h
	g++ -c -g Split.cpp -o Split.o

IOHandler.o : IOHandler.cpp IOHandler.h
	g++ -c -g IOHandler.cpp -o IOHandler.o

InferAnc.o : InferAnc.cpp InferAnc.h
	g++ -c -g InferAnc.cpp -o InferAnc.o

install:  inferanc
	mv inferanc ../bin
clean:
	rm -f *.o inferanc



