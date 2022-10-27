CXX=g++
CC=g++
CFLG= -O3 -Wall -o3

TEST: Fif.o TestFIF.o FFT.o interp.o 
	$(CXX) TestFIF.o Fif.o FFT.o interp.o -o test -lfftw3l -lm -g #-fPIC
	
TestFIF.o: TestFIF.c Fif.h
	$(CXX) -c TestFIF.c #-fPIC

Fif.o: Fif.c Fif.h
	$(CXX) -c Fif.c #-fPIC
	
FFT.o: FFT.c FFT.h
	$(CXX) -c FFT.c -lfftw3l -lm #-fPIC

interp.o: interp.c interp.h
	$(CXX) -c interp.c -lm #-fPIC

clean:
	rm *.o
