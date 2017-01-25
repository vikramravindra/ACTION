OBJECTS = CoinedPValues.o
SYSLIBS = -lm 
LFLAGS = -static
CC = g++  
CFLAGS = -O3     -I .

#PGROUP = (0==TMP_Zaykin,1==Fisher,2==Fisher TAU=1)


Example_CoinedPValues_v2 :  Example_CoinedPValues_v2.cpp  $(OBJECTS) CoinedPValues.mak
	    $(CC) -Wall $(LFLAGS) $(CFLAGS) Example_CoinedPValues_v2.cpp $(OBJECTS) -o Example_CoinedPValues_v2 $(SYSLIBS)

Example_CoinedPValues_v1 :  Example_CoinedPValues_v1.cpp  $(OBJECTS) CoinedPValues.mak
	    $(CC) -Wall $(LFLAGS) $(CFLAGS) Example_CoinedPValues_v1.cpp $(OBJECTS) -o Example_CoinedPValues_v1 $(SYSLIBS)

CoinedPValues.o : CoinedPValues.h CoinedPValues.cpp CoinedPValues.mak
		$(CC) -Wall $(CFLAGS)  -c  CoinedPValues.cpp

 