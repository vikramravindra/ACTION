g++ -Wall -c CoinedPValues.cpp
g++ -Wall Example_CoinedPValues_v1.cpp CoinedPValues.o -o Example_CoinedPValues_v1 -lm
g++ -Wall Example_CoinedPValues_v2.cpp CoinedPValues.o -o Example_CoinedPValues_v2  -lm