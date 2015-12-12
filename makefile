keyboard : keyboard.cpp lodepng.cpp
	g++ keyboard.cpp lodepng.cpp -o keyboard -lboost_system -lboost_filesystem -lyaml-cpp -ggdb -Wall -DDEBUG
