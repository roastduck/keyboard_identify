plainio : analyze.cpp plainio.cpp lodepng.cpp
	g++ analyze.cpp plainio.cpp lodepng.cpp -o keyboard -O2

yamlio : analyze.cpp yamlio.cpp lodepng.cpp
	g++ analyze.cpp yamlio.cpp lodepng.cpp -o keyboard -lboost_system -lboost_filesystem -lyaml-cpp -O2
