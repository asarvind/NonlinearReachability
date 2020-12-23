incPath1 = /usr/local/include/
incPath2 = /usr/include/
incPath3 = src/

executable = src/pywrite/simulate

all: src/ui.py
	sudo rm -rf src/pywrite
	python src/ui.py

compile10: src/reachability.cpp
	g++-10 -std=c++11 -fopenmp -I $(incPath1) -I $(incPath2) -I $(incPath3) src/main.cpp -o $(executable)

compile9: src/reachability.cpp
	g++-9 -std=c++11 -fopenmp -I $(incPath1) -I $(incPath2) -I $(incPath3) src/main.cpp -o $(executable)

rand: src/test.cpp
	g++-10 -std=c++11 -fopenmp -I $(incPath1) -I $(incPath2) -I $(incPath3) src/test.cpp -o $(executable)

execute: $(executable)
	time ./$(executable)

