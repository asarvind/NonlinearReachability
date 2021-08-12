incPath1 = /usr/local/include/
incPath2 = /usr/include/
incPath3 = src/
incPath4 = /opt/homebrew/include/
incPath5 = /opt/homebrew/include/

executable = src/pywrite/simulate

all: src/ui.py
	sudo rm -rf src/pywrite
	python src/ui.py

compile11: src/reachability.cpp
	g++-11 -std=c++20 -fopenmp -I $(incPath4) -I $(incPath5) -I $(incPath3) src/main.cpp -o $(executable)

compile9: src/reachability.cpp
	g++-9 -std=c++11 -fopenmp -I $(incPath1) -I $(incPath2) -I $(incPath3) src/main.cpp -o $(executable)

rand10: src/randmain.cpp
	g++-10 -std=c++11 -fopenmp -I $(incPath1) -I $(incPath2) -I $(incPath3) src/randmain.cpp -o $(executable)

rand9: src/randmain.cpp
	g++-9 -std=c++11 -fopenmp -I $(incPath1) -I $(incPath2) -I $(incPath3) src/randmain.cpp -o $(executable)

execute: $(executable)
	time ./$(executable)

