incPath1 = /usr/local/include/
incPath2 = /usr/include/
incPath3 = src/

executable = simulate

all: src/ui.py
	rm -rf src/pywrite
	python src/ui.py

# compile10: src/reachability.cpp
# 	g++-10 -std=c++11 -fopenmp -I $(incPath1) -I $(incPath2) -I $(incPath3) main.cpp -o $(executable)
