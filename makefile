CXX=g++-mp-4.9
CXX=clang++
CXX=g++
CC=gcc
ASA_VER := ASA_CG-3.0
#IPOPT_INCLUDE := -DHAVE_IPOPT -I$(HOME)/local/include/coin
#IPOPT_LIBS := -L${HOME}/local/lib -lipopt -llapack -lm -lcoinmumps -lpthread -lblas -lcoinhsl -lcoinmetis
INCLUDE := -I. -I./eigen -I/usr/include/eigen3 -I/opt/local/include/eigen3 -I$(ASA_VER)/ $(IPOPT_INCLUDE)
CFLAGS := -O3 -g
CXXFLAGS := -O0 -g -Wall -Wextra -pedantic-errors -std=c++11 $(INCLUDE) #-fopenmp
CXXFLAGSTEST := -g -Wall -Wextra -pedantic-errors -std=c++11 $(INCLUDE) -Igtest-1.7.0/include #-fopenmp 
DEPS := src/*.h src/*.cpp

.PHONY: main test

main: src/main.cpp $(DEPS) src/lbfgs.o  $(ASA_VER)/asa_cg.o 
	$(CXX) $(CXXFLAGS) -o main src/main.cpp  src/lbfgs.o $(ASA_VER)/asa_cg.o $(IPOPT_LIBS)

test: src/unittests.cpp $(DEPS) src/lbfgs.o $(ASA_VER)/asa_cg.o 
	$(CXX) $(CXXFLAGSTEST) -o test src/unittests.cpp libgtest.a -lpthread src/lbfgs.o $(ASA_VER)/asa_cg.o 

testdual: src/testdual.cpp unsupported/Eigen/* unsupported/Eigen/src/*
	$(CXX) $(CXXFLAGSTEST) -o testdual src/testdual.cpp libgtest.a -lpthread

clean:
	rm main test $(ASA_VER)/asa_cg.o src/lbfgs.o

# google-testing-framework
install-gtest:
	rm -f gtest-1.7.0.zip
	rm -fR gtest-1.7.0
	wget -O gtest-1.7.0.zip https://googletest.googlecode.com/files/gtest-1.7.0.zip
	unzip gtest-1.7.0.zip
	g++ -Igtest-1.7.0/include -Igtest-1.7.0 -c "gtest-1.7.0/src/gtest-all.cc" 
	ar -rv libgtest.a gtest-all.o
	rm -f gtest-1.7.0.zip

# eigen library
install-eigen:
	wget -c http://bitbucket.org/eigen/eigen/get/3.2.2.tar.bz2 -O eigen-3.2.2.tar.bz2
	bunzip2 eigen-3.2.2.tar.bz2
	tar xvf eigen-3.2.2.tar
	mv eigen-eigen-* eigen 
	rm -Rf eigen-3.2.2.tar
	rm -Rf eigen-3.2.2.tar.bz2

# gpl'd ASA_CG code
install-asa:
	wget http://users.clas.ufl.edu/hager/papers/CG/Archive/$(ASA_VER).tar.gz
	tar xvf $(ASA_VER).tar.gz
	rm -Rf $(ASA_VER).tar.gz

install: install-gtest install-eigen install-asa
	ln -s . CppNumericalSolvers
	echo "installed."
