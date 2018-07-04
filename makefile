CXX=g++-mp-4.9
CXX=clang++
CXX=g++
CC=gcc
ASA_VER := ASA_CG-3.0
GTESTVERSION=1.8.0
EIGENVERSION=3.2.2
EIGENVERSION=3.3.4
INCLUDES=   -I Eigen
#IPOPT_INCLUDE := -DHAVE_IPOPT -I$(HOME)/local/include/coin
#IPOPT_LIBS := -L${HOME}/local/lib -lipopt -llapack -lm -lcoinmumps -lpthread -lblas -lcoinhsl -lcoinmetis
INCLUDE := -I. -IEigen -I$(ASA_VER)/ $(IPOPT_INCLUDE)
CFLAGS := -O3 -g
CXXFLAGS := -O0 -g -Wall -Wextra -pedantic-errors -std=c++11 $(INCLUDE)
CXXFLAGSTEST := -g -Wall -Wextra -pedantic-errors -std=c++11 $(INCLUDE) \
	-I googletest-release-$(GTESTVERSION)/googletest/include
DEPS := src/*.h src/*.cpp

.PHONY: main test testdual

main: src/main.cpp $(DEPS) src/lbfgs.o  $(ASA_VER)/asa_cg.o 
	$(CXX) $(CXXFLAGS) -o main src/main.cpp  src/lbfgs.o $(ASA_VER)/asa_cg.o $(IPOPT_LIBS)

test: src/unittests.cpp $(DEPS) src/lbfgs.o $(ASA_VER)/asa_cg.o 
	$(CXX) $(CXXFLAGSTEST) -o test src/unittests.cpp libgtest.a -lpthread src/lbfgs.o $(ASA_VER)/asa_cg.o 

testdual: src/testdual.cpp unsupported/Eigen/* unsupported/Eigen/src/*
	$(CXX) $(CXXFLAGSTEST) -o testdual src/testdual.cpp libgtest.a -lpthread

clean:
	rm main test $(ASA_VER)/asa_cg.o src/lbfgs.o

# google-testing-framework
libgtest.a:
	-rm -f googletest-release-$(GTESTVERSION).zip
	rm -fR googletest-release-$(GTESTVERSION)
	wget -O googletest-release-$(GTESTVERSION).zip https://github.com/google/googletest/archive/release-$(GTESTVERSION).zip
	unzip googletest-release-$(GTESTVERSION).zip
	g++ -Igoogletest-release-$(GTESTVERSION)/googletest/include -Igoogletest-release-$(GTESTVERSION)/googletest -c "googletest-release-$(GTESTVERSION)/googletest/src/gtest-all.cc" 
	ar -rv libgtest.a gtest-all.o
	rm -f googletest-release-$(GTESTVERSION).zip

# eigen library
install-eigen:
	-rm eigen-$(EIGENVERSION).tar.bz2
	wget -c http://bitbucket.org/eigen/eigen/get/$(EIGENVERSION).tar.bz2 -O eigen-$(EIGENVERSION).tar.bz2
	tar xjvf eigen-$(EIGENVERSION).tar.bz2
	-rm -r Eigen
	mv eigen-eigen-* Eigen 
	rm -Rf eigen-$(EIGENVERSION).tar
	rm -Rf eigen-$(EIGENVERSION).tar.bz2

# gpl'd ASA_CG code
install-asa:
	wget http://users.clas.ufl.edu/hager/papers/CG/Archive/$(ASA_VER).tar.gz
	tar xvf $(ASA_VER).tar.gz
	rm -Rf $(ASA_VER).tar.gz

install: libgtest.a install-eigen install-asa
	-ln -s . CppNumericalSolvers
	echo "installed."
