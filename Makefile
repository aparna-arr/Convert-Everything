CXX := g++
CXXFLAGS := -Wall -std=c++11 -pthread -g
OBJECTS := Main.o Utils.o Files.o

convert-everything: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $@

Main.o: Main.cpp Utils.h Files.h
	$(CXX) $(CXXFLAGS) -c Main.cpp

Utils.o: Utils.cpp Utils.h
	$(CXX) $(CXXFLAGS) -c Utils.cpp 

Files.o: Files.cpp Files.h
	$(CXX) $(CXXFLAGS) -c Files.cpp

clean: 
	rm -f $(OBJECTS) convert-everything
