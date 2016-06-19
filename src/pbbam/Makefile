CFLAGS=-fPIC -g -O2 -Wall -I../../inst/include/ --std=c++11 -Wno-deprecated-register -Wno-unused-local-typedefs

# Should be equivalent to your list of C files, if you don't build selectively
#SRC=$(wildcard *.cpp)
CPP_FILES := $(wildcard *.cpp)
OBJ_FILES := $(CPP_FILES:.cpp=.o)

XML_FILES := $(wildcard pugixml/*.cpp)
XML_OBJ_FILES := $(CPP_FILES:.cpp=.o)

#test: $(SRC)
#	gcc -o test $(SRC) $(CFLAGS) $(LIBS)

# $@ - The name of the target
# $^ - All dependencies
# %.o - wildcard matching,
# See http://www.cprogramming.com/tutorial/makefiles_continued.html

libpbbam.a: $(OBJ_FILES) pugixml.o
	$(AR) -r $@ $(OBJ_FILES) pugixml.o


%.o: %.cpp
	g++ $(CFLAGS) $(LIBS) -c -o $@ $<

pugixml.o: $(XML_FILES)
	echo $(XML_FILES)
	g++ $(CFLAGS) -c $< -o $@

clean:
	rm *.o
	rm *.a
