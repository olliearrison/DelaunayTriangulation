APP_NAME=triangle

OBJS=triangle.o

CXX = g++
CXXFLAGS = -Wall -O3 -std=c++17 -m64 -I. -fopenmp -Wno-unknown-pragmas

SANITIZE_FLAGS = -fsanitize=address,undefined -g -O1 -fno-omit-frame-pointer

all: $(APP_NAME)

$(APP_NAME): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

sanitize: CXXFLAGS := $(CXXFLAGS) $(SANITIZE_FLAGS)
sanitize: clean $(APP_NAME)

%.o: %.cpp %.h
	$(CXX) $(CXXFLAGS) -c $<

clean:
	/bin/rm -rf *~ *.o $(APP_NAME) *.class