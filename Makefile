APP_NAME=triangle

OBJS=triangle.o

CXX = g++

BASE_FLAGS = -Wall -O3 -std=c++17 -m64 -I. -Wno-unknown-pragmas
OPENMP_FLAGS = -fopenmp
CXXFLAGS = $(BASE_FLAGS) $(OPENMP_FLAGS)

SANITIZE_FLAGS = -fsanitize=address,undefined -g -O1 -fno-omit-frame-pointer

all: $(APP_NAME)

$(APP_NAME): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

%.o: %.cpp %.h
	$(CXX) $(CXXFLAGS) -c $<

sanitize: CXXFLAGS := $(BASE_FLAGS) $(SANITIZE_FLAGS)
sanitize: clean $(APP_NAME)

clean:
	/bin/rm -rf *~ *.o $(APP_NAME) *.class