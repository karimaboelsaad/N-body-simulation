CXX := c++
CXXFLAGS := -std=c++20 -O2 -Wall -Wextra -Wpedantic -Iinclude

TARGET := nbody
SOURCES := $(wildcard src/*.cpp)
HEADERS := $(wildcard include/*.hpp)

$(TARGET): $(SOURCES) $(HEADERS)
	$(CXX) $(CXXFLAGS) $(SOURCES) -o $(TARGET)

.PHONY: clean
clean:
	$(RM) $(TARGET)
