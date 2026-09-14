CXX := c++
CXXFLAGS := -std=c++20 -O2 -Wall -Wextra -Wpedantic -Iinclude

TARGET := nbody
TEST_TARGET := nbody_tests
SOURCES := $(wildcard src/*.cpp)
TEST_SOURCES := tests/tests.cpp src/gravity.cpp src/quadtree.cpp src/integrators.cpp src/scenarios.cpp
HEADERS := $(wildcard include/*.hpp)

$(TARGET): $(SOURCES) $(HEADERS)
	$(CXX) $(CXXFLAGS) $(SOURCES) -o $(TARGET)

$(TEST_TARGET): $(TEST_SOURCES) $(HEADERS)
	$(CXX) $(CXXFLAGS) $(TEST_SOURCES) -o $(TEST_TARGET)

.PHONY: test clean
test: $(TEST_TARGET)
	./$(TEST_TARGET)

clean:
	$(RM) $(TARGET) $(TEST_TARGET)
