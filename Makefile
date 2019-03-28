CXXFLAGS = -O3

LogicomeProfiler: main.cpp logicome_profiler.cpp

	$(CXX) $(CXXFLAGS) -o LogicomeProfiler main.cpp logicome_profiler.cpp
