#include "logicome_profiler.h"

int main(int argc, char* argv[]){
  if(argc != 5){
    cout << "The number of argument is invalid." << endl;
    exit(1);
  }

  LogicomeProfiler temp_logicome_profiler;
  temp_logicome_profiler.SetParameters(argv);
  temp_logicome_profiler.Run();
}
