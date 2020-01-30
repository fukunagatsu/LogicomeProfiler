#ifndef LOGICOME_PROFILER_H
#define LOGICOME_PROFILER_H

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include <map>
#include <algorithm>
#include <stdlib.h>
#include <string.h>

using namespace std;

class Node{
 public:
  Node(int x, int y, int z, int w, double f){
    a = x; b = y; c = z; l = w; pv = f;
  }
  Node(){}
  int a;  int b;  int c;  int l;  double pv;
};

class LogicomeProfiler{
 public:
  LogicomeProfiler(){
    _input_data_file_name = "";
    _output_file_name = "";
    _number_of_samples = 0;
    _number_of_elements = 0;
    _criteria = 0;
    _alpha = 0.0;
  }
  void SetParameters(char* argv[]);
  void Run();

 private:
  void ReadData();
  void CountElements();
  void CalcLogFact();
  void MakeLogicOccurrenceTensor();
  void FirstScreeningForFWER();
  void SecondScreeningForFWER();
  void CalcPVForFDR();
  void CalcProbSecondCondition(int a,int b,int c,double logic,double threshold,double& prob2,double& prob3);
  int CountLogicOccurrence(int a, int b,int l);
  int CountLogicOccurrence(int a, int b, int c, int l);
  
  double CalcProb(int x, int y, int z, int N, double threshold);
  double CalcCombination(int x, int y);
  double max3(double x, double y, double z);

  vector<vector<vector<int> > > _logic_occurrence_tensor;
  vector<vector<bool> > _dataset;
  vector<int> _count_vector;
  vector<Node> _result_vector;
  vector<vector<Node> > _fdr_result_vector;
  vector<map<string,double> > _fdr_result_map;
  vector<double> _log_fact;
  vector<string> _name_vector;
  string _input_data_file_name;
  string _output_file_name;
  int _number_of_samples;
  int _number_of_elements;
  double _alpha;
  int _iteration_number;
  int _criteria;
};

#endif
