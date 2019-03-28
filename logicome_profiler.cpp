#include "logicome_profiler.h"

bool compare_node( const Node& left, const Node& right ) {
  return(left.pv < right.pv);
}

void LogicomeProfiler::Run(){
  ReadData();
  CountElements();
  MakeLogicOccurrenceTensor();
  CalcLogFact();
  FirstScreening();
  SecondScreening();
}

void LogicomeProfiler::SecondScreening(){
  int N = _number_of_samples;
  long long int ne = _number_of_elements;
  long long int max_size_count = 4*ne*(ne-1)*(ne-2); 
  double bf_threshold = _alpha/max_size_count;
  double threshold = _alpha/_result_vector.size();
  ofstream ofs(_output_file_name.c_str());
  ofs << "1st_significance_level" << bf_threshold << " 2nd_significance_level " << threshold << endl;
  ofs << "Logic_formula Pv_of_condition_a Pv_of_condition_b Pv_of_condition_c p(c|l(a,b)) p(l(a,b)|c) Logic_type"<< endl;
  
  for(vector<Node>::iterator it = _result_vector.begin(); it != _result_vector.end(); ++it){    
    int a  = it->a;
    int b  = it->b;
    int c  = it->c;
    int logic = it->l;
    
    int x = _logic_occurrence_tensor[a][b][logic];
    int y = _count_vector[c];
    int z = CountLogicOccurrence(a,b,c,logic);
    double prob2;
    double prob3;
    double es1 = (double)z/x;
    double es2 = (double)z/y;
    
    if(logic == 1){
      int w1 = _logic_occurrence_tensor[a][c][1];
      int v1 = _count_vector[a];
      prob2 = CalcProb(x,w1,z,v1,threshold);
      
      int w2 = _logic_occurrence_tensor[b][c][1];
      int v2 = _count_vector[b];
      prob3 = CalcProb(x,w2,z,v2,threshold);
      
    }else if(logic == 2){
      
      int r = _logic_occurrence_tensor[a][b][1];
      int s1 = _logic_occurrence_tensor[a][c][5];
      int t = r - CountLogicOccurrence(a,b,c,1);
      int u1 =  _count_vector[a];
      prob2 = CalcProb(r,s1,t,u1,threshold);
      
      int s2 = _logic_occurrence_tensor[b][c][5];
      int u2 =  _count_vector[b];
      prob3 = CalcProb(r,s2,t,u2,threshold);
      
    }else if(logic == 3){
      
      int r = _logic_occurrence_tensor[a][b][4];
      int s1 = _logic_occurrence_tensor[a][c][4];
      int t = r - CountLogicOccurrence(a,b,c,4);
      int u1 =  _count_vector[a];
      prob2 = CalcProb(r,s1,t,N-u1,threshold);
      
      int s2 = _logic_occurrence_tensor[b][c][4];
      int u2 =  _count_vector[b];
      prob3 = CalcProb(r,s2,t,N-u2,threshold);
      
    }else if(logic == 4){
      
      int w1 = _logic_occurrence_tensor[c][a][5];
      int v1 = _count_vector[a];
      prob2 = CalcProb(x,w1,z,N-v1,threshold);
      
      int w2 = _logic_occurrence_tensor[c][b][5];
      int v2 = _count_vector[b];
      prob3 = CalcProb(x,w2,z,N-v2,threshold);
    }else if(logic == 5){
      
      int w1 = _logic_occurrence_tensor[a][c][1];
      int v1 = _count_vector[a];
      prob2 = CalcProb(x,w1,z,v1,threshold);
      
      int w2 = _logic_occurrence_tensor[c][b][5];
      int v2 = _count_vector[b];
      prob3 = CalcProb(x,w2,z,N-v2,threshold);
      
    }else if(logic == 6){
      int r = _logic_occurrence_tensor[a][b][5];
      int s1 = _logic_occurrence_tensor[a][c][5];
      int t = r - CountLogicOccurrence(a,b,c,5);
      int u1 =  _count_vector[a];
      prob2 = CalcProb(r,s1,t,u1,threshold);
      
      int s2 = _logic_occurrence_tensor[b][c][4];
      int u2 =  _count_vector[b];
      prob3 = CalcProb(r,s2,t,N-u2,threshold);
    }

    if(prob2 < threshold && prob3 < threshold){
      ofs << _name_vector[c] << "=" << _name_vector[a] << "," << _name_vector[b] << " " << it->pv << " "<< prob2 << " "<< prob3 << " " << es1 << " " << es2 << " " << logic << endl;
    }
  }
}

void LogicomeProfiler::FirstScreening(){
  int N = _number_of_samples;
  
  long long int count = 0;
  long long int ne = _number_of_elements;
  long long int max_size_count = 4*ne*(ne-1)*(ne-2); 
  double bf_threshold = _alpha/max_size_count;
  _result_vector.reserve(10000);

  for(int c = 0; c < _number_of_elements; c++){
    for(int a = 0; a < _number_of_elements; a++){
      if(a!=c){
	for(int b = 0; b < _number_of_elements; b++){
	  if(b!=a && b!=c){
	    for(int l = 1; l <=6; l++){
	      if(l>=5 || a < b){
		count++;
		int x = _logic_occurrence_tensor[a][b][l];
		int y = _count_vector[c];
		int z = CountLogicOccurrence(a,b,c,l);
		double prob = CalcProb(x,y,z,N,bf_threshold);
		
		if(prob < bf_threshold){
		  Node new_node(a,b,c,l,prob);
		  _result_vector.push_back(new_node);
		}

		if(count%1000000 == 0){
		  cout << count<< "/" << max_size_count << " " << _result_vector.size() << endl;
		}
	      }
	    }
	  }
	}
      }
    }
  }
}

void LogicomeProfiler::CalcLogFact(){
  int N = _number_of_samples;
  _log_fact.resize(N+1,0.0);

  double temp = 0.0;
  for(int i = 1; i<=N;i++){
    temp += log(i);
    _log_fact[i] = temp;
  }
}

double LogicomeProfiler::CalcCombination(int x, int y){
  return(_log_fact[x] - (_log_fact[y] + _log_fact[x-y]));
}

double LogicomeProfiler::CalcProb(int x, int y, int z, int N, double threshold){
  int z_max = x > y ? y : x;
  double sum = 0.0;
  double temp = CalcCombination(N,y);
  for(int k = z; k<=z_max; k++){
    double prob = -temp;
    prob = prob + CalcCombination(x,k)+CalcCombination(N-x,y-k);
    prob = exp(prob);
    sum += prob;
    if(sum > threshold){
      return(sum);
    }
  }
  return(sum);
}

int LogicomeProfiler::CountLogicOccurrence(int a, int b, int l){
  int x = 0;
  for(int i = 0; i < _number_of_samples;i++){
    switch(l){
    case 1:
      if(_dataset[i][a]==true && _dataset[i][b]==true){
	x += 1;
      }
      break;
    case 2:
      if(!(_dataset[i][a]==true && _dataset[i][b]==true)){
	x += 1;
      }
      break;
    case 3:
      if(_dataset[i][a]==true || _dataset[i][b]==true){
	x += 1;
      }
      break;
    case 4:
      if(!(_dataset[i][a]==true || _dataset[i][b]==true)){
	x += 1;
      }
      break;
    case 5:
      if(_dataset[i][a]==true && _dataset[i][b]!=true){
	x += 1;
      }
      break;
    case 6:
      if(_dataset[i][a] !=true || _dataset[i][b]==true){
	x += 1;
      }
      break;
    }
  }
  return x;
}

void LogicomeProfiler::MakeLogicOccurrenceTensor(){
  
  
  for(int i = 0; i < _number_of_elements; i++){
    vector<vector<int> > temp_matrix; temp_matrix.reserve(_number_of_elements);
    for(int j = 0; j < _number_of_elements;j++){
      vector<int> temp_vector; temp_vector.reserve(7);
      temp_vector.push_back(0);
      
      for(int k = 1;k <= 6;k++){
	int count = CountLogicOccurrence(i,j,k);
	temp_vector.push_back(count);
      }
      temp_matrix.push_back(temp_vector);
    }
    _logic_occurrence_tensor.push_back(temp_matrix);
  }
}

int LogicomeProfiler::CountLogicOccurrence(int a, int b, int c, int l){
  int z  = 0;
  for(int i = 0; i < _number_of_samples;i++){
    switch(l){
    case 1:
      if(_dataset[i][a]==true && _dataset[i][b]==true && _dataset[i][c]==true){
	z += 1;
      }
      break;
    case 2:
      if(!(_dataset[i][a]==true && _dataset[i][b]==true) && _dataset[i][c]==true){
	z += 1;
      }
      break;
    case 3:
      if((_dataset[i][a]==true || _dataset[i][b]==true) && _dataset[i][c]==true){
	z += 1;
      }
      break;
    case 4:
      if((!(_dataset[i][a]==true || _dataset[i][b]==true)) && _dataset[i][c]==true){
	z += 1;
      }
      break;
    case 5:
      if(_dataset[i][a]==true && _dataset[i][b]!=true && _dataset[i][c]==true){
	z += 1;
      }
      break;
    case 6:
      if((_dataset[i][a]!=true || _dataset[i][b]==true) && _dataset[i][c]==true){
	z += 1;
      }
      break;
    }
  }
  return z;
}

void LogicomeProfiler::CountElements(){
  _count_vector.reserve(_number_of_elements);
  for(int i = 0; i < _number_of_elements; i++){
    int count = 0;
    for(int j = 0; j < _number_of_samples;j++){
      if(_dataset[j][i]==true){
	count += 1;
      }
    }
    _count_vector.push_back(count);
  }
}

void LogicomeProfiler::ReadData(){
  ifstream fp;
  fp.open(_input_data_file_name.c_str(), ios::in);
  if (!fp) {
    cout << "Cannot open " + _input_data_file_name << endl;
    exit(1);
  }
  char buf[1000000];
  fp.getline(buf, 1000000);
  char* tp;
  tp = strtok(buf, " ");
  _number_of_samples = atoi(tp);
  tp = strtok(NULL, " ");
  _number_of_elements = atoi(tp);
  _name_vector.reserve(_number_of_elements);
  fp.getline(buf, 1000000);
  
  _dataset.resize(_number_of_samples, vector<bool>(_number_of_elements, 0));
  for (int i = 0; i < _number_of_elements; i++){
    fp.getline(buf, 1000000);
    tp = strtok(buf, " ");
    _name_vector.push_back(string(tp));
    for(int j = 0; j < _number_of_samples; j++){
      tp = strtok(NULL, " ");
      _dataset[j][i] = atoi(tp) == 1 ? true : false;
    }
  }
  fp.close();
}

void LogicomeProfiler::SetParameters(char* argv[]){
  _input_data_file_name = argv[1];
  _output_file_name = argv[2];
}
