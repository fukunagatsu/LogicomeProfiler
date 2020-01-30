#include "logicome_profiler.h"

bool compare_node( const Node& left, const Node& right ) {
  return(left.pv > right.pv);
}

void LogicomeProfiler::Run(){
  ReadData();
  CountElements();
  MakeLogicOccurrenceTensor();
  CalcLogFact();
  if(_criteria == 0){
    FirstScreeningForFWER();
    SecondScreeningForFWER();
  }else if(_criteria == 1){
    CalcPVForFDR();
  }
}

void LogicomeProfiler::CalcPVForFDR(){
  int N = _number_of_samples;
  long long int count = 0;
  long long int ne = _number_of_elements;
  long long int max_size_count = 4*ne*(ne-1)*(ne-2);
  long long int ignored_count[3] = {0};
  double threshold = _alpha;
  _fdr_result_vector.reserve(3);
  
  for(int i = 0; i < 3; i++){
    vector<Node> temp_vector; temp_vector.reserve(10000);
    _fdr_result_vector.push_back(temp_vector);
  }
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
		double prob[3] = {0.0}; 
		prob[0] = CalcProb(x,y,z,N,threshold);
		CalcProbSecondCondition(a,b,c,l,threshold,prob[1],prob[2]);

		for(int i = 0; i < 3; i++){
		  double temp = (double)max_size_count*log(max_size_count)/(max_size_count - ignored_count[i]);		 
		  if(prob[i]*temp < threshold){
		    Node new_node(a,b,c,l,prob[i]);
		    _fdr_result_vector[i].push_back(new_node);
		    
		  }else{
		    ignored_count[i]++;
		  }
		}

		if(count%1000000 == 0){
		  cout << count<< "/" << max_size_count << " " << _fdr_result_vector[0].size() << " " << _fdr_result_vector[1].size()
		       << " " << _fdr_result_vector[2].size()<< endl;
		}
	      }
	    }
	  }
	}
      }
    }
  }
  
  for(int i = 0 ; i < 3; i++){
    sort(_fdr_result_vector[i].begin(), _fdr_result_vector[i].end(), compare_node);
    
    long long int denominator = _fdr_result_vector[i].size();
    double prev_qv = 100;
    for(vector<Node>::iterator it = _fdr_result_vector[i].begin(); it != _fdr_result_vector[i].end(); ++it){
      it->pv = min(prev_qv, it->pv*((double)max_size_count*log(max_size_count)/denominator));
      prev_qv = it->pv;
      denominator--;
    }

    bool break_flag = false;
    for(vector<Node>::iterator it = _fdr_result_vector[i].begin(); it != _fdr_result_vector[i].end(); ++it){
      if(it->pv < threshold){
	_fdr_result_vector[i].erase(_fdr_result_vector[i].begin(),it);
	_fdr_result_vector[i].shrink_to_fit();
	break_flag = true;
	break;
      }
    }
    if(!break_flag){
      _fdr_result_vector[i].clear();
      _fdr_result_vector[i].shrink_to_fit();
    }
    cout << _fdr_result_vector[i].size() << endl;
  }
  _fdr_result_map.resize(2);
  for(int i = 1; i < 3; i++){
    for(vector<Node>::iterator it = _fdr_result_vector[i].begin(); it != _fdr_result_vector[i].end(); ++it){
      string temp_string = to_string(it->a)+" "+to_string(it->b)+" "+to_string(it->c)+" "+to_string(it->l);
      _fdr_result_map[i-1].insert(make_pair(temp_string,it->pv));
    }
    _fdr_result_vector[i].clear();
    _fdr_result_vector[i].shrink_to_fit();
  }
  

  ofstream ofs(_output_file_name.c_str());
  ofs << "significance_level " << threshold << endl;
  ofs << "Logic_formula Qv_of_condition_a Qv_of_condition_b Qv_of_condition_c Logic_type"<< endl;
  for(vector<Node>::iterator it = _fdr_result_vector[0].begin(); it != _fdr_result_vector[0].end(); ++it){
    string temp_string = to_string(it->a)+" "+to_string(it->b)+" "+to_string(it->c)+" "+to_string(it->l);
    if(_fdr_result_map[0].find(temp_string) != _fdr_result_map[0].end()
       && _fdr_result_map[1].find(temp_string) != _fdr_result_map[1].end()){
      
      ofs << _name_vector[it->c] << "=" << _name_vector[it->a] << "," << _name_vector[it->b] << " " << it->pv << " "<< _fdr_result_map[0][temp_string] << " "<< _fdr_result_map[1][temp_string] << " " << it->l << endl;
    }
  }
  ofs.close();
}

void LogicomeProfiler::CalcProbSecondCondition(int a,int b,int c,double logic,double threshold,double& prob2,double& prob3){
  int x = _logic_occurrence_tensor[a][b][logic];
  int y = _count_vector[c];
  int z = CountLogicOccurrence(a,b,c,logic);
   int N = _number_of_samples;
   
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
}

void LogicomeProfiler::SecondScreeningForFWER(){
  int N = _number_of_samples;
  long long int ne = _number_of_elements;
  long long int max_size_count = 4*ne*(ne-1)*(ne-2); 
  double threshold = _alpha/max_size_count;
  ofstream ofs(_output_file_name.c_str());
  ofs << "significance_level " << threshold << endl;
  ofs << "Logic_formula Pv_of_condition_a Pv_of_condition_b Pv_of_condition_c Logic_type"<< endl;
  
  for(vector<Node>::iterator it = _result_vector.begin(); it != _result_vector.end(); ++it){    
    int a  = it->a;
    int b  = it->b;
    int c  = it->c;
    int logic = it->l;
    double prob2 = 0.0;
    double prob3 = 0.0;
    CalcProbSecondCondition(a,b,c,logic,threshold,prob2,prob3);

    if(prob2 < threshold && prob3 < threshold){
      ofs << _name_vector[c] << "=" << _name_vector[a] << "," << _name_vector[b] << " " << it->pv << " "<< prob2 << " "<< prob3 << " " << logic << endl;
    }
  }
}

void LogicomeProfiler::FirstScreeningForFWER(){
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
  _criteria = atoi(argv[3]);
  _alpha = atof(argv[4]);
}
