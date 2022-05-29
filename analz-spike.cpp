#include <vector>
#include <string>
#include <set>
#include <cmath>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <iterator>
#include <map>
#include <sstream>
using namespace std;

string replace_space(string s){
  for(int i=0;i<s.size();i++)
    if(s[i] == ' ') s[i] = '-';

  return s;
}


struct GISAID {
  string Gene_name;
  string Isolate_name;
  string Date;
  string Isolate_ID;
  string History;
  string Type;
  string Host;
  string Originating_lab;
  string Submitting_lab;
  string Submitter;
  string Location;

  string sequence;

  GISAID(vector<string> &elems, string & seq){

    Gene_name = elems[0];
    Isolate_name = elems[1];
    Date = elems[2];
    Isolate_ID = elems[3];
    History = elems[4];
    Type = elems[5];
    Host = elems[6];
    Originating_lab = elems[7];
    Submitting_lab = elems[8];
    Submitter = elems[9];
    Location = elems[10];
    sequence = seq;
  }

};




vector<string> split(const string &s, char delim) {
    vector<string> elems;
    stringstream ss(s);
    string number;
    while(getline(ss, number, delim)) {
      elems.push_back(number);
    }
    return elems;
}



int main(int argc, char ** argv){
  //SEQUENCES
  vector<string> desc;
  //  vector<string> seq;


//   ifstream in; in.open(argv[1]);  
//   istream_iterator<string> iter(in),end;
//   int nline = 0;
//   while(iter!=end){
    
//     //seq.push_back(*iter); ++iter;
//     if(nline%2==0){
//       cout<<nline<<" "<<*iter<<"\n";;
//     }

//     ++nline;
//     ++iter;
//   }

  vector<GISAID> gisaid;
  ifstream in; in.open(argv[1]);
  int min = atoi(argv[2]);
  int max = atoi(argv[3]);
  int nline = 0;
  vector<string> tok;
  string seq;
  while(!in.eof()){
    string line;
    getline (in,line);
   
    if(nline%2==0){
      tok = split(line,'|');
    } else {
      seq = line;
      
    }
    ++nline;
    if(nline%2==0){
      if(tok.size()==11){
	GISAID g(tok,seq);
	gisaid.push_back(g);
	//	cerr<<nline<<" "<<tok[1]<<" "<<seq<<"\r";
      }
    }

    if(nline%10000==0) 
      cerr<<"Lines ... "<<nline<<"\r";

    //    if(gisaid.size()==1000000){
    //      cerr<< "DEBUG END ..."<<endl;
    //      break;
    //    }
  }
  cerr<<">> Processed "<<gisaid.size()<<" sequences"<<endl;

  cerr<<">> Searching ..."<<endl;

  //  map<string,
  for(int i=0;i<gisaid.size();i++){
    if(gisaid[i].Host!="Human") continue;
    std::size_t xpos = gisaid[i].sequence.find("X");
    if(xpos==string::npos&&gisaid[i].sequence.size()>min&&gisaid[i].sequence.size()<max){
      //      set<string>::iterator seq_pos = uniq_seq.find(gisaid[i].sequence);
      cout<<replace_space(gisaid[i].Isolate_name)<<"\t"
	  <<gisaid[i].Host<<"\t"
	  <<gisaid[i].Date<<"\t"
	  <<gisaid[i].sequence<<"\n";
      //      uniq_seq.insert(gisaid[i].sequence);
      
    }

  }
  //  cerr<<">> Found "<<uniq_seq.size()<<" unique sequnces"<<endl;
  cerr <<"DONE ..."<<endl;

}

