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

string ref_spike = 
  "MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFS"
  "NVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIV"
  "NNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLE"
  "GKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQT"
  "LLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETK"
  "CTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISN"
  "CVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIAD"
  "YNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPC"
  "NGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVN"
  "FNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITP"
  "GTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSY"
  "ECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTI"
  "SVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQE"
  "VFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDC"
  "LGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAM"
  "QMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALN"
  "TLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRA"
  "SANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPA"
  "ICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDP"
  "LQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDL"
  "QELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDD"
  "SEPVLKGVKLHYT";

double entropy(map<char,int> &ss){
  double N = 0;
  for(auto &x : ss) N+=x.second;

  double E = 0;
  for(auto &x : ss){
    double p = x.second/N;
    E += p*log(p);
  }

  return E;
}


double nmut (map<char,int> &ss){
  double N = 0;
  for(auto &x : ss) N+=x.second;
  double nmax= 0;
  for(auto &x : ss){
    if(x.second>nmax) nmax = x.second;
  }
  
  return (N-nmax)/N;

}


double nmut (map<char,int> &ss, char ref){
  double N = 0;
  for(auto &x : ss) N+=x.second;
  double nmut= 0;
  for(auto &x : ss){
    if(x.first!=ref) nmut += x.second ;
  }
  
  return (nmut)/N;

}

char cons (map<char,int> &ss){
  char c;
  int nmax = 0;
  for(auto &x : ss){
    if(x.second>nmax){
      nmax = x.second;
      c = x.first;
    }
  }
  
  return c;

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


struct SPIKE {
  string ID;
  string Organism;
  string date;
  string sequence;
  
};


struct BLAST {
  string query;
  string subj;
};

int main(int argc, char ** argv){
  cerr<<ref_spike<<endl;
  cerr<<ref_spike.size()<<endl;
  vector<SPIKE> spike;
  map<string,BLAST> blast;

  {
    ifstream in; in.open(argv[1]);  
    istream_iterator<string> iter(in),end;
    int nline = 0;

    while(iter!=end){
      SPIKE s;
      s.ID = *iter; ++iter;
      s.Organism = *iter; ++iter;
      s.date = *iter; ++iter;
      s.sequence = *iter; ++iter;
      int pos = s.sequence.find("*");
      if(pos != string::npos){
	s.sequence = string(s.sequence.begin(),s.sequence.begin()+pos);
      }
      spike.push_back(s);
      ++nline;
      if(nline%10000==0)
	cerr<<">> line "<<nline<<"\r";
      //DEBUG
      //      if(nline == 1000000)
      //	break;
    }

  cerr<<">> READ "<<spike.size()<<" sequences"<<endl;
  }

  
  {
    ifstream in; in.open(argv[2]);  
    istream_iterator<string> iter(in),end;
    int nline = 0;
    while(iter!=end){
      string s = *iter; ++iter;
      string q = *iter; ++iter;
      string h = *iter; ++iter;
      //SUBJECT IS REF-SPIKE
      BLAST b; b.query = q; b.subj = h;

      blast[s] = b;
    }

    cerr <<"READ "<<blast.size() << " blasts"<<endl;
  }

  int range_beg = atoi(argv[3]);
  int range_end = atoi(argv[4]);

  

  //NOW do stats ...
  vector<map<char,int> > mutations;
  mutations.resize(ref_spike.size());
  for(int i = 0; i< spike.size();i++){
    //find
    map<string,BLAST>::iterator pos = blast.find(spike[i].sequence);

    if(pos == blast.end()){
      //cout<<"<"<<spike[i].sequence<<">"<<endl;
      //      cerr<<"> skipping "<<i<<endl;
      continue;
    }

    int na = 0;
    int nins = 0;
    for(int j=0;j<pos->second.subj.size();j++){
      if(pos->second.subj[j] != '-'){
	na = na + 1;
      }
      else {
	nins = nins + 1;
      }
    }
    if(na != ref_spike.size()){
      cerr<<i+1<<" NA "<<na<<" "<<nins<<" ("<<ref_spike.size()<<")"<<endl;
      cerr<<pos->second.query<<"\n";
      cerr<<pos->second.subj<<"\n\n";
      continue;
    }


    string t_date;
    for(int a=0;a<spike[i].date.size();a++){
      if(spike[i].date[a] != '-')
	t_date += spike[i].date[a];
    }
    
    int n_date = atoi(t_date.c_str());

    if(n_date<range_beg) continue;
    if(n_date>range_end) continue;
    
    //cout<<"DATE "<<n_date<<" "<<spike[i].date<<endl;

    //cout<<t_date<<" "<<atoi(t_date.c_str())<<endl;



    //    if(nins> 0){
    //      cout<<"NA "<<na<<" "<<nins<<endl;
    //    }

    na = 0;
    for(int j=0;j<pos->second.subj.size();j++){
      mutations[na][pos->second.query[j]] += 1;
      if(pos->second.subj[j] != '-'){
	na = na + 1;
      }
    }

    if(na!= ref_spike.size()){
      cerr<<"NA "<<na<<" error...."<<endl;
      continue;
    }


  }

  {
    //Proces ....
    cout<<"POSITION\tENTROPY\tNMUT\tNMUTR\tCONS\tREF\tDIST\n";
    for(int i=0;i<mutations.size();i++){

      map<char,int> & m = mutations[i];
      double ent = entropy(m);
      double nm  = nmut(m);
      double nmr = nmut(m,ref_spike[i]);
      char c = cons(m);

      cout<<i+1<<"\t"<<ent<<"\t"<<nm<<"\t"<<nmr<<"\t"<<c<<"\t"<<ref_spike[i]<<"\t";

     
      map<char,int>::iterator iter = m.begin();
      while(iter!=m.end()){
	cout<<iter->first<<" "<<iter->second<<" ";
	++iter;
      }
      cout<<"\n";
      

    }
  }


}
