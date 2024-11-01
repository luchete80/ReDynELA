#include "lsdynaReader.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

bool isComment(const string &line){
  
  bool ret = false;
  if (line.find('$') == 0)
    ret = true;
  return ret;
}
// str.erase(
  // std::remove_if(str.begin(), str.end(), [](char chr){ return chr == '&' || chr == ' ';}),
  // str.end());

int findNextCommandLine(const int &curr_line, std::vector <std::string> line) {
  string first;
  int i = curr_line;
  bool end = false;
  while (!end) {
    i++;
    if (i==line.size()) end = true; 
      first = line[i].substr(0,1);
      //cout << "first "<<first <<", "<<line[i]<<endl;
      //Can a line start with spaces??
      //if (first.find("$")!=string::npos)
      if (first=="*"){
        //cout << "FOUND * at line "<<i<<endl;
        end = true;
      }
  }// while end
  
        //cout << "* NOT FOUND "<<i<<endl;
  return i;
}

void lsdynaReader::removeComments(){
  std::vector<std::string>::iterator it;
  int i=0;
  for (it = m_line.begin();it!=m_line.end();){
    //cout << *it<<endl;
    string first = it->substr(0,1);
    //cout << first <<endl;
    if (first=="$"){
      m_line.erase(it);
      m_line_count--;
    } else {
      it ++;
      //cout << "NOT $ FOUND "<<first<<", LIN: "<<*it<<endl;
    }
    i++;
  }
}

void lsdynaReader::readNodes() {
  bool end = false;
  int i = 0;
  while (!end){

      if (m_line[i].find("*NODE") != std::string::npos){
        cout << "Node command found at line "<<m_line_count<<endl;
        
        m_node_count = findNextCommandLine(i,m_line) - i;
        cout << "Node count: "<<m_node_count<<endl;
        end = true;
      }
    if (i==m_line_count) {
      end = true;
      cout << "NODES not defined "<<endl;
    }
    i++;
  } 
  
  for (i=0;i<m_node_count;i++){
    int id;
    ls_node nod();
  }

}  //line


lsdynaReader::lsdynaReader(const char *fname){
  string line;
  m_line_count = 0;
  int start, end;
  char dl = ' ';
  ifstream file(fname);
  if (file.is_open()) {
    while (getline(file, line)) {
      m_line.push_back(line);
      m_line_count++;
    }
  }
        file.close();
  cout << "Line count: "<< m_line_count << endl;
  
  removeComments();
  cout << "Line count w/o comments: "<< m_line_count << endl;
  
  readNodes();
  
  std::vector<std::string>::iterator it = m_line.begin();
  for (int i=0;i<m_line_count;i++) {
          start = end = 0;
          //cout << m_line[i]<<endl;
          // start = line.find_first_not_of(dl, end);
          if (!isComment(line)){
            if (m_line[i].find('*') != std::string::npos){
              //cout << "COMMAND FOUND at "<< line.find("*")<<", line "<<m_line_count<<endl;
                
            } 
            if (m_line[i].find("*NODE") != std::string::npos){
              cout << "Node command found at line "<<m_line_count<<endl;
              
              cout << "Node count: "<<findNextCommandLine(i,m_line) - i<<endl;
            }
          } else {
            // cout << "Commented line "<<m_line_count<<endl;
          }
          // if (start != string::npos){
            
          // }
          
          
    it++;

    }  //line

}