#ifndef _LS_DYNA_READER_
#define _LS_DYNA_READER_

#include <string>
#include <vector>

class Keyword {
public:


};

struct ls_node {
  void ls_nodeconst (const int &id_, const double &x, const double &y, const double &z){
    m_id = id_;
    m_x[0]=x;m_x[1]=y;m_x[2]=z;
  }
  int m_id;
  double m_x[3];
};

struct ls_element {
  int id;
  std::vector <int> node;
};

class lsdynaReader{
public:  
  lsdynaReader(){}
  lsdynaReader(const char *);
protected:
  int m_line_count;
  int m_node_count;
  std::vector <std::string> m_line;
  void readNodes();
  void removeComments();
  
  std::vector < ls_node    > m_node;
  std::vector < ls_element > m_elem;
};


#endif
