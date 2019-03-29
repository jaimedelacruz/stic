#include <string>
#include <cstring>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>

#include "rh.h"
#include "readAtomFile.h"

using namespace std;

// ------------------------------------------------------------------------------ //

char **matrix_char2(int Nrow, int Ncol)
{
  char *theMatrix, **Matrix;
  int   typeSize = sizeof(char), pointerSize = sizeof(char *);

  theMatrix = (char *)  calloc(Nrow * Ncol, typeSize);
  Matrix    = (char **) malloc(Nrow * pointerSize);
  for (int i = 0;  i < Nrow;  i++, theMatrix += Ncol)
    Matrix[i] = theMatrix;

  return Matrix;
}

// ------------------------------------------------------------------------------ //

inline std::string removeSpaces(std::string input){
  input.erase(std::remove(input.begin(),input.end(),' '),input.end());
  return input;
}

// ------------------------------------------------------------------------------ //

inline std::vector<std::string> strsplit(std::string &var, std::string token, bool rmspaces){
    
  std::vector<std::string> res;
  std::stringstream ss(var);
    
  std::string tmp;
  if(rmspaces)while(std::getline(ss, tmp, *(token.c_str()))) res.push_back(removeSpaces(tmp));
  else while(std::getline(ss, tmp, *(token.c_str()))) res.push_back(tmp);
  
  return res;
}

// ------------------------------------------------------------------------------ //

inline bool file_exists_bool(const std::string &name) {
  std::ifstream f(name.c_str(),std::ifstream::in);
  if (f.good()) {
    f.close();
    return true;
  } else {
    f.close();
    return false;
  }
}

// ------------------------------------------------------------------------------ //

inline std::string cleanLine(std::string var,  std::string token, int rmspaces)
{
  std::vector<std::string> res = strsplit(var, token, rmspaces);
  return res[0];
}

// ------------------------------------------------------------------------------ //

char** readAtomFile(char* const fname_in, const char com, int remove_empty)
{
  static std::string const nu = "\0";
  std::string const fname = string(fname_in);  
  bool exists = file_exists_bool(fname);
  if(!exists) fprintf(stderr,"error: readAtomFile: cannot open file [%s]\n", fname.c_str());
  
  std::ifstream in(fname.c_str(), std::ios::in | std::ios::binary);
  std::vector<std::string> res;
  std::string line, com1(" ");
  com1.assign(1, com);
  
  if(in){
    while(std::getline(in, line)){
      if(remove_empty) if((line == "") || (line == " ")) continue;
      if(line[0] != com) res.push_back(cleanLine(line,com1,false)+nu);
    }
  }

  
  
  // --- ok now we have the textfile in a C++ vector of string
  // --- Let's copy it to a C-style 2D array

  const size_t iim = res.size();
  size_t linemax = 0;
  for(size_t ii=0; ii<iim; ++ii)
    linemax = std::max<size_t>(linemax, res[ii].size());

  char **tmp = (char**)matrix_char2(iim, linemax+1);

  for(size_t ii = 0; ii<iim; ++ii){
    strcpy(tmp[ii], res[ii].c_str());
  }
    
  return tmp;
}

// ------------------------------------------------------------------------------ //

