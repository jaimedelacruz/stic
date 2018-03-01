//#include <sys/types.h>
#ifndef MMEM_H
#define MMEM_H

#include <stdint.h>
#include <algorithm>

namespace mmem{
/*
  Memory allocation functions -> templates overview
 */ 
template <class T,class U> T **mat2d(U nx1, U nx2);
template <class T,class U> T ***mat3d(U nx1, U nx2, U nx3);
template <class T,class U> T ****mat4d(U nx1, U nx2, U nx3, U nx4);
template <class T,class U> T **var2dim(T *data, U nx1, U nx2);
template <class T,class U> T ***var3dim(T *data, U nx1, U nx2, U nx3);
template <class T,class U> T ****var4dim(T *data, U nx1, U nx2, U nx3, U nx4);
//
template <class T,class U> void del_mat2d(T **p);
template <class T,class U> void del_mat3d(T ***p);
template <class T,class U> void del_mat4d(T ****p);
//

//
/*
  Definitions
 */
template <class T,class U> T **mat2d(U nx1, U nx2){
  T **p;
  p = new T* [nx1];
  p[0] = new T [nx1 * nx2];
  for(U x1=1;x1<=nx1-1;++x1) p[x1] = p[x1-1] + nx2;
  return p;
}
//
template <class T,class U> T ***mat3d(U nx1, U nx2, U nx3){
  //
  T ***p;
  p=new T** [nx1];
  p[0]=new T* [nx1*nx2];
  p[0][0]=new T [nx1*nx2*nx3];
  for(U x2=1;x2<=(nx2-1);++x2) p[0][x2]=p[0][x2-1]+nx3;
  for(U x1=1;x1<=(nx1-1);++x1){
    p[x1]=p[x1-1]+nx2;
    p[x1][0]=p[x1-1][0]+nx2*nx3;
    for(U x2=1;x2<=(nx2-1);++x2) p[x1][x2]=p[x1][x2-1]+nx3;
  }
  return p;
}
//
template <class T,class U> T ****mat4d(U nx1, U nx2, U nx3, U nx4){
  T ****p;
  p = new T*** [nx1];
  p[0] = new T** [nx1 * nx2];
  p[0][0] = new T* [nx1 * nx2 * nx3];
  p[0][0][0] = new T [nx1 * nx2 * nx3 * nx4];
  //
  for(U x3=1;x3<=nx3-1;++x3) p[0][0][x3]=p[0][0][x3-1] + nx4;
  for(U x2=1;x2<=nx2-1;++x2){
    p[0][x2] = p[0][x2-1] + nx3;
    p[0][x2][0] = p[0][x2-1][0] + (nx3 * nx4);
    for(U x3=1;x3<=nx3-1;++x3) p[0][x2][x3] = p[0][x2][x3-1] + nx4;
  }
  for(U x1=1;x1<=nx1-1;++x1) {
    p[x1] = p[x1-1] + nx2;
    p[x1][0]=p[x1-1][0]+nx2*nx3;
    p[x1][0][0]=p[x1-1][0][0] + (nx2 * nx3 * nx4);
    for(U x3=1;x3<=nx3-1;++x3) p[x1][0][x3]=p[x1][0][x3-1] + nx4;
    for(U x2=1;x2<=nx2-1;++x2){
      p[x1][x2]=p[x1][x2-1] + nx3;
      p[x1][x2][0]=p[x1][x2-1][0] + (nx3*nx4);
      for(U x3=1;x3<=nx3-1;++x3) p[x1][x2][x3]=p[x1][x2][x3-1] + nx4;
    }
  }
  return p;
}
//
template <class T,class U> T **var2dim(T *data, U nx1, U nx2){
  T **p;
  p=new T* [nx1];
  p[0] = data;
  for(U x1=1;x1<=nx1-1;++x1) p[x1]=p[x1-1] + nx2;
  return p;
}
//
template <class T,class U> T ***var3dim(T *data, U nx1, U nx2, U nx3){
  T ***p;
  p=new T** [nx1];
  p[0]=new T* [nx1 * nx2];
  p[0][0]=data;
  for(U x2=1;x2<=nx2-1;++x2) p[0][x2]=p[0][x2-1] + nx3;
  for(U x1=1;x1<=nx1-1;++x1){
    p[x1]=p[x1-1] + nx2;
    p[x1][0]=p[x1-1][0] + nx2 * nx3;
    for(U x2=1;x2<=nx2-1;++x2) p[x1][x2]=p[x1][x2-1] + nx3;
  }
  return p;
}
//
template <class T,class U> T ****var4dim(T *data, U nx1, U nx2, U nx3, U nx4){
  T ****p;
  p=new T*** [nx1];
  p[0]=new T** [nx1*nx2];
  p[0][0]=new T* [nx1*nx2*nx3];
  p[0][0][0]=data;
  for(U x3=1;x3<=nx3-1;++x3) p[0][0][x3]=p[0][0][x3-1] + nx4;
  for(U x2=1;x2<=nx2-1;++x2){
    p[0][x2]=p[0][x2-1] + nx3;
    p[0][x2][0]=p[0][x2-1][0] + (nx3 * nx4);
    for(U x3=1;x3<=nx3-1;++x3) p[0][x2][x3]=p[0][x2][x3-1] + nx4;
  }
  for(U x1=1;x1<=nx1-1;++x1) {
    p[x1]=p[x1-1] + nx2;
    p[x1][0]=p[x1-1][0] + (nx2 * nx3);
    p[x1][0][0]=p[x1-1][0][0] + (nx2 * nx3 * nx4);
    for(U x3=1;x3<=nx3-1;++x3) p[x1][0][x3]=p[x1][0][x3-1] + nx4;
    for(U x2=1;x2<=nx2-1;++x2){
      p[x1][x2]=p[x1][x2-1] + nx3;
      p[x1][x2][0]=p[x1][x2-1][0] + (nx3 * nx4);
      for(U x3=1;x3<=nx3-1;++x3) p[x1][x2][x3]=p[x1][x2][x3-1] + nx4;
    }
  }
  return p;
}
//
template <class T,class U> void del_mat2d(T **p){
  delete[] (p[0]);
  delete[] (p);
}
//
template <class T,class U> void del_mat3d(T ***p){
  delete[] (p[0][0]);
  delete[] (p[0]);
  delete[] (p);
}
//
template <class T,class U> void del_mat4d(T ****p){
  delete[] (p[0][0][0]);
  delete[] (p[0][0]);
  delete[] (p[0]);
  delete[] (p);
}
};
#endif
