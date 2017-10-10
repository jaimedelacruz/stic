#include <iostream>
#include <string>
#include <algorithm>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <fstream>      // std::ifstream
#include <sstream>
#include "cmemt.h"
#include "input.h"
#include "physical_consts.h"
//
using namespace std;
//
std::string removeSpaces(std::string input){
  input.erase(std::remove(input.begin(),input.end(),' '),input.end());
  return input;
}
//
std::vector<std::string> strsplit(std::string &var, std::string token, bool rmspaces){
  
  std::vector<std::string> res;
  std::stringstream ss(var);

  std::string tmp;
  if(rmspaces) while(std::getline(ss, tmp, *(token.c_str()))) res.push_back(removeSpaces(tmp));
  else while(std::getline(ss, tmp, *(token.c_str()))) res.push_back(tmp);
  
  return res;
}
//

std::vector<double> fill_lambdas(iput_t &input, bool air){

  /* --- Copy wavelength array and check which lines are computed at each wav --- */
  
  int nlambda = 0;
  for(auto &it: input.regions){
    
    it.off = nlambda;
    nlambda += it.nw;	
    it.wav.resize(it.nw);
    it.nu.resize(it.nw);
    for(int ii = 0; ii<it.nw; ii++){
      // Compute wavelength array in vacuum
      it.wav[ii] = inv_convl( convl(it.w0) + it.dw * double(ii)); // lambda
      it.nu[ii] = phyc::CC  / (it.wav[ii] * 1.e-8);  // nu (Freq. in s^-1)
    }
    
  }
  
  /* --- Create array of all lambdas --- */
  
  std::vector<double> lambda;
  lambda.resize(nlambda);

  int kk = 0;
  for(auto &it: input.regions){
    for(int ii = 0; ii<it.nw; ii++){
      if(air)
	lambda[kk++] = it.wav[ii];
      else
	lambda[kk++] = convl(it.wav[ii]);
    }
  }

  return(lambda);
}
//
iput_t read_input(std::string filename, bool verbose){
  
  iput_t input = {};
  
  // Default values
  input.mu = 1.0; // default
  input.npack = 1; // default
  input.mode = 1; // default
  input.nInv = 1; // default
  input.chi2_thres = 0.0; // default
  input.max_inv_iter = 40; // default
  input.master_threads = 1; // default
  input.sparse_threshold = 0.70; //
  input.wavelet_order = 4;
  input.wavelet_type = "daub"; 
  input.dpar = 1.e-2; // Default
  input.nw_tot = 0;
  input.nodes.nnodes = 0;
  input.nodes.depth_t = 0;
  input.verbose = false; // default
  memset(&input.nodes.toinv[0],0, 7*sizeof(int));
  input.solver = 0;
  input.centder = 0;
  input.init_step = -1.0;
  input.thydro = 1;
  input.nodes.nnodes = 0;
  input.dint = 0;
  input.keep_nne = 0;
  input.marquardt_damping = -1.0;
  input.svd_thres = 1.0e-5;
  input.svd_split = 1;
  input.regularize = 0.0;
  input.random_first = 0;
  input.depth_model = 0;
  
  memset(input.nodes.regul_type, 0, 6*sizeof(int));
  
  // Open File and read
  std::ifstream in(filename, std::ios::in | std::ios::binary);
  if (in){
    std::string iline;
    while (std::getline(in, iline)) {
      // Split string at "=" and remove empty spaces
      std::string::size_type n;
      n = iline.find("=");
      std::string key = removeSpaces(iline.substr(0, n));
      std::string field = removeSpaces(iline.substr(n+1));
      bool set = false;
      
      // Fill in input structure
      if     (key == "input_model"){
	input.imodel = field;
	set = true;
      }
      else if(key == "input_profiles"){
	input.iprof = field;
	set = true;	
      }
      else if(key == "output_model"){
	input.omodel = field;
	set = true;
      }
      else if(key == "output_atmos"){
	input.oatmos = field;
	set = true;
      }
      else if(key == "output_profiles"){
	input.oprof = field;
	set = true;
      }
      else if(key == "mu"){
	input.mu = atof(field.c_str());
	set = true;
      }
      else if(key == "marquardt_damping"){
	input.marquardt_damping = atof(field.c_str());
	set = true;
      }
      else if(key == "svd_thres"){
	input.svd_thres = atof(field.c_str());
	set = true;
      }
      else if(key == "mpi_pack"){
	input.npack = atoi(field.c_str());
	set = true;
      }
      else if(key == "depth_model"){
	input.depth_model = atoi(field.c_str());
	set = true;
      }
      else if(key == "rt_solver"){
	input.solver = atoi(field.c_str());
	set = true;
      }
      else if(key == "centered_derivatives"){
	input.centder = atoi(field.c_str());
	set = true;
      }
      else if(key == "recompute_hydro"){
	input.thydro = atoi(field.c_str());
	set = true;
      }
      else if(key == "depth_interpolation"){
	input.dint = atoi(field.c_str());
	set = true;
      }
      else if(key == "regularize"){
	input.regularize = atof(field.c_str());
	set = true;
      }
      else if(key == "depth_t"){
	input.nodes.depth_t = atoi(field.c_str());
	set = true;
      }
      else if(key == "regularization_type"){
	std::vector<std::string> param = strsplit(field,",");
	int dum = (int)param.size();
	if(dum > 7) dum = 7;
	for(int ii = 0; ii<dum; ii++) input.nodes.regul_type[ii] = std::stoi(param[ii]);
	set = true;
      }
      else if(key == "mode"){
	input.mode = atoi(field.c_str());
	set = true;
      }
      else if(key == "keep_nne"){
	input.keep_nne = atoi(field.c_str());
	set = true;
      }
      else if(key == "randomize_first"){
	input.random_first = atoi(field.c_str());
	set = true;
      }
      else if(key == "randomize_inversions"){
	input.nInv = atoi(field.c_str());
	set = true;
      }
      else if(key == "chi2_threshold"){
	input.chi2_thres = atof(field.c_str());
	set = true;
      }
      else if(key == "init_step"){
	input.init_step = atof(field.c_str());
	set = true;
      }
      else if(key == "atmosphere_type"){
	input.atmos_type = field;
	set = true;
	input.atmos_len = field.size();
      }
      else if(key == "instrument"){
	input.instrument = field;
	set = true;
	input.inst_len = field.size();
      }
      else if(key == "max_inv_iter"){
	input.max_inv_iter = atoi(field.c_str());
	set = true;
      }
      else if(key == "master_threads"){
	input.master_threads = atoi(field.c_str());
	set = true;
      }
      else if(key == "svd_split_singular"){
	input.svd_split = atoi(field.c_str());
	set = true;
      }
      else if(key == "wavelet_order"){
	input.wavelet_order = atoi(field.c_str());
	set = true;
      }
      else if(key == "verbose"){
	if(field == "true") input.verbose = true;
	else input.verbose = false;
	set = true;
      }
      else if(key == "sparse_threshold"){
	input.sparse_threshold = atof(field.c_str());
	set = true;
      }
      else if(key == "parameter_perturbation"){
	input.dpar = atof(field.c_str());
	set = true;
      }
      else if(key == "wavelet_type"){
	input.wavelet_type =field;
	set = true;
      }
      else if(key == "abundance_file"){
	input.abfile = field;
	set = true;
	input.ab_len = field.size();
      }
      else if(key == "lines"){
	std::vector<std::string> param = strsplit(field,",");
	for(auto &it: param) input.ilines.push_back(it); 
	set = true;
      }
      else if(key == "region"){
	std::vector<std::string> param = strsplit(field,",");
	region_t reg = {};
	std::string::size_type sz = 0;    
	reg.w0 = inv_convl(std::stod(param[0])); // convert to vacuum;
	reg.dw = std::stod(param[1]);
	reg.nw = std::stoi(param[2]);
	reg.cscal = std::stod(param[3]);
	//
	if(param.size() == 6){
	  reg.inst = param[4];
	  reg.ifile = param[5];
	}else{
	  reg.inst = "none";
	  reg.ifile = "none";
	}

	//
	input.nw_tot += reg.nw;
	input.regions.push_back(reg);
	set = true;
      }
      else if(key == "nodes_temp"){
	std::vector<std::string> param = strsplit(field,",");
	if(param.size() == 1){
	  double dum = std::stoi(param[0]);
	  int nn = 1;
	  if(int(dum+0.5) == 0) nn = 0;
	  input.nodes.temp.resize(nn, dum);
	} else if(param.size() > 1){
	  input.nodes.temp.resize(param.size());
	  for(int ii = 0; ii<(int)param.size(); ii++) input.nodes.temp[ii] = std::stod(param[ii]);
	}
	set = true;
      }
      else if(key == "nodes_vlos"){
	std::vector<std::string> param = strsplit(field,",");
	if(param.size() == 1) {
	  double dum = std::stoi(param[0]);
	  int nn = 1;
	  if(int(dum+0.5) == 0) nn = 0;
	  input.nodes.v.resize(nn, dum);
	} 
	else if(param.size() > 1){
	  input.nodes.v.resize(param.size());
	  for(int ii = 0; ii<(int)param.size(); ii++) input.nodes.v[ii] = std::stod(param[ii]);
	}
	set = true;
      }
      else if(key == "nodes_blong"){
	std::vector<std::string> param = strsplit(field,",");
	if(param.size() == 1) {
	  double dum = std::stoi(param[0]);
	  int nn = 1;
	  if(int(dum+0.5) == 0) nn = 0;
	  input.nodes.bl.resize(nn, dum);
	} 
	else if(param.size() > 1){
	  input.nodes.bl.resize(param.size());
	  for(int ii = 0; ii<(int)param.size(); ii++) input.nodes.bl[ii] = std::stod(param[ii]);
	}
	set = true;
      }
      else if(key == "nodes_bhor"){
	std::vector<std::string> param = strsplit(field,",");
	if(param.size() == 1) {
	  double dum = std::stoi(param[0]);
	  int nn = 1;
	  if(int(dum+0.5) == 0) nn = 0;
	  input.nodes.bh.resize(nn, dum);
	} 
	else if(param.size() > 1){
	  input.nodes.bh.resize(param.size());
	  for(int ii = 0; ii<(int)param.size(); ii++) input.nodes.bh[ii] = std::stod(param[ii]);
	}
	set = true;
      }
      else if(key == "nodes_azi"){
	std::vector<std::string> param = strsplit(field,",");
	if(param.size() == 1){
	  double dum = std::stoi(param[0]);
	  int nn = 1;
	  if(int(dum+0.5) == 0) nn = 0;
	  input.nodes.azi.resize(nn, dum);
	} 
	else if(param.size() > 1){
	  input.nodes.azi.resize(param.size());
	  for(int ii = 0; ii<(int)param.size(); ii++) input.nodes.azi[ii] = std::stod(param[ii]);
	}
	set = true;
      }
      else if(key == "invert_pgas_boundary"){
	std::vector<std::string> param = strsplit(field,",");
	if(param.size() == 1){
	  int dum = std::stoi(param[0]);
	  if(dum > 0) input.nodes.toinv[6] = 1;
	} 
	set = true;
      }
      else if(key == "nodes_vturb"){
	std::vector<std::string> param = strsplit(field,",");
	if(param.size() == 1){
	  double dum = std::stoi(param[0]);
	  int nn = 1;
	  if(int(dum+0.5) == 0) nn = 0;
	  input.nodes.vturb.resize(nn, dum);
	} 
	else if(param.size() > 1){
	  input.nodes.vturb.resize(param.size());
	  for(int ii = 0; ii<(int)param.size(); ii++) input.nodes.vturb[ii] = std::stod(param[ii]);
	}
	set = true;
      }
      else if(key == "") set = false;
      else{
	if(verbose) std::cout << "read_input: ignoring line: " << iline<<std::endl;
      }
      if(set && verbose) std::cout << "read_input: " << key << " -> "<<field << std::endl;
    }
  }

  if(input.nodes.temp.size() == 1) input.nodes.nnodes += input.nodes.temp[0];
  else input.nodes.nnodes += input.nodes.temp.size();
  
  if(input.nodes.vturb.size() == 1) input.nodes.nnodes += input.nodes.vturb[0];
  else input.nodes.nnodes += input.nodes.vturb.size();
  
  if(input.nodes.v.size() == 1) input.nodes.nnodes += input.nodes.v[0];
  else input.nodes.nnodes += input.nodes.v.size();
  
  if(input.nodes.bl.size() == 1) input.nodes.nnodes += input.nodes.bl[0];
  else input.nodes.nnodes += input.nodes.bl.size();
  
  if(input.nodes.bh.size() == 1) input.nodes.nnodes += input.nodes.bh[0];
  else input.nodes.nnodes += input.nodes.bh.size();
  
  if(input.nodes.azi.size() == 1) input.nodes.nnodes += input.nodes.azi[0];
  else input.nodes.nnodes += input.nodes.azi.size();

  if(input.nodes.toinv[6] > 0) input.nodes.nnodes += 1;
  
  if(verbose && (input.nodes.nnodes > 0)) {
    std::cout << "read_input: total number of nodes = "<<input.nodes.nnodes<<std::endl;
    std::cout << std::endl;
  }
  in.close();

  
  return input;
}
//
void read_lines(std::string filename, iput_t &input, bool verbose){
  std::string inam = "read_lines: ";
  //
  std::ifstream in(filename, std::ios::in | std::ios::binary);
  if (!in){
    //exit(0);
    return;
  }
  
  // 
  // Extract info from the string array
  //
  std::vector<line_t> line;
  std::string iline;
  
  //
  while(std::getline(in, iline)){
    if(iline == "" || iline[0] == '#') continue;

    // Convert to sstream
    std::stringstream stream(iline);
    line_t tmp = {0};
    std::string strtmp;
      
    // Check label
    stream >> strtmp;
    strcpy(tmp.label, (removeSpaces(strtmp)).c_str());
      
    bool cont = false;
    for(auto &it: input.ilines) if((it == strtmp) || (it == "all")) cont = true;
      
    if(!cont) continue;
      
    // element
    stream >> strtmp;
    strcpy(tmp.elem, strtmp.c_str());
      
    // ion
    stream >> tmp.ion;
      
    // anum
    stream >> tmp.anum;
    tmp.amass = phyc::AMASS[tmp.anum-1];
      
    // wav
    stream >> tmp.w0;
    tmp.w0 = inv_convl(tmp.w0); // Convert to vacuum
    tmp.nu0 = (phyc::CC * 1.0E8)  / tmp.w0; // in Hz (s^-1)
      
    // log gf
    stream >> tmp.gf;
    tmp.gf = pow(10.0, tmp.gf); // log(gf) -> gf
      
    // j_low & Jup
    stream >> tmp.Jlow;
    stream >> tmp.Jup;
      
    // g_low & gup
    stream >> tmp.Glow;
    stream >> tmp.Gup;
      
    // e_low & ionization energy
    stream >> tmp.e_low;
    tmp.e_up = 1.E8 / tmp.w0 * phyc::CM1_TO_EV + tmp.e_low;
    if(tmp.ion == 1) tmp.eion = phyc::EION1[tmp.anum-1] - tmp.e_low;
    else if(tmp.ion == 2) tmp.eion = phyc::EION2[tmp.anum-1] - tmp.e_low;
    else tmp.eion = 0;
    // Convert to ERGs
    tmp.e_up  *= phyc::EV;
    tmp.e_low *= phyc::EV;
    tmp.eion  *= phyc::EV;
      
    // gamma rad
    stream >> tmp.g_rad;
    if(tmp.g_rad == 0.0) tmp.g_rad = 0.22 / (tmp.w0*tmp.w0*1.e-16); // Gray (2005) eq. 11.13
    else tmp.g_rad = pow(10.0,tmp.g_rad);

    // gamma Stark (Vald gives log10(gamma_4/nne at 10000. [K]).
    stream >> tmp.g_str;
    if(tmp.g_str != 0.0) tmp.g_str = pow(10.0, tmp.g_str);
    else tmp.g_str = 0.0;
    
    // vDW & Barklem
    stream >> tmp.g_vdw;
      
    // INIT Barklem constants here
    if(tmp.g_vdw > 20.0){
      tmp.b_sig = int(tmp.g_vdw);
      tmp.b_alp = tmp.g_vdw - int(tmp.g_vdw);
    } else{
      tmp.b_sig = 0.0;
      tmp.b_alp = 0.0;
    }
      
    // width
    stream >> tmp.width;
      

    tmp.firsttime = 1;
      
    input.lines.push_back(tmp);
      
    if(verbose){
      if(input.lines.size() == 1) cout<<inam<<" Using input lines:"<<endl;
      fprintf(stdout, "%-15s %3s %2d %3d %10.4f %8.3f %5.2f %5.2f %5.2f %5.2f %8.3f %8.3f %8.3f\n",
	      tmp.label, tmp.elem, tmp.ion, tmp.anum, convl(tmp.w0), log10(tmp.gf),
	      tmp.Jlow, tmp.Jup, tmp.Glow, tmp.Gup, tmp.g_rad, tmp.g_str, tmp.g_vdw);
	
	
    }
    iline.clear();
  } // while(getline)
  
  in.close();
  
  if(input.lines.size() == 0){
    cerr << inam << "ERROR, input file does not contain lines present in "<< filename<<endl;
    exit(0);
  }
  
}


/* --- convl & inv_convl, adapted from MULTI, Carlsson 1986 --- */
double convl(double alamb){
  if(alamb < 2000.) return alamb;
  else return alamb/(1.0+2.735182e-4+131.4182/alamb/alamb+ 
		     2.76249e8/alamb/alamb/alamb/alamb);
}

double inv_convl(double lambda_air){
  // Init values
  double lambda_vacuum;
  if(lambda_air > 2000) lambda_vacuum = lambda_air / 1.00029;
  else lambda_vacuum = lambda_air;

  // Iterate
  double error = 1.0;
  while(error > 1.e-6){
    error = lambda_air - convl(lambda_vacuum);
    lambda_vacuum = lambda_vacuum + error / 1.0029;
  }
  return lambda_vacuum;
}

// int set_nodes(nodes_t &n, double min, double max, bool verbose){
//   int nn = n.nnodes;
//   string inam = "set_nodes: ";

//   /* --- resize array for the node type --- */
//   n.ntype.resize(nn, none_node);
//   int k = 0;
//   n.tosend = 0;
  
//   /* --- temp --- */
//   nn = n.temp.size();
//   if(nn > 0) n.temp_off = k;
//   else n.temp_off = -1;
  
//   if(nn == 1){
//     nn = n.temp[0];
//     n.temp.resize(nn);
//     equidist(n.temp, min, max);
//   }
//   for(int ii = 0; ii<nn; ii++) n.ntype[k++] = temp_node;
//   if(verbose) cout << inam << "Temp -> "<< formatVect<double>(n.temp)<<endl;

//   if(n.temp.size() > 0) n.toinv[0] = 1;
//   else{
//     n.toinv[0] = 0;
//     n.tosend += 1;
//   }
  
//   /* --- vlos --- */
//   nn = n.v.size();
//   if(nn > 0) n.v_off = k;
//   else n.v_off = -1;

//   if(nn == 1){
//     nn = n.v[0];
//     n.v.resize(nn);
//     equidist(n.v, min, max);
//   }
//   for(int ii = 0; ii<nn; ii++) n.ntype[k++] = v_node;
//   if(verbose) cout << inam << "Vlos -> "<< formatVect<double>(n.v)<<endl;
//   if(n.v.size() > 0) n.toinv[1] = 1;
//   else{
//     n.toinv[1] = 0;
//     n.tosend += 1;
//   }
  
//   /* --- vturb --- */
//   nn = n.vturb.size();
//   if(nn > 0) n.vturb_off = k;
//   else n.vturb_off = -1;

//   if(nn == 1){
//     nn = n.vturb[0];
//     n.vturb.resize(nn);
//     equidist(n.vturb, min, max);
//   }
//   for(int ii = 0; ii<nn; ii++) n.ntype[k++] = vturb_node;
//   if(verbose) cout << inam << "Vturb -> "<< formatVect<double>(n.vturb)<<endl;
//   if(n.vturb.size() > 0) n.toinv[2] = 1;
//   else{
//     n.toinv[2] = 0;
//     n.tosend += 1;
//   }
  
//   /* --- blong --- */
//   nn = n.bl.size();
//   if(nn > 0) n.bl_off = k;
//   else n.bl_off = -1;
  
//   if(nn == 1){
//     nn = n.bl[0];
//     n.bl.resize(nn);
//     equidist(n.bl, min, max);
//   }
//   for(int ii = 0; ii<nn; ii++) n.ntype[k++] = bl_node;
//   if(verbose) cout << inam << "Blong -> "<< formatVect<double>(n.bl)<<endl;
//   if(n.bl.size() > 0) n.toinv[3] = 1;
//   else{
//     n.toinv[3] = 0;
//     n.tosend += 1;
//   }
  
//   /* --- bhor --- */
//   nn = n.bh.size();
//   if(nn > 0) n.bh_off = k;
//   else n.bh_off = -1;
  
//   if(nn == 1){
//     nn = n.bh[0];
//     n.bh.resize(nn);
//     equidist(n.bh, min, max);
//   }
//   for(int ii = 0; ii<nn; ii++) n.ntype[k++] = bh_node;
//   if(verbose) cout << inam << "Bhor -> "<< formatVect<double>(n.bh)<<endl;
//   if(n.bh.size() > 0) n.toinv[4] = 1;
//   else{
//     n.toinv[4] = 0;
//     n.tosend += 1;
//   }
  
//   /* --- azi --- */
//   nn = n.azi.size();
//   if(nn > 0) n.azi_off = k;
//   else n.azi_off = -1;
  
//   if(nn == 1){
//     nn = n.azi[0];
//     n.azi.resize(nn);
//     equidist(n.azi, min, max);
//   }
//   for(int ii = 0; ii<nn; ii++) n.ntype[k++] = azi_node;
//   if(verbose) cout << inam << "Azi -> "<< formatVect<double>(n.azi)<<endl;
//   if(n.azi.size() > 0) n.toinv[5] = 1;
//   else{
//     n.toinv[5] = 0;
//     n.tosend += 1;
//   }
  
  
//   return n.nnodes;
  
// }

int set_nodes(nodes_t &n, vector<double> &itau, int dint, bool verbose){
  int nn = n.nnodes;
  string inam = "set_nodes: ";

  /* --- resize array for the node type --- */
  
  n.ntype.resize(nn, none_node);
  int k = 0;
  n.tosend = 0;
  
  /* --- temp --- */
  nn = n.temp.size();
  if(nn > 0) n.temp_off = k;
  else n.temp_off = -1;
  
  if(nn == 1){
    nn = n.temp[0];
    n.temp.resize(nn);
    equidist(n.temp, itau, ((nn >= 3) && (dint == 3))?true:false);
  }else if(nn > 1) for(int ii = 0; ii<nn; ii++) n.temp[ii] = nodeLocation(itau, n.temp[ii]);
 
  for(int ii = 0; ii<nn; ii++) n.ntype[k++] = temp_node;
  if(verbose) cout << inam << "Temp -> "<< formatVect<double>(n.temp)<<endl;

  if(n.temp.size() > 0) n.toinv[0] = 1;
  else{
    n.toinv[0] = 0;
    n.tosend += 1;
  }
  
  /* --- vlos --- */
  nn = n.v.size();
  if(nn > 0) n.v_off = k;
  else n.v_off = -1;

  if(nn == 1){
    nn = n.v[0];
    n.v.resize(nn);
    equidist(n.v, itau, ((nn >= 3) && (dint == 3))?true:false);
  }else if(nn > 1) for(int ii = 0; ii<nn; ii++) n.v[ii] = nodeLocation(itau, n.v[ii]);
  
  for(int ii = 0; ii<nn; ii++) n.ntype[k++] = v_node;
  if(verbose) cout << inam << "Vlos -> "<< formatVect<double>(n.v)<<endl;
  if(n.v.size() > 0) n.toinv[1] = 1;
  else{
    n.toinv[1] = 0;
    n.tosend += 1;
  }
  
  /* --- vturb --- */
  nn = n.vturb.size();
  if(nn > 0) n.vturb_off = k;
  else n.vturb_off = -1;

  if(nn == 1){
    nn = n.vturb[0];
    n.vturb.resize(nn);
    equidist(n.vturb, itau, ((nn >= 3) && (dint == 3))?true:false);
  }else if(nn > 1) for(int ii = 0; ii<nn; ii++) n.vturb[ii] = nodeLocation(itau, n.vturb[ii]);
  
  for(int ii = 0; ii<nn; ii++) n.ntype[k++] = vturb_node;
  if(verbose) cout << inam << "Vturb -> "<< formatVect<double>(n.vturb)<<endl;
  if(n.vturb.size() > 0) n.toinv[2] = 1;
  else{
    n.toinv[2] = 0;
    n.tosend += 1;
  }
  
  /* --- blong --- */
  nn = n.bl.size();
  if(nn > 0) n.bl_off = k;
  else n.bl_off = -1;
  
  if(nn == 1){
    nn = n.bl[0];
    n.bl.resize(nn);
    equidist(n.bl, itau, ((nn >=3) && (dint == 3))?true:false);
  }else if(nn > 1) for(int ii = 0; ii<nn; ii++) n.bl[ii] = nodeLocation(itau, n.bl[ii]);
  
  for(int ii = 0; ii<nn; ii++) n.ntype[k++] = bl_node;
  if(verbose) cout << inam << "Blong -> "<< formatVect<double>(n.bl)<<endl;
  if(n.bl.size() > 0) n.toinv[3] = 1;
  else{
    n.toinv[3] = 0;
    n.tosend += 1;
  }
  
  /* --- bhor --- */
  nn = n.bh.size();
  if(nn > 0) n.bh_off = k;
  else n.bh_off = -1;
  
  if(nn == 1){
    nn = n.bh[0];
    n.bh.resize(nn);
    equidist(n.bh, itau, ((nn >=3) && (dint == 3))?true:false);
  }else if(nn > 1) for(int ii = 0; ii<nn; ii++) n.bh[ii] = nodeLocation(itau, n.bh[ii]);

  
  for(int ii = 0; ii<nn; ii++) n.ntype[k++] = bh_node;
  if(verbose) cout << inam << "Bhor -> "<< formatVect<double>(n.bh)<<endl;
  if(n.bh.size() > 0) n.toinv[4] = 1;
  else{
    n.toinv[4] = 0;
    n.tosend += 1;
  }
  
  /* --- azi --- */
  nn = n.azi.size();
  if(nn > 0) n.azi_off = k;
  else n.azi_off = -1;
  
  if(nn == 1){
    nn = n.azi[0];
    n.azi.resize(nn);
    equidist(n.azi, itau, ((nn >=3) && (dint == 3))?true:false);
  }else if(nn > 1) for(int ii = 0; ii<nn; ii++) n.azi[ii] = nodeLocation(itau, n.azi[ii]);
  
  for(int ii = 0; ii<nn; ii++) n.ntype[k++] = azi_node;
  if(verbose) cout << inam << "Azi -> "<< formatVect<double>(n.azi)<<endl;
  if(n.azi.size() > 0) n.toinv[5] = 1;
  else{
    n.toinv[5] = 0;
    n.tosend += 1;
  }


  
  /* --- Pgas boundary --- */
  if(n.toinv[6] > 0){
    n.pgas_off = k;
    n.ntype[k] = pgas_node;
  }
  
  return n.nnodes;
  
}

void equidist(vector<double> &var, double min, double max){
  int n = var.size();

  if(n == 1) var[0] = 0.0;
  else{
    double dr = (max - min) / ((double)n - 1.0);
    
    var[0] = min;
    for(int ii = 0; ii<n; ii++) var[ii] = min + (double)ii * dr;
  }
  var[n-1] = max;
}

void equidist(vector<double> &var, vector<double> &itau, bool cent_grid){
  int n = (int)var.size();
  int n1 = (int)itau.size();
  double mmi = itau[0];
  double mma = itau[0];
  
  for(auto &it: itau){
    if(it > mma) mma = it;
    if(it < mmi) mmi = it;
  }
  

  if(n == 1) var[0] = 0.0;
  else{
    double dx = (mma-mmi) / (n-1.0);
    if(cent_grid){
      double ddx = 1.e10, odx = 0.0;
      while(fabs(ddx) > 1.e-5){
	odx = dx;
	mmi = itau[0]    + dx/2;
	mma = itau[n1-1] - dx/2;
	dx = (mma-mmi) / (n-1.0);
	ddx = dx - odx;
      }
    }
    
    for(int ii = 1; ii<n-1; ii++){
      double xloc = dx*ii + mmi;
      var[ii] = nodeLocation(itau, xloc);
      
    } // ii
    var[0] = mmi;
    var[n-1] = mma;
  } // else

}

double nodeLocation(vector<double> &itau, double iloc){

  /* --- get closest location to iloc in itau --- */

  int idx = -1, ndep = (int)itau.size();
  double imin = 1.e10;
  
  for(int zz=0; zz < ndep; zz++){
    double idif = abs(itau[zz] - iloc);
    if(idif < imin){
      imin = idif;
      idx = zz;
    }
  }

  return itau[idx];
}
