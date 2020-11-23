#include "csu.h"

// channel HSU method funcions
// really very simple but here to be developed

// constructor
csu::csu(const int& id_, const int& p_idx_,
	 const double& area_, const double& timestep_):
  // set properties
  id{id_}, p_idx{p_idx_}, area{area_}, timestep{timestep_}
{
  // initialise local (non)constants
  // local constants used in the evaluation
  total_volume = 0.0;
  n_add = 0.0;
}

int csu::get_id(){return id;}

double csu::get_vol(){return total_volume;}

double csu::get_inflow(){return total_volume/(n_add*timestep);}

void csu::initialise(){total_volume = 0.0;n_add = 0;}

void csu::add_inflow(std::vector<double>& obs, std::vector<double>& il_sf_rec,
		     std::vector<double>& il_sz_rec, double& ipa){
  ipa = area*timestep*obs[p_idx];
  total_volume += ipa + il_sf_rec[id] + il_sz_rec[id];
  n_add += 1;
}

 
