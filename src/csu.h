#include <vector>
// channel HSU header
// really very simple but here to be developed
class csu {
private:
  // properties of the HSU
  const int id, p_idx; 
  const double area;
  // time step
  const double timestep;
  // local variables computed from properties
  double total_volume;
  int n_add; // number fo timesteps added
public:
  // constructor
  csu(const int& id_, const int& p_idx_,
      const double& area_, const double& timestep_);
  // return id
  int get_id();
  // return volume of water stored
  double get_vol();
  // return volume as average flow rate
  double get_inflow();
  // initialise to empty
  void initialise();
  // add volumes from precipiation, surface and saturated zone inputs
  void add_inflow(std::vector<double>& obs, std::vector<double>& il_sf_rec,
		  std::vector<double>& il_sz_rec, double& ipa);
};

