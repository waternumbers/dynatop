#include "helpers.h"

std::vector<hru> makeHRUs(Rcpp::DataFrame mdl,
			  std::vector<double> &s_sf, std::vector<double> &s_rz, std::vector<double> &s_uz, std::vector<double> &s_sz, std::vector<double> &q_sz, std::vector<double> &q_sf,
			  std::vector<double> &precip, std::vector<double> &pet,
			  std::vector<double> &q_sf_in, std::vector<double> &q_sz_in, std::vector<double> &q_sf_in_prev, std::vector<double> &q_sz_in_prev, 
			  std::vector<double> &r_sf_rz, std::vector<double> &r_rz_uz, std::vector<double> &r_uz_sz, std::vector<double> &e_a,
			  double const &Dt, double const &vtol, double const &qtol, int const &max_it){
  int nhru = mdl.nrows(); // number of HRUs
  std::vector<hru> hrus;
  Rcpp::NumericVector s_bar = mdl["s_bar"]; // average gradient
  Rcpp::NumericVector area = mdl["area"]; // surface area (plan)
  Rcpp::NumericVector width = mdl["width"]; // contour length of outflow
  Rcpp::List sf_desc = mdl["sf_cpp"]; // descriptions and parameters for surface store
  Rcpp::List sp_desc = mdl["sp_cpp"]; // descriptions and parameters for spill store
  Rcpp::List rz_desc = mdl["sf_cpp"]; // descriptions and parameters for root zone
  Rcpp::List uz_desc = mdl["sf_cpp"]; // descriptions and parameters for unsaturated zone
  Rcpp::List sz_desc = mdl["sz_cpp"]; // descriptions and parameters for saturated zone
  for(int ii=0; ii<nhru; ++ii){
    Rcpp::List tmp_sf = sf_desc[ii];
    Rcpp::List tmp_sz = sz_desc[ii];
    hrus.push_back( hru(s_sf[ii],s_rz[ii],s_uz[ii],s_sz[ii], // states passed - references
			q_sf[ii],q_sz[ii], q_sf_in[ii], q_sz_in[ii], // external fluxes - references
			precip[ii], pet[ii], e_a[ii], // external fluxes - references
			r_sf_sp[ii], r_sf_rz[ii], r_sp_rz[ii], r_rz_uz[ii], r_uz_sz[ii], // vertical fluxes - references
			tmp_sf["type_int"], tmp_sf["param"], // surface type and parameters passed explicitly
			tmp_sp["type_int"], tmp_sp["param"], // spill type and parameters passed explicitly
			tmp_rz["type_int"], tmp_rz["param"], // root zone type and parameters passed explicitly
			tmp_uz["type_int"], tmp_uz["param"], // unsaturated zone type and parameters passed explicitly
			tmp_sz["type_int"], tmp_sz["param"], // saturated zone type and parameters passed explicitly
			Dt, vtol,qtol,max_it) ); // simulation parameters passed as references
  }
  return( hrus );
}


outFlux::outFlux(std::vector<int> out_idx_, std::vector<int> idx_, std::vector<int> flux_type_):
  out_idx(out_idx_), idx(idx_), flux_type(flux_type_)
{};

void outFlux::apply(std::vector<double> &out, double &sc,
		    std::vector<double> &precip, std::vector<double> &pet, std::vector<double> &a_e,
		    std::vector<double> &q_sf, std::vector<double> &q_sf_in,
		    std::vector<double> &q_sz, std::vector<double> &q_sz_in,
		    std::vector<double> &s_sf, std::vector<double> &s_rz,
		    std::vector<double> &s_uz, std::vector<double> &s_sz,
		    std::vector<double> &r_sf_rz, std::vector<double> &r_rz_uz, std::vector<double> &r_uz_sz){
  for(uint ii=0; ii<out_idx.size(); ++ii){
    int &oi = out_idx[ii];
    int &i = idx[ii];
    int &fl = flux_type[ii];
    switch (fl) {
    case 1: // precip
      out[oi] += precip[i] / sc;
      break;
    case 2: // pet
      out[oi] += pet[i] / sc;
      break;
    case 3: // aet
      out[oi] += a_e[i] / sc;
      break;
    case 4: // q_sf
      out[oi] += q_sf[i] / sc;
      break;
    case 5: // q_sf_in
      out[oi] += q_sf_in[i] / sc;
      break;
    case 6: // q_sz
      out[oi] += q_sz[i] / sc;
      break;
    case 7: // q_sz_in
      out[oi] += q_sz_in[i] / sc;
      break;
    case 8: // s_sf
      out[oi] = s_sf[i];
      break;
    case 9: // s_rz
      out[oi] = s_rz[i];
      break;
    case 10: // s_uz
      out[oi] = s_uz[i];
      break;
    case 11: // s_sz
      out[oi] = s_sz[i];
      break;
    case 12: // r_sf_rz
      out[oi] += r_sf_rz[i] / sc;
      break;
    case 13: // r_rz_uz
      out[oi] += r_rz_uz[i] / sc;
      break;
    case 14: // q_uz_sz
      out[oi] += r_uz_sz[i] / sc;
      break;
    }
  }
}
