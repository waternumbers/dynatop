#include "helpers.h"

std::vector<hru> makeHRUs(Rcpp::List mdl){
  int nhru = mdl.size(); // number of HRUs
  std::vector<hru> hrus;
  
  for(int ii=0; ii<nhru; ++ii){
    Rcpp::List m = mdl[ii];
    Rcpp::NumericVector svec = m["states"];
    Rcpp::NumericVector pvec = m["properties"];
    Rcpp::List sf_list = m["sf"];
    Rcpp::List rz_list = m["rz"];
    Rcpp::List uz_list = m["uz"];
    Rcpp::List sz_list = m["sz"];
    Rcpp::List pcp_list = m["precip"];
    Rcpp::List pet_list = m["pet"];
    Rcpp::List q_sf_list = m["sf_flow_direction"];
    Rcpp::List q_sz_list = m["sz_flow_direction"];
    int id = m["id"];

    //svec = svec * pvec["area"];

    // all passed explicity, not by reference
    hrus.push_back( hru( id, //m["id"], // id passed explicitly
			 Rcpp::as<std::vector<double>>(svec),
			 Rcpp::as<std::vector<double>>(pvec),
			 Rcpp::as<int>(sf_list["type"]), Rcpp::as<std::vector<double>>(sf_list["parameters"]), // surface type and parameters
			 rz_list["parameters"], // root zone type and parameters passed explicitly
			 uz_list["parameters"], // unsaturated zone type and parameters passed explicitly
			 Rcpp::as<int>(sz_list["type"]), Rcpp::as<std::vector<double>>(sz_list["parameters"]), // saturated zone type and parameters passed explicitly
			 Rcpp::as<std::vector<int>>(pcp_list["idx"]), Rcpp::as<std::vector<double>>(pcp_list["fraction"]), // precipiataion inputs
			 Rcpp::as<std::vector<int>>(pet_list["idx"]), Rcpp::as<std::vector<double>>(pet_list["fraction"]), // pet inputs
			 Rcpp::as<std::vector<int>>(q_sf_list["id"]), Rcpp::as<std::vector<double>>(q_sf_list["fraction"]), // surface zone redistribution
			 Rcpp::as<std::vector<int>>(q_sz_list["id"]), Rcpp::as<std::vector<double>>(q_sz_list["fraction"]) // saturated zone redistribution
			 )
		    );
    
  }
  
  return(hrus);
}

Rcpp::List makeStateList(std::vector<hru> &hrus){
  int nhru = hrus.size(); // number of HRUs
  Rcpp::List state_list;
  for(int ii=0; ii<nhru; ++ii){
    Rcpp::NumericVector s = Rcpp::NumericVector::create(
							Rcpp::Named("s_sf", hrus[ii].s_sf / hrus[ii].area),
							Rcpp::Named("s_rz", hrus[ii].s_rz / hrus[ii].area),
							Rcpp::Named("s_uz", hrus[ii].s_uz / hrus[ii].area),
							Rcpp::Named("s_sz", hrus[ii].s_sz / hrus[ii].area),
							Rcpp::Named("q_sf", hrus[ii].q_sf / hrus[ii].area),
							Rcpp::Named("q_sz", hrus[ii].q_sz / hrus[ii].area));
    Rcpp::List L = Rcpp::List::create(Rcpp::Named("id") = hrus[ii].id , Rcpp::Named("states") = s);
    state_list.push_back(L);
  }
  return( state_list );
}



outFlux::outFlux(std::vector<int> out_idx_, std::vector<int> idx_, std::vector<int> flux_type_, std::vector<double> scale_, double nstep):
  out_idx(out_idx_), idx(idx_), flux_type(flux_type_), scale(scale_)
{
  for(long unsigned int ii=0; ii<scale.size(); ++ii){
    scale[ii] = scale[ii] / nstep;
  }
}

void outFlux::apply(std::vector<hru> &hrus, std::vector<double> &out){
  for(long unsigned int ii=0; ii<out_idx.size(); ++ii){
    
    int &oi = out_idx[ii];
    int &i = idx[ii];
    int &fl = flux_type[ii];
    double &sc = scale[ii];
    
    switch (fl) {
    case 1: // precip
      out[oi] += hrus[i].precip * sc;
      break;
    case 2: // pet
      out[oi] += hrus[i].pet * sc;
      break;
    case 3: // aet
      out[oi] += hrus[i].aet * sc;
      break;
    case 4: // q_sf
      out[oi] += hrus[i].q_sf * sc;
      break;
    case 5: // q_sf_in
      //Rcpp::Rcout << "q_sf_in " << i << " " <<  hrus[i].q_sf_in << std::endl;
      out[oi] += hrus[i].q_sf_in * sc;
      break;
    case 6: // q_sz
      out[oi] += hrus[i].q_sz * sc;
      break;
    case 7: // q_sz_in
      //Rcpp::Rcout << "q_sz_in " << i << " " <<  hrus[i].q_sz_in << std::endl;
      out[oi] += hrus[i].q_sz_in  * sc;
      break;
    case 8: // s_sf
      out[oi] = hrus[i].s_sf  * sc;
      break;
    case 9: // s_rz
      out[oi] = hrus[i].s_rz  * sc;
      break;
    case 10: // s_uz
      out[oi] = hrus[i].s_uz * sc;
      break;
    case 11: // s_sz
      out[oi] = hrus[i].s_sz * sc;
      break;
    case 12: // v_sf_rz
      out[oi] += hrus[i].v_sf_rz * sc;
      break;
    case 13: // v_rz_uz
      out[oi] += hrus[i].v_rz_uz * sc;
      break;
    case 14: // v_uz_sz
      out[oi] += hrus[i].v_uz_sz * sc;
      break;
    }
  }
}
