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
    Rcpp::List q_sf_list = m["sf_links"];
    Rcpp::List q_sz_list = m["sz_links"];
    int id = m["id"];
    
    // hrus.push_back( hru( id, //m["id"], // id passed explicitly
    // 			 svec["s_sf"], svec["s_rz"], svec["s_uz"], svec["s_sz"], // states passed explicitly
    // 			 pvec["area"], pvec["s_bar"], pvec["width"], // properties passed explicity
    // 			 sf_list["type"], sf_list["parameters"], // surface type and parameters passed explicitly
    // 			 rz_list["parameters"], // root zone type and parameters passed explicitly
    // 			 uz_list["parameters"], // unsaturated zone type and parameters passed explicitly
    // 			 sz_list["type"], sz_list["parameters"], // saturated zone type and parameters passed explicitly
    // 			 pcp_list["id"], pcp_list["fraction"], // precipiataion inputs
    // 			 pet_list["id"], pet_list["fraction"], // pet inputs
    // 			 q_sf_list["id"], q_sf_list["fraction"], // surface zone redistribution
    // 			 q_sz_list["id"], q_sz_list["fraction"] // saturated zone redistribution
    // 			 )
    // 		    );
    hrus.push_back( hru( id, //m["id"], // id passed explicitly
			 Rcpp::as<double>(svec["s_sf"]), Rcpp::as<double>(svec["s_rz"]), Rcpp::as<double>(svec["s_uz"]), Rcpp::as<double>(svec["s_sz"]), // states passed explicitly
			 Rcpp::as<double>(pvec["area"]), Rcpp::as<double>(pvec["s_bar"]), Rcpp::as<double>(pvec["width"]), // properties passed explicity
			 Rcpp::as<int>(sf_list["type"]), Rcpp::as<std::vector<double>>(sf_list["parameters"]), // surface type and parameters passed explicitly
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
};

Rcpp::List makeStateList(std::vector<hru> &hrus){
  int nhru = hrus.size(); // number of HRUs
  Rcpp::List state_list;
  for(int ii=0; ii<nhru; ++ii){
    Rcpp::NumericVector s = Rcpp::NumericVector::create(
							Rcpp::Named("s_sf", hrus[ii].s_sf),
							Rcpp::Named("s_rz", hrus[ii].s_rz),
							Rcpp::Named("s_uz", hrus[ii].s_uz),
							Rcpp::Named("s_sz", hrus[ii].s_sz));
    Rcpp::List L = Rcpp::List::create(Rcpp::Named("id") = hrus[ii].id , Rcpp::Named("states") = s);
    state_list.push_back(L);
  }
  return( state_list );
};



outFlux::outFlux(std::vector<int> out_idx_, std::vector<int> idx_, std::vector<int> flux_type_):
  out_idx(out_idx_), idx(idx_), flux_type(flux_type_)
{};

void outFlux::apply(std::vector<hru> &hrus, std::vector<double> &out, double nstep){
  for(uint ii=0; ii<out_idx.size(); ++ii){
    int &oi = out_idx[ii];
    int &i = idx[ii];
    int &fl = flux_type[ii];
    switch (fl) {
    case 1: // precip
      out[oi] += hrus[i].precip * hrus[i].area / nstep;
      break;
    case 2: // pet
      out[oi] += hrus[i].pet  * hrus[i].area / nstep;
      break;
    case 3: // aet
      out[oi] += hrus[i].aet  * hrus[i].area / nstep;
      break;
    case 4: // q_sf
      out[oi] += hrus[i].q_sf * hrus[i].area / nstep;
      break;
    case 5: // q_sf_in
      out[oi] += hrus[i].q_sf_in * hrus[i].area / nstep;
      break;
    case 6: // q_sz
      out[oi] += hrus[i].q_sz  * hrus[i].area/ nstep;
      break;
    case 7: // q_sz_in
      out[oi] += hrus[i].q_sz_in * hrus[i].area / nstep;
      break;
    case 8: // s_sf
      out[oi] = hrus[i].s_sf * hrus[i].area / nstep;
      break;
    case 9: // s_rz
      out[oi] = hrus[i].s_rz * hrus[i].area / nstep;
      break;
    case 10: // s_uz
      out[oi] = hrus[i].s_uz  * hrus[i].area / nstep;
      break;
    case 11: // s_sz
      out[oi] = hrus[i].s_sz  * hrus[i].area / nstep;
      break;
    case 12: // r_sf_rz
      out[oi] += hrus[i].r_sf_rz  * hrus[i].area / nstep;
      break;
    case 13: // r_rz_uz
      out[oi] += hrus[i].r_rz_uz  * hrus[i].area / nstep;
      break;
    case 14: // r_uz_sz
      out[oi] += hrus[i].r_uz_sz  * hrus[i].area / nstep;
      break;
    }
  }
}
