#include <Rcpp.h>
using namespace Rcpp;

//' cpp wrapper function for computation of redistribution matrices
//' 
//' @param dem Digital elevation model
//' @param grad average downslope gradient
//' @param land_area Area of land surface in pixel
//' @param channel_area Area of channel surface in pixel
//' @param channel_id UID of channel present in the pixel
//' @param hillslope_id hillslope class of each pixel
//' @param offset - difference between cell index of adjacent cells and current cell index - clockwise from top left
//' @param dx distance between cell centres - from top left in clockwise direction
//' @param cl contour length - from top left in a clockwise direction. The 9th value is used for cells split beteen land and channel
//' @param max_index maximum value of the hillslope and channel id's
//' @return list of hillslope and channel properties
//'

// function for external input
void f_ex(List h, NumericVector obs){
  //Rcout << "start ex" << "\n";
  
  NumericVector out = h["output"];
  for(int j=0;j<out.length();j++){
    out[j] = 0;
  }
  CharacterVector in_series = h["series"];
  CharacterVector in_var = in_series.names();
  for(int j=0;j<in_series.length();j++){
    //Rcout << j << " " << in_series[j] << " " << in_var[j] << "\n";
    std::string lbl = as<std::string>(in_series[j]);
    std::string var = as<std::string>(in_var[j]);
    out[var] = obs[lbl];
  }
  //Rcout << "leave ex" << "\n";
};

// function for mux
void f_mux(List h, List hru){
  CharacterVector lbl_in = h["uptree"];
  NumericVector weight = h["weight"];
  NumericVector out = h["output"];
  for(int j=0;j<out.length();j++){
    out[j] = 0;
  }
  for(int j=0;j<lbl_in.length();j++){
    //Rcout << j << " " << lbl_in[j]<< "\n";
    std::string lbl = as<std::string>(lbl_in[j]);
    List h_in = hru[lbl];
    NumericVector h_input = h_in["output"];
    NumericVector w_in = weight[j]*h_input;
    out = out + w_in;
  }
  // h["output"] = out;
  //Rcout << "leave mux" << "\n";
 };

//function for channel
void f_channel(List h, NumericVector input, List opt){
  NumericVector output = h["output"];
  double ts = opt["time_step"];
  NumericVector prop = h["prop"];
  output["qch"] = (input["lsz"] + input["lex"] + ts * prop["area"] * input["precip"])/(3600*ts);
  //h["output"]=input;
  //Rcout << "leave channel" << "\n";
};

double fode(double a,double b, double x0, double t){
  double x(1);
  if( b < 1e-10 ){b = 1e-10;}
  x = x0*exp(-b*t) + (a/b)*exp(-b*t);
  return x;
}

//function for hillslope
void f_hillslope(List h, NumericVector input, List opt){
  NumericVector prop = h["prop"];
  NumericVector state = h["state"];
  NumericVector param = h["param"];
  NumericVector output = h["output"];
  double ts = opt["time_step"];
  double omega = opt["omega"];
  double theta = opt["theta"];
  

  // solve surface excess
  double ex = fode(input["lex"]/(prop["area"]*ts),
		   1/param["tex"],
		   state["ex"],ts);
  double lex = ex + (input["lex"]/prop["area"]) - state["ex"];
  state["ex"] = ex;
  output["lex"] = lex*prop["area"];

  // solve the root zone
  double q_ex_rz = min(state["ex"]/ts, param["qex_max"]);
  double tilde_rz = fode(input["precip"] + q_ex_rz,
			 input["pet"]/param["srz_max"],
			 state["rz"],ts);
  state["rz"] = min(tilde_rz,param["srz_max"]);
  double q_rz_ex = 0 ;
  double q_rz_uz = tilde_rz - state["rz"];
  if( state["sz"] <= 0 ){
    q_rz_ex = q_rz_uz;
    q_rz_uz = 0;
  }

  // unsaturated zone
  double tilde_uz = fode( q_rz_uz/ts,
			  1/(state["sz"]*param["td"]),
			  state["uz"],ts);
  double q_uz_sz = state["uz"] + q_rz_uz - tilde_uz;
  state["uz"] = tilde_uz;

  // saturated zone
  double qbar = ( state["lsz"] + state["lsz_in"] + input["lsz"]
		  + max(state["lsz"], input["lsz"]) )/4;
  double cbar = (qbar * prop["delta_x"])/param["m"];
  double lambda = omega + theta*cbar*ts/prop["delta_x"];
  double lambda_prime = omega + (1-theta)*cbar*ts/prop["delta_x"];
  double k = lambda_prime*state["lsz"] +
    (1-lambda_prime)*min(state["lsz_in"],state["lsz_max"]) +
    cbar*q_uz_sz/prop["delta_x"];
  double lsz = (k - lambda_prime*min(input["lsz"],state["lsz_max"]))/lambda;
  lsz = min(lsz,state["lsz_max"]);

  // mass ballance on saturated zone
  double tilde_sz = state["sz"] + ts*(state["lsz"]+lsz)/2 - ts*(input["lsz"]+state["lsz_in"])/2;
  state["sz"] = max(tilde_sz,0);
  double q_sz_ex = state["sz"] - tilde_sz;
  state["lsz"] = lsz;
  state["lsz_in"] = input["lsz"];
  output["lsz"] = state["lsz"]*prop["area"];

  // correct stores
  state["ex"] = state["ex"] + q_rz_ex + q_sz_ex;
  if( state["sz"] <= 0 ){
    state["ex"] = state["ex"] + state["uz"];
    state["uz"] = 0;
  }
  
  //h["output"]=input;
  //Rcout << "leave hillslope" << "\n";
};


// [[Rcpp::export]]
void rcpp_dynatop(List hru, NumericVector obs, List opt){

  CharacterVector type_list = "mux";
  //CharacterVector::create(Named("mux","mux"),Named("leaf","leaf"));
  
  Rcout << "started" << "\n";
  int n = hru.size();
  int cnt = 0;
  
  Rcout << "n=" << n << "\n";
  for( int i=0; i<n; ++i){
    List h = hru[i];
    std::string type = h["type"];

    if( type == "ex" ){
      f_ex(h,obs);
      //NumericVector out = h["output"];
      //out["precip"] = 9.2;
    }
    
    // mux operator to join data together
    if( type == "mux" ){
      //Rcout << type << "\n";
      f_mux(h,hru);
    }

    // hillslope type for a dynamic topmodel hillslope
    if( type == "hillslope" ){
      std::string lbl = as<std::string>(h["uptree"]);
      List h_in = hru[lbl];
      NumericVector input = h_in["output"];
      f_hillslope(h,input,opt);
      //Rcout << type << "\n";
    }

    if( type == "channel" ){
      std::string lbl = as<std::string>(h["uptree"]);
      List h_in = hru[lbl];
      NumericVector input = h_in["output"];
      f_channel(h,input,opt);
      //Rcout << type << "\n";
    }
    
    cnt =  cnt + 1;
    //Rcout << type << "\n";

    // user interupt
    Rcpp::checkUserInterrupt();
  }

  Rcout << cnt << "\n";
  //return hru;

}

