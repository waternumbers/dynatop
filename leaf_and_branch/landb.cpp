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
// [[Rcpp::export]]
List rcpp_dynatop(List hru){
  
  Rcout << "started" << "\n";
  int n = hru.size();

  for( int i=0; i<n; i++){
    List h = hru[i];
    CharacterVector type = h["type"];

    
    Rcout << type << "\n";
  }

  return hru;

}
