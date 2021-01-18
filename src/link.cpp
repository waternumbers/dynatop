// a test for passing data to the hillslope class
#include "link.h"

flink::flink(double& frc_, double& q_sf_out_, double& q_sz_out_):
  frc(frc_), q_sf_out(q_sf_out_), q_sz_out(q_sz_out_){}

void flink::eval(double& q_sf_in, double& q_sz_in){
  q_sf_out += frc*q_sf_in;
  q_sz_out += frc*q_sz_in;
}
