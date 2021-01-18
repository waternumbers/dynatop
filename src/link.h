class flink {
private:
  // fraction and destination
  double &frc, &q_sf_out, &q_sz_out;
public:
  // constructor
  flink(double& frc_, double& q_sf_out_, double& q_sz_out_);
  void eval(double& q_sf_in, double& q_sz_in);
};
