class flink {
private:
  // fraction and destination
  double &frc, &q_from, &q_to;
  
public:
  // constructor
  flink(double& frc_, double& q_from_, double& q_to_);
  void eval(); //double& q_sf_in, double& q_sz_in);
};
