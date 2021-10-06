#include "car.h"

namespace {
#define DIM 8
#define EDIM 8
#define MEDIM 8
typedef void (*Hfun)(double *, double *, double *);

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}
const static double MAHA_THRESH_25 = 3.8414588206941227;
const static double MAHA_THRESH_24 = 5.991464547107981;
const static double MAHA_THRESH_30 = 3.8414588206941227;
const static double MAHA_THRESH_26 = 3.8414588206941227;
const static double MAHA_THRESH_27 = 3.8414588206941227;
const static double MAHA_THRESH_29 = 3.8414588206941227;
const static double MAHA_THRESH_28 = 5.991464547107981;

/******************************************************************************
 *                      Code generated with sympy 1.7.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_6401802909237175111) {
   out_6401802909237175111[0] = delta_x[0] + nom_x[0];
   out_6401802909237175111[1] = delta_x[1] + nom_x[1];
   out_6401802909237175111[2] = delta_x[2] + nom_x[2];
   out_6401802909237175111[3] = delta_x[3] + nom_x[3];
   out_6401802909237175111[4] = delta_x[4] + nom_x[4];
   out_6401802909237175111[5] = delta_x[5] + nom_x[5];
   out_6401802909237175111[6] = delta_x[6] + nom_x[6];
   out_6401802909237175111[7] = delta_x[7] + nom_x[7];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_5339040875024180657) {
   out_5339040875024180657[0] = -nom_x[0] + true_x[0];
   out_5339040875024180657[1] = -nom_x[1] + true_x[1];
   out_5339040875024180657[2] = -nom_x[2] + true_x[2];
   out_5339040875024180657[3] = -nom_x[3] + true_x[3];
   out_5339040875024180657[4] = -nom_x[4] + true_x[4];
   out_5339040875024180657[5] = -nom_x[5] + true_x[5];
   out_5339040875024180657[6] = -nom_x[6] + true_x[6];
   out_5339040875024180657[7] = -nom_x[7] + true_x[7];
}
void H_mod_fun(double *state, double *out_7553745406291488644) {
   out_7553745406291488644[0] = 1.0;
   out_7553745406291488644[1] = 0.0;
   out_7553745406291488644[2] = 0.0;
   out_7553745406291488644[3] = 0.0;
   out_7553745406291488644[4] = 0.0;
   out_7553745406291488644[5] = 0.0;
   out_7553745406291488644[6] = 0.0;
   out_7553745406291488644[7] = 0.0;
   out_7553745406291488644[8] = 0.0;
   out_7553745406291488644[9] = 1.0;
   out_7553745406291488644[10] = 0.0;
   out_7553745406291488644[11] = 0.0;
   out_7553745406291488644[12] = 0.0;
   out_7553745406291488644[13] = 0.0;
   out_7553745406291488644[14] = 0.0;
   out_7553745406291488644[15] = 0.0;
   out_7553745406291488644[16] = 0.0;
   out_7553745406291488644[17] = 0.0;
   out_7553745406291488644[18] = 1.0;
   out_7553745406291488644[19] = 0.0;
   out_7553745406291488644[20] = 0.0;
   out_7553745406291488644[21] = 0.0;
   out_7553745406291488644[22] = 0.0;
   out_7553745406291488644[23] = 0.0;
   out_7553745406291488644[24] = 0.0;
   out_7553745406291488644[25] = 0.0;
   out_7553745406291488644[26] = 0.0;
   out_7553745406291488644[27] = 1.0;
   out_7553745406291488644[28] = 0.0;
   out_7553745406291488644[29] = 0.0;
   out_7553745406291488644[30] = 0.0;
   out_7553745406291488644[31] = 0.0;
   out_7553745406291488644[32] = 0.0;
   out_7553745406291488644[33] = 0.0;
   out_7553745406291488644[34] = 0.0;
   out_7553745406291488644[35] = 0.0;
   out_7553745406291488644[36] = 1.0;
   out_7553745406291488644[37] = 0.0;
   out_7553745406291488644[38] = 0.0;
   out_7553745406291488644[39] = 0.0;
   out_7553745406291488644[40] = 0.0;
   out_7553745406291488644[41] = 0.0;
   out_7553745406291488644[42] = 0.0;
   out_7553745406291488644[43] = 0.0;
   out_7553745406291488644[44] = 0.0;
   out_7553745406291488644[45] = 1.0;
   out_7553745406291488644[46] = 0.0;
   out_7553745406291488644[47] = 0.0;
   out_7553745406291488644[48] = 0.0;
   out_7553745406291488644[49] = 0.0;
   out_7553745406291488644[50] = 0.0;
   out_7553745406291488644[51] = 0.0;
   out_7553745406291488644[52] = 0.0;
   out_7553745406291488644[53] = 0.0;
   out_7553745406291488644[54] = 1.0;
   out_7553745406291488644[55] = 0.0;
   out_7553745406291488644[56] = 0.0;
   out_7553745406291488644[57] = 0.0;
   out_7553745406291488644[58] = 0.0;
   out_7553745406291488644[59] = 0.0;
   out_7553745406291488644[60] = 0.0;
   out_7553745406291488644[61] = 0.0;
   out_7553745406291488644[62] = 0.0;
   out_7553745406291488644[63] = 1.0;
}
void f_fun(double *state, double dt, double *out_5122005430319732036) {
   out_5122005430319732036[0] = state[0];
   out_5122005430319732036[1] = state[1];
   out_5122005430319732036[2] = state[2];
   out_5122005430319732036[3] = state[3];
   out_5122005430319732036[4] = state[4];
   out_5122005430319732036[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_5122005430319732036[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_5122005430319732036[7] = state[7];
}
void F_fun(double *state, double dt, double *out_3354449505010087220) {
   out_3354449505010087220[0] = 1;
   out_3354449505010087220[1] = 0;
   out_3354449505010087220[2] = 0;
   out_3354449505010087220[3] = 0;
   out_3354449505010087220[4] = 0;
   out_3354449505010087220[5] = 0;
   out_3354449505010087220[6] = 0;
   out_3354449505010087220[7] = 0;
   out_3354449505010087220[8] = 0;
   out_3354449505010087220[9] = 1;
   out_3354449505010087220[10] = 0;
   out_3354449505010087220[11] = 0;
   out_3354449505010087220[12] = 0;
   out_3354449505010087220[13] = 0;
   out_3354449505010087220[14] = 0;
   out_3354449505010087220[15] = 0;
   out_3354449505010087220[16] = 0;
   out_3354449505010087220[17] = 0;
   out_3354449505010087220[18] = 1;
   out_3354449505010087220[19] = 0;
   out_3354449505010087220[20] = 0;
   out_3354449505010087220[21] = 0;
   out_3354449505010087220[22] = 0;
   out_3354449505010087220[23] = 0;
   out_3354449505010087220[24] = 0;
   out_3354449505010087220[25] = 0;
   out_3354449505010087220[26] = 0;
   out_3354449505010087220[27] = 1;
   out_3354449505010087220[28] = 0;
   out_3354449505010087220[29] = 0;
   out_3354449505010087220[30] = 0;
   out_3354449505010087220[31] = 0;
   out_3354449505010087220[32] = 0;
   out_3354449505010087220[33] = 0;
   out_3354449505010087220[34] = 0;
   out_3354449505010087220[35] = 0;
   out_3354449505010087220[36] = 1;
   out_3354449505010087220[37] = 0;
   out_3354449505010087220[38] = 0;
   out_3354449505010087220[39] = 0;
   out_3354449505010087220[40] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_3354449505010087220[41] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_3354449505010087220[42] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_3354449505010087220[43] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_3354449505010087220[44] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_3354449505010087220[45] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_3354449505010087220[46] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_3354449505010087220[47] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_3354449505010087220[48] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_3354449505010087220[49] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_3354449505010087220[50] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_3354449505010087220[51] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_3354449505010087220[52] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_3354449505010087220[53] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_3354449505010087220[54] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_3354449505010087220[55] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_3354449505010087220[56] = 0;
   out_3354449505010087220[57] = 0;
   out_3354449505010087220[58] = 0;
   out_3354449505010087220[59] = 0;
   out_3354449505010087220[60] = 0;
   out_3354449505010087220[61] = 0;
   out_3354449505010087220[62] = 0;
   out_3354449505010087220[63] = 1;
}
void h_25(double *state, double *unused, double *out_2630864984075806183) {
   out_2630864984075806183[0] = state[6];
}
void H_25(double *state, double *unused, double *out_7302634080082881807) {
   out_7302634080082881807[0] = 0;
   out_7302634080082881807[1] = 0;
   out_7302634080082881807[2] = 0;
   out_7302634080082881807[3] = 0;
   out_7302634080082881807[4] = 0;
   out_7302634080082881807[5] = 0;
   out_7302634080082881807[6] = 1;
   out_7302634080082881807[7] = 0;
}
void h_24(double *state, double *unused, double *out_1291237442698331645) {
   out_1291237442698331645[0] = state[4];
   out_1291237442698331645[1] = state[5];
}
void H_24(double *state, double *unused, double *out_7184311969083448776) {
   out_7184311969083448776[0] = 0;
   out_7184311969083448776[1] = 0;
   out_7184311969083448776[2] = 0;
   out_7184311969083448776[3] = 0;
   out_7184311969083448776[4] = 1;
   out_7184311969083448776[5] = 0;
   out_7184311969083448776[6] = 0;
   out_7184311969083448776[7] = 0;
   out_7184311969083448776[8] = 0;
   out_7184311969083448776[9] = 0;
   out_7184311969083448776[10] = 0;
   out_7184311969083448776[11] = 0;
   out_7184311969083448776[12] = 0;
   out_7184311969083448776[13] = 1;
   out_7184311969083448776[14] = 0;
   out_7184311969083448776[15] = 0;
}
void h_30(double *state, double *unused, double *out_2906059046360312072) {
   out_2906059046360312072[0] = state[4];
}
void H_30(double *state, double *unused, double *out_2311264830866301473) {
   out_2311264830866301473[0] = 0;
   out_2311264830866301473[1] = 0;
   out_2311264830866301473[2] = 0;
   out_2311264830866301473[3] = 0;
   out_2311264830866301473[4] = 1;
   out_2311264830866301473[5] = 0;
   out_2311264830866301473[6] = 0;
   out_2311264830866301473[7] = 0;
}
void h_26(double *state, double *unused, double *out_2972226441198403016) {
   out_2972226441198403016[0] = state[7];
}
void H_26(double *state, double *unused, double *out_8608869852037722952) {
   out_8608869852037722952[0] = 0;
   out_8608869852037722952[1] = 0;
   out_8608869852037722952[2] = 0;
   out_8608869852037722952[3] = 0;
   out_8608869852037722952[4] = 0;
   out_8608869852037722952[5] = 0;
   out_8608869852037722952[6] = 0;
   out_8608869852037722952[7] = 1;
}
void h_27(double *state, double *unused, double *out_7632538251736651828) {
   out_7632538251736651828[0] = state[3];
}
void H_27(double *state, double *unused, double *out_3598846818702926785) {
   out_3598846818702926785[0] = 0;
   out_3598846818702926785[1] = 0;
   out_3598846818702926785[2] = 0;
   out_3598846818702926785[3] = 1;
   out_3598846818702926785[4] = 0;
   out_3598846818702926785[5] = 0;
   out_3598846818702926785[6] = 0;
   out_3598846818702926785[7] = 0;
}
void h_29(double *state, double *unused, double *out_7095864345180785352) {
   out_7095864345180785352[0] = state[1];
}
void H_29(double *state, double *unused, double *out_6135063599995133426) {
   out_6135063599995133426[0] = 0;
   out_6135063599995133426[1] = 1;
   out_6135063599995133426[2] = 0;
   out_6135063599995133426[3] = 0;
   out_6135063599995133426[4] = 0;
   out_6135063599995133426[5] = 0;
   out_6135063599995133426[6] = 0;
   out_6135063599995133426[7] = 0;
}
void h_28(double *state, double *unused, double *out_4462485337195563497) {
   out_4462485337195563497[0] = state[5];
   out_4462485337195563497[1] = state[6];
}
void H_28(double *state, double *unused, double *out_4244599626389607974) {
   out_4244599626389607974[0] = 0;
   out_4244599626389607974[1] = 0;
   out_4244599626389607974[2] = 0;
   out_4244599626389607974[3] = 0;
   out_4244599626389607974[4] = 0;
   out_4244599626389607974[5] = 1;
   out_4244599626389607974[6] = 0;
   out_4244599626389607974[7] = 0;
   out_4244599626389607974[8] = 0;
   out_4244599626389607974[9] = 0;
   out_4244599626389607974[10] = 0;
   out_4244599626389607974[11] = 0;
   out_4244599626389607974[12] = 0;
   out_4244599626389607974[13] = 0;
   out_4244599626389607974[14] = 1;
   out_4244599626389607974[15] = 0;
}
#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}




}
extern "C" {

void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
}
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<2, 3, 0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
}
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
}
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
}
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
}
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
}
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<2, 3, 0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
}
void car_err_fun(double *nom_x, double *delta_x, double *out_6401802909237175111) {
  err_fun(nom_x, delta_x, out_6401802909237175111);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_5339040875024180657) {
  inv_err_fun(nom_x, true_x, out_5339040875024180657);
}
void car_H_mod_fun(double *state, double *out_7553745406291488644) {
  H_mod_fun(state, out_7553745406291488644);
}
void car_f_fun(double *state, double dt, double *out_5122005430319732036) {
  f_fun(state,  dt, out_5122005430319732036);
}
void car_F_fun(double *state, double dt, double *out_3354449505010087220) {
  F_fun(state,  dt, out_3354449505010087220);
}
void car_h_25(double *state, double *unused, double *out_2630864984075806183) {
  h_25(state, unused, out_2630864984075806183);
}
void car_H_25(double *state, double *unused, double *out_7302634080082881807) {
  H_25(state, unused, out_7302634080082881807);
}
void car_h_24(double *state, double *unused, double *out_1291237442698331645) {
  h_24(state, unused, out_1291237442698331645);
}
void car_H_24(double *state, double *unused, double *out_7184311969083448776) {
  H_24(state, unused, out_7184311969083448776);
}
void car_h_30(double *state, double *unused, double *out_2906059046360312072) {
  h_30(state, unused, out_2906059046360312072);
}
void car_H_30(double *state, double *unused, double *out_2311264830866301473) {
  H_30(state, unused, out_2311264830866301473);
}
void car_h_26(double *state, double *unused, double *out_2972226441198403016) {
  h_26(state, unused, out_2972226441198403016);
}
void car_H_26(double *state, double *unused, double *out_8608869852037722952) {
  H_26(state, unused, out_8608869852037722952);
}
void car_h_27(double *state, double *unused, double *out_7632538251736651828) {
  h_27(state, unused, out_7632538251736651828);
}
void car_H_27(double *state, double *unused, double *out_3598846818702926785) {
  H_27(state, unused, out_3598846818702926785);
}
void car_h_29(double *state, double *unused, double *out_7095864345180785352) {
  h_29(state, unused, out_7095864345180785352);
}
void car_H_29(double *state, double *unused, double *out_6135063599995133426) {
  H_29(state, unused, out_6135063599995133426);
}
void car_h_28(double *state, double *unused, double *out_4462485337195563497) {
  h_28(state, unused, out_4462485337195563497);
}
void car_H_28(double *state, double *unused, double *out_4244599626389607974) {
  H_28(state, unused, out_4244599626389607974);
}
void car_predict(double *in_x, double *in_P, double *in_Q, double dt) {
  predict(in_x, in_P, in_Q, dt);
}
void car_set_mass(double x) {
  set_mass(x);
}
void car_set_rotational_inertia(double x) {
  set_rotational_inertia(x);
}
void car_set_center_to_front(double x) {
  set_center_to_front(x);
}
void car_set_center_to_rear(double x) {
  set_center_to_rear(x);
}
void car_set_stiffness_front(double x) {
  set_stiffness_front(x);
}
void car_set_stiffness_rear(double x) {
  set_stiffness_rear(x);
}
}

const EKF car = {
  .name = "car",
  .kinds = { 25, 24, 30, 26, 27, 29, 28 },
  .feature_kinds = {  },
  .f_fun = car_f_fun,
  .F_fun = car_F_fun,
  .err_fun = car_err_fun,
  .inv_err_fun = car_inv_err_fun,
  .H_mod_fun = car_H_mod_fun,
  .predict = car_predict,
  .hs = {
    { 25, car_h_25 },
    { 24, car_h_24 },
    { 30, car_h_30 },
    { 26, car_h_26 },
    { 27, car_h_27 },
    { 29, car_h_29 },
    { 28, car_h_28 },
  },
  .Hs = {
    { 25, car_H_25 },
    { 24, car_H_24 },
    { 30, car_H_30 },
    { 26, car_H_26 },
    { 27, car_H_27 },
    { 29, car_H_29 },
    { 28, car_H_28 },
  },
  .updates = {
    { 25, car_update_25 },
    { 24, car_update_24 },
    { 30, car_update_30 },
    { 26, car_update_26 },
    { 27, car_update_27 },
    { 29, car_update_29 },
    { 28, car_update_28 },
  },
  .Hes = {
  },
  .sets = {
    { "mass", car_set_mass },
    { "rotational_inertia", car_set_rotational_inertia },
    { "center_to_front", car_set_center_to_front },
    { "center_to_rear", car_set_center_to_rear },
    { "stiffness_front", car_set_stiffness_front },
    { "stiffness_rear", car_set_stiffness_rear },
  },
  .extra_routines = {
  },
};

ekf_init(car);
