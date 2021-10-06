#pragma once
#include "rednose/helpers/common_ekf.h"
extern "C" {
void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_err_fun(double *nom_x, double *delta_x, double *out_6401802909237175111);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_5339040875024180657);
void car_H_mod_fun(double *state, double *out_7553745406291488644);
void car_f_fun(double *state, double dt, double *out_5122005430319732036);
void car_F_fun(double *state, double dt, double *out_3354449505010087220);
void car_h_25(double *state, double *unused, double *out_2630864984075806183);
void car_H_25(double *state, double *unused, double *out_7302634080082881807);
void car_h_24(double *state, double *unused, double *out_1291237442698331645);
void car_H_24(double *state, double *unused, double *out_7184311969083448776);
void car_h_30(double *state, double *unused, double *out_2906059046360312072);
void car_H_30(double *state, double *unused, double *out_2311264830866301473);
void car_h_26(double *state, double *unused, double *out_2972226441198403016);
void car_H_26(double *state, double *unused, double *out_8608869852037722952);
void car_h_27(double *state, double *unused, double *out_7632538251736651828);
void car_H_27(double *state, double *unused, double *out_3598846818702926785);
void car_h_29(double *state, double *unused, double *out_7095864345180785352);
void car_H_29(double *state, double *unused, double *out_6135063599995133426);
void car_h_28(double *state, double *unused, double *out_4462485337195563497);
void car_H_28(double *state, double *unused, double *out_4244599626389607974);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}