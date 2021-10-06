#pragma once
#include "rednose/helpers/common_ekf.h"
extern "C" {
void live_update_3(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_4(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_9(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_10(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_12(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_32(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_13(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_14(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_19(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_H(double *in_vec, double *out_3486321576378820949);
void live_err_fun(double *nom_x, double *delta_x, double *out_6774737359060317795);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_8658247340784801970);
void live_H_mod_fun(double *state, double *out_6591325931598470811);
void live_f_fun(double *state, double dt, double *out_5129652178799571317);
void live_F_fun(double *state, double dt, double *out_2159681337421382206);
void live_h_3(double *state, double *unused, double *out_3516931370966504480);
void live_H_3(double *state, double *unused, double *out_4431687741703305444);
void live_h_4(double *state, double *unused, double *out_8262644902824909893);
void live_H_4(double *state, double *unused, double *out_6985753214206701721);
void live_h_9(double *state, double *unused, double *out_8786092844168326250);
void live_H_9(double *state, double *unused, double *out_2092456026997137100);
void live_h_10(double *state, double *unused, double *out_2234905382731837231);
void live_H_10(double *state, double *unused, double *out_46647814378935336);
void live_h_12(double *state, double *unused, double *out_7315926699410965140);
void live_H_12(double *state, double *unused, double *out_2308406059128783501);
void live_h_31(double *state, double *unused, double *out_3522737920654910723);
void live_H_31(double *state, double *unused, double *out_7631167739922692019);
void live_h_32(double *state, double *unused, double *out_7304326213810782349);
void live_H_32(double *state, double *unused, double *out_7162882495176947108);
void live_h_13(double *state, double *unused, double *out_1120575925540488616);
void live_H_13(double *state, double *unused, double *out_4492631502119963767);
void live_h_14(double *state, double *unused, double *out_8786092844168326250);
void live_H_14(double *state, double *unused, double *out_2092456026997137100);
void live_h_19(double *state, double *unused, double *out_4125313563859337844);
void live_H_19(double *state, double *unused, double *out_6738581207089662209);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}