/* -*- coding:utf-8-unix mode:c++ -*- */
#ifndef FOOT_H_
#define FOOT_H_
#include <iostream>
#include <cmath>
#include <hrpUtil/Eigen3d.h>
// #include "hrpsys/util/Hrpsys.h"
#include "PreviewController.h"

namespace rats
{
  class foot_guided_control_base
  {
  private:
    void calc_u(const double remain_t, const double ref_dcm, const double ref_zmp);
    void calc_u(const double remain_t, const double ref_dcm, const double ref_ccm, const double ref_zmp);
    void calc_x_k();
    void calc_hgain_list(const std::size_t N);
    void calc_refzmp_gain_list(const std::size_t N);
    double calc_hgain(const std::size_t N);
  protected:
    Eigen::Matrix<double, 2, 2> A; // state matrix
    Eigen::Matrix<double, 2, 1> b; // input matrix
    Eigen::Matrix<double, 2, 2> Phi; // convert matrix : pos, vel => dcm, ccm
    Eigen::Matrix<double, 2, 2> Phi_inv; // convert matrix : dcm, ccm => pos, vel
    Eigen::Matrix<double, 2, 1> x_k; // pos, vel
    Eigen::Matrix<double, 2, 1> w_k; // dcm, ccm
    hrp::dvector hgain_list; // hgain_list
    hrp::dvector refzmp_gain_list; // refzmp_gain_list (const)
    double u_k; // zmp
    double dt, g, t_max;
    double dz, xi, h, h_;
    size_t N_max;
  public:
    foot_guided_control_base(const double _dt,  const double _dz,
                             const double _g = DEFAULT_GRAVITATIONAL_ACCELERATION,
                             const double _t_max = 1.6)
      : dt(_dt), g(_g), t_max(_t_max), N_max(round(t_max / dt)),
        x_k(Eigen::Matrix<double, 2, 1>::Zero()), u_k(0.0)
    {
      set_param(_dz, N_max);
    }
    void update_x_k(const double remain_t, const double ref_dcm, const double ref_zmp);
    void update_x_k(const double remain_t, const double ref_dcm, const double ref_ccm, const double ref_zmp);
    void set_param(const double _dz, const size_t N_max);
    void get_pos (double& ret) { ret = x_k(0); }
    void get_vel (double& ret) { ret = x_k(1); }
    void get_zmp (double& ret) { ret = u_k; }
    double get_hgain(size_t i) { return hgain_list(i); }
    double get_refzmp_gain(size_t i) { return refzmp_gain_list(i); }
    double get_xi() { return xi; }
  };
}
#endif /*FOOT_H_*/
