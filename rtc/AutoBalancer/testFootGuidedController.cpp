/* -*- coding:utf-8-unix; mode:c++; -*- */

#include "matplotlibcpp.h"
#include "PreviewController.h"
#include "FootGuidedController.h"

using namespace hrp;
using namespace rats;
namespace plt = matplotlibcpp;

struct gait_parameter /* for test */
{
  double tm;
  hrp::Vector3 ref_zmp;
  gait_parameter (const double _tm, const hrp::Vector3& _ref_zmp)
    : tm(_tm), ref_zmp(_ref_zmp) {};
};

#include<cstdio>

int main(int argc, char* argv[])
{
  /* this is c++ version example of test-preview-filter1-modified in euslib/jsk/preview.l*/
  bool use_gnuplot = true;
  if (argc >= 2) {
      if ( std::string(argv[1])== "--use-gnuplot" ) {
          use_gnuplot = (std::string(argv[2])=="true");
      }
  }

  double dt = 0.004, max_tm = 8.0;
  std::queue<hrp::Vector3> ref_zmp_list;
  std::queue<hrp::Vector3> ref_dcm_list;
  std::deque<double> tm_list;
  double switch_tm_list[4] = { 2.00, 4.00, 6.00, 8.00 };
  for (size_t i = 0; i < 4; i ++)
    // switch_tm_list[i] += 0.04;
    switch_tm_list[i] += 0.016; // TODO ????
  for (size_t i = 0; i < static_cast<size_t>(round(max_tm / dt)); i++) {
    double tmp_tm = i * dt;
    tm_list.push_back(tmp_tm);
    hrp::Vector3 v;
    hrp::Vector3 v_;
    if (tmp_tm == 0){
      v << 0, 0.0, 0;
      v_ << 0, -0.02, 0;
    } else if (tmp_tm < 2) {
      v << 0, -0.02, 0;
      v_ << -0.1, 0.02, 0;
    } else if (tmp_tm < 4) {
      v << -0.1, 0.02, 0;
      v_ << 0.1, 0.00, 0;
    } else if (tmp_tm < 6) {
      v << 0.1, 0.0, 0;
      v_ << 0.0, -0.01, 0.0;
    } else {
      v << 0.0, -0.01, 0.0;
      v_ << 0.0, -0.01, 0.0;
    }
    ref_zmp_list.push(v);
    ref_dcm_list.push(v_);
  }

  foot_guided_control_base fgx(dt, 0.8);
  foot_guided_control_base fgy(dt, 0.8);
  preview_dynamics_filter<preview_control> pc(dt, 0.8, ref_zmp_list.front());
  preview_dynamics_filter<extended_preview_control> df(dt, 0.8, ref_zmp_list.front());
  std::string fname("/tmp/plot.dat");
  FILE* fp = fopen(fname.c_str(), "w");
  double cart_zmp[3], refzmp[3];
  bool r = true;
  size_t index = 3;
  while (r) {
    hrp::Vector3 p, x;
    std::vector<hrp::Vector3> qdata;
    r = df.update(p, x, qdata, ref_zmp_list.front(), qdata, !ref_zmp_list.empty());
    if (r) {
      index++;
      df.get_cart_zmp(cart_zmp);
      df.get_current_refzmp(refzmp);
      double remain_t = 0;
      for (size_t i = 0; i < 4; i++){
          if (tm_list[index] < switch_tm_list[i]) {
              remain_t = switch_tm_list[i] - tm_list[index];
              break;
          }
      }
      hrp::Vector3 ref_dcm = ref_dcm_list.front();
      fgx.update_x_k(round(remain_t / dt), ref_dcm[0], refzmp[0]);
      fgy.update_x_k(round(remain_t / dt), ref_dcm[1], refzmp[1]);
      ref_dcm_list.pop();
      double pos[2];
      double zmp[2];
      fgx.get_pos(pos[0]);
      fgx.get_zmp(zmp[0]);
      fgy.get_pos(pos[1]);
      fgy.get_zmp(zmp[1]);
      fprintf(fp, "%f %f %f %f %f %f %f %f %f %f %f\n",
              tm_list[index],
              cart_zmp[0], /* zmpx ;; this zmp is "zmp as a table-cart model" */
              x[0], /* cogy */
              refzmp[0], /* refzmpx */
              pos[0], // pos (foot guided control)
              zmp[0], // input (from foot guided control)
              cart_zmp[1], /* zmpy ;; this zmp is "zmp as a table-cart model" */
              x[1], /* cogy */
              refzmp[1], /* refzmpy */
              pos[1], // pos (foot guided control)
              zmp[1] // input zmp (foot guided control)
              );
    } else if ( !ref_zmp_list.empty() ) r = true;
    if (!ref_zmp_list.empty()) ref_zmp_list.pop();
  }
  fclose(fp);
  FILE* gp[3];
  std::string titles[2] = {"X", "Y"};
  for (size_t ii = 0; ii < 2; ii++) {
    gp[ii] = popen("gnuplot", "w");
    fprintf(gp[ii], "set title \"%s\"\n", titles[ii].c_str());
    fprintf(gp[ii], "plot \"%s\" using 1:%zu with lines title \"cart-table zmp\"\n", fname.c_str(), ( ii * 5 + 2));
    fprintf(gp[ii], "replot \"%s\" using 1:%zu with lines title \"cog\"\n", fname.c_str(), ( ii * 5 + 3));
    fprintf(gp[ii], "replot \"%s\" using 1:%zu with lines title \"refzmp\"\n", fname.c_str(), ( ii * 5 + 4));
    fprintf(gp[ii], "replot \"%s\" using 1:%zu with lines title \"cog of foot guided control\"\n", fname.c_str(), ( ii * 5 + 5));
    fprintf(gp[ii], "replot \"%s\" using 1:%zu with lines title \"zmp of foot guided control\"\n", fname.c_str(), ( ii * 5 + 6));
    fflush(gp[ii]);
  }
  double tmp;
  std::cin >> tmp;
  for (size_t j = 0; j < 2; j++) pclose(gp[j]);
  return 0;
}
