/* -*- coding:utf-8-unix; mode:c++; -*- */

#include "../SequencePlayer/interpolator.h"
#include "PreviewController.h"
#include "WrenchDistributor.h"

using namespace hrp;
using namespace rats;

// calculate time
#include <cstdio>
#include <sys/time.h>
#include <unistd.h>
#include <stdarg.h>

struct __mtimer__ {
  double start;
  char msg[100];
  bool is_sec_print;
  __mtimer__(const char* format, ...)
  __attribute__((format(printf, 2, 3)))
  {
    va_list args;
    va_start(args, format);
    vsnprintf(msg, sizeof(msg), format, args);
    va_end(args);

    start = sec();
  }
  ~__mtimer__() {
    fprintf(stderr, "[%s] => %.6f [ms]\n", msg, 1000*(sec() - start));
  }
  double sec() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec * 1e-6;
  }
  operator bool() { return false; }
};

#define mtimer(...) if(__mtimer__ __b__ = __mtimer__(__VA_ARGS__));else

// walk
void decide_contact_states(const size_t index, std::map<std::string, EndEffectorParam>& tmp_eeparam_map, std::map<std::string, EndEffectorParam>& eeparam_map)
{
    if ( index < 25 ) {
        tmp_eeparam_map["rleg"] = eeparam_map["rleg"];
        tmp_eeparam_map["lleg"] = eeparam_map["lleg"];
    } else if ( index < 100 ){
        tmp_eeparam_map["lleg"] = eeparam_map["lleg"];
    } else if ( index < 150 ){
        eeparam_map["rleg"].setEEParam(hrp::Vector3(0.15,-0.1,0.0), hrp::Matrix33::Identity());
        tmp_eeparam_map["rleg"] = eeparam_map["rleg"];
        tmp_eeparam_map["lleg"] = eeparam_map["lleg"];
    } else if (index < 225 ){
        tmp_eeparam_map["rleg"] = eeparam_map["rleg"];
    } else {
        eeparam_map["lleg"].setEEParam(hrp::Vector3(0.15,0.1,0.0), hrp::Matrix33::Identity());
        tmp_eeparam_map["rleg"] = eeparam_map["rleg"];
        tmp_eeparam_map["lleg"] = eeparam_map["lleg"];
    }
};

// shuffle motion
void decide_contact_states2(const double dt, const size_t index, std::map<std::string, EndEffectorParam>& tmp_eeparam_map, std::map<std::string, EndEffectorParam>& eeparam_map, interpolator* ip)
{
    double x,v,a;
    static double start = 0.0, goal = 0.15, goalv = 0.0;
    if ( index == 0 ){
        ip->set(&start);
        ip->setGoal(&goal, &goalv, 125*dt);
        ip->get(&x, &v, &a, true);
        eeparam_map["rleg"].move_vec = hrp::Vector3(1,0,0);
        eeparam_map["rleg"].setEEParam(hrp::Vector3(x,-0.1,0.0), hrp::Matrix33::Identity());
    } else if ( index < 125 ) {
        ip->get(&x, &v, &a, true);
        eeparam_map["rleg"].move_vec = hrp::Vector3(1,0,0);
        eeparam_map["rleg"].setEEParam(hrp::Vector3(x,-0.1,0.0), hrp::Matrix33::Identity());
    } else if ( index == 125 ) {
        ip->set(&start);
        ip->setGoal(&goal, &goalv, 125*dt);
        ip->get(&x, &v, &a, true);
        eeparam_map["rleg"].move_vec = hrp::Vector3(0,0,0);
        eeparam_map["rleg"].setEEParam(hrp::Vector3(goal,-0.1,0.0), hrp::Matrix33::Identity());
        eeparam_map["lleg"].move_vec = hrp::Vector3(1,0,0);
        eeparam_map["lleg"].setEEParam(hrp::Vector3(x,0.1,0.0), hrp::Matrix33::Identity());
    } else {
        ip->get(&x, &v, &a, true);
        eeparam_map["lleg"].move_vec = hrp::Vector3(1,0,0);
        eeparam_map["lleg"].setEEParam(hrp::Vector3(x,0.1,0.0), hrp::Matrix33::Identity());
    }
    tmp_eeparam_map["rleg"] = eeparam_map["rleg"];
    tmp_eeparam_map["lleg"] = eeparam_map["lleg"];
};

int main(int argc, char* argv[])
{
    //
    bool use_gnuplot = true;
    double dt = 0.04, max_tm = 10.0;
    std::queue<hrp::Vector3> ref_zmp_list;
    std::deque<double> tm_list;
    //
    double mass = 56;
    double g = 9.8066;
    double z_c = 0.8;
    WrenchDistributor wd(mass, g);
    //
    for (size_t i = 0; i < static_cast<size_t>(round(max_tm / dt)); i++) {
        double tmp_tm = i * dt;
        tm_list.push_back(tmp_tm);
        hrp::Vector3 v;
        double step_time = 5;
        double dspr = 0.2;
        if (tmp_tm < step_time * dspr) { // dsp
            v << 0, 0.1 * tmp_tm / (step_time * dspr) , 0;
        } else if (tmp_tm < step_time * (1 - dspr)) { // ssp
            v << 0, 0.1, 0;
        } else if (tmp_tm < step_time * (1 + dspr)) { //dsp
            v << 0.15 - 0.15 * (step_time * (1 + dspr) - tmp_tm) / ( 2 * step_time * dspr) , -0.1 + 0.2 * (step_time * (1 + dspr) - tmp_tm) / ( 2 * step_time * dspr), 0;
        } else if (tmp_tm < step_time * (2 - dspr)) {
            v << 0.15, -0.1, 0;
        } else {
            v << 0.15, (10 - tmp_tm) / ( step_time * dspr ) * -0.1, 0;
        }
        ref_zmp_list.push(v);
    }
    //
    interpolator* ip = new interpolator(1, dt, interpolator::HOFFARBIB);
    preview_dynamics_filter<extended_preview_control> df(dt, z_c, ref_zmp_list.front());
    std::string fname("/tmp/plot.dat");
    FILE* fp = fopen(fname.c_str(), "w");
    double cart_zmp[3], refzmp[3];
    bool r = true;
    size_t index = 0;
    // for eeparam
    std::map<std::string, EndEffectorParam> eeparam_map;
    EndEffectorParam rleg(hrp::Vector3(0, -0.1, 0), hrp::Matrix33::Identity());
    EndEffectorParam lleg(hrp::Vector3(0,  0.1, 0), hrp::Matrix33::Identity());
    eeparam_map.insert(std::pair<std::string, EndEffectorParam>("rleg", rleg));
    eeparam_map.insert(std::pair<std::string, EndEffectorParam>("lleg", lleg));
    while (r) {
        hrp::Vector3 p, x; // x cog
        std::vector<hrp::Vector3> qdata;
        r = df.update(p, x, qdata, ref_zmp_list.front(), qdata, !ref_zmp_list.empty());
        if (r) {
            index++;
            df.get_cart_zmp(cart_zmp);
            df.get_current_refzmp(refzmp);
            // wrench distributor
            std::map<std::string, EndEffectorParam> tmp_eeparam_map;
            double ratio;
            //// reference dynamics
            double omega2 = g / z_c;
            hrp::Vector3 com_pos = hrp::Vector3(x[0],x[1],0.8);
            hrp::Vector3 ref_linear_momentum_rate = mass * hrp::Vector3(omega2*(x[0]-cart_zmp[0]),omega2*(x[1]-cart_zmp[1]),0);
            hrp::Vector3 ref_angular_momentum_rate = hrp::Vector3(0,0,0);
            // decide_contact_states(index, tmp_eeparam_map, eeparam_map);
            decide_contact_states2(dt, index, tmp_eeparam_map, eeparam_map, ip);
            mtimer("time"){
                wd.DistributeWrench(com_pos, ref_linear_momentum_rate, ref_angular_momentum_rate, tmp_eeparam_map);
            }
            for ( std::map<std::string, EndEffectorParam>::iterator it = eeparam_map.begin(); it != eeparam_map.end(); it++ )
                it->second.wrench = hrp::dvector::Zero(it->second.state_dim);
            for ( std::map<std::string, EndEffectorParam>::iterator it = tmp_eeparam_map.begin(); it != tmp_eeparam_map.end(); it++ )
                eeparam_map[it->first] = it->second;
            // wd.printResult(eeparam_map);
            std::cerr << "at " << index * dt << " [s]"<< std::endl;
            // print
            fprintf(fp, "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",
                    tm_list[index],
                    cart_zmp[0], /* zmpx ;; this zmp is "zmp as a table-cart model" */
                    x[0], /* cogy */
                    refzmp[0], /* refzmpx */
                    cart_zmp[1], /* zmpy ;; this zmp is "zmp as a table-cart model" */
                    x[1], /* cogy */
                    refzmp[1], /* refzmpy */
                    ref_linear_momentum_rate[0] / mass, // reference cog_acc x
                    wd.linear_momentum_rate[0] / mass, // cog_acc x
                    ref_linear_momentum_rate[1] / mass, // reference cog_acc y
                    wd.linear_momentum_rate[1] / mass, // cog_acc y
                    ref_linear_momentum_rate[2] / mass, // reference cog_acc y
                    wd.linear_momentum_rate[2] / mass, // cog_acc y
                    ref_angular_momentum_rate[0], // reference angular momentum rate x
                    wd.angular_momentum_rate[0], // angular momentum rate x
                    ref_angular_momentum_rate[1], // reference angular momentum rate y
                    wd.angular_momentum_rate[1], // angular momentum rate y
                    ref_angular_momentum_rate[2], // reference angular momentum rate z
                    wd.angular_momentum_rate[2], // angular momentum rate z
                    eeparam_map["rleg"].wrench[0], // rleg force x
                    eeparam_map["lleg"].wrench[0], // lleg force x
                    eeparam_map["rleg"].wrench[1], // rleg force y
                    eeparam_map["lleg"].wrench[1], // lleg force y
                    eeparam_map["rleg"].wrench[2], // rleg force z
                    eeparam_map["lleg"].wrench[2] // lleg force z
                    );
        } else if ( !ref_zmp_list.empty() ) r = true;
        if (!ref_zmp_list.empty()) ref_zmp_list.pop();
    }
    fclose(fp);
    if (use_gnuplot) {
        FILE* gp[11];
        std::string titles[] = {"com pos X", "com pos Y", "com acc X", "com acc Y", "com acc Z", "angular momentum rate x", "angular momentum rate y", "angular momentum rate z", "force x", "force y", "force z"};
        for (size_t ii = 0; ii < 2; ii++) {
            gp[ii] = popen("gnuplot", "w");
            fprintf(gp[ii], "set title \"%s\"\n", titles[ii].c_str());
            fprintf(gp[ii], "plot \"%s\" using 1:%zu with lines title \"cart-table zmp\"\n", fname.c_str(), ( ii * 3 + 2));
            fprintf(gp[ii], "replot \"%s\" using 1:%zu with lines title \"cog\"\n", fname.c_str(), ( ii * 3 + 3));
            fprintf(gp[ii], "replot \"%s\" using 1:%zu with lines title \"refzmp\"\n", fname.c_str(), ( ii * 3 + 4));
            fflush(gp[ii]);
        }
        for (size_t ii = 0; ii < 3; ii++) {
            gp[ii] = popen("gnuplot", "w");
            fprintf(gp[ii], "set title \"%s\"\n", titles[2+ii].c_str());
            fprintf(gp[ii], "plot \"%s\" using 1:%zu with lines title \"reference cog acc\"\n", fname.c_str(), ( ii * 2 + 8 ));
            fprintf(gp[ii], "replot \"%s\" using 1:%zu with lines title \"cog acc\"\n", fname.c_str(), ( ii * 2 + 9 ));
            fflush(gp[ii]);
        }
        for (size_t ii = 0; ii < 3; ii++) {
            gp[ii] = popen("gnuplot", "w");
            fprintf(gp[ii], "set title \"%s\"\n", titles[5+ii].c_str());
            fprintf(gp[ii], "plot \"%s\" using 1:%zu with lines title \"reference angular momentum rate\"\n", fname.c_str(), ( ii * 2 + 14 ));
            fprintf(gp[ii], "replot \"%s\" using 1:%zu with lines title \"angular momentum rate\"\n", fname.c_str(), ( ii * 2 + 15 ));
            fflush(gp[ii]);
        }
        for (size_t ii = 0; ii < 3; ii++) {
            gp[ii] = popen("gnuplot", "w");
            fprintf(gp[ii], "set title \"%s\"\n", titles[8+ii].c_str());
            fprintf(gp[ii], "plot \"%s\" using 1:%zu with lines title \"rleg force\"\n", fname.c_str(), ( ii * 2 + 20 ));
            fprintf(gp[ii], "replot \"%s\" using 1:%zu with lines title \"lleg force\"\n", fname.c_str(), ( ii * 2 + 21 ));
            fflush(gp[ii]);
        }
        double tmp;
        std::cin >> tmp;
        for (size_t j = 0; j < 8; j++) pclose(gp[j]);
    }

    // // Initialization
    // double mass = 56;
    // double g = 9.8066;
    // WrenchDistributor wd(mass, g);

    // // Proc One Tic
    // //// End Effector
    // size_t ee_num = 2;
    // std::vector<EndEffectorParam> eeparam_vec;
    // eeparam_vec.resize(ee_num);
    // //// Constraints
    // eeparam_vec[0].setEEParam(hrp::Vector3(0,-0.105,0), hrp::Matrix33::Identity());
    // eeparam_vec[1].setEEParam(hrp::Vector3(0,0.105,0), hrp::Matrix33::Identity());
    // // eeparam_vec[2].setEEParam(hrp::Vector3(0.2,-0.105,0.5), hrp::Matrix33::Identity());
    // // eeparam_vec[3].setEEParam(hrp::Vector3(0.2,0.105,0.5), hrp::Matrix33::Identity());
    // for ( std::vector<EndEffectorParam>::iterator it = eeparam_vec.begin(); it != eeparam_vec.end(); it++ ){
    //     it->calcConstraintsMatrix();
    // }
    // ///// Reference Dynamics
    // hrp::Vector3 com_pos = hrp::Vector3(0.2,0,0.7);
    // hrp::Vector3 ref_lmr = hrp::Vector3(0.1,0,0);
    // hrp::Vector3 ref_amr = hrp::Vector3(0,0,0);
    // //// Wegiht
    // hrp::dvector weight_vector(1e-5*hrp::dvector::Ones(6*ee_num+6));
    // weight_vector.head(6) = hrp::dvector::Ones(6);
    // ////
    // wd.setParam(com_pos, ref_lmr, ref_amr, ee_num, eeparam_vec, weight_vector);
    // mtimer("time") {
    //     wd.solveWrenchQP();
    // }
    // wd.printResult();
}
