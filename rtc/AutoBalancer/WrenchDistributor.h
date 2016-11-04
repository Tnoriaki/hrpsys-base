// -*- C++ -*-
/*!
 * @file  WrenchDistributor.h
 * @brief Wrench distribution
 * @date  $Date$
 *
 * $Id$
 */

#ifndef WRENCH_DISTRIBUTOR_H
#define WRENCH_DISTRIBUTOR_H

#include "../Stabilizer/ZMPDistributor.h"

// #ifdef USE_QPOASES
#include <qpOASES.hpp>
using namespace qpOASES;
// #endif

// #ifdef USE_QPOASES

// enum ee_name {RLEG, LLEG, RARM, LARM};
// enum state_type {UNILATERAL, ATTACHED};
// enum contact_type {FLOAT, STATIC, KINETIC};

class EndEffectorParam
{
    public:
    size_t state_dim;
    size_t c_dim;
    double weight; // standard 1
    hrp::Vector3 pos;
    hrp::Matrix33 rot;
    hrp::Vector3 e_vec; // selection vector (0,0,1,0,0,0) => unilateral / selection vector (0,0,0,0,0,0) => attached
    hrp::Vector3 mu_vec; // (mu_s, mu_k, mu_r)
    hrp::Vector3 move_vec; // ex) fix : (0,0,0) / sliding : (1 0 0),(0.5,0.5),(0,1,0) / float (0,0,1),(1,1,1)...
    hrp::dvector support_polygon_vec; // only rectangle
    hrp::dmatrix Cmat;
    hrp::dvector wrench;
    EndEffectorParam() : state_dim(6), c_dim(0), weight(1), e_vec(hrp::Vector3(0,0,1)), mu_vec(hrp::Vector3(0.3,0.1,0.03)), move_vec(hrp::Vector3(0,0,0)) {
        support_polygon_vec.resize(4);
        support_polygon_vec << 0.1,0.1,0.05,0.05;
        wrench = hrp::dvector::Zero(6);
    };
    EndEffectorParam(const size_t _state_dim) : state_dim(_state_dim), c_dim(0), weight(1),
                                                e_vec(hrp::Vector3(0,0,1)), mu_vec(hrp::Vector3(0.3,0.1,0.03)), move_vec(hrp::Vector3(0,0,0)) {
        support_polygon_vec.resize(4);
        support_polygon_vec << 0.1,0.1,0.05,0.05;
        wrench = hrp::dvector::Zero(6);
    };
    EndEffectorParam(const hrp::Vector3& _pos, const hrp::Matrix33& _rot) : state_dim(6), c_dim(0), weight(1),
                                                                                  pos(_pos), rot(_rot),
                                                                                  e_vec(hrp::Vector3(0,0,1)), mu_vec(hrp::Vector3(0.3,0.1,0.03)),
                                                                                  move_vec(hrp::Vector3(0,0,0))
    {
        support_polygon_vec.resize(4);
        support_polygon_vec << 0.1,0.1,0.05,0.05;
        wrench = hrp::dvector::Zero(6);
    };
    EndEffectorParam(const hrp::Vector3& _pos, const hrp::Matrix33& _rot, const size_t _state_dim) : state_dim(_state_dim), c_dim(0), weight(1),
                                                                                                           pos(_pos), rot(_rot),
                                                                                                           e_vec(hrp::Vector3(0,0,1)), mu_vec(hrp::Vector3(0.3,0.1,0.03)),
                                                                                                           move_vec(hrp::Vector3(0,0,0))
    {
        support_polygon_vec.resize(4);
        support_polygon_vec << 0.1,0.1,0.05,0.05;
        wrench = hrp::dvector::Zero(6);
    };
    void setEEParam(const hrp::Vector3& _pos, const hrp::Matrix33& _rot, const double _weight = 1.0){
        pos = _pos;
        rot = _rot;
        weight = _weight;
    }
    void setCCParam(const hrp::Vector3& _e_vec, const hrp::Vector3& _mu_vec, const hrp::dvector _support_polygon_vec){
        e_vec = _e_vec;
        mu_vec = _mu_vec;
        support_polygon_vec = _support_polygon_vec;
    }
    void calcStateConstraintsMatrix(hrp::dmatrix& C, hrp::Vector3& e_vec);
    void calcFrictionConstraintsMatrix(hrp::dmatrix& C, hrp::Vector3& mu_vec);
    void calcMomentumConstraintsMatrix(hrp::dmatrix& C, hrp::dvector& support_polygon_vec);
    virtual void calcConstraintsMatrix();
};

class WrenchDistributor : public EndEffectorParam
{
    private:
    double mass;
    double gravitational_acceleration;
    public:
    size_t ee_num;
    size_t states_dim;
    size_t constraints_dim;
    hrp::dvector weight_vector;
    hrp::dvector wrenches;
    hrp::dmatrix Hmat;
    hrp::dmatrix gvec;
    hrp::dmatrix Amat;
    hrp::dmatrix Phimat;
    hrp::dmatrix Ximat;
    hrp::Vector3 ref_cog;
    hrp::Vector3 ref_linear_momentum_rate;
    hrp::Vector3 ref_angular_momentum_rate;
    hrp::Vector3 linear_momentum_rate;
    hrp::Vector3 angular_momentum_rate;

    WrenchDistributor(const double _mass, const double _gravitational_acceleration)
        : mass(_mass), gravitational_acceleration(_gravitational_acceleration)
    {};
    void DistributeWrench(const hrp::Vector3& _ref_cog, const hrp::Vector3& _ref_linear_momentum_rate, const hrp::Vector3& _ref_angular_momentum_rate, std::map<std::string, EndEffectorParam>& _eeparam_map){
        ee_num = _eeparam_map.size();
        ref_cog = _ref_cog;
        ref_linear_momentum_rate = _ref_linear_momentum_rate;
        ref_angular_momentum_rate = _ref_angular_momentum_rate;
        calcAugmentedConstraintsMatrix(_eeparam_map);
        calcEvaluationFunctionMatrix(_eeparam_map);
        solveWrenchQP();
        calcMomentumRate();
        size_t index = 0;
        for ( std::map<std::string, EndEffectorParam>::iterator it = _eeparam_map.begin(); it != _eeparam_map.end(); it++ ){
            it->second.wrench.segment(0,it->second.state_dim) = wrenches.segment(index, it->second.state_dim);
            index += it->second.state_dim;
        }
    }
    void calcAugmentedConstraintsMatrix(std::map<std::string, EndEffectorParam>& eeparam_map);
    void calcEvaluationFunctionMatrix(const std::map<std::string, EndEffectorParam>& eeparam_map);
    void solveWrenchQP ();
    void calcMomentumRate();
    void printResult (const std::map<std::string, EndEffectorParam>& eeparam_map);
};

// #endif
#endif // WRENCH_DISTRIBUTOR_H
