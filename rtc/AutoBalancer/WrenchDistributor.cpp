/* -*- coding:utf-8-unix; mode:c++; -*- */

#include "WrenchDistributor.h"


void EndEffectorParam::calcStateConstraintsMatrix(hrp::dmatrix& C, hrp::Vector3& _e_vec)
{
    C = hrp::dmatrix::Zero(1,state_dim);
    C.block(0,0,1,3) << _e_vec(0), _e_vec(1), _e_vec(2);
}

void EndEffectorParam::calcFrictionConstraintsMatrix(hrp::dmatrix& C, hrp::Vector3& _mu_vec, hrp::Vector3& _move_vec)
{
    C = hrp::dmatrix::Zero(4, state_dim);
    C.block(0,0,4,2) <<
        1, 0,
        -1,0,
        0, 1,
        0,-1;
    if ( _move_vec.norm() != 0 ) _move_vec.normalize();
    for ( size_t i = 0; i < 2; i++ ) { // x or y
        if ( _move_vec(i) == 0 ) { // static
            C(2*i,2) = C(2*i+1,2) = _mu_vec(0);
        } else { // sliding
            C(2*i,2) = -_mu_vec(1)*(-_move_vec(i));
            C(2*i+1,2) = _mu_vec(1)*(-_move_vec(i));
        }
    }
}

void EndEffectorParam::calcMomentumConstraintsMatrix(hrp::dmatrix& C, hrp::dvector& _support_polygon_vec, hrp::Vector3& _mu_vec)
{
    // spv = (dx+,dx-,dy+,dy-)
    if ( state_dim == 6 ){ // TODO (line contact)
        C = hrp::dmatrix::Zero(6 ,state_dim);
        // Support Polygon Constraints (taux, tauy)
        C.block(0,0,4,state_dim) <<
            0,0,_support_polygon_vec(2),-1,0,0,
            0,0,_support_polygon_vec(3), 1,0,0,
            0,0,_support_polygon_vec(0),0,-1,0,
            0,0,_support_polygon_vec(1),0,1,0;
        // Rotation Slip Suppression (tauz)
        C.block(4,0,2,state_dim) <<
            0,0,_mu_vec(2),0,0,-1,
            0,0,_mu_vec(2),0,0,1;
    } else {
        C.resize(0,0); // TODO
    }
}

void EndEffectorParam::calcConstraintsMatrix()
{
    hrp::dmatrix C_state, C_friction, C_moment;
    calcStateConstraintsMatrix(C_state, e_vec);
    calcFrictionConstraintsMatrix(C_friction, mu_vec, move_vec);
    calcMomentumConstraintsMatrix(C_moment, support_polygon_vec, mu_vec);
    size_t cs_dim = C_state.rows();
    size_t cf_dim = C_friction.rows();
    size_t cm_dim = C_moment.rows();
    c_dim = cs_dim+cf_dim+cm_dim;
    Cmat = hrp::dmatrix::Zero(c_dim, state_dim);
    Cmat.block(0,0,cs_dim,state_dim) = C_state;
    Cmat.block(cs_dim,0,cf_dim,state_dim) = C_friction;
    if (state_dim == 6) Cmat.block(cs_dim+cf_dim,0,cm_dim,state_dim) = C_moment;
    // Local Frame => World Frame
    hrp::dmatrix Rmat(hrp::dmatrix::Zero(state_dim, state_dim));
    Rmat.block(0,0,3,3) = rot;
    if (state_dim == 6) Rmat.block(3,3,3,3) = rot; //TODO
    Cmat = Cmat * Rmat;
}

void ObjectParam::calcConstraintsMatrix(std::map<std::string, EndEffectorParam>& eeparam_map)
{
    states_dim = 0;
    // assuming eeparam_map[]->calcConstraintsMatrix
    for ( std::map<std::string, EndEffectorParam>::iterator it = eeparam_map.begin(); it != eeparam_map.end(); it++ ){
        states_dim += it->second.state_dim;
    }
    Cmat = hrp::dmatrix::Zero(11, states_dim); // TODO
    for ( std::vector<std::string>::iterator it = object_contact_eename_vec.begin(); it != object_contact_eename_vec.end(); it++ ){
        size_t index = eeparam_map[*it].index;
        size_t tmp_state_dim = eeparam_map[*it].state_dim;
        hrp::dmatrix tmpC_state, tmpC_friction;
        eeparam_map[*it].calcStateConstraintsMatrix(tmpC_state, e_vec);
        eeparam_map[*it].calcFrictionConstraintsMatrix(tmpC_friction, mu_vec, move_vec);
        Cmat.block(0, index, 1, tmp_state_dim) = tmpC_state * rot;
        Cmat.block(1, index, 4, tmp_state_dim) = tmpC_friction * rot;
    }
    hrp::dmatrix tmpC_moment;
    calcAugmentedMomentumConstraintsMatrix(tmpC_moment, support_polygon_vec, mu_vec, eeparam_map);
    Cmat.block(5, 0, 6, states_dim) = tmpC_moment;
}

void ObjectParam::calcAugmentedMomentumConstraintsMatrix(hrp::dmatrix& C, hrp::dvector& _support_polygon_vec, hrp::Vector3& _mu_vec, std::map<std::string, EndEffectorParam>& eeparam_map)
{
    C = hrp::dmatrix::Zero(6, states_dim);
    for ( std::vector<std::string>::iterator it = object_contact_eename_vec.begin(); it != object_contact_eename_vec.end(); it++ ){
        size_t index = eeparam_map[*it].index;
        size_t tmp_state_dim = eeparam_map[*it].state_dim;
        hrp::dmatrix tmpC_moment(hrp::dmatrix::Zero(6, 6));
        hrp::dmatrix crossM(hrp::hat(eeparam_map[*it].pos - pos));
        tmpC_moment <<
            0,0,_support_polygon_vec(2),-1,0,0,
            0,0,_support_polygon_vec(3), 1,0,0,
            0,0,_support_polygon_vec(0),0,-1,0,
            0,0,_support_polygon_vec(1),0,1,0,
            0,0,_mu_vec(2),0,0,-1,
            0,0,_mu_vec(2),0,0,1;
        for ( size_t i = 0; i < 3; i++ ){
            tmpC_moment.block(2*i, 0, 1, 3) -= crossM.block(i,0,1,3);
            tmpC_moment.block(2*i+1, 0, 1, 3) += crossM.block(i,0,1,3);
        }
        C.block(0, index, 6, tmp_state_dim) = tmpC_moment.block(0, 0, 6, tmp_state_dim);
    }
}

void WrenchDistributor::calcAugmentedConstraintsMatrix(std::map<std::string, EndEffectorParam>& eeparam_map)
{
    constraints_dim = 0;
    states_dim = 0;
    for ( std::map<std::string, EndEffectorParam>::iterator it = eeparam_map.begin(); it != eeparam_map.end(); it++ ){
        it->second.calcConstraintsMatrix();
        constraints_dim += it->second.c_dim;
        states_dim += it->second.state_dim;
    }
    weight_vector = hrp::dvector::Ones(6+states_dim);
    Amat = hrp::dmatrix::Zero(constraints_dim, states_dim);
    size_t current_dim = 0;
    size_t index = 0;
    for ( std::map<std::string, EndEffectorParam>::iterator it = eeparam_map.begin(); it != eeparam_map.end(); it++ ){
        weight_vector.segment( 6 + index , 3 ) << 1e-5, 1e-5, 1e-5;
        weight_vector.segment( 6 + index , it->second.state_dim ) *= it->second.weight;
        Amat.block(current_dim, index, it->second.c_dim, it->second.state_dim) = it->second.Cmat;
        current_dim += it->second.c_dim;
        it->second.index = index;
        index += it->second.state_dim;
    }
}

void WrenchDistributor::calcAugmentedConstraintsMatrix(std::map<std::string, EndEffectorParam>& eeparam_map, ObjectParam& oparam)
{
    constraints_dim = 0;
    states_dim = 0;
    for ( std::map<std::string, EndEffectorParam>::iterator it = eeparam_map.begin(); it != eeparam_map.end(); it++ ){
        it->second.calcConstraintsMatrix();
        constraints_dim += it->second.c_dim;
        states_dim += it->second.state_dim;
    }
    weight_vector = hrp::dvector::Ones(6+states_dim);
    constraints_dim += 11; // TODO
    Amat = hrp::dmatrix::Zero(constraints_dim, states_dim);
    size_t current_dim = 0;
    size_t index = 0;
    for ( std::map<std::string, EndEffectorParam>::iterator it = eeparam_map.begin(); it != eeparam_map.end(); it++ ){
        weight_vector.segment( 6 + index , 3 ) << 1e-5, 1e-5, 1e-5;
        weight_vector.segment( 6 + index , it->second.state_dim ) *= it->second.weight;
        Amat.block(current_dim, index, it->second.c_dim, it->second.state_dim) = it->second.Cmat;
        current_dim += it->second.c_dim;
        it->second.index = index;
        index += it->second.state_dim;
    }
    oparam.calcConstraintsMatrix(eeparam_map);
    Amat.block(current_dim, 0, oparam.c_dim, states_dim) = oparam.Cmat;
}

void WrenchDistributor::calcEvaluationFunctionMatrix(const std::map<std::string, EndEffectorParam>& eeparam_map)
{
    size_t Phirows = 6 + states_dim;
    size_t Phicols = states_dim;
    hrp::dmatrix Wmat(hrp::dmatrix::Zero(Phirows, Phirows));
    Ximat = hrp::dmatrix::Zero(Phirows, 1);
    Phimat = hrp::dmatrix::Zero(Phirows, Phicols);
    Wmat = weight_vector.asDiagonal();
    Ximat.block(0,0,3,1) = ref_linear_momentum_rate + hrp::Vector3(0,0,mass*gravitational_acceleration);
    Ximat.block(3,0,3,1) = ref_angular_momentum_rate;
    for ( std::map<std::string, EndEffectorParam>::const_iterator it = eeparam_map.begin(); it != eeparam_map.end(); it++ ){
        Phimat.block(0, it->second.index, 3,3) = hrp::dmatrix::Identity(3,3);
        Phimat.block(3, it->second.index, 3,3) << hrp::hat(it->second.pos - ref_cog);
        if ( it->second.state_dim == 6 ) Phimat.block(3, it->second.index + 3,3, 3) = hrp::dmatrix::Identity(3,3); // TODO
    }
    Phimat.block(6,0,Phirows-6,Phicols) = hrp::dmatrix::Identity(Phirows-6,Phicols);
    Hmat = Phimat.transpose()*Wmat*Phimat;
    gvec = - Phimat.transpose()*Wmat*Ximat;
}

void WrenchDistributor::solveWrenchQP ()
{
    size_t e_dim = 0;
    real_t* H = new real_t[states_dim*states_dim];
    real_t* g = new real_t[states_dim];
    real_t* A = new real_t[(constraints_dim+e_dim)*states_dim];
    real_t* lb = new real_t[states_dim];
    real_t* ub = new real_t[states_dim];
    real_t* lbA = new real_t[constraints_dim+e_dim];
    real_t* ubA = new real_t[constraints_dim+e_dim];
    for (size_t i = 0; i < states_dim; i++) {
        for (size_t j = 0; j < states_dim; j++) {
            H[i*states_dim+j] = Hmat(i,j);
        }
        g[i] = gvec(i);
        lb[i] = -1e10;
        ub[i] = 1e10;
    }
    for (size_t i = 0; i < constraints_dim; i++) {
        for (size_t j = 0; j < states_dim; j++){
            A[i*states_dim+j] = Amat(i,j);
        }
        lbA[i] = 0;
        ubA[i] = 1e10;
    }
    for (size_t i = 0; i < e_dim; i++){
        size_t index = i + constraints_dim;
        for (size_t j = 0; j < states_dim; j++){
            A[index*states_dim+j] = Phimat(i,j);
        }
        lbA[index] = Ximat(i);
        ubA[index] = Ximat(i);
    }
    QProblem example( states_dim, constraints_dim + e_dim );
    //
    Options options;
    options.printLevel = PL_LOW;
    options.initialStatusBounds = ST_INACTIVE;
    options.numRefinementSteps = 1;
    options.enableCholeskyRefactorisation = 1;
    //
    example.setOptions( options );
    int nWSR = 1000;
    example.init(H,g,A,lb,ub,lbA,ubA,nWSR);
    real_t* xOpt = new real_t[states_dim];
    example.getPrimalSolution( xOpt );
    size_t count = 0;
    wrenches.resize(states_dim);
    for ( size_t i = 0; i < states_dim; i++ )
        wrenches(i) = xOpt[i];
    delete[] H;
    delete[] g;
    delete[] A;
    delete[] lb;
    delete[] ub;
    delete[] lbA;
    delete[] ubA;
    delete[] xOpt;
};

void WrenchDistributor::calcMomentumRate ()
{
    linear_momentum_rate = (Phimat * wrenches).segment(0,3) - hrp::Vector3(0,0,mass * gravitational_acceleration);
    angular_momentum_rate = (Phimat * wrenches).segment(3,3);
}

void WrenchDistributor::printResult (const std::map<std::string, EndEffectorParam>& eeparam_map)
{
    std::cerr << "------------------------------------------" << std::endl;
    std::cerr << "CoM" << std::endl;
    std::cerr << ref_cog.transpose() << std::endl;
    std::cerr << "------------------------------------------" << std::endl;
    std::cerr << "reference linear momentum" << std::endl;
    std::cerr << ref_linear_momentum_rate.transpose() << std::endl;
    std::cerr << "reference angular momentum" << std::endl;
    std::cerr << ref_angular_momentum_rate.transpose() << std::endl;
    std::cerr << "------------------------------------------" << std::endl;
    std::cerr << "linear momentum" << std::endl;
    std::cerr << linear_momentum_rate.transpose() << std::endl;
    std::cerr << "angluar momentum" << std::endl;
    std::cerr << angular_momentum_rate.transpose() << std::endl;
    std::cerr << "wrenches" << std::endl;
    for ( std::map<std::string, EndEffectorParam>::const_iterator it = eeparam_map.begin(); it != eeparam_map.end(); it++ )
        std::cerr << it->second.wrench.transpose() << std::endl;
}
