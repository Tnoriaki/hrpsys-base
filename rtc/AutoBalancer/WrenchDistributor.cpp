/* -*- coding:utf-8-unix; mode:c++; -*- */

#include "WrenchDistributor.h"


void EndEffectorParam::calcStateConstraintsMatrix(hrp::dmatrix& C)
{
    C = hrp::dmatrix::Zero(1,state_dim);
    C.block(0,0,1,3) << e_vec(0), e_vec(1), e_vec(2);
}

void EndEffectorParam::calcFrictionConstraintsMatrix(hrp::dmatrix& C)
{
    if ( move_vec.norm() == 0 ){ // static
        C = hrp::dmatrix::Zero(4, state_dim);
        C.block(0,0,1,3) <<  1, 0, mu_vec(0);
        C.block(1,0,1,3) << -1, 0, mu_vec(0);
        C.block(2,0,1,3) <<  0, 1, mu_vec(0);
        C.block(3,0,1,3) <<  0,-1, mu_vec(0);
    } else if ( move_vec(2) == 0 ) { // sliding
        move_vec.normalize();
        C = hrp::dmatrix::Zero(4, state_dim);
        C.block(0,0,1,3) <<  1, 0, -mu_vec(1)*(-move_vec(0));
        C.block(1,0,1,3) << -1, 0,  mu_vec(1)*(-move_vec(0));
        C.block(2,0,1,3) <<  0, 1, -mu_vec(1)*(-move_vec(1));
        C.block(3,0,1,3) <<  0,-1,  mu_vec(1)*(-move_vec(1));
    } else { // float
        C = hrp::dmatrix::Zero(2 * state_dim, state_dim);
        C.block(0,0,state_dim,state_dim) = hrp::dmatrix::Identity(state_dim,state_dim);
        C.block(state_dim,0,state_dim,state_dim) = -hrp::dmatrix::Identity(state_dim,state_dim);
    }
}

void EndEffectorParam::calcMomentumConstraintsMatrix(hrp::dmatrix& C)
{
    // spv = (dx+,dx-,dy+,dy-)
    if ( state_dim == 6 ){ // TODO (line contact)
        C = hrp::dmatrix::Zero(6 ,state_dim);
        // Support Polygon Constraints (taux, tauy)
        C.block(0,0,4,state_dim) <<
            0,0,support_polygon_vec(2),-1,0,0,
            0,0,support_polygon_vec(3), 1,0,0,
            0,0,support_polygon_vec(0),0,-1,0,
            0,0,support_polygon_vec(1),0,1,0;
        // Rotation Slip Suppression (tauz)
        C.block(4,0,2,state_dim) <<
            0,0,mu_vec(2),0,0,1,
            0,0,mu_vec(2),0,0,-1;
    } else {
        C.resize(0,0); // TODO
    }
}

void EndEffectorParam::calcConstraintsMatrix()
{
    hrp::dmatrix C_state, C_friction, C_moment;
    calcStateConstraintsMatrix(C_state);
    calcFrictionConstraintsMatrix(C_friction);
    calcMomentumConstraintsMatrix(C_moment);
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
    weight_vector.head(6) = 1e5*hrp::dvector::Ones(6);
    Amat = hrp::dmatrix::Zero(constraints_dim, states_dim);
    size_t current_dim = 0;
    size_t index = 0;
    for ( std::map<std::string, EndEffectorParam>::iterator it = eeparam_map.begin(); it != eeparam_map.end(); it++ ){
        weight_vector.segment( 6 + index , 3 ) << 1e-5, 1e-5, 1e-5;
        weight_vector.segment( 6 + index , it->second.state_dim ) *= it->second.weight;
        Amat.block(current_dim, index, it->second.c_dim, it->second.state_dim) = it->second.Cmat;
        current_dim += it->second.c_dim;
        index += it->second.state_dim;
    }
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
    size_t index = 0;
    for ( std::map<std::string, EndEffectorParam>::const_iterator it = eeparam_map.begin(); it != eeparam_map.end(); it++ ){
        Phimat.block(0, index, 3,3) = hrp::dmatrix::Identity(3,3);
        Phimat.block(3, index, 3,3) << hrp::hat(it->second.pos - ref_cog);
        if ( it->second.state_dim == 6 ) Phimat.block(3, index + 3,3, 3) = hrp::dmatrix::Identity(3,3); // TODO
        index += it->second.state_dim;
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
    int nWSR = 100;
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
