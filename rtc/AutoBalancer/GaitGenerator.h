/* -*- mode:c++ -*- */
#ifndef GAITGENERATOR_H
#define GAITGENERATOR_H
#include "PreviewController.h"
#include "../ImpedanceController/RatsMatrix.h"
#include "interpolator.h"
#include <vector>
#include <queue>
#include <boost/assign.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/shared_ptr.hpp>

namespace rats
{
    void cycloid_midpoint (hrp::Vector3& ret,
                           const double ratio, const hrp::Vector3& start,
                           const hrp::Vector3& goal, const double height,
                           const double default_top_ratio = 0.5);
    void multi_mid_coords (coordinates& mid_coords, const std::vector<coordinates>& cs);

    enum orbit_type {SHUFFLING, CYCLOID, RECTANGLE, STAIR, CYCLOIDDELAY, CYCLOIDDELAYKICK, CROSS};
    enum leg_type {RLEG, LLEG, RARM, LARM, BOTH, ALL};

    struct step_node
    {
        leg_type l_r;
        coordinates worldcoords;
        double step_height, step_time, toe_angle, heel_angle;
        step_node () : l_r(RLEG), worldcoords(coordinates()),
                       step_height(), step_time(),
                       toe_angle(), heel_angle(){};
        step_node (const leg_type _l_r, const coordinates& _worldcoords,
                   const double _step_height, const double _step_time,
                   const double _toe_angle, const double _heel_angle)
            : l_r(_l_r), worldcoords(_worldcoords),
              step_height(_step_height), step_time(_step_time),
              toe_angle(_toe_angle), heel_angle(_heel_angle) {};
        step_node (const std::string& _l_r, const coordinates& _worldcoords,
                   const double _step_height, const double _step_time,
                   const double _toe_angle, const double _heel_angle)
            : l_r((_l_r == "rleg") ? RLEG :
                  (_l_r == "rarm") ? RARM :
                  (_l_r == "larm") ? LARM :
                  LLEG), worldcoords(_worldcoords),
              step_height(_step_height), step_time(_step_time),
              toe_angle(_toe_angle), heel_angle(_heel_angle) {};
        friend std::ostream &operator<<(std::ostream &os, const step_node &sn)
        {
            os << "footstep" << std::endl;
            os << "  name = [" << ((sn.l_r==LLEG)?std::string("lleg"):
                                   (sn.l_r==RARM)?std::string("rarm"):
                                   (sn.l_r==LARM)?std::string("larm"):
                                   std::string("rleg")) << "]" << std::endl;
            os << "  pos =";
            os << (sn.worldcoords.pos).format(Eigen::IOFormat(Eigen::StreamPrecision, 0, ", ", ", ", "", "", " [", "]")) << std::endl;
            os << "  rot =";
            os << (sn.worldcoords.rot).format(Eigen::IOFormat(Eigen::StreamPrecision, 0, ", ", "", " [", "]")) << std::endl;
            os << "  step_height = " << sn.step_height << "[m], step_time = " << sn.step_time << "[s], "
               << "toe_angle = " << sn.toe_angle << "[deg], heel_angle = " << sn.heel_angle << "[deg]";
            return os;
        };
    };

    /* footstep parameter */
    struct footstep_parameter
    {
        /* translate pos is translate position of a leg from default foot_midcoords
         *   vector -> (list rleg-pos[mm] lleg-pos[mm] )
         */
        std::vector<hrp::Vector3> leg_default_translate_pos;
        /* stride params indicate max stride ( [mm], [mm], [deg] ) */
        double stride_fwd_x, stride_y, stride_theta, stride_bwd_x;
        footstep_parameter (const std::vector<hrp::Vector3>& _leg_pos,
                            const double _stride_fwd_x, const double _stride_y, const double _stride_theta, const double _stride_bwd_x)
            : leg_default_translate_pos(_leg_pos),
              stride_fwd_x(_stride_fwd_x), stride_y(_stride_y), stride_theta(_stride_theta), stride_bwd_x(_stride_bwd_x)  {};
    };

    /* velocity parameter for velocity mode */
    struct velocity_mode_parameter
    {
        /* velocity is [mm/s], [mm/s], [deg/s] */
        double velocity_x, velocity_y, velocity_theta;
        void set (const double _vx, const double _vy, const double _vth)
        {
            velocity_x = _vx;
            velocity_y = _vy;
            velocity_theta = _vth;
        };
        velocity_mode_parameter ()
            :velocity_x(0), velocity_y(0), velocity_theta(0) {};
    };

    /* Phase name of toe heel contact.
     *   SOLE0 : Sole contact (start). Foot angle is 0.
     *   SOLE2TOE : Transition of foot angle (0 -> toe_angle).
     *   TOE2SOLE : Transition of foot angle (toe_angle -> 0).
     *   SOLE1 : Foot_angle is 0.
     *   SOLE2HEEL : Transition of foot angle (0 -> -1 * heel_angle).
     *   HEEL2SOLE : Transition of foot angle (-1 * heel_angle -> 0).
     *   SOLE2 : Sole contact (end). Foot angle is 0.
     */
    enum toe_heel_phase {SOLE0, SOLE2TOE, TOE2SOLE, SOLE1, SOLE2HEEL, HEEL2SOLE, SOLE2, NUM_TH_PHASES};

    /* Manager for toe heel phase. */
    class toe_heel_phase_counter
    {
        double toe_heel_phase_ratio[NUM_TH_PHASES];
        size_t toe_heel_phase_count[NUM_TH_PHASES], one_step_count;
        bool calc_toe_heel_phase_count_from_raio ()
        {
            double ratio_sum = 0.0;
            for (size_t i = 0; i < NUM_TH_PHASES; i++) {
                ratio_sum += toe_heel_phase_ratio[i];
                toe_heel_phase_count[i] = static_cast<size_t>(one_step_count * ratio_sum);
            }
        };
    public:
        toe_heel_phase_counter () : one_step_count(0)
        {
            toe_heel_phase_ratio[SOLE0] = 0.05;
            toe_heel_phase_ratio[SOLE2TOE] = 0.25;
            toe_heel_phase_ratio[TOE2SOLE] = 0.2;
            toe_heel_phase_ratio[SOLE1] = 0.0;
            toe_heel_phase_ratio[SOLE2HEEL] = 0.2;
            toe_heel_phase_ratio[HEEL2SOLE] = 0.25;
            toe_heel_phase_ratio[SOLE2] = 0.05;
        };
        bool check_toe_heel_phase_ratio_validity (const std::vector<double>& ratio)
        {
            bool ret = true;
            // Check size
            if (ratio.size() != NUM_TH_PHASES) {
                ret = false;
            }
            // Check sum == 1.0
            double sum_ratio = 0.0;
            for (int i = 0; i < NUM_TH_PHASES; i++) sum_ratio += ratio[i];
            if (std::fabs(sum_ratio-1.0) > 1e-3) {
                ret = false;
            }
            if (!ret) {
                std::cerr << "toe_heel_phase_ratio is not set, "
                          << ", required length = " << NUM_TH_PHASES << " != input length " << ratio.size()
                          << ", sum_ratio = " << sum_ratio << " is not 1.0."
                          << std::endl;
            } else {
                std::cerr << "toe_heel_phase_ratio is successfully set." << std::endl;
            }
            return ret;
        };
        // setter
        void set_one_step_count (const size_t _count)
        {
            one_step_count = _count;
            calc_toe_heel_phase_count_from_raio();
        };
        bool set_toe_heel_phase_ratio (const std::vector<double>& ratio)
        {
            if (check_toe_heel_phase_ratio_validity(ratio)) {
                for (size_t i = 0; i < NUM_TH_PHASES; i++) toe_heel_phase_ratio[i] = ratio[i];
                return true;
            } else {
                return false;
            }
        };
        // getter
        void get_toe_heel_phase_ratio (std::vector<double>& ratio) const
        {
            for (size_t i = 0; i < NUM_TH_PHASES; i++) ratio[i] = toe_heel_phase_ratio[i];
        };
        int get_NUM_TH_PHASES () const { return NUM_TH_PHASES; };
        // functions for checking phase and calculating phase time and ratio
        bool is_phase_starting (const size_t current_count, const toe_heel_phase _phase) const
        {
            return (current_count == toe_heel_phase_count[_phase]);
        };
        bool is_between_phases (const size_t current_count, const toe_heel_phase phase0, const toe_heel_phase phase1) const
        {
            return (toe_heel_phase_count[phase0] <= current_count) && (current_count < toe_heel_phase_count[phase1]);
        };
        bool is_between_phases (const size_t current_count, const toe_heel_phase phase1) const
        {
            return (current_count < toe_heel_phase_count[phase1]);
        };
        bool is_no_SOLE1_phase () const { return toe_heel_phase_count[TOE2SOLE] == toe_heel_phase_count[SOLE1]; };
        double calc_phase_period (const toe_heel_phase start_phase, const toe_heel_phase goal_phase, const double _dt) const
        {
            return _dt * (toe_heel_phase_count[goal_phase]-toe_heel_phase_count[start_phase]);
        };
        double calc_phase_ratio (const size_t current_count, const toe_heel_phase start_phase, const toe_heel_phase goal_phase) const
        {
            return static_cast<double>(current_count-toe_heel_phase_count[start_phase]) / (toe_heel_phase_count[goal_phase]-toe_heel_phase_count[start_phase]);
        };
        double calc_phase_ratio (const size_t current_count, const toe_heel_phase goal_phase) const
        {
            return static_cast<double>(current_count) / (toe_heel_phase_count[goal_phase]);
        };
    };

    /* refzmp_generator to generate current refzmp from footstep_node_list */
    class refzmp_generator
    {
#ifdef HAVE_MAIN
    public:
#endif
      std::vector<hrp::Vector3> refzmp_cur_list;
      std::map<leg_type, double> zmp_weight_map;
      std::vector< std::vector<hrp::Vector3> > foot_x_axises_list; // Swing foot x axis list according to refzmp_cur_list
      std::vector< std::vector<leg_type> > swing_leg_types_list; // Swing leg list according to refzmp_cur_list
      std::vector<size_t> step_count_list; // Swing leg list according to refzmp_cur_list
      std::vector<hrp::Vector3> default_zmp_offsets; /* list of RLEG and LLEG */
      size_t refzmp_index, refzmp_count, one_step_count;
      double toe_zmp_offset_x, heel_zmp_offset_x; // [m]
      double dt;
      toe_heel_phase_counter* thp_ptr;
      bool use_toe_heel_transition;
      boost::shared_ptr<interpolator> zmp_weight_interpolator;
      std::vector<size_t> same_footstep_index_list;
      void calc_current_refzmp (hrp::Vector3& ret, std::vector<hrp::Vector3>& swing_foot_zmp_offsets, const double default_double_support_ratio_before, const double default_double_support_ratio_after, const double default_double_support_static_ratio_before, const double default_double_support_static_ratio_after);
      const bool is_start_double_support_phase () const { return refzmp_index == 0; };
      const bool is_second_phase () const { return refzmp_index == 1; };
      const bool is_second_last_phase () const { return refzmp_index == refzmp_cur_list.size()-2; };
      const bool is_end_double_support_phase () const { return refzmp_index == refzmp_cur_list.size() - 1; };
#ifndef HAVE_MAIN
    public:
#endif
      refzmp_generator(toe_heel_phase_counter* _thp_ptr, const double _dt)
        : refzmp_cur_list(), foot_x_axises_list(), swing_leg_types_list(), step_count_list(), default_zmp_offsets(),
          refzmp_index(0), refzmp_count(0), one_step_count(0),
          toe_zmp_offset_x(0), heel_zmp_offset_x(0), dt(_dt),
          thp_ptr(_thp_ptr), use_toe_heel_transition(false), same_footstep_index_list()
      {
          default_zmp_offsets.push_back(hrp::Vector3::Zero());
          default_zmp_offsets.push_back(hrp::Vector3::Zero());
          default_zmp_offsets.push_back(hrp::Vector3::Zero());
          default_zmp_offsets.push_back(hrp::Vector3::Zero());
          double zmp_weight_initial_value[4] = {1.0, 1.0, 0.1, 0.1};
          zmp_weight_map = boost::assign::map_list_of<leg_type, double>(RLEG, zmp_weight_initial_value[0])(LLEG, zmp_weight_initial_value[1])(RARM, zmp_weight_initial_value[2])(LARM, zmp_weight_initial_value[3]);
          zmp_weight_interpolator = boost::shared_ptr<interpolator>(new interpolator(4, dt));
          zmp_weight_interpolator->set(zmp_weight_initial_value); /* set initial value */
          zmp_weight_interpolator->setName("GaitGenerator zmp_weight_interpolator");
      };
      ~refzmp_generator()
      {
          thp_ptr = NULL;
      };
      void remove_refzmp_cur_list_over_length (const size_t len)
      {
        while ( refzmp_cur_list.size() > len) refzmp_cur_list.pop_back();
        while ( foot_x_axises_list.size() > len) foot_x_axises_list.pop_back();
        while ( swing_leg_types_list.size() > len) swing_leg_types_list.pop_back();
        while ( step_count_list.size() > len) step_count_list.pop_back();
      };
      void reset (const size_t _refzmp_count)
      {
        set_indices(0);
        one_step_count = _refzmp_count;
        set_refzmp_count(_refzmp_count);
        refzmp_cur_list.clear();
        foot_x_axises_list.clear();
        swing_leg_types_list.clear();
        step_count_list.clear();
        same_footstep_index_list.clear();
      };
      void push_refzmp_from_footstep_nodes_for_dual (const std::vector<step_node>& fns,
                                                     const std::vector<step_node>& _support_leg_steps,
                                                     const std::vector<step_node>& _swing_leg_steps);
      void push_refzmp_from_footstep_nodes_for_single (const std::vector<step_node>& fns, const std::vector<step_node>& _support_leg_stepsconst);
      void update_refzmp (const std::vector< std::vector<step_node> >& fnsl);
      void push_same_footstep_index (const size_t _same_footstep_index){
        same_footstep_index_list.push_back(_same_footstep_index);
      }
      // setter
      void set_indices (const size_t idx) { refzmp_index = idx; };
      void set_refzmp_count(const size_t _refzmp_count) { refzmp_count = _refzmp_count; };
      void set_default_zmp_offsets(const std::vector<hrp::Vector3>& tmp) { default_zmp_offsets = tmp; };
      void set_toe_zmp_offset_x (const double _off) { toe_zmp_offset_x = _off; };
      void set_heel_zmp_offset_x (const double _off) { heel_zmp_offset_x = _off; };
      void set_use_toe_heel_transition (const double _u) { use_toe_heel_transition = _u; };
      void set_zmp_weight_map (const std::map<leg_type, double> _map) {
          double zmp_weight_array[4] = {_map.find(RLEG)->second, _map.find(LLEG)->second, _map.find(RARM)->second, _map.find(LARM)->second};
          if (zmp_weight_interpolator->isEmpty()) {
              zmp_weight_interpolator->clear();
              double zmp_weight_initial_value[4] = {zmp_weight_map[RLEG], zmp_weight_map[LLEG], zmp_weight_map[RARM], zmp_weight_map[LARM]};
              zmp_weight_interpolator->set(zmp_weight_initial_value);
              zmp_weight_interpolator->go(zmp_weight_array, 2.0, true);
          } else {
              std::cerr << "zmp_weight_map cannot be set because interpolating." << std::endl;
          }
      };
      // getter
      bool get_current_refzmp (hrp::Vector3& rzmp, std::vector<hrp::Vector3>& swing_foot_zmp_offsets, const double default_double_support_ratio_before, const double default_double_support_ratio_after, const double default_double_support_static_ratio_before, const double default_double_support_static_ratio_after)
      {
        if (refzmp_cur_list.size() > refzmp_index ) calc_current_refzmp(rzmp, swing_foot_zmp_offsets, default_double_support_ratio_before, default_double_support_ratio_after, default_double_support_static_ratio_before, default_double_support_static_ratio_after);
        return refzmp_cur_list.size() > refzmp_index;
      };
      const hrp::Vector3& get_refzmp_cur () const { return refzmp_cur_list.front(); };
      const hrp::Vector3& get_default_zmp_offset (const leg_type lt) const { return default_zmp_offsets[lt]; };
      double get_toe_zmp_offset_x () const { return toe_zmp_offset_x; };
      double get_heel_zmp_offset_x () const { return heel_zmp_offset_x; };
      bool get_use_toe_heel_transition () const { return use_toe_heel_transition; };
      const std::map<leg_type, double> get_zmp_weight_map () const { return zmp_weight_map; };
      void proc_zmp_weight_map_interpolation () {
          if (!zmp_weight_interpolator->isEmpty()) {
              double zmp_weight_output[4];
              zmp_weight_interpolator->get(zmp_weight_output, true);
              zmp_weight_map = boost::assign::map_list_of<leg_type, double>(RLEG, zmp_weight_output[0])(LLEG, zmp_weight_output[1])(RARM, zmp_weight_output[2])(LARM, zmp_weight_output[3]);
          }
      };
    };


    class delay_hoffarbib_trajectory_generator
    {
    private:
      hrp::Vector3 pos, vel, acc; // [m], [m/s], [m/s^2]
      hrp::Vector3 skate_acc;
      double dt; // [s]
      // Implement hoffarbib to configure remain_time;
      void hoffarbib_interpolation (const double tmp_remain_time, const hrp::Vector3& tmp_goal)
      {
        hrp::Vector3 jerk = (-9.0/ tmp_remain_time) * acc +
          (-36.0 / (tmp_remain_time * tmp_remain_time)) * vel +
          (60.0 / (tmp_remain_time * tmp_remain_time * tmp_remain_time)) * (tmp_goal - pos);
        acc = acc + dt * jerk;
        vel = vel + dt * acc;
        pos = pos + dt * vel;
      };
      void hoffarbib_interpolation_all (const double tmp_remain_time, const hrp::Vector3& tmp_goal, const hrp::Vector3& tmp_goal_vel, const hrp::Vector3& tmp_goal_acc)
      {
        hrp::Vector3 jerk = (-9.0/ tmp_remain_time) * (acc - tmp_goal_acc / 3.0) +
          (-36.0 / (tmp_remain_time * tmp_remain_time)) * (tmp_goal_vel * 2.0 / 3.0 + vel) +
          (60.0 / (tmp_remain_time * tmp_remain_time * tmp_remain_time)) * (tmp_goal - pos);
        acc = acc + dt * jerk;
        vel = vel + dt * acc;
        pos = pos + dt * vel;
      };
    protected:
      double time_offset; // [s]
      double final_distance_weight;
      size_t one_step_count, current_count, double_support_count_before, double_support_count_after; // time/dt
      hrp::Vector3 take_off_vel; // [m/s^2]
      hrp::Vector3 landing_vel; // [m/s^2]
      virtual hrp::Vector3 interpolate_antecedent_path (const hrp::Vector3& start, const hrp::Vector3& goal, const double height, const double tmp_ratio) = 0;
    public:
      delay_hoffarbib_trajectory_generator () : time_offset(0.35), final_distance_weight(1.0), one_step_count(0), current_count(0), double_support_count_before(0), double_support_count_after(0), take_off_vel(hrp::Vector3(-0.3,0,0)), landing_vel(hrp::Vector3(0.0,0,0)) {};
      ~delay_hoffarbib_trajectory_generator() { };
      void set_dt (const double _dt) { dt = _dt; };
      void set_swing_trajectory_delay_time_offset (const double _time_offset) { time_offset = _time_offset; };
      void set_swing_trajectory_final_distance_weight (const double _final_distance_weight) { final_distance_weight = _final_distance_weight; };
      void set_swing_leg_take_off_vel(const hrp::Vector3 _take_off_vel) { take_off_vel = _take_off_vel; };
      void set_swing_leg_landing_vel(const hrp::Vector3 _landing_vel) { landing_vel = _landing_vel; };
      void reset (const size_t _one_step_len, const double default_double_support_ratio_before, const double default_double_support_ratio_after)
      {
        one_step_count = _one_step_len;
        current_count = 0;
        double_support_count_before = (default_double_support_ratio_before*one_step_count);
        double_support_count_after = (default_double_support_ratio_after*one_step_count);
      };
      void reset_all (const double _dt, const size_t _one_step_len,
                      const double default_double_support_ratio_before, const double default_double_support_ratio_after,
                      const double _time_offset, const double _final_distance_weight)
      {
          set_dt(_dt);
          reset (_one_step_len, default_double_support_ratio_before, default_double_support_ratio_after);
          set_swing_trajectory_delay_time_offset(_time_offset);
          set_swing_trajectory_final_distance_weight(_final_distance_weight);
      };
      void get_trajectory_point (hrp::Vector3& ret, const hrp::Vector3& start, const hrp::Vector3& goal, const double height)
      {
        if ( double_support_count_before <= current_count && current_count < one_step_count - double_support_count_after ) { // swing phase
          size_t swing_remain_count = one_step_count - current_count - double_support_count_after;
          size_t swing_one_step_count = one_step_count - double_support_count_before - double_support_count_after;
          if (swing_remain_count*dt > time_offset) { // antecedent path is still interpolating
            hoffarbib_interpolation (time_offset, interpolate_antecedent_path(start, goal, height, ((swing_one_step_count - swing_remain_count) / (swing_one_step_count - time_offset/dt))));
          } else if (swing_remain_count > 0) { // antecedent path already reached to goal
            hoffarbib_interpolation (swing_remain_count*dt, goal);
          } else {
            pos = goal;
          }
        } else if ( current_count < double_support_count_before ) { // first double support phase
          pos = start;
          vel = hrp::Vector3::Zero();
          acc = hrp::Vector3::Zero();
        } else { // last double support phase
          pos = goal;
          vel = hrp::Vector3::Zero();
          acc = hrp::Vector3::Zero();
        }
        ret = pos;
        current_count++;
      };
      void get_trajectory_point_for_skate (hrp::Vector3& ret, const hrp::Vector3& start, const hrp::Vector3& goal, const double height, const bool use_start_vel, const bool use_goal_vel)
      {
        /// tmp
        landing_vel = hrp::Vector3::Zero();
        /* landing_vel = hrp::Vector3(-0.1,0,0); */
        /// set paramter
        hrp::Vector3 take_off_dist = take_off_vel * double_support_count_before * dt / 2.0;
        hrp::Vector3 landing_dist = - landing_vel * double_support_count_after * dt / 2.0;
        hrp::Vector3 start_acc = hrp::Vector3::Zero();
        hrp::Vector3 start_vel = hrp::Vector3::Zero();
        hrp::Vector3 goal_acc = hrp::Vector3::Zero();
        hrp::Vector3 goal_vel = hrp::Vector3::Zero();
        double double_support_count_before_for_skate = double_support_count_before;
        double double_support_count_after_for_skate = double_support_count_after;
        if ( use_start_vel ) {
            double_support_count_before_for_skate *= 0.5;
            start_acc = 1.5 * (take_off_vel - landing_vel) / (double_support_count_before * dt);
            start_vel = 0.5 * take_off_vel + 7.0 * landing_vel / 16.0;
            take_off_dist -= ( 3.0 *  take_off_vel + 5.0 * landing_vel) * double_support_count_before * dt / 32.0;
        }
        if ( use_goal_vel ) {
            double_support_count_after_for_skate *= 0.5;
            goal_acc = 1.5 * (take_off_vel - landing_vel) / (double_support_count_after * dt);
            goal_vel = 0.5 * take_off_vel + 7.0 * landing_vel / 16.0;
            landing_dist = - ( 3.0 * take_off_vel + 5.0 * landing_vel ) * double_support_count_after * dt /32.0;
        }
        /// update
        if ( current_count == 0 ){
          acc = start_acc;
          vel = start_vel + acc * dt;
          pos = start + vel * dt;
          skate_acc = acc;
        } else if ( current_count == one_step_count ){
          acc = goal_acc;
          vel = goal_vel;
          pos = goal;
          skate_acc = acc;
        } else if ( double_support_count_before_for_skate <= current_count && current_count < one_step_count - double_support_count_after_for_skate ) { // swing phase
          size_t swing_remain_count = one_step_count - current_count - double_support_count_after_for_skate;
          size_t swing_one_step_count = one_step_count - double_support_count_before_for_skate - double_support_count_after_for_skate;
          if (swing_remain_count*dt > time_offset) { // antecedent path is still interpolating
            hoffarbib_interpolation_all (time_offset, interpolate_antecedent_path(start + take_off_dist, goal + landing_dist, height, ((swing_one_step_count - swing_remain_count) / (swing_one_step_count - time_offset/dt))), landing_vel, hrp::Vector3::Zero());
            skate_acc = hrp::Vector3::Zero();
          } else if (swing_remain_count > 0) { // antecedent path already reached to goal
            hoffarbib_interpolation_all (swing_remain_count*dt, goal + landing_dist, landing_vel, hrp::Vector3::Zero());
            skate_acc = hrp::Vector3::Zero();
          } else { // swing_remain_count = 0
            pos += vel * dt;
            skate_acc = hrp::Vector3::Zero();
          }
        } else if ( current_count < double_support_count_before_for_skate ) { // first double support phase
          hoffarbib_interpolation_all ((double_support_count_before_for_skate  - current_count) * dt, start + take_off_dist, take_off_vel, hrp::Vector3::Zero());
          skate_acc = acc;
        } else { // last double support phase
          hoffarbib_interpolation_all((one_step_count - current_count + 1) * dt, goal, goal_vel, goal_acc);
          skate_acc = acc;
        }
        ret = pos;
        current_count++;
      };
      double get_swing_trajectory_delay_time_offset () const { return time_offset; };
      double get_swing_trajectory_final_distance_weight () const { return final_distance_weight; };
      hrp::Vector3 get_swing_leg_take_off_vel () const { return take_off_vel; };
      hrp::Vector3 get_swing_leg_landing_vel () const { return landing_vel; };
      hrp::Vector3 get_skate_acc () const { return skate_acc; };
      // interpolate path vector
      //   tmp_ratio : ratio value [0, 1]
      //   org_point_vec : vector of via points
      //   e.g., move tmp_ratio from 0 to 1 => move point from org_point_vec.front() to org_point_vec.back()
      hrp::Vector3 interpolate_antecedent_path_base (const double tmp_ratio, const std::vector<hrp::Vector3> org_point_vec)
      {
        std::vector<hrp::Vector3> point_vec;
        std::vector<double> distance_vec;
        double total_path_length = 0;
        point_vec.push_back(org_point_vec.front());
        // remove distance-zero points
        for (size_t i = 0; i < org_point_vec.size()-1; i++) {
          double tmp_distance = (org_point_vec[i+1]-org_point_vec[i]).norm();
          if (i==org_point_vec.size()-2) tmp_distance*=final_distance_weight;
          if ( tmp_distance > 1e-5 ) {
            point_vec.push_back(org_point_vec[i+1]);
            distance_vec.push_back(tmp_distance);
            total_path_length += tmp_distance;
          }
        }
        if ( total_path_length < 1e-5 ) { // if total path is zero, return goal point.
          return org_point_vec.back();
        }
        // point_vec        : [p0, p1, ..., pN-1, pN]
        // distance_vec     : [  d0, ...,     dN-1  ]
        // sum_distance_vec : [l0, l1, ..., lN-1, lN] <= lj = \Sum_{i=0}^{j-1} di
        std::vector<double> sum_distance_vec;
        sum_distance_vec.push_back(0);
        double tmp_dist = 0;
        for (size_t i = 0; i < distance_vec.size(); i++) {
          sum_distance_vec.push_back(tmp_dist + distance_vec[i]);
          tmp_dist += distance_vec[i];
        }
        // select current segment in which 'tmp_ratio' is included
        double current_length = tmp_ratio * total_path_length;
        for (size_t i = 0; i < sum_distance_vec.size(); i++) {
          if ( (sum_distance_vec[i] <= current_length) && (current_length <= sum_distance_vec[i+1]) ) {
            double tmpr = ((current_length - sum_distance_vec[i]) / distance_vec[i]);
            return ((1-tmpr) * point_vec[i] + tmpr * point_vec[1+i]);
          }
        }
        // if illegal tmp-ratio
        if (current_length < 0) return org_point_vec.front();
        else org_point_vec.back();
      };
    };

    class rectangle_delay_hoffarbib_trajectory_generator : public delay_hoffarbib_trajectory_generator
    {
      hrp::Vector3 interpolate_antecedent_path (const hrp::Vector3& start, const hrp::Vector3& goal, const double height, const double tmp_ratio)
      {
        std::vector<hrp::Vector3> rectangle_path;
        double max_height = std::max(start(2), goal(2))+height;
        rectangle_path.push_back(start);
        rectangle_path.push_back(hrp::Vector3(start(0), start(1), max_height));
        rectangle_path.push_back(hrp::Vector3(goal(0), goal(1), max_height));
        rectangle_path.push_back(goal);
        return interpolate_antecedent_path_base(tmp_ratio, rectangle_path);
      };
    };

    class stair_delay_hoffarbib_trajectory_generator : public delay_hoffarbib_trajectory_generator
    {
      hrp::Vector3 way_point_offset;
      hrp::Vector3 interpolate_antecedent_path (const hrp::Vector3& start, const hrp::Vector3& goal, const double height, const double tmp_ratio)
      {
        std::vector<hrp::Vector3> path;
        double max_height = std::max(start(2), goal(2))+height;
        hrp::Vector3 diff_vec = goal - start;
        diff_vec(2) = 0.0; // projection on horizontal plane
        path.push_back(start);
        // currently way_point_offset(1) is not used.
        //if ( diff_vec.norm() > 1e-4 && (goal(2) - start(2)) > way_point_offset(2) ) {
        if ( diff_vec.norm() > 1e-4 && (goal(2) - start(2)) > 0.02) {
          path.push_back(hrp::Vector3(start+-1*way_point_offset(0)*diff_vec.normalized()+hrp::Vector3(0,0,way_point_offset(2)+max_height-start(2))));
        }
        path.push_back(hrp::Vector3(start(0), start(1), max_height));
        path.push_back(hrp::Vector3(goal(0), goal(1), max_height));
        //if ( diff_vec.norm() > 1e-4 && (start(2) - goal(2)) > way_point_offset(2) ) {
        if ( diff_vec.norm() > 1e-4 && (start(2) - goal(2)) > 0.02) {
          path.push_back(hrp::Vector3(goal+way_point_offset(0)*diff_vec.normalized()+hrp::Vector3(0,0,way_point_offset(2)+max_height-goal(2))));
        }
        // if (height > 20 * 1e-3) {
        //   path.push_back(hrp::Vector3(goal(0), goal(1), 20*1e-3+goal(2)));
        // }
        path.push_back(goal);
        return interpolate_antecedent_path_base(tmp_ratio, path);
      };
    public:
      stair_delay_hoffarbib_trajectory_generator () : delay_hoffarbib_trajectory_generator(), way_point_offset(hrp::Vector3(0.03, 0.0, 0.0)) {};
      ~stair_delay_hoffarbib_trajectory_generator () {};
      void set_stair_trajectory_way_point_offset (const hrp::Vector3 _offset) { way_point_offset = _offset; };
      hrp::Vector3 get_stair_trajectory_way_point_offset() const { return way_point_offset; };
    };

    class cycloid_delay_hoffarbib_trajectory_generator : public delay_hoffarbib_trajectory_generator
    {
      hrp::Vector3 interpolate_antecedent_path (const hrp::Vector3& start, const hrp::Vector3& goal, const double height, const double tmp_ratio)
      {
        std::vector<hrp::Vector3> cycloid_path;
        hrp::Vector3 tmpv, via_goal(goal);
        double ratio = 0.4;
        via_goal(2) += ratio*height;
        double tmpheight = ((start(2)+goal(2))/2.0+height-(start(2)+via_goal(2))/2.0);
        cycloid_path.push_back(start);
        cycloid_midpoint(tmpv, 0.2, start, via_goal, tmpheight);
        cycloid_path.push_back(tmpv);
        cycloid_midpoint(tmpv, 0.4, start, via_goal, tmpheight);
        cycloid_path.push_back(tmpv);
        cycloid_midpoint(tmpv, 0.6, start, via_goal, tmpheight);
        cycloid_path.push_back(tmpv);
        cycloid_midpoint(tmpv, 0.8, start, via_goal, tmpheight);
        cycloid_path.push_back(tmpv);
        cycloid_path.push_back(via_goal);
        cycloid_path.push_back(goal);
        return interpolate_antecedent_path_base(tmp_ratio, cycloid_path);
      };
    };

    class cycloid_delay_kick_hoffarbib_trajectory_generator : public delay_hoffarbib_trajectory_generator
    {
    private:
      hrp::Matrix33 start_rot;
      hrp::Vector3 kick_point_offset;
    public:
      cycloid_delay_kick_hoffarbib_trajectory_generator() : delay_hoffarbib_trajectory_generator(), kick_point_offset(hrp::Vector3(-0.01, 0.0, 0.03)) {};
      void set_cycloid_delay_kick_point_offset (const hrp::Vector3 _offset) { kick_point_offset = _offset; };
      void set_start_rot (const hrp::Matrix33 _offset) { start_rot = _offset; };
      hrp::Vector3 get_cycloid_delay_kick_point_offset () const { return kick_point_offset; };
      hrp::Vector3 interpolate_antecedent_path (const hrp::Vector3& start, const hrp::Vector3& goal, const double height, const double tmp_ratio)
      {
        std::vector<hrp::Vector3> cycloid_path;
        hrp::Vector3 tmpv, via_goal(goal);
        double ratio = 0.4;
        via_goal(2) += ratio*height;
        double tmpheight = ((start(2)+goal(2))/2.0+height-(start(2)+via_goal(2))/2.0);
        /* cycloid_path.push_back(start); */
        if(height > 1e-4){
            cycloid_path.push_back(start + start_rot * kick_point_offset);
            cycloid_midpoint(tmpv, 0.2, start + start_rot * kick_point_offset, via_goal, tmpheight);
            cycloid_path.push_back(tmpv);
            cycloid_midpoint(tmpv, 0.4, start + start_rot * kick_point_offset, via_goal, tmpheight);
            cycloid_path.push_back(tmpv);
            cycloid_midpoint(tmpv, 0.6, start + start_rot * kick_point_offset, via_goal, tmpheight);
            cycloid_path.push_back(tmpv);
            cycloid_midpoint(tmpv, 0.8, start + start_rot * kick_point_offset, via_goal, tmpheight);
            cycloid_path.push_back(tmpv);
        }
        cycloid_path.push_back(via_goal);
        cycloid_path.push_back(goal);
        return interpolate_antecedent_path_base(tmp_ratio, cycloid_path);
      };
    };

    class cross_delay_hoffarbib_trajectory_generator : public delay_hoffarbib_trajectory_generator
    {
    private:
      leg_type swing_leg;
    public:
      cross_delay_hoffarbib_trajectory_generator () : delay_hoffarbib_trajectory_generator(), way_point_offset(hrp::Vector3(0.04, 0.15, 0.0)) {};
      ~cross_delay_hoffarbib_trajectory_generator () {};
      void set_swing_leg (leg_type _lr) { swing_leg = _lr; };
      hrp::Vector3 way_point_offset;
      hrp::Vector3 interpolate_antecedent_path (const hrp::Vector3& start, const hrp::Vector3& goal, const double height, const double tmp_ratio)
      {
        std::vector<hrp::Vector3> path;
        double max_height = std::max(start(2), goal(2))+height;
        hrp::Vector3 diff_vec = goal - start;
        diff_vec(2) = 0.0; // projection on horizontal plane
        path.push_back(start);
        if ( swing_leg == LLEG ) { // swing_leg is left
            path.push_back(hrp::Vector3(start+-1*way_point_offset(0)*diff_vec.normalized()+hrp::Vector3(0,way_point_offset(1),way_point_offset(2)+max_height-start(2))));
            path.push_back(hrp::Vector3(goal+way_point_offset(0)*diff_vec.normalized()+hrp::Vector3(0,way_point_offset(1),way_point_offset(2)+max_height-goal(2))));
        } else { // swing_leg is right
            path.push_back(hrp::Vector3(start+-1*way_point_offset(0)*diff_vec.normalized()+hrp::Vector3(0,-way_point_offset(1),way_point_offset(2)+max_height-start(2))));
            path.push_back(hrp::Vector3(goal+way_point_offset(0)*diff_vec.normalized()+hrp::Vector3(0,-way_point_offset(1),way_point_offset(2)+max_height-goal(2))));
        }
        if (height > 30 * 1e-3) {
          path.push_back(hrp::Vector3(goal(0), goal(1), 30*1e-3+goal(2)));
        }
        path.push_back(goal);
        return interpolate_antecedent_path_base(tmp_ratio, path);
      };
    };

    /* leg_coords_generator to generate current swing_leg_coords and support_leg_coords from footstep_node_list */
    class leg_coords_generator
    {
#ifdef HAVE_MAIN
    public:
#endif
      std::vector< std::vector<step_node> > swing_leg_dst_steps_list, support_leg_steps_list;
      // Support leg coordinates.
      std::vector<step_node> support_leg_steps;
      // Swing leg coordinates is interpolated from swing_leg_src_coords to swing_leg_dst_coords during swing phase.
      std::vector<step_node> swing_leg_steps, swing_leg_src_steps, swing_leg_dst_steps;
      double default_step_height, default_top_ratio, current_step_height, swing_ratio, swing_rot_ratio, foot_midcoords_ratio, dt, current_toe_angle, current_heel_angle;
      double time_offset, final_distance_weight;
      std::vector<double> current_swing_time;
      // Index for current footstep. footstep_index should be [0,footstep_node_list.size()]. Current footstep is footstep_node_list[footstep_index].
      size_t footstep_index;
      // one_step_count is total counter num of current steps (= step_time/dt). lcg_count is counter for lcg. During one step, lcg_count decreases from one_step_count to 0.
      size_t lcg_count, one_step_count, next_one_step_count;
      // Current support leg
      std::vector<leg_type> support_leg_types, swing_leg_types;
      orbit_type default_orbit_type;
      bool is_swing_phase;
      // Foot trajectory generators
      std::vector<rectangle_delay_hoffarbib_trajectory_generator> rdtg;
      stair_delay_hoffarbib_trajectory_generator sdtg;
      std::vector<cycloid_delay_hoffarbib_trajectory_generator> cdtg;
      cycloid_delay_kick_hoffarbib_trajectory_generator cdktg;
      cross_delay_hoffarbib_trajectory_generator crdtg;
      toe_heel_phase_counter* thp_ptr;
      interpolator* foot_ratio_interpolator;
      interpolator* swing_foot_rot_ratio_interpolator;
      // Parameters for toe-heel contact
      interpolator* toe_heel_interpolator;
      double toe_pos_offset_x, heel_pos_offset_x, toe_angle, heel_angle, foot_dif_rot_angle;
      bool use_toe_joint;
      void calc_current_swing_leg_steps (std::vector<step_node>& rets, const double step_height, const double _current_toe_angle, const double _current_heel_angle, const double footstep_num);
      double calc_interpolated_toe_heel_angle (const toe_heel_phase start_phase, const toe_heel_phase goal_phase, const double start, const double goal);
      void modif_foot_coords_for_toe_heel_phase (coordinates& org_coords, const double _current_toe_angle, const double _current_heel_angle);
      void cycloid_midcoords (coordinates& ret, const coordinates& start,
                              const coordinates& goal, const double height) const;
      void rectangle_midcoords (coordinates& ret, const coordinates& start,
                                const coordinates& goal, const double height, const size_t swing_trajectory_generator_idx);
      void stair_midcoords (coordinates& ret, const coordinates& start,
                            const coordinates& goal, const double height);
      void cycloid_delay_midcoords (coordinates& ret, const coordinates& start,
                                    const coordinates& goal, const double height, const size_t swing_trajectory_generator_idx);
      void cycloid_delay_kick_midcoords (coordinates& ret, const coordinates& start,
                                         const coordinates& goal, const double height, const double footstep_num);
      void cross_delay_midcoords (coordinates& ret, const coordinates& start,
                                  const coordinates& goal, const double height, leg_type lr);
      void calc_ratio_from_double_support_ratio (const double default_double_support_ratio_before, const double default_double_support_ratio_after);
#ifndef HAVE_MAIN
    public:
#endif
      leg_coords_generator(const double _dt, toe_heel_phase_counter* _thp_ptr)
        : support_leg_steps(), swing_leg_steps(), swing_leg_src_steps(), swing_leg_dst_steps(),
          default_step_height(0.05), default_top_ratio(0.5), current_step_height(0.0), swing_ratio(0), swing_rot_ratio(0), foot_midcoords_ratio(0), dt(_dt),
          current_toe_angle(0), current_heel_angle(0),
          time_offset(0.35), final_distance_weight(1.0),
          footstep_index(0), lcg_count(0), default_orbit_type(CYCLOID),
          rdtg(), cdtg(),
          thp_ptr(_thp_ptr),
          foot_ratio_interpolator(NULL), swing_foot_rot_ratio_interpolator(NULL), toe_heel_interpolator(NULL),
          toe_pos_offset_x(0.0), heel_pos_offset_x(0.0), toe_angle(0.0), heel_angle(0.0), foot_dif_rot_angle(0.0), use_toe_joint(false)
      {
        support_leg_types = boost::assign::list_of<leg_type>(RLEG);
        swing_leg_types = boost::assign::list_of<leg_type>(LLEG);
        current_swing_time = boost::assign::list_of<double>(0.0)(0.0)(0.0)(0.0);
        sdtg.set_dt(dt);
        cdktg.set_dt(dt);
        crdtg.set_dt(dt);
        if (foot_ratio_interpolator == NULL) foot_ratio_interpolator = new interpolator(1, dt);
        if (swing_foot_rot_ratio_interpolator == NULL) swing_foot_rot_ratio_interpolator = new interpolator(1, dt);
        //if (foot_ratio_interpolator == NULL) foot_ratio_interpolator = new interpolator(1, dt, interpolator::LINEAR);
        if (toe_heel_interpolator == NULL) toe_heel_interpolator = new interpolator(1, dt);
        foot_ratio_interpolator->setName("GaitGenerator foot_ratio_interpolator");
        swing_foot_rot_ratio_interpolator->setName("GaitGenerator swing_foot_rot_ratio_interpolator");
        toe_heel_interpolator->setName("GaitGenerator toe_heel_interpolator");
      };
      ~leg_coords_generator()
      {
        if (foot_ratio_interpolator != NULL) {
            delete foot_ratio_interpolator;
            foot_ratio_interpolator = NULL;
        }
        if (swing_foot_rot_ratio_interpolator != NULL) {
            delete swing_foot_rot_ratio_interpolator;
            swing_foot_rot_ratio_interpolator = NULL;
        }
        if (toe_heel_interpolator != NULL) {
            delete toe_heel_interpolator;
            toe_heel_interpolator = NULL;
        }
        thp_ptr = NULL;
      };
      void set_default_step_height (const double _tmp) { default_step_height = _tmp; };
      void set_default_top_ratio (const double _tmp) { default_top_ratio = _tmp; };
      void set_default_orbit_type (const orbit_type _tmp) { default_orbit_type = _tmp; };
      void set_swing_trajectory_delay_time_offset (const double _time_offset)
      {
        sdtg.set_swing_trajectory_delay_time_offset(_time_offset);
        cdktg.set_swing_trajectory_delay_time_offset(_time_offset);
        crdtg.set_swing_trajectory_delay_time_offset(_time_offset);
        time_offset = _time_offset;
      };
      void set_swing_trajectory_final_distance_weight (const double _final_distance_weight)
      {
        sdtg.set_swing_trajectory_final_distance_weight(_final_distance_weight);
        cdktg.set_swing_trajectory_final_distance_weight(_final_distance_weight);
        crdtg.set_swing_trajectory_final_distance_weight(_final_distance_weight);
        final_distance_weight = _final_distance_weight;
      };
      void set_swing_leg_take_off_vel (const hrp::Vector3 _take_off_vel) { cdktg.set_swing_leg_take_off_vel(_take_off_vel); };
      void set_swing_leg_landing_vel (const hrp::Vector3 _landing_vel) { cdktg.set_swing_leg_landing_vel(_landing_vel); };
      void set_stair_trajectory_way_point_offset (const hrp::Vector3 _offset) { sdtg.set_stair_trajectory_way_point_offset(_offset); };
      void set_cycloid_delay_kick_point_offset (const hrp::Vector3 _offset) { cdktg.set_cycloid_delay_kick_point_offset(_offset); };
      void set_toe_pos_offset_x (const double _offx) { toe_pos_offset_x = _offx; };
      void set_heel_pos_offset_x (const double _offx) { heel_pos_offset_x = _offx; };
      void set_toe_angle (const double _angle) { toe_angle = _angle; };
      void set_heel_angle (const double _angle) { heel_angle = _angle; };
      void set_use_toe_joint (const bool ut) { use_toe_joint = ut; };
      void set_swing_support_steps_list (const std::vector< std::vector<step_node> >& fnsl)
      {
          std::vector<step_node> prev_support_leg_steps = support_leg_steps_list.front();
          support_leg_steps_list.clear();
          swing_leg_dst_steps_list.clear();
          support_leg_steps_list.push_back(prev_support_leg_steps);
          swing_leg_dst_steps_list = fnsl;
          for (size_t i = 0; i < fnsl.size(); i++) {
              if (i > 0) {
                  if (is_same_footstep_nodes(fnsl.at(i), fnsl.at(i-1))) {
                      support_leg_steps_list.push_back(support_leg_steps_list.back());
                  } else {
                      /* current support leg steps = prev swing leg dst steps + (prev support leg steps without current swing leg names) */
                      std::vector<step_node> tmp_support_leg_steps = swing_leg_dst_steps_list.at(i-1);
                      std::copy(support_leg_steps_list.back().begin(),
                                support_leg_steps_list.back().end(),
                                std::back_inserter(tmp_support_leg_steps));
                      for (size_t j = 0; j < swing_leg_dst_steps_list.at(i).size(); j++) {
                          std::vector<step_node>::iterator it = std::remove_if(tmp_support_leg_steps.begin(),
                                                                               tmp_support_leg_steps.end(),
                                                                               (&boost::lambda::_1->* &step_node::l_r == swing_leg_dst_steps_list.at(i).at(j).l_r));
                          tmp_support_leg_steps.erase(it, tmp_support_leg_steps.end());
                      }
                      support_leg_steps_list.push_back(tmp_support_leg_steps);
              }
            }
          }
      };
      void reset(const size_t _one_step_count, const size_t _next_one_step_count,
                 const std::vector<step_node>& _swing_leg_dst_steps,
                 const std::vector<step_node>& _swing_leg_src_steps,
                 const std::vector<step_node>& _support_leg_steps,
                 const double default_double_support_ratio_before,
                 const double default_double_support_ratio_after)
      {
        support_leg_steps_list.clear();
        swing_leg_dst_steps_list.clear();
        /* swing_leg_steps.clear(); */
        swing_leg_dst_steps = _swing_leg_dst_steps;
        swing_leg_src_steps = _swing_leg_src_steps;
        support_leg_steps = _support_leg_steps;
        support_leg_steps_list.push_back(support_leg_steps);
        one_step_count = lcg_count = _one_step_count;
        next_one_step_count = _next_one_step_count;
        thp_ptr->set_one_step_count(one_step_count);
        footstep_index = 0;
        current_step_height = 0.0;
        switch (default_orbit_type) {
        case RECTANGLE:
            rdtg.clear();
            for (size_t i = 0; i < swing_leg_dst_steps.size(); i++) {
                rdtg.push_back(rectangle_delay_hoffarbib_trajectory_generator());
                rdtg.back().reset_all(dt, one_step_count,
                                      default_double_support_ratio_before, default_double_support_ratio_after,
                                      time_offset, final_distance_weight);
            }
            break;
        case STAIR:
            sdtg.reset(one_step_count, default_double_support_ratio_before, default_double_support_ratio_after);
            break;
        case CYCLOIDDELAY:
            cdtg.clear();
            for (size_t i = 0; i < swing_leg_dst_steps.size(); i++) {
                cdtg.push_back(cycloid_delay_hoffarbib_trajectory_generator());
                cdtg.back().reset_all(dt, one_step_count,
                                      default_double_support_ratio_before, default_double_support_ratio_after,
                                      time_offset, final_distance_weight);
            }
            break;
        case CYCLOIDDELAYKICK:
            cdktg.reset(one_step_count, default_double_support_ratio_before, default_double_support_ratio_after);
            break;
        case CROSS:
            crdtg.reset(one_step_count, default_double_support_ratio_before, default_double_support_ratio_after);
            break;
        default:
            break;
        }
        reset_foot_ratio_interpolator();
      };
      void reset_foot_ratio_interpolator ()
      {
        double tmp_ratio = 0.0;
        foot_ratio_interpolator->clear();
        foot_ratio_interpolator->set(&tmp_ratio);
        tmp_ratio = 1.0;
        //foot_ratio_interpolator->go(&tmp_ratio, dt*one_step_count, true);
        foot_ratio_interpolator->setGoal(&tmp_ratio, dt*one_step_count, true);
        foot_ratio_interpolator->sync();
      };
      void clear_interpolators ( ) {
        double tmp;
        while (!swing_foot_rot_ratio_interpolator->isEmpty()) {
            swing_foot_rot_ratio_interpolator->get(&tmp, true);
        }
        while (!foot_ratio_interpolator->isEmpty()) {
            foot_ratio_interpolator->get(&tmp, true);
        }
        while (!toe_heel_interpolator->isEmpty()) {
            toe_heel_interpolator->get(&tmp, true);
        }
      };
      bool is_same_footstep_nodes(const std::vector<step_node>& fns_1, const std::vector<step_node>& fns_2) const;
      void update_leg_steps (const std::vector< std::vector<step_node> >& fnsl, const double default_double_support_ratio_before, const double default_double_support_ratio_after);
      size_t get_footstep_index() const { return footstep_index; };
      size_t get_lcg_count() const { return lcg_count; };
      size_t get_one_step_count() const { return one_step_count; };
      double get_current_swing_time(const size_t idx) const { return current_swing_time.at(idx); };
      const std::vector<step_node>& get_swing_leg_steps() const { return swing_leg_steps; };
      const std::vector<step_node>& get_support_leg_steps() const { return support_leg_steps; };
      const std::vector<step_node>& get_swing_leg_src_steps() const { return swing_leg_src_steps; };
      const std::vector<step_node>& get_swing_leg_dst_steps() const { return swing_leg_dst_steps; };
      const std::vector<step_node>& get_swing_leg_dst_steps_idx(const size_t idx) const { return swing_leg_dst_steps_list[idx]; };
      const std::vector<step_node>& get_support_leg_steps_idx(const size_t idx) const { return support_leg_steps_list[idx]; };
      const std::vector<leg_type>& get_support_leg_types() const { return support_leg_types;};
      const std::vector<leg_type>& get_swing_leg_types() const { return swing_leg_types;};
      double get_default_step_height () const { return default_step_height;};
      void get_swing_support_mid_coords(coordinates& ret) const
      {
        std::vector<coordinates> swg_coords, sup_coords;
        for (std::vector<step_node>::const_iterator it_src = swing_leg_src_steps.begin(), it_dst = swing_leg_dst_steps.begin();
             it_src != swing_leg_src_steps.end() && it_dst != swing_leg_dst_steps.end();
             it_src++, it_dst++) {
            coordinates tmp;
            mid_coords(tmp, foot_midcoords_ratio, it_src->worldcoords, it_dst->worldcoords);
            if (it_src->l_r == RLEG or it_src->l_r == LLEG) swg_coords.push_back(tmp);
        }
        for (std::vector<step_node>::const_iterator it = support_leg_steps.begin(); it != support_leg_steps.end(); it++) {
            if (it->l_r == RLEG or it->l_r == LLEG) sup_coords.push_back(it->worldcoords);
        }
        coordinates tmp_swg_mid, tmp_sup_mid;
        if (swg_coords.size() > 0) multi_mid_coords(tmp_swg_mid, swg_coords);
        if (sup_coords.size() > 0) multi_mid_coords(tmp_sup_mid, sup_coords);
        mid_coords(ret, static_cast<double>(sup_coords.size()) / (swg_coords.size() + sup_coords.size()), tmp_swg_mid, tmp_sup_mid);
      };
      std::vector<leg_type> get_current_support_states () const
      {
          if ( is_swing_phase ) {
              return get_support_leg_types();
          } else {
              std::vector<leg_type> tmp_sup_types = get_support_leg_types();
              std::vector<leg_type> tmp_swg_types = get_swing_leg_types();
              std::copy(tmp_swg_types.begin(), tmp_swg_types.end(), std::back_inserter(tmp_sup_types));
              return tmp_sup_types;
          }
      };
      orbit_type get_default_orbit_type () const { return default_orbit_type; };
      double get_swing_trajectory_delay_time_offset () const { return time_offset; };
      double get_swing_trajectory_final_distance_weight () const { return final_distance_weight; };
      hrp::Vector3 get_stair_trajectory_way_point_offset () const { return sdtg.get_stair_trajectory_way_point_offset(); };
      hrp::Vector3 get_cycloid_delay_kick_point_offset () const { return cdktg.get_cycloid_delay_kick_point_offset() ; };
      hrp::Vector3 get_swing_leg_take_off_vel () const { return cdktg.get_swing_leg_take_off_vel() ; };
      hrp::Vector3 get_swing_leg_landing_vel () const { return cdktg.get_swing_leg_landing_vel() ; };
      hrp::Vector3 get_skate_acc () const { return cdktg.get_skate_acc() ; };
      double get_toe_pos_offset_x () const { return toe_pos_offset_x; };
      double get_heel_pos_offset_x () const { return heel_pos_offset_x; };
      double get_toe_angle () const { return toe_angle; };
      double get_heel_angle () const { return heel_angle; };
      double get_foot_dif_rot_angle () const { return foot_dif_rot_angle; };
      bool get_use_toe_joint () const { return use_toe_joint; };
    };

  class gait_generator
  {

  public:
#ifndef HAVE_MAIN
  private:
#endif

    enum velocity_mode_flag { VEL_IDLING, VEL_DOING, VEL_ENDING };
    enum emergency_flag { IDLING, EMERGENCY_STOP, STOPPING };

    /* member variables for gait_generator */
    // Footstep list to be executed
    //   First and last footstep are used for double support phase.
    std::vector< std::vector<step_node> > footstep_nodes_list;
    // Footstep list for overwriting future footstep queue
    std::vector< std::vector<step_node> > overwrite_footstep_nodes_list;
    toe_heel_phase_counter thp;
    refzmp_generator rg;
    leg_coords_generator lcg;
    footstep_parameter footstep_param;
    velocity_mode_parameter vel_param, offset_vel_param;
    hrp::Vector3 cog, refzmp, prev_que_rzmp; /* cog by calculating proc_one_tick */
    std::vector<hrp::Vector3> swing_foot_zmp_offsets, prev_que_sfzos;
    double dt; /* control loop [s] */
    std::vector<std::string> all_limbs;
    double default_step_time;
    double default_double_support_ratio_before, default_double_support_ratio_after, default_double_support_static_ratio_before, default_double_support_static_ratio_after;
    double default_double_support_ratio_swing_before; /*first double support time for leg coords generator */
    double default_double_support_ratio_swing_after; /*last double support time for leg coords generator */
    double gravitational_acceleration;
    size_t finalize_count, optional_go_pos_finalize_footstep_num;
    // overwrite_footstep_index is used for footstep overwriting.
    //   When overwrite_footstep_index == get_overwritable_index(), overwrite footsteps after overwrite_footstep_index.
    size_t overwrite_footstep_index;
    // overwritable_footstep_index_offset is used for emergency stop and velocity mode.
    //   overwritable footstep index is "footstep_index + overwritable_footstep_index_offset", which is obtained by get_overwritable_index().
    size_t overwritable_footstep_index_offset;
    velocity_mode_flag velocity_mode_flg;
    emergency_flag emergency_flg;
    bool use_inside_step_limitation;
    std::map<leg_type, std::string> leg_type_map;
    coordinates initial_foot_mid_coords;
    bool solved;

    /* preview controller parameters */
    //preview_dynamics_filter<preview_control>* preview_controller_ptr;
    preview_dynamics_filter<extended_preview_control>* preview_controller_ptr;

    void append_go_pos_step_nodes (const coordinates& _ref_coords,
                                   const std::vector<leg_type>& lts)
    {
        append_go_pos_step_nodes(_ref_coords, lts, footstep_nodes_list);
    };

    void append_go_pos_step_nodes (const coordinates& _ref_coords,
                                   const std::vector<leg_type>& lts,
                                   std::vector< std::vector<step_node> >& _footstep_nodes_list) const
    {
      std::vector<step_node> sns;
      for (size_t i = 0; i < lts.size(); i++) {
          sns.push_back(step_node(lts.at(i), _ref_coords,
                                  lcg.get_default_step_height(), default_step_time,
                                  lcg.get_toe_angle(), lcg.get_heel_angle()));
          sns.at(i).worldcoords.pos += sns.at(i).worldcoords.rot * footstep_param.leg_default_translate_pos[lts.at(i)];
      }
      _footstep_nodes_list.push_back(sns);
    };
    void overwrite_refzmp_queue(const std::vector< std::vector<step_node> >& fnsl);
    void calc_ref_coords_trans_vector_velocity_mode (coordinates& ref_coords, hrp::Vector3& trans, double& dth, const std::vector<step_node>& sup_fns, const velocity_mode_parameter& cur_vel_param) const;
    void calc_next_coords_velocity_mode (std::vector< std::vector<step_node> >& ret_list, const size_t idx, const size_t future_step_num = 3);
    void append_footstep_list_velocity_mode ();
    void append_footstep_list_velocity_mode (std::vector< std::vector<step_node> >& _footstep_nodes_list, const velocity_mode_parameter& cur_vel_param) const;

#ifndef HAVE_MAIN
    /* inhibit copy constructor and copy insertion not by implementing */
    gait_generator (const gait_generator& _p);
    gait_generator &operator=(const gait_generator &_p);
  public:
#endif
    gait_generator (double _dt,
                    /* arguments for footstep_parameter */
                    const std::vector<hrp::Vector3>& _leg_pos, std::vector<std::string> _all_limbs,
                    const double _stride_fwd_x, const double _stride_y, const double _stride_theta, const double _stride_bwd_x)
        : footstep_nodes_list(), overwrite_footstep_nodes_list(), thp(), rg(&thp, _dt), lcg(_dt, &thp), all_limbs(_all_limbs),
        footstep_param(_leg_pos, _stride_fwd_x, _stride_y, _stride_theta, _stride_bwd_x),
        vel_param(), offset_vel_param(), cog(hrp::Vector3::Zero()), refzmp(hrp::Vector3::Zero()), prev_que_rzmp(hrp::Vector3::Zero()),
        dt(_dt), default_step_time(1.0), default_double_support_ratio_before(0.1), default_double_support_ratio_after(0.1), default_double_support_static_ratio_before(0.0), default_double_support_static_ratio_after(0.0), default_double_support_ratio_swing_before(0.1), default_double_support_ratio_swing_after(0.1), gravitational_acceleration(DEFAULT_GRAVITATIONAL_ACCELERATION),
        finalize_count(0), optional_go_pos_finalize_footstep_num(0), overwrite_footstep_index(0), overwritable_footstep_index_offset(1),
        velocity_mode_flg(VEL_IDLING), emergency_flg(IDLING),
        use_inside_step_limitation(true),
        preview_controller_ptr(NULL) {
        swing_foot_zmp_offsets = boost::assign::list_of<hrp::Vector3>(hrp::Vector3::Zero());
        prev_que_sfzos = boost::assign::list_of<hrp::Vector3>(hrp::Vector3::Zero());
        leg_type_map = boost::assign::map_list_of<leg_type, std::string>(RLEG, "rleg")(LLEG, "lleg")(RARM, "rarm")(LARM, "larm");
    };
    ~gait_generator () {
      if ( preview_controller_ptr != NULL ) {
        delete preview_controller_ptr;
        preview_controller_ptr = NULL;
      }
    };
    void initialize_gait_parameter (const hrp::Vector3& cog,
                                    const std::vector<step_node>& initial_support_leg_steps,
                                    const std::vector<step_node>& initial_swing_leg_dst_steps,
                                    const double delay = 1.6);
    bool proc_one_tick ();
    void append_footstep_nodes (const std::vector<std::string>& _legs, const std::vector<coordinates>& _fss)
    {
        std::vector<step_node> tmp_sns;
        for (size_t i = 0; i < _legs.size(); i++) {
            tmp_sns.push_back(step_node(_legs[i], _fss[i], lcg.get_default_step_height(), default_step_time, lcg.get_toe_angle(), lcg.get_heel_angle()));
        }
        footstep_nodes_list.push_back(tmp_sns);
    };
    void append_footstep_nodes (const std::vector<std::string>& _legs, const std::vector<coordinates>& _fss, const double _step_height, const double _step_time, const double _toe_angle, const double _heel_angle)
    {
        std::vector<step_node> tmp_sns;
        for (size_t i = 0; i < _legs.size(); i++) {
            tmp_sns.push_back(step_node(_legs[i], _fss[i], _step_height, _step_time, _toe_angle, _heel_angle));
        }
        footstep_nodes_list.push_back(tmp_sns);
    };
    void clear_footstep_nodes_list () {
        footstep_nodes_list.clear();
        overwrite_footstep_nodes_list.clear();
        overwrite_footstep_index = 0;
    };
    bool go_pos_param_2_footstep_nodes_list (const double goal_x, const double goal_y, const double goal_theta, /* [mm] [mm] [deg] */
                                             const std::vector<coordinates>& initial_support_legs_coords, coordinates start_ref_coords,
                                             const std::vector<leg_type>& initial_support_legs,
                                             const bool is_initialize = true) {
        std::vector< std::vector<step_node> > new_footstep_nodes_list;
        return go_pos_param_2_footstep_nodes_list (goal_x, goal_y, goal_theta,
                                                   initial_support_legs_coords, start_ref_coords,
                                                   initial_support_legs,
                                                   new_footstep_nodes_list,
                                                   is_initialize);
    };
    bool go_pos_param_2_footstep_nodes_list (const double goal_x, const double goal_y, const double goal_theta, /* [mm] [mm] [deg] */
                                             const std::vector<coordinates>& initial_support_legs_coords, coordinates start_ref_coords,
                                             const std::vector<leg_type>& initial_support_legs,
                                             std::vector< std::vector<step_node> >& new_footstep_nodes_list,
                                             const bool is_initialize = true);
    void go_pos_param_2_footstep_nodes_list_core (const double goal_x, const double goal_y, const double goal_theta, /* [mm] [mm] [deg] */
                                                  const std::vector<coordinates>& initial_support_legs_coords, coordinates start_ref_coords,
                                                  const std::vector<leg_type>& initial_support_legs,
                                                  std::vector< std::vector<step_node> >& new_footstep_nodes_list,
                                                  const bool is_initialize, const size_t overwritable_fs_index) const;
    void go_single_step_param_2_footstep_nodes_list (const double goal_x, const double goal_y, const double goal_z, const double goal_theta, /* [mm] [mm] [mm] [deg] */
                                               const std::string& tmp_swing_leg,
                                               const coordinates& _support_leg_coords);
    void initialize_velocity_mode (const coordinates& _ref_coords,
				   const double vel_x, const double vel_y, const double vel_theta, /* [mm/s] [mm/s] [deg/s] */
                                   const std::vector<leg_type>& current_legs);
    void finalize_velocity_mode ();
    void append_finalize_footstep ()
    {
      append_finalize_footstep(footstep_nodes_list);
    };
    void append_finalize_footstep (std::vector< std::vector<step_node> >& _footstep_nodes_list) const
    {
      std::vector<step_node> sns = _footstep_nodes_list[_footstep_nodes_list.size()-2];
      for (size_t i = 0; i < sns.size(); i++) {
          sns.at(i).step_height = sns.at(i).toe_angle = sns.at(i).heel_angle = 0.0;
      }
      _footstep_nodes_list.push_back(sns);
    };
    void emergency_stop ()
    {
      if (!footstep_nodes_list.empty()) {
        velocity_mode_flg = VEL_IDLING;
        emergency_flg = EMERGENCY_STOP;
      }
    };
    /* parameter setting */
    void set_default_step_time (const double _default_step_time) { default_step_time = _default_step_time; };
    void set_default_double_support_ratio_before (const double _default_double_support_ratio_before) { default_double_support_ratio_before = _default_double_support_ratio_before; };
    void set_default_double_support_ratio_after (const double _default_double_support_ratio_after) { default_double_support_ratio_after = _default_double_support_ratio_after; };
    void set_default_double_support_static_ratio_before (const double _default_double_support_static_ratio_before) { default_double_support_static_ratio_before = _default_double_support_static_ratio_before; };
    void set_default_double_support_static_ratio_after (const double _default_double_support_static_ratio_after) { default_double_support_static_ratio_after = _default_double_support_static_ratio_after; };
    void set_default_double_support_ratio_swing_before (const double _default_double_support_ratio_swing_before) { default_double_support_ratio_swing_before = _default_double_support_ratio_swing_before; };
    void set_default_double_support_ratio_swing_after (const double _default_double_support_ratio_swing_after) { default_double_support_ratio_swing_after = _default_double_support_ratio_swing_after; };
    void set_default_zmp_offsets(const std::vector<hrp::Vector3>& tmp) { rg.set_default_zmp_offsets(tmp); };
    void set_toe_zmp_offset_x (const double _off) { rg.set_toe_zmp_offset_x(_off); };
    void set_heel_zmp_offset_x (const double _off) { rg.set_heel_zmp_offset_x(_off); };
    void set_use_toe_heel_transition (const double _u) { rg.set_use_toe_heel_transition(_u); };
    void set_zmp_weight_map (const std::map<leg_type, double> _map) { rg.set_zmp_weight_map(_map); };
    void set_default_step_height(const double _tmp) { lcg.set_default_step_height(_tmp); };
    void set_default_top_ratio(const double _tmp) { lcg.set_default_top_ratio(_tmp); };
    void set_velocity_param (const double vel_x, const double vel_y, const double vel_theta) /* [mm/s] [mm/s] [deg/s] */
    {
      vel_param.set(vel_x, vel_y, vel_theta);
    };
    void set_offset_velocity_param (const double vel_x, const double vel_y, const double vel_theta) /* [mm/s] [mm/s] [deg/s] */
    {
      offset_vel_param.set(vel_x, vel_y, vel_theta);
    };
    void set_stride_parameters (const double _stride_fwd_x, const double _stride_y, const double _stride_theta, const double _stride_bwd_x)
    {
      footstep_param.stride_fwd_x = _stride_fwd_x;
      footstep_param.stride_y = _stride_y;
      footstep_param.stride_theta = _stride_theta;
      footstep_param.stride_bwd_x = _stride_bwd_x;
    };
    void set_use_inside_step_limitation(const bool uu) { use_inside_step_limitation = uu; };
    void set_default_orbit_type (const orbit_type type) { lcg.set_default_orbit_type(type); };
    void set_swing_trajectory_delay_time_offset (const double _time_offset) { lcg.set_swing_trajectory_delay_time_offset(_time_offset); };
    void set_swing_trajectory_final_distance_weight (const double _final_distance_weight) { lcg.set_swing_trajectory_final_distance_weight(_final_distance_weight); };
    void set_swing_leg_take_off_vel (const hrp::Vector3 _take_off_vel) { lcg.set_swing_leg_take_off_vel(_take_off_vel); };
    void set_swing_leg_landing_vel (const hrp::Vector3 _landing_vel) { lcg.set_swing_leg_landing_vel(_landing_vel); };
    void set_stair_trajectory_way_point_offset (const hrp::Vector3 _offset) { lcg.set_stair_trajectory_way_point_offset(_offset); };
    void set_cycloid_delay_kick_point_offset (const hrp::Vector3 _offset) { lcg.set_cycloid_delay_kick_point_offset(_offset); };
    void set_gravitational_acceleration (const double ga) { gravitational_acceleration = ga; };
    void set_toe_pos_offset_x (const double _offx) { lcg.set_toe_pos_offset_x(_offx); };
    void set_heel_pos_offset_x (const double _offx) { lcg.set_heel_pos_offset_x(_offx); };
    void set_toe_angle (const double _angle) { lcg.set_toe_angle(_angle); };
    void set_heel_angle (const double _angle) { lcg.set_heel_angle(_angle); };
    bool set_toe_heel_phase_ratio (const std::vector<double>& ratio) { return thp.set_toe_heel_phase_ratio(ratio); };
    void set_use_toe_joint (const bool ut) { lcg.set_use_toe_joint(ut); };
    void set_leg_default_translate_pos (const std::vector<hrp::Vector3>& off) { footstep_param.leg_default_translate_pos = off;};
    void set_optional_go_pos_finalize_footstep_num (const size_t num) { optional_go_pos_finalize_footstep_num = num; };
    void set_all_limbs (const std::vector<std::string>& _all_limbs) { all_limbs = _all_limbs; };
    void set_overwritable_footstep_index_offset (const size_t _of) { overwritable_footstep_index_offset = _of;};
    void set_foot_steps_list (const std::vector< std::vector<step_node> >& fnsl)
    {
        clear_footstep_nodes_list();
        footstep_nodes_list = fnsl;
        append_finalize_footstep();
        print_footstep_nodes_list();
    };
    void set_overwrite_foot_steps_list (const std::vector< std::vector<step_node> >& fnsl)
    {
        overwrite_footstep_nodes_list.clear();
        overwrite_footstep_nodes_list = fnsl;
        append_finalize_footstep(overwrite_footstep_nodes_list);
        print_footstep_nodes_list(overwrite_footstep_nodes_list);
    };
    /* Get overwritable footstep index. For example, if overwritable_footstep_index_offset = 1, overwrite next footstep. If overwritable_footstep_index_offset = 0, overwrite current swinging footstep. */
    size_t get_overwritable_index () const
    {
        return lcg.get_footstep_index()+overwritable_footstep_index_offset;
    };
    bool set_overwrite_foot_step_index (const size_t idx)
    {
        if (idx >= get_overwritable_index()) {
            overwrite_footstep_index = idx;
            return true;
        } else {
            return false;
        }
    };
    bool get_footstep_nodes_by_index (std::vector<step_node>& csl, const size_t idx) const
    {
        if (footstep_nodes_list.size()-1 >= idx) {
            csl = footstep_nodes_list.at(idx);
            return true;
        } else {
            return false;
        }
    };
    void print_footstep_nodes_list (const std::vector< std::vector<step_node> > _footstep_nodes_list) const
    {
        for (size_t i = 0; i < _footstep_nodes_list.size(); i++) {
            for (size_t j = 0; j < _footstep_nodes_list.at(i).size(); j++) {
                std::cerr << "[" << i << "] " << _footstep_nodes_list.at(i).at(j) << std::endl;
            }
        }
    };
    void print_footstep_nodes_list () const
    {
      print_footstep_nodes_list(footstep_nodes_list);
    };
    /* parameter getting */
    const hrp::Vector3& get_cog () const { return cog; };
    hrp::Vector3 get_cog_vel () const {
      double refcog_vel[3];
      preview_controller_ptr->get_refcog_vel(refcog_vel);
      return hrp::Vector3(refcog_vel[0], refcog_vel[1], refcog_vel[2]);
    };
    hrp::Vector3 get_cog_acc () const {
      double refcog_acc[3];
      preview_controller_ptr->get_refcog_acc(refcog_acc);
      return hrp::Vector3(refcog_acc[0], refcog_acc[1], refcog_acc[2]);
    };
    const hrp::Vector3& get_refzmp () {
      //for skate
      /* hrp::Vector3 modif_refzmp = refzmp; */
      /* modif_refzmp(0) += lcg.get_skate_acc()(0) * cog(2) / gravitational_acceleration; */
      refzmp(0) += lcg.get_skate_acc()(0) * cog(2) / gravitational_acceleration;
      return refzmp;
    };
    hrp::Vector3 get_cart_zmp () const
    {
        double czmp[3];
        preview_controller_ptr->get_cart_zmp(czmp);
        return hrp::Vector3(czmp[0], czmp[1], czmp[2]);
    };
    std::vector<std::string> convert_leg_types_to_names (const std::vector<leg_type>& lts) const {
      std::vector<std::string> ret;
      for (std::vector<leg_type>::const_iterator it = lts.begin(); it != lts.end(); it++) {
          ret.push_back(leg_type_map.find(*it)->second);
      }
      return ret;
    };
    const std::vector<hrp::Vector3>& get_swing_foot_zmp_offsets () const { return swing_foot_zmp_offsets;};
    std::vector<hrp::Vector3> get_support_foot_zmp_offsets () const {
      std::vector<hrp::Vector3> ret;
      for (size_t i = 0; i < lcg.get_support_leg_types().size(); i++) {
          ret.push_back(rg.get_default_zmp_offset(lcg.get_support_leg_types().at(i)));
      }
      return ret;
    };
    double get_toe_zmp_offset_x () const { return rg.get_toe_zmp_offset_x(); };
    double get_heel_zmp_offset_x () const { return rg.get_heel_zmp_offset_x(); };
    bool get_use_toe_heel_transition () const { return rg.get_use_toe_heel_transition(); };
    const std::map<leg_type, double> get_zmp_weight_map () const { return rg.get_zmp_weight_map(); };
    void proc_zmp_weight_map_interpolation () { return rg.proc_zmp_weight_map_interpolation(); };
    std::vector<std::string> get_footstep_front_leg_names () const {
      std::vector<leg_type> lts;
      for (size_t i = 0; i < footstep_nodes_list[0].size(); i++) {
          lts.push_back(footstep_nodes_list[0].at(i).l_r);
      }
      return convert_leg_types_to_names(lts);
    };
    std::vector<std::string> get_footstep_back_leg_names () const {
      std::vector<leg_type> lts;
      for (size_t i = 0; i < footstep_nodes_list.back().size(); i++) {
          lts.push_back(footstep_nodes_list.back().at(i).l_r);
      }
      return convert_leg_types_to_names(lts);
    };
    const std::vector<std::string> get_support_leg_names() const { return convert_leg_types_to_names(lcg.get_support_leg_types());};
    const std::vector<std::string> get_swing_leg_names() const { return convert_leg_types_to_names(lcg.get_swing_leg_types());};
    const std::vector<step_node>& get_swing_leg_steps() const { return lcg.get_swing_leg_steps(); };
    const std::vector<step_node>& get_support_leg_steps() const { return lcg.get_support_leg_steps(); };
    const std::vector<step_node>& get_swing_leg_src_steps() const { return lcg.get_swing_leg_src_steps(); };
    const std::vector<step_node>& get_swing_leg_dst_steps() const { return lcg.get_swing_leg_dst_steps(); };
    const coordinates get_dst_foot_midcoords() const /* get foot_midcoords calculated from swing_leg_dst_coords */
    {
      coordinates tmp(lcg.get_swing_leg_dst_steps().front().worldcoords);
      tmp.pos += tmp.rot * hrp::Vector3(-1*footstep_param.leg_default_translate_pos[lcg.get_swing_leg_dst_steps().front().l_r]);
      return tmp;
    };
    void get_swing_support_mid_coords(coordinates& ret) const { lcg.get_swing_support_mid_coords(ret); };
    void get_stride_parameters (double& _stride_fwd_x, double& _stride_y, double& _stride_theta, double& _stride_bwd_x) const
    {
      _stride_fwd_x = footstep_param.stride_fwd_x;
      _stride_y = footstep_param.stride_y;
      _stride_theta = footstep_param.stride_theta;
      _stride_bwd_x = footstep_param.stride_bwd_x;
    };
    size_t get_footstep_index() const { return lcg.get_footstep_index(); };
    size_t get_lcg_count() const { return lcg.get_lcg_count(); };
    size_t get_one_step_count() const { return lcg.get_one_step_count(); };
    double get_current_swing_time(const size_t idx) const { return lcg.get_current_swing_time(idx); };
    std::vector<leg_type> get_current_support_states() const { return lcg.get_current_support_states();};
    double get_default_step_time () const { return default_step_time; };
    double get_default_step_height () const { return lcg.get_default_step_height(); };
    double get_default_double_support_ratio_before () const { return default_double_support_ratio_before; };
    double get_default_double_support_ratio_after () const { return default_double_support_ratio_after; };
    double get_default_double_support_static_ratio_before () const { return default_double_support_static_ratio_before; };
    double get_default_double_support_static_ratio_after () const { return default_double_support_static_ratio_after; };
    double get_default_double_support_ratio_swing_before () const {return default_double_support_ratio_swing_before; };
    double get_default_double_support_ratio_swing_after () const {return default_double_support_ratio_swing_after; };
    std::vector< std::vector<step_node> > get_remaining_footstep_nodes_list () const
    {
        std::vector< std::vector<step_node> > fsnl;
        size_t fsl_size = (footstep_nodes_list.size()>lcg.get_footstep_index() ? footstep_nodes_list.size()-lcg.get_footstep_index() : 0);
        // The rest of fsl are swing dst coords from now.
        for (size_t i = 0; i < fsl_size; i++) {
            fsnl.push_back(footstep_nodes_list[i+lcg.get_footstep_index()]);
        }
        return fsnl;
    };
    orbit_type get_default_orbit_type () const { return lcg.get_default_orbit_type(); };
    double get_swing_trajectory_delay_time_offset () const { return lcg.get_swing_trajectory_delay_time_offset(); };
    double get_swing_trajectory_final_distance_weight () const { return lcg.get_swing_trajectory_final_distance_weight(); };
    hrp::Vector3 get_stair_trajectory_way_point_offset () const { return lcg.get_stair_trajectory_way_point_offset(); };
    hrp::Vector3 get_cycloid_delay_kick_point_offset () const { return lcg.get_cycloid_delay_kick_point_offset(); };
    hrp::Vector3 get_swing_leg_take_off_vel() const {return lcg.get_swing_leg_take_off_vel(); };
    hrp::Vector3 get_swing_leg_landing_vel() const {return lcg.get_swing_leg_landing_vel(); };
    hrp::Vector3 get_skate_acc() const { return lcg.get_skate_acc();};
    double get_gravitational_acceleration () const { return gravitational_acceleration; } ;
    double get_toe_pos_offset_x () const { return lcg.get_toe_pos_offset_x(); };
    double get_heel_pos_offset_x () const { return lcg.get_heel_pos_offset_x(); };
    double get_toe_angle () const { return lcg.get_toe_angle(); };
    double get_heel_angle () const { return lcg.get_heel_angle(); };
    double get_foot_dif_rot_angle () const { return lcg.get_foot_dif_rot_angle(); };
    void get_toe_heel_phase_ratio (std::vector<double>& ratio) const { thp.get_toe_heel_phase_ratio(ratio); };
    int get_NUM_TH_PHASES () const { return thp.get_NUM_TH_PHASES(); };
    bool get_use_toe_joint () const { return lcg.get_use_toe_joint(); };
    void get_leg_default_translate_pos (std::vector<hrp::Vector3>& off) const { off = footstep_param.leg_default_translate_pos; };
    size_t get_overwritable_footstep_index_offset () const { return overwritable_footstep_index_offset; };
    const std::vector<leg_type> calc_counter_leg_types_from_footstep_nodes (const std::vector<step_node>& fns, std::vector<std::string> _all_limbs) const;
    const std::map<leg_type, std::string> get_leg_type_map () const { return leg_type_map; };
    size_t get_optional_go_pos_finalize_footstep_num () const { return optional_go_pos_finalize_footstep_num; };
    bool is_finalizing (const double tm) const { return ((preview_controller_ptr->get_delay()*2 - default_step_time/dt)-finalize_count) <= (tm/dt)-1; };
    size_t get_overwrite_check_timing () const { return static_cast<size_t>(footstep_nodes_list[lcg.get_footstep_index()][0].step_time/dt * 0.5) - 1;}; // Almost middle of step time
    void print_param (const std::string& print_str = "") const
    {
        double stride_fwd_x, stride_y, stride_th, stride_bwd_x;
        get_stride_parameters(stride_fwd_x, stride_y, stride_th, stride_bwd_x);
        std::cerr << "[" << print_str << "]   stride_parameter = " << stride_fwd_x << "[m], " << stride_y << "[m], " << stride_th << "[deg], " << stride_bwd_x << "[m]" << std::endl;
        std::cerr << "[" << print_str << "]   leg_default_translate_pos = ";
        for (size_t i = 0; i < footstep_param.leg_default_translate_pos.size(); i++) {
            std::cerr << footstep_param.leg_default_translate_pos[i].format(Eigen::IOFormat(Eigen::StreamPrecision, 0, ", ", ", ", "", "", "    [", "]"));
        }
        std::cerr << std::endl;
        std::cerr << "[" << print_str << "]   default_step_time = " << get_default_step_time() << "[s]" << std::endl;
        std::cerr << "[" << print_str << "]   default_step_height = " << get_default_step_height() << "[m]" << std::endl;
        std::cerr << "[" << print_str << "]   default_double_support_ratio_before = " << get_default_double_support_ratio_before() << ", default_double_support_ratio_before = " << get_default_double_support_ratio_after() << ", default_double_support_static_ratio_before = " << get_default_double_support_static_ratio_before() << ", default_double_support_static_ratio_after = " << get_default_double_support_static_ratio_after() << ", default_double_support_static_ratio_swing_before = " << get_default_double_support_ratio_swing_before() << ", default_double_support_static_ratio_swing_after = " << get_default_double_support_ratio_swing_after() << std::endl;
        std::cerr << "[" << print_str << "]   default_orbit_type = ";
        if (get_default_orbit_type() == SHUFFLING) {
            std::cerr << "SHUFFLING" << std::endl;
        } else if (get_default_orbit_type() == CYCLOID) {
            std::cerr << "CYCLOID" << std::endl;
        } else if (get_default_orbit_type() == RECTANGLE) {
            std::cerr << "RECTANGLE" << std::endl;
        } else if (get_default_orbit_type() == STAIR) {
            std::cerr << "STAIR" << std::endl;
        } else if (get_default_orbit_type() == CYCLOIDDELAY) {
            std::cerr << "CYCLOIDDELAY" << std::endl;
        } else if (get_default_orbit_type() == CYCLOIDDELAYKICK) {
            std::cerr << "CYCLOIDDELAYKICK" << std::endl;
        } else if (get_default_orbit_type() == CROSS) {
            std::cerr << "CROSS" << std::endl;
        }
        std::cerr << "[" << print_str << "]   swing_trajectory_delay_time_offset = " << get_swing_trajectory_delay_time_offset() << "[s], swing_trajectory_final_distance_weight = " << get_swing_trajectory_final_distance_weight() << std::endl;
        hrp::Vector3 tmpv;
        tmpv = get_stair_trajectory_way_point_offset();
        std::cerr << "[" << print_str << "]   stair_trajectory_way_point_offset = " << tmpv.format(Eigen::IOFormat(Eigen::StreamPrecision, 0, ", ", ", ", "", "", "    [", "]")) << "[m]" << std::endl;
        tmpv = get_cycloid_delay_kick_point_offset();
        std::cerr << "[" << print_str << "]   cycloid_delay_kick_point_offset = " << tmpv.format(Eigen::IOFormat(Eigen::StreamPrecision, 0, ", ", ", ", "", "", "    [", "]")) << "[m]" << std::endl;
        tmpv = get_swing_leg_take_off_vel();
        std::cerr << "[" << print_str << "]   swing_leg_take_off_vel = " << tmpv.format(Eigen::IOFormat(Eigen::StreamPrecision, 0, ", ", ", ", "", "", "    [", "]")) << "[m]" << std::endl;
        tmpv = get_swing_leg_landing_vel();
        std::cerr << "[" << print_str << "]   swing_leg_landing_vel = " << tmpv.format(Eigen::IOFormat(Eigen::StreamPrecision, 0, ", ", ", ", "", "", "    [", "]")) << "[m]" << std::endl;
        std::cerr << "[" << print_str << "]   gravitational_acceleration = " << get_gravitational_acceleration() << "[m/s^2]" << std::endl;
        std::cerr << "[" << print_str << "]   toe_pos_offset_x = " << get_toe_pos_offset_x() << "[mm], heel_pos_offset_x = " << get_heel_pos_offset_x() << "[mm]" << std::endl;
        std::cerr << "[" << print_str << "]   toe_zmp_offset_x = " << get_toe_zmp_offset_x() << "[mm], heel_zmp_offset_x = " << get_heel_zmp_offset_x() << "[mm]" << std::endl;
        std::cerr << "[" << print_str << "]   toe_angle = " << get_toe_angle() << "[deg]" << std::endl;
        std::cerr << "[" << print_str << "]   heel_angle = " << get_heel_angle() << "[deg]" << std::endl;
        std::cerr << "[" << print_str << "]   use_toe_joint = " << (get_use_toe_joint()?"true":"false") << ", use_toe_heel_transition = " << (get_use_toe_heel_transition()?"true":"false") << std::endl;
        std::vector<double> tmp_ratio(get_NUM_TH_PHASES(), 0.0);
        get_toe_heel_phase_ratio(tmp_ratio);
        std::cerr << "[" << print_str << "]   toe_heel_phase_ratio = [";
        for (int i = 0; i < get_NUM_TH_PHASES(); i++) std::cerr << tmp_ratio[i] << " ";
        std::cerr << "]" << std::endl;
        std::cerr << "[" << print_str << "]   optional_go_pos_finalize_footstep_num = " << optional_go_pos_finalize_footstep_num << ", overwritable_footstep_index_offset = " << overwritable_footstep_index_offset << std::endl;
    };
  };
}
#endif /* GAITGENERATOR_H */
