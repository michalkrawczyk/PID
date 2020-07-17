//
// Created by Student235325 on 15.07.2020.
//

#ifndef PID_PID_HPP
#define PID_PID_HPP

#if defined(__cplusplus) && (__cplusplus >= 201103L)
#include <cstdint>
#include <limits>

namespace pid {
    enum class AntiWindup_t : uint8_t {
        NONE,
        BACK_CALCULATION,
        CLAMPING,
        SATURATION
    };

    class PID {
    public:
        explicit PID(const double &p, const double &i, const double &d);

        explicit PID(const double &p, const double &d);

        explicit PID(const double &p);

        void idealOrPararell(bool is_ideal);    //serial won't be included
        double calculate(const double &set_point, const double &p_v, const double &dt = 1.0);

        void setAntiWindupMode(const AntiWindup_t &mode);

        void setBoundries(const double &upper_bound, const double &lower_bound = 0.0);


    private:
        double _kp, _ki, _kd;

        double _up_bound{std::numeric_limits<double>::max() / 2};
        double _low_bound{std::numeric_limits<double>::lowest() / 2};
        double _integral{0.0}, _last_error{0.0};
        //pointer to set_point

        bool _is_ideal{false};
        AntiWindup_t _a_w_mode{AntiWindup_t::NONE};

        double Integrate(double &error, const double &dt);

        double Derive(double &error, const double &dt);

        double calculateCV(double &error, const double &dt);
    };
}
#endif //defined(__cplusplus) && (__cplusplus >= 201103L)

#endif //PID_PID_HPP
