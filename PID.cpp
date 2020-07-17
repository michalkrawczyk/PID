#if defined(__cplusplus) && (__cplusplus >= 201103L)
#include "PID.hpp"
#include <algorithm>

using namespace pid;

PID::PID(const double &p, const double &i, const double &d):
_kp(p),
_ki(i),
_kd(d)
{}

PID::PID(const double &p, const double &d):
PID(p, 0.0, d)
{}

PID::PID(const double &p):
PID(p, 0.0, 0.0)
{}


void PID::setAntiWindupMode(const AntiWindup_t &mode)
{
    _a_w_mode = mode;
}


void PID::idealOrPararell(bool is_ideal)
{
    _is_ideal = is_ideal;
}

void PID::setBoundries(const double &upper_bound, const double &lower_bound)
{
    _up_bound = upper_bound;
    _low_bound = lower_bound;
}

double PID::calculate(const double &set_point, const double &p_v, const double &dt)
{
    double err = set_point - p_v; //calculate error

    return calculateCV(err, dt);
}

double PID::Integrate(double &error, const double &dt)
{
    return _integral + error * dt;
}

double PID::Derive(double &error, const double &dt)
{
    return (error - _last_error) / dt;
}

double PID::calculateCV(double &error, const double &dt)
{
    double c_v;
    double temp_integral(Integrate(error, dt));

    double d_out(_kd * Derive(error, dt));
    double i_out(_ki * temp_integral);

    if(_is_ideal)
    {
        c_v = _kp * (error + i_out + d_out);
    }
    else    //parallel
    {
        c_v = _kp * error + i_out + d_out;
    }

    switch(_a_w_mode)
    {
        case AntiWindup_t::CLAMPING:
        {
            if (c_v > _up_bound || c_v < _low_bound)
            {
                c_v -= _is_ideal ? _kp * i_out : i_out;
                std::swap(_last_error, error);
                return c_v;
            }

            break;
        }
        case AntiWindup_t::BACK_CALCULATION:
        {
            //TODO: Back Calculation
            break;
        }
        case AntiWindup_t::SATURATION:
        {
            if (c_v > _up_bound)
            {
                c_v = _up_bound;
            }
            else if(c_v < _low_bound)
            {
                c_v = _low_bound;
            }
            break;
        }
        default:
            break;
    }


    std::swap(_last_error, error);
    std::swap(_integral, temp_integral);

    return c_v;
}
#endif