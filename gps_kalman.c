//
// Created by liaoy on 2021/2/25.
//

#include "gps_kalman.h"
#include <stdlib.h>

static void kalman2_init(kalman2_state *state, const double init_x[2], double init_p[2][2])
{
    state->x[0]    = init_x[0];
    state->x[1]    = init_x[1];
    state->p[0][0] = init_p[0][0];
    state->p[0][1] = init_p[0][1];
    state->p[1][0] = init_p[1][0];
    state->p[1][1] = init_p[1][1];
    //state->A       = {{1, 0.1}, {0, 1}};
    state->A[0][0] = 1;
    state->A[0][1] = 0.1;
    state->A[1][0] = 0;
    state->A[1][1] = 1;
    //state->H       = {1,0};
    state->H[0]    = 1;
    state->H[1]    = 0;
    //state->q       = {{10e-6,0}, {0,10e-6}};
    state->q[0]    = 10e-7;
    state->q[1]    = 10e-7; /* estimated error covariance */
    state->r       = 10e-7; /* measure noise covariance */
}

static double kalman2_filter(kalman2_state *state, double z_measure)
{
    /* Step1: Predict */
    state->x[0] = state->A[0][0] * state->x[0] + state->A[0][1] * state->x[1];
    state->x[1] = state->A[1][0] * state->x[0] + state->A[1][1] * state->x[1];
    /* p(n|n-1)=A^2*p(n-1|n-1)+q */
    state->p[0][0] = state->A[0][0] * state->p[0][0] + state->A[0][1] * state->p[1][0] + state->q[0];
    state->p[0][1] = state->A[0][0] * state->p[0][1] + state->A[1][1] * state->p[1][1];
    state->p[1][0] = state->A[1][0] * state->p[0][0] + state->A[0][1] * state->p[1][0];
    state->p[1][1] = state->A[1][0] * state->p[0][1] + state->A[1][1] * state->p[1][1] + state->q[1];

    /* Step2: Measurement */
    /* gain = p * H^T * [r + H * p * H^T]^(-1), H^T means transpose. */
    double temp0 = state->p[0][0] * state->H[0] + state->p[0][1] * state->H[1];
    double temp1 = state->p[1][0] * state->H[0] + state->p[1][1] * state->H[1];
    double temp  = state->r + state->H[0] * temp0 + state->H[1] * temp1;
    state->gain[0] = temp0 / temp;
    state->gain[1] = temp1 / temp;
    /* x(n|n) = x(n|n-1) + gain(n) * [z_measure - H(n)*x(n|n-1)]*/
    temp = state->H[0] * state->x[0] + state->H[1] * state->x[1];
    state->x[0] = state->x[0] + state->gain[0] * (z_measure - temp);
    state->x[1] = state->x[1] + state->gain[1] * (z_measure - temp);

    /* Update @p: p(n|n) = [I - gain * H] * p(n|n-1) */
    state->p[0][0] = (1 - state->gain[0] * state->H[0]) * state->p[0][0];
    state->p[0][1] = (1 - state->gain[0] * state->H[1]) * state->p[0][1];
    state->p[1][0] = (1 - state->gain[1] * state->H[0]) * state->p[1][0];
    state->p[1][1] = (1 - state->gain[1] * state->H[1]) * state->p[1][1];

    return state->x[0];
}

static void kalman2_set_param(kalman2_state *state, const double q[2], double r)
{
    state->q[0] = q[0];
    state->q[1] = q[1];
    state->r = r;
}

gps_filter_t *gps_init(void)
{
    gps_filter_t *gps_kalman_filter = (gps_filter_t *) malloc(sizeof(gps_filter_t));
    gps_kalman_filter->is_first = true;
    return gps_kalman_filter;
}

bool gps_filter(gps_filter_t *gps_kalman_filter, double in_longitude, double in_latitude, double *out_longitude,
                double *out_latitude)
{
    if (gps_kalman_filter == NULL)
    {
        *out_longitude = in_longitude;
        *out_latitude = in_latitude;
        return false;
    }
    if (gps_kalman_filter->is_first)
    {
        double init_longitude[2] = {in_longitude, 0};
        double init_longitude_p[2][2] = {{0, 0},
                                         {0, 0}};
        kalman2_init(&gps_kalman_filter->longitude_filter, init_longitude, init_longitude_p);

        double init_latitude[2] = {in_latitude, 0};
        double init_latitude_p[2][2] = {{0, 0},
                                        {0, 0}};
        kalman2_init(&gps_kalman_filter->latitude_filter, init_latitude, init_latitude_p);
        gps_kalman_filter->is_first = false;
    }
    double process_noise_err[2] = {0.00001, 0.00001};
    double measurement_err = (in_longitude - gps_kalman_filter->longitude_filter.x[0]) *
                             (in_longitude - gps_kalman_filter->longitude_filter.x[0]);
    measurement_err += (in_latitude - gps_kalman_filter->latitude_filter.x[0]) *
                       (in_latitude - gps_kalman_filter->latitude_filter.x[0]);
    measurement_err *= 1300;
    kalman2_set_param(&gps_kalman_filter->longitude_filter, process_noise_err, measurement_err);
    kalman2_set_param(&gps_kalman_filter->latitude_filter, process_noise_err, measurement_err);

    *out_longitude = kalman2_filter(&gps_kalman_filter->longitude_filter, in_longitude);
    *out_latitude = kalman2_filter(&gps_kalman_filter->latitude_filter, in_latitude);

    return true;
}

void gps_de_init(gps_filter_t *gps_kalman_filter)
{
    if (gps_kalman_filter != NULL)
    {
        free(gps_kalman_filter);
        gps_kalman_filter = NULL;
    }
}