//
// Created by liaoy on 2021/2/25.
//

#ifndef __GPS_KALMAN_H
#define __GPS_KALMAN_H

#include <stdbool.h>

/* 2 Dimension */
typedef struct {
    double x[2];     /* state: [0]-angle [1]-difference of angle, 2x1 */
    double A[2][2];  /* X(n)=A*X(n-1)+U(n),U(n)~N(0,q), 2x2 */
    double H[2];     /* Z(n)=H*X(n)+W(n),W(n)~N(0,r), 1x2   */
    double q[2];     /* process(predict) noise covariance,2x1 [q0,0; 0,q1] */
    double r;        /* measure noise variance */
    double p[2][2];  /* estimated error covariance,2x2 [p0 p1; p2 p3] */
    double gain[2];  /* 2x1 */
} kalman2_state;

typedef struct
{
    kalman2_state longitude_filter;
    kalman2_state latitude_filter;
    bool is_first;
} gps_filter_t;

gps_filter_t *gps_init(void);

bool gps_filter(gps_filter_t *gps_kalman_filter, double in_longitude, double in_latitude, double *out_longitude,
                double *out_latitude);

void gps_de_init(gps_filter_t *gps_kalman_filter);

#endif //__GPS_KALMAN_H
