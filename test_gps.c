#include <stdio.h>
#include "gps_kalman.h"

int main()
{
    double in_longitude, in_latitude, out_longitude, out_latitude;
    FILE *fp = fopen("../data/data1.txt", "r");
    FILE *fp_out = fopen("../data/data1_out.txt", "w+");
    if (fp == NULL)
    {
        return -1;
    }
    gps_filter_t *gps = gps_init();
    while (!feof(fp))
    {
        fscanf(fp, "%lf %lf", &in_longitude, &in_latitude);
        gps_filter(gps, in_longitude, in_latitude, &out_longitude, &out_latitude);
        fprintf(fp_out, "%lf %lf %lf %lf\n", in_longitude, in_latitude, out_longitude, out_latitude);
    }
    gps_de_init(gps);
    fclose(fp);
    fclose(fp_out);
    return 0;
}