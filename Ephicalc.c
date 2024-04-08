#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define n 6




// constatns.

const double mu = 398600.4418e+9; //- earth gravity
const double J = 1082625.75e-9; //- Ephi zonal harmonic decomposition of the geopotential in a series of spherical functions;
const double a = 6378136;       //- the equatorial radius of the earth;
const double w = 7.292115e-5;   // Rad/s is the angular velocity of the earth's rotation;



double len(double x[n]) {
    double len = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
    return len;
}




void function(double coordinate_before[n], double* intermidiate) {
    double r = len(coordinate_before);
    double Ak = -mu / (r * r * r);
    double Ck = 1 - 5 * coordinate_before[2] * coordinate_before[2] / (r * r);
    double Bk = (-3 / 2) * J * mu * a * a / pow(r, 5) * Ck;
    
    
    for (int i = 0; i<(int)n/2; ++i)
    {
        intermidiate[i] = coordinate_before[i + (int)n / 2];
    }
    

    intermidiate[3] = Ak * coordinate_before[0] + Bk * coordinate_before[0] + w * w * coordinate_before[0] + 2 * w * coordinate_before[4];
    intermidiate[4] = Ak * coordinate_before[1] + Bk * coordinate_before[1] + w * w * coordinate_before[1] - 2 * w * coordinate_before[3];
    intermidiate[5] = Ak * coordinate_before[2] - (3./ 2.) * J * mu * a * a / pow(r, 5) * coordinate_before[2] * (3 - 5 * coordinate_before[2] * coordinate_before[2] / (r * r));
    

}



void transformation(double coordinate_before[n], double m, double h, double* k)
{
    for (int i = 0; i < n; ++i)
    {
        k[i] = (k[i] * 0.5*h)+ coordinate_before[i];
    }
    
}


double tu(double te, double ti) {
    return (ti - te) * 3600;
}

void Runge(const double coordinate_before[n], const int step, double te, double ti,double * rungeg) {
    double h =(tu(te, ti) / (double)step);
    double  k1[n], k2[n], k3[n], k4[n];
    double  k1_1[n], k2_1[n], k3_1[n];
    
    double m[2] = { 0.5,1 };
    function(coordinate_before, k1);
    memcpy(k1_1, k1, sizeof(k1[0]) * n);
    transformation(coordinate_before, m[0], h, k1);
    function(k1, k2);
    memcpy(k2_1, k2, sizeof(k2[0]) * n);
    transformation(k1, m[0], h, k2);
    function(k2, k3);
    memcpy(k3_1, k3, sizeof(k3[0]) * n);
    transformation(k2, m[1], h, k3);
    function(k3, k4);
   
    for (int i = 0; i < n; i++) {
        rungeg[i] = coordinate_before[i] + (k1_1[i] + 2 * k2_1[i] + 2 * k3_1[i] + k4[i]) * h / 6;
    }
    
    
    
}
void printE(double x[n]) {
    printf("(x:%+2.lf y:%+.2lf z:%+.2lf Vx:%+.2lf Vy:%+.2lf Vz:%+.2lf)\n", x[0], x[1], x[2], x[3], x[4], x[5]);
}

void PredictionEphemeris(double te, double ti, double coordinate_before[n], double  coordinate_after[n], int step) {
   
    Runge(coordinate_before, step, te, ti, coordinate_after);
    for (int i = 0; i < step; i++)
    {
        Runge(coordinate_after, step, te, ti, coordinate_after);
    }
    
}


int main(int argc, char* argv[])
{
 
    double te = 6.25;
    double ti = 6.5;
    int step = 100;
    
    double coordinate_after[n] =  {0};
    double coordinate_before[n] = { -14081.752701e+3, 18358.958252e+3,10861.302124e+3,-1.02576358e+3,1.08672147e+3, -3.15732343e+3 };
    
    PredictionEphemeris(te, ti, coordinate_before, coordinate_after,step);
    printE(coordinate_after);

    



    return 0;
}
