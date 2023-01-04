#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define n 6
//*************************************************************************
//coordinats of target
//-------------------------------------------------------------------------

typedef struct Ephi {
    double x = -14081.752701e+3; // x coordinat, m
    double y = 18358.958252e+3; // y coordinat, m
    double z = 10861.302124e+3; // z coordinat, m
    double x_vel = -1.02576358e+3; // component x of velocity, m/s 
    double y_vel = 1.08672147e+3; // component y of velocity, m/s
    double z_vel = -3.15732343e+3; // component z of velocity, ms

}Ephi;

//-------------------------------------------------------------------------

//*************************************************************************
// constatns.

const double mu = 398600.4418e+9; //- earth gravity
const double J = 1082625.75e-9; //- Ephi zonal harmonic decomposition of the geopotential in a series of spherical functions;
const double a = 6378136;       //- the equatorial radius of the earth;
const double w = 7.292115e-5;   // Rad/s is the angular velocity of the earth's rotation;
//-------------------------------------------------------------------------

//double *XYZ(double *x,double h,double t){
//    double *arr;
//    arr=(double*)malloc(6*sizeof (double));
//    double S=
//
//}
//*************************************************************************
//-------------------------------------------------------------------------
// Title: len
//-------------------------------------------------------------------------
// Input parameters:
//          *x={x,y,z}-coordinates;
//
//-------------------------------------------------------------------------

//return value: len-length of Ephi;
//-------------------------------------------------------------------------
double len(double* x) {
    double len = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
    return len;
}
//-------------------------------------------------------------------------



//*************************************************************************
//-------------------------------------------------------------------------
// Title: f
//-------------------------------------------------------------------------
// Input parameters:
//          *x={x,y,z}-coordinates;
//          *y={Vx,Vy,Vz}-velocity;
//-------------------------------------------------------------------------
//Output parameters:
//                    arr[0] -dx/dt;
//                    arr[1] - dy/dt;
//                    arr[2] - dz/dt;
//                    arr[3] - dVx/dt;
//                    arr[4] - dVy/dt;
//                    arr[5] - dVz/dt;

//return value: arr - array of solved equations;
//-------------------------------------------------------------------------


double* f(double* x) {
    double r = len(x);
    double Ak = -mu / (r * r * r);
    double Ck = 1 - 5 * x[2] * x[2] / (r * r);
    double Bk = (-3 / 2) * J * mu * a * a / pow(r, 5) * Ck;
    double* arr;
    arr = (double*)malloc(6 * sizeof(double));
    arr[0] = x[3];
    arr[1] = x[4];
    arr[2] = x[5];

    arr[3] = Ak * x[0] + Bk * x[0] + w * w * x[0] + 2 * w * x[4];
    arr[4] = Ak * x[1] + Bk * x[1] + w * w * x[1] - 2 * w * x[3];
    arr[5] = Ak * x[2] - (3 / 2) * J * mu * a * a / pow(r, 5) * x[2] * (3 - 5 * x[2] * x[2] / (r * r));
    return arr;

}



//*************************************************************************
//-------------------------------------------------------------------------
// Title: multyM
//-------------------------------------------------------------------------
// Input parameters:
//          *x={x,y,z};
//          k-const;
//-------------------------------------------------------------------------
//Output parameters:
//                    arr[i]=x[i]*k
//
//return value: arr - multiplying a Ephi by a constant;
//-------------------------------------------------------------------------

double* multyM(double* x, double k) {
    double* arr;
    arr = (double*)malloc(6 * sizeof(double));
    for (int i = 0;i < n;i++) {
        arr[i] = x[i] * k;
    }
    return arr;
}
//-------------------------------------------------------------------------




//*************************************************************************
//-------------------------------------------------------------------------
// Title: addM
//-------------------------------------------------------------------------
// Input parameters:
//          *x={x,y,z}-Ephi;
//          *y={x,y,z}-Ephi;
//-------------------------------------------------------------------------
//Output parameters:
//                    arr[i]=x[i]+y[i];
//
//return value: arr - adding arrays;
//-------------------------------------------------------------------------
double* addM(double* x, double* y) {
    double* arr;
    arr = (double*)malloc(n * sizeof(double));
    for (int i = 0;i < n;i++) {
        arr[i] = x[i] + y[i];
    }
    return arr;
}

//-------------------------------------------------------------------------




double tu(double te, double ti) {
    return (ti - te) * 3600;
}

double* Runge(double* x, double step, double te, double ti) {
    double h = tu(te, ti) / step;
    double* k1, * k2, * k3, * k4;
    double* rungeg;
    rungeg = (double*)malloc(n * sizeof(double));
    k1 = f(x);

    // printf("k1:%lf s %lf s %lf\n",k1[0],k1[1],k1[2]);
    k2 = f(addM(x, multyM(multyM(k1, 0.5), h)));
    // printf("k2:%lf s %lf s %lf\n",k2[0],k2[1],k2[2]);
    k3 = f(addM(x, multyM(multyM(k2, 0.5), h)));
    // printf("k3:%lf s %lf s %lf\n",k3[0],k3[1],k3[2]);
    k4 = f(addM(x, multyM(multyM(k3, 1), h)));
    //printf("k4:%lf s %lf s %lf\n",k4[0],k4[1],k4[2]);
    for (int i = 0;i < n;i++) {
        rungeg[i] = x[i] + (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) * h / 6;
    }
    // printf("r:%lf s %lf s %lf\n",rungeg[0],rungeg[1],rungeg[2]);
    return rungeg;
}

void PredictionEphemeris(double te, double ti, Ephi &x) {
    double arr[6];
    arr[0] = x.x;
    arr[1] = x.y;
    arr[2] = x.z;
    arr[3] = x.x_vel;
    arr[4] = x.y_vel;
    arr[5] = x.z_vel;

    double step = 100;
   
   

    double* q[100];

    q[0] = Runge(arr, step, te, ti);
    //printf("(x:%+2.lf y:%+.2lf z:%+.2lf)\n",q[0][0],q[0][1],q[0][2]);

    for (int i = 1;i < step;i++) {



        q[i] = Runge(q[i - 1], step, te, ti);


    }
    
        x.x= q[99][0];
        x.y = q[99][1];
        x.z = q[99][2];
        x.x_vel= q[99][3];
        x.y_vel= q[99][4];
        x.z_vel= q[99][5];

    

}
void printE(Ephi x) {
    printf("(x:%+2.lf y:%+.2lf z:%+.2lf Vx:%+.2lf Vy:%+.2lf Vz:%+.2lf)\n", x.x,x.y,x.z,x.x_vel,x.y_vel,x.z_vel);
}
int main(int argc, char* argv[])
{
 int step = 50;
    double te = 6.25;
    double ti = 6.5;
    Ephi z;
    z.x = -14081.752701e+3;
    z.y = 18358.958252e+3;
    z.z = 10861.302124e+3;
    z.x_vel = -1.02576358e+3;
    z.y_vel = 1.08672147e+3;
    z.z_vel = -3.15732343e+3;
    //double x[6] = { -14081.752701e+3,18358.958252e+3,10861.302124e+3,-1.02576358e+3,1.08672147e+3,-3.15732343e+3 };
    PredictionEphemeris( te, ti, z);

    printE(z);



    return 0;
}
