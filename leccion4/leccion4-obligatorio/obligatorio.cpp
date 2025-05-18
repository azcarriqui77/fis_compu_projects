/* COHETE A LA LUNA: PROBLEMA DE LOS TRES CUERPOS CON ALGORITMO DE RUNGE-KUTTA DE CUARTO ORDEN */

#include <iostream>
#include <cmath>
#include <fstream>
#include <stdio.h>
#include <cstdlib>

#define PI 3.14159265

using namespace std;

double rprima (double r, double phi, double t);

int main (void)
{
    int i,j,k,n;
    double G=6.67E-11, MT=5.9736E24, ML=0.07349E24, dTL=3.844E8, omega=2.6617E-6, RT=6.378160E6, RL=1.7374E6;
    double t,h,theta,v,delta,mu;
    double r,phi,pr,pphi,rL,phiL;
    double k1[4],k2[4],k3[4],k4[4];
    double H;

    ofstream fichero;
    ofstream constante;

    fichero.open("/Users/alzorrcarri/Documents/cphys/leccion4-obligatorio/data.dat");
    constante.open("/Users/alzorrcarri/Documents/cphys/leccion4-obligatorio/const.dat");

    /* Parámetros */
    n=1000000; h=1; t=0;
    delta=G*MT/pow(dTL,3); mu=ML/MT;

    /* Condiciones iniciales del cohete sobre la superficie terrestre. Reescalamos magnitudes. El cohete lo lanzamos desde el polo norte terrestre de forma tangencial. */
    r=RT; r=r/dTL;
    phi=PI/2;
    theta=2.99*PI/8.0;
    v=sqrt(2*G*MT/RT); v=0.991*v/dTL; //Fracción de la velocidad de escape
    pr=v*cos(theta-phi);
    pphi=r*v*sin(theta-phi);
    rL=dTL; rL=rL/dTL;
    phiL=0;

    /* Escribimos en el fichero las posiciones iniciales del cohete, la luna y la Tierra. */
    fichero << 0 << ',' << 0 << '\n';
    fichero << rL*cos(phiL) << ',' << rL*sin(phiL) << '\n';
    fichero << r*cos(phi) << ',' << r*sin(phi) << '\n';
    fichero << '\n';

    for ( i = 1; i <= n; i++)
    {
        H=0;
        /* MÉTODO DE RUNGE-KUTTA DE ORDEN CUATRO */
        /* Calculamos k1 */
        k1[0]=h*pr;
        k1[1]=h*pphi/(r*r);
        k1[2]=h*pphi*pphi/(r*r*r)-h*delta*(1/(r*r)+mu*(r-cos(phi-omega*t))/pow(rprima(r,phi,t),3));
        k1[3]=-h*delta*mu*r*sin(phi-omega*t)/pow(rprima(r,phi,t),3);

        /* Calculamos k2 */
        k2[0]=h*(pr+k1[2]/2);
        k2[1]=h*(pphi+k1[3]/2)/pow(r+k1[0]/2,2);
        k2[2]=h*pow(pphi+k1[3]/2,2)/pow(r+k1[0]/2,3)-h*delta*(1/pow(r+k1[0]/2,2)+mu*(r+k1[0]/2-cos(phi+k1[1]/2-omega*(t+h/2)))/pow(rprima(r+k1[0]/2,phi+k1[1]/2,t+h/2),3));
        k2[3]=-h*delta*mu*(r+k1[0]/2)*sin(phi+k1[1]/2-omega*(t+h/2))/pow(rprima(r+k1[0]/2,phi+k1[1]/2,t+h/2),3);

        /* Calculamos k3 */
        k3[0]=h*(pr+k2[2]/2);
        k3[1]=h*(pphi+k2[3]/2)/pow(r+k2[0]/2,2);
        k3[2]=h*pow(pphi+k2[3]/2,2)/pow(r+k2[0]/2,3)-h*delta*(1/pow(r+k2[0]/2,2)+mu*(r+k2[0]/2-cos(phi+k2[1]/2-omega*(t+h/2)))/pow(rprima(r+k2[0]/2,phi+k2[1]/2,t+h/2),3));
        k3[3]=-h*delta*mu*(r+k2[0]/2)*sin(phi+k2[1]/2-omega*(t+h/2))/pow(rprima(r+k2[0]/2,phi+k2[1]/2,t+h/2),3);

        /* Calculamos k4 */
        k4[0]=h*(pr+k3[2]);
        k4[1]=h*(pphi+k3[3])/pow(r+k3[0],2);
        k4[2]=h*pow(pphi+k3[3],2)/pow(r+k3[0],3)-h*delta*(1/pow(r+k3[0],2)+mu*(r+k3[0]-cos(phi+k3[1]-omega*(t+h)))/pow(rprima(r+k3[0],phi+k3[1],t+h),3));
        k4[3]=-h*delta*mu*(r+k3[0])*sin(phi+k3[1]-omega*(t+h))/pow(rprima(r+k3[0],phi+k3[1],t+h),3);

        /* Actualizamos el valor de las variables de movimiento del cohete y de la Luna */
        r=r+(k1[0]+2*k2[0]+2*k3[0]+k4[0])/6.0;
        phi=phi+(k1[1]+2*k2[1]+2*k3[1]+k4[1])/6.0;
        pr=pr+(k1[2]+2*k2[2]+2*k3[2]+k4[2])/6.0;
        pphi=pphi+(k1[3]+2*k2[3]+2*k3[3]+k4[3])/6.0;
        t=t+h;
        phiL=omega*t;

        /* Escribimos en el fichero los valores de posición del cohete y la Luna cada 100 pasos. */
        if (i%(n/100)==0)
        {
            fichero << 0 << ',' << 0 << '\n';
            fichero << rL*cos(phiL) << ',' << rL*sin(phiL) << '\n';
            fichero << r*cos(phi) << ',' << r*sin(phi) << '\n';
            fichero << '\n';
        }
        H=dTL*dTL*((pr*pr+pphi*pphi/(r*r))/2.0-omega*pphi)-G*(MT/r+ML/sqrt(1+r*r-2*r*cos(phi-omega*t)))/dTL;
        constante << t << '\t' << H << '\n';
        
    }

    fichero.close();
    constante.close();
    return 0;
}

double rprima (double r, double phi, double t)
{
    double aux,dTL=3.844E8,omega=2.6617E-6;
    aux=sqrt(1+r*r+-2*r*cos(phi-omega*t));
    return(aux);
}