/* SCRIPT DEL EJERCICIO OBLIGATORIO DE LA LECCIÓN 3 DE FÍSICA COMPUTACIONAL */
/* SE RESUELVE COMPUTACIONALMENTE LA ECUACIÓN DE SCHRÖDINGER PARA UN POTENCIAL 
   ESCALÓN USANDO COMO FUNCIÓN DE ONDA INICIAL UNA GAUSSIANA */

#include <iostream>
#include <cmath>
#include <fstream>
#include <stdio.h>
#include <cstdlib>
#include "complex.h"

#define N 1000  
#define PI 3.14159265

using namespace std;

int main (void)
{
    int i,j,nciclos;
    double lambda,s,x0,sigma,k0,norma,t;
    double V[N+1];

    fcomplex phi[N+1];
    fcomplex chi[N+1];
    fcomplex A0[N+1];
    fcomplex alpha[N];
    fcomplex beta[N];
    fcomplex b[N+1];

    ofstream fichero1;
    ofstream fichero2;
    ofstream fichero3;

    fichero1.open("/Users/alzorrcarri/Documents/cphys/leccion3-obligatorio/funciondeonda.dat");
    fichero2.open("/Users/alzorrcarri/Documents/cphys/leccion3-obligatorio/norma.dat");
    fichero3.open("/Users/alzorrcarri/Documents/cphys/leccion3-obligatorio/probabilidad.dat");

    nciclos=20; //Empezamos con nciclos bajos, y subimos hasta N/4 como máximo

    /* Parámetros de la función de onda gaussiana */
    x0=N/4.0;
    sigma=N/16.0;
    k0=2*PI*nciclos/N;

    /* Espaciado temporal */
    s=1.0/(4.0*k0*k0);
    t=0;

    /* Establecemos el potencial */
    cout << "Valor de la altura de potencial:";
    cin >> lambda;
    for ( j = 0; j <= N; j++)
    {
        if ((0.4*N <= j) && (j<= 0.6*N))
        {
            V[j]=lambda*k0*k0;
        } else V[j]=0;
    }

    /* Establecemos la función de onda en el instante inicial t=0 y calculamos su norma */
    norma=0.0;
    phi[0].r=0; phi[N].r=0;
    phi[0].i=0; phi[N].i=0;
    for ( j = 1; j <= N-1; j++)
    {
        phi[j]=Cgauss(k0*j,exp((-1)*pow((j-x0),2)/(2*sigma*sigma)));
        norma=norma+pow(Cabs(phi[j]),2);
    }

    /* Normalizamos la función de onda inicial */
    for ( j = 0; j <= N; j++)
    {
        phi[j]=RCmul(1.0/sqrt(norma),phi[j]);
    }

    /* Escribimos la función de onda inicial en el fichero */
    for ( j = 0; j <= N; j++)
    {
        fichero1 << j << ',' << phi[j].r << ',' << phi[j].i << ',' << V[j] << '\n';
        fichero3 << j << ',' << pow(Cabs(phi[j]),2) << ',' << V[j] << '\n';
    }
    fichero1 << '\n';
    fichero3 << '\n';
    
    /* Generamos las componentes del vector A0, que no dependen del tiempo */
    for ( j = 0; j <= N+1; j++)
    {
        A0[j].r=-2.0-V[j];
        A0[j].i=2.0/s;
    }
    
    /* Generamos las componentes de alpha, que no dependen del tiempo */
    alpha[N-1]=Complex(0.0,0.0);
    for ( j = 2; j <= N; j++)
    {
        alpha[N-j]=Cadd(A0[j],alpha[N-j+1]);
        alpha[N-j]=Cdiv(Complex(-1.0,0.0),alpha[N-j]);
    }

    /* PROCESO ITERATIVO EN EL TIEMPO i */
    for ( i = 1; i <= 5000; i++)
    {
        t=t+s;
        /* Calculo el vector b */
        for ( j = 0; j <= N; j++)
        {
            b[j]=phi[j];
            b[j]=Cmul(Complex(0.0,1.0),b[j]);
            b[j]=RCmul(4.0/s,b[j]);
        }

        /* Calculo el vector beta */
        beta[N-1].r=0;
        beta[N-1].i=0;
        for ( j = 2; j <= N; j++)
        {
            beta[N-j]=Csub(b[N-j+1],beta[N-j+1]);
            beta[N-j]=Cmul(Cdiv(Complex(1.0,0.0),Cadd(A0[j],alpha[N-j+1])),beta[N-j]);
        }
        
        /* Calculamos el vector chi */
        chi[0].r=0; chi[0].i=0;
        chi[N].r=0; chi[N].i=0;
        for ( j = 1; j <= N-1; j++)
        {
            chi[j]=Cadd(Cmul(alpha[j-1],chi[j-1]),beta[j-1]);
        }
        
        /* Calculamos el vector phi para el nuevo instante de tiempo y su norma */
        norma=0.0;
        for ( j = 1; j <= N; j++)
        {
            phi[j]=Csub(chi[j],phi[j]);
            norma=norma+pow(Cabs(phi[j]),2);
        }

        /* Escribimos la función de onda en cada instante de tiempo en un fichero y la norma */
        for ( j = 0; j <= N; j++)
        {
            fichero1 << j << ',' << phi[j].r << ',' << phi[j].i << ',' << V[j] << '\n';
            fichero3 << j << ',' << pow(Cabs(phi[j]),2) << ',' << V[j] << '\n';
        }
        fichero1 << '\n';
        fichero2 << t << '\t' << norma << '\n';
        fichero3 << '\n';
        
        
    }
    
fichero1.close();
fichero2.close();
fichero3.close();
return 0;

}