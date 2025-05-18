/* SIMULACIÓN DEL MODELO DE ISING CON EL MÉTODO DE MONTECARLO */
/* EJERCICIO OBLIGATORIO DE LA LECCIÓN 2 DE FÍSICA COMPTACIONAL */

#include <iostream>
#include <cmath>
#include <fstream>
#include <stdio.h>
#include "gsl_rng.h"

using namespace std;

gsl_rng *tau;

int main (void)
{
    /* Definimos inicialmente algunas variables necesarias.*/
    int i,j,k,random,pasos,dim;
    int aleatorio;
    double energia,p,temp,aux;
    ofstream fichero;
    
    /* Pedimos los valores iniciales de la dimensión de la cuadrícula, la temperatura del sistema y el número de pasos Montecarlo.*/
    cout << "Introducir la dimensión de la cuadrícula:";
    cin >> dim;
    
    cout << "Introducir la temperatura del sistema en Kelvin:";
    cin >> temp;

    cout << "Introducir el número de pasos Montecarlo:";
    cin >> pasos;

    cout << "¿Configuración inicial aleatoria (1) o no aleatoria(0)?";
    cin >> aleatorio;

    pasos=pasos*dim*dim;

    fichero.open("datos_ising.dat");

    short int matriz[dim+2][dim+2];

    extern gsl_rng *tau; //Puntero al estado del número aleatorio
    int semilla=135254; //Semilla del generador de números aleatorios

    tau=gsl_rng_alloc(gsl_rng_taus); //Inicializamos el puntero
    gsl_rng_set(tau,semilla); //Inicializamos la semilla
    
    /* Generamos una configuración aleatoria o no aleatoria de la matriz de espines */
    for ( i = 1; i <= dim; i++)
    {
        for ( j = 1; j <= dim; j++)
        {
            random=gsl_rng_uniform_int(tau,2);
            if (random==1 && aleatorio==1)
            {
                matriz[i][j]=-1;
            }
            else matriz[i][j]=1;
        }
        
    }
    
    /* Aplicamos las condiciones periódicas a la matriz de espines for columnas y por filas */

    for ( i = 1; i <= dim; i++)
    {
        matriz[i][0]=matriz[i][dim];
        matriz[i][dim+1]=matriz[i][1];
    }

    for ( j = 1; j <= dim; j++)
    {
        matriz[0][j]=matriz[dim][j];
        matriz[dim+1][j]=matriz[1][j];
    }
    
    /* Copiamos la matriz en el fichero */

    for ( i = 1; i <= dim; i++)
    {
        for ( j = 1; j <= dim-1; j++)
        {
            fichero << matriz[i][j] << ',';
            cout << matriz[i][j] << ',';
        }
        fichero << matriz[i][dim+1] << '\n';
        cout << matriz[i][dim+1] << '\n';
        
    }
    fichero << '\n';


    /* Hacemos los pasos Montecarlo, en los que se hace el intento de cambio de los espines */
    for ( k = 1; k <= pasos; k++)
    {
        i=gsl_rng_uniform_int(tau,dim+1);
        j=gsl_rng_uniform_int(tau,dim+1);
        energia=2*matriz[i][j]*(matriz[i+1][j]+matriz[i-1][j]+matriz[i][j+1]+matriz[i][j-1]);
        p=exp(-energia/temp);
        if (p > 1.0)
        {
            p=1.0;
        }
        
        aux=gsl_rng_uniform(tau);
        if (aux < p)
        {
            matriz[i][j]=(-1)*matriz[i][j];
            /* Si la casilla está en el borde de la matriz, aplicamos las condiciones de periodicidad */
            if (i == 1)
            {
                matriz[dim+1][j]=(-1)*matriz[dim+1][j];
            }
            if (i == dim)
            {
                matriz[0][j]=(-1)*matriz[0][j];
            }
            if (j == 1)
            {
                matriz[i][dim+1]=(-1)*matriz[i][dim+1];
            }
            if (j == dim)
            {
                matriz[i][0]=(-1)*matriz[i][0];
            }   
        }
        /* Escribimos la matriz en el fichero solo en cada paso Montecarlo*/
        if (k%(dim*dim) == 1)
        {
            for ( i = 1; i <= dim; i++)
            {
                for ( j = 1; j <= dim-1; j++)
                {
                    fichero << matriz[i][j] << ',';
                    cout << matriz[i][j] << ',';
                }
                fichero << matriz[i][dim] << endl;
                cout << matriz[i][dim] << endl; 
            }
            fichero << endl;
            
        }
        
    }
    
    fichero.close();
    return 0;
}