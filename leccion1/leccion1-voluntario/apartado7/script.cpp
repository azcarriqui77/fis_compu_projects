/* SIMULACIÓN DE LA DINÁMICA MOLECULAR DE UN GAS CON UN POTENCIAL DE LENNARD-JONES */
/* EJERCICIO VOLUNTARIO 1 LECCIÓN 1 SISTEMA SOLAR */

#include <iostream>
#include <cmath>
#include <fstream>
#include <stdlib.h>
#include <time.h>


using namespace std;

/* Definimos unas funciones necesarias en nuestro programa */
double distancia (double x1, double y1, double x2, double y2, double dif_coord[]);
double posicion (double x);
double fuerza (double x1, double y1, double x2, double y2, double dif_coord[]);

int main (void)
{
    int n; n=16;
    double x[n], y[n], vx[n], vy[n], ax[n], ay[n], wx[n], wy[n], dif_coord[2];
    double t,h,energia,T,V,f,r,r6,temp,aux;
    int i,j,k;
    ofstream fichero_desviacion;
    ofstream fichero_temperatura;

    srand(time(NULL));

    fichero_desviacion.open("/Users/alzorrcarri/Documents/cphys/leccion1-voluntario/apartado7/desviacion_conincr.dat");
    fichero_temperatura.open("/Users/alzorrcarri/Documents/cphys/leccion1-voluntario/apartado7/temperatura.dat");


    /* Damos valores iniciales a las posiciones dentro de la cuadrícula 4x4 */

    cout << "X" << '\t' << "Y" << '\t' << "V_X" << '\t' << "V_Y" << endl;
   for ( i = 0; i < 4; i++)
   {
      for ( j = 0; j < 4; j++)
      {
         x[4*i+j]=i+0.5;
         y[4*i+j]=j+0.5;
      }
      
   }
        
      /* Consideramos velocidades iniciales nulas. */
   for ( i = 0; i <= n-1; i++)
    {
      vx[i] = 0;
      vy[i] = 0;
      cout << x[i] << '\t' << y[i] << '\t' << vx[i] << '\t' << vy[i] << endl; // Vemos si todo es correcto.
    }
   
    /* Iniciamos el proceso iterativo de cálculo de las aceleraciones, posiciones, vector auxiliar, nuevas aceleraciones
    y velocidades según el algoritmo de Verlet. */
   h=0.002;
   t=0;
   for ( i = 0; i <= n-1; i++)
    {
        ax[i]=0;
        ay[i]=0;
    }

   /* En primer lugar el cálculo de la aceleración a t=0. */
   for ( i = 0; i <= n-1; i++)
   {
      for ( j = 0; j <= n-1; j++)
      {
         if (j!=i)
         {
            f=fuerza(x[i],y[i],x[j],y[j],dif_coord);
            ax[i]=ax[i]+f*dif_coord[0];
            ay[i]=ay[i]+f*dif_coord[1];
         }   
      }
   }

   /* ALGORTIMO DE VERLET */
   temp=0;
   for ( k = 0; k <= 150000; k++)
   {
      aux=0;
      for ( i = 0; i <= n-1; i++)
      {
         aux=aux+vx[i]*vx[i]+vy[i]*vy[i];
      }
      aux=aux/(2*16.0);
      temp=temp+aux;

      /* Estudiamos la separación entre las partículas 5 y 6 y la temperatura en un instante de tiempo. */
      fichero_desviacion << t << '\t' << pow(distancia(x[5],y[5],x[6],y[6],dif_coord),2) << '\n';

      if (k%100==0)
      {
         fichero_temperatura << t << '\t' << temp/100.0 << '\n';
         temp=0;
      }
      
      

      /* Actualizamos el tiempo y calculamos las nuevas posiciones y velocidades de los planetas. */
      t=t+h;

      /* Calentamos el sistema en ciertos instantes de tiempo. */
      if (k*h==60.0 || k*h==120.0 || k*h==180.0 || k*h==240.0 || k*h==300.0)
      {
         for ( i = 0; i < n-1; i++)
         {
            vx[i]=1.2*vx[i];
            vy[i]=1.2*vy[i];
         }
         
      }
      

      for (i = 0; i <= n-1; i++)
      {
         x[i]=x[i]+h*vx[i]+h*h*ax[i]/2;
         x[i]=posicion(x[i]);
         y[i]=y[i]+h*vy[i]+h*h*ay[i]/2;
         y[i]=posicion(y[i]);
      }
      for (i = 0; i <= n-1; i++)
      {
         wx[i]=vx[i]+h*ax[i]/2;
         wy[i]=vy[i]+h*ay[i]/2;
         ax[i]=0; ay[i]=0;
      }

      for ( i = 0; i <= n-1; i++)
      {
         for ( j = 0; j <= n-1; j++)
         {
            if (j!=i)
            {
               f=fuerza(x[i],y[i],x[j],y[j],dif_coord);
               ax[i]=ax[i]+f*dif_coord[0];
               ay[i]=ay[i]+f*dif_coord[1];
            }
         }
      }
      for ( i = 0; i <= n-1; i++)
      {
         vx[i]=wx[i]+h*ax[i]/2;
         vy[i]=wy[i]+h*ay[i]/2;
      }
   }

    fichero_desviacion.close();
    fichero_temperatura.close();
    return 0;
}

/* Definición de las funciones empleadas */
/* La función distancia calcula la distancia entre dos partículas teniendo en cuenta las condiciones de contorno periódicas. */
/* Usamos que los arrays siempre se pasan por referencia en las funciones. */
double distancia (double x1, double y1, double x2, double y2, double dif_coord[])
{
    double difx, dify, difx2, dify2;
    double d1,d2,d3,d4;
    double r;

    difx=x1-x2;
    dify=y1-y2;
    if (x1<=x2)
    {
      difx2=x1-x2+4;
    } else difx2=x1-x2-4;
    if (y1<=y2)
    {
      dify2=y1-y2+4;
    } else dify2=y1-y2-4;

    // Calculamos las cuatro posibles distancias y escogemos la mínima.
    d1=difx*difx+dify*dify;
    d2=difx*difx+dify2*dify2;
    d3=difx2*difx2+dify*dify;
    d4=difx2*difx2+dify2*dify2;

    if (d1<=d2 && d1<=d3 && d1<=d4)
    {
      r=sqrt(d1);
      dif_coord[0]=difx;
      dif_coord[1]=dify;
    }
    if (d2<=d1 && d2<=d3 && d2<=d4)
    {
      r=sqrt(d2);
      dif_coord[0]=difx;
      dif_coord[1]=dify2;
    }
    if (d3<=d1 && d3<=d2 && d3<=d4)
    {
      r=sqrt(d3);
      dif_coord[0]=difx2;
      dif_coord[1]=dify;
    }
    if (d4<=d1 && d4<=d2 && d4<=d3)
    {
      r=sqrt(d4);
      dif_coord[0]=difx2;
      dif_coord[1]=dify2;
    }
    return r;   
    
}

/* La función posición comprubea que todas las partículas queden dentro de la caja, y aquéllas que están fuera, las mete dentro
siguiendo las condiciones de contorno periódicas. */
double posicion (double x)
{
    if (x>4.0 || x<0.0)
    {
        x=x-4*floor(x/4.0);
    }
    return x;
}

/* La función fuerza devuelve el valor del módulo de la fuerza entre dos partículas según la distancia que las separa.*/
double fuerza (double x1, double y1, double x2, double y2, double dif_coord[])
{
   double r,r6,f;
   r=distancia(x1,y1,x2,y2,dif_coord);
   r6=pow(r,6);
   f=24*(2.0/r6-1)/(r6*r*r);
   return f;
}
