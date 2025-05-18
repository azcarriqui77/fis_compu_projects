/* PROGRAMA QUE CALCULA LAS POSICIONES Y VELOCIDADES EN VARIOS INSTANTES DE TIEMPO
 DE LOS PLANETAS DEL SISTEMA SOLAR (Mercurio, Venus, ..., Neptuno) EN UN PLANO BIDIMENSIONAL 
 SEGÚN UNOS PARÁMETROS INICIALES DE POSICIÓN, MASA Y VELOCIDAD PARA CADA PLANETA Y LA MASA DEL SOL */

 #include <iostream>
 #include <cmath>
 #include <fstream>

 using namespace std;

 int main (void)
 {
    double m[9], x[9], y[9], vx[9], vy[9], ax[9], ay[9], wx[9], wy[9];
    double t,h,energia,l;
    double periodo[8];
    double Ms, c, G;
    int i,j,k;
    double modulo, epsilon;
    ofstream fichero_planetas;
    ofstream fichero_energia;
    ofstream fichero_periodos;

    Ms=1.99e+30;
    c=1.496e+11;
    G=6.67e-11;
    epsilon=0.1;

    /* Definimos las masas de los planetas, sus posiciones iniciales y sus velocidades iniciales*/
    /* Inicialmente consideramos los planetas sobre el eje X, luego su velocidad inicial solo tendrá
    componente en el eje Y*/
    /* Hacemos ya el reescalamiento*/
    /* El índice 0 de los vectores está reservada para el Sol y el resto de índices para los planetas en orden creciente de órbita (Mercurio-1, Venus-2...)*/
    /* Perturbo ligeramente la posición inicial de un planeta y cambio de signo la velocidad inicial de otro.*/
    m[0]=1; x[0]=0; y[0]=0; vx[0]=0; vy[0]=0; ax[0]=0; ay[0]=0;
    m[1]=0.330e+24/Ms; x[1]=57.9e+9/c; y[1]=0; vx[1]=0; vy[1]=47.9e+3/sqrt(G*Ms/c); ax[1]=0; ay[1]=0; periodo[0]=0;
    m[2]=4.87e+24/Ms; x[2]=108.2e+9/c; y[2]=0; vx[2]=0; vy[2]=35.0e+3/sqrt(G*Ms/c); ax[2]=0; ay[2]=0; periodo[1]=0;
    m[3]=5.97e+24/Ms; x[3]=149.6e+9/c; y[3]=0; vx[3]=0; vy[3]=29.8e+3/sqrt(G*Ms/c); ax[3]=0; ay[3]=0; periodo[2]=0;
    m[4]=0.642e+24/Ms; x[4]=227.9e+9/c; y[4]=0; vx[4]=0; vy[4]=24.1e+3/sqrt(G*Ms/c); ax[4]=0; ay[4]=0; periodo[3]=0;
    m[5]=1899e+24/Ms; x[5]=778.6e+9/c; y[5]=0; vx[5]=0; vy[5]=13.1e+3/sqrt(G*Ms/c); ax[5]=0; ay[5]=0; periodo[4]=0;
    m[6]=568e+24/Ms; x[6]=1433.5e+9/c; y[6]=0; vx[6]=0; vy[6]=9.7e+3/sqrt(G*Ms/c); ax[6]=0; ay[6]=0; periodo[5]=0;
    m[7]=86.8e+24/Ms; x[7]=2872.5e+9/c; y[7]=0; vx[7]=0; vy[7]=6.8e+3/sqrt(G*Ms/c); ax[7]=0; ay[7]=0; periodo[6]=0;
    m[8]=102e+24/Ms; x[8]=4495e+9/c; y[8]=0; vx[8]=0; vy[8]=5.4e+3/sqrt(G*Ms/c); ax[8]=0; ay[8]=0; periodo[7]=0;

   // x[1]=x[1]+epsilon; // vy[5]=-vy[5];

    /* Para comprobar que todo está correcto, vamos a escribir estos vectores iniciales en pantalla*/
   cout << "Masa  Posición x  Posición y  Velocidad x Velocidad y" << endl;
   for ( i = 0; i <= 8; i++)
   {
      cout << m[i] << '\t' << x[i] << '\t' << y[i] << '\t' << vx[i] << '\t' << vy[i] << endl;
   }

    /* Abrimos el fichero donde vamos a escribir los valores de las posiciones de los planetas
    según el formato descrito en el script de la animación gif proporcionado en clase*/
    /* Abrimos otro fichero para escribir los valores de energía y momento angular en función del tiempo y el último para escribir los períodos de los planetas*/
    fichero_planetas.open("planets_data.dat");
    fichero_energia.open("energy_time.dat");
    fichero_periodos.open("periodos.dat");

    /* Hacemos un formato bonito para los documentos*/
    fichero_energia << "Tiempo astronómico" << '\t' << "Energía" << '\t' << "Momento angular" << endl;
    fichero_periodos << "Planetas" << '\t' << "Tiempo astronómico" << '\t' << "Tiempo (días)" << endl;

   /* Iniciamos el proceso iterativo de cálculo de las aceleraciones, posiciones, vector auxiliar, nuevas aceleraciones
   y velocidades según el algoritmo de Verlet.*/
   h=0.01; /* En el tiempo gravitatorio 1 día equivale a 0.17 segundos aproximadamente. */
   t=0;

   /* En primer lugar el cálculo de la aceleración a t=0.*/
   for ( i = 1; i <= 8; i++)
   {
      for ( j = 0; j <= 8; j++)
      {
         if (j!=i)
         {
            modulo=pow(x[i]-x[j],2)+pow(y[i]-y[j],2);
            modulo=pow(modulo,1.5);
            ax[i]=ax[i]-m[j]*(x[i]-x[j])/modulo;
            ay[i]=ay[i]-m[j]*(y[i]-y[j])/modulo;
         }   
      }
   }

   /* ALGORTIMO DE VERLET */

   for ( k = 0; k < 110000; k++)
   {
      /* Escribimos las posiciones de los planetas y la energía del sistema en los ficheros correspondientes antes de pasar al siguiente tiempo.*/
      for ( i = 0; i <= 8; i++)
      {
         fichero_planetas << x[i] << ',' << y[i] << endl;
      }
      fichero_planetas << endl;

      energia=0; l=0;
      for ( i = 1; i <=8; i++)
      {
         energia=energia+m[i]*(vx[i]*vx[i]+vy[i]*vy[i])/2;
         l=l+m[i]*(x[i]*vy[i]-y[i]*vx[i]);
         for ( j = 0; j <=8; j++)
         {
            if (j!=i)
            {
               energia=energia-m[i]*m[j]*pow(pow(x[i]-x[j],2)+pow(y[i]-y[j],2),-0.5);
            }
         } 
      }
      fichero_energia << t << '\t' << energia << '\t' << l << endl;

      /* Para calcular los períodos de los planetas, uso tres condiciones: [y muy próximo a cero] y [x positivo] y [k mayor que un cierto número de 
      iteraciones (para evitar el instante inicial)]. En ese momento, escribe en el fichero de periodos el tiempo que ha tardado cada planeta en dar una
      órbita.*/

      for ( i = 0; i <=7 ; i++)
      {
         if ((periodo[i]==0) && (abs(y[i+1]) <= 0.01) && (x[i+1] >= 0) && (k>100))
         {
            fichero_periodos << "Planeta número " << i+1 << '\t' << t << '\t' << t/(sqrt(G*Ms/pow(c,3))*3600*24) << endl;
            periodo[i]=1;
         }
         
      }
      

      /* Actualizamos el tiempo y calculamos las nuevas posiciones y velocidades de los planetas. */
      t=t+h;

      for (i = 1; i <= 8; i++)
      {
         x[i]=x[i]+h*vx[i]+h*h*ax[i]/2;
         y[i]=y[i]+h*vy[i]+h*h*ay[i]/2;
      }
      for (i = 1; i <=8; i++)
      {
         wx[i]=vx[i]+h*ax[i]/2;
         wy[i]=vy[i]+h*ay[i]/2;
         ax[i]=0; ay[i]=0;
      }

      for ( i = 1; i <=8; i++)
      {
         for ( j = 0; j <= 8; j++)
         {
            if (j!=i)
            {
               modulo=pow(x[i]-x[j],2)+pow(y[i]-y[j],2);
               modulo=pow(modulo,1.5);
               ax[i]=ax[i]-m[j]*(x[i]-x[j])/modulo;
               ay[i]=ay[i]-m[j]*(y[i]-y[j])/modulo;
            }
         }
      }
      for ( i = 1; i <= 8; i++)
      {
         vx[i]=wx[i]+h*ax[i]/2;
         vy[i]=wy[i]+h*ay[i]/2;
      }
   }

   fichero_planetas.close();
   fichero_energia.close();
   fichero_periodos.close();
   return 0;
 }
