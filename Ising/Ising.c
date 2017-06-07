#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <math.h>

# define J 0.0
# define B 1.0
# define k 1.0

# define RED   "\x1b[31m"
# define BLUE  "\x1b[34m"
# define RESET "\x1b[0m"

float flip(int *lattice, int n, float T, int *M);          // Me genera la nueva red segun metropolis monte-carlo para un cambio de spin.
float metropolis(int *lattice, int n, float T, int *M);    // Me genra lo pasos de Metropolis.
float energy_change(int *lattice, int n, int i, int j);    // Me da el cambio de enrgia al cambiar un spin.
float energy(int *lattice, int n);                         // Me da la energia totaldel sistema.
float valor_medio(int *c, int b);                          // calcula el valor medio del vector c de tamaño b.
float varianza(int *c, int b);                             // Calcula la varianza del vector c de tamaño b.
float correlation(int *c, int b, int tau);                 // Me da la correlacion del vector c y tamaño b en salto tau.
int magnetization(int *lattice, int n);	                   // Me da la magnetizcion del sistema.
int neighbors(int *lattice, int n, int i, int j);		   // Me da el valor de la red para la posicon i-j.
int fill_lattice(int *lattice, int n, float p);            // Crea la red de n*n con valores -1 / 1 con probabilidad p.
int print_lattice(int *lattice, int n);					   // Imprime la red en pantalla.
int pick_site(int *lattice, int n);                        // Me genera un numero rando entre 0 y n-1.


int main(){
	
	/* Parametros */
	float p = 0.5;              // Probabilidad de llenado.
	int   n = 128;		        // Tamaño red.
	int   term  = 100000;		// Pasos termlizacion.
	int   niter = term*100;      // Pasos luego de termalizar.


	/* *************************************/ // Termalizacion //
/*
	float T = 2.26; // Temperatura. Tc = 2.26
	int   i;
	FILE* fichero;
    fichero = fopen("E_termalizacion.txt", "wt");
    FILE* fichero2;
    fichero2 = fopen("M_termalizacion.txt", "wt");
	
	// Genero la red.
	srand(time(NULL));
	int *lattice;
	lattice=(int *)malloc(n*n*sizeof(int));
	fill_lattice(lattice,n,p);

	// Monte Carlo Termalización.
	float Et = energy(lattice,n); 
	int M  = magnetization(lattice, n);
	print_lattice(lattice,n);
	for (i=0; i<term; i++){
		fprintf(fichero, "%i    %f\n",i, Et);
		fprintf(fichero2,"%i    %i\n",i, M);
    	float Dt = metropolis(lattice, n, T, &M);
    	Et = Et+Dt; // Nueva energia.
	}
	int m = magnetization(lattice, n);
	float e = energy(lattice,n);
	printf("Magnetizacion = %i /  %i\n", m, M);
	printf("Energia       = %f /  %f\n", e, Et);
	print_lattice(lattice,n);
	fclose(fichero);
	free(lattice);


	/* *************************************/ // M/E vs T //
/*
	int puntos = 100, j, i;
	float T_vec[puntos], T;
	FILE* fichero3;
    fichero3 = fopen("MvsT.txt", "wt");

    FILE* fichero4;
    fichero4 = fopen("EvsT.txt", "wt");

    // Genero el vector de temperaturas.
	for (i=0; i<puntos; i++){ 
		T_vec[i] = puntos - i;
	}
	for (i=0; i<puntos; i++){
		T_vec[i] = T_vec[i]/puntos*4; // t entre 0 y 4.
	}

	// Genero la red.
	srand(time(NULL)); 
	int *lattice;
	lattice=(int *)malloc(n*n*sizeof(int));
	fill_lattice(lattice,n,p);

	int *m;
	m=(int *)malloc(niter*sizeof(int));

	int *e;
	e=(int *)malloc(niter*sizeof(int));

	for (j=0; j<puntos; j++){
		printf("cargando %i porciento\n", (int)(100*j/puntos));
		T = T_vec[j];
		int M = magnetization(lattice,n);
		float Et = energy(lattice,n);
		for (i=0; i<term; i++){
			float Dt = metropolis(lattice, n, T, &M);
    		Et = Et+Dt; // Nueva energia.
		}
		for (i=0; i<niter; i++){
			float Dt = metropolis(lattice, n, T, &M);
    		Et = Et+Dt; // Nueva energia.
			*(m+i) = M;
			*(e+i) = Et;
		}

		double mm = 0;
		double ee = 0;
		for (i=0; i<niter; i++){
			mm = mm + (double)*(m+i);
			ee = ee + (double)*(e+i);
		}	
		mm = mm/(double)niter;
		ee = ee/(double)niter;

		fprintf(fichero3,"%f    %f\n", T, mm);	
		fprintf(fichero4,"%f    %f\n", T, ee);	

	}

	fclose(fichero3);
	fclose(fichero4);
	free(lattice);
	free(m);
	free(e);


	/* *************************************/ // Correlacion //
/*
	float T = 4.5; // Temperatura.
    int i, tau;

	FILE* fichero4;
	FILE* fichero5;
    fichero4 = fopen("Corelation_E.txt", "wt");
    fichero5 = fopen("Corelation_M.txt", "wt");

    // Genero la red.
    srand(time(NULL));
	int *lattice;
	lattice=(int *)malloc(n*n*sizeof(int));
	fill_lattice(lattice,n,p);

	// Termalizo.
	int   Mt = magnetization(lattice,n);
	for (i=0; i<term; i++){
    	metropolis(lattice, n, T, &Mt);
	}

	int *m;
	m=(int *)malloc(niter*sizeof(int));
	int *e;
	e=(int *)malloc(niter*sizeof(int));

	int M = magnetization(lattice,n);
	float Et = energy(lattice,n);
	for (i=0; i<term; i++){
		float Dt = metropolis(lattice, n, T, &M);
		Et = Et+Dt; // Nueva energia.
	}
	for (i=0; i<niter; i++){
		float Dt = metropolis(lattice, n, T, &M);
		Et = Et+Dt; // Nueva energia.
		*(m+i) = M;
		*(e+i) = (int)Et;
	}
	print_lattice(lattice,n);


	// Correlacion
	for (tau=0; tau<(int)(niter/2); tau = tau + (int)((niter/2)/100)){
		printf("cargando correlacion %i porciento\n", (int)(100*tau/(niter/2)));
		float mm = correlation(m,niter,tau);
		float ee = correlation(e,niter,tau);
		fprintf(fichero4,"%i    %f\n",tau, ee);
		fprintf(fichero5,"%i    %f\n",tau, mm);

	}

	fclose(fichero4);
	fclose(fichero5);
	free(lattice);
	free(e);
	free(m);


	/* *************************************/ // TCor vs T //
/*
    int i, j, tau;
    int puntos = 100;
    float T_vec[puntos], T;
	FILE* fichero4;
    fichero4 = fopen("Cor_E_vs_T.txt", "wt");

    // Genero el vector de temperaturas.
	for (i=0; i<puntos; i++){ 
		T_vec[i] = puntos - i;
	}
	for (i=0; i<puntos; i++){
		T_vec[i] = T_vec[i]/puntos*2+1.5; // t entre 1,5 y 3,5.
	}

    // Genero la red.
    srand(time(NULL));
	int *lattice;
	lattice=(int *)malloc(n*n*sizeof(int));
	fill_lattice(lattice,n,p);

	int *m;
	m=(int *)malloc(niter*sizeof(int));
	int *e;
	e=(int *)malloc(niter*sizeof(int));

	for (j=0; j<puntos; j++){

		printf("cargando %i porciento\n", (int)(100*j/puntos));
		T = T_vec[j];

		// Termalizo.
		int   Mt = magnetization(lattice,n);
		for (i=0; i<term; i++){
	    	metropolis(lattice, n, T, &Mt);
		}

		int M = magnetization(lattice,n);
		float Et = energy(lattice,n);
		for (i=0; i<term; i++){
			float Dt = metropolis(lattice, n, T, &M);
			Et = Et+Dt; // Nueva energia.
		}
		for (i=0; i<niter; i++){
			float Dt = metropolis(lattice, n, T, &M);
			Et = Et+Dt; // Nueva energia.
			*(m+i) = M;
			*(e+i) = (int)Et;
		}


		// Correlacion
		for (tau=0; tau<(int)(niter); tau = tau + (int)((niter)/1000)){
			//float mm = correlation(m,niter,tau);
			float ee = correlation(e,niter,tau);
			if (ee<=0){
				fprintf(fichero4,"%f    %i\n",T, tau);
				break;
			}

		}
	}

	fclose(fichero4);
	free(lattice);
	free(e);
	free(m);


	/* *************************************/ // Hist M //
/*
	int salto = 5000;
	int j, i, puntos = 500000;
	float T = 2.26;
	FILE* fichero6;
    fichero6 = fopen("Hist_M.txt", "wt");
    FILE* fichero7;
    fichero7 = fopen("time_M.txt", "wt");

	// Genero la red.
	srand(time(NULL)); 
	int *lattice;
	lattice=(int *)malloc(n*n*sizeof(int));

	int *m;
	m=(int *)malloc(puntos*sizeof(int));

	fill_lattice(lattice,n,p);
	int M = magnetization(lattice,n);
	for (i=0; i<term; i++){
		metropolis(lattice, n, T, &M);
	}
	for (i=0; i<puntos; i++){
		printf("cargando %i porciento\n", (int)(100*i/puntos));
		for (j=0; j<salto; j++){
			metropolis(lattice, n, T, &M);
		}
		*(m+i) = M;
	}

	int *hist;
	hist=(int *)malloc((2*n*n+1)*sizeof(int));
	for (i=0; i<(2*n*n+1); i++){
		*(hist+i) = 0;
	}

	for (i=0; i<puntos; i++){
		*(hist + n*n + *(m+i)) = *(hist + n*n + *(m+i)) + 1;
	}

	for (i=0; i<(2*n*n+1); i++){
		fprintf(fichero6,"%i    %i\n",i-n*n, *(hist+i));
		*(hist+i) = 0;
	}

	for (i=0; i<puntos; i++){
		fprintf(fichero7,"%i    %i\n",i, *(m+i));
	}

	fclose(fichero6);
	fclose(fichero7);
	free(lattice);
	free(hist);
	free(m);


	/* *************************************/ // Cv //

	int puntos = 100, j, i, c;
	int salto  = 10*n*n, var_p = 5000;
	float T_vec[puntos], T;
	FILE* fichero3;
    fichero3 = fopen("varMvsT_J=0.txt", "wt");

    FILE* fichero4;
    fichero4 = fopen("varEvsT_J=0.txt", "wt");

    // Genero el vector de temperaturas.
	for (i=0; i<puntos; i++){ 
		T_vec[i] = puntos - i;
	}
	for (i=0; i<puntos; i++){
		T_vec[i] = T_vec[i]/puntos*1+2; // t entre 2 y 3.
	}

	// Genero la red.
	srand(time(NULL)); 
	int *lattice;
	lattice=(int *)malloc(n*n*sizeof(int));
	fill_lattice(lattice,n,p);

	int *m;
	m=(int *)malloc(var_p*sizeof(int));

	int *e;
	e=(int *)malloc(var_p*sizeof(int));

	for (j=0; j<puntos; j++){
		printf("cargando %i porciento\n", (int)(100*j/puntos));
		T = T_vec[j];
		int M = magnetization(lattice,n);
		float Et = energy(lattice,n);
		for (i=0; i<term; i++){
			float Dt = metropolis(lattice, n, T, &M);
    		Et = Et+Dt; // Nueva energia.
		}
		for (i=0; i<var_p; i++){
			for (c=0; c<salto; c++){
				float Dt = metropolis(lattice, n, T, &M);
	    		Et = Et+Dt; // Nueva energia.
	    	}
			*(m+i) = M;
			*(e+i) = Et;
		}

		float mm = varianza(m,var_p);
		float ee = varianza(e,var_p);

		fprintf(fichero3,"%f    %f\n", T, mm);	
		fprintf(fichero4,"%f    %f\n", T, ee);	

	}

	fclose(fichero3);
	fclose(fichero4);
	free(lattice);
	free(m);
	free(e);


	/* *************************************/ // Fin //

	return 0;

}


/* Funciones auxiliares */


int fill_lattice(int *lattice, int n, float p){ // Crea la red de n*n con valores -1 / 1 con probabilidad p.

	float R;
	int i;	

	for (i=0; i<n*n; i++){
		R = (float)rand()/(float)RAND_MAX;
		if (R<=p){
			*(lattice+i) =  1;
		}else{
			*(lattice+i) = -1;
		}
	}

	return 0;

}


int print_lattice(int *lattice, int n){ // Imprime la red en pantalla.

	int i;	

	if(n<=40){  // Porque sino queda feo en pantalla si es mas grande.
		printf("\n");
		for (i=0; i<n*n; i++){
	 		if ((i+1)%n == 0){
	 			if(*(lattice+i)==-1){
	 				printf(RED "%2i" RESET "\n", *(lattice+i));
	 			}else{
	 				printf(BLUE "%2i" RESET "\n", *(lattice+i));
	 			}			
	 		}else{
	 			if(*(lattice+i)==-1){
	 				printf(RED "%2i" RESET, *(lattice+i));
	 			}else{
	 				printf(BLUE "%2i" RESET, *(lattice+i));
	 			}	
	 		}
		}
	}

	printf("\n");

	return 0;

}


int pick_site(int *lattice, int n){ // Me genera un numero rando entre 0 y n-1.

	float rnd;
	int i;
	while(1){
		rnd = (float)rand()/(float)RAND_MAX; 
		i = (int)(n*rnd);
		if (i!=n){
			break;
		}
	}
	
	return i;

}


float flip(int *lattice, int n, float T, int *M){ // Me genera la nueva red segun metropolis monte-carlo para un cambio de spin.

	int i = pick_site(lattice,n); //selcciono un spin.
	int j = pick_site(lattice,n); //selcciono un spin.
	float R = (float)rand()/(float)RAND_MAX;

	float dt = energy_change(lattice,n,i,j);

	float p = exp(-dt/(k*T)); if(p>=1){p = 1;}

	float Dt = 0;
	if (R<=p){
		*(lattice+i*n+j) = *(lattice+i*n+j) *-1;
		Dt = dt;
		*M = *M + *(lattice+i*n+j) *2;
	}

	return Dt;

}


float energy_change(int *lattice, int n, int i, int j){ // Me da el cambio de enrgia al cambiar un spin.

	int c,u,d,l,r;

	c = neighbors(lattice,n,i,j);
	u = neighbors(lattice,n,i-1,j);
	d = neighbors(lattice,n,i+1,j);
	l = neighbors(lattice,n,i,j-1);
	r = neighbors(lattice,n,i,j+1);

	// Segundos vecinos.
	int ul,ur,dl,dr;
	ul = neighbors(lattice,n,i-1,j-1);
	ur = neighbors(lattice,n,i-1,j+1);
	dl = neighbors(lattice,n,i+1,j-1);
	dr = neighbors(lattice,n,i+1,j+1);
	float E0 = -J*c*(u+d+l+r)+(J/sqrt(2))*c*(ul+ur+dl+dr)-B*c;
	float E1 =  J*c*(u+d+l+r)-(J/sqrt(2))*c*(ul+ur+dl+dr)+B*c;	

	//float E0 = -J*c*(u+d+l+r)-B*c;
	//float E1 =  J*c*(u+d+l+r)+B*c;

	float dt = E1-E0;

	return dt;

}


float energy(int *lattice, int n){  // Me da la energia total del sistema.

	int i,j;
	int c,d,r;
	float E = 0;
	for(i=0; i<n; i++){
		for(j=0; j<n; j++){	
	
		c = neighbors(lattice,n,i,j);
		d = neighbors(lattice,n,i+1,j);
		r = neighbors(lattice,n,i,j+1);	

		// Segundos vecinos.
		int ur,dr;
		ur = neighbors(lattice,n,i-1,j+1);
		dr = neighbors(lattice,n,i+1,j+1);
		float E0 = -J*c*(d+r)+(J/sqrt(2))*c*(ur+dr)-B*c;

		//float E0 = -J*c*(d+r)-B*c;
		E = E + E0;

		}
	}

	return E;

}

	
float metropolis(int *lattice, int n, float T, int *M){ // Me genra lo pasos de Metropolis.

	float Dt = flip(lattice,n,T,M);

	return Dt;

}


int neighbors(int *lattice, int n, int i, int j){ // Me da el valor de la red para la posicon i-j.

	int f = (i+n)%n;
	int c = (j+n)%n;
	int v = *(lattice+f*n+c);

	return v;

}


int magnetization(int *lattice, int n){ // Me da la magnetizcion del sistema.  

	int i,j;
	int M = 0;

	for(i=0; i<n; i++){
		for(j=0; j<n; j++){	
			M = M + *(lattice+i*n+j);
		}
	}

	return M;

}


float valor_medio(int *c, int b){ // Calcula el valor medio del vector c de tamaño b.
	
	float suma=0, media;
	int i;
	for (i=0; i<b; i++){
		suma = suma + (float)*(c+i);
	}
	media = suma/(float)b;

	return media;

}


float varianza(int *c, int b){ // Calcula la varianza del vector c de tamaño b.
	
	float suma=0, media, var=0;
	int i;

	for (i=0; i<b; i++){
		suma = suma + (float)*(c+i);
	}
	media = suma/(float)b;


	for (i=0; i<b; i++){
		var = var + (float)(((float)*(c+i)-media)*((float)*(c+i)-media));
	}
	var = var/(float)b;

	return var;

}


float correlation(int *c, int b, int tau){ // Me da la correlacion del vector c y tamaño b en salto tau.

	int i;

	float mean = valor_medio(c,b);
	float var = varianza(c,b);
	float R = 0;
	for(i=0; i<b-tau; i++){
		R = R + (float)(((float)*(c+i)-mean)*((float)*(c+i+tau)-mean));
	}
	R = R/(float)(b-tau);
	R = R/var;

	return R;

}