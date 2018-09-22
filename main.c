#include <stdio.h>
#include <stdlib.h>
#include <math.h>

 double x;
 double alpha;
 double Fx[4];
 double A;
 double B;
 double C;
 double D;
 double vec_aux[4] = {0,0,0,0};
 double PolDum[3] = {0,0,0};
 double PolDdois[2] = {0,0};
 double raizesPol[2] = {0,0};
 double intervalos[5] = {0,0,0,0,80};
 double raizlinlin;
 double p[4];
 double q[3];
 int termina = 0;
 double a;
void avaliaA()
{
    A = -0.001;
}

void avaliaB(double alpha)
{
    B = 0.105 - 5/(3014.01*pow(cos(alpha),2));
}

void avaliaC(double alpha)
{
    C = -3 + tan(alpha) - 500/(3014.01*pow(cos(alpha),2));
}

void avaliaD(double alpha)
{
    D = 50*tan(alpha) - 12500/(3014.01*pow(cos(alpha),2));
}


void avaliaCoef(double alpha)
{
    avaliaA();
    avaliaB(alpha);
    avaliaC(alpha);
    avaliaD(alpha);
}


/*
    assumimos que esta função recebe os coeficiente de um polinômio de terceiro grau.
    isso significa que deverá receber um vetor de tamanho 4
*/
void coef_d1(double coefPol[])
{
    PolDum[0] = 3*coefPol[0];
    PolDum[1] = 2*coefPol[1];
    PolDum[2] = coefPol[2];
}

/*
    assumimos que esta função recebe os coeficiente de um polinômio de segundo grau.
    isso significa que deverá receber um vetor de tamanho 3
*/
void coef_d2(double coefPol[])
{
    PolDdois[0] = 6*coefPol[0];
    PolDdois[1] = 2*coefPol[1];
}

double termo(double a_im1,double a_i1, double z)
{
    double b = a_im1*z + a_i1;
    return b;
}

/*
LEMBRAR DE ADICIONAR X COMO ARGUMENTO PARA O ALGORITMO DE HORNER.
ISTO ESTAVA ATRAPALHANDO A RESOLVER OUTRAS COISAS.
SERÁ NECESSÁRIO FAZER ADAPTAÇÕES.
*/
double horner(double coef[], int tam, double z)
{
    vec_aux[0] = coef[0];
    int i;
    for (i=1; i<tam; i++)
    {
        vec_aux[i] = termo(vec_aux[i-1],coef[i],z);
    }
    return vec_aux[tam-1];
}

void RaizFlin()
{
    double A1 = PolDum[0];
    double B1 = PolDum[1];
    double C1 = PolDum[2];
    double delta = B1*B1 - 4*A1*C1;
    if (delta < 0)
    {
        /*printf("delta menor que zero");*/
        return;
    }
    raizesPol[0] = (-B1 + pow(delta,0.5))/(2*A1);
    raizesPol[1] = (-B1 - pow(delta,0.5))/(2*A1);
    if ((raizesPol[0] < 0) || raizesPol[0] > 80)
    {
        raizesPol[0] = 0;
    }
    if ((raizesPol[1] < 0) || raizesPol[1] > 80)
    {
        raizesPol[1] = 0;
    }
    return;
}

void RaizFlinlin()
{
    raizlinlin = -(PolDdois[1]/PolDdois[0]);
    if ((raizlinlin < 0) || raizlinlin > 80)
    {
        raizlinlin = 0;
    }
    return;
}

double alcMax(double alpha)
{
    return -50 + ((54.9*54.9)/10)*sin(2*alpha);
}

void InsertionSort(double original[], int lenght)
{
	int i, j;
	double atual;
	for (i = 1; i < lenght; i++)
	{
		atual = original[i];
		j = i - 1;
		while ((j >= 0) && (atual < original[j]))
		{
			original[j + 1] = original[j];
            j = j - 1;
		}
		original[j + 1] = atual;
	}

	return;
}

double NewtonIte(double pol[])
{
    double F = horner(pol,4,x);
    double Flin = horner(vec_aux,3,x);
    return (x - (F/Flin));
}

void Newton(double pol[])
{
    while ((horner(pol,4,x) > 0.0001) || horner(pol,4,x) < -0.0001)
    {
        /*printf("x eh: %f\n",x);*/
        x = NewtonIte(pol);
    }
    return;
}

double avaliaGalfa(double a)
{
    double G = ((100*tan(a)) - (16.5892/( cos(a)*cos(a) )) - 12.5);
    return G;
}

double avaliaGlin(double a)
{
    double Glin = ((100 - 33.1784*tan(a))/( (cos(a))*(cos(a)) ));
    return Glin;
}
double AlfaNewtonIte(double a)
{
    return (a - (avaliaGalfa(a)/avaliaGlin(a)));
}

void AlfaNewton()
{
    while ((avaliaGalfa(a) > 0.0001) || avaliaGalfa(a) < -0.0001)
    {
        /*printf("a eh: %f\n",a);*/
        a = AlfaNewtonIte(a);
    }
    return;
}

void resolveAlvo()
{
    /*
    Calularemos agora o valor de alfa tal que x =50 e y = 12.5 (alvo)
    */

    /*
    */

    /*
        a seguir estão definido os intervalos feitos para checar
    */
    double Ialfa[5] = {0,0.180475,1.06997,1.25044,1.5707};
    if (( ((avaliaGalfa(Ialfa[0]+0.0005))*(avaliaGalfa(Ialfa[1]-0.005))) < 0))
        {
            /* printf ("f(r0) = %f e f(r1) = %f\n",(avaliaGalfa(Ialfa[0]+0.0005)),(avaliaGalfa(Ialfa[1]-0.005)));*/
            /* printf ("f(r0)*f(r1) = %f\n", ((avaliaGalfa(Ialfa[0]+0.0005))*(avaliaGalfa(Ialfa[1]-0.005))) );*/
            /* printf("primeiro if de 0 ate r1\n");*/
            a = (Ialfa[0]+0.0005);
            if ( ((Ialfa[0]+0.0005) < AlfaNewtonIte(a)) && ( AlfaNewtonIte(a) < (Ialfa[1]-0.0005) ) )
            {
                /* printf ("deve ser diferente de %f\n",a);*/
                AlfaNewton();
                printf("Um candidato a alfa eh %f\n", a);
            }
            else
            {
                a = (Ialfa[1]-0.0005);
                /* printf ("deve ser diferente de %f\n",a);*/
                AlfaNewton();
                printf("Um candidato a alfa eh %f\n", a);
            }
        }
    if (( ((avaliaGalfa(Ialfa[1]+0.0005))*(avaliaGalfa(Ialfa[2]-0.005))) < 0))
        {
            /* printf ("f(r1) = %f e f(r2) = %f\n",(avaliaGalfa(Ialfa[1]+0.0005)),(avaliaGalfa(Ialfa[2]-0.005)));*/
            /* printf ("f(r1)*f(r2) = %f\n", ((avaliaGalfa(Ialfa[1]+0.0005))*(avaliaGalfa(Ialfa[2]-0.005))) );*/
            /* printf("segundo if de r1 ate r2\n");*/
            a = (Ialfa[1]+0.0005);
            if ( ((Ialfa[1]+0.0005) < AlfaNewtonIte(a)) && ( AlfaNewtonIte(a) < (Ialfa[2]-0.0005) ) )
            {
                /* printf ("deve ser diferente de %f\n",a);*/
                AlfaNewton();
                printf("Um candidato a alfa eh %f\n", a);
            }
            else
            {
                a = (Ialfa[2]-0.0005);
                /* printf ("deve ser diferente de %f\n",a);*/
                AlfaNewton();
                printf("Um candidato a alfa eh %f\n", a);
            }
        }
        if (( ((avaliaGalfa(Ialfa[2]+0.0005))*(avaliaGalfa(Ialfa[3]-0.005))) < 0))
        {
            /* printf ("f(r2) = %f e f(r3) = %f\n",(avaliaGalfa(Ialfa[2]+0.0005)),(avaliaGalfa(Ialfa[3]-0.005)));*/
            /* printf ("f(r2)*f(r3) = %f\n", ((avaliaGalfa(Ialfa[2]+0.0005))*(avaliaGalfa(Ialfa[3]-0.005))) );*/
            /* printf("terceiro if de r2 ate r3\n");*/
            a = (Ialfa[2]+0.0005);
            if ( ((Ialfa[2]+0.0005) < AlfaNewtonIte(a)) && ( AlfaNewtonIte(a) < (Ialfa[3]-0.0005) ) )
            {
                /* printf ("deve ser diferente de %f\n",a);*/
                AlfaNewton();
                printf("Um candidato a alfa eh %f\n", a);
            }
            else
            {
                a = (Ialfa[3]-0.0005);
                /* printf ("deve ser diferente de %f\n",a);*/
                AlfaNewton();
                printf("Um candidato a alfa eh %f\n", a);
            }
        }
        if (( ((avaliaGalfa(Ialfa[3]+0.0005))*(avaliaGalfa(Ialfa[4]-0.005))) < 0))
        {
            /* printf ("f(r3) = %f e f(r3) = %f\n",(avaliaGalfa(Ialfa[3]+0.0005)),(avaliaGalfa(Ialfa[4]-0.005)));*/
            /* printf ("f(r3)*f(r4) = %f\n", ((avaliaGalfa(Ialfa[3]+0.0005))*(avaliaGalfa(Ialfa[4]-0.005))) );*/
            /* printf("quarto if de r3 ate r4\n");*/
            a = (Ialfa[3]+0.0005);
            if ( ((Ialfa[3]+0.0005) < AlfaNewtonIte(a)) && ( AlfaNewtonIte(a) < (Ialfa[4]-0.0005) ) )
            {
               /* printf ("deve ser diferente de %f\n",a);*/
                AlfaNewton();
               printf("Um candidato a alfa eh %f\n", a);
            }
            else
            {
                a = (Ialfa[4]-0.0005);
                /* printf ("deve ser diferente de %f\n",a);*/
                AlfaNewton();
                printf("Um candidato a alfa eh %f\n", a);
            }
        }
        return;
}

void resolveAlfa(double alpha)
{
        /*printf("alcance maior que zero \n");*/
        vec_aux[0] = A;
        vec_aux[1] = B;
        vec_aux[2] = C;
        vec_aux[3] = D;
        /*printf("Os coeficientes da F(x) sao %f %f %f %f\n",vec_aux[0],vec_aux[1],vec_aux[2],vec_aux[3]);*/
        coef_d1(vec_aux);
        RaizFlin();
        intervalos[1] = raizesPol[0];
        intervalos[2] = raizesPol[1];
        /*printf("raizes de F' sao: %f %f \n",intervalos[1],intervalos[2]);*/
        coef_d2(vec_aux);
        RaizFlinlin();
        intervalos[3] = raizlinlin;
        /*printf("A raiz de F'' eh %f\n",intervalos[3]);*/
        InsertionSort(intervalos,5);
        /*printf("Os intervalos sao: [%f, %f, %f, %f, %f ] \n",intervalos[0],intervalos[1],intervalos[2],intervalos[3],intervalos[4]);*/
        Fx[0] = A;
        Fx[1] = B;
        Fx[2] = C;
        Fx[3] = D;
        if ((horner(Fx,4,intervalos[0]))*(horner(Fx,4,intervalos[1]-0.0005)) < 0 )
        {
            termina = 1;
            x = intervalos[0];
            /*printf("primeiro if 0 a r1\n");*/
            if ( (0 < NewtonIte(Fx)) && (NewtonIte(Fx) < (intervalos[1]-0.0005)) )
            {
                x = 0;
                /*printf ("deve ser diferente de %f\n",x);*/
                Newton(Fx);
                printf("O projetil atinge o solo da montanha em x = %f\n", x);
            }
            else
            {
                x = intervalos[1]-0.0005;
                /*printf ("deve ser diferente de %f\n",x);*/
                Newton(Fx);
                printf("O projetil atinge o solo da montanha em x = %f\n", x);
            }
        }
        if ((((horner(Fx,4,intervalos[1]+0.0005))*(horner(Fx,4,intervalos[2]-0.0005))) < 0 ) && !termina)
        {
            /*printf("segundo if de r1 ate r2\n");*/
            termina = 1;
            if ( ((intervalos[1]+0.0005) < NewtonIte(Fx)) && (NewtonIte(Fx) < (intervalos[2]-0.0005)) )
            {
                /*printf("O vetor Fx eh: %f %f %f %f\n",Fx[0], Fx[1], Fx[2], Fx[3]);*/
                x = (intervalos[1]+0.0005);
                /*printf ("deve ser diferente de %f\n",x);*/
                Newton(Fx);
                printf("O projetil atinge o solo da montanha em x = %f\n", x);
            }
            else
            {
                x = (intervalos[2]-0.0005);
                /*printf ("deve ser diferente de %f\n",x);*/
                Newton(Fx);
                printf("O projetil atinge o solo da montanha em x = %f\n", x);
            }
        }
        if ((((horner(Fx,4,intervalos[2]+0.0005))*(horner(Fx,4,intervalos[3]-0.0005))) < 0 ) && !termina)
        {
            /*printf("terceiro if de r2 ate r3");*/
            termina = 1;
            if ( ((intervalos[2]+0.0005) < NewtonIte(Fx)) && (NewtonIte(Fx) < (intervalos[3]-0.0005)) )
            {
                x = (intervalos[2]+0.0005);
                /*printf ("deve ser diferente de %f\n",x);*/
                Newton(Fx);
                printf("O projetil atinge o solo da montanha em x = %f\n", x);
            }
            else
            {
                x = (intervalos[3]-0.0005);
                /*printf ("deve ser diferente de %f\n",x);*/
                Newton(Fx);
                printf("O projetil atinge o solo da montanha em x = %f\n", x);
            }
        }
        if ((((horner(Fx,4,intervalos[3]+0.0005))*(horner(Fx,4,intervalos[4]-0.0005))) < 0 ) && !termina)
        {
            /*printf("quarto if de r3 ate r4\n");*/
            termina = 1;
            if ( ((intervalos[3]+0.0005) < NewtonIte(Fx)) && (NewtonIte(Fx) < (intervalos[4]-0.0005)) )
            {
                x = (intervalos[3]+0.0005);
                /*printf ("deve ser diferente de %f\n",x);*/
                Newton(Fx);
                printf("O projetil atinge o solo da montanha em x = %f\n", x);
            }
            else
            {
                x = (intervalos[4]-0.0005);
                /*printf ("deve ser diferente de %f\n",x);*/
                Newton(Fx);
                printf("O projetil atinge o solo da montanha em x = %f\n", x);
            }
        }
}


int main()
{
    resolveAlvo();
    float angulo;
    printf("Digite um angulo(em radianos) de lancamento maior que 0 e menor que pi/2:\n");
    scanf("%f",&angulo);
    /*printf("seu angulo eh: %f\n",angulo);*/
    alpha = angulo;
    double alcance = alcMax(alpha);
    /*printf("O alcance eh: %f\n",alcance);*/

    if (alcance <= 0)
    {
        printf("O projetil atinge o solo em (x,y) = (%f,0)\n",alcance);
    }
    avaliaCoef(alpha);
    if (alcance > 0)
    {
        resolveAlfa(alpha);
    }
    return 0;
}
