#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
const double MIN = 0;
const double MAX = 5000;  //Максимальные U и I
const int ERROR = -1;
const int nmeas = 12;
const int len1 = 20;
const int len2 = 30;
const int len3 = 50;
const double S = 0.00107;
const double amp_err = 1.875;
const double maxI = 400;
const double tol = 3e-8;
double R_len1 = 0, R_len2 = 0, R_len3 = 0;  //ВНИМАНИЕ: глобальные переменные
                                            //Функции
int Read_Data(double *I, double *U_len1, double *U_len2, double *U_len3, char read_file[]);
double Least_Square(double U[], double I[]);
void Calculate(double I[], double U_len1[], double U_len2[], double U_len3[], \
               double *p_len1, double *p_len2, double *p_len3, \
			   double *sigma_len1, double *sigma_len2, double *sigma_len3);
double Rand_Error(double U[], double I[], double R);
int Write_Data (double *p_len1, double *p_len2, double *p_len3, \
                double *sigma_len1, double *sigma_len2, double *sigma_len3, char write_file[]); 
void MyFree (double **arr);
int Check_Data(double *I, double *U_len1, double *U_len2, double *U_len3, int line);

//-------------------------------------------------------------------------------------------------------------------------------------	

int main()
{
	double p_len1 = 0, p_len2 = 0, p_len3 = 0, sigma_len1 = 0, sigma_len2 = 0, sigma_len3 = 0;
	char read_file[] = {"data.txt"};
	char write_file[] = {"result.txt"};
	double *I = calloc(nmeas, sizeof(double));
	double *U_len1 = calloc(nmeas, sizeof(double));
	double *U_len2 = calloc(nmeas, sizeof(double));
	double *U_len3 = calloc(nmeas, sizeof(double));
	if (I == NULL || U_len1 == NULL || U_len2 == NULL || U_len3 == NULL)
		{printf("Can't create array.\n"); return ERROR;}
	if (Read_Data(I, U_len1, U_len2, U_len3, read_file) == ERROR)
		{printf("Function Read_Data failed.\n"); return ERROR;}
	Calculate(I, U_len1, U_len2, U_len3, &p_len1, &p_len2, &p_len3, &sigma_len1, &sigma_len2, &sigma_len3);
	if (Write_Data(&p_len1, &p_len2, &p_len3, &sigma_len1, &sigma_len2, &sigma_len3, write_file) == ERROR)
		{printf("Function Write_Data failed.\n"); return ERROR;}
	MyFree(&I);
	MyFree(&U_len1);
	MyFree(&U_len2);
	MyFree(&U_len3);
	return 0;
}

//-------------------------------------------------------------------------------------------------------------------------------------	

int Read_Data(double *I, double *U_len1, double *U_len2, double *U_len3, char read_file[])
{
	FILE* read = fopen(read_file, "r");
	if (!read) {printf("Cannot open data.txt.\n"); return ERROR;}
	fscanf (read, "%*s %*s %*s %*s %*s %*s %*s");  //Пропуск заголовков
	int line = 0;
	for(;;) {
		if (EOF == fscanf (read, "%lg %lg %lg %lg", &I[line], &U_len1[line], &U_len2[line], &U_len3[line])) break;
		if (Check_Data(I, U_len1, U_len2, U_len3, line) == ERROR)
			{printf("Wrong data in line %d.\n", line+1); return ERROR;}
		line++;
	}
	fclose (read);
	return line;
}

//-------------------------------------------------------------------------------------------------------------------------------------

double Least_Square(double U[], double I[])      //Расчет среднего сопротивления
{
	double num = 0, denom = 0, R = 0;
	for (int i=0; i<nmeas; ++i) {
		num += U[i]*I[i];
		denom += I[i]*I[i];
	}
	return R = num/denom;
}

//-------------------------------------------------------------------------------------------------------------------------------------

double Rand_Error(double U[], double I[], double R)     //Расчет случайной погрешности
{
	double num = 0, denom = 0;
	for (int i=0; i<nmeas; ++i) {
		num += U[i]*U[i];
		denom += I[i]*I[i];
	}
	return 1/sqrt(nmeas)*sqrt(num/denom - R*R);
}
	
//-------------------------------------------------------------------------------------------------------------------------------------

void Calculate(double I[], double U_len1[], double U_len2[], double U_len3[], \
               double *p_len1, double *p_len2, double *p_len3, \
			   double *sigma_len1, double *sigma_len2, double *sigma_len3)
{
	double rand_len1 = 0, rand_len2 = 0, rand_len3 = 0, syst_len1 = 0, syst_len2 = 0, syst_len3 = 0;
	R_len1 = Least_Square(U_len1, I);
	R_len2 = Least_Square(U_len2, I);
	R_len3 = Least_Square(U_len3, I);
	
	*p_len1 = (R_len1*S)/len1;
	*p_len2 = (R_len2*S)/len2;            //Расчет удельного сопротивления
	*p_len3 = (R_len3*S)/len3;
	
	rand_len1 = Rand_Error (U_len1, I, R_len1);
	rand_len2 = Rand_Error (U_len2, I, R_len2);
	rand_len3 = Rand_Error (U_len3, I, R_len3);
	
	syst_len1 = R_len1*amp_err/maxI;
	syst_len2 = R_len2*amp_err/maxI;      //Расчет систематической погрешности
	syst_len3 = R_len3*amp_err/maxI;
	
	*sigma_len1 = sqrt(rand_len1*rand_len1 + syst_len1*syst_len1);
	*sigma_len2 = sqrt(rand_len2*rand_len2 + syst_len2*syst_len2);   //Расчет итоговой погрешности
	*sigma_len3 = sqrt(rand_len3*rand_len3 + syst_len3*syst_len3);
}

//-------------------------------------------------------------------------------------------------------------------------------------

int Write_Data (double *p_len1, double *p_len2, double *p_len3, \
                double *sigma_len1, double *sigma_len2, double *sigma_len3, char write_file[]) 
{
	FILE* result = fopen(write_file, "w");
	if (!result) {printf("Cannot open result.txt.\n"); return ERROR;}
	fprintf (result,"\t  Среднее сопротивление\tУдельное сопротивление\n");
	fprintf (result,"l=20см\t\t%.4g±%.3g\t\t%.4g\n", R_len1, *sigma_len1, *p_len1*10000);
	fprintf (result,"l=30см\t\t%.4g±%.3g\t\t%.4g\n", R_len2, *sigma_len2, *p_len2*10000);
	fprintf (result,"l=50см\t\t%.4g±%.3g\t\t%.4g\n", R_len3, *sigma_len3, *p_len3*10000);
	printf ("Mission completed!\n");
	fclose(result);
	return 0;
}

//-------------------------------------------------------------------------------------------------------------------------------------

void MyFree (double **arr)
{
	free(*arr);
	*arr = NULL;
}

//-------------------------------------------------------------------------------------------------------------------------------------

int Check_Data(double *I, double *U_len1, double *U_len2, double *U_len3, int line)
{
	if (I[line] <= (MIN-tol) || I[line] > (MAX+tol) || U_len1[line] <= (MIN-tol) || U_len1[line] > (MAX+tol) \
		|| U_len2[line] <= (MIN-tol) || U_len2[line] > (MAX+tol) || U_len3[line] <= (MIN-tol) || U_len3[line] > (MAX+tol))
		return ERROR;
		else return 0;
}

//Пасхалка
