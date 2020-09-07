/*
   Calculating the value of pi using reduction : Serial Implementation
Author : Omkar Damle.
Date : August 2016.
 */

#include<stdio.h>
#include<math.h>
#include<omp.h>
#include <limits.h>
#include <stdbool.h>
#include<time.h>
#include<string.h>
#include<stdlib.h>

//  Using the MONOTONIC clock 
#define CLK CLOCK_MONOTONIC
int N,k,NoCenter;
double **data;
double **centers;
int *output;
/* Function to compute the difference between two points in time */
struct timespec diff(struct timespec start, struct timespec end);

/* 
   Function to computes the difference between two time instances

   Taken from - http://www.guyrutenberg.com/2007/09/22/profiling-code-using-clock_gettime/ 

   Further reading:
http://stackoverflow.com/questions/6749621/how-to-create-a-high-resolution-timer-in-linux-to-measure-program-performance
http://stackoverflow.com/questions/3523442/difference-between-clock-realtime-and-clock-monotonic
 */
struct timespec diff(struct timespec start, struct timespec end){
	struct timespec temp;
	if((end.tv_nsec-start.tv_nsec)<0){
		temp.tv_sec = end.tv_sec-start.tv_sec-1;
		temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
	}
	else{
		temp.tv_sec = end.tv_sec-start.tv_sec;
		temp.tv_nsec = end.tv_nsec-start.tv_nsec;
	}
	return temp;
}

	void getData(int n){
		FILE* f=fopen("input.csv","r");
		if(f==NULL){
			printf("ERROR OPENING FILE\n");
			exit(1);
		}
	
		int i,j;
		char temp[1],tempRead[1];
		while(1){
			fscanf(f,"%c",temp);
			if(*temp=='\n'){
				N++;
				break;
			}
			else if(*temp==',') k++;		
		}
		k++;
	
//		while(1){
//			if(fscanf(f,"%c",temp)==EOF)
//				break;
//			else if(*temp=='\n'){
//				N++;
//			}
//		}
		N=n;
		fclose(f);
		f=fopen("input.csv","r");	
		
		data=malloc(sizeof(double *)*N);
		for(i=0;i<N;i++)
			data[i]=malloc(sizeof(double *)*k);
	
		centers=malloc(sizeof(double *)*NoCenter);
		for(i=0;i<NoCenter;i++)
			centers[i]=malloc(sizeof(double *)*k);
	
		output=malloc(sizeof(double *)*N);	
	
		double val;
		for(i=0;i<N;i++){
			for(j=0;j<k;j++){
				fscanf(f,"%lf",&val);
				fscanf(f,"%c",tempRead);
				data[i][j]=val;
				
			}
		}
		fclose(f);
	}
	
	double FindDistance(double a[],double b[]){
		double ans=0;
		int i;
		for(i=0;i<k;i++){
			ans+=((a[i]-b[i])*(a[i]-b[i]));
		}
		return sqrt(ans);
	}
	
	void ChooseCenter(){
		int i,j,num;
		bool temp[N];
		for(i=0;i<N;i++) temp[i]=false;
		i=0;
		while(i != NoCenter){
			num=rand()%N;	// Generate random numbers from 0 to n-1
			if(temp[num] == false){
				temp[num]=true;
				for(j=0;j<k;j++) centers[i][j]=data[num][j];			
				i++;
			} 	 			
		}
	}
	
	double Dunn_index()
	{
		double  Dis;		
		double  MinClusterDis=INT_MAX;
		int i,j;
		for(i=0;i<NoCenter;i++){
			for(j=i+1;j<NoCenter;j++)
			{
				Dis= FindDistance(centers[i],centers[j]);	
				if(Dis < MinClusterDis)
					MinClusterDis=Dis;		
			}
		}
		double  MaxSize=INT_MIN;		
		for(i=0;i<NoCenter;i++){
			int c=0;
			double  Size=0;
			for(j=0;j<N;j++){
				if(output[j]==i+1){
					Size+=(FindDistance(centers[i],data[j]));
					c++;
				}
			}
			Size/=c;
			if(Size > MaxSize)
				MaxSize=Size;
		}
		return (MinClusterDis/MaxSize);
	}

int main(int argc, char* argv[])
{
	struct timespec start_e2e, end_e2e, start_alg, end_alg, e2e, alg;
	/* Should start before anything else */
	clock_gettime(CLK, &start_e2e);

	/* Check if enough command-line arguments are taken in. */
	if(argc < 3){
		printf( "Usage: %s n p \n", argv[0] );
		return -1;
	}

	int n=atoi(argv[1]);	/* size of input array */
	int p=atoi(argv[2]);	/* number of processors*/
	char *problem_name = "matrix_multiplication";
	char *approach_name = "block";
//	char buffer[10];
//	FILE* inputFile;
	FILE* outputFile;
	//	inputFile = fopen(argv[3],"r");

	char outputFileName[50];		
	sprintf(outputFileName,"output/%s_%s_%s_%s_output.txt",problem_name,approach_name,argv[1],argv[2]);

	
	
	NoCenter=30;
	getData(n);
	ChooseCenter();
	int i,j,It=10;
	clock_gettime(CLK, &start_alg);	/* Start the algo timer */

	/*----------------------Core algorithm starts here----------------------------------------------*/
	while(It--){
		for(i=0;i<N;i++){
			double Min_Dis=INT_MAX;
			for(j=0;j<NoCenter;j++){
				double Dis=FindDistance(centers[j],data[i]);
				if(Min_Dis>Dis){
					Min_Dis=Dis;
					output[i]=j+1;
				}
			}
			for(j=0;j<k;j++){
				centers[output[i]-1][j]=(centers[output[i]-1][j]+data[i][j])/2.0;
			}
			//printf("YES");
		}
//		double dun=Dunn_index();
//		printf("Dun_Index:-%f\n",dun);
	}
	/*----------------------Core algorithm finished--------------------------------------------------*/

	clock_gettime(CLK, &end_alg);	/* End the algo timer */
	/* Ensure that only the algorithm is present between these two
	   timers. Further, the whole algorithm should be present. */

	
	
	/* Should end before anything else (printing comes later) */
	clock_gettime(CLK, &end_e2e);
	e2e = diff(start_e2e, end_e2e);
	alg = diff(start_alg, end_alg);

	outputFile = fopen(outputFileName,"w");

	//fprintf(outputFile,"%.8f\n",n);		

	/* problem_name,approach_name,n,p,e2e_sec,e2e_nsec,alg_sec,alg_nsec
	   Change problem_name to whatever problem you've been assigned
	   Change approach_name to whatever approach has been assigned
	   p should be 0 for serial codes!! 
	 */
	printf("%s,%s,%d,%d,%d,%ld,%d,%ld\n", problem_name, approach_name, n, p, e2e.tv_sec, e2e.tv_nsec, alg.tv_sec, alg.tv_nsec);

	return 0;

}

