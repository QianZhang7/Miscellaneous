#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

void swap(int * v, int i, int j)
{
  int t = v[i];
  v[i] = v[j];
  v[j] = t;
}


/* (bubble) sort array v; array is of length n */
void bubblesort(int * v, unsigned long n)
{
  long i, j;
  for (i = n-2; i >= 0; i--)
    for (j = 0; j <= i; j++)
      if (v[j] > v[j+1])
        swap(v, j, j+1);
}
void numgen(int * data, unsigned long c,int dist,int me)
{
	int range=1024;
	int quarter=range/4;
	int q;
	unsigned long i;
	srand(me+1);
	if (dist==0)
	{
		for (i=0;i<c;i++)
			////data uniform distribution ;
			data[i]=rand()%range;
			//printf("%d\n",data[i]);
				
	}
	else
	{
		for (i=0;i<c;i++)
		{
			////data quarter distribution
			q=rand()%8; ///quarter selection
			if (q<1){
				data[i]=rand()%quarter;
			}
			else
			{	data[i]=rand()%(range-quarter)+quarter; 
				//cnt++;
			}
		}
	}

}
int main(int argc, char ** argv)
{
	unsigned long total=256,c;
	int dist=0,swlen=16;
	int *data=NULL;
	int *tempdata=NULL;
	int *first=NULL,*last=NULL;
	int *top=NULL,*buttom=NULL;
	MPI_Status stat; 
	int i,j,me,np;
	double t1,t2;
	if (argc>2)
	{
		total=atoi(argv[1]);
		dist=atoi(argv[2]);
	}	
	MPI_Init(&argc,&argv);                 // Initialize
	MPI_Comm_size(MPI_COMM_WORLD, &np);  // Get # processors
	MPI_Comm_rank(MPI_COMM_WORLD, &me);      // Get my rank (id)
	c=total/np;	//size of array on each cluster
	if (2*swlen>c) swlen=c/2;

	data=(int*)malloc(c*sizeof(int));
	tempdata=(int*)malloc(c*sizeof(int)); // temporary receive buffer
	first=(int*)malloc(np*sizeof(int)); //first and last element in each node
	last=(int*)malloc(np*sizeof(int));
	   
	// Generate data array on each cluster
	numgen(data,c,dist,me);
	//debug: print data
/*	
	printf("-----%d:Initial-----\n ",me);
	for (i=0;i<c;i++)
 	{
		printf("%d\t",data[i]);
		tempdata[i]=-1;
	}
	printf("\n");
*/	
	t1=MPI_Wtime();
 	//sort locally
	bubblesort (data,c);
	unsigned long subc=c/np; //chunk size of subarray
	//printf("chunk size is %d \n",subc);
	i=0;
	MPI_Barrier(MPI_COMM_WORLD);	
	// scatter array to all nodes
	for (i=0;i<np;i++)
	{
		MPI_Scatter(data,subc,MPI_INT,tempdata+i*subc,subc,MPI_INT,i,MPI_COMM_WORLD);
	}
	data=tempdata;
	//sort locally
	bubblesort(data,c);
	int finish=0;
	tempdata=data-swlen;	
	MPI_Barrier(MPI_COMM_WORLD);

	while(finish<np-1)
	{
		if (me!=np-1)
			MPI_Send(data+c-swlen,swlen,MPI_INT,me+1,0,MPI_COMM_WORLD);
		if (me!=0)
			MPI_Recv(tempdata,swlen,MPI_INT,me-1,0,MPI_COMM_WORLD,&stat);
		if (me!=0)
			MPI_Send(data,swlen,MPI_INT,me-1,0,MPI_COMM_WORLD);
		if (me!=np-1)
			MPI_Recv(data+c,swlen,MPI_INT,me+1,0,MPI_COMM_WORLD,&stat);
		if (me!=0)
			bubblesort(tempdata,2*swlen);
		if (me!=np-1)
			bubblesort(data+c-swlen,2*swlen);	
		bubblesort(data,c);
		MPI_Barrier(MPI_COMM_WORLD);

	//debug
	/*	
		printf("%d: ",me);
		for (i=0;i<c+2*swlen;i++)
		{
			printf("%d\t",tempdata[i]);
		}
		printf("\n");
	*/
		MPI_Gather(data,1,MPI_INT,first,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Gather(data+c-1,1,MPI_INT,last,1,MPI_INT,0,MPI_COMM_WORLD);
		
		if (me==0)
		{
			finish=0;
			for(i=1;i<np;i++)
			{	
				//printf("first: %d last: %d\n",first[i],last[i-1]);
				if (first[i]>=last[i-1])
						finish++;
			}
		}
		
	MPI_Bcast(&finish,1,MPI_INT,0,MPI_COMM_WORLD);
	//printf("%d  finish: %d\n",me,finish);
	MPI_Barrier(MPI_COMM_WORLD);
	}
	t2=MPI_Wtime();
	//debug: print data
	/*printf("=====%d:Sorted===== ",me);
	for (i=0;i<c;i++)
 	{
		printf("%d\t",data[i]);
		tempdata[i]=-1;
	}
	printf("\n");
	*/
	if(me==0)
	printf("\nSorting Success! Time: %f\n",t2-t1);
	free(data);
	free(tempdata);
	free(first);
	free(last);
	MPI_Finalize();
	return 0;
}
