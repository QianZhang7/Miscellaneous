// qsort.c - Parallel sorting algorithm based on quicksort
// compile: mpicc -Wall -O -o qsort qsort.c
// run:     mpirun -np num_procs qsort in_file out_file

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>


// swap entries in array v at positions i and j; used by quicksort
void swap(int * v, int i, int j)
{
  int t = v[i];
  v[i] = v[j];
  v[j] = t;
}


// (quick) sort slice of array v; slice starts at s and is of length n
void quicksort(int * v, int s, int n)
{
  int x, p, i;
  // base case?
  if (n <= 1)
    return;
  // pick pivot and swap with first element
  x = v[s + n/2];
  swap(v, s, s + n/2);
  // partition slice starting at s+1
  p = s;
  for (i = s+1; i < s+n; i++)
    if (v[i] < x) {
      p++;
      swap(v, i, p);
    }
  // swap pivot into place
  swap(v, s, p);
  // recurse into partition
  quicksort(v, s, p-s);
  quicksort(v, p+1, s+n-p-1);
}


// merge two sorted arrays v1, v2 of lengths n1, n2, respectively
int * merge(int * v1, int n1, int * v2, int n2)
{
  int * result = (int *)malloc((n1 + n2) * sizeof(int));
  int i = 0;
  int j = 0;
  int k;
  for (k = 0; k < n1 + n2; k++) {
    if (i >= n1) {
      result[k] = v2[j];
      j++;
    }
    else if (j >= n2) {
      result[k] = v1[i];
      i++;
    }
    else if (v1[i] < v2[j]) { // indices in bounds as i < n1 && j < n2
      result[k] = v1[i];
      i++;
    }
    else { // v2[j] <= v1[i]
      result[k] = v2[j];
      j++;
    }
  }
  return result;
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
  int dist;
  long n;
  int * data = NULL;
  long c, s;
  int * chunk;
  int o;
  int * other;
  int step;
  int p, id;
  MPI_Status status;
  double elapsed_time;
  FILE * file = NULL;
  int i;
  if (argc>2)
	{
		n=atoi(argv[1]);
		dist=atoi(argv[2]);
	}	
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);
  c=n/p;	//size of array on each cluster
  data=(int*)malloc(c*sizeof(int));
  numgen(data,c,dist,id);

  // start the timer
  MPI_Barrier(MPI_COMM_WORLD);
  elapsed_time = - MPI_Wtime();

  // broadcast size
  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

  // compute chunk size
  c = n/p; if (n%p) c++;

  // scatter data
  chunk = (int *)malloc(c * sizeof(int));
  MPI_Scatter(data, c, MPI_INT, chunk, c, MPI_INT, 0, MPI_COMM_WORLD);
  free(data);
  data = NULL;

  // compute size of own chunk and sort it
  if (n >= c * (id+1))
    s = c;
  else
    s = n - c * id;
  quicksort(chunk, 0, s);

  // up to log_2 p merge steps
  for (step = 1; step < p; step = 2*step) {
    if (id % (2*step)) {
      // id is no multiple of 2*step: send chunk to id-step and exit loop
      MPI_Send(chunk, s, MPI_INT, id-step, 0, MPI_COMM_WORLD);
      break;
    }
    // id is multiple of 2*step: merge in chunk from id+step (if it exists)
    if (id+step < p) {
      // compute size of chunk to be received
      if (n >= c * (id+2*step))
        o = c * step;
      else
        o = n - c * (id+step);
      // receive other chunk
      other = (int *)malloc(o * sizeof(int));
      MPI_Recv(other, o, MPI_INT, id+step, 0, MPI_COMM_WORLD, &status);
      // merge and free memory
      data = merge(chunk, s, other, o);
      free(chunk);
      free(other);
      chunk = data;
      s = s + o;
    }
  }

  // stop the timer
  elapsed_time += MPI_Wtime();

    // printf("Quicksort %d ints on %d procs: %f secs\n", n, p, elapsed_time);
	if(id==0)
    printf("%2d %f\n",  p, elapsed_time);
 

  MPI_Finalize();
  return 0;
}
