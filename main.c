#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <unistd.h>

#define TIMEUNIT 0.1

typedef struct body body;

struct body
{
	double x, y, mass;
	double forceX, forceY;
	double velocityX, velocityY;
};

struct quadTree_t
{
        struct quadTree_t * quadrents[4];
	struct quadTree_t * parent;
	double x, y, mass, centerOfMassX, centerOfMassY;
	double xmin, xmax, ymin, ymax;
	int bodyLabel;
	short externFlag;
};

typedef struct quadTree_t quadTree;

//externFlag = 0 means the node is internal
void quadTreeInit(quadTree * par, quadTree * work)
{
        work->parent = par;
	work->quadrents[0] = NULL;
	work->quadrents[1] = NULL;
	work->quadrents[2] = NULL;
	work->quadrents[3] = NULL;
	work->x = 0;
	work->y = 0;
	work->ymin = 0;
	work->ymax = 0;
	work->xmin = 0;
	work->xmax = 0;
	work->mass = 0;
	work->bodyLabel = -1;
	work->externFlag = 0;
	work->centerOfMassX = 0;
	work->centerOfMassY = 0;
}

char *colormap[] = {
		    "31", "32", "33", "34", "35", "36", "37", "30;1", "31;1", "32;1", "33;1", "34;1", "35;1", "36;1", "37;1"
};

int n_colors = 15;

char *bodycolor(int index) {
  return colormap[index % n_colors];
}

void printBody(struct body *b, int index) {
  printf("Body \e[%sm%3d\e[0m - Location: <%+e, %+e> Velocity: <%+e, %+e>\e[K\n",
	 bodycolor(index), index, b->x, b->y, b->velocityX, b->velocityY);
}

//given two masses and two coordinates calculates their center of mass
double centMass(double mass1, double mass2, double coord1, double coord2)
{
	double mass = mass1 + mass2;
	double centerOfMass = (coord1*mass1) + (coord2*mass2);
	centerOfMass /= mass;
	return centerOfMass;
}

double midPoint(double coord1, double coord2)
{
	double retval = coord1 + coord2;
	retval /= 2;
	return retval;
}

void treeInsert(quadTree * root, body * bodies, int index)
{
	//if internal
	if(root->externFlag == 0)
	{
		body * temp = &(bodies[index]);
		root->centerOfMassX = centMass(root->mass, temp->mass, root->centerOfMassX, temp->x);
		root->centerOfMassY = centMass(root->mass, temp->mass, root->centerOfMassY, temp->y);
		root->mass += temp->mass;
		int quadrent;
		if(temp->x < root->x)
		{
			if(temp->y < root->y)
			{
				quadrent = 2;
			}
			else
			{
				quadrent = 1;
			}
		}
		else
		{
			if(temp->y < root->y)
			{
				quadrent = 3;
			}
			else
			{
				quadrent = 0;
			}
		}
		//does quadrent exist?
		if(root->quadrents[quadrent] == NULL)
		{
			quadTree * work = malloc(sizeof(quadTree));
			root->quadrents[quadrent] = work;
			quadTreeInit(root, work);
			work->externFlag = 1;
			work->bodyLabel = index;
			work->centerOfMassX = temp->x;
			work->centerOfMassY = temp->y;
			work->mass = temp->mass;
			switch(quadrent)
			{
				case 0: 
					work->xmin = root->x;
					work->xmax = root->xmax;
					work->ymin = root->y;
					work->ymax = root->ymax;
					break;
				case 1:
					work->xmin = root->xmin;
					work->xmax = root->x;
					work->ymin = root->y;
					work->ymax = root->ymax;
					break;
				case 2:
					work->xmin = root->xmin;
					work->xmax = root->x;
					work->ymin = root->ymin;
					work->ymax = root->y;
					break;
				case 3:
					work->xmin = root->x;
					work->xmax = root->xmax;
					work->ymin = root->ymin;
					work->ymax = root->y;
					break;
			}
			work->x = midPoint(work->xmax, work->xmin);
			work->y = midPoint(work->ymax, work->ymin);
		}
		else
		{
			treeInsert(root->quadrents[quadrent], bodies, index);
		}
	}
	//else external
	else
	{
		int temp = root->bodyLabel;
		root->externFlag = 0;
		root->centerOfMassX = 0;
		root->centerOfMassY = 0;
		root->mass = 0;
		root->bodyLabel = -1;
		treeInsert(root, bodies, temp);
		treeInsert(root, bodies, index);
	}
}

void barnsHutConstruct(body * bodies, int size, quadTree * root)
{
	int i;
	for(i = 0; i < size; i++)
	{
		treeInsert(root, bodies, i);
	}
}

void generateBodies(body * bodies, int size, double xdim, double ydim, long maxvel)
{
	int i;
	for(i = 0; i < size; i++)
	{
		bodies[i].x = rand() % (int) xdim;
		bodies[i].y = rand() % (int) ydim;
		bodies[i].mass = rand();
		bodies[i].velocityX = rand() % (int) maxvel;
		bodies[i].velocityY = rand() % (int) maxvel;
		bodies[i].forceX = 0;
		bodies[i].forceY = 0;
	}
}

double dist1D(double coord1, double coord2)
{
	double retval = coord2-coord1;
	retval *= retval;
	retval = sqrt(retval);
	return retval;
}

double distance(double coord1x, double coord2x, double coord1y, double coord2y)
{
	double retval = coord2x-coord1x;
	double temp = coord2y - coord1y;
	temp *= temp;
	retval *= retval;
	retval += temp;
	retval = sqrt(retval);
	return retval;
}

void calculateForces(body * bodies, int size, quadTree * root, int index, double threshold)
{
	double G = 0.0000000000667;
	if(root->externFlag != 0 && root->bodyLabel != index)
	{
		double dist = distance(root->centerOfMassX, bodies[index].x, root->centerOfMassY, bodies[index].y);
		double force = (G * root->mass * bodies[index].mass) / (dist*dist);
		bodies[index].forceX += force * dist1D(root->centerOfMassX, bodies[index].x) / dist;
		bodies[index].forceY += force * dist1D(root->centerOfMassY, bodies[index].y) / dist;
	}
	else if(((root->xmax - root->xmin) / distance(root->centerOfMassX, bodies[index].x, root->centerOfMassY, bodies[index].y)) < threshold)
	{
		double dist = distance(root->centerOfMassX, bodies[index].x, root->centerOfMassY, bodies[index].y);
		double force = (G * root->mass * bodies[index].mass) / (dist*dist);
		bodies[index].forceX += force * dist1D(root->centerOfMassX, bodies[index].x) / dist;
		bodies[index].forceY += force * dist1D(root->centerOfMassY, bodies[index].y) / dist;

	}
	else
	{
		if(root->quadrents[0] != NULL)
		{
			calculateForces(bodies, size, root->quadrents[0], index, threshold);
		}
		if(root->quadrents[1] != NULL)
		{
			calculateForces(bodies, size, root->quadrents[1], index, threshold);
		}
		if(root->quadrents[2] != NULL)
		{
			calculateForces(bodies, size, root->quadrents[2], index, threshold);
		}
		if(root->quadrents[3] != NULL)
		{
			calculateForces(bodies, size, root->quadrents[3], index, threshold);
		}
	}
}

void resetForces(body * bodies, int size)
{
	int i;
	for(i = 0; i < size; i++)
	{
		bodies[i].forceX = 0;
		bodies[i].forceY = 0;
	}
}

void updateVelocityPosition(body * bodies, int index)
{
	bodies[index].velocityX += TIMEUNIT * bodies[index].forceX / bodies[index].mass;
	bodies[index].velocityY += TIMEUNIT * bodies[index].forceY / bodies[index].mass;
    bodies[index].x += TIMEUNIT * bodies[index].velocityX;
	bodies[index].y += TIMEUNIT * bodies[index].velocityY;
}

void wrapPositions(body *bodies, quadTree *root, int size) {
  int i;
  for(i = 0; i < size; i++) {
    if(bodies[i].x > root->xmax) {
      bodies[i].x = root->xmax;
      bodies[i].velocityX = -bodies[i].velocityX;
    }
    if(bodies[i].x < root->xmin) {
      bodies[i].x = root->xmin;
      bodies[i].velocityX = -bodies[i].velocityX;
    }
    if(bodies[i].y > root->ymax) {
      bodies[i].y = root->ymax;
      bodies[i].velocityY = -bodies[i].velocityY;
    }
    if(bodies[i].y < root->ymin) {
      bodies[i].y = root->ymin;
      bodies[i].velocityY = -bodies[i].velocityY;
    }
  }
}

void timeStep(body * bodies, int size, quadTree * root, double threshold, int startIndex, int endIndex)
{
	int i;
	for(i = startIndex; i < endIndex; i++)
	{
		calculateForces(bodies, size, root, i, threshold);
		updateVelocityPosition(bodies,i);
	}
	resetForces(bodies, size);
	wrapPositions(bodies, root, size);
}

void pack(double * move, body * bodies, int size)
{
	int i;
	for(i = 0; i < size; i++)
	{
		move[(7*i)] = bodies[i].x;
		move[(7*i)+1] = bodies[i].y;
		move[(7*i)+2] = bodies[i].mass;
		move[(7*i)+3] = bodies[i].forceX;
		move[(7*i)+4] = bodies[i].forceY;
		move[(7*i)+5] = bodies[i].velocityX;
		move[(7*i)+6] = bodies[i].velocityY;
	}
}

void unpack(double * move, body * bodies, int size)
{
	int i;
	for(i = 0; i < size; i++)
	{
		bodies[i].x = move[(7*i)];
		bodies[i].y = move[(7*i)+1];
		bodies[i].mass = move[(7*i)+2];
		bodies[i].forceX = move[(7*i)+3];
		bodies[i].forceY = move[(7*i)+4];
		bodies[i].velocityX = move[(7*i)+5];
		bodies[i].velocityY = move[(7*i)+6];
	}
}

void unpackRank(double * move, body * bodies, int size, int rangeMin, int rangeMax)
{
	int i;
	for(i = rangeMin; i < rangeMax; i++)
	{
		bodies[i].x = move[(7*i)];
		bodies[i].y = move[(7*i)+1];
		bodies[i].mass = move[(7*i)+2];
		bodies[i].forceX = move[(7*i)+3];
		bodies[i].forceY = move[(7*i)+4];
		bodies[i].velocityX = move[(7*i)+5];
		bodies[i].velocityY = move[(7*i)+6];
	}
}

void freeTree(quadTree * root)
{
	if(root->quadrents[0] != NULL)
	{
		freeTree(root->quadrents[0]);
	}
	if(root->quadrents[1] != NULL)
	{
		freeTree(root->quadrents[1]);	
	}
	if(root->quadrents[2] != NULL)
	{
		freeTree(root->quadrents[2]);
	}
	if(root->quadrents[3] != NULL)
	{
		freeTree(root->quadrents[3]);
	}
	free(root);
}

/* argv contians the following in the following order
executable name
accuracy threshold
simulation space y diminsion max
simulation space x dimintion max
Number of bodys to spawn for simulation
Number of time steps to simulate
Max Initial body velocity
*/
int main(int argc, char* argv[])
{
	if(argc != 7)
	{
		printf("wrong number of arguments. usage: ./program [accuracy threshold] [y diminsion] [x diminsion] [number of bodies to simulate] [number of time steps to simulate for] [max initial body velocity]\n");
		printf("warning user input is not validated, invalid arguments will result in a crash.\n");
		printf("suggested threshold value 0.5, value must be within the range [0,1]\n");
		printf("note max velocity and timesteps both expect longs, and size expects an int, all other parameters are doubles.\n");
	}
	double threshold = strtod(argv[1], NULL);
	double ydim = strtod(argv[2], NULL);
	double xdim = strtod(argv[3], NULL);
	int size = atoi(argv[4]);
	long timesteps = atol(argv[5]);
	long maxvel = atol(argv[6]);
	
	MPI_Init(NULL, NULL);

	//init cores
	int cores;
	MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Comm_size(comm, &cores);
	int rank;

	//sets rank = to what core
	MPI_Comm_rank(comm, &rank);
	MPI_Status status;
	//determain waht part of the structure each node will work on
	int range = size/cores;
	int rangeMin = rank * range;
	int rangeMax = rangeMin + range;
	//determine if there are leftover values that the root needs to do
	int remainderFlag = 0;
	if(size % cores != 0)
	{
		remainderFlag = 1;
	}
	
	int excessRangemin;
	if(rank == 0 && remainderFlag == 1)
	{
		excessRangemin = ((cores - 1) * range) + range;
	}
	
	srand(time(NULL));
	
	quadTree world;
	body * bodies;
	bodies = malloc(sizeof(body) * size);
	if(rank == 0)
	{
		generateBodies(bodies,size,xdim,ydim,maxvel);
	}
	double * move = malloc(sizeof(double)*7*size);
	int i;
	printf("\e[2J");
	//do the time steps
	for(i = 0; i < timesteps; i++)
	{
		if(rank == 0)
		{
			pack(move,bodies,size);
			int n;
			for(n=0;n<cores;n++)
			{
				MPI_Send(move, 7 * size, MPI_DOUBLE, n, 0, comm);
			}
		}
		else
		{
			MPI_Recv(move, 7 * size, MPI_DOUBLE, 0, 0, comm, &status);
			unpack(move, bodies, size);
		}
		quadTreeInit(NULL, &world);
		world.xmin = 0;
		world.ymin = 0;
		world.xmax = xdim;
		world.ymax = ydim;
		world.x = midPoint(world.xmax, world.xmin);
		world.y = midPoint(world.ymax, world.ymin);
		barnsHutConstruct(bodies,size,&world);
		
		//insert a print or write out function here to display the results
		timeStep(bodies, size, &world, threshold,rangeMin, rangeMax);	
		if(rank == 0 && remainderFlag == 1)
		{
			timeStep(bodies,size,&world, threshold, excessRangemin, size);
		}
		
		if(rank == 0)
		{
			int n;
			for(n = 0; n < cores; n++)
			{
				MPI_Recv(move, 7 * size, MPI_DOUBLE, n, 0, comm, &status);
				int tempMin = n * range;
				int tempMax = tempMin + range;
				unpackRank(move, bodies, size, tempMin, tempMax);
			}
			printf("\e[1;1HTime: %d\n+", i);
			for(n = 0; n < world.xmax + 1; n++) printf("-");
			printf("+\n");
			for(n = 0; n < world.ymax + 1; n++) {
			  int m;
			  printf("|");
			  for(m = 0; m < world.xmax + 1; m++) printf(" ");
			  printf("|\n");
			}
			printf("+");
			for(n = 0; n < world.xmax + 1; n++) printf("-");
			printf("+\n");
			for(n = 0; n < size; n++) {
			  printf("\e[%dA\e[%dG\e[%sm%d\e[0m\e[%dB", (int) bodies[n].y + 2, (int) bodies[n].x + 1, bodycolor(n), n, (int) bodies[n].y + 2);
			}
			fflush(stdout);
			printf("\n");
			for(n = 0; n < size; n++) {
			  printBody(&bodies[n], n);
			}
			sleep(1);
		}
		else
		{
			pack(move, bodies, size);
			MPI_Send(move, 7 * size, MPI_DOUBLE, 0, 0, comm);
		}
	}
	
	free(bodies);
	free(move);
	
	MPI_Finalize();
	return 0;
}
