#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <inttypes.h>

/////////////////
/*Global Values*/
/////////////////

#define uint unsigned int
#define PI 3.14159265359
#define GOLDEN_RATIO 0.61803398875
#define GOLDEN_ANGLE_DEGREES 137.5077640500378546463487
#define GOLDEN_ANGLE_RADS 2.39996322972865332


//////////////////////////
/*Global Data Structures*/
//////////////////////////

typedef struct Particle {

	float latitude;
	float longitude;

	float height;	//delta from sea level

	bool removed;

}Particle;




////////////////////////////////////
/*Functional Function Doing Things*/
////////////////////////////////////

/**
* Checks to see if an uint is a power of two
*
* @param particle - The integer to check
* @return bool - Is it a power of two
*/
bool IsPower2(uint x) {

	return x && !(x & (x - 1));

}


/**
* returns the amount of particles to simulate for the current rank. Mallocs a new array
* and copies over the previous one if it exists
*
* @param rankID - The cuurent rank
* @param rankCount - The total rank count
* @param particleCount - The particle count to distribute
* @return particles - The particles to simulate
* @return count - The number of particles particles to simulate
* @return index - The index of the first particle for this rank
*/
inline void ParticlestoSimulate(uint rankID, uint rankCount, uint particleCount, Particle*& particles, uint& count, uint& index) {

	count = particleCount / rankCount;

	uint remainder = particleCount % rankCount;
	if (rankID < remainder) {

		++count;
		index = count*rankID;

	}
	else {

		index = (count + 1)*rankID;

	}

	Particle* temp = (Particle*)malloc(sizeof(Particle)*count);

	if (particles) {
		//TODO: Copy over particles that are not removed 
		//memcpy(temp, particles, initialCount*sizeof(struct Particle));         
		free(particles);
	}

	particles = temp;

}

/**
* Setup the Particle.
*
* @param particle - The particle obvs
*/

void InitParticle(Particle& particle, uint particleID, uint particleCount) {

	particle.height = 0.0f;
	particle.removed = false;

	particle.longitude = GOLDEN_ANGLE_RADS*particleID;

	//TODO: is this needed?

	/*lon /= 2 * PI;
	lon -= floor(lon);
	lon *= 2 * PI;*/

	if (particle.longitude > PI) {
		particle.longitude -= 2 * PI;
	}

	particle.latitude = asin(-1 + 2 * particleID / (float)particleCount);

}

//////////////////////
/*Mainly the Program*/
//////////////////////

int main(int argc, char **argv)
{

	////////////////////////////
	/*Setup Global Information*/
	////////////////////////////


	if (argc != 5) {
		std::cout << "Incorrect argument count.Usage:" << std::endl
			<< "Particle Count" << std::endl
			<< "Ticks" << std::endl
			<< "MeshSize" << std::endl
			<< "Nearest Neighbors" << std::endl;
	}

	//inputs
	uint initialParticleCount = strtoumax(argv[1], NULL, 10);
	uint simulationTicks = strtoumax(argv[2], NULL, 10);
	uint meshSize = strtoumax(argv[3], NULL, 10); //must be a power of two
	uint nearestNeighbors = strtoumax(argv[4], NULL, 10);

	//error checking on inputs
	if ((initialParticleCount == UINTMAX_MAX && errno == ERANGE) || initialParticleCount < 0) {
		fprintf(stderr, "Incorrect particle count paramenter.\n");
	}
	if ((simulationTicks == UINTMAX_MAX && errno == ERANGE) || simulationTicks < 0) {
		fprintf(stderr, "Incorrect ticks paramenter.\n");
	}
	if ((meshSize == UINTMAX_MAX && errno == ERANGE) || meshSize < 0 || IsPower2(meshSize)) {
		fprintf(stderr, "Incorrect mesh size paramenter. (Must be power of 2)\n");
	}
	if ((nearestNeighbors == UINTMAX_MAX && errno == ERANGE) || nearestNeighbors < 1) {
		fprintf(stderr, "Incorrect nearest neighbors paramenter.\n");
	}


	//init MPI
	int rankCount, ID;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &rankCount);
	MPI_Comm_rank(MPI_COMM_WORLD, &ID);


	//////////////////////////
	/*Setup Rank Information*/
	//////////////////////////
	uint currentParticleCount;
	uint particleOffset; //the offset of the local particles into the global particle count
	uint particlestoSimulate = 0; //initial count needed =0
	Particle* particles;  //the local array of particles being simulated

	ParticlestoSimulate(ID, rankCount, initialParticleCount, particles, particlestoSimulate, particleOffset);


	/////////////////////////
	/*Initialize Simulation*/
	/////////////////////////

	for (int i = 0; i < particlestoSimulate; ++i) {
		InitParticle(particles[i], particleOffset + i, initialParticleCount);
	}


	////////////////////
	/*Start Simulation*/
	////////////////////

	currentParticleCount = initialParticleCount;

	for (int i = 0; i < simulationTicks; ++i) {
		//TODO: move particles (e.g. Integration. Please do not use Euler. Do like some Verlet integration as a minimum)

		//TODO: create acceleration structure (or send all particles to all other ranks, so inefficient)
		//IF ACCELERATION:	Assign morton code to each particle
		//					sort all the morton codes in parallel (radix sort is fastest and really easy)
		//					k nearest neighbors will be the particles above and below in the array

		//TODO: update particles using k nearest neighbors

		//TODO: create and remove particles

		//update rank information
		ParticlestoSimulate(ID, rankCount, currentParticleCount, particles, particlestoSimulate, particleOffset);

	}
	//////////////////
	/*End Simulation*/
	//////////////////

	//TODO: create mesh (probably icosphere) indices/vertices

	//TODO: update mesh heights based on particles

	//TODO: write mesh to .obj format. meshSize will be the number of faces (icosphere is power of two)

	MPI_Finalize();

	return 0;
}
