#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <inttypes.h>
#include <limits>

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

class Particle {

public:
	//x,y,z position of the particle on a sphere Radius: 1
	float x;
	float y;
	float z;

	float height;	//delta from sea level

	bool removed;

	uint64_t mortonCode;

	void UpdateCode() {

		//map the floating point values to .. [0,UINT_MAX]
		int max = std::numeric_limits<int>::max();
		uint xMap = max*x + max;
		uint yMap = max*y + max;
		uint zMap = max*z + max;

		mortonCode = CalculateMortonCode(xMap, yMap, zMap);
	}

private:

	//compact the bits
	inline uint64_t CompactBits(const uint a) {
		const uint64_t masks[6] = {
			0x1fffff,
			0x1f00000000ffff,
			0x1f0000ff0000ff,
			0x100f00f00f00f00f,
			0x10c30c30c30c30c3,
			0x1249249249249249 };

		uint64_t x = a;
		x = x & masks[0];
		x = (x | x << 32) & masks[1];
		x = (x | x << 16) & masks[2];
		x = (x | x << 8)  & masks[3];
		x = (x | x << 4)  & masks[4];
		x = (x | x << 2)  & masks[5];
		return x;
	}


	//interleave the bits
	inline uint64_t CalculateMortonCode(const uint x, const uint y, const uint z) {
		return CompactBits(x) | (CompactBits(y) << 1) | (CompactBits(z) << 2);
	}

};




////////////////////////////////////
/*Functional Functions Doing Things*/
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

	Particle* temp = new Particle[count];

	if (particles) {
		//TODO: Copy over particles that are not removed into the new pointer
		//memcpy(temp, particles, initialCount*sizeof(struct Particle));         
		delete particles;
	}

	particles = temp;

}

/**
* Setup the Particle.
*
* @param particle - The particle obvs
*/

void InitParticle(Particle& particle, uint particleID, uint particleCount) {

	//calc lat and long in radians
	float longitude = GOLDEN_ANGLE_RADS*particleID;

	//TODO: is this needed?

	/*lon /= 2 * PI;
	lon -= floor(lon);
	lon *= 2 * PI;*/

	if (longitude > PI) {
		longitude -= 2 * PI;
	}

	float latitude = asin(-1 + 2 * particleID / (float)particleCount);

	particle.height = 0.0f;
	particle.removed = false;

	particle.x = cos(latitude) * cos(longitude);
	particle.y = cos(latitude) * sin(longitude);
	particle.z = sin(latitude);

}

//////////////////////
/*Mainly the Program*/
//////////////////////

int main(int argc, char **argv)
{

	////////////////////////////
	/*Setup Global Information*/
	////////////////////////////

	Particle* globalParticles;

	//init MPI
	int rankCount, ID;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &rankCount);
	MPI_Comm_rank(MPI_COMM_WORLD, &ID);

	//input gets called by all mpi ranks anywho
	if (ID == 0) {
		if (argc != 5) {
			std::cout << "Incorrect argument count.Usage:" << std::endl
				<< "Particle Count" << std::endl
				<< "Ticks" << std::endl
				<< "MeshSize" << std::endl
				<< "Nearest Neighbors" << std::endl;
		}
	}

	//inputs done for all ranks
	uint initialParticleCount = strtoumax(argv[1], NULL, 10);
	uint simulationTicks = strtoumax(argv[2], NULL, 10);
	uint meshSize = strtoumax(argv[3], NULL, 10); //must be a power of two
	uint nearestNeighbors = strtoumax(argv[4], NULL, 10);

	if (ID == 0) {
		//error checking on inputs
		if ((initialParticleCount == UINTMAX_MAX && errno == ERANGE) || initialParticleCount < 0) {
			std::cout << "Incorrect particle count paramenter." << std::endl;
			return 1;
		}
		if ((simulationTicks == UINTMAX_MAX && errno == ERANGE) || simulationTicks < 0) {
			std::cout << "Incorrect ticks paramenter." << std::endl;
			return 1;
		}
		if ((meshSize == UINTMAX_MAX && errno == ERANGE) || meshSize < 0 || !IsPower2(meshSize)) {
			std::cout << "Incorrect mesh size paramenter. (Must be power of 2)" << std::endl;
			return 1;
		}
		if ((nearestNeighbors == UINTMAX_MAX && errno == ERANGE) || nearestNeighbors < 1) {
			std::cout << "Incorrect nearest neighbors paramenter." << std::endl;
			return 1;
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);

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


		/////////////////////////////////////
		/*Create the Acceleration Structure*/
		/////////////////////////////////////

		//update all the paricles morton codes
		for (int i = 0; i < particlestoSimulate; ++i) {
			particles[i].UpdateCode();
		}

		//sort all the particles in the system by morton code
		//TODO: GLOBAl SORT
		//Now ok to call KNearest for this timestep

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
	MPI_File file;
	MPI_Status status;

	MPI_File_open(MPI_COMM_WORLD, (char *)"planet.obj", MPI_MODE_CREATE | MPI_MODE_WRONLY,
		MPI_INFO_NULL, &file);
	//MPI_File_seek(file, _, MPI_SEEK_SET);
	//MPI_File_write(file, _, _, _, &status);
	MPI_File_close(&file);



	MPI_Finalize();

	return 0;
}
