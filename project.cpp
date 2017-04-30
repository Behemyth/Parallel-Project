#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <cerrno>
#include <inttypes.h>
#include <limits>
#include <map>
#include <vector>
#include <iomanip> 
#include <sstream> 
#include <iostream>
#include <algorithm>
#include <functional>
#include <queue>
#include <cassert>
#include "clcg4.h"
#include <algorithm>

/////////////////
/*Global Values*/
/////////////////

#define uint unsigned int
#define byte uint8_t
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

	uint plateID;
	uint currentRank;
	uint64_t mortonCode;

	bool operator < (const Particle& p) const
	{
		return (mortonCode < p.mortonCode);
	}

	bool operator > (const Particle& p) const
	{
		return (mortonCode > p.mortonCode);
	}

};

class Vertex {

public:
	float x;
	float y;
	float z;

	Vertex(float x_, float y_, float z_) {
		x = x_;
		y = y_;
		z = z_;
	}

	// A default constructor is needed to call resize()
	// Otherwise, this should not be used
	Vertex() {
		x = 0;
		y = 0;
		z = 0;
	}
};

class Face {
public:
	int v1;
	int v2;
	int v3;

	Face(int v1_, int v2_, int v3_) {
		v1 = v1_;
		v2 = v2_;
		v3 = v3_;
	}

	// A default constructor is needed to call resize()
	// Otherwise, this should not be used
	Face() {
		v1 = 1;
		v2 = 2;
		v3 = 3;
	}
};

class Plate {
public:
	std::vector<Particle> particles;
	int plateID;
};

///////////////
/*Global Data*/
///////////////

//Holds the vertices of the final icosphere object
std::vector<Vertex> vertices;
//Caches the vertices in a time efficient manner for easy retrieval
//  Maps ID based on the parents points to location in vertices
std::map<uint64_t, int> vertexCache;
//Holds the faces of the final icosphere object
std::vector<Face> faces;


////////////////////////////////////
/*Functional Functions Doing Things*/
////////////////////////////////////


//compact the bits
inline uint64_t CompactBits(const uint a) {

	uint64_t x = a;
	x = x & 0x1fffff;
	x = (x | x << 32) & 0x1f00000000ffff;
	x = (x | x << 16) & 0x1f0000ff0000ff;
	x = (x | x << 8) & 0x100f00f00f00f00f;
	x = (x | x << 4) & 0x10c30c30c30c30c3;
	x = (x | x << 2) & 0x1249249249249249;
	return x;
}


//interleave the bits
inline uint64_t CalculateMortonCode(const uint x, const uint y, const uint z) {
	return CompactBits(x) | (CompactBits(y) << 1) | (CompactBits(z) << 2);
}

//unit box [-1,1]
uint64_t CreateCode(float x, float y, float z) {

	//map the floating point values to .. [0,UINT_MAX]
	int max = std::numeric_limits<int>::max();
	uint xMap = int(max*x) + max;
	uint yMap = int(max*y) + max;
	uint zMap = int(max*z) + max;

	return CalculateMortonCode(xMap, yMap, zMap);
}

/**
* Checks to see if an uint is a power of two
*
* @param x - The integer to check
* @return bool - Is it a power of two
*/
bool IsPower2(uint x) {

	return x && !(x & (x - 1));

}

/**
* Checks to see if an unit is a legal number of
* faces for an isophere
*
* @param x - The integer to check
* @return bool - Is it a legal number of faces
*/
bool IsLegalIcosphereFaceNumber(uint x) {
	x = x / 20;
	if (x == 0)
		return false;
	while (x != 1)
	{
		if (x % 4 != 0)
			return false;
		x = x / 4;
	}
	return true;
}

/**
* Determines how may levels of recursion will be necessary to get the desired
* number of faces on the icosphere.
* Should only be used after ensuring x is legal with IsLegalIcosphereFaceNumber
*
* @param x - the number of faces in the desired icosphere
* @return the levels of recursion necessary
*/
int IcosphereLevel(uint x) {
	x = x / 20;
	return log(x) / log(4);
}

/**
* Determine the number of faces on the mesh based on the given level of sphere recursion
*
* @param sphereLevel - the levels of recursion when generating the icosphere
* @return - the mesh size given that number of recursions
*/
int meshSize(uint sphereLevel) {
	int mesh = pow(sphereLevel, 4);
	return (mesh * 20);
}


/**
*Odd-Even Parallel sort
* Sorts the particles by Morton ID, after calculating the Morton ID
*the global particles will be updated by this function, will cordinate with all other ranks
*
* @param data - The global particles to sort by Morton ID
* @param size - The particles size
* @param localOffset - the global offset into the local simulation
* @param localSize - the local size to simulate
*/
void Sort(std::vector<Particle>& data, uint size, uint localOffset, uint localSize, uint rankCount, uint rankID) {

	//update all the local paricles morton codes
	for (int i = 0; i < localSize; ++i) {
		data[i + localOffset].mortonCode = CreateCode(data[i + localOffset].x, data[i + localOffset].y, data[i + localOffset].z);
	}

	//sort the local array
	std::sort(data.begin() + localOffset, data.begin() + localOffset + localSize);

	uint inSize = sizeof(Particle)*size;
	uint outSize = sizeof(Particle)*localSize;


	int *byteCount = new int[rankCount];
	int *count = new int[rankCount];
	int *countLeft = new int[rankCount];
	int *displacementBytes = new int[rankCount];
	int *displacement = new int[rankCount];

	//calc the rank statistics for variable particle count
	uint localCount = size / rankCount;
	uint remainder = size % rankCount;
	uint localCountBytes = localCount * sizeof(Particle);

	for (int i = 0; i < rankCount; i++) {
		if (i < remainder) {
			byteCount[i] = localCountBytes + sizeof(Particle);
			count[i] = localCount + 1;
		}
		else {
			byteCount[i] = localCountBytes;
			count[i] = localCount;
		}
		countLeft[i] = count[i];
	}

	displacementBytes[0] = 0;
	displacement[0] = 0;
	for (int i = 1; i < rankCount; i++) {
		displacementBytes[i] = displacementBytes[i - 1] + byteCount[i - 1];
		displacement[i] = displacement[i - 1] + count[i - 1];
	}

	/*assert(localSize * sizeof(Particle) == byteCount[rankID]);
	assert(localOffset * sizeof(Particle) == displacementBytes[rankID]);*/

	//gather all rank mortons on to rank 0
	MPI_Gatherv(data.data() + localOffset, outSize, MPI_BYTE, data.data(), byteCount, displacementBytes, MPI_BYTE, 0, MPI_COMM_WORLD);


	if (rankID == 0) {
		//merge all sorted arrays with min-heap sort
		std::priority_queue<Particle, std::vector<Particle>, std::greater<Particle> > least;
		std::vector<Particle> finalData(data.size());

		uint itr = 0;

		//add least from each bin
		for (int i = 0; i < rankCount; i++) {
			if (count[i] > 0) {
				least.push(data[displacement[i]]);
				countLeft[i]--;
			}
		}

		while (itr != size) {
			Particle temp = least.top();
			least.pop();
			finalData[itr] = temp;
			if (countLeft[temp.currentRank] > 0) {
				least.push(data[displacement[temp.currentRank] + (count[temp.currentRank] - countLeft[temp.currentRank])]);
				countLeft[temp.currentRank]--;
			}

			++itr;

		}

		//copy sorted
		data = finalData;
	}

	//scatter new data
	MPI_Bcast(data.data(), inSize, MPI_BYTE, 0, MPI_COMM_WORLD);

	delete byteCount;
	delete displacementBytes;
	delete displacement;
	delete count;
	delete countLeft;

}

/**
* returns the amount of particles to simulate for the current rank.
*
* @param rankID - The current rank
* @param rankCount - The total rank count
* @param particleCount - The particle count to distribute
* @return count - The number of particles particles to simulate
* @return index - The index of the first particle for this rank
*/
inline void ParticlestoSimulate(const uint rankID, const uint rankCount, const uint particleCount, uint& count, uint& index) {

	count = particleCount / rankCount;
	uint remainder = particleCount % rankCount;

	if (rankID < remainder) {

		++count;
		index = count*rankID;

	}
	else {

		index = (count*rankID) + remainder;

	}

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

	particle.x = cos(latitude) * cos(longitude);
	particle.y = cos(latitude) * sin(longitude);
	particle.z = sin(latitude);

	particle.plateID = particleID;
	particle.height = 0.3f;
	particle.mortonCode = 0;
}

/**
* Add a vertex to a list of vertices, adjusting it so it sits on the unit sphere
*
* @param v - The vertex to be added
* @param vertices - The list of vertices to add to
* @return the location of the new vertex
*/
int addVertex(Vertex v)
{
	double length = sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
	vertices.push_back(Vertex(v.x / length, v.y / length, v.z / length));
	return vertices.size() - 1;
}

/**
* Finds the midpoint between two points
* Used to form the next segment of the isosphere
*
* @param v1 - The location of the first vertex in vertices
* @param v2 - The location of the second vertex in vertices
* @return the location in vertices of the new vertex
*/
int findMidpoint(int v1, int v2)
{
	// determines if the midpoint already exists, by first
	// getting the unique ID of these two vertices and then
	// checking the cache
	uint64_t smaller = (v1 < v2) ? v1 : v2;
	uint64_t greater = (v1 < v2) ? v2 : v1;
	uint64_t key = (smaller << 32) + greater;

	if (vertexCache.count(key) > 0)
	{
		return vertexCache[key];
	}

	// if this midpoint does not exist yet, calculuate its location
	Vertex vertex1 = vertices[v1];
	Vertex vertex2 = vertices[v2];
	Vertex midpoint = Vertex((vertex1.x + vertex2.x) / 2.0,
		(vertex1.y + vertex2.y) / 2.0,
		(vertex1.z + vertex2.z) / 2.0);
	int newV = addVertex(midpoint);

	// store the new midpoint in the cache
	vertexCache[key] = newV;
	return newV;
}

/**
* Return the k particles closest to this one, not including the current particle
*
* @param k - the amount of neighbors to return
* @param particles - the global particle array
* @param size - current size of the particle array
* @param position - the index of the particle to check
* @return the nearest neighbors
*/
std::vector<Particle> getNearestNeighbors(uint k, std::vector<Particle>& particles, uint size, uint position) {
	std::priority_queue<Particle, std::vector<Particle>, std::less<Particle> > least;
	int temp = 0;
	for (int test = position; temp < k; temp++) {

		//offset by 1
		if (test - temp - 1 >= 0) {
			least.push(particles[test - temp - 1]);
		}

		//offset by 1 
		if (test + temp + 1 < size) {
			least.push(particles[test + temp + 1]);
		}
	}
	while (least.size() > k) {
		least.pop();
	}


	std::vector<Particle> parts;

	while (least.size() > 0) {
		parts.push_back(least.top());
		least.pop();
	}

	return parts;
}

/**
* Return the k particles closest to an arbitrary point in space
*
* @param k - the amount of neighbors to return
* @param particles - the global particl array
* @param size - current size of the particle array
* @param mortonCode - the code of the particle to check
* @return the nearest neighbors
*/
std::vector<Particle> getNearestNeighbors(uint k, std::vector<Particle>& particles, uint size, uint64_t mortonCode) {
	std::priority_queue<Particle, std::vector<Particle>, std::less<Particle> > least;
	int temp = 0;

	//binary search
	uint left = 0;
	uint right = size - 1;
	uint position = 0;

	while (left <= right) {

		int middle = (left + right) / 2;

		if (left == right - 1) {
			position = left;
			break;
		}
		else if (particles[middle].mortonCode > mortonCode) {
			right = middle - 1;
		}
		else {
			left = middle + 1;
		}

	}

	for (int test = position; temp < k; temp++) {

		//not offset because morton takes the lowest bound
		if (test - temp >= 0) {
			least.push(particles[test - temp]);
		}

		//offset by 1 
		if (test + temp + 1 < size) {
			least.push(particles[test + temp + 1]);
		}
	}
	while (least.size() > k) {
		least.pop();
	}


	std::vector<Particle> parts;

	while (least.size() > 0) {
		parts.push_back(least.top());
		least.pop();
	}
	return parts;
}

int getSmallestPlate(std::map<int, std::vector<int>>& plates) {
	std::map<int, std::vector<int>>::iterator itr = plates.begin();
	int size = itr->second.size();
	int id = itr->first;
	for( ; itr != plates.end(); itr++) {
		if(itr->second.size() < size) {
			size = itr->second.size();
			id = itr->first;
		}
	}
	return id;
}

void printVector(std::vector<int> v) {
	for(int i = 0; i < v.size(); i++) {
		std::cout << v[i] << "  ";
	}
}

void printMap(std::map<int, std::vector<int>> m) {
	std::map<int, std::vector<int>>::iterator itr = m.begin();
	for(; itr != m.end(); itr++) {
		std::cout << itr->first << ":  ";
		printVector(itr->second);
		std::cout << std::endl;
	}
}

/**
* Return the particles weights based on squared distance
*
* @param nearest - the particles to weight
* @param x,y,z - the position in space
* @return the weights [0,1] for each particle
*/
std::vector<float> getWeights(std::vector<Particle>& nearest, float x, float y, float z) {
	std::vector<float> weights(nearest.size());

	float furthest = 0.0f;
	for (int i = 0; i < nearest.size(); i++) {

		float xSqr = (x - nearest[i].x) * (x - nearest[i].x);
		float ySqr = (y - nearest[i].y) * (y - nearest[i].y);
		float zSqr = (z - nearest[i].z) * (z - nearest[i].z);

		float dist2 = xSqr + ySqr + zSqr;

		if (furthest < dist2) {
			furthest = dist2;
		}

		weights[i] = dist2;
	}

	for (int i = 0; i < nearest.size(); i++) {
		weights[i] = (1.0f - (weights[i] / furthest)) / nearest.size();
	}

	return weights;
}

//////////////////////
/*Mainly the Program*/
//////////////////////

int main(int argc, char **argv)
{

	////////////////////////////
	/*Setup Global Information*/
	////////////////////////////

	std::vector<Particle> particles;
	std::map<int, Plate> plates;


	/********** Initialize MPI **********/
	int rankCount, ID;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &rankCount);
	MPI_Comm_rank(MPI_COMM_WORLD, &ID);

	// Init 16,384 RNG streams - each rank has an independent stream
	InitDefault();
	MPI_Barrier(MPI_COMM_WORLD);

	if (ID == 0) {
		//input gets called by all mpi ranks anywho
		if (argc != 6) {
			std::cout << "Incorrect argument count.Usage:" << std::endl
				<< "Particle Count" << std::endl
				<< "Ticks" << std::endl
				<< "Sphere Level" << std::endl
				<< "Nearest Neighbors" << std::endl
				<< "Number of Plates" << std::endl;
		}
	}

	//inputs done for all ranks
	uint initialParticleCount = strtoumax(argv[1], NULL, 10);
	uint simulationTicks = strtoumax(argv[2], NULL, 10);
	uint sphereLevel = strtoumax(argv[3], NULL, 10);
	uint nearestNeighbors = strtoumax(argv[4], NULL, 10);
	uint numberOfPlates = strtoumax(argv[5], NULL, 10);
	double startTime;

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
		if ((sphereLevel == UINTMAX_MAX && errno == ERANGE) || sphereLevel < 0) {
			std::cout << "Incorrect sphere level paramenter." << std::endl;
			return 1;
		}
		if ((nearestNeighbors == UINTMAX_MAX && errno == ERANGE) || nearestNeighbors < 1) {
			std::cout << "Incorrect nearest neighbors paramenter." << std::endl;
			return 1;
		}
		if ((numberOfPlates == UINTMAX_MAX && errno == ERANGE) || numberOfPlates < 1) {
			std::cout << "Incorrect number of plates paramenter." << std::endl;
			return 1;
		}

		//start time
		if (ID == 0) {
			startTime = MPI_Wtime();
		}
	}

	//create the particle datatype

	MPI_Barrier(MPI_COMM_WORLD);

	//////////////////////////
	/*Setup Rank Information*/
	//////////////////////////
	uint currentParticleCount;
	uint particleOffset = 0; //the offset of the local particles into the global particle count
	uint particlestoSimulate = 0;

	ParticlestoSimulate(ID, rankCount, initialParticleCount, particlestoSimulate, particleOffset);


	/////////////////////////
	/*Initialize Simulation*/
	/////////////////////////
	//(Plate ID, temp, and any gosh darn variable this simulation would be cool with)

	particles.resize(initialParticleCount);

	//init particles sets the ID to the global Particle ID
	for (int i = 0; i < particlestoSimulate; ++i) {
		InitParticle(particles[particleOffset + i], particleOffset + i, initialParticleCount);
		particles[particleOffset + i].currentRank = ID;
	}


	//Sort data (updates the global array)
	Sort(particles, initialParticleCount, particleOffset, particlestoSimulate, rankCount, ID);

	//ACTUAL plate assigning
	if(ID == 0) {
		//First, get vectors for every plate
		std::map<int, std::vector<int>> plateByParticleId; // maps plateID to particles (in term of their location in particles)
		for (int p = 0; p < particles.size(); p++) {
			Particle particle = particles[p];
			plateByParticleId[particle.plateID].push_back(p);
		}
		while(plateByParticleId.size() > numberOfPlates) {
			int smallestPlateID = getSmallestPlate(plateByParticleId);
			std::vector<int> smallestPlate = plateByParticleId[smallestPlateID];
			int closestPlate;
			for(int i = 0; i < smallestPlate.size(); i++) {
				Particle p = particles[smallestPlate[i]];
				std::vector<Particle> neighbors = getNearestNeighbors(nearestNeighbors, particles, initialParticleCount, (uint)smallestPlate[i]);
				std::random_shuffle(neighbors.begin(), neighbors.end()); //ensure the nearest neighbor is randomly chosen
				closestPlate = -1;
				for(int n = 0; n < neighbors.size(); n++) {
					Particle neighbor = neighbors[n];
					if(neighbor.plateID != p.plateID) {
						closestPlate = neighbor.plateID;
						break;
					}
				}
				if(closestPlate > -1) {
					break;
				}
			}
			if(closestPlate == -1) {
				std::cout << "Error: No nearest neighbor" << std::endl;
			}
			//std::vector<int> parentPlate = plateByParticleId[closestPlate];
			plateByParticleId[closestPlate].insert(plateByParticleId[closestPlate].end(), plateByParticleId[smallestPlateID].begin(), plateByParticleId[smallestPlateID].end());
			for(int i = 0; i < plateByParticleId[closestPlate].size(); i++) {
				particles[plateByParticleId[closestPlate][i]].plateID = closestPlate;
			}
			plateByParticleId.erase(smallestPlateID);
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	// Send out the new particle and plate vectors
	int p_size;
	if (ID == 0) {
		p_size = particles.size() * sizeof(Particle);
	}
	// First, send out the amount of data that will be sent
	MPI_Bcast(&p_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
	// Get space ready for the broadcast
	if (ID != 0) {
		particles.resize(p_size / sizeof(Particle));
	}
	// Then, send out the data
	MPI_Bcast(&particles[0], p_size, MPI_BYTE, 0, MPI_COMM_WORLD);

	//Fill the plates map based on the recieved particles
	for(int p = 0; p < particles.size(); p++) {
		Particle particle = particles[p];
		plates[particle.plateID].plateID = particle.plateID;
	}

	/*//DEBUG:
	if(ID == 0) {
		std::cout << plates.size() << std::endl;
		for(std::map<int, Plate>::iterator itr = plates.begin(); itr != plates.end(); itr++) {
			Plate plate = itr->second;
			std::cout << plate.plateID << std::endl;
		}
	}*/

	////////////////////
	/*Start Simulation*/
	////////////////////

	currentParticleCount = initialParticleCount;

	for (int i = 0; i < simulationTicks; ++i) {
		//TODO: move particles (e.g. Integration. Please do not use Euler. Do like some Verlet integration as a minimum)


		/////////////////////////////////////
		/*Create the Acceleration Structure*/
		/////////////////////////////////////

		//sort all the particles in the system by morton code
		Sort(particles, currentParticleCount, particleOffset, particlestoSimulate, rankCount, ID);
		//Now ok to call KNearest for this timestep

		//TODO: update particles using k nearest neighbors

		//TODO: create and remove particles

		//update rank information
		ParticlestoSimulate(ID, rankCount, currentParticleCount, particlestoSimulate, particleOffset);

		//update the current rank processor
		for (int i = 0; i < particlestoSimulate; ++i) {
			particles[particleOffset + i].currentRank = ID;
		}
	}
	//////////////////
	/*End Simulation*/
	//////////////////

	MPI_Barrier(MPI_COMM_WORLD);


	//TODO: update mesh heights based on K-N particles

////////////////////////////////
/*Create the initial icosphere*/
////////////////////////////////

	MPI_File file;
	MPI_Status status;

	if (ID == 0) {
		////////////////////////////////
		/*Create the initial icosphere*/
		////////////////////////////////

		int t = (1.0 + sqrt(5.0)) / 2.0;

		// Create the initial vertices of the icohedron
		addVertex(Vertex(-1, t, 0));
		addVertex(Vertex(1, t, 0));
		addVertex(Vertex(-1, -t, 0));
		addVertex(Vertex(1, -t, 0));

		addVertex(Vertex(0, -1, t));
		addVertex(Vertex(0, 1, t));
		addVertex(Vertex(0, -1, -t));
		addVertex(Vertex(0, 1, -t));

		addVertex(Vertex(t, 0, -1));
		addVertex(Vertex(t, 0, 1));
		addVertex(Vertex(-t, 0, -1));
		addVertex(Vertex(-t, 0, 1));

		// Create the initial faces of the icohedron
		// Faces around point 0
		faces.push_back(Face(0, 11, 5));
		faces.push_back(Face(0, 5, 1));
		faces.push_back(Face(0, 1, 7));
		faces.push_back(Face(0, 7, 10));
		faces.push_back(Face(0, 10, 11));

		// Adjacent faces 
		faces.push_back(Face(1, 5, 9));
		faces.push_back(Face(5, 11, 4));
		faces.push_back(Face(11, 10, 2));
		faces.push_back(Face(10, 7, 6));
		faces.push_back(Face(7, 1, 8));

		// Faces around point 3
		faces.push_back(Face(3, 9, 4));
		faces.push_back(Face(3, 4, 2));
		faces.push_back(Face(3, 2, 6));
		faces.push_back(Face(3, 6, 8));
		faces.push_back(Face(3, 8, 9));

		// Adjacent faces 
		faces.push_back(Face(4, 9, 5));
		faces.push_back(Face(2, 4, 11));
		faces.push_back(Face(6, 2, 10));
		faces.push_back(Face(8, 6, 7));
		faces.push_back(Face(9, 8, 1));

		// Now, begin splitting up the triangle faces to form an isohedron of a desired number of faces
		int levels = sphereLevel;

		for (int i = 0; i < levels; i++)
		{
			std::vector<Face> newFaces;
			for (int j = 0; j < faces.size(); j++)
			{
				Face face = faces[j];
				// Split this face into four new faces
				int a = findMidpoint(face.v1, face.v2);
				int b = findMidpoint(face.v2, face.v3);
				int c = findMidpoint(face.v3, face.v1);

				newFaces.push_back(Face(face.v1, a, c));
				newFaces.push_back(Face(face.v2, b, a));
				newFaces.push_back(Face(face.v3, c, b));
				newFaces.push_back(Face(a, b, c));
			}
			faces = newFaces;
		}
	}
	// Send the newly generated vertex and face arrays to the other ranks for writing
	int v_size;
	int f_size;
	if (ID == 0) {
		v_size = vertices.size() * sizeof(Vertex);
		f_size = faces.size() * sizeof(Face);
	}

	// First, send out the amount of data that will be sent
	MPI_Bcast(&v_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&f_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
	// Get space ready for the broadcast
	if (ID != 0) {
		vertices.resize(v_size / sizeof(Vertex));
		faces.resize(f_size / sizeof(Face));
	}

	// Then, send out the data
	MPI_Bcast(&vertices[0], v_size, MPI_BYTE,
		0, MPI_COMM_WORLD);
	MPI_Bcast(&faces[0], f_size, MPI_BYTE,
		0, MPI_COMM_WORLD);


	int verticesToWrite = vertices.size() / rankCount;
	int extraVertices = vertices.size() % rankCount;
	if (ID < extraVertices) {
		verticesToWrite += 1;
	}
	int facesToWrite = faces.size() / rankCount;
	int extraFaces = faces.size() % rankCount;
	if (ID < extraFaces) {
		facesToWrite += 1;
	}

	// Write the vertices to the obj file

	// Open the file, deleting it if it exists already
	int exists = MPI_File_open(MPI_COMM_WORLD,
		(char *)"planet.obj",
		MPI_MODE_CREATE | MPI_MODE_EXCL | MPI_MODE_WRONLY,
		MPI_INFO_NULL,
		&file);
	if (exists != MPI_SUCCESS) {
		if (ID == 0) {
			MPI_File_delete((char *)"planet.obj", MPI_INFO_NULL);
		}
		MPI_File_open(MPI_COMM_WORLD,
			(char *)"planet.obj",
			MPI_MODE_CREATE | MPI_MODE_WRONLY,
			MPI_INFO_NULL,
			&file);
	}

	//	MPI_INFO_NULL, &file);

	std::stringstream stream;
	int start = ID * verticesToWrite + (ID >= extraVertices ? extraVertices : 0);
	int end = start + verticesToWrite;
	int vertexBytesPerLine = 1 				// 'v'
		+ (13 * 3) + 1 						// Numbers  
		+ 4 								// spaces
		+ 1;								// newline
	MPI_Offset offset = vertexBytesPerLine * start;
	MPI_File_seek(file, offset, MPI_SEEK_SET);
	// Debug call
	if (end > vertices.size()) {
		std::cout << "Error in getting write range (vertices)" << std::endl;
	}
	for (int v = start; v < end; v++) {

		//modify vertices based on nearest neighbors
		uint64_t mort = CreateCode(vertices[v].x, vertices[v].y, vertices[v].z);
		std::vector<Particle> nearest = getNearestNeighbors(nearestNeighbors, particles, currentParticleCount, mort);
		std::vector<float> weights = getWeights(nearest, vertices[v].x, vertices[v].y, vertices[v].z);

		float avgDelta = 0.0f;
		for (int l = 0; l < weights.size(); l++) {
			avgDelta += (weights[l] * (nearest[l].height)) / weights.size();
		}
		float avgHeight = avgDelta + 1.0f;

		vertices[v].x *= avgHeight; 
		vertices[v].y *= avgHeight; 
		vertices[v].z *= avgHeight;

		stream << std::fixed << std::setprecision(10) << std::setw(13) << vertices[v].x;
		std::string x = stream.str();
		stream.str(std::string());

		stream << std::fixed << std::setprecision(10) << std::setw(13) << vertices[v].y;
		std::string y = stream.str();
		stream.str(std::string());

		stream << std::fixed << std::setprecision(10) << std::setw(13) << vertices[v].z;
		std::string z = stream.str();
		stream.str(std::string());

		std::string line = "v " + x + " " + y + " " + z + " 1" + "\n";

		MPI_File_write(file, (void*)line.c_str(), line.size(), MPI_CHAR, &status);
	}
	stream.str(std::string());

	//Write the faces to the obj file
	start = ID * facesToWrite + (ID >= extraFaces ? extraFaces : 0);
	end = start + facesToWrite;
	int faceBytesPerLine = 1 			// 'v'
		+ (16 * 3) 						// Numbers  
		+ 3 							// spaces
		+ 1;							// newline
	offset = (faceBytesPerLine * start) + (vertexBytesPerLine * vertices.size());
	MPI_File_seek(file, offset, MPI_SEEK_SET);
	// Debug call
	if (end > faces.size()) {
		std::cout << "Error in getting write range (faces)" << std::endl;
	}
	for (int f = start; f < end; f++) {
		stream << std::setw(16) << (faces[f].v1 + 1);
		std::string v1 = stream.str();
		stream.str(std::string());

		stream << std::setw(16) << (faces[f].v2 + 1);
		std::string v2 = stream.str();
		stream.str(std::string());

		stream << std::setw(16) << (faces[f].v3 + 1);
		std::string v3 = stream.str();
		stream.str(std::string());

		std::string line = "f " + v1 + " " + v2 + " " + v3 + "\n";

		MPI_File_write(file, (void*)line.c_str(), line.size(), MPI_CHAR, &status);
	}

	MPI_File_close(&file);

	//end time
	if (ID == 0) {
		double endTime = MPI_Wtime();

		double calcTime = endTime - startTime;

		printf("The simulation took %f seconds.\n", calcTime);
	}

	MPI_Finalize();

	return 0;
}
