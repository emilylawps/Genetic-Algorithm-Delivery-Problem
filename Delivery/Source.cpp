#include <iostream>
#include <ctime>
#include <algorithm>
#include <fstream>
using namespace std;

const int POPSIZE = 10;
const int GENENUM = 5;
const int DISTANCES[6][6] =	{ 
	{ 0, 1, 2, 5, 2, 4 },
	{ 1, 0, 2, 4, 1, 3 },
	{ 2, 2, 0, 1, 4, 3 },
	{ 5, 4, 1, 0, 5, 2 },
	{ 2, 1, 4, 5, 0, 3 },
	{ 4, 3, 3, 2, 3, 0 }
};
const int PARENTSNUM = 2;
const double CROSS = 0.9;
const double MUTATE = 0.9;
const int SURVIVAL = 2;
const int GENERATION = 20;

int chromosome[POPSIZE][GENENUM];
int totalDist[POPSIZE];
int accDist;
double fitness[POPSIZE];
int parents[PARENTSNUM][GENENUM];
int parentCands[PARENTSNUM][2];
int parentIndex[PARENTSNUM];
int p1, p2;
int child[PARENTSNUM][GENENUM];
int newChromo[POPSIZE][GENENUM];
double survivalFitness[SURVIVAL];
int bestChromosome[GENENUM];
double aveFitness;
double bestFitness;
ofstream avefit;
ofstream bestfit;
ofstream bestchromo;


void initializePopulation() {

	int x[GENENUM] = { 2, 3, 4, 5, 6 };		//city 1 is always the first city, hence no need to include in chromosome
	int temp;

	for (int i = 0; i < POPSIZE; i++) {
		for (int k = 1; k < GENENUM; k++) {
			int r = rand() % k;				//select random index
			temp = x[k];					//swap chromosome indexes to generate random permutation chromosome arrays
			x[k] = x[r];
			x[r] = temp;
		}
		for (int j = 0; j < GENENUM; j++) {
			chromosome[i][j] = x[j];
		}
	}
}

void evaluateChromosome() {

	for (int i = 0; i < POPSIZE; i++) {
		accDist = 0;
		accDist += DISTANCES[0][chromosome[i][0] - 1];		//calculate the distance between city 1 with the first city in chromosome
		for (int j = 0; j < GENENUM - 1; j++) {
			accDist += DISTANCES[chromosome[i][j] - 1][chromosome[i][j + 1] - 1];		//calculate the distance between each cities in the chromosome
		}
		totalDist[i] = accDist;
		fitness[i] = 1.0 / totalDist[i];
	}
}

void printFitness() {
	cout << "Chromosomes\tTotal Distance\t  Fitness\n";
	for (int i = 0; i < POPSIZE; i++) {

		for (int j = 0; j < GENENUM; j++) {
			cout << chromosome[i][j] << " ";
		}
		cout << "\t  " << totalDist[i] << "\t\t  " << fitness[i] << endl;
	}
	cout << endl;
}

void printPopulation() {

	for (int i = 0; i < POPSIZE; i++) {
		cout << "\tChromosome " << i + 1 << ": ";
		for (int j = 0; j < GENENUM; j++) {
			cout << chromosome[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}

void parentSelection() {
	int i;

	do {										//loop if the parents selected are the same
		i = 0;
		while (i < PARENTSNUM) {				//loop while parents < 2
			do {								//loop if numbers are the same
				p1 = rand() % POPSIZE;
				p2 = rand() % POPSIZE;
			} while (p1 == p2);

			if (fitness[p1] >= fitness[p2]) {				//select the chromosome with higher fitness value as one of the parents
				for (int j = 0; j < GENENUM; j++) {
					parents[i][j] = chromosome[p1][j];
				}
				parentIndex[i] = p1;
			}
			else {
				for (int j = 0; j < GENENUM; j++) {
					parents[i][j] = chromosome[p2][j];
				}
				parentIndex[i] = p2;
			}
			parentCands[i][0] = p1;			//copy the index of parents candidates
			parentCands[i][1] = p2;
			i++;
		}
	} while (parentIndex[0] == parentIndex[1]);

	cout << "--------------------PARENT SELECTION--------------------" << endl;
	for (int i = 0; i < PARENTSNUM; i++) {
		cout << "Parent candicate pair " << i + 1 << ": " << parentCands[i][0] + 1 << " " << parentCands[i][1] + 1 << endl;
	}
	cout << endl;
}

void printParent() {
	
	cout << "Selected parents: " << endl;

	for (int i = 0; i < PARENTSNUM; i++) {
		cout << "\tChromosome " << parentIndex[i] + 1 << ": ";
		for (int j = 0; j < GENENUM; j++) {
			cout << parents[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}

void crossover() {
	cout << "------------------CROSSOVER & MUTATION------------------" << endl;
	float point = float(rand() % 11) / 10;
	cout << "Crossover Probability: " << CROSS << "\t" << "Random number: " << point << endl;

	if (point < CROSS) {
		//random 2 numbers for starting and ending index
		int cutting1, cutting2, start, end;
		do {
			cutting1 = rand() % GENENUM;
			cutting2 = rand() % GENENUM;
		} while (cutting1 == cutting2);

		//set min index as starting and max index as ending
		start = min(cutting1, cutting2);
		end = max(cutting1, cutting2);

		cout << "Crossover starting and ending point: " << start + 1 << " & " << end + 1 << endl;

		//copy the selected genes into child
		for (int i = 0; i < PARENTSNUM; i++) {
			
			for (int x = 0; x < GENENUM; x++) {
				child[i][x] = 0;		//make sure the child chromosomes does not store the previous value to avoid possible error
			}
			for (int j = start; j <= end; j++) {
				child[i][j] = parents[i][j];		//copy selected genes
			}
		}

		int p;
		//if ending is not the last gene, start copying parent genes after the ending index
		if (end < GENENUM - 1) {

			//compare the parent genes with selected gene in child to avoid repetition of gene
			int s = 1;
			for (int r = 0; r < PARENTSNUM; r++) {
				int q = end + 1;	//to count the index of child
				for (int i = end + 1; i < GENENUM; i++) {
					p = start;

					while (parents[r][i] != child[s][p] && p < GENENUM) {		//while current parent gene not same with the child and not reach the last index, keep comparing all child gene
						p++;
					}
					if (p >= GENENUM) {			//if parent genes != with all the present child genes (ady compared all of the genes till p>=GENENUM), copy the gene value to child 
						if (q < GENENUM) {		//if current child index not exceed last gene yet then continue copying gene to the current index
							child[s][q] = parents[r][i];
							q++;
						}
						else {		//update index of child to the first gene if current index exceed last gene (all genes after ending index had been filled up)
							q = 0;
							child[s][q] = parents[r][i];
							q++;
						}
					}
				}
				//when the gene in parent to copy has exceeded last gene, change index n copy from the first gene
				for (int j = 0; j < GENENUM; j++) {
					p = start;
					while (parents[r][j] != child[s][p] && p < GENENUM) {		//while current parent gene not same with the child and not reach the last index, keep comparing all child gene
						p++;
					}
					if (p >= GENENUM) {
						if (q < GENENUM) {
							child[s][q] = parents[r][j];
							q++;
						}
						else {
							q = 0;
							child[s][q] = parents[r][j];
							q++;
						}
					}
				}
				s--;	//repeat the same for another child;
			}

		}
		else {			//if ending index is the last gene
			int s = 1;
			for (int r = 0; r < PARENTSNUM; r++) {
				int q = 0;
				for (int i = 0; i < GENENUM; i++) {
					p = start;

					while (parents[r][i] != child[s][p] && p <= GENENUM) {		//while current parent gene not same with the child and not reach the ending index, keep comparing all child gene
						p++;
					}
					if (p >= GENENUM) {
						child[s][q] = parents[r][i];
						q++;
					}
				}
				s--;	//repeat the same for another child;
			}
		}
	}
	else {
		for (int i = 0; i < PARENTSNUM; i++) {
			for (int j = 0; j < GENENUM; j++) {
				child[i][j] = parents[i][j];
			}
		}
	}
	
}
//swap mutation
void mutation() {
	float point = float(rand() % 11) / 10;
	cout << "Mutation Probability: " << MUTATE << "\t" << "Random number: " << point << endl;

	if (point < MUTATE) {
		int swap1, swap2, temp;
		do {
			swap1 = rand() % GENENUM;
			swap2 = rand() % GENENUM;
		} while (swap1 == swap2);

		cout << "Mutation swapping points: " << swap1 + 1 << " & " << swap2 + 1 << endl;

		for (int i = 0; i < PARENTSNUM; i++) {
			temp = child[i][swap1];

			child[i][swap1] = child[i][swap2];
			child[i][swap2] = temp;
		}
	}
}

void printChild() {
	for (int i = 0; i < PARENTSNUM; i++) {
		cout << "\tChild " << i + 1 << ": ";
		for (int j = 0; j < GENENUM; j++) {
			cout << child[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}

void survivalSelection(int count) {

	for (int i = 0; i < SURVIVAL; i++) {
		accDist = 0;
		accDist += DISTANCES[0][child[i][0] - 1];
		for (int j = 0; j < GENENUM - 1; j++) {
			accDist += DISTANCES[child[i][j] - 1][child[i][j + 1] - 1];
		}
		totalDist[i] = accDist;
		survivalFitness[i] = 1.0 / totalDist[i];
	}

	if (survivalFitness[0] > survivalFitness[1]) {
		for (int j = 0; j < GENENUM; j++) {
			newChromo[count][j] = child[0][j];
		}
	}
	else {
		for (int j = 0; j < GENENUM; j++) {
			newChromo[count][j] = child[1][j];
		}
	}

	if (fitness[parentIndex[0]] > fitness[parentIndex[1]]) {
		for (int j = 0; j < GENENUM; j++) {
			newChromo[count + 1][j] = parents[0][j];
		}
	}
	else {
		for (int j = 0; j < GENENUM; j++) {
			newChromo[count + 1][j] = parents[1][j];
		}
	}
}

void printNewChromo(int count) {
	cout << count << " chromosomes in new generation: \n";
	for (int i = 0; i < count; i++) {
		cout << "\t";
		for (int j = 0; j < GENENUM; j++) {
			cout << newChromo[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}

void copyChromo() {
	for (int i = 0; i < POPSIZE; i++) {
		for (int j = 0; j < GENENUM; j++) {
			chromosome[i][j] = newChromo[i][j];
		}
	}
}

void calculateAverageFitness() {
	double totalFitness = 0;

	for (int i = 0; i < POPSIZE; i++) {
		totalFitness += fitness[i];
	}
	aveFitness = totalFitness / POPSIZE;
	cout << "Average fitness: " << aveFitness << endl;

	avefit << aveFitness << endl;
}

void recordBestFitness() {
	bestFitness = 0;
	for (int i = 0; i < POPSIZE; i++) {
		if (fitness[i] > bestFitness) {
			bestFitness = fitness[i];

			for (int j = 0; j < GENENUM; j++) {
				bestChromosome[j] = chromosome[i][j];
			}
		}
	}

	cout << "Best Chromosome: ";
	for (int i = 0; i < GENENUM; i++) {
		cout << bestChromosome[i] << " ";
		bestchromo << bestChromosome[i];
	}
	bestchromo << endl;
	cout << "\t Best Fitness: " << bestFitness << endl;
	bestfit << bestFitness << endl;

}

void main() {
	srand(time(NULL));

	avefit.open("Average Fitness.csv");
	bestchromo.open("Best Chromosome.csv");
	bestfit.open("Best Fitness.csv");
	
	int count;
	initializePopulation();

	for (int g = 0; g < GENERATION; g++) {
		count = 0;
		cout << "\n*****************************************************" << endl;
		cout << "Generation " << g+1 << " of chromosome: " << endl;

		printPopulation();
		evaluateChromosome();
		printFitness();

		while (count < POPSIZE) {
			parentSelection();
			printParent();
			crossover();
			printChild();
			mutation();
			printChild();
			survivalSelection(count);
			count += 2;
			printNewChromo(count);
		}
		copyChromo();
		calculateAverageFitness();
		recordBestFitness();
	}
	avefit.close();
	bestchromo.close();
	bestfit.close();
}

