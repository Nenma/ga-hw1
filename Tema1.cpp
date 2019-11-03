// Tema1.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
#include <time.h>
#include <math.h>
#include <random>
#include <chrono>

using namespace std;

constexpr auto PI = 3.141592;
constexpr auto ITERATIONS = 30;
int PRECISION = 3; // default precision
typedef chrono::high_resolution_clock Clock;

template<typename T>
void printVector(vector<T> vector) {
	for (int i = 0; i < vector.size(); ++i) {
		cout << vector[i] << " ";
	}
	cout << endl;
}

// 1. Rastrigin's function (f : [-5.12, 5.12])
double Rastrigin(vector<double> x) {
	double res = 0;
	for (unsigned int i = 0; i < x.size(); ++i) {
		res += x[i] * x[i] - 10 * cos(2 * PI * x[i]);
	}
	res += (double)10 * x.size();
	return res;
}

// 2. De Jong's (Sphere) function (f : [-5.12, 5.12])
double Sphere(vector<double> x) {
	double res = 0;
	for (unsigned int i = 0; i < x.size(); ++i) {
		res += x[i] * x[i];
	}
	return res;
}

// 3. Schwefel's function (f : [-500, 500])
double Schwefel(vector<double> x) {
	double res = 0;
	for (unsigned int i = 0; i < x.size(); ++i) {
		res += -x[i] * sin(sqrt(abs(x[i])));
	}
	return res;
}

// 4. Michalewicz's function (f : [0, pi])
double Michalewicz(vector<double> x) {
	double res = 0;
	for (unsigned int i = 0; i < x.size(); ++i) {
		res += sin(x[i]) * pow(sin(((i + 1) * x[i] * x[i]) / PI), 20);
	}
	return -res;
}


double randomDouble(double lowerBound, double upperBound) {
	uniform_real_distribution<double> unif(lowerBound, upperBound);
	random_device rd;
	return unif(rd);
}

int randomInt(int lowerBound, int upperBound) {
	uniform_int_distribution<int> unif(lowerBound, upperBound);
	random_device rd;
	return unif(rd);
}

vector<bool> generateRandomBitstring(int size) {
	vector<bool> bitstring;
	for (int i = 0; i < size; ++i) {
		bitstring.push_back(randomInt(0, 1));
	}
	return bitstring;
}

int decimal(vector<bool> bitstring) {
	int x = 0;
	for (int i = 0; i < bitstring.size(); ++i) {
		x *= 2;
		x += bitstring[i];
	}
	return x;
}

vector<double> decode(vector<bool> bitstring, double lowerBound, double upperBound, int dim) {
	vector<double> candidate;
	int N = (upperBound - lowerBound) * pow(10, PRECISION);
	int size = ceil(log2(N));
	for (int i = 0; i < size * dim; i += size) {
		vector<bool> temporary;
		for (int j = 0; j < size; ++j) {
			temporary.push_back(bitstring[j + i]);
		}
		double component = lowerBound + decimal(temporary) * (upperBound - lowerBound) / (pow(2, size) - 1);
		candidate.push_back(component);
	}
	return candidate;
}

vector<bool> neighbour(vector<bool> bitstring, int position) {
	vector<bool> neighbour = bitstring;
	neighbour[position] = 1 - bitstring[position];
	return neighbour;
}

vector<bool> bestImprovement(vector<bool> bitstring, double f(vector<double>), double lowerBound, double upperBound, int dim) {
	vector<bool> best = bitstring;
	for (int i = 0; i < bitstring.size(); ++i) {
		vector<bool> neigh = neighbour(bitstring, i);
		if (f(decode(neigh, lowerBound, upperBound, dim)) < f(decode(best, lowerBound, upperBound, dim))) {
			best = neigh;
		}
	}
	return best;
}

vector<bool> firstImprovement(vector<bool> bitstring, double f(vector<double>), double lowerBound, double upperBound, int dim) {
	vector<bool> first = bitstring;
	for (int i = 0; i < bitstring.size(); ++i) {
		vector<bool> neigh = neighbour(bitstring, i);
		if (f(decode(neigh, lowerBound, upperBound, dim)) < f(decode(first, lowerBound, upperBound, dim))) {
			return neigh;
		}
	}
	return first;
}

void hillClimbing(vector<bool> improvement(vector<bool>, double(vector<double>), double, double, int), double f(vector<double>), double lowerBound, double upperBound, int dim) {
	double solMin = DBL_MAX;
	double solMax = DBL_MIN;
	double solAvg = 0;

	int timeMin = INT_MAX;
	int timeMax = INT_MIN;
	int timeAvg = 0;

	auto begin = Clock::now();

	for (int t = 0; t < ITERATIONS; ++t) {
		auto start = Clock::now();

		int size = ceil(log2((upperBound - lowerBound) * pow(10, PRECISION)));
		vector<bool> bitstring = generateRandomBitstring(size * dim);
		while (true) {

			vector<bool> neigh = improvement(bitstring, f, lowerBound, upperBound, dim);
			
			if (f(decode(neigh, lowerBound, upperBound, dim)) < f(decode(bitstring, lowerBound, upperBound, dim))) {
				bitstring = neigh;
			}
			else {
				auto finish = Clock::now();

				cout << t + 1 << ". ";
				printVector(decode(bitstring, lowerBound, upperBound, dim));
				double sol = f(decode(bitstring, lowerBound, upperBound, dim));
				int time = chrono::duration_cast<chrono::milliseconds>(finish - start).count();
				cout << sol << " " << time << "ms" << endl << endl;

				if (sol < solMin) solMin = sol;
				if (sol > solMax) solMax = sol;
				solAvg += sol;

				if (time < timeMin) timeMin = time;
				if (time > timeMax) timeMax = time;
				timeAvg += time;

				break;
			}
		}
	}

	auto end = Clock::now();

	std::cout << "Min: " << solMin << endl
		<< "Max: " << solMax << endl
		<< "Average: " << solAvg / ITERATIONS << endl;

	std::cout << endl;

	std::cout << "Min time: " << timeMin << "ms" << endl
		<< "Max time: " << timeMax << "ms" << endl
		<< "Average: " << timeAvg / ITERATIONS << "ms" << endl;

	std::cout << endl;

	std::cout << "Total time: " << chrono::duration_cast<chrono::milliseconds>(end - begin).count() << "ms" << endl;
}

void simulatedAnnealing(double f(vector<double>), double lowerBound, double upperBound, int dim) {
	double solMin = DBL_MAX;
	double solMax = DBL_MIN;
	double solAvg = 0;

	int timeMin = INT_MAX;
	int timeMax = INT_MIN;
	int timeAvg = 0;

	auto begin = Clock::now();
	
	for (int t = 0; t < ITERATIONS; ++t) {
		int size = ceil(log2((upperBound - lowerBound) * pow(10, PRECISION)));
		vector<bool> bitstring = generateRandomBitstring(size * dim);
		double T = ITERATIONS;

		auto start = Clock::now();

		while (T > 1) {
			vector<bool> neigh = neighbour(bitstring, randomInt(0, bitstring.size() - 1));
			
			if (f(decode(bitstring, lowerBound, upperBound, dim)) < f(decode(neigh, lowerBound, upperBound, dim))) {
				bitstring = neigh;
			}
			else if (randomDouble(0, 1) < exp(-abs(f(decode(bitstring, lowerBound, upperBound, dim)) - f(decode(neigh, lowerBound, upperBound, dim)))) / T) {
				bitstring = neigh;
			}

			T *= 0.9;
		}

		auto finish = Clock::now();

		cout << t + 1 << ". ";
		printVector(decode(bitstring, lowerBound, upperBound, dim));
		double sol = f(decode(bitstring, lowerBound, upperBound, dim));
		int time = chrono::duration_cast<chrono::milliseconds>(finish - start).count();
		cout << sol << " " << time << "ms" << endl << endl;

		if (sol < solMin) solMin = sol;
		if (sol > solMax) solMax = sol;
		solAvg += sol;

		if (time < timeMin) timeMin = time;
		if (time > timeMax) timeMax = time;
		timeAvg += time;

	}

	auto end = Clock::now();

	std::cout << "Min: " << solMin << endl
		<< "Max: " << solMax << endl
		<< "Average: " << solAvg / ITERATIONS << endl;

	std::cout << endl;

	std::cout << "Min time: " << timeMin << "ms" << endl
		<< "Max time: " << timeMax << "ms" << endl
		<< "Average: " << timeAvg / ITERATIONS << "ms" << endl;

	std::cout << endl;

	std::cout << "Total time: " << chrono::duration_cast<chrono::milliseconds>(end - begin).count() << "ms" << endl;
}

int main()
{
	int function, dimension;
	
	std::cout << "Please specify the function you want to be tested and the number of dimensions." << endl
		<< "Your options are: " << endl << endl
		<< "1 - Rastrigin's function" << endl
		<< "2 - De Jong's (Sphere) function" << endl
		<< "3 - Schwefel's function" << endl
		<< "4 - Michalewicz's function" << endl << endl;
	
	std::cout << "Your function option: "; cin >> function;
	while (function != 1 && function != 2 && function != 3 && function != 4) {
		std::cout << "That is not a valid option. Try again: "; cin >> function;
	}

	std::cout << "Your dimension option: "; cin >> dimension;
	std::cout << endl << endl;

	if (function == 1) {

		if (dimension == 10) {
			PRECISION = 2;
		}
		else if (dimension == 30) {
			PRECISION = 1;
		}

		std::cout << "=======BEST  IMPROVEMENT=======" << endl;
		hillClimbing(bestImprovement, Rastrigin, -5.12, 5.12, dimension);

		std::cout << endl << endl;

		std::cout << "=======FIRST IMPROVEMENT=======" << endl;
		hillClimbing(firstImprovement, Rastrigin, -5.12, 5.12, dimension);

		std::cout << endl << endl;

		std::cout << "======SIMULATED ANNEALING======" << endl;
		simulatedAnnealing(Rastrigin, -5.12, 5.12, dimension);
	}
	else if (function == 2) {

		if (dimension == 10) {
			PRECISION = 2;
		}
		else if (dimension == 30) {
			PRECISION = 1;
		}

		std::cout << "=======BEST  IMPROVEMENT=======" << endl;
		hillClimbing(bestImprovement, Sphere, -5.12, 5.12, dimension);

		std::cout << endl << endl;

		std::cout << "=======FIRST IMPROVEMENT=======" << endl;
		hillClimbing(firstImprovement, Sphere, -5.12, 5.12, dimension);

		std::cout << endl << endl;

		std::cout << "======SIMULATED ANNEALING======" << endl;
		simulatedAnnealing(Sphere, -5.12, 5.12, dimension);
	}
	else if (function == 3) {

		if (dimension == 10) {
			PRECISION = 2;
		}
		else if (dimension == 30) {
			PRECISION = 1;
		}

		std::cout << "=======BEST  IMPROVEMENT=======" << endl;
		hillClimbing(bestImprovement, Schwefel, -500, 500, dimension);

		std::cout << endl << endl;

		std::cout << "=======FIRST IMPROVEMENT=======" << endl;
		hillClimbing(firstImprovement, Schwefel, -500, 500, dimension);

		std::cout << endl << endl;

		std::cout << "======SIMULATED ANNEALING======" << endl;
		simulatedAnnealing(Schwefel, -500, 500, dimension);
	}
	else {

		if (dimension == 10) {
			PRECISION = 2;
		}
		else if (dimension == 30) {
			PRECISION = 1;
		}

		std::cout << "=======BEST  IMPROVEMENT=======" << endl;
		hillClimbing(bestImprovement, Michalewicz, 0, PI, dimension);

		std::cout << endl << endl;

		std::cout << "=======FIRST IMPROVEMENT=======" << endl;
		hillClimbing(firstImprovement, Michalewicz, 0, PI, dimension);

		std::cout << endl << endl;

		std::cout << "======SIMULATED ANNEALING======" << endl;
		simulatedAnnealing(Michalewicz, 0, PI, dimension);
	}

	
	return 0;
}