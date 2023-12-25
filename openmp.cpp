#include "pch.h"
#include <iostream>
#include <omp.h>
#include <stdio.h>
#include <time.h>     
#include <vector>
#include <stdlib.h>
#include <limits>
#include <cmath>
#include <windows.h>     //раскоментировать только для последней задачи

// Task 1
typedef unsigned int Element;
typedef std::vector<Element> Array;

unsigned int getMinElement(const Array& array, int threadsCount) {
	omp_set_num_threads(threadsCount);
	Array minElements(
		threadsCount,
		std::numeric_limits<Element>::max()
	);
	int size(array.size());
#pragma omp parallel for schedule(dynamic, 1000)
	for (int i = 0; i < size; ++i) {
		int threadNum(omp_get_thread_num());
		if (array[i] <= minElements[threadNum])
			minElements[threadNum] = array[i];
	}
	Element minElement(minElements[0]);
	for (int i = 1; i < minElements.size(); ++i)
		if (minElements[i] < minElement)
			minElement = minElements[i];
	return minElement;
}

int main() {
	Array a;
	for (unsigned int i = 0; i < 1000000; ++i)
		a.push_back(rand());
	int n;
	for (int i = 0; i < 5; i++) {
		clock_t startTime(clock());
		unsigned int minElement(getMinElement(a, pow(2, i)));
		clock_t arrayProcessTime(clock() - startTime);
		std::cout << "Threads: "
			<< pow(2, i)
			<< "\nTime: "
			<< double(arrayProcessTime) / CLOCKS_PER_SEC 
			<< std::endl;
	}
	return 0;
}


//Task 2
int main() { 
	int i, n, chunk;
	float a[100000], b[100000], result;
	n = 100000; chunk = 1000;
	
	for (i = 0; i < n; i++) {
		a[i] = rand(); b[i] = rand();
	}
	for (int k = 0; k < 5; k++) {
		result = 0.0;
		omp_set_num_threads(pow(2, k));
		clock_t startTime(clock());
#pragma omp parallel for schedule(dynamic, 1000) reduction(+:result)
		for (i = 0; i < n; i++)
			result = result + (a[i] * b[i]);
		clock_t arrayProcessTime(clock() - startTime);
		std::cout << "Threads: "
			<< pow(2, k)
			<< " Time: "
			<< double(arrayProcessTime) / CLOCKS_PER_SEC
			<< std::endl;
	}
	return 0;
}


//Task 3
double func(double x)
{
	return exp(-x * x);
}
int main(int argc, char **argv)
{
	const double a = -4.0;
	const double b = 4.0;
	const int n = 1000000;
	for (int i = 0; i < 5; i++) {
		clock_t starttime(clock());
		omp_set_num_threads(pow(2, i));
		double h = (b - a) / n;
		double s = 0.0;

#pragma omp parallel for reduction(+:s)
		for (int i = 0; i < n; i++)
			s += func(a + h * i);
		s *= h;
		clock_t arrayprocesstime(clock() - starttime);
		std::cout << "threads: "
			<< pow(2, i)
			<< " time: "
			<< double(arrayprocesstime) / CLOCKS_PER_SEC
			<< std::endl;
	}
	return 0;
}

//Task 4
#define N 500
int main() { 
	int i, j, n, m;
	int min, max;
	int a[N][N], res[N];
	n = N; 
	m = N;
	srand(100);
	for (i = 0; i < n; i++) {
		res[i] = -1;
		for (j = 0; j < m; j++)
			a[i][j] = rand();
	}
	for (int k = 0; k < 5; k++) {
		clock_t startTime(clock());
		omp_set_num_threads(pow(2, k));
#pragma omp parallel for shared(res) schedule(dynamic, 1)
		for (i = 0; i < n; i++)
		{

			res[i] = a[i][0];
#pragma omp parallel for
			for (j = 0; j < m; j++)
				if (a[i][j] < res[i])
					res[i] = a[i][j];
		}
		max = res[0];
		for (i = 0; i < n; i++)
			if (res[i] > max)
				max = res[i];
		clock_t arrayProcessTime(clock() - startTime);
		std::cout << "threads: "
			<< pow(2, k)
			<< " time: "
			<< double(arrayProcessTime) / CLOCKS_PER_SEC
			<< std::endl;
	}
	return 0;
}

//Task 11
int main() { 
	int i, n, chunk;
	float a[100000];
	n = 100000; chunk = 5;

	std::cout << "static\n";
	for (i = 0; i < n; i++) {
		a[i] = 0 + rand() % (1 - 0);;
	}
	for (int k = 0; k < 5; k++) {
		omp_set_num_threads(pow(2, k));
		clock_t startTime(clock());
#pragma omp parallel for  schedule(static) 
		for (i = 0; i < n; i++)
			Sleep(a[i]);
		clock_t arrayProcessTime(clock() - startTime);
		std::cout << "Threads: "
			<< pow(2, k)
			<< " Time: "
			<< double(arrayProcessTime) / CLOCKS_PER_SEC
			<< std::endl;
	}
	std::cout << "dynamic\n";
	for (int k = 0; k < 5; k++) {
		omp_set_num_threads(pow(2, k));
		clock_t startTime(clock());
#pragma omp parallel for  schedule(dynamic) 
		for (i = 0; i < n; i++)
			Sleep(a[i]);
		clock_t arrayProcessTime(clock() - startTime);
		std::cout << "Threads: "
			<< pow(2, k)
			<< " Time: "
			<< double(arrayProcessTime) / CLOCKS_PER_SEC
			<< std::endl;
	}
	std::cout << "dynamic with chunk \n";
	for (int k = 0; k < 5; k++) {
		omp_set_num_threads(pow(2, k));
		clock_t startTime(clock());
#pragma omp parallel for  schedule(dynamic, 1000) 
		for (i = 0; i < n; i++)
			Sleep(a[i]);
		clock_t arrayProcessTime(clock() - startTime);
		std::cout << "Threads: "
			<< pow(2, k)
			<< " Time: "
			<< double(arrayProcessTime) / CLOCKS_PER_SEC
			<< std::endl;
	}
	std::cout << "guided\n";
	for (int k = 0; k < 5; k++) {
		omp_set_num_threads(pow(2, k));
		clock_t startTime(clock());
#pragma omp parallel for  schedule(guided) 
		for (i = 0; i < n; i++)
			Sleep(a[i]);
		clock_t arrayProcessTime(clock() - startTime);
		std::cout << "Threads: "
			<< pow(2, k)
			<< " Time: "
			<< double(arrayProcessTime) / CLOCKS_PER_SEC
			<< std::endl;
	}
	return 0;
}

