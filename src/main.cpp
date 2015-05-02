//============================================================================
// Name        : main.cpp
// Author      : Sávio S. Dias
// Version     : 1.0
// Copyright   :
// Description : A GRASP-based heuristic for the UrApHMP
//============================================================================

#include <iostream>
#include <vector>
#include <cstdio>
#include <ctime>

#include "UrApHMP.h"
#include "solution.h"
#include "grasp.h"

using namespace std;



int main() {
	int n;

	scanf("%d", &n);
	double X = 3.0, alpha_1 = 0.75, delta = 2.0;
	uraphmp instance(n, X, alpha_1, delta);

	vector< vector< double> > aux;
	for(int i = 0; i < n; i++){
		vector< double > aux2;
		for(int j = 0; j < n; j++){
			double temp;
			scanf("%lf", &temp);
			aux2.push_back(temp);
		}
		aux.push_back(aux2);
	}
	instance.set_traffics(aux);

	aux.clear();
	for(int i = 0; i < n; i++){
		vector< double > aux2;
		for(int j = 0; j < n; j++){
			double temp;
			scanf("%lf", &temp);
			aux2.push_back(temp);
		}
		aux.push_back(aux2);
	}
	instance.set_distances(aux);
	aux.clear();

	int max_iterations = 10;
	int p = 4;
	int r = 2;
	double alpha_2 = 0.8;
	srand(time(NULL));
	grasp GRASP(instance, max_iterations, p, r, alpha_2);
	GRASP.execute();
//	GRASP.set_best();
//	GRASP.get_best().show_data();
//	GRASP.greedy_randomized_construction().show_data();

	return 0;
}