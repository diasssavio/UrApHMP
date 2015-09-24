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
#include "FWChrono.h"

using namespace std;



int main() {
	FWChrono timer;
	timer.start();
//	double timer_1 = cpuTime();

	srand(time(NULL));

	int n;

	scanf("%d", &n);
	double X = 1.0, alpha_1 = 0.2, delta = 1.0;
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

	int max_iterations = 0.1 * n;
	int p = 5;
	int r = 3;
	double alpha_2 = 0.2;
	grasp GRASP(instance, max_iterations, p, r, alpha_2, timer);
	solution result = GRASP.execute();
//	double timer_2 = cpuTime();
	timer.stop();
	printf("Tempo de execução: %.2lf", timer.getStopTime());
//	printf("Tempo de execução: %.2lf", (timer_2-timer_1));
	result.show_data();

	printf("\nIT_LOG:\n");
	vector< pair< double, int> > it_log = GRASP.get_it_log();
	vector< double > times = GRASP.get_times();
	vector< int > path = GRASP.get_path();
	double min_time = 0.0;
	for(int i = 0; i < it_log.size(); i++){
		printf("#%d:\t%.2lf\t%.2lf\n", it_log[i].second, it_log[i].first, times[i]);

		if(it_log[i].first == result.get_total_cost() && min_time == 0.0)
			min_time = times[i];
	}

	printf("\nPATH RELINKING USAGE:\n");
	for(int i = 0; i < path.size(); i++)
		printf("%d\t", path[i]);

	printf("\nMIN_TIME: %.2lf\n", min_time);

	return 0;
}
