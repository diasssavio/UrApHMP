/*
 * grasp.h
 *
 *  Created on: Apr 30, 2015
 *      Author: Sávio S. Dias
 */

#ifndef GRASP_H_
#define GRASP_H_

#include "UrApHMP.h"
#include "solution.h"

#include <utility>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <limits>

using namespace std;

class grasp {
private:
	// GRASP Parameters
	size_t max_iterations; // Max number of iterations
	double alpha; // alpha [0,1]
//	time_t seed; // Seed for random selection

	// Input instance
	uraphmp instance;

	// Current Solution
	int p; // Number of hubs to be allocated
	int r; // Number of hubs to be assigned to each node
	solution best;
	vector< solution > neighbors;

	// Private methods
	static bool my_comparison( pair< double, int >, pair< double, int > );

public:
	// Constructors & Destructors
	grasp( uraphmp&, size_t, int, int, double );
	virtual ~grasp();

	// Setters
	void set_max_iterations( size_t );
	void set_instance( uraphmp& );
	void set_best( const solution& );

	// Getters
	size_t get_max_iterations();
	uraphmp& get_instance();
	solution& get_best();

	// Useful Methods
	solution greedy_randomized_construction();
	solution local_search_n1( solution& ); // TODO Implement the local search testing new neighborhoods
	solution local_search_n1_min( solution& );
	solution local_search_n2( solution& );

	vector<solution> neighborhood1( solution& );
	vector<solution> neighborhood1_min( solution& );
	vector<solution> neighborhood2( solution& );
	// TODO Create a exchange route neighborhood Na

	solution& execute();
};

#endif /* GRASP_H_ */
