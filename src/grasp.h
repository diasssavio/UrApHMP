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
#include "FWChrono.h"

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
	double alpha; // alpha in [0,1]

	// Input instance
	uraphmp instance;

	// Current Solution
	int p; // Number of hubs to be allocated
	int r; // Number of hubs to be assigned to each node
	solution best;

	// Preprocessed mesh of probable hubs
	vector< int > mesh;

	// Iterations
	vector< pair< double, int > > it_log;
	vector< double > times;
	vector< int > path;
	FWChrono timer;

	// Private methods
	static bool my_comparison( pair< double, int >, pair< double, int > );

public:
	// Constructors & Destructors
	grasp( uraphmp&, size_t, int, int, double, FWChrono& );
	virtual ~grasp();

	// Setters
	void set_max_iterations( size_t );
	void set_instance( uraphmp& );
	void set_best( const solution& );

	// Getters
	size_t get_max_iterations();
	uraphmp& get_instance();
	solution& get_best();
	vector< pair< double, int > >& get_it_log();
	vector< double >& get_times();
	vector< int >& get_path();

	// Useful Methods
	void preprocessing();
	solution greedy_randomized_construction();
	solution local_search_n1( solution& );
	solution local_search_rn1( solution& );
	solution local_search_c2n1( solution& );
	solution local_search_na( solution& );
//	solution local_search_n2( solution& );
//	solution local_search_rn2( solution& );

	vector< solution > neighborhood1( solution& );
	vector< solution > r_neighborhood1( solution& );
	vector< solution > closest2_n1( solution& );
	vector< solution > neighborhood_a( solution& );
//	vector< solution > neighborhood2( solution& );
//	vector< solution > r_neighborhood2( solution& );

	solution path_relinking( solution&, solution& );

	// TODO 5.Use a Hash Table to avoid calculus of same solutions on the neighbors **
	// 		Check whether the same hubs were allocated in a solution so to avoid the calculation
	solution& execute();
};

#endif /* GRASP_H_ */
