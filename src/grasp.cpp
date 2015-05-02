/*
 * grasp.cpp
 *
 *  Created on: Apr 30, 2015
 *      Author: Sávio S. Dias
 */

#include "grasp.h"

grasp::grasp( uraphmp& instance, size_t max_iterations, int p, int r, double alpha ) : max_iterations(max_iterations), p(p), r(r), alpha(alpha) {
	this->set_instance(instance);
}

grasp::~grasp() {

}

void grasp::set_max_iterations( size_t max_iterations ){
	this->max_iterations = max_iterations;
}

void grasp::set_instance( uraphmp& instance ){
	this->instance = instance;
}

void grasp::set_best( const solution& sol ){
	this->best = sol;
}

size_t grasp::get_max_iterations(){
	return this->max_iterations;
}

uraphmp& grasp::get_instance(){
	return this->instance;
}

solution& grasp::get_best(){
	return this->best;
}

bool grasp::my_comparison( pair< double, int > p1, pair< double, int > p2 ){
	return (p1.first < p2.first);
}

solution grasp::greedy_randomized_construction(){
//	srand(this->seed);

	vector< vector< double > > traffics = instance.get_traffics();
	vector< vector< double > > distances = instance.get_distances();

	solution sol(instance, p, r);

	// Step 1 - Locating the hubs - alloc_hubs
	vector< int > hubs;
	for(int p = 0; p < sol.get_p(); p++){
		// Calculating the g(h) for each node unselected as hub
		vector< pair< double, int > > g;
		for(int h = 0; h < instance.get_n(); h++){
			// Adaptive part of the procedure
			if(find(hubs.begin(), hubs.end(), h) != hubs.end()) continue;

			// Creating the e(j,h)
			vector< double > e;
			for(int j = 0; j < instance.get_n(); j++){
				// Adaptive part of the procedure
				if(find(hubs.begin(), hubs.end(), j) != hubs.end()) continue;

				double aux = 0.0;
				for(int i = 0; i < instance.get_n(); i++)
					aux += traffics[j][i];
				aux *= distances[j][h];
				e.push_back(aux);
			}

			// Sorting the elements of e(j,h)
			sort(e.begin(), e.end());

			// Chosing the r elements to compose g(h)
			double sum = 0.0;
			for(int j = 0; j < sol.get_r(); j++)
				sum += e[j];
			g.push_back(make_pair(sum, h));
		}

		// Creating the RCL based on g(h)
		double g_min = (*min_element(g.begin(), g.end(), my_comparison)).first;
		double g_max = (*max_element(g.begin(), g.end(), my_comparison)).first;
		vector< int > RCL;
		for(int i = 0; i < g.size(); i++)
			if(g[i].first <= (g_min + this->alpha*(g_max - g_min)))
				RCL.push_back(g[i].second);

		// Selecting randomly the hub from RCL
		hubs.push_back(RCL[ rand() % RCL.size() ]);
	}
	sol.set_alloc_hubs(hubs);

	// Step 2 - Assignment of the r hubs - assigned_hubs
	sol.assign_hubs();

	// Step 3 - Route traffics
	sol.route_traffics();

	return sol;
}

solution grasp::local_search(solution& p_sol){
	// TODO Implement a VND local search with a 4 neighborhoods structure

	// Evaluating the neighbors
	vector< solution > neighbors = neighborhood2(p_sol);
	solution improved = *min_element(neighbors.begin(), neighbors.end(), solution::my_sol_comparison);

	// Checking the best elements in the Neighborhood 1 & 2
	if(improved.get_total_cost() < p_sol.get_total_cost())
		return improved;
	/*else{
		neighbors = neighborhood2(p_sol);
		improved = *min_element(neighbors.begin(), neighbors.end(), solution::my_sol_comparison);
		if(improved.get_total_cost() < p_sol.get_total_cost())
			return improved;
	}*/
	return p_sol;
}

vector<solution> grasp::neighborhood1( solution& p_sol ){
	/**
	 * Generate a Neighborhood based on a one-point exchange of the worst evaluated
	 * Hub in the p_sol (partial solution) given.
	 */

	vector< double > temp = p_sol.get_hubs_cost();
	vector< double >::iterator temp_it = max_element(temp.begin(), temp.end());
	int h = temp_it - temp.begin();

	// Generating Neighborhood
	vector< solution > neighbors;
	for(int i = 0; i < instance.get_n(); i++){
		if(p_sol.is_hub(i)) continue;
		vector< int > hubs(p_sol.get_alloc_hubs());
		hubs[h] = i;
		solution s1(instance, p, r);
		s1.set_alloc_hubs(hubs);
		s1.assign_hubs();
		s1.route_traffics();
		neighbors.push_back(s1);
	}
	/*for(int i = 0; i < neighbors.size(); i++){
		if(p_sol.is_hub(i)) continue;
		neighbors[i].show_data();
		printf("\nneighbor #%d: %.2lf\n", i, neighbors[i].get_total_cost());
	}*/
	return neighbors;
}

vector<solution> grasp::neighborhood2( solution& p_sol ){
	/**
	 * Generate a Neighborhood based on a two-point exchange of the two-worst evaluated
	 * Hubs in the p_sol (partial solution) given.
	 */

	// Finding the two expensive hubs in the solution
	vector< double > temp(p_sol.get_hubs_cost());
	vector< double >::iterator it1 = max_element(temp.begin(), temp.end());
	int h1 = it1 - temp.begin();
	temp.erase(it1);
	it1 = max_element(temp.begin(), temp.end());
	int h2 = it1 - temp.begin();

	// Generating the neighborhood
	vector<solution> neighbors;
	for(int i = 0; i < instance.get_n(); i++){
		if(p_sol.is_hub(i)) continue;
		for(int j = 0; j < instance.get_n(); j++){
			if(p_sol.is_hub(j)) continue;
			vector< int > hubs(p_sol.get_alloc_hubs());
			hubs[h1] = i;
			hubs[h2] = j;
			solution s1(instance, p, r);
			s1.set_alloc_hubs(hubs);
			neighbors.push_back(s1);
		}
	}
	for(int i = 0; i < neighbors.size(); i++){
		neighbors[i].assign_hubs();
		neighbors[i].route_traffics();
	}
	/*for(int i = 0; i < neighbors.size(); i++){
//		if(p_sol.is_hub(i)) continue;
		neighbors[i].show_data();
//		printf("\nneighbor #%d: %.2lf\n", i, neighbors[i].get_total_cost());
	}*/

	return neighbors;
}

solution& grasp::execute(){
	for(size_t i = 0; i < max_iterations; i++){
		solution initial = greedy_randomized_construction();
		printf("\niteration #%d:\n", i);
//		printf("Initial Solution:\n");
//		initial.show_data();

		solution improved = local_search(initial);
//		printf("Improved Solution:\n");
//		improved.show_data();

		if(i != 0){
			if(improved.get_total_cost() < best.get_total_cost())
				set_best(improved);
		}else best = improved;

		printf("Best Solution:\n"); best.show_data();
	}
	return best;
}
