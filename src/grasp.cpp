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
	vector< vector< double > > traffics = instance.get_traffics();
	vector< vector< double > > distances = instance.get_distances();

	solution sol(instance, p, r);

	// Step 1 - Locating the hubs - alloc_hubs
	vector< int > hubs;
	for(int p = 0; p < sol.get_p(); p++){
		// Calculating the g(h) for each node unselected as hub
		vector< pair< double, int > > g;
		for(int h = 0; h < mesh.size(); h++){
			// Adaptive part of the procedure
			if(find(hubs.begin(), hubs.end(), mesh[h]) != hubs.end()) continue;

			// Creating the e(j,h)
			vector< double > e;
			for(int j = 0; j < instance.get_n(); j++){
				// Adaptive part of the procedure
				if(find(hubs.begin(), hubs.end(), j) != hubs.end()) continue;

				double aux = 0.0;
				for(int i = 0; i < instance.get_n(); i++)
					aux += traffics[j][i];
				aux *= distances[j][ mesh[h] ];
				e.push_back(aux);
			}

			// Sorting the elements of e(j,h)
			sort(e.begin(), e.end());

			// Chosing the k = |n/p| elements to compose g(h)
			int k = instance.get_n()/sol.get_p();
			double sum = 0.0;
			for(int j = 0; j < k; j++)
				sum += e[j];
			g.push_back(make_pair(sum, mesh[h]));
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

solution grasp::local_search_n1(solution& p_sol){
	// Evaluating the neighbors
	vector< solution > neighbors = neighborhood1(p_sol);
	solution improved = *min_element(neighbors.begin(), neighbors.end(), solution::my_sol_comparison);

	// Checking the best elements in the Neighborhood 1
	if(improved.get_total_cost() < p_sol.get_total_cost())
		return improved;
	return p_sol;
}

solution grasp::local_search_rn1( solution& p_sol ){
	vector< solution > neighbors = r_neighborhood1(p_sol);
	solution improved = *min_element(neighbors.begin(), neighbors.end(), solution::my_sol_comparison);

	// Checking the best elements in the Neighborhood 1
	if(improved.get_total_cost() < p_sol.get_total_cost())
		return improved;
	return p_sol;
}

solution grasp::local_search_n2(solution& p_sol){
	// Evaluating the neighbors
	vector< solution > neighbors = neighborhood2(p_sol);
	solution improved = *min_element(neighbors.begin(), neighbors.end(), solution::my_sol_comparison);

	// Checking the best elements in the Neighborhood 1
	if(improved.get_total_cost() < p_sol.get_total_cost())
		return improved;
	return p_sol;
}

solution grasp::local_search_rn2(solution& p_sol){
	// Evaluating the neighbors
	vector< solution > neighbors = r_neighborhood2(p_sol);
	solution improved = *min_element(neighbors.begin(), neighbors.end(), solution::my_sol_comparison);

	// Checking the best elements in the Neighborhood 1
	if(improved.get_total_cost() < p_sol.get_total_cost())
		return improved;
	return p_sol;
}

solution grasp::local_search_na( solution& p_sol ){
	vector< solution > neighbors = neighborhood_a(p_sol);
	solution improved = *min_element(neighbors.begin(), neighbors.end(), solution::my_sol_comparison);

	// Checking the best elements in the Neighborhood 1
	if(improved.get_total_cost() < p_sol.get_total_cost())
		return improved;
	return p_sol;
}

vector<solution> grasp::neighborhood1( solution& p_sol ){
	/**
	 * Generate a Neighborhood based on a one-point exchange of the worst evaluated
	 * Hub in the p_sol (partial solution) given.
	 */

	vector< double > temp = p_sol.get_hubs_cost();
	vector< double >::iterator temp_it = max_element(temp.begin(), temp.end());
	int h = temp_it - temp.begin(); // hub to be exchanged

	// Generating Neighborhood
	vector< solution > neighbors;
	for(int i = 0; i < mesh.size(); i++){
		if(p_sol.is_hub(mesh[i])) continue;
		vector< int > hubs(p_sol.get_alloc_hubs());
		hubs[h] = mesh[i];

		solution s1(instance, p, r);
		s1.set_alloc_hubs(hubs);
		s1.assign_hubs();
//		s1.set_assigned_hubs(p_sol.get_assigned_hubs());
//		s1.assign_partial_hubs(h, i); // Changing assigned_hubs in one position
		s1.route_traffics();
		neighbors.push_back(s1);
	}
	/*printf("\nneighbor #\t Total cost:\t Hubs Cost:\t Razão\n");
	for(int i = 0; i < neighbors.size(); i++){
		if(p_sol.is_hub(i)) continue;
		neighbors[i].show_data();
		printf("\n#%d: %.2lf\n", i, neighbors[i].get_total_cost());
	}*/
	return neighbors;
}

vector<solution> grasp::r_neighborhood1( solution& p_sol ){
	/**
	 * Generate a Neighborhood based on a one-point exchange of a random Hub
	 * in p_sol (partial solution) given.
	 */

	// Generating Neighborhood
	vector< solution > neighbors;
	for(int i = 0; i < mesh.size(); i++){
		if(p_sol.is_hub(mesh[i])) continue;
		int h = rand() % p; // hub to be exchanged
		vector< int > hubs(p_sol.get_alloc_hubs());
		hubs[h] = mesh[i];

		solution s1(instance, p, r);
		s1.set_alloc_hubs(hubs);
		s1.assign_hubs();
		s1.route_traffics();
		neighbors.push_back(s1);
	}

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
	for(int i = 0; i < mesh.size(); i++){
		if(p_sol.is_hub(mesh[i])) continue;
		for(int j = 0; j < mesh.size(); j++){
			if(p_sol.is_hub(mesh[j])) continue;
			vector< int > hubs(p_sol.get_alloc_hubs());
			hubs[h1] = mesh[i];
			hubs[h2] = mesh[j];
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

vector<solution> grasp::r_neighborhood2( solution& p_sol ){
	/**
	 * Generate a Neighborhood based on a two-point exchange of two random Hubs
	 * in p_sol (partial solution) given.
	 */

	// Generating Neighborhood
	vector< solution > neighbors;
	for(int h = 0; h < p; h++){
		vector< int > hubs(p_sol.get_alloc_hubs());
		hubs[h] = rand() % mesh.size();

		solution s1(instance, p, r);
		s1.set_alloc_hubs(hubs);
		s1.assign_hubs();
		s1.route_traffics();
		neighbors.push_back(s1);
	}

	return neighbors;
}

vector< solution > grasp::neighborhood_a( solution& p_sol ){
	// Generating Neighborhood
	vector< solution > neighbors;
	vector< int > hubs(p_sol.get_alloc_hubs());
	for(int i = 0; i < instance.get_n(); i++){
		vector< int > assigned_hubs(p_sol.get_assigned_hubs(i));
		vector< int > to_assign;
		for(int j = 0; j < hubs.size(); j++)
			if(find(assigned_hubs.begin(), assigned_hubs.end(), hubs[j]) == assigned_hubs.end())
				to_assign.push_back(hubs[j]);
		int x = rand() % to_assign.size();
		int h = rand() % r;

		solution s1(instance, p, r);
		s1.set_alloc_hubs(hubs);
		s1.set_assigned_hubs(p_sol.get_assigned_hubs());
		s1.set_assigned_hub(i, h, to_assign[x]);
		s1.route_traffics();
		neighbors.push_back(s1);
	}

	return neighbors;
}

void grasp::preprocessing(){
	vector< vector< double > > traffics = instance.get_traffics();
	vector< vector< double > > distances = instance.get_distances();
	int n = instance.get_n();

	vector< pair< double, int > > aux;
	for(int h = 0; h < n; h++){
		// Creating the e(j,h)
		vector< double > e;
		for(int j = 0; j < n; j++){
			if(h == j) continue;

			double aux = 0.0;
			for(int i = 0; i < n; i++)
				aux += traffics[j][i];
			aux *= distances[j][h];
			e.push_back(aux);
		}

		// Sorting the elements of e(j,h)
		sort(e.begin(), e.end());

		// Chosing the k = |n/p| elements to compose g(h)
		int k = n/p;
		double sum = 0.0;
		for(int j = 0; j < k; j++)
			sum += e[j];
		aux.push_back(make_pair(sum, h));
	}

	// Eliminating the x percent most expensive
	vector< int > mesh;
	int x = (0.1 * n);
	sort(aux.begin(), aux.end(), my_comparison);
	for(int i = 0; i < aux.size() - x; i++)
		mesh.push_back(aux[i].second);

	this->mesh = mesh;
}

solution& grasp::execute(){
	// Pre-processing
	preprocessing();

	// processing
	vector< double > partial_improment;
	int t_max = 0.1 * max_iterations;
	double mean = 0.0;
	for(int i = 0; i < max_iterations; i++){
		solution initial = greedy_randomized_construction();
//		printf("\niteration #%d:\n", i);
//		printf("Initial Solution:\n");
//		initial.show_data();

		// Filtering mechanism
		if(i == t_max - 1){
			for(int j = 0; j < t_max; j++)
				mean += partial_improment[j];
			mean /= instance.get_n();
		}
		if(i >= t_max){
			double delta_pi = (initial.get_total_cost() - best.get_total_cost()) / initial.get_total_cost();
			if(delta_pi >= mean) continue;
		}

		solution improved = local_search_rn1(initial);
		improved = local_search_na(improved);
//		if(improved.get_total_cost() == initial.get_total_cost())
//			improved = local_search_rn2(initial);
//		printf("Improved Solution:\n");
//		improved.show_data();

		// Partial improvement for the filtering mechanism
		if(i < t_max){
			double _aux = (initial.get_total_cost() - improved.get_total_cost()) / initial.get_total_cost();
			partial_improment.push_back(_aux);
		}

		// Acceptance criteria
		if(i != 0){
			if(improved.get_total_cost() < best.get_total_cost())
				set_best(improved);
		}else best = improved;

//		printf("Best Solution:\n"); best.show_data();
	}

	// Post-processing
//	set_best(local_search_n2(best));

	return best;
}
