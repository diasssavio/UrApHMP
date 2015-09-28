/*
 * ils.cpp
 *
 *  Created on: Apr 20, 2015
 *      Author: Sávio S. Dias
 */

#include "../include/ils.h"

ils::ils( uraphmp& instance, size_t max_iterations, int p, int r, FWChrono& timer ) : max_iterations(max_iterations), p(p), r(r) {
	this->set_instance(instance);
	this->timer = timer;
}

ils::~ils() { }

void ils::set_c2n1(const vector<solution>& c2n1) { this->c2n1 = c2n1; }

void ils::set_na(const vector<solution>& na) { this->na = na; }

void ils::set_rn1(const vector<solution>& rn1) { this->rn1 = rn1; }

const solution& ils::get_best() const { return best; }

const uraphmp& ils::get_instance() const { return instance; }

size_t ils::get_max_iterations() const { return max_iterations; }

const vector<pair<double, int> >& ils::get_it_log() const { return it_log; }

const vector<double>& ils::get_times() const { return times; }

void ils::set_best(const solution& best) { this->best = best; }

void ils::set_instance(const uraphmp& instance) { this->instance = instance; }

void ils::set_max_iterations(size_t maxIterations) { max_iterations = maxIterations; }

bool ils::my_comparison( pair< double, int > p1, pair< double, int > p2 ){
	return (p1.first < p2.first);
}

void ils::preprocessing(){
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
	for(unsigned i = 0; i < aux.size() - x; i++)
		mesh.push_back(aux[i].second);

	this->mesh = mesh;
}

solution ils::constructor(){
	vector< vector< double > > traffics = instance.get_traffics();
	vector< vector< double > > distances = instance.get_distances();

	solution sol(instance, p, r);

	// Step 1 - Locating the hubs - alloc_hubs
	vector< int > hubs;
	for(int p = 0; p < sol.get_p(); p++){
		// Calculating the g(h) for each node unselected as hub - The Candidate List
		vector< pair< double, int > > g;
		for(unsigned h = 0; h < mesh.size(); h++){
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

			// Choosing the k = |n/p| elements to compose g(h)
			int k = instance.get_n()/sol.get_p();
			double sum = 0.0;
			for(int j = 0; j < k; j++)
				sum += e[j];
			g.push_back(make_pair(sum, mesh[h]));
		}

		// Creating the RCL based on g(h) -- GRASP
		/*double g_min = (*min_element(g.begin(), g.end(), my_comparison)).first;
		double g_max = (*max_element(g.begin(), g.end(), my_comparison)).first;
		vector< int > RCL;
		for(unsigned i = 0; i < g.size(); i++)
			if(g[i].first <= (g_min + this->alpha*(g_max - g_min)))
				RCL.push_back(g[i].second);

		// Selecting randomly the hub from RCL
		hubs.push_back(RCL[ rand() % RCL.size() ]);*/

		pair< double, int > g_min = *min_element(g.begin(), g.end(), my_comparison);
		hubs.push_back(g_min.second);
	}
	sol.set_alloc_hubs(hubs);

	// Step 2 - Assignment of the r hubs - assigned_hubs
	sol.assign_hubs();

	// Step 3 - Route traffics
	sol.route_traffics();

	return sol;
}

solution ils::local_search_rn1( solution& p_sol ){
	// Checking the best elements in the Neighborhood 1
	solution partial = p_sol;
	solution improved;
	bool is_improved = true;
	while(is_improved){
		improved = r_neighborhood1(partial);
//		improved = neighbors[ neighbors.size() - 1 ];
		if(improved.get_total_cost() < partial.get_total_cost())
			partial = improved;
		else is_improved = false;
	}

	return partial;
}

solution ils::local_search_c2n1( solution& p_sol ){
	// TODO 1.2.Test the use w/ RN1

	// Checking the best elements in the Adjacent points
	solution partial = p_sol;
	solution improved;
	bool is_improved = true;
	while(is_improved){
		improved = closest2_n1(partial);
//		improved = neighbors[ neighbors.size() - 1 ];
		if(improved.get_total_cost() < partial.get_total_cost())
			partial = improved;
		else is_improved = false;
	}

	return partial;
}

solution ils::local_search_na( solution& p_sol ){
	// Checking the best solution in the Neighborhood A

	solution partial = p_sol;
	solution improved;
	bool is_improved = true;
	while(is_improved){
		improved = neighborhood_a(partial);
//		improved = neighbors[ neighbors.size() - 1 ];
		if(improved.get_total_cost() < partial.get_total_cost())
			partial = improved;
		else is_improved = false;
	}

	return improved;
}

solution& ils::r_neighborhood1( solution& p_sol ){
	/**
	 * Generate a Neighborhood based on a one-point exchange of a random Hub
	 * in p_sol (partial solution) given.
	 */

	// Generating Neighborhood
	vector< solution > neighbors;
	for(unsigned i = 0; i < mesh.size(); i++){
		if(p_sol.is_hub(mesh[i])) continue;
		int h = rand() % p; // hub to be exchanged
		vector< int > hubs(p_sol.get_alloc_hubs());
		hubs[h] = mesh[i];

		solution s1(instance, p, r);
		s1.set_alloc_hubs(hubs);
		s1.assign_hubs();
		s1.route_traffics();
		neighbors.push_back(s1);
//		if(s1.get_cost() < p_sol.get_cost()) break;
	}
	set_rn1(neighbors);

	return *min_element(rn1.begin(), rn1.end(), solution::my_sol_comparison);
}

solution& ils::closest2_n1( solution& p_sol ){
	/**
	 * Generate a Neighborhood based on a one-point exchange of the two closest points of Hubs
	 * in p_sol (partial solution) given.
	 */
	// TODO 1.1.Check the generation of previously calculated solutions as neighbors
	//		There's many intersection of neighbors

	vector< vector< double > > distances = instance.get_distances();

	// Generating Neighborhood
	vector< solution > neighbors; // Analyzing the 2 closest indexes points for each hub
	vector< int > alloc_hubs = p_sol.get_alloc_hubs();
	for(int i = 0; i < p; i++){
		// Finding the first two non-hub nodes
		unsigned j = 0;
		while(p_sol.is_hub(mesh[j]))
			j++;
		unsigned k = j + 1;
		while(p_sol.is_hub(mesh[k]))
			k++;
		pair< unsigned, double > min1, min2;
		if(distances[alloc_hubs[i]][mesh[k]] < distances[alloc_hubs[i]][mesh[j]]){
			min1 = make_pair(mesh[k], distances[alloc_hubs[i]][mesh[k]]);
			min2 = make_pair(mesh[j], distances[alloc_hubs[i]][mesh[j]]);
		} else {
			min2 = make_pair(mesh[k], distances[alloc_hubs[i]][mesh[k]]);
			min1 = make_pair(mesh[j], distances[alloc_hubs[i]][mesh[j]]);
		}
		// Getting the two Euclidian closest non-hub nodes
		for(unsigned j = k + 1; j < mesh.size(); j++){
			if(alloc_hubs[i] == mesh[j] || p_sol.is_hub(mesh[j])) continue; // main diagonal & hub nodes avoidance
			if(distances[alloc_hubs[i]][mesh[j]] < min1.second){
				min2 = min1;
				min1 = make_pair(mesh[j], distances[alloc_hubs[i]][mesh[j]]);
			}
		}

		// Making the solution objects
		vector< int > hubs1(alloc_hubs);
		vector< int > hubs2(alloc_hubs);
		hubs1[i] = min1.first;
		hubs2[i] = min2.first;

		solution s1(instance, p, r);
		solution s2(instance, p, r);
		s1.set_alloc_hubs(hubs1);
		s1.assign_hubs();
		s1.route_traffics();
		s2.set_alloc_hubs(hubs2);
		s2.assign_hubs();
		s2.route_traffics();

		neighbors.push_back(s1);
		neighbors.push_back(s2);
	}

	set_c2n1(neighbors);

	return *min_element(c2n1.begin(), c2n1.end(), solution::my_sol_comparison);
}

solution& ils::neighborhood_a( solution& p_sol ){
	// Generating Neighborhood A

	vector< solution > neighbors;
	vector< int > hubs(p_sol.get_alloc_hubs());
	for(int i = 0; i < instance.get_n(); i++){
		vector< int > assigned_hubs(p_sol.get_assigned_hubs(i));

		// Finding the symmetric difference
		vector< int > to_assign;
		for(unsigned j = 0; j < hubs.size(); j++)
			if(find(assigned_hubs.begin(), assigned_hubs.end(), hubs[j]) == assigned_hubs.end())
				to_assign.push_back(hubs[j]);

		// Generating the neighbors
		for(int j = 0; j < this->r; j++)
			for(unsigned k = 0; k < to_assign.size(); k++){
				solution s1(instance, p, r);
				s1.set_alloc_hubs(hubs);
				s1.set_assigned_hubs(p_sol.get_assigned_hubs());
				s1.set_assigned_hub(i, j, to_assign[k]);
				s1.route_traffics();
				neighbors.push_back(s1);
//				if(s1.get_cost() < p_sol.get_cost())
//					return neighbors;
			}
	}

	set_na(neighbors);

	return *min_element(na.begin(), na.end(), solution::my_sol_comparison);
}

solution& ils::execute(){
	// Pre-processing
	preprocessing();

	// Processing

	// Constructing initial solution
	solution initial = constructor();
	solution improved = initial;

	unsigned i = 1, k = 0;
	bool first = true;
	while(i < max_iterations){
		// Local Search
//		improved = local_search_c2n1(improved);
		improved = local_search_rn1(improved);

		// Acceptance criterion & VND
		if(!first){
			if(improved.get_total_cost() < best.get_total_cost()) {
				set_best(improved);
				i = 1;
			} else { // Testing the LS_a only when LS_h doesn't improve the solution
				improved = local_search_na(improved);
				if(improved.get_total_cost() < best.get_total_cost()){
					set_best(improved);
					i = 1;
				}else i++;
			}
		} else {
			best = improved;
			first = false;
		}

		// Shaking phase
		improved = rn1[rand() % rn1.size()];

		// Saving the execution logs
		it_log.push_back(make_pair(best.get_total_cost(), k++));
		times.push_back(((double) timer.getMilliSpan() / 1000));
	}


	// Post-processing
//	set_best(local_search_rn1(best));

	return best;
}
