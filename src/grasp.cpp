/*
 * grasp.cpp
 *
 *  Created on: Apr 30, 2015
 *      Author: Sávio S. Dias
 */

#include "grasp.h"

grasp::grasp( uraphmp& instance, size_t max_iterations, int p, int r, double alpha, FWChrono& timer ) : max_iterations(max_iterations), p(p), r(r), alpha(alpha) {
	this->set_instance(instance);
	this->timer = timer;
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

vector< pair< double, int > >& grasp::get_it_log(){
	return this->it_log;
}

vector< double >& grasp::get_times(){
	return this->times;
}

vector< int >& grasp::get_path(){
	return this->path;
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
		for(unsigned i = 0; i < g.size(); i++)
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
	// Checking the best elements in the Neighborhood 1
	solution partial = p_sol;
	solution improved;
	vector< solution > neighbors;
	do{
		neighbors = r_neighborhood1(partial);
		improved = *min_element(neighbors.begin(), neighbors.end(), solution::my_sol_comparison);
		//		improved = neighbors[ neighbors.size() - 1 ];
		if(improved.get_total_cost() < partial.get_total_cost())
			partial = improved;
	}while(improved.get_cost() == partial.get_cost());

	return partial;
}

solution grasp::local_search_c2n1( solution& p_sol ){
	// TODO 2.Test the use w/ RN1

	// Checking the best elements in the Adjacent points
	solution partial = p_sol;
	solution improved;
	vector< solution > neighbors;
	do{
		neighbors = closest2_n1(partial);
		improved = *min_element(neighbors.begin(), neighbors.end(), solution::my_sol_comparison);
//		improved = neighbors[ neighbors.size() - 1 ];
		if(improved.get_total_cost() < partial.get_total_cost())
			partial = improved;
	}while(improved.get_cost() == partial.get_cost());

	return partial;
}

solution grasp::local_search_na( solution& p_sol ){
	// Checking the best elements in the Neighborhood 1
	solution partial = p_sol;
	solution improved;
	vector< solution > neighbors;
	do{
		neighbors = neighborhood_a(partial);
		//		improved = *min_element(neighbors.begin(), neighbors.end(), solution::my_sol_comparison);
		improved = neighbors[ neighbors.size() - 1 ];
		if(improved.get_total_cost() < partial.get_total_cost())
			partial = improved;
	}while(improved.get_cost() == partial.get_cost());

	return partial;
}

//solution grasp::local_search_n2(solution& p_sol){
//	// Evaluating the neighbors
//	vector< solution > neighbors = neighborhood2(p_sol);
//	solution improved = *min_element(neighbors.begin(), neighbors.end(), solution::my_sol_comparison);
//
//	// Checking the best elements in the Neighborhood 1
//	if(improved.get_total_cost() < p_sol.get_total_cost())
//		return improved;
//	return p_sol;
//}
//
//solution grasp::local_search_rn2(solution& p_sol){
//	// Evaluating the neighbors
//	vector< solution > neighbors = r_neighborhood2(p_sol);
//	solution improved = *min_element(neighbors.begin(), neighbors.end(), solution::my_sol_comparison);
//
//	// Checking the best elements in the Neighborhood 1
//	if(improved.get_total_cost() < p_sol.get_total_cost())
//		return improved;
//	return p_sol;
//}

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
	for(unsigned i = 0; i < mesh.size(); i++){
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

	return neighbors;
}

vector< solution > grasp::closest2_n1( solution& p_sol ){
	/**
	 * Generate a Neighborhood based on a one-point exchange of the two closest points of Hubs
	 * in p_sol (partial solution) given.
	 */
	// TODO 1.Check the generation of previously calculated solutions as neighbors
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

//		if(alloc_hubs[i] - 1 > 0 && !p_sol.is_hub(alloc_hubs[i] - 1)){
//			vector< int > hubs(alloc_hubs);
//			hubs[i] = alloc_hubs[i] - 1;
//
//			solution s1(instance, p, r);
//			s1.set_alloc_hubs(hubs);
//			s1.assign_hubs();
//			s1.route_traffics();
//			neighbors.push_back(s1);
//		}
//		if((alloc_hubs[i] + 1) < (instance.get_n() - 1) && !p_sol.is_hub(alloc_hubs[i] + 1)){
//			vector< int > hubs(alloc_hubs);
//			hubs[i] = alloc_hubs[i] + 1;
//
//			solution s1(instance, p, r);
//			s1.set_alloc_hubs(hubs);
//			s1.assign_hubs();
//			s1.route_traffics();
//			neighbors.push_back(s1);
//		}
	}

	return neighbors;
}

vector< solution > grasp::neighborhood_a( solution& p_sol ){
	// Generating Neighborhood
	// TODO 3.Implement a full N-A, where all exchanges in positions are made, not only at random way
	// TODO 4.Implement the partial evaluation of the solution for this neighborhood
	vector< solution > neighbors;
	vector< int > hubs(p_sol.get_alloc_hubs());
	for(int i = 0; i < instance.get_n(); i++){
		vector< int > assigned_hubs(p_sol.get_assigned_hubs(i));
		vector< int > to_assign;
		for(unsigned j = 0; j < hubs.size(); j++)
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
		if(s1.get_cost() < p_sol.get_cost()) break;
	}

	return neighbors;
}

//vector<solution> grasp::neighborhood2( solution& p_sol ){
//	/**
//	 * Generate a Neighborhood based on a two-point exchange of the two-worst evaluated
//	 * Hubs in the p_sol (partial solution) given.
//	 */
//
//	// Finding the two expensive hubs in the solution
//	vector< double > temp(p_sol.get_hubs_cost());
//	vector< double >::iterator it1 = max_element(temp.begin(), temp.end());
//	int h1 = it1 - temp.begin();
//	temp.erase(it1);
//	it1 = max_element(temp.begin(), temp.end());
//	int h2 = it1 - temp.begin();
//
//	// Generating the neighborhood
//	vector<solution> neighbors;
//	for(unsigned i = 0; i < mesh.size(); i++){
//		if(p_sol.is_hub(mesh[i])) continue;
//		for(unsigned j = 0; j < mesh.size(); j++){
//			if(p_sol.is_hub(mesh[j])) continue;
//			vector< int > hubs(p_sol.get_alloc_hubs());
//			hubs[h1] = mesh[i];
//			hubs[h2] = mesh[j];
//			solution s1(instance, p, r);
//			s1.set_alloc_hubs(hubs);
//			neighbors.push_back(s1);
//		}
//	}
//	for(unsigned i = 0; i < neighbors.size(); i++){
//		neighbors[i].assign_hubs();
//		neighbors[i].route_traffics();
//	}
//	/*for(int i = 0; i < neighbors.size(); i++){
////		if(p_sol.is_hub(i)) continue;
//		neighbors[i].show_data();
////		printf("\nneighbor #%d: %.2lf\n", i, neighbors[i].get_total_cost());
//	}*/
//
//	return neighbors;
//}
//
//vector<solution> grasp::r_neighborhood2( solution& p_sol ){
//	/**
//	 * Generate a Neighborhood based on a two-point exchange of two random Hubs
//	 * in p_sol (partial solution) given.
//	 */
//
//	// Generating Neighborhood
//	vector< solution > neighbors;
//	for(int h = 0; h < p; h++){
//		vector< int > hubs(p_sol.get_alloc_hubs());
//		hubs[h] = rand() % mesh.size();
//
//		solution s1(instance, p, r);
//		s1.set_alloc_hubs(hubs);
//		s1.assign_hubs();
//		s1.route_traffics();
//		neighbors.push_back(s1);
//	}
//
//	return neighbors;
//}

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
	for(unsigned i = 0; i < aux.size() - x; i++)
		mesh.push_back(aux[i].second);

	this->mesh = mesh;
}

solution grasp::path_relinking(solution& origin, solution& destination){
	vector<int> hubs_dif;
	vector<int> hubs_o = origin.get_alloc_hubs();
	for(int i = 0; i < p; i++){
		if( destination.is_hub(hubs_o[i]) ) continue;
		hubs_dif.push_back(hubs_o[i]);
	}

	if(hubs_dif.size() > 1){
		int k = 0;
		vector<int> hubs_d(destination.get_alloc_hubs());
		for(int i = 0; i < p; i++){
			if( origin.is_hub(hubs_d[i]) ) continue;

			solution s1(instance, p, r);
			hubs_d[i] = hubs_dif[k];
			k++;

			s1.set_alloc_hubs(hubs_d);
			s1.assign_hubs();
			s1.route_traffics();
			solution improved = local_search_rn1(s1);

			if( improved.get_cost() < destination.get_cost() ) return improved;
		}
	}

	return destination;
}

solution& grasp::execute(){
	// Pre-processing
	preprocessing();

	// processing
	vector< double > partial_improment;
	//	int t_max = 0.1 * max_iterations;
	//	double mean = 0.0;
	unsigned i = 1, k = 0;
	bool first = true;
	while(i < max_iterations){
		solution initial = greedy_randomized_construction();
		solution improved = local_search_c2n1(initial);
//		solution improved = local_search_rn1(initial);

		// Acceptance criterion
		if(!first){
			if(improved.get_total_cost() < best.get_total_cost())
			{
				set_best(improved);
				i = 1;
			}/*else{
				improved = path_relinking(improved, best);
				if(improved.get_total_cost() < best.get_total_cost()){
					set_best(improved);
					i = 1;
					path.push_back(k);
				}else i++;
			}*/else i++;
		}else{
			best = improved;
			first = false;
		}

		it_log.push_back(make_pair(best.get_total_cost(), k++));
		times.push_back(((double) timer.getMilliSpan() / 1000));
	}


	// Post-processing
	//	set_best(local_search_rn1(best));

	return best;
}
