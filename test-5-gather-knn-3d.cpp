#include <cmath>
#include <cstdlib>
#include <list>
#include <vector>

#include <sys/time.h>

#include <gsl/gsl_statistics.h>

// -------- k-d tree search --------
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
//#include <CGAL/Orthogonal_incremental_neighbor_search.h>
#include <CGAL/Search_traits_3.h>

#include "common.hpp"

/*-------------------------------- IMPORTED VARIABLES --------------------------------*/

extern int dim_x, dim_y, dim_z;
extern Position3d *position_matrix_3d;
extern int position_counter;

extern bool *row_is_sorted;
extern bool *col_is_sorted;

extern bool verbose;
extern double tcomp, tread, twrit;

extern float diagonal_angle;
//extern float diagonal_width;

/*-------------------------------- LOCAL TYPES --------------------------------*/

typedef CGAL::Simple_cartesian<float> Kernel;
typedef Kernel::Point_3 Point_3;
typedef CGAL::Search_traits_3<Kernel> Traits;
typedef CGAL::Orthogonal_k_neighbor_search<Traits> Neighbor_search;
//typedef CGAL::Orthogonal_incremental_neighbor_search<TreeTraits> Incremental_neighbor_search;

/*-------------------------------- LOCAL VARIABLES --------------------------------*/

float *range_kd_exact;
float *range_kd_approximate;
float *range_nm;

Neighbor_search::Tree *cgal_tree = NULL;

int *candidate;

/*-------------------------------- KD-TREE FUNCTIONS --------------------------------*/

void setup_kd()
{
	cgal_tree = new Neighbor_search::Tree();
	//Incremental_neighbor_search::Tree tree;
	for (int i = 0; i < position_counter; i++) {
		cgal_tree->insert(Point_3(position_matrix_3d[i].x, position_matrix_3d[i].y, position_matrix_3d[i].z));
	}
	//tree.statistics(std::cout);
	//std::cout << "tree=" << tree.size() << "\n";
}

void query_kd_exact(int k)
{
	for (int i = 0; i < position_counter; i++) {
		Point_3 query(position_matrix_3d[i].x, position_matrix_3d[i].y, position_matrix_3d[i].z);
		//std::cout << std::setw(4) << i << ":kd  " << position_matrix[i].x << ',' << position_matrix[i].y << " -> ";

		Neighbor_search search(*cgal_tree, query, k + 1);
		//Incremental_neighbor_search search(tree, query);

		Neighbor_search::iterator it = search.begin();
		// skip self and other neighbors
//		if (i % 10000 == 0) {
//			std::cout << "kd #" << i << ' ';
//			for (int j = 1; j <= k; j++) {
//				std::cout << std::setw(5) << it->second << ' ';
//				it++;
//			}
//			std::cout << std::setw(5) << it->second << '\n';
//		}
//		else {
			for (int j = 1; j <= k; j++) {
				it++;
			}
//		}
		range_kd_exact[i] = it->second; // store squared distance of nearest point

		//while (it != search.end()) {
		//nearest_kd[i].x = it->first.x();
		//nearest_kd[i].y = it->first.y();
		//std::cout << std::setw(5) << it->first.x() << ',' << std::setw(5) << it->first.y();
		//std::cout << " d=" << std::setw(5) << it->second << '\n';
		//it++;
		//}
		//std::cout << "\n";
		//search.statistics(std::cout);
	}
}

void query_kd_approximate(int k, float eps)
{
	for (int i = 0; i < position_counter; i++) {
		Point_3 query(position_matrix_3d[i].x, position_matrix_3d[i].y, position_matrix_3d[i].z);

		Neighbor_search search(*cgal_tree, query, k + 1, eps);

		Neighbor_search::iterator it = search.begin();
		for (int j = 1; j <= k; j++) {
			it++;
		}
		range_kd_approximate[i] = it->second; // store squared distance of nearest point
	}
}

void finish_kd()
{
	delete cgal_tree;
	cgal_tree = NULL;
}

/*-------------------------------- NEIGHBORHOOD MATRIX FUNCTIONS --------------------------------*/

void query_nm(int k, int n_size)
{
	if (k == 1) {
		// single nearest neighbor case
		for (int z = 0; z < dim_z; z++) {
			for (int y = 0; y < dim_y; y++) {
				for (int x = 0; x < dim_x; x++) {
					//std::cout << std::setw(4) << position_matrix[i].id << ":nm  " << position_matrix[i].x << ',' << position_matrix[i].y << " -> ";
					int i = x + y * dim_x + z * dim_x * dim_y;
					Position3d& current = position_matrix_3d[i];
					get_hard_3d_neighborhood(x, y, z, n_size, candidate);
					float min_dist = FLT_MAX;
					int j = 0;
					while (candidate[j] != -1) {
						Position3d& neighbor = position_matrix_3d[candidate[j]];
						float dx = current.x - neighbor.x;
						float dy = current.y - neighbor.y;
						float dz = current.z - neighbor.z;
						float sd = dx * dx + dy * dy + dz * dz;
						if (sd < min_dist) {
							min_dist = sd;
						}
						j++;
					}
					range_nm[i] = min_dist;
					//std::cout << "nm #" << position_matrix_3d[i].id << ' ' << std::setw(5) << min_dist << '\n';
				}
			}
		}
	}
	else {
		// general k-nearest neighbors case
//		float nearest[k];
//		for (int y = 0; y < dim_y; y++) {
//			for (int x = 0; x < dim_x; x++) {
//				int i = x + y * dim_x;
//				Position& current = position_matrix[i];
//				int *candidate = (*func)(x, y);
//				for (int m = 0; m < k; m++) {
//					nearest[m] = FLT_MAX;
//				}
//				int j = 0;
//				while (candidate[j] != -1) {
//					Position& neighbor = position_matrix[candidate[j]];
//					float dx = current.x - neighbor.x;
//					float dy = current.y - neighbor.y;
//					float sd = dx * dx + dy * dy;
//					// insert distance if near enough
//					if (sd < nearest[k - 1]) {
//						int m;
//						for (m = k - 1; m >= 1; m--) {
//							if (sd < nearest[m - 1]) {
//								nearest[m] = nearest[m - 1];
//							}
//							else {
//								break;
//							}
//						}
//						nearest[m] = sd;
//					}
//					j++;
//				}
//				range_nm[i] = nearest[k - 1];
////				if (position_matrix[i].id % 10000 == 0) {
////					std::cout << "nm #" << position_matrix[i].id << ' ';
////					for (int n = 0; n < k; n++) {
////						std::cout << std::setw(5) << nearest[n] << ' ';
////					}
////					std::cout << '\n';
////				}
//			}
//		}
	}
}

/*-------------------------------- TEST FUNCTIONS --------------------------------*/

void evaluate_kd_misses(double& count, double& error_sum)
{
	count = error_sum = 0;
	for (int i = 0; i < position_counter; i++) {
		if (range_kd_approximate[i] != range_kd_exact[i]) {
			//std::cout << '#' << pos << " not matched: kd=" << nearest_kd[pos] << " <> nm=" << nearest_nm[i] << '\n';
			count += 1;
			error_sum += fabsf(sqrtf(range_kd_approximate[i]) - sqrtf(range_kd_exact[i]));
		}
	}
	//std::cout << "not matched: " << c << '\n';
}

void evaluate_nm_misses(double& count, double& error_sum)
{
	count = error_sum = 0;
	for (int i = 0; i < position_counter; i++) {
		int pos = position_matrix_3d[i].id;
		if (range_nm[i] != range_kd_exact[pos]) {
			//std::cout << '#' << pos << " not matched: kd=" << nearest_kd[pos] << " <> nm=" << nearest_nm[i] << '\n';
			count += 1;
			error_sum += fabsf(sqrtf(range_nm[i]) - sqrtf(range_kd_exact[pos]));
		}
	}
	//std::cout << "not matched: " << c << '\n';
}

void test_knn(const char *dist_name,  function_ptr dist_func,
		      const char *setup_name, function_ptr setup_func,
		      int runs, int k, float eps)
{
	double setup_kde[runs], setup_kda[runs], query_kde[runs], query_kda[runs];
	double misse_kda[runs], error_kda[runs];

	double setup_nm[runs];
	double query_nm1[runs], query_nm2[runs], query_nm3[runs], query_nm4[runs];
	double misse_nm1[runs], misse_nm2[runs], misse_nm3[runs], misse_nm4[runs];
	double error_nm1[runs], error_nm2[runs], error_nm3[runs], error_nm4[runs];

	double smean, ssdev, qmean, qsdev, mmean, emean;
	struct timeval t0, t1;

	for (int r = 1; r <= runs; r++) {
		random_seed(r);
		tcomp = tread = twrit = 0;

		(*dist_func)();

		/*---------------- evaluate exact kd-tree ----------------*/
		gettimeofday(&t0, 0);
		setup_kd();
		gettimeofday(&t1, 0);
		setup_kde[r - 1] = (t1.tv_sec - t0.tv_sec) * 1000000.0 + (t1.tv_usec - t0.tv_usec);

		gettimeofday(&t0, 0);
		query_kd_exact(k);
		gettimeofday(&t1, 0);
		query_kde[r - 1] = (t1.tv_sec - t0.tv_sec) * 1000000.0 + (t1.tv_usec - t0.tv_usec);

		finish_kd();

		/*---------------- evaluate approximate kd-tree ----------------*/
		if (eps != 0) {
			gettimeofday(&t0, 0);
			setup_kd();
			gettimeofday(&t1, 0);
			setup_kda[r - 1] = (t1.tv_sec - t0.tv_sec) * 1000000.0 + (t1.tv_usec - t0.tv_usec);

			gettimeofday(&t0, 0);
			query_kd_approximate(k, eps);
			gettimeofday(&t1, 0);
			query_kda[r - 1] = (t1.tv_sec - t0.tv_sec) * 1000000.0 + (t1.tv_usec - t0.tv_usec);
			evaluate_kd_misses(misse_kda[r - 1], error_kda[r - 1]);

			finish_kd();
		}

		/*---------------- evaluate neighborhood matrix ----------------*/
		gettimeofday(&t0, 0);
		(*setup_func)();
		gettimeofday(&t1, 0);
		setup_nm[r - 1] = (t1.tv_sec - t0.tv_sec) * 1000000.0 + (t1.tv_usec - t0.tv_usec);

		gettimeofday(&t0, 0);
		query_nm(k, 1);
		gettimeofday(&t1, 0);
		query_nm1[r - 1] = (t1.tv_sec - t0.tv_sec) * 1000000.0 + (t1.tv_usec - t0.tv_usec);
		evaluate_nm_misses(misse_nm1[r - 1], error_nm1[r - 1]);

		gettimeofday(&t0, 0);
		query_nm(k, 2);
		gettimeofday(&t1, 0);
		query_nm2[r - 1] = (t1.tv_sec - t0.tv_sec) * 1000000.0 + (t1.tv_usec - t0.tv_usec);
		evaluate_nm_misses(misse_nm2[r - 1], error_nm2[r - 1]);

		gettimeofday(&t0, 0);
		query_nm(k, 3);
		gettimeofday(&t1, 0);
		query_nm3[r - 1] = (t1.tv_sec - t0.tv_sec) * 1000000.0 + (t1.tv_usec - t0.tv_usec);
		evaluate_nm_misses(misse_nm3[r - 1], error_nm3[r - 1]);

		gettimeofday(&t0, 0);
		query_nm(k, 4);
		gettimeofday(&t1, 0);
		query_nm4[r - 1] = (t1.tv_sec - t0.tv_sec) * 1000000.0 + (t1.tv_usec - t0.tv_usec);
		evaluate_nm_misses(misse_nm4[r - 1], error_nm4[r - 1]);
	}

	smean = gsl_stats_mean(setup_kde, 1, runs);
	ssdev = gsl_stats_sd(setup_kde, 1, runs);
	qmean = gsl_stats_mean(query_kde, 1, runs);
	qsdev = gsl_stats_sd(query_kde, 1, runs);
	std::cout << dist_name  << "  kd-tree   kd-exac  ";
	std::cout << std::setprecision(0) << std::fixed << ' ';
	std::cout << std::setw(10) << smean << "  " << std::setw(10) << ssdev << "  ";
	std::cout << std::setw(10) << qmean << "  " << std::setw(10) << qsdev << "  ";
	std::cout << std::setw(10) << smean + qmean << '\n';

	if (eps != 0) {
		smean = gsl_stats_mean(setup_kda, 1, runs);
		ssdev = gsl_stats_sd(setup_kda, 1, runs);
		qmean = gsl_stats_mean(query_kda, 1, runs);
		qsdev = gsl_stats_sd(query_kda, 1, runs);
		mmean = gsl_stats_mean(misse_kda, 1, runs);
		emean = gsl_stats_mean(error_kda, 1, runs);
		std::cout << dist_name  << "  kd-tree   kd-appr  ";
		std::cout << std::setprecision(0) << std::fixed << ' ';
		std::cout << std::setw(10) << smean << "  " << std::setw(10) << ssdev << "  ";
		std::cout << std::setw(10) << qmean << "  " << std::setw(10) << qsdev << "  ";
		std::cout << std::setw(10) << smean + qmean << "  ";
		std::cout << std::setw(10) << mmean << "  ";
		std::cout << std::setprecision(7) << std::setw(10) << 100 * mmean / position_counter << "%  ";
		std::cout << std::setw(10) << emean << "  " << std::setw(10) << emean / mmean << std::setprecision(4) << "  eps=" << eps << '\n';
	}

	smean = gsl_stats_mean(setup_nm, 1, runs);
	ssdev = gsl_stats_sd(setup_nm, 1, runs);

	qmean = gsl_stats_mean(query_nm1, 1, runs);
	qsdev = gsl_stats_sd(query_nm1, 1, runs);
	mmean = gsl_stats_mean(misse_nm1, 1, runs);
	emean = gsl_stats_mean(error_nm1, 1, runs);
	std::cout << dist_name  << "  " << setup_name << "  hard_26  ";
	std::cout << std::setprecision(0) << std::fixed << ' ';
	std::cout << std::setw(10) << smean << "  " << std::setw(10) << ssdev << "  ";
	std::cout << std::setw(10) << qmean << "  " << std::setw(10) << qsdev << "  ";
	std::cout << std::setw(10) << smean + qmean << "  ";
	std::cout << std::setw(10) << mmean << "  ";
	std::cout << std::setprecision(7) << std::setw(10) << 100 * mmean / position_counter << "%  ";
	std::cout << std::setw(10) << emean << "  " << std::setw(10) << emean / mmean << '\n';

	qmean = gsl_stats_mean(query_nm2, 1, runs);
	qsdev = gsl_stats_sd(query_nm2, 1, runs);
	mmean = gsl_stats_mean(misse_nm2, 1, runs);
	emean = gsl_stats_mean(error_nm2, 1, runs);
	std::cout << dist_name  << "  " << setup_name << "  hard_124 ";
	std::cout << std::setprecision(0) << std::fixed << ' ';
	std::cout << std::setw(10) << smean << "  " << std::setw(10) << ssdev << "  ";
	std::cout << std::setw(10) << qmean << "  " << std::setw(10) << qsdev << "  ";
	std::cout << std::setw(10) << smean + qmean << "  ";
	std::cout << std::setw(10) << mmean << "  ";
	std::cout << std::setprecision(7) << std::setw(10) << 100 * mmean / position_counter << "%  ";
	std::cout << std::setw(10) << emean << "  " << std::setw(10) << emean / mmean << '\n';

	qmean = gsl_stats_mean(query_nm3, 1, runs);
	qsdev = gsl_stats_sd(query_nm3, 1, runs);
	mmean = gsl_stats_mean(misse_nm3, 1, runs);
	emean = gsl_stats_mean(error_nm3, 1, runs);
	std::cout << dist_name  << "  " << setup_name << "  hard_342 ";
	std::cout << std::setprecision(0) << std::fixed << ' ';
	std::cout << std::setw(10) << smean << "  " << std::setw(10) << ssdev << "  ";
	std::cout << std::setw(10) << qmean << "  " << std::setw(10) << qsdev << "  ";
	std::cout << std::setw(10) << smean + qmean << "  ";
	std::cout << std::setw(10) << mmean << "  ";
	std::cout << std::setprecision(7) << std::setw(10) << 100 * mmean / position_counter << "%  ";
	std::cout << std::setw(10) << emean << "  " << std::setw(10) << emean / mmean << '\n';

	qmean = gsl_stats_mean(query_nm4, 1, runs);
	qsdev = gsl_stats_sd(query_nm4, 1, runs);
	mmean = gsl_stats_mean(misse_nm4, 1, runs);
	emean = gsl_stats_mean(error_nm4, 1, runs);
	std::cout << dist_name  << "  " << setup_name << "  hard_728 ";
	std::cout << std::setprecision(0) << std::fixed << ' ';
	std::cout << std::setw(10) << smean << "  " << std::setw(10) << ssdev << "  ";
	std::cout << std::setw(10) << qmean << "  " << std::setw(10) << qsdev << "  ";
	std::cout << std::setw(10) << smean + qmean << "  ";
	std::cout << std::setw(10) << mmean << "  ";
	std::cout << std::setprecision(7) << std::setw(10) << 100 * mmean / position_counter << "%  ";
	std::cout << std::setw(10) << emean << "  " << std::setw(10) << emean / mmean << '\n';
}

int main()
{
	int side = 64;
	int runs = 10;
	int k = 1;
	float eps = 0.0;

	verbose = false;
	dim_x = dim_y = dim_z = side;
	position_matrix_3d = new Position3d[dim_x * dim_y * dim_z];
	position_counter = dim_x * dim_y * dim_z;

    //row_is_sorted = new bool[dim_y];
    //col_is_sorted = new bool[dim_x];

	range_kd_exact       = new float[dim_x * dim_y * dim_z];
	range_kd_approximate = new float[dim_x * dim_y * dim_z];
	range_nm             = new float[dim_x * dim_y * dim_z];

	candidate = new int[9 * 9 * 9]; // FIXME: hard-coded maximum neighborhood size

	std::cout << "\nside=" << side << " runs=" << runs << " k=" << k << "\n\n";
	std::cout << "distrib   setup     query     prepr_mean  prepr_sdev  query_mean  query_sdev  total_mean   miss_mean        miss%   error_sum  error_mean\n";
                //uniform   sp_shell  neigh_4        90156        1035        3441         140       93598       99702  39.8808000%  68.9624678   0.0006917

	const char *sort_name;
	function_ptr sort_func;

	//test_knn("uniform ", &setup_uniform, "p_quicks", &partial_quicksort, runs);
	//test_knn("uniform ", &setup_uniform, "f_oddevn", &odd_even_full_sort, runs);

	sort_name  = "full_ins";
	sort_func = &full_insertion_3d;
	test_knn("uniform ", &setup_uniform_3d,         sort_name, sort_func, runs, k, eps);

//	sort_name  = "quick_sk";
//	sort_func = &full_quick_insertion_skip;
//	test_knn("uniform ", &setup_uniform,         sort_name, sort_func, runs, k, eps);
//
//	sort_name  = "sp_sh_sk";
//	sort_func = &full_spatial_shell_sort_skip;
//	test_knn("uniform ", &setup_uniform,         sort_name, sort_func, runs, k, eps);

//	test_knn("uniform ", &setup_uniform,         sort_name, sort_func, runs, k, eps);
//
//	test_knn("gradient", &setup_gradient,        sort_name, sort_func, runs, k, eps);
//	test_knn("gaussian", &setup_gaussian,        sort_name, sort_func, runs, k, eps);
//	test_knn("shuffle ", &setup_shuffle,         sort_name, sort_func, runs, k, eps);
//	test_knn("inside  ", &setup_inside_circle,   sort_name, sort_func, runs, k, eps);
//	test_knn("outside ", &setup_outside_circle,  sort_name, sort_func, runs, k, eps);
//	test_knn("clusters", &setup_random_clusters, sort_name, sort_func, runs, k, eps);
//	test_knn("holes   ", &setup_random_holes,    sort_name, sort_func, runs, k, eps);
//	test_knn("streets ", &setup_streets,         sort_name, sort_func, runs, k, eps);
//	test_knn("tilted  ", &setup_tilted,          sort_name, sort_func, runs, k, eps);
//
	diagonal_angle = 20;
//	test_knn("diagon20", &setup_diagonal,        sort_name, sort_func, runs, k, eps);
//	diagonal_angle = 45;
//	test_knn("diagon45", &setup_diagonal,        sort_name, sort_func, runs, k, eps);
//
//	test_knn("cross   ", &setup_cross,           sort_name, sort_func, runs, k, eps);
//	test_knn("triangle", &setup_triangle,        sort_name, sort_func, runs, k, eps);

	return 0;
}
