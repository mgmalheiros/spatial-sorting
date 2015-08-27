#include <cmath>
#include <cstdlib>
#include <list>
#include <vector>
#include <fstream>

#include <time.h>

#include <gsl/gsl_statistics.h>

// CGAL NOTES:
// - kd-incremental neighbor search takes more time (and not less)

// -------- k-d tree search --------
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
//#include <CGAL/Orthogonal_incremental_neighbor_search.h>
#include <CGAL/Search_traits_2.h>

// -------- nanoflann k-d tree --------
#include "nanoflann.hpp"

#include "common.hpp"

/*-------------------------------- IMPORTED VARIABLES --------------------------------*/

extern int dim_x, dim_y;
extern Position *position_matrix;
extern int position_counter;

extern bool *row_is_sorted;
extern bool *col_is_sorted;

extern bool verbose;
extern double tcomp, tread, twrit;

extern float diagonal_angle;
//extern float diagonal_width;

/*-------------------------------- LOCAL TYPES --------------------------------*/

// -------- cgal types --------

typedef CGAL::Simple_cartesian<float> Kernel;
typedef Kernel::Point_2 Point_2;
typedef CGAL::Search_traits_2<Kernel> Traits;
typedef CGAL::Orthogonal_k_neighbor_search<Traits> Neighbor_search;
//typedef CGAL::Orthogonal_incremental_neighbor_search<TreeTraits> Incremental_neighbor_search;

// -------- nanoflann types --------

struct Adaptor
{
	// Must return the number of data points
	inline size_t kdtree_get_point_count() const { return position_counter; }

	// Returns the distance between the vector "p1[0:size-1]" and the data point with index "idx_p2" stored in the class:
	inline float kdtree_distance(const float *p1, const size_t idx_p2, UNUSED size_t size) const
	{
		//const float d0=p1[0]-pts[idx_p2].x;
		//const float d1=p1[1]-pts[idx_p2].y;
		//const float d2=p1[2]-pts[idx_p2].z;
		const float d0 = p1[0] - position_matrix[idx_p2].x;
		const float d1 = p1[1] - position_matrix[idx_p2].y;
		return d0 * d0 + d1 * d1;
	}

	// Returns the dim'th component of the idx'th point in the class:
	// Since this is inlined and the "dim" argument is typically an immediate value, the
	//  "if/else's" are actually solved at compile time.
	inline float kdtree_get_pt(const size_t idx, int dim) const
	{
		if (dim==0) return position_matrix[idx].x;
		else return position_matrix[idx].y;
	}

	// Optional bounding-box computation: return false to default to a standard bbox computation loop.
	//   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
	//   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
	template <class BBOX>
	bool kdtree_get_bbox(UNUSED BBOX &bb) const { return false; }
};

// construct a kd-tree index:
typedef nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<float, Adaptor>, Adaptor, 2 /* dim */> NanoTree;

/*-------------------------------- LOCAL VARIABLES --------------------------------*/

float *nearest_exact;
float *nearest_trial;

Neighbor_search::Tree *cgal_tree = NULL;

NanoTree *nano_tree = NULL;

/*-------------------------------- CGAL KD-TREE FUNCTIONS --------------------------------*/

void setup_kd()
{
	cgal_tree = new Neighbor_search::Tree();
	//Incremental_neighbor_search::Tree cgal_tree;
	for (int i = 0; i < position_counter; i++) {
		cgal_tree->insert(Point_2(position_matrix[i].x, position_matrix[i].y));
	}
	//cgal_tree.statistics(std::cout);
	//std::cout << "cgal_tree=" << cgal_tree.size() << "\n";
}

void query_kd_exact(int k)
{
	for (int i = 0; i < position_counter; i++) {
		Point_2 query(position_matrix[i].x, position_matrix[i].y);
		//std::cout << std::setw(4) << i << ":kd  " << position_matrix[i].x << ',' << position_matrix[i].y << " -> ";

		Neighbor_search search(*cgal_tree, query, k + 1);
		//Incremental_neighbor_search search(cgal_tree, query);

		Neighbor_search::iterator it = search.begin();
		it++; // skip self
		for (int j = 0; j < k; j++) {
			nearest_exact[i * k + j] = it->second; // store squared distance
			it++;
		}
	}
}

void query_kd_approximate(int k, float eps)
{
	for (int i = 0; i < position_counter; i++) {
		Point_2 query(position_matrix[i].x, position_matrix[i].y);

		Neighbor_search search(*cgal_tree, query, k + 1, eps);

		Neighbor_search::iterator it = search.begin();
		it++; // skip self
		for (int j = 0; j < k; j++) {
			nearest_trial[i * k + j] = it->second; // store squared distance
			it++;
		}
	}
}

void finish_kd()
{
	delete cgal_tree;
	cgal_tree = NULL;
}

/*-------------------------------- NANOFLANN KD-TREE FUNCTIONS --------------------------------*/

void setup_nanoflann()
{
	Adaptor dataset;

	nano_tree = new NanoTree(2 /*dim*/, dataset, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
	nano_tree->buildIndex();
}

void query_nanoflann(int k)
{
	float query_pt[2];
	std::vector<size_t> ret_index(k + 1);
	std::vector<float> out_dist_sqr(k + 1);

	// ----------------------------------------------------------------
	// knnSearch():  Perform a search for the N closest points
	// ----------------------------------------------------------------
	for (int i = 0; i < position_counter; i++) {
		query_pt[0] = position_matrix[i].x;
		query_pt[1] = position_matrix[i].y;
		nano_tree->knnSearch(&query_pt[0], k + 1, &ret_index[0], &out_dist_sqr[0]);

		for (int j = 0; j < k; j++) {
			nearest_trial[i * k + j] = out_dist_sqr[j + 1]; // store squared distance of nearest point
		}

	}

	// ----------------------------------------------------------------
	// radiusSearch():  Perform a search for the N closest points
	// ----------------------------------------------------------------
//	const float search_radius = static_cast<float>(0.000001); // NOTE: squared distance
//	std::vector<std::pair<size_t,float> >   ret_matches;
//
//	nanoflann::SearchParams params;
//	//params.sorted = false;
//
//	const size_t nMatches = index.radiusSearch(&query_pt[0],search_radius, ret_matches, params);
//
//	std::cout << "radiusSearch(): radius=" << search_radius << " -> " << nMatches << " matches\n";
//	for (size_t i=0;i<nMatches;i++)
//		std::cout << "idx["<< i << "]=" << ret_matches[i].first << " dist["<< i << "]=" << ret_matches[i].second << std::endl;
//	std::cout << "\n";
}

void finish_nanoflann()
{
	delete nano_tree;
	nano_tree = NULL;
}

/*-------------------------------- NEIGHBORHOOD MATRIX FUNCTIONS --------------------------------*/

void query_nm(neighbor_ptr func, int k)
{
	if (k == 1) {
		// single nearest neighbor case
		for (int y = 0; y < dim_y; y++) {
			for (int x = 0; x < dim_x; x++) {
				//std::cout << std::setw(4) << position_matrix[i].id << ":nm  " << position_matrix[i].x << ',' << position_matrix[i].y << " -> ";
				int i = x + y * dim_x;
				Position& current = position_matrix[i];
				int *candidate = (*func)(x, y);
				float min_dist = FLT_MAX;
				int j = 0;
				while (candidate[j] != -1) {
					Position& neighbor = position_matrix[candidate[j]];
					float dx = current.x - neighbor.x;
					float dy = current.y - neighbor.y;
					float sd = dx * dx + dy * dy;
					if (sd < min_dist) {
						min_dist = sd;
					}
					j++;
				}
				nearest_trial[current.id] = min_dist;
				//std::cout << "nm #" << position_matrix[i].id << ' ' << std::setw(5) << min_dist << '\n';
			}
		}
	}
	else {
		// general k-nearest neighbors case
		for (int y = 0; y < dim_y; y++) {
			for (int x = 0; x < dim_x; x++) {
				int i = x + y * dim_x;
				Position& current = position_matrix[i];
				float *nearest = &nearest_trial[current.id * k];
				int *candidate = (*func)(x, y);
				for (int m = 0; m < k; m++) {
					nearest[m] = FLT_MAX;
				}
				int j = 0;
				while (candidate[j] != -1) {
					Position& neighbor = position_matrix[candidate[j]];
					float dx = current.x - neighbor.x;
					float dy = current.y - neighbor.y;
					float sd = dx * dx + dy * dy;
					// insert distance if near enough
					if (sd < nearest[k - 1]) {
						int m;
						for (m = k - 1; m >= 1; m--) {
							if (sd < nearest[m - 1]) {
								nearest[m] = nearest[m - 1];
							}
							else {
								break;
							}
						}
						nearest[m] = sd;
					}
					j++;
				}
			}
		}
	}
}

/*-------------------------------- TEST FUNCTIONS --------------------------------*/

void evaluate_misses(int k, double& pos_count, double& total_count)
{
	pos_count = total_count = 0;
	for (int i = 0; i < position_counter; i++) {
		float *exact = &nearest_exact[i * k];
		float *trial = &nearest_trial[i * k];
		int count = 0;
		for (int j = 0; j < k; j++) {
			if (*exact == *trial) {
				exact++;
				trial++;
			}
			else {
				count++;
				exact++;
			}
		}
		if (count > 0) {
			pos_count++;
			total_count += count;
		}
	}
}

void test_knn(const char *dist_name,  function_ptr dist_func,
		      const char *setup_name, function_ptr setup_func,
		      int runs, int k, float eps, const char *type, bool nano)
{
	double setup_kde[runs], setup_kda[runs], query_kde[runs], query_kda[runs];
	double pmiss_kda[runs], tmiss_kda[runs];

	double setup_nf[runs], query_nf[runs], pmiss_nf[runs], tmiss_nf[runs];

	double setup_nm[runs];
	double query_nm8[runs], query_nm24[runs], query_nm48[runs], query_nm80[runs], query_nm120[runs];
	double pmiss_nm8[runs], pmiss_nm24[runs], pmiss_nm48[runs], pmiss_nm80[runs], pmiss_nm120[runs];
	double tmiss_nm8[runs], tmiss_nm24[runs], tmiss_nm48[runs], tmiss_nm80[runs], tmiss_nm120[runs];

	double smean, ssdev, qmean, qsdev, pmmean, tmmean;
	struct timespec t0, t1;

	neighbor_ptr n8_func, n24_func, n48_func, n80_func, n120_func;
	if (strcmp(type, "hard") == 0) {
		n8_func  = &get_hard_8_neighborhood;
		n24_func = &get_hard_24_neighborhood;
		n48_func = &get_hard_48_neighborhood;
		n80_func = &get_hard_80_neighborhood;
		n120_func = &get_hard_120_neighborhood;
	}
	else {
		n8_func  = &get_soft_8_neighborhood;
		n24_func = &get_soft_24_neighborhood;
		n48_func = &get_soft_48_neighborhood;
		n80_func = &get_soft_80_neighborhood;
		n120_func = &get_soft_120_neighborhood;
	}

	for (int r = 1; r <= runs; r++) {
		random_seed(r);
		tcomp = tread = twrit = 0;

		(*dist_func)();

		/*---------------- evaluate exact kd-cgal_tree ----------------*/
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t0);
		setup_kd();
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t1);
		setup_kde[r - 1] = (t1.tv_sec - t0.tv_sec) * 1000.0 + (t1.tv_nsec - t0.tv_nsec) * 0.000001;

		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t0);
		query_kd_exact(k);
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t1);
		query_kde[r - 1] = (t1.tv_sec - t0.tv_sec) * 1000.0 + (t1.tv_nsec - t0.tv_nsec) * 0.000001;

		finish_kd();

		/*---------------- evaluate approximate kd-cgal_tree ----------------*/
		if (eps != 0) {
			clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t0);
			setup_kd();
			clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t1);
			setup_kda[r - 1] = (t1.tv_sec - t0.tv_sec) * 1000.0 + (t1.tv_nsec - t0.tv_nsec) * 0.000001;

			clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t0);
			query_kd_approximate(k, eps);
			clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t1);
			query_kda[r - 1] = (t1.tv_sec - t0.tv_sec) * 1000.0 + (t1.tv_nsec - t0.tv_nsec) * 0.000001;
			evaluate_misses(k, pmiss_kda[r - 1], tmiss_kda[r - 1]);

			finish_kd();
		}

		if (nano) {
			/*---------------- evaluate nanoflann kd-cgal_tree ----------------*/
			clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t0);
			setup_nanoflann();
			clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t1);
			setup_nf[r - 1] = (t1.tv_sec - t0.tv_sec) * 1000.0 + (t1.tv_nsec - t0.tv_nsec) * 0.000001;

			clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t0);
			query_nanoflann(k);
			clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t1);
			query_nf[r - 1] = (t1.tv_sec - t0.tv_sec) * 1000.0 + (t1.tv_nsec - t0.tv_nsec) * 0.000001;
			evaluate_misses(k, pmiss_nf[r - 1], tmiss_nf[r - 1]);

			finish_nanoflann();
		}

		/*---------------- evaluate neighborhood matrix ----------------*/
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t0);
		(*setup_func)();
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t1);
		setup_nm[r - 1] = (t1.tv_sec - t0.tv_sec) * 1000.0 + (t1.tv_nsec - t0.tv_nsec) * 0.000001;

		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t0);
		query_nm(n8_func, k);
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t1);
		query_nm8[r - 1] = (t1.tv_sec - t0.tv_sec) * 1000.0 + (t1.tv_nsec - t0.tv_nsec) * 0.000001;
		evaluate_misses(k, pmiss_nm8[r - 1], tmiss_nm8[r - 1]);

		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t0);
		query_nm(n24_func, k);
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t1);
		query_nm24[r - 1] = (t1.tv_sec - t0.tv_sec) * 1000.0 + (t1.tv_nsec - t0.tv_nsec) * 0.000001;
		evaluate_misses(k, pmiss_nm24[r - 1], tmiss_nm24[r - 1]);

		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t0);
		query_nm(n48_func, k);
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t1);
		query_nm48[r - 1] = (t1.tv_sec - t0.tv_sec) * 1000.0 + (t1.tv_nsec - t0.tv_nsec) * 0.000001;
		evaluate_misses(k, pmiss_nm48[r - 1], tmiss_nm48[r - 1]);

		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t0);
		query_nm(n80_func, k);
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t1);
		query_nm80[r - 1] = (t1.tv_sec - t0.tv_sec) * 1000.0 + (t1.tv_nsec - t0.tv_nsec) * 0.000001;
		evaluate_misses(k, pmiss_nm80[r - 1], tmiss_nm80[r - 1]);

		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t0);
		query_nm(n120_func, k);
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t1);
		query_nm120[r - 1] = (t1.tv_sec - t0.tv_sec) * 1000.0 + (t1.tv_nsec - t0.tv_nsec) * 0.000001;
		evaluate_misses(k, pmiss_nm120[r - 1], tmiss_nm120[r - 1]);
	}

	smean = gsl_stats_mean(setup_kde, 1, runs);
	ssdev = gsl_stats_sd(setup_kde, 1, runs);
	qmean = gsl_stats_mean(query_kde, 1, runs);
	qsdev = gsl_stats_sd(query_kde, 1, runs);
	std::cout << dist_name  << "  cgal      exact    ";

	std::cout << std::setprecision(2) << std::fixed << ' ';
	std::cout << std::setw(10) << smean << "  " << std::setw(10) << ssdev << "  ";
	std::cout << std::setw(10) << qmean << "  " << std::setw(10) << qsdev << "  ";
	std::cout << std::setw(10) << smean + qmean << "           0           0           0           0\n";

	if (eps != 0) {
		smean = gsl_stats_mean(setup_kda, 1, runs);
		ssdev = gsl_stats_sd(setup_kda, 1, runs);
		qmean = gsl_stats_mean(query_kda, 1, runs);
		qsdev = gsl_stats_sd(query_kda, 1, runs);
		pmmean = gsl_stats_mean(pmiss_kda, 1, runs);
		tmmean = gsl_stats_mean(tmiss_kda, 1, runs);
		std::cout << dist_name  << "  cgal      approx   ";
		std::cout << std::setprecision(2) << std::fixed << ' ';
		std::cout << std::setw(10) << smean << "  " << std::setw(10) << ssdev << "  ";
		std::cout << std::setw(10) << qmean << "  " << std::setw(10) << qsdev << "  ";
		std::cout << std::setw(10) << smean + qmean << "  ";
		std::cout << std::setprecision(0) << std::setw(10) << pmmean << "  ";
		std::cout << std::setprecision(4) << std::setw(10) << 100 * pmmean / position_counter << "  ";
		std::cout << std::setprecision(0) << std::setw(10) << tmmean << "  ";
		std::cout << std::setprecision(4) << std::setw(10) << 100 * tmmean / (position_counter * k) << '\n';
	}

	if (nano) {
		// nanoflann k-d tree
		smean = gsl_stats_mean(setup_nf, 1, runs);
		ssdev = gsl_stats_sd(setup_nf, 1, runs);
		qmean = gsl_stats_mean(query_nf, 1, runs);
		qsdev = gsl_stats_sd(query_nf, 1, runs);
		pmmean = gsl_stats_mean(pmiss_nf, 1, runs);
		tmmean = gsl_stats_mean(tmiss_nf, 1, runs);
		std::cout << dist_name  << "  nanoflann          ";
		std::cout << std::setprecision(2) << std::fixed << ' ';
		std::cout << std::setw(10) << smean << "  " << std::setw(10) << ssdev << "  ";
		std::cout << std::setw(10) << qmean << "  " << std::setw(10) << qsdev << "  ";
		std::cout << std::setw(10) << smean + qmean << "  ";
		std::cout << std::setprecision(0) << std::setw(10) << pmmean << "  ";
		std::cout << std::setprecision(4) << std::setw(10) << 100 * pmmean / position_counter << "  ";
		std::cout << std::setprecision(0) << std::setw(10) << tmmean << "  ";
		std::cout << std::setprecision(4) << std::setw(10) << 100 * tmmean / (position_counter * k) << '\n';
	}

	// neighborhood matrix
	smean = gsl_stats_mean(setup_nm, 1, runs);
	ssdev = gsl_stats_sd(setup_nm, 1, runs);

	qmean = gsl_stats_mean(query_nm8, 1, runs);
	qsdev = gsl_stats_sd(query_nm8, 1, runs);
	pmmean = gsl_stats_mean(pmiss_nm8, 1, runs);
	tmmean = gsl_stats_mean(tmiss_nm8, 1, runs);
	std::cout << dist_name  << "  " << setup_name << "  " << type << "_8   ";
	std::cout << std::setprecision(2) << std::fixed << ' ';
	std::cout << std::setw(10) << smean << "  " << std::setw(10) << ssdev << "  ";
	std::cout << std::setw(10) << qmean << "  " << std::setw(10) << qsdev << "  ";
	std::cout << std::setw(10) << smean + qmean << "  ";
	std::cout << std::setprecision(0) << std::setw(10) << pmmean << "  ";
	std::cout << std::setprecision(4) << std::setw(10) << 100 * pmmean / position_counter << "  ";
	std::cout << std::setprecision(0) << std::setw(10) << tmmean << "  ";
	std::cout << std::setprecision(4) << std::setw(10) << 100 * tmmean / (position_counter * k) << '\n';

	qmean = gsl_stats_mean(query_nm24, 1, runs);
	qsdev = gsl_stats_sd(query_nm24, 1, runs);
	pmmean = gsl_stats_mean(pmiss_nm24, 1, runs);
	tmmean = gsl_stats_mean(tmiss_nm24, 1, runs);
	std::cout << dist_name  << "  " << setup_name << "  " << type << "_24  ";
	std::cout << std::setprecision(2) << std::fixed << ' ';
	std::cout << std::setw(10) << smean << "  " << std::setw(10) << ssdev << "  ";
	std::cout << std::setw(10) << qmean << "  " << std::setw(10) << qsdev << "  ";
	std::cout << std::setw(10) << smean + qmean << "  ";
	std::cout << std::setprecision(0) << std::setw(10) << pmmean << "  ";
	std::cout << std::setprecision(4) << std::setw(10) << 100 * pmmean / position_counter << "  ";
	std::cout << std::setprecision(0) << std::setw(10) << tmmean << "  ";
	std::cout << std::setprecision(4) << std::setw(10) << 100 * tmmean / (position_counter * k) << '\n';

	qmean = gsl_stats_mean(query_nm48, 1, runs);
	qsdev = gsl_stats_sd(query_nm48, 1, runs);
	pmmean = gsl_stats_mean(pmiss_nm48, 1, runs);
	tmmean = gsl_stats_mean(tmiss_nm48, 1, runs);
	std::cout << dist_name  << "  " << setup_name << "  " << type << "_48  ";
	std::cout << std::setprecision(2) << std::fixed << ' ';
	std::cout << std::setw(10) << smean << "  " << std::setw(10) << ssdev << "  ";
	std::cout << std::setw(10) << qmean << "  " << std::setw(10) << qsdev << "  ";
	std::cout << std::setw(10) << smean + qmean << "  ";
	std::cout << std::setprecision(0) << std::setw(10) << pmmean << "  ";
	std::cout << std::setprecision(4) << std::setw(10) << 100 * pmmean / position_counter << "  ";
	std::cout << std::setprecision(0) << std::setw(10) << tmmean << "  ";
	std::cout << std::setprecision(4) << std::setw(10) << 100 * tmmean / (position_counter * k) << '\n';

	qmean = gsl_stats_mean(query_nm80, 1, runs);
	qsdev = gsl_stats_sd(query_nm80, 1, runs);
	pmmean = gsl_stats_mean(pmiss_nm80, 1, runs);
	tmmean = gsl_stats_mean(tmiss_nm80, 1, runs);
	std::cout << dist_name  << "  " << setup_name << "  " << type << "_80  ";
	std::cout << std::setprecision(2) << std::fixed << ' ';
	std::cout << std::setw(10) << smean << "  " << std::setw(10) << ssdev << "  ";
	std::cout << std::setw(10) << qmean << "  " << std::setw(10) << qsdev << "  ";
	std::cout << std::setw(10) << smean + qmean << "  ";
	std::cout << std::setprecision(0) << std::setw(10) << pmmean << "  ";
	std::cout << std::setprecision(4) << std::setw(10) << 100 * pmmean / position_counter << "  ";
	std::cout << std::setprecision(0) << std::setw(10) << tmmean << "  ";
	std::cout << std::setprecision(4) << std::setw(10) << 100 * tmmean / (position_counter * k) << '\n';

	qmean = gsl_stats_mean(query_nm120, 1, runs);
	qsdev = gsl_stats_sd(query_nm120, 1, runs);
	pmmean = gsl_stats_mean(pmiss_nm120, 1, runs);
	tmmean = gsl_stats_mean(tmiss_nm120, 1, runs);
	std::cout << dist_name  << "  " << setup_name << "  " << type << "_120 ";
	std::cout << std::setprecision(2) << std::fixed << ' ';
	std::cout << std::setw(10) << smean << "  " << std::setw(10) << ssdev << "  ";
	std::cout << std::setw(10) << qmean << "  " << std::setw(10) << qsdev << "  ";
	std::cout << std::setw(10) << smean + qmean << "  ";
	std::cout << std::setprecision(0) << std::setw(10) << pmmean << "  ";
	std::cout << std::setprecision(4) << std::setw(10) << 100 * pmmean / position_counter << "  ";
	std::cout << std::setprecision(0) << std::setw(10) << tmmean << "  ";
	std::cout << std::setprecision(4) << std::setw(10) << 100 * tmmean / (position_counter * k) << '\n';
}

int main()
{
	int side = 1000;
	int runs = 10;
	int k = 1;
	float eps = 0.0; //1.0; //0.5;
	const char *ne_type = "hard";
	//const char *ne_type = "soft";

	verbose = false;
	dim_x = dim_y = side;
	position_matrix = new Position[dim_x * dim_y];
	position_counter = dim_x * dim_y;

    row_is_sorted = new bool[dim_y];
    col_is_sorted = new bool[dim_x];

    nearest_exact = new float[dim_x * dim_y * k];
    nearest_trial  = new float[dim_x * dim_y * k];

	std::cout << "\nside=" << side << " runs=" << runs << " k=" << k << " eps=" << eps << "\n\n";
	std::cout << "distrib   setup     query     prepr_mean  prepr_sdev  query_mean  query_sdev  total_mean   pos_miss#   pos_miss% total_miss# total_miss%\n";
    //            uniform   cgal      exact           2.75        -nan      903.18        -nan      905.93           0           0           0           0

//	const char *sort_name;
//	function_ptr sort_func;

//	sort_name = "sp_sh_sk";
//	sort_func = &full_spatial_shell_sort_skip;
//	sort_name = "sqi     ";
//	sort_func = &full_quick_insertion_skip;

	//test_knn("uniform ", &setup_uniform,         "part_qk ", &partial_quicksort,            runs, k, eps, ne_type, false);
	//test_knn("uniform ", &setup_uniform,         "soe     ", &odd_even_full_sort,           runs, k, eps, ne_type, false);
	test_knn("uniform ", &setup_uniform,         "sqi     ", &full_quick_insertion_skip,    runs, k, eps, ne_type, true);
	test_knn("uniform ", &setup_uniform,         "ssi     ", &full_spatial_shell_sort_skip, runs, k, eps, ne_type, false);
	test_knn("uniform ", &setup_uniform,         "simple  ", &simple_sort,                  runs, k, eps, ne_type, false);

//	test_knn("uniform ", &setup_uniform,         sort_name, sort_func, runs, k, eps, ne_type, false);
//	test_knn("compact ", &setup_compact,         sort_name, sort_func, runs, k, eps, ne_type, false);
//	test_knn("gradient", &setup_gradient,        sort_name, sort_func, runs, k, eps, ne_type, false);
//	test_knn("gaussian", &setup_gaussian,        sort_name, sort_func, runs, k, eps, ne_type, false);
//	test_knn("shuffle ", &setup_shuffle,         sort_name, sort_func, runs, k, eps, ne_type, false);
//	test_knn("inside  ", &setup_inside_circle,   sort_name, sort_func, runs, k, eps, ne_type, false);
//	test_knn("outside ", &setup_outside_circle,  sort_name, sort_func, runs, k, eps, ne_type, false);
//	test_knn("clusters", &setup_random_clusters, sort_name, sort_func, runs, k, eps, ne_type, false);
//	test_knn("holes   ", &setup_random_holes,    sort_name, sort_func, runs, k, eps, ne_type, false);
//	test_knn("streets ", &setup_streets,         sort_name, sort_func, runs, k, eps, ne_type, false);
//	test_knn("tilted  ", &setup_tilted,          sort_name, sort_func, runs, k, eps, ne_type, false);
//	test_knn("triangle", &setup_triangle,        sort_name, sort_func, runs, k, eps, ne_type, false);

	diagonal_angle = 20;
//	test_knn("diagon20", &setup_diagonal,        sort_name, sort_func, runs, k, eps, ne_type, false);
//	diagonal_angle = 45;
//	test_knn("diagon45", &setup_diagonal,        sort_name, sort_func, runs, k, eps, ne_type, false);
//	test_knn("cross   ", &setup_cross,           sort_name, sort_func, runs, k, eps, ne_type, false);

	return 0;
}

