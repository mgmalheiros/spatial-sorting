#include <list>
#include <cmath>

#include <time.h>

#include <gsl/gsl_statistics.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Fuzzy_sphere.h>

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

typedef CGAL::Simple_cartesian<float> Kernel;
typedef Kernel::Point_2 Point_2;
typedef CGAL::Search_traits_2<Kernel> Traits;
typedef CGAL::Kd_tree<Traits> Tree;
typedef CGAL::Fuzzy_sphere<Traits> Fuzzy_circle;

/*-------------------------------- LOCAL VARIABLES --------------------------------*/

int *range_kd_exact;
int *range_kd_approximate;
int *range_nm;

Tree *cgal_tree = NULL;

/*-------------------------------- KD-TREE FUNCTIONS --------------------------------*/

void setup_kd()
{
	cgal_tree = new Tree();
	for (int i = 0; i < position_counter; i++) {
		cgal_tree->insert(Point_2(position_matrix[i].x, position_matrix[i].y));
	}
	//tree->statistics(std::cout);
	//std::cout << "tree=" << tree->size() << "\n";
}

void query_kd_exact(float range)
{
	for (int i = 0; i < position_counter; i++) {
		Point_2 center(position_matrix[i].x, position_matrix[i].y);
		Fuzzy_circle circle(center, range); // define exact circular range query (fuzziness=0)

		std::list<Point_2> result;
		cgal_tree->search(std::back_inserter(result), circle);
		range_kd_exact[i] = result.size() - 1; // skip self

//		std::cout << "kd #" << i << ' ';
//		int c = 0;
//		for (std::list<Point_2>::iterator it = result.begin(); it != result.end(); it++) {
//			if (*it == center) {
//				// std::cout << "self ";
//				continue;
//			}
//			// std::cout << std::setw(5) << it->x() << ',' << std::setw(5) << it->y() << ' ';
//			c++;
//		}
//	    std::cout << " c=" << c << '\n';
	}
}

void query_kd_approximate(float range, float eps)
{
	for (int i = 0; i < position_counter; i++) {
		Point_2 center(position_matrix[i].x, position_matrix[i].y);
		Fuzzy_circle circle(center, range, eps); // define approximate circular range query with given fuzziness

		std::list<Point_2> result;
		cgal_tree->search(std::back_inserter(result), circle);
		range_kd_approximate[i] = result.size() - 1; // skip self
	}
}

void finish_kd()
{
	delete cgal_tree;
	cgal_tree = NULL;
}

/*-------------------------------- NEIGHBORHOOD MATRIX FUNCTIONS --------------------------------*/

void query_nm(neighbor_ptr func, float range)
{
	range = range * range;
	for (int y = 0; y < dim_y; y++) {
		for (int x = 0; x < dim_x; x++) {
			//std::cout << std::setw(4) << position_matrix[i].id << ":nm  " << position_matrix[i].x << ',' << position_matrix[i].y << " -> ";
			int i = x + y * dim_x;
			Position& current = position_matrix[i];
			int *candidate = (*func)(x, y);
			int c = 0, j = 0;
			while (candidate[j] != -1) {
				Position& neighbor = position_matrix[candidate[j]];
				float dx = current.x - neighbor.x;
				float dy = current.y - neighbor.y;
				float sd = dx * dx + dy * dy;
				if (sd <= range) {
					c++;
				}
				j++;
			}
			range_nm[i] = c;
			//std::cout << "nm #" << position_matrix[i].id << " c=" << c << '\n';
		}
	}
}

/*-------------------------------- TEST FUNCTIONS --------------------------------*/

int evaluate_kd_misses()
{
	int count = 0;
	for (int i = 0; i < position_counter; i++) {
		if (range_kd_approximate[i] != range_kd_exact[i]) {
			//std::cout << '#' << pos << " not matched: kd=" << nearest_kd[pos] << " <> nm=" << nearest_nm[i] << '\n';
			count += 1;
		}
	}
	//std::cout << "not matched: " << c << '\n';
	return count;
}

int evaluate_nm_misses()
{
	int count = 0;
	for (int i = 0; i < position_counter; i++) {
		int pos = position_matrix[i].id;
		if (range_nm[i] != range_kd_exact[pos]) {
			//std::cout << '#' << pos << " not matched: kd=" << nearest_kd[pos] << " <> nm=" << nearest_nm[i] << '\n';
			count += 1;
		}
	}
	//std::cout << "not matched: " << c << '\n';
	return count;
}

long int count_all_neighbors(int *matrix)
{
	long int count = 0;
	for (int i = 0; i < position_counter; i++) {
		count += matrix[i];
	}
	return count;
}

void test_knn(const char *dist_name,  function_ptr dist_func,
		      const char *setup_name, function_ptr setup_func,
		      int runs, float range, float eps, const char *type)
{
	double setupkde[runs], querykde[runs], setupkda[runs], querykda[runs], setupnm[runs];
	double missea[runs], totkde[runs], totkda[runs];
	double query4[runs], query8[runs], query24[runs], query48[runs], query80[runs];
	double misse4[runs], misse8[runs], misse24[runs], misse48[runs], misse80[runs];
	double totnm4[runs], totnm8[runs], totnm24[runs], totnm48[runs], totnm80[runs];
	double smean, ssdev, qmean, qsdev, mmean, tmean, tsdev;
	double total_neighbors;
	struct timespec t0, t1;

	neighbor_ptr n4_func, n8_func, n24_func, n48_func, n80_func;
	if (strcmp(type, "hard") == 0) {
		n4_func  = &get_hard_4_neighborhood;
		n8_func  = &get_hard_8_neighborhood;
		n24_func = &get_hard_24_neighborhood;
		n48_func = &get_hard_48_neighborhood;
		n80_func = &get_hard_80_neighborhood;
	}
	else {
		n4_func  = &get_soft_4_neighborhood;
		n8_func  = &get_soft_8_neighborhood;
		n24_func = &get_soft_24_neighborhood;
		n48_func = &get_soft_48_neighborhood;
		n80_func = &get_soft_80_neighborhood;
	}

	for (int r = 1; r <= runs; r++) {
		random_seed(r);
		tcomp = tread = twrit = 0;

		(*dist_func)();

		/*---------------- evaluate exact kd-tree ----------------*/
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t0);
		setup_kd();
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t1);
		setupkde[r - 1] = (t1.tv_sec - t0.tv_sec) * 1000.0 + (t1.tv_nsec - t0.tv_nsec) * 0.000001;

		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t0);
		query_kd_exact(range);
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t1);
		querykde[r - 1] = (t1.tv_sec - t0.tv_sec) * 1000.0 + (t1.tv_nsec - t0.tv_nsec) * 0.000001;
		totkde[r - 1] = count_all_neighbors(range_kd_exact);

		finish_kd();

		/*---------------- evaluate approximate kd-tree ----------------*/
		if (eps != 0) {
			clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t0);
			setup_kd();
			clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t1);
			setupkda[r - 1] = (t1.tv_sec - t0.tv_sec) * 1000.0 + (t1.tv_nsec - t0.tv_nsec) * 0.000001;

			clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t0);
			query_kd_approximate(range, eps);
			clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t1);
			querykda[r - 1] = (t1.tv_sec - t0.tv_sec) * 1000.0 + (t1.tv_nsec - t0.tv_nsec) * 0.000001;
			missea[r - 1] = evaluate_kd_misses();
			totkda[r - 1] = count_all_neighbors(range_kd_approximate);

			finish_kd();
		}

		/*---------------- evaluate neighborhood matrix ----------------*/
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t0);
		(*setup_func)();
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t1);
		setupnm[r - 1] = (t1.tv_sec - t0.tv_sec) * 1000.0 + (t1.tv_nsec - t0.tv_nsec) * 0.000001;

		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t0);
		query_nm(n4_func, range);
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t1);
		query4[r - 1] = (t1.tv_sec - t0.tv_sec) * 1000.0 + (t1.tv_nsec - t0.tv_nsec) * 0.000001;
		misse4[r - 1] = evaluate_nm_misses();
		totnm4[r - 1] = count_all_neighbors(range_nm);

		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t0);
		query_nm(n8_func, range);
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t1);
		query8[r - 1] = (t1.tv_sec - t0.tv_sec) * 1000.0 + (t1.tv_nsec - t0.tv_nsec) * 0.000001;
		misse8[r - 1] = evaluate_nm_misses();
		totnm8[r - 1] = count_all_neighbors(range_nm);

		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t0);
		query_nm(n24_func, range);
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t1);
		query24[r - 1] = (t1.tv_sec - t0.tv_sec) * 1000.0 + (t1.tv_nsec - t0.tv_nsec) * 0.000001;
		misse24[r - 1] = evaluate_nm_misses();
		totnm24[r - 1] = count_all_neighbors(range_nm);

		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t0);
		query_nm(n48_func, range);
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t1);
		query48[r - 1] = (t1.tv_sec - t0.tv_sec) * 1000.0 + (t1.tv_nsec - t0.tv_nsec) * 0.000001;
		misse48[r - 1] = evaluate_nm_misses();
		totnm48[r - 1] = count_all_neighbors(range_nm);

		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t0);
		query_nm(n80_func, range);
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t1);
		query80[r - 1] = (t1.tv_sec - t0.tv_sec) * 1000.0 + (t1.tv_nsec - t0.tv_nsec) * 0.000001;
		misse80[r - 1] = evaluate_nm_misses();
		totnm80[r - 1] = count_all_neighbors(range_nm);
	}

	smean = gsl_stats_mean(setupkde, 1, runs);
	ssdev = gsl_stats_sd(setupkde, 1, runs);
	qmean = gsl_stats_mean(querykde, 1, runs);
	qsdev = gsl_stats_sd(querykde, 1, runs);
	tmean = gsl_stats_mean(totkde, 1, runs);
	tsdev = gsl_stats_sd(totkde, 1, runs);
	std::cout << dist_name  << "  kd-tree   kd-exac  ";
	std::cout << std::setprecision(2) << std::fixed << ' ';
	std::cout << std::setw(10) << smean << "  " << std::setw(10) << ssdev << "  ";
	std::cout << std::setw(10) << qmean << "  " << std::setw(10) << qsdev << "  ";
	std::cout << std::setw(10) << smean + qmean << "           0           0  ";
	std::cout << std::setprecision(4) << std::setw(10) << tmean / position_counter << "  ";
	std::cout << std::setprecision(4) << std::setw(10) << tsdev / position_counter << "            0\n";

	total_neighbors = tmean;

	if (eps != 0) {
		smean = gsl_stats_mean(setupkda, 1, runs);
		ssdev = gsl_stats_sd(setupkda, 1, runs);
		qmean = gsl_stats_mean(querykda, 1, runs);
		qsdev = gsl_stats_sd(querykda, 1, runs);
		mmean = gsl_stats_mean(missea, 1, runs);
		tmean = gsl_stats_mean(totkda, 1, runs);
		tsdev = gsl_stats_sd(totkda, 1, runs);
		std::cout << dist_name  << "  kd-tree   kd-appr  ";
		std::cout << std::setprecision(2) << std::fixed << ' ';
		std::cout << std::setw(10) << smean << "  " << std::setw(10) << ssdev << "  ";
		std::cout << std::setw(10) << qmean << "  " << std::setw(10) << qsdev << "  ";
		std::cout << std::setw(10) << smean + qmean << "  ";
		std::cout << std::setprecision(0) << std::setw(10) << mmean << "  ";
		std::cout << std::setprecision(4) << std::setw(10) << 100 * mmean / position_counter << "  ";
		std::cout << std::setprecision(4) << std::setw(10) << tmean / position_counter << "  ";
		std::cout << std::setprecision(4) << std::setw(10) << tsdev / position_counter << "   ";
		std::cout << std::setprecision(4) << std::setw(10) << 100 * (1 - tmean / total_neighbors) << '\n';
	}

	smean = gsl_stats_mean(setupnm, 1, runs);
	ssdev = gsl_stats_sd(setupnm, 1, runs);

	qmean = gsl_stats_mean(query4, 1, runs);
	qsdev = gsl_stats_sd(query4, 1, runs);
	mmean = gsl_stats_mean(misse4, 1, runs);
	tmean = gsl_stats_mean(totnm4, 1, runs);
	tsdev = gsl_stats_sd(totnm4, 1, runs);
	std::cout << dist_name  << "  " << setup_name << "  " << type << "_4   ";
	std::cout << std::setprecision(2) << std::fixed << ' ';
	std::cout << std::setw(10) << smean << "  " << std::setw(10) << ssdev << "  ";
	std::cout << std::setw(10) << qmean << "  " << std::setw(10) << qsdev << "  ";
	std::cout << std::setw(10) << smean + qmean << "  ";
	std::cout << std::setprecision(0) << std::setw(10) << mmean << "  ";
	std::cout << std::setprecision(4) << std::setw(10) << 100 * mmean / position_counter << "  ";
	std::cout << std::setprecision(4) << std::setw(10) << tmean / position_counter << "  ";
	std::cout << std::setprecision(4) << std::setw(10) << tsdev / position_counter << "   ";
	std::cout << std::setprecision(4) << std::setw(10) << 100 * (1 - tmean / total_neighbors) << "\n";

	qmean = gsl_stats_mean(query8, 1, runs);
	qsdev = gsl_stats_sd(query8, 1, runs);
	mmean = gsl_stats_mean(misse8, 1, runs);
	tmean = gsl_stats_mean(totnm8, 1, runs);
	tsdev = gsl_stats_sd(totnm8, 1, runs);
	std::cout << dist_name  << "  " << setup_name << "  " << type << "_8   ";
	std::cout << std::setprecision(2) << std::fixed << ' ';
	std::cout << std::setw(10) << smean << "  " << std::setw(10) << ssdev << "  ";
	std::cout << std::setw(10) << qmean << "  " << std::setw(10) << qsdev << "  ";
	std::cout << std::setw(10) << smean + qmean << "  ";
	std::cout << std::setprecision(0) << std::setw(10) << mmean << "  ";
	std::cout << std::setprecision(4) << std::setw(10) << 100 * mmean / position_counter << "  ";
	std::cout << std::setprecision(4) << std::setw(10) << tmean / position_counter << "  ";
	std::cout << std::setprecision(4) << std::setw(10) << tsdev / position_counter << "   ";
	std::cout << std::setprecision(4) << std::setw(10) << 100 * (1 - tmean / total_neighbors) << "\n";

	qmean = gsl_stats_mean(query24, 1, runs);
	qsdev = gsl_stats_sd(query24, 1, runs);
	mmean = gsl_stats_mean(misse24, 1, runs);
	tmean = gsl_stats_mean(totnm24, 1, runs);
	tsdev = gsl_stats_sd(totnm24, 1, runs);
	std::cout << dist_name  << "  " << setup_name << "  " << type << "_24  ";
	std::cout << std::setprecision(2) << std::fixed << ' ';
	std::cout << std::setw(10) << smean << "  " << std::setw(10) << ssdev << "  ";
	std::cout << std::setw(10) << qmean << "  " << std::setw(10) << qsdev << "  ";
	std::cout << std::setw(10) << smean + qmean << "  ";
	std::cout << std::setprecision(0) << std::setw(10) << mmean << "  ";
	std::cout << std::setprecision(4) << std::setw(10) << 100 * mmean / position_counter << "  ";
	std::cout << std::setprecision(4) << std::setw(10) << tmean / position_counter << "  ";
	std::cout << std::setprecision(4) << std::setw(10) << tsdev / position_counter << "   ";
	std::cout << std::setprecision(4) << std::setw(10) << 100 * (1 - tmean / total_neighbors) << "\n";

	qmean = gsl_stats_mean(query48, 1, runs);
	qsdev = gsl_stats_sd(query48, 1, runs);
	mmean = gsl_stats_mean(misse48, 1, runs);
	tmean = gsl_stats_mean(totnm48, 1, runs);
	tsdev = gsl_stats_sd(totnm48, 1, runs);
	std::cout << dist_name  << "  " << setup_name << "  " << type << "_48  ";
	std::cout << std::setprecision(2) << std::fixed << ' ';
	std::cout << std::setw(10) << smean << "  " << std::setw(10) << ssdev << "  ";
	std::cout << std::setw(10) << qmean << "  " << std::setw(10) << qsdev << "  ";
	std::cout << std::setw(10) << smean + qmean << "  ";
	std::cout << std::setprecision(0) << std::setw(10) << mmean << "  ";
	std::cout << std::setprecision(4) << std::setw(10) << 100 * mmean / position_counter << "  ";
	std::cout << std::setprecision(4) << std::setw(10) << tmean / position_counter << "  ";
	std::cout << std::setprecision(4) << std::setw(10) << tsdev / position_counter << "   ";
	std::cout << std::setprecision(4) << std::setw(10) << 100 * (1 - tmean / total_neighbors) << "\n";

	qmean = gsl_stats_mean(query80, 1, runs);
	qsdev = gsl_stats_sd(query80, 1, runs);
	mmean = gsl_stats_mean(misse80, 1, runs);
	tmean = gsl_stats_mean(totnm80, 1, runs);
	tsdev = gsl_stats_sd(totnm80, 1, runs);
	std::cout << dist_name  << "  " << setup_name << "  " << type << "_80  ";
	std::cout << std::setprecision(2) << std::fixed << ' ';
	std::cout << std::setw(10) << smean << "  " << std::setw(10) << ssdev << "  ";
	std::cout << std::setw(10) << qmean << "  " << std::setw(10) << qsdev << "  ";
	std::cout << std::setw(10) << smean + qmean << "  ";
	std::cout << std::setprecision(0) << std::setw(10) << mmean << "  ";
	std::cout << std::setprecision(4) << std::setw(10) << 100 * mmean / position_counter << "  ";
	std::cout << std::setprecision(4) << std::setw(10) << tmean / position_counter << "  ";
	std::cout << std::setprecision(4) << std::setw(10) << tsdev / position_counter << "   ";
	std::cout << std::setprecision(4) << std::setw(10) << 100 * (1 - tmean / total_neighbors) << "\n";
}

int main()
{
	int side = 1000;
	int runs = 2;
	float range = 1.0 / side;
	float eps   = 0.2 * range;
	const char *type = "hard";
	//const char *type = "soft";

	verbose = false;
	dim_x = dim_y = side;
	position_matrix = new Position[dim_x * dim_y];
	position_counter = dim_x * dim_y;

    row_is_sorted = new bool[dim_y];
    col_is_sorted = new bool[dim_x];

    range_kd_exact       = new int[dim_x * dim_y];
	range_kd_approximate = new int[dim_x * dim_y];
	range_nm             = new int[dim_x * dim_y];

	std::cout << "\nside=" << side << " runs=" << runs << " range=" << range << " eps=" << eps << "\n\n";
	std::cout << "distrib   setup     query     prepr_mean  prepr_sdev  query_mean  query_sdev  total_mean      miss_#      miss_%   neig_mean   neig_sdev  total_miss%\n";
                //uniform   sp_shell  hard_80        88408        2303       61550         193      149957          94   0.0374400      3.1361      0.0060   0.0124981

	const char *sort_name  = "sp_shell";
	function_ptr sort_func = &full_spatial_shell_sort_skip;

	test_knn("uniform ", &setup_uniform,         sort_name, sort_func, runs, range, eps, type);
//	test_knn("gradient", &setup_gradient,        sort_name, sort_func, runs, range, eps, type);
//	test_knn("gaussian", &setup_gaussian,        sort_name, sort_func, runs, range, eps, type);
//	test_knn("shuffle ", &setup_shuffle,         sort_name, sort_func, runs, range, eps, type);
//	test_knn("inside  ", &setup_inside_circle,   sort_name, sort_func, runs, range, eps, type);
//	test_knn("outside ", &setup_outside_circle,  sort_name, sort_func, runs, range, eps, type);
//	test_knn("clusters", &setup_random_clusters, sort_name, sort_func, runs, range, eps, type);
//	test_knn("holes   ", &setup_random_holes,    sort_name, sort_func, runs, range, eps, type);
//	test_knn("streets ", &setup_streets,         sort_name, sort_func, runs, range, eps, type);
//	test_knn("tilted  ", &setup_tilted,          sort_name, sort_func, runs, range, eps, type);
//	test_knn("triangle", &setup_triangle,        sort_name, sort_func, runs, range, eps, type);

	diagonal_angle = 20;
//	test_knn("diagon20", &setup_diagonal,        sort_name, sort_func, runs, range, eps, type);
//	diagonal_angle = 45;
//	test_knn("diagon45", &setup_diagonal,        sort_name, sort_func, runs, range, eps, type);
//	test_knn("cross   ", &setup_cross,           sort_name, sort_func, runs, range, eps, type);

	return 0;
}
