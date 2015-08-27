/*-------------------------------- INCLUDES --------------------------------*/

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>

#include <time.h>

#include <gsl/gsl_statistics.h>

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

/*-------------------------------- LOCAL VARIABLES --------------------------------*/

double *comp, *msec;

/*------------------------ FULL SORT FUNCTIONS ------------------------*/

typedef void (*function_ptr)();

void single_sort(const char *dist_name, function_ptr dist_func, const char *sort_name, function_ptr sort_func, int runs, int side)
{
	double cmean, csdev, tmean, tsdev;
	struct timespec t0, t1;

	int progress = (runs <= 10) ? 1 : runs / 10;

	std::cout << std::setw(4) << side << "  ";
	std::cout << dist_name << "  " << sort_name << "  " << std::flush;
	for (int i = 1; i <= runs; i++) {
		random_seed(i);
		tcomp = tread = twrit = 0;
		(*dist_func)();

		if ((i % progress) == 0) {
			std::cout << '.' << std::flush;
		}
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t0);
		(*sort_func)();
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t1);

		comp[i - 1] = tcomp;
		msec[i - 1] = (t1.tv_sec - t0.tv_sec) * 1000.0 + (t1.tv_nsec - t0.tv_nsec) * 0.000001;
		if (! is_spatially_sorted()) {
			std::cout << "result is not sorted!\n";
		}
	}
	cmean = gsl_stats_mean(comp, 1, runs);
	csdev = gsl_stats_sd(comp, 1, runs);
	tmean = gsl_stats_mean(msec, 1, runs);
	tsdev = gsl_stats_sd(msec, 1, runs);
	std::cout << std::setprecision(0) << std::fixed;
	std::cout << std::setw(12) << cmean << ' ' << std::setw(12) << csdev << ' ';
	std::cout << std::setprecision(2) << std::fixed;
	std::cout << std::setw(12) << tmean << ' ' << std::setw(12) << tsdev << '\n';
}

void all_sorts(const char *dist_name, function_ptr dist_func, int runs, int side)
{
	//single_sort(dist_name, dist_func, "rep_odd_even     ", &odd_even_full_sort, runs, side);
	//single_sort(dist_name, dist_func, "2x_quick+oe      ", &full_quick_odd_even, runs, side);
	//single_sort(dist_name, dist_func, "2x_quick+in      ", &full_quick_insertion, runs, side);
	//single_sort(dist_name, dist_func, "2x_quick+in_sk   ", &full_quick_insertion_skip, runs, side);
	//single_sort(dist_name, dist_func, "sp_shell+oe      ", &full_shell_odd_even, runs, side);
	//single_sort(dist_name, dist_func, "sp_shell+in      ", &full_spatial_shell_sort, runs, side);
	single_sort(dist_name, dist_func, "sp_shell+in_sk   ", &full_spatial_shell_sort_skip, runs, side);
}

extern float diagonal_angle;
//extern float diagonal_width;

int main()
{
	int runs = 10;

	comp = new double[runs];
	msec = new double[runs];
	verbose = false;

	std::cout << "\nruns=" << runs << "\n\n";
	std::cout << "side  distrib   sort                            mean_comp    sdev_comp    mean_msec    sdev_msec\n";

	for (int side = 200; side <= 1200; side += 20) {
		dim_x = dim_y = side;
		position_matrix = new Position[dim_x * dim_y];
		position_counter = dim_x * dim_y;

                row_is_sorted = new bool[dim_x];
                col_is_sorted = new bool[dim_y];

		all_sorts("uniform ", &setup_uniform, runs, side);
//		all_sorts("gaussian", &setup_gaussian, runs, side);
//		all_sorts("shuffle ", &setup_shuffle, runs, side);
//		all_sorts("inside  ", &setup_inside_circle, runs, side);
//		all_sorts("outside ", &setup_outside_circle, runs, side);
//		all_sorts("clusters", &setup_random_clusters, runs, side);
//		all_sorts("holes   ", &setup_random_holes, runs, side);
//		all_sorts("streets ", &setup_streets, runs, side);
//		all_sorts("tilted  ", &setup_tilted, runs, side);

//		diagonal_angle = 20;
//		all_sorts("diagon20", &setup_diagonal, runs, side);
//		diagonal_angle = 45;
//		all_sorts("diagon45", &setup_diagonal, runs, side);

//		all_sorts("cross   ", &setup_cross, runs, side);
//		all_sorts("triangle", &setup_triangle, runs, side);

		delete[] position_matrix;
                delete[] row_is_sorted;
                delete[] col_is_sorted;
	}
	return 0;
}
