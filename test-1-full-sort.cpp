/*-------------------------------- INCLUDES --------------------------------*/

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <iomanip>

#include <time.h>

#include <gsl/gsl_statistics.h>

#include <openssl/md5.h>

#include "common.hpp"

/*-------------------------------- IMPORTED VARIABLES --------------------------------*/

extern int dim_x, dim_y;
extern Position *position_matrix;
extern int position_counter;

extern bool *row_is_sorted;
extern bool *col_is_sorted;

extern int verbose;
extern double tcomp, tread, twrit, tpass;

extern function_ptr prev_pass_func;
extern function_ptr post_pass_func;

extern float diagonal_angle;
//extern float diagonal_width;

/*-------------------------------- LOCAL VARIABLES --------------------------------*/

double *comp, *msec, *pass;
int *last_row_for_id;
int *last_col_for_id;

/*------------------------ FULL SORT FUNCTIONS ------------------------*/

void save_row_and_col_per_id()
{
	for (int row = 0; row < dim_y; row++) {
		for (int col = 0; col < dim_x; col++) {
			int index = col + row * dim_x;
			int id = position_matrix[index].id;
			last_row_for_id[id] = row;
			last_col_for_id[id] = col;
		}
	}
}

void calculate_total_travel()
{
	double euclidean = 0;
	double max_dx, max_dy, sum_dx, sum_dy;
	max_dx = max_dy = sum_dx = sum_dy = 0;
	for (int row = 0; row < dim_y; row++) {
		for (int col = 0; col < dim_x; col++) {
			int index = col + row * dim_x;
			int id = position_matrix[index].id;
			int old_row = last_row_for_id[id];
			int old_col = last_col_for_id[id];
			euclidean += sqrt((double) (row - old_row) * (row - old_row) + (col - old_col) * (col - old_col));
			if (abs(col - old_col) > max_dx) {
				max_dx = abs(col - old_col);
			}
			if (abs(row - old_row) > max_dy) {
				max_dy = abs(row - old_row);
			}
			sum_dx += abs(col - old_col);
			sum_dy += abs(row - old_row);
		}
	}
	std::cout << std::setprecision(0) << std::fixed;
	std::cout << "     euc=" << euclidean << "  man=" << sum_dx + sum_dy;
	std::cout << "  sum_dx=" << sum_dx << "  sum_dy=" << sum_dy << "  max_dx=" << max_dx << "  max_dy=" << max_dy <<'\n';
}

int get_hash()
{
	unsigned char *data = (unsigned char *) position_matrix;
	unsigned long  size = sizeof(Position) * position_counter;
	unsigned char *hash = MD5(data, size, NULL);
	return (hash[0] << 8) | hash[1];
}

void update(float percent, float offset)
{
	int n = (int) (position_counter * percent / 100);
	for (int j = 0; j < n; j++) {
		int i = (int) random_range(0, position_counter - 1);
		float x = position_matrix[i].x + random_range(- offset, offset);
		float y = position_matrix[i].y + random_range(- offset, offset);
		position_matrix[i].x = (x <= 0) ? 0 : (x <= 1) ? x : 1;
		position_matrix[i].y = (y <= 0) ? 0 : (y <= 1) ? y : 1;
	}
}

void single_sort(const char *dist_name, function_ptr dist_func, const char *sort_name, function_ptr sort_func, int runs /*, float percent, float offset*/)
{
	double cmean, csdev, tmean, tsdev, pmean, psdev;
	struct timespec t0, t1;
	int hash = 0;
	//std::cout << "percent=" << percent << "  offset=" << offset << "  ";

	int progress = (runs <= 10) ? 1 : runs / 10;
	std::cout << std::setprecision(0) << std::fixed;

	if (! verbose) {
		std::cout << dist_name << "  " << sort_name << "  " << std::flush;
	}

	for (int i = 1; i <= runs; i++) {
		random_seed(i);
		tcomp = tread = twrit = tpass = 0;
		(*dist_func)();

		if (verbose) {
			if (verbose == 1) {
				std::cout << "seed=" << i << "  ";
			}
			else {
				std::cout << "---------------- seed=" << i << '\n';
			}
			(*sort_func)();
			msec[i - 1] = 0;
		}
		else {
			if ((i % progress) == 0) {
				std::cout << '.' << std::flush;
			}
			//(*sort_func)();
			//update(percent, offset);
			clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t0);
			(*sort_func)();
			clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t1);
			msec[i - 1] = (t1.tv_sec - t0.tv_sec) * 1000.0 + (t1.tv_nsec - t0.tv_nsec) * 0.000001;
		}

		comp[i - 1] = tcomp;
		pass[i - 1] = tpass;
		if (i == 1) {
			hash = get_hash();
		}
		else {
			hash ^= get_hash();
		}
		if (! is_spatially_sorted()) {
			std::cout << "result is not sorted!\n";
		}
	}

	if (verbose) {
		std::cout << dist_name << "  " << sort_name << "            ";
	}
	else {
		if (progress == 1) {
			for (int i = runs; i < 10; i++) {
				std::cout << ' ';
			}
		}
	}
	cmean = gsl_stats_mean(comp, 1, runs);
	csdev = gsl_stats_sd(comp, 1, runs);
	tmean = gsl_stats_mean(msec, 1, runs);
	tsdev = gsl_stats_sd(msec, 1, runs);
	pmean = gsl_stats_mean(pass, 1, runs);
	psdev = gsl_stats_sd(pass, 1, runs);
	std::cout << std::setw(12) << cmean << ' ' << std::setw(12) << csdev << ' ';
	std::cout << std::setprecision(2) << std::fixed;
	std::cout << std::setw(12) << tmean << ' ' << std::setw(12) << tsdev << ' ';
	std::cout << std::setw(12) << pmean << ' ' << std::setw(12) << psdev << "  ";
	std::cout << std::hex << hash << std::dec << '\n';
}

void all_sorts(const char *dist_name, function_ptr dist_func, int runs)
{
//	single_sort(dist_name, dist_func, "rep_insertion_skp", &full_insertion_skip, runs,          1, 0.2);
//	single_sort(dist_name, dist_func, "rep_insertion_skp", &full_insertion_skip, runs,          2, 0.2);
//	single_sort(dist_name, dist_func, "rep_insertion_skp", &full_insertion_skip, runs,         10, 0.2);
//	single_sort(dist_name, dist_func, "rep_insertion_skp", &full_insertion_skip, runs,         25, 0.2);
//	single_sort(dist_name, dist_func, "rep_insertion_skp", &full_insertion_skip, runs,         50, 0.2);
//	single_sort(dist_name, dist_func, "sp_shell+in_skp  ", &full_spatial_shell_sort_skip, runs,  1, 0.2);
//	single_sort(dist_name, dist_func, "sp_shell+in_skp  ", &full_spatial_shell_sort_skip, runs,  2, 0.2);
//	single_sort(dist_name, dist_func, "sp_shell+in_skp  ", &full_spatial_shell_sort_skip, runs, 10, 0.2);
//	single_sort(dist_name, dist_func, "sp_shell+in_skp  ", &full_spatial_shell_sort_skip, runs, 25, 0.2);
//	single_sort(dist_name, dist_func, "sp_shell+in_skp  ", &full_spatial_shell_sort_skip, runs, 50, 0.2);

	//single_sort(dist_name, dist_func, "simple           ", &simple_sort, runs, 10, 0.2);
	//single_sort(dist_name, dist_func, "rep_odd_even     ", &odd_even_full_sort, runs);
	//single_sort(dist_name, dist_func, "rep_insertion    ", &full_insertion, runs);
	//single_sort(dist_name, dist_func, "rep_insertion_skp", &full_insertion_skip, runs);
	//single_sort(dist_name, dist_func, "rep_quick_sort   ", &full_quicksort, runs);
	//single_sort(dist_name, dist_func, "2x_quick+oe      ", &full_quick_odd_even, runs);
	//single_sort(dist_name, dist_func, "2x_quick+in      ", &full_quick_insertion, runs);
	single_sort(dist_name, dist_func, "2x_quick+in_skp  ", &full_quick_insertion_skip, runs);
	//single_sort(dist_name, dist_func, "sp_shell+oe      ", &full_shell_odd_even, runs);
	//single_sort(dist_name, dist_func, "sp_shell+in      ", &full_spatial_shell_sort, runs);
	single_sort(dist_name, dist_func, "sp_shell+in_skp  ", &full_spatial_shell_sort_skip, runs);
}

int main()
{
	int side = 1000;
	int runs = 10;
	verbose = 0;

	comp = new double[runs];
	msec = new double[runs];
	pass = new double[runs];
	dim_x = dim_y = side;
	position_matrix = new Position[dim_x * dim_y];
	position_counter = dim_x * dim_y;

	row_is_sorted = new bool[dim_y];
	col_is_sorted = new bool[dim_x];

	if (verbose == 2) {
		last_row_for_id = new int[dim_x * dim_y];
		last_col_for_id = new int[dim_x * dim_y];
		prev_pass_func = save_row_and_col_per_id;
		post_pass_func = calculate_total_travel;
	}

	std::cout << "\nside=" << side << " runs=" << runs << "\n\n";
	std::cout << "distrib   sort                            mean_comp    sdev_comp    mean_msec    sdev_msec    mean_pass    sdev_pass  hash\n";
	//            uniform   rep_odd_even                     17625600       184640            0            0          135            1  128c

	all_sorts("uniform ", &setup_uniform, runs);

	all_sorts("gradient", &setup_gradient, runs);
	all_sorts("gaussian", &setup_gaussian, runs);
	all_sorts("shuffle ", &setup_shuffle, runs);
//	all_sorts("inside  ", &setup_inside_circle, runs);
//	all_sorts("outside ", &setup_outside_circle, runs);
//	all_sorts("clusters", &setup_random_clusters, runs);
//	all_sorts("holes   ", &setup_random_holes, runs);
//	all_sorts("streets ", &setup_streets, runs);
//	all_sorts("tilted  ", &setup_tilted, runs);
//	all_sorts("triangle", &setup_triangle, runs);

	diagonal_angle = 20;
//	all_sorts("diagon20", &setup_diagonal, runs);
//	diagonal_angle = 45;
//	all_sorts("diagon45", &setup_diagonal, runs);
//	all_sorts("cross   ", &setup_cross, runs);

	return 0;
}

