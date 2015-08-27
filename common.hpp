#ifndef COMMON_HPP
#define COMMON_HPP

/*-------------------------------- MACRO DEFINITIONS --------------------------------*/

#define X_COMPARE_BOTH(A,B) ((A.x > B.x) || (A.x == B.x && A.y > B.y))
#define Y_COMPARE_BOTH(A,B) ((A.y > B.y) || (A.y == B.y && A.x > B.x))

#define X_COMPARE_ONLY(A,B) (A.x > B.x)
#define Y_COMPARE_ONLY(A,B) (A.y > B.y)
#define Z_COMPARE_ONLY(A,B) (A.z > B.z)

#define X_COMPARE X_COMPARE_ONLY
#define Y_COMPARE Y_COMPARE_ONLY
#define Z_COMPARE Z_COMPARE_ONLY

/*-------------------------------- TYPE DEFINITIONS --------------------------------*/

struct Position {
	float x, y;
	int id;
};

struct Position3d {
	float x, y, z;
	int id;
};

typedef void (*function_ptr)();
typedef int *(*neighbor_ptr)(int, int);

/*-------------------------------- ALGORITHMS MODULE FUNCTIONS --------------------------------*/

bool are_all_rows_sorted();
bool are_all_columns_sorted();
bool is_spatially_sorted();

bool odd_even_partial_sort();
void odd_even_full_sort();

void linear_selection_sort_on_each_row();
void linear_selection_sort_on_each_col();
void full_selection();

int  linear_insertion_sort_on_each_row();
int  linear_insertion_sort_on_each_col();
int  linear_insertion_sort_on_each_row_3d();
int  linear_insertion_sort_on_each_col_3d();
int  linear_insertion_sort_on_each_sta_3d();
int  linear_insertion_sort_on_each_row_skip();
int  linear_insertion_sort_on_each_col_skip();
void full_insertion();
void full_insertion_3d();
void full_insertion_skip();

void linear_shell_sort_on_each_row();
void linear_shell_sort_on_each_col();
void spatial_shell_sort();
void full_spatial_shell_sort();
void full_spatial_shell_sort_skip();
void full_shell_odd_even();

void linear_quicksort_on_each_row();
void linear_quicksort_on_each_col();
void partial_quicksort();
void full_quicksort();
void full_quick_odd_even();
void full_quick_insertion();
void full_quick_insertion_skip();
void interleaved_quicksort();

void simple_sort();

int *get_hard_4_neighborhood(int x, int y);
int *get_hard_8_neighborhood(int x, int y);
int *get_hard_24_neighborhood(int x, int y);
int *get_hard_48_neighborhood(int x, int y);
int *get_hard_80_neighborhood(int x, int y);
int *get_hard_120_neighborhood(int x, int y);

void get_hard_3d_neighborhood(int x, int y, int z, int n_siz, int *candidate);

int *get_soft_4_neighborhood(int x, int y);
int *get_soft_8_neighborhood(int x, int y);
int *get_soft_24_neighborhood(int x, int y);
int *get_soft_48_neighborhood(int x, int y);
int *get_soft_80_neighborhood(int x, int y);
int *get_soft_120_neighborhood(int x, int y);

float get_distance(Position& a, Position& b);
float get_squared_distance(Position& a, Position& b);

/*-------------------------------- DISTRIBUTIONS MODULE FUNCTIONS --------------------------------*/

void  random_init();
void  random_seed(int seed);
void  random_done();
float random_range(float min, float max);

void setup_uniform();
void setup_uniform_3d();
void setup_shuffle();
void setup_gaussian();
void setup_inside_circle();
void setup_outside_circle();
void setup_random_clusters();
void setup_random_holes();
void setup_streets();
void setup_tilted();

void setup_diagonal();
void setup_cross();
void setup_triangle();
void setup_gradient();
void setup_compact();

#endif // COMMON_HPP

