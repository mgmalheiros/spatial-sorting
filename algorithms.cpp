/*-------------------------------- INCLUDES --------------------------------*/

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <iomanip>

#include "common.hpp"

/*-------------------------------- EXPORTED VARIABLES --------------------------------*/

int dim_x, dim_y, dim_z;
Position *position_matrix;
Position3d *position_matrix_3d;
int position_counter;

bool *row_is_sorted;
bool *col_is_sorted;

int verbose = 0;
double tcomp, tread, twrit, tpass;

function_ptr prev_pass_func = NULL;
function_ptr post_pass_func = NULL;

/*-------------------------------- LOCAL VARIABLES --------------------------------*/

long int qcomp, qread, qwrit;

int neighborhood_4[5];
int neighborhood_8[9];
int neighborhood_24[25];
int neighborhood_48[49];
int neighborhood_80[81];
int neighborhood_120[121];

/*------------------------ SPATIAL SORT VALIDATION ------------------------*/

bool are_all_rows_sorted()
{
    // validate each row
    for (int row = 0; row < dim_y; row++) {
    	for (int col = 0; col < dim_x - 1; col++) {
        	int curr = col + row * dim_x;
            int neig = curr + 1;
            Position& current  = position_matrix[curr];
            Position& neighbor = position_matrix[neig];
            if (X_COMPARE(current, neighbor)) {
                return false;
            }
        }
    }
    return true;
}

bool are_all_columns_sorted()
{
    // validate each column
    for (int col = 0; col < dim_x; col++) {
    	for (int row = 0; row < dim_y - 1; row++) {
            int curr = col + row * dim_x;
            int neig = curr + dim_x;
            Position& current  = position_matrix[curr];
            Position& neighbor = position_matrix[neig];
            if (Y_COMPARE(current, neighbor)) {
            	return false;
            }
        }
    }
    return true;
}

bool is_spatially_sorted()
{
	return are_all_rows_sorted() && are_all_columns_sorted();
}

/*------------------------ ODD EVEN SORT ------------------------*/

bool odd_even_partial_sort()
{
	long int comp = 0, read = 0, writ = 0;

    // compare starting at an odd column -- somehow, column-wise seems to be faster
    for (int col = 1; col < dim_x - 1; col += 2) {
        for (int row = 0; row < dim_y; row++) {
        	int curr = col + row * dim_x;
            int neig = curr + 1;
            Position& current  = position_matrix[curr]; read++;
            Position& neighbor = position_matrix[neig]; read++;
            comp++;
            if (X_COMPARE(current, neighbor)) {
                Position temp = current;
                position_matrix[curr] = neighbor; writ++;
                position_matrix[neig] = temp;     writ++;
            }
        }
    }
    tpass++;

    // compare starting at an even column
    for (int col = 0; col < dim_x - 1; col += 2) {
        for (int row = 0; row < dim_y; row++) {
            int curr = col + row * dim_x;
            int neig = curr + 1;
            Position& current  = position_matrix[curr]; read++;
            Position& neighbor = position_matrix[neig]; read++;
            comp++;
            if (X_COMPARE(current, neighbor)) {
                Position temp = current;
                position_matrix[curr] = neighbor; writ++;
                position_matrix[neig] = temp;     writ++;
            }
        }
    }
    tpass++;

    // compare starting at an odd row
    for (int row = 1; row < dim_y - 1; row += 2) {
    	for (int col = 0; col < dim_x; col++) {
    		int curr = col + row * dim_x;
    		int neig = curr + dim_x;
    		Position& current  = position_matrix[curr]; read++;
    		Position& neighbor = position_matrix[neig]; read++;
    		comp++;
            if (Y_COMPARE(current, neighbor)) {
    			Position temp = current;
    			position_matrix[curr] = neighbor; writ++;
    			position_matrix[neig] = temp;     writ++;
    		}
    	}
    }
    tpass++;

    // compare starting at an even row
    for (int row = 0; row < dim_y - 1; row += 2) {
    	for (int col = 0; col < dim_x; col++) {
    		int curr = col + row * dim_x;
    		int neig = curr + dim_x;
    		Position& current  = position_matrix[curr]; read++;
    		Position& neighbor = position_matrix[neig]; read++;
    		comp++;
            if (Y_COMPARE(current, neighbor)) {
    			Position temp = current;
    			position_matrix[curr] = neighbor; writ++;
    			position_matrix[neig] = temp;     writ++;
    		}
    	}
    }
    tpass++;

    if (verbose >= 2) {
    	std::cout << "oep  comp=" << comp << " read=" << read << " writ=" << writ << '\n';
    }
    tcomp += comp; tread += read; twrit += writ;
    return writ == 0;
}

void odd_even_full_sort()
{
    while (odd_even_partial_sort() == false) {
    }
    if (verbose >= 1) {
    	std::cout << "oef  tcomp=" << tcomp << " read=" << tread << " twrit=" << twrit << " tpass=" << tpass << '\n';
    }
}

/*------------------------ SELECTION SORT ------------------------*/

void linear_selection_sort_on_each_row()
{
	long int comp = 0, read = 0, writ = 0;
	for (int row = 0; row < dim_y; row++) {
		for (int j = 0; j < dim_x - 1; j++) {
			int min = j + row * dim_x;
			for (int i = j + 1; i < dim_x; i++) {
				int cur = i + row * dim_x;
				Position& current = position_matrix[cur]; read++;
				Position& minimum = position_matrix[min]; read++;
				comp++;
				if (X_COMPARE(minimum, current)) {
					min = cur;
				}
			}
			int pos = j + row * dim_x;
			if (min != pos) {
				Position temp = position_matrix[pos];
				position_matrix[pos] = position_matrix[min];
				position_matrix[min] = temp;
				read += 2;
				writ += 2;
			}
		}
	}
    if (verbose >= 2) {
    	std::cout << "lsr  comp=" << comp << " read=" << read << " writ=" << writ << '\n';
    }
    tcomp += comp; tread += read; twrit += writ; tpass++;
}

void linear_selection_sort_on_each_col()
{
	long int comp = 0, read = 0, writ = 0;
	for (int col = 0; col < dim_x; col++) {
		for (int j = 0; j < dim_y - 1; j++) {
			int min = col + j * dim_x;
			for (int i = j + 1; i < dim_y; i++) {
				int cur = col + i * dim_x;
				Position& current = position_matrix[cur]; read++;
				Position& minimum = position_matrix[min]; read++;
				comp++;
				if (Y_COMPARE(minimum, current)) {
					min = cur;
				}
			}
			int pos = col + j * dim_x;
			if (min != pos) {
				Position temp = position_matrix[pos];
				position_matrix[pos] = position_matrix[min];
				position_matrix[min] = temp;
				read += 2;
				writ += 2;
			}
		}
	}
    if (verbose >= 2) {
    	std::cout << "lsc  comp=" << comp << " read=" << read << " writ=" << writ << '\n';
    }
    tcomp += comp; tread += read; twrit += writ; tpass++;
}

void full_selection()
{
	double old_twrit = 0;
	while (true) {
		linear_selection_sort_on_each_row();
		linear_selection_sort_on_each_col();
		if (old_twrit == twrit) {
			break;
		}
		old_twrit = twrit;
	}
    if (verbose >= 1) {
    	std::cout << "fse  tcomp=" << tcomp << " read=" << tread << " twrit=" << twrit << " tpass=" << tpass << '\n';
    }
}

// the X full linear sort should calculate matrix positions differently (like column-major in memory)
//void full_linear_selection_sort_on_y()
//{
//	for (int j = 0; j < position_counter - 1; j++) {
//		int min = j;
//		for (int i = j + 1; i < position_counter; i++) {
//			Position& current = position_matrix[i];
//			Position& minimum = position_matrix[min];
//			if (Y_COMPARE(minimum, current)) {
//				min = i;
//			}
//		}
//		if (min != j) {
//			Position temp = position_matrix[j];
//			position_matrix[j] = position_matrix[min];
//			position_matrix[min] = temp;
//		}
//	}
//}

/*------------------------ INSERTION SORT ------------------------*/

int linear_insertion_sort_on_each_row()
{
	long int comp = 0, read = 0, writ = 0;
	for (int row = 0; row < dim_y; row++) {
		for (int i = 1; i < dim_x; i++) {
			Position temp = position_matrix[i + row * dim_x]; read++;
			int j;
			for (j = i; j >= 1; j--) {
				Position curr = position_matrix[(j - 1) + row * dim_x]; read++;
				comp++;
				if (X_COMPARE(curr, temp)) {
					position_matrix[j + row * dim_x] = curr; writ++;
				}
				else {
					break;
				}
			}
			if (j != i) {
				position_matrix[j + row * dim_x] = temp; writ++;
			}
		}
	}
    if (verbose >= 2) {
    	std::cout << "lir  comp=" << comp << " read=" << read << " writ=" << writ << '\n';
    }
    tcomp += comp; tread += read; twrit += writ; tpass++;
    return writ;
}

int linear_insertion_sort_on_each_col()
{
	long int comp = 0, read = 0, writ = 0;
	for (int col = 0; col < dim_x; col++) {
		for (int i = 1; i < dim_y; i++) {
			Position temp = position_matrix[col + i * dim_x]; read++;
			int j;
			for (j = i; j >= 1; j--) {
				Position curr = position_matrix[col + (j - 1) * dim_x]; read++;
				comp++;
				if (Y_COMPARE(curr, temp)) {
					position_matrix[col + j * dim_x] = curr; writ++;
				}
				else {
					break;
				}
			}
			if (j != i) {
				position_matrix[col + j * dim_x] = temp; writ++;
			}
		}
	}
    if (verbose >= 2) {
    	std::cout << "lic  comp=" << comp << " read=" << read << " writ=" << writ << '\n';
    }
    tcomp += comp; tread += read; twrit += writ; tpass++;
    return writ;
}

int linear_insertion_sort_on_each_row_3d()
{
	long int comp = 0, read = 0, writ = 0;
	for (int layer = 0; layer < dim_z; layer++) {
		for (int row = 0; row < dim_y; row++) {
			for (int i = 1; i < dim_x; i++) {
				Position3d temp = position_matrix_3d[i + row * dim_x + layer * dim_x * dim_y]; read++;
				int j;
				for (j = i; j >= 1; j--) {
					Position3d curr = position_matrix_3d[(j - 1) + row * dim_x + layer * dim_x * dim_y]; read++;
					comp++;
					if (X_COMPARE(curr, temp)) {
						position_matrix_3d[j + row * dim_x + layer * dim_x * dim_y] = curr; writ++;
					}
					else {
						break;
					}
				}
				if (j != i) {
					position_matrix_3d[j + row * dim_x + layer * dim_x * dim_y] = temp; writ++;
				}
			}
		}
	}
    if (verbose >= 2) {
    	std::cout << "lir3d  comp=" << comp << " read=" << read << " writ=" << writ << '\n';
    }
    tcomp += comp; tread += read; twrit += writ; tpass++;
    return writ;
}

int linear_insertion_sort_on_each_col_3d()
{
	long int comp = 0, read = 0, writ = 0;
	for (int layer = 0; layer < dim_z; layer++) {
		for (int col = 0; col < dim_x; col++) {
			for (int i = 1; i < dim_y; i++) {
				Position3d temp = position_matrix_3d[col + i * dim_x + layer * dim_x * dim_y]; read++;
				int j;
				for (j = i; j >= 1; j--) {
					Position3d curr = position_matrix_3d[col + (j - 1) * dim_x + layer * dim_x * dim_y]; read++;
					comp++;
					if (Y_COMPARE(curr, temp)) {
						position_matrix_3d[col + j * dim_x + layer * dim_x * dim_y] = curr; writ++;
					}
					else {
						break;
					}
				}
				if (j != i) {
					position_matrix_3d[col + j * dim_x + layer * dim_x * dim_y] = temp; writ++;
				}
			}
		}
	}
    if (verbose >= 2) {
    	std::cout << "lic3d  comp=" << comp << " read=" << read << " writ=" << writ << '\n';
    }
    tcomp += comp; tread += read; twrit += writ; tpass++;
    return writ;
}

int linear_insertion_sort_on_each_sta_3d()
{
	long int comp = 0, read = 0, writ = 0;
	for (int row = 0; row < dim_y; row++) {
		for (int col = 0; col < dim_x; col++) {
			for (int i = 1; i < dim_z; i++) {
				Position3d temp = position_matrix_3d[col + row * dim_x + i * dim_x * dim_y]; read++;
				int j;
				for (j = i; j >= 1; j--) {
					Position3d curr = position_matrix_3d[col + row * dim_x + (j - 1) * dim_x * dim_y]; read++;
					comp++;
					if (Z_COMPARE(curr, temp)) {
						position_matrix_3d[col + row * dim_x + j * dim_x * dim_y] = curr; writ++;
					}
					else {
						break;
					}
				}
				if (j != i) {
					position_matrix_3d[col + row * dim_x + j * dim_x * dim_y] = temp; writ++;
				}
			}
		}
	}
    if (verbose >= 2) {
    	std::cout << "lis3d  comp=" << comp << " read=" << read << " writ=" << writ << '\n';
    }
    tcomp += comp; tread += read; twrit += writ; tpass++;
    return writ;
}

int linear_insertion_sort_on_each_row_skip()
{
	long int comp = 0, read = 0, writ = 0;
	for (int row = 0; row < dim_y; row++) {
		if (row_is_sorted[row]) {
			continue;
		}
		for (int i = 1; i < dim_x; i++) {
			Position temp = position_matrix[i + row * dim_x]; read++;
			int j;
			for (j = i; j >= 1; j--) {
				Position curr = position_matrix[(j - 1) + row * dim_x]; read++;
				comp++;
				if (X_COMPARE(curr, temp)) {
					position_matrix[j + row * dim_x] = curr; writ++;
					col_is_sorted[j] = false;
				}
				else {
					break;
				}
			}
			if (j != i) {
				position_matrix[j + row * dim_x] = temp; writ++;
				col_is_sorted[j] = false;
			}
		}
		row_is_sorted[row] = true;
	}
    if (verbose >= 2) {
    	std::cout << "lirsk comp=" << comp << " read=" << read << " writ=" << writ << '\n';
    }
    tcomp += comp; tread += read; twrit += writ; tpass++;
    return writ;
}

int linear_insertion_sort_on_each_col_skip()
{
	long int comp = 0, read = 0, writ = 0;
	for (int col = 0; col < dim_x; col++) {
		if (col_is_sorted[col]) {
			continue;
		}
		for (int i = 1; i < dim_y; i++) {
			Position temp = position_matrix[col + i * dim_x]; read++;
			int j;
			for (j = i; j >= 1; j--) {
				Position curr = position_matrix[col + (j - 1) * dim_x]; read++;
				comp++;
				if (Y_COMPARE(curr, temp)) {
					position_matrix[col + j * dim_x] = curr; writ++;
					row_is_sorted[j] = false;
				}
				else {
					break;
				}
			}
			if (j != i) {
				position_matrix[col + j * dim_x] = temp; writ++;
				row_is_sorted[j] = false;
			}
		}
		col_is_sorted[col] = true;
	}
    if (verbose >= 2) {
    	std::cout << "licsk  comp=" << comp << " read=" << read << " writ=" << writ << '\n';
    }
    tcomp += comp; tread += read; twrit += writ; tpass++;
    return writ;
}

void full_insertion()
{
	int writes = 1;
	while (writes != 0) {
		if (prev_pass_func) { (*prev_pass_func)(); }
		writes = linear_insertion_sort_on_each_row();
		if (post_pass_func) { (*post_pass_func)(); }

		if (writes == 0) {
			break;
		}

		if (prev_pass_func) { (*prev_pass_func)(); }
		writes = linear_insertion_sort_on_each_col();
		if (post_pass_func) { (*post_pass_func)(); }
	}
    if (verbose >= 1) {
    	std::cout << "fin  tcomp=" << tcomp << " read=" << tread << " twrit=" << twrit << " tpass=" << tpass << '\n';
    }
}

void full_insertion_3d()
{
	int writes = 1;
	while (writes != 0) {
		if (prev_pass_func) { (*prev_pass_func)(); }
		writes = linear_insertion_sort_on_each_row_3d();
		if (post_pass_func) { (*post_pass_func)(); }

		if (writes == 0) {
			break;
		}

		if (prev_pass_func) { (*prev_pass_func)(); }
		writes = linear_insertion_sort_on_each_col_3d();
		if (post_pass_func) { (*post_pass_func)(); }

		if (writes == 0) {
			break;
		}

		if (prev_pass_func) { (*prev_pass_func)(); }
		writes = linear_insertion_sort_on_each_sta_3d();
		if (post_pass_func) { (*post_pass_func)(); }
	}
    if (verbose >= 1) {
    	std::cout << "fin3d  tcomp=" << tcomp << " read=" << tread << " twrit=" << twrit << " tpass=" << tpass << '\n';
    }
}

void full_insertion_skip()
{
	int writes = 1;

	memset(row_is_sorted, 0, sizeof(bool) * dim_y);
	memset(col_is_sorted, 0, sizeof(bool) * dim_x);

	linear_insertion_sort_on_each_row();
	writes = linear_insertion_sort_on_each_col();

	while (writes != 0) {
		if (prev_pass_func) { (*prev_pass_func)(); }
		writes = linear_insertion_sort_on_each_row_skip();
		if (post_pass_func) { (*post_pass_func)(); }

		if (writes == 0) {
			break;
		}

		if (prev_pass_func) { (*prev_pass_func)(); }
		writes = linear_insertion_sort_on_each_col_skip();
		if (post_pass_func) { (*post_pass_func)(); }
	}
    if (verbose >= 1) {
    	std::cout << "finsk  tcomp=" << tcomp << " read=" << tread << " twrit=" << twrit << " tpass=" << tpass << '\n';
    }
}

/*------------------------ SHELL SORT ------------------------*/

void linear_shell_sort_on_each_row()
{
	int gaps[] = {1750, 701, 301, 132, 57, 23, 10, 4, 1};
	long int comp = 0, read = 0, writ = 0;

	for (int g = 0; g < (int) (sizeof(gaps) / sizeof(int)); g++) {
		int gap = gaps[g];
		if (gap >= dim_x) {
			continue;
		}

		for (int row = 0; row < dim_y; row++) {
			for (int col = gap; col < dim_x; col++) {
				Position temp = position_matrix[col + row * dim_x]; read++;
				int j;
				for (j = col; j >= gap; j -= gap) {
					Position curr = position_matrix[(j - gap) + row * dim_x]; read++;
					comp++;
					if (X_COMPARE(curr, temp)) {
						position_matrix[j + row * dim_x] = curr; writ++;
					}
					else {
						break;
					}
				}
				if (j != col) {
					position_matrix[j + row * dim_x] = temp; writ++;
				}
			}
	   	}
	}
    if (verbose >= 2) {
    	std::cout << "lhr  comp=" << comp << " read=" << read << " writ=" << writ << '\n';
    }
    tcomp += comp; tread += read; twrit += writ; tpass++;
}

void linear_shell_sort_on_each_col()
{
	int gaps[] = {1750, 701, 301, 132, 57, 23, 10, 4, 1};
	long int comp = 0, read = 0, writ = 0;

	for (int g = 0; g < (int) (sizeof(gaps) / sizeof(int)); g++) {
		int gap = gaps[g];
		if (gap >= dim_y) {
			continue;
		}

		for (int col = 0; col < dim_x; col++) {
			for (int row = gap; row < dim_y; row++) {
				Position temp = position_matrix[col + row * dim_x]; read++;
				int j;
				for (j = row; j >= gap; j -= gap) {
					Position curr = position_matrix[col + (j - gap) * dim_x]; read++;
					comp++;
					if (Y_COMPARE(curr, temp)) {
						position_matrix[col + j * dim_x] = curr; writ++;
					}
					else {
						break;
					}
				}
				if (j != row) {
					position_matrix[col + j * dim_x] = temp; writ++;
				}
			}
		}
	}
    if (verbose >= 2) {
    	std::cout << "lhc  comp=" << comp << " read=" << read << " writ=" << writ << '\n';
    }
    tcomp += comp; tread += read; twrit += writ; tpass++;
}

void spatial_shell_sort()
{
	int gaps[] = {1750, 701, 301, 132, 57, 23, 10, 4, 1};
	long int comp = 0, read = 0, writ = 0;

	for (int g = 0; g < (int) (sizeof(gaps) / sizeof(int)); g++) {
		int gap = gaps[g];

	   	if (gap < dim_x) {
	   		// sort rows using current gap
	   		//std::cout << "row gap=" << gap << " ";

	   		// for each row
	   		for (int row = 0; row < dim_y; row++) {
	   			// do a gapped insertion sort on this row
	   			for (int col = gap; col < dim_x; col++) {
	   				// add a[i] to the elements that have been gap sorted
	   				// save a[i] in temp and make a hole at position i
	   				Position temp = position_matrix[col + row * dim_x]; read++;
	   				// shift earlier gap-sorted elements up until the correct location for a[i] is found
	   				int j;
	   				for (j = col; j >= gap; j -= gap) {
	   					Position curr = position_matrix[(j - gap) + row * dim_x]; read++;
	   					comp++;
	   					if (X_COMPARE(curr, temp)) {
	   						position_matrix[j + row * dim_x] = curr; writ++;
	   					}
	   					else {
	   						break;
	   					}
	   				}
	   				// put temp (the original a[i]) in its correct location, if changed
	   				if (j != col) {
	   					position_matrix[j + row * dim_x] = temp; writ++;
	   				}
	   			}
	   		}
	   		tpass++;
	   	}
	   	else {
	   		//std::cout << "skipped row gap=" << gap << '\n';
	   	}

	   	if (gap < dim_y) {
	   		// sort columns using current gap
	   		//std::cout << "col gap=" << gap << " ";

	   		// for each column
	   		for (int col = 0; col < dim_x; col++) {
	   			// do a gapped insertion sort on this column
	   			for (int row = gap; row < dim_y; row++) {
	   				// add a[i] to the elements that have been gap sorted
	   				// save a[i] in temp and make a hole at position i
	   				Position temp = position_matrix[col + row * dim_x]; read++;
	   				// shift earlier gap-sorted elements up until the correct location for a[i] is found
	   				int j;
	   				for (j = row; j >= gap; j -= gap) {
	   					Position curr = position_matrix[col + (j - gap) * dim_x]; read++;
	   					comp++;
	   					if (Y_COMPARE(curr, temp)) {
	   						position_matrix[col + j * dim_x] = curr; writ++;
	   					}
	   					else {
	   						break;
	   					}
	   				}
	   				// put temp (the original a[i]) in its correct location, if changed
	   				if (j != row) {
	   					position_matrix[col + j * dim_x] = temp; writ++;
	   				}
	   			}
	   		}
	   		tpass++;
	   	}
	   	else {
	   		//std::cout << "skipped col gap=" << gap << '\n';
	   	}
	}
    if (verbose >= 2) {
    	std::cout << "psh  comp=" << comp << " read=" << read << " writ=" << writ << '\n';
    }
    tcomp += comp; tread += read; twrit += writ;
}

void full_spatial_shell_sort()
{
	if (prev_pass_func) { (*prev_pass_func)(); }
	spatial_shell_sort();
	if (post_pass_func) { (*post_pass_func)(); }

	full_insertion();
	if (verbose >= 1) {
    	std::cout << "fsh  tcomp=" << tcomp << " tread=" << tread << " twrit=" << twrit << " tpass=" << tpass << '\n';
    }
}

void full_spatial_shell_sort_skip()
{
	if (prev_pass_func) { (*prev_pass_func)(); }
	spatial_shell_sort();
	if (post_pass_func) { (*post_pass_func)(); }

	full_insertion_skip();
	if (verbose >= 1) {
    	std::cout << "fsh  tcomp=" << tcomp << " tread=" << tread << " twrit=" << twrit << " tpass=" << tpass << '\n';
    }
}

void full_shell_odd_even()
{
	if (prev_pass_func) { (*prev_pass_func)(); }
	spatial_shell_sort();
	if (post_pass_func) { (*post_pass_func)(); }

	bool is_sorted = false;
    while (! is_sorted) {
    	if (prev_pass_func) { (*prev_pass_func)(); }
    	is_sorted = odd_even_partial_sort();
		if (post_pass_func) { (*post_pass_func)(); }
    }
    if (verbose >= 1) {
    	std::cout << "fshoe  tcomp=" << tcomp << " tread=" << tread << " twrit=" << twrit << " tpass=" << tpass << '\n';
    }
}

/*------------------------ QUICKSORT ------------------------*/

int partition_on_row(int row, int left, int right, int pivot_index)
{
	Position pivot = position_matrix[pivot_index + row * dim_x]; qread++;
	position_matrix[pivot_index + row * dim_x] = position_matrix[right + row * dim_x]; qwrit++; qread++;

	int store_index = left;
    for (int i = left; i < right; i++) {
    	Position curr = position_matrix[i + row * dim_x]; qread++;
    	qcomp++;
    	if (X_COMPARE(pivot, curr)) {
    		Position store = position_matrix[store_index + row * dim_x]; qread++;
    		position_matrix[store_index + row * dim_x] = curr; qwrit++;
    		position_matrix[i + row * dim_x] = store; qwrit++;
    		store_index++;
    	}
    }
    position_matrix[right + row * dim_x] = position_matrix[store_index + row * dim_x]; qwrit++; qread++;
    position_matrix[store_index + row * dim_x] = pivot; qwrit++;
    return store_index;
}

void quicksort_on_row(int row, int left, int right)
{
    if (left < right) {
    	int pivot_index = (left + right) / 2;
        int new_pivot_index = partition_on_row(row, left, right, pivot_index);
        quicksort_on_row(row, left, new_pivot_index - 1);
        quicksort_on_row(row, new_pivot_index + 1, right);
    }
}

void linear_quicksort_on_each_row()
{
	qcomp = qread = qwrit = 0;
	for (int row = 0; row < dim_y; row++) {
		quicksort_on_row(row, 0, dim_x - 1);
	}
    if (verbose >= 2) {
    	std::cout << "lqr  comp=" << qcomp << " read=" << qread << " writ=" << qwrit << '\n';
    }
    tcomp += qcomp; tread += qread; twrit += qwrit; tpass++;
}

int partition_on_col(int col, int left, int right, int pivot_index)
{
	Position pivot = position_matrix[col + pivot_index * dim_x]; qread++;
	position_matrix[col + pivot_index * dim_x] = position_matrix[col + right * dim_x]; qwrit++; qread++;

	int store_index = left;
    for (int i = left; i < right; i++) {
    	Position curr = position_matrix[col + i * dim_x]; qread++;
    	qcomp++;
    	if (Y_COMPARE(pivot, curr)) {
    		Position store = position_matrix[col + store_index * dim_x]; qread++;
    		position_matrix[col + store_index * dim_x] = curr; qwrit++;
    		position_matrix[col + i * dim_x] = store; qwrit++;
    		store_index++;
    	}
    }
    position_matrix[col + right * dim_x] = position_matrix[col + store_index * dim_x]; qwrit++; qread++;
    position_matrix[col + store_index * dim_x] = pivot; qwrit++;
    return store_index;
}

void quicksort_on_col(int col, int left, int right)
{
    if (left < right) {
    	int pivot_index = (left + right) / 2;
        int new_pivot_index = partition_on_col(col, left, right, pivot_index);
        quicksort_on_col(col, left, new_pivot_index - 1);
        quicksort_on_col(col, new_pivot_index + 1, right);
    }
}

void linear_quicksort_on_each_col()
{
	qcomp = qread = qwrit = 0;
	for (int col = 0; col < dim_x; col++) {
		quicksort_on_col(col, 0, dim_y - 1);
	}
    if (verbose >= 2) {
    	std::cout << "lqc  comp=" << qcomp << " read=" << qread << " writ=" << qwrit << '\n';
    }
    tcomp += qcomp; tread += qread; twrit += qwrit; tpass++;
}

void partial_quicksort()
{
	linear_quicksort_on_each_row();
	linear_quicksort_on_each_col();
    if (verbose >= 1) {
    	std::cout << "pqui  tcomp=" << tcomp << " tread=" << tread << " twrit=" << twrit << " tpass=" << tpass << '\n';
    }
}

void full_quicksort()
{
	// run quick sort on each row and column until it converges
	do {
		linear_quicksort_on_each_row();
		linear_quicksort_on_each_col();
	} while (! is_spatially_sorted());
    if (verbose >= 1) {
    	std::cout << "fqu  tcomp=" << tcomp << " tread=" << tread << " twrit=" << twrit << " tpass=" << tpass << '\n';
    }
}

void full_quick_odd_even()
{
	if (prev_pass_func) { (*prev_pass_func)(); }
	linear_quicksort_on_each_row();
	//if (post_pass_func) { (*post_pass_func)(); }

	//if (prev_pass_func) { (*prev_pass_func)(); }
	linear_quicksort_on_each_col();
	if (post_pass_func) { (*post_pass_func)(); }

	bool is_sorted = false;
    while (! is_sorted) {
    	if (prev_pass_func) { (*prev_pass_func)(); }
    	is_sorted = odd_even_partial_sort();
		if (post_pass_func) { (*post_pass_func)(); }
    }
    if (verbose >= 1) {
    	std::cout << "fquoe  tcomp=" << tcomp << " tread=" << tread << " twrit=" << twrit << " tpass=" << tpass << '\n';
    }
}

void full_quick_insertion()
{
	if (prev_pass_func) { (*prev_pass_func)(); }
	linear_quicksort_on_each_row();
	//if (post_pass_func) { (*post_pass_func)(); }

	//if (prev_pass_func) { (*prev_pass_func)(); }
	linear_quicksort_on_each_col();
	if (post_pass_func) { (*post_pass_func)(); }

	full_insertion();
    if (verbose >= 1) {
    	std::cout << "fquins  tcomp=" << tcomp << " tread=" << tread << " twrit=" << twrit << " tpass=" << tpass << '\n';
    }
}

void full_quick_insertion_skip()
{
	if (prev_pass_func) { (*prev_pass_func)(); }
	linear_quicksort_on_each_row();
	//if (post_pass_func) { (*post_pass_func)(); }

	//if (prev_pass_func) { (*prev_pass_func)(); }
	linear_quicksort_on_each_col();
	if (post_pass_func) { (*post_pass_func)(); }

	full_insertion_skip();
    if (verbose >= 1) {
    	std::cout << "fquins  tcomp=" << tcomp << " tread=" << tread << " twrit=" << twrit << " tpass=" << tpass << '\n';
    }
}

void interleaved_quicksort()
{
	// alternate rows and cols -> bad results
	qcomp = qread = qwrit = 0;
	for (int i = 0; i < dim_x; i++) {
		quicksort_on_row(i, 0, dim_x - 1);
		quicksort_on_col(i, 0, dim_y - 1);
	}
    if (verbose >= 1) {
    	std::cout << "comp=" << qcomp << " read=" << qread << " writ=" << qwrit << '\n';
    }
    tcomp += qcomp; tread += qread; twrit += qwrit; tpass++;
}

/*------------------------ SIMPLE SORT ------------------------*/

int simple_global_partition_for_y(int left, int right, int pivot_index)
{
	Position pivot = position_matrix[pivot_index]; qread++;
	position_matrix[pivot_index] = position_matrix[right]; qwrit++; qread++;

	int store_index = left;
    for (int i = left; i < right; i++) {
    	Position curr = position_matrix[i]; qread++;
    	qcomp++;
    	if (Y_COMPARE(pivot, curr)) {
    		Position store = position_matrix[store_index]; qread++;
    		position_matrix[store_index] = curr; qwrit++;
    		position_matrix[i] = store; qwrit++;
    		store_index++;
    	}
    }
    position_matrix[right] = position_matrix[store_index]; qwrit++; qread++;
    position_matrix[store_index] = pivot; qwrit++;
    return store_index;
}

void simple_global_quicksort_for_y(int left, int right)
{
    if (left < right) {
    	int pivot_index = (left + right) / 2;
        int new_pivot_index = simple_global_partition_for_y(left, right, pivot_index);
        simple_global_quicksort_for_y(left, new_pivot_index - 1);
        simple_global_quicksort_for_y(new_pivot_index + 1, right);
    }
}

int simple_partition_on_row_for_x(int row, int left, int right, int pivot_index)
{
	Position pivot = position_matrix[pivot_index + row * dim_x]; qread++;
	position_matrix[pivot_index + row * dim_x] = position_matrix[right + row * dim_x]; qwrit++; qread++;

	int store_index = left;
    for (int i = left; i < right; i++) {
    	Position curr = position_matrix[i + row * dim_x]; qread++;
    	qcomp++;
    	if (X_COMPARE(pivot, curr)) {
    		Position store = position_matrix[store_index + row * dim_x]; qread++;
    		position_matrix[store_index + row * dim_x] = curr; qwrit++;
    		position_matrix[i + row * dim_x] = store; qwrit++;
    		store_index++;
    	}
    }
    position_matrix[right + row * dim_x] = position_matrix[store_index + row * dim_x]; qwrit++; qread++;
    position_matrix[store_index + row * dim_x] = pivot; qwrit++;
    return store_index;
}

void simple_quicksort_on_row_for_x(int row, int left, int right)
{
    if (left < right) {
    	int pivot_index = (left + right) / 2;
        int new_pivot_index = simple_partition_on_row_for_x(row, left, right, pivot_index);
        simple_quicksort_on_row_for_x(row, left, new_pivot_index - 1);
        simple_quicksort_on_row_for_x(row, new_pivot_index + 1, right);
    }
}

void simple_sort()
{
	qcomp = qread = qwrit = 0;
	simple_global_quicksort_for_y(0, position_counter - 1);
    //if (verbose >= 2) {
    //std::cout << "global y  comp=" << qcomp << " read=" << qread << " writ=" << qwrit << '\n';
    tcomp += qcomp; tread += qread; twrit += qwrit; tpass++;
    //}

	qcomp = qread = qwrit = 0;
	for (int row = 0; row < dim_y; row++) {
		quicksort_on_row(row, 0, dim_x - 1);
	}
    //if (verbose >= 2) {
    //std::cout << "row x     comp=" << qcomp << " read=" << qread << " writ=" << qwrit << '\n';
    tcomp += qcomp; tread += qread; twrit += qwrit; tpass++;
    //}
    //std::cout << "simple sort tcomp=" << tcomp << " tread=" << tread << " twrit=" << twrit << '\n';

    //if (is_spatially_sorted()) {
        //std::cout << "is sorted!\n";
    //}
    //else {
        //std::cout << "is not sorted...\n";
    //}
}

/*------------------------ NEIGHBORHOOD FUNCTIONS ------------------------*/

int *get_hard_4_neighborhood(int x, int y)
{
	int nx, ny;
	int *p = neighborhood_4;
	// bottom
	nx = x;
	ny = y - 1;
	if (nx >= 0 && nx < dim_x && ny >= 0 && ny < dim_y) {
		*p++ = nx + ny * dim_x;
	}
	// top
	//nx = x;
	ny = y + 1;
	if (nx >= 0 && nx < dim_x && ny >= 0 && ny < dim_y) {
		*p++ = nx + ny * dim_x;
	}
	// left
	nx = x - 1;
	ny = y;
	if (nx >= 0 && nx < dim_x && ny >= 0 && ny < dim_y) {
		*p++ = nx + ny * dim_x;
	}
	// right
	nx = x + 1;
	//ny = y;
	if (nx >= 0 && nx < dim_x && ny >= 0 && ny < dim_y) {
		*p++ = nx + ny * dim_x;
	}
	*p = -1; // mark list end
	return neighborhood_4;
}

int *get_hard_8_neighborhood(int x, int y)
{
	int self = x + y * dim_x;
	int min_x = x - 1, max_x = x + 1;
	int min_y = y - 1, max_y = y + 1;
	if (min_x < 0) {
		min_x = 0;
	}
	else if (max_x >= dim_x) {
		max_x = dim_x - 1;
	}
	if (min_y < 0) {
		min_y = 0;
	}
	else if (max_y >= dim_y) {
		max_y = dim_y - 1;
	}
	int *p = neighborhood_8;
	for (int j = min_y; j <= max_y; j++) {
		for (int i = min_x; i <= max_x; i++) {
			int pos = i + j * dim_x;
			if (pos != self) {
				*p++ = pos;
			}
		}
	}
	*p = -1; // mark list end
	return neighborhood_8;
}

int *get_hard_24_neighborhood(int x, int y)
{
	int self = x + y * dim_x;
	int min_x = x - 2, max_x = x + 2;
	int min_y = y - 2, max_y = y + 2;
	if (min_x < 0) {
		min_x = 0;
	}
	else if (max_x >= dim_x) {
		max_x = dim_x - 1;
	}
	if (min_y < 0) {
		min_y = 0;
	}
	else if (max_y >= dim_y) {
		max_y = dim_y - 1;
	}
	int *p = neighborhood_24;
	for (int j = min_y; j <= max_y; j++) {
		for (int i = min_x; i <= max_x; i++) {
			int pos = i + j * dim_x;
			if (pos != self) {
				*p++ = pos;
			}
		}
	}
	*p = -1; // mark list end
	return neighborhood_24;
}

int *get_hard_48_neighborhood(int x, int y)
{
	int self = x + y * dim_x;
	int min_x = x - 3, max_x = x + 3;
	int min_y = y - 3, max_y = y + 3;
	if (min_x < 0) {
		min_x = 0;
	}
	else if (max_x >= dim_x) {
		max_x = dim_x - 1;
	}
	if (min_y < 0) {
		min_y = 0;
	}
	else if (max_y >= dim_y) {
		max_y = dim_y - 1;
	}
	int *p = neighborhood_48;
	for (int j = min_y; j <= max_y; j++) {
		for (int i = min_x; i <= max_x; i++) {
			int pos = i + j * dim_x;
			if (pos != self) {
				*p++ = pos;
			}
		}
	}
	*p = -1; // mark list end
	return neighborhood_48;
}

int *get_hard_80_neighborhood(int x, int y)
{
	int self = x + y * dim_x;
	int min_x = x - 4, max_x = x + 4;
	int min_y = y - 4, max_y = y + 4;
	if (min_x < 0) {
		min_x = 0;
	}
	else if (max_x >= dim_x) {
		max_x = dim_x - 1;
	}
	if (min_y < 0) {
		min_y = 0;
	}
	else if (max_y >= dim_y) {
		max_y = dim_y - 1;
	}
	int *p = neighborhood_80;
	for (int j = min_y; j <= max_y; j++) {
		for (int i = min_x; i <= max_x; i++) {
			int pos = i + j * dim_x;
			if (pos != self) {
				*p++ = pos;
			}
		}
	}
	*p = -1; // mark list end
	return neighborhood_80;
}

int *get_hard_120_neighborhood(int x, int y)
{
	int self = x + y * dim_x;
	int min_x = x - 5, max_x = x + 5;
	int min_y = y - 5, max_y = y + 5;
	if (min_x < 0) {
		min_x = 0;
	}
	else if (max_x >= dim_x) {
		max_x = dim_x - 1;
	}
	if (min_y < 0) {
		min_y = 0;
	}
	else if (max_y >= dim_y) {
		max_y = dim_y - 1;
	}
	int *p = neighborhood_120;
	for (int j = min_y; j <= max_y; j++) {
		for (int i = min_x; i <= max_x; i++) {
			int pos = i + j * dim_x;
			if (pos != self) {
				*p++ = pos;
			}
		}
	}
	*p = -1; // mark list end
	return neighborhood_120;
}

void get_hard_3d_neighborhood(int x, int y, int z, int n_size, int *candidate)
{
	int self = x + y * dim_x + z * dim_x * dim_y;
	int min_x = x - n_size, max_x = x + n_size;
	int min_y = y - n_size, max_y = y + n_size;
	int min_z = z - n_size, max_z = z + n_size;
	if (min_x < 0) {
		min_x = 0;
	}
	else if (max_x >= dim_x) {
		max_x = dim_x - 1;
	}
	if (min_y < 0) {
		min_y = 0;
	}
	else if (max_y >= dim_y) {
		max_y = dim_y - 1;
	}
	if (min_z < 0) {
		min_z = 0;
	}
	else if (max_z >= dim_z) {
		max_z = dim_z - 1;
	}
	for (int k = min_z; k <= max_z; k++) {
		for (int j = min_y; j <= max_y; j++) {
			for (int i = min_x; i <= max_x; i++) {
				int pos = i + j * dim_x + k * dim_x * dim_y;
				if (pos != self) {
					*candidate++ = pos;
				}
			}
		}
	}
	*candidate = -1; // mark list end
}

// keep neighborhood fully inside matrix

int *get_soft_4_neighborhood(int x, int y)
{
	int self = x + y * dim_x;
	if (x < 1) {
		x = 1;
	}
	else if (x >= dim_x - 1) {
		x = dim_x - 1 - 1;
	}
	if (y < 1) {
		y = 1;
	}
	else if (y >= dim_y - 1) {
		y = dim_y - 1 - 1;
	}
	int *p = neighborhood_4;
	int pos;
	// bottom
	pos = x + (y - 1) * dim_x;
	if (pos != self) {
		*p++ = pos;
	}
	// top
	pos = x + (y + 1) * dim_x;
	if (pos != self) {
		*p++ = pos;
	}
	// left
	pos = (x - 1) + y * dim_x;
	if (pos != self) {
		*p++ = pos;
	}
	// right
	pos = (x + 1) + y * dim_x;
	if (pos != self) {
		*p++ = pos;
	}
	*p = -1; // mark list end
	return neighborhood_4;
}

int *get_soft_8_neighborhood(int x, int y)
{
	int self = x + y * dim_x;
	if (x < 1) {
		x = 1;
	}
	else if (x >= dim_x - 1) {
		x = dim_x - 2;
	}
	if (y < 1) {
		y = 1;
	}
	else if (y >= dim_y - 1) {
		y = dim_y - 2;
	}
	int *p = neighborhood_8;
	for (int j = -1; j <= 1; j++) {
		for (int i = -1; i <= 1; i++) {
			int pos = (x + i) + (y + j) * dim_x;
			if (pos != self) {
				*p++ = pos;
			}
		}
	}
	*p = -1; // mark list end
	return neighborhood_8;
}

int *get_soft_24_neighborhood(int x, int y)
{
	int self = x + y * dim_x;
	if (x < 2) {
		x = 2;
	}
	else if (x >= dim_x - 2) {
		x = dim_x - 3;
	}
	if (y < 2) {
		y = 2;
	}
	else if (y >= dim_y - 2) {
		y = dim_y - 3;
	}
	int *p = neighborhood_24;
	for (int j = -2; j <= 2; j++) {
		for (int i = -2; i <= 2; i++) {
			int pos = (x + i) + (y + j) * dim_x;
			if (pos != self) {
				*p++ = pos;
			}
		}
	}
	*p = -1; // mark list end
	return neighborhood_24;
}

int *get_soft_48_neighborhood(int x, int y)
{
	int self = x + y * dim_x;
	if (x < 3) {
		x = 3;
	}
	else if (x >= dim_x - 3) {
		x = dim_x - 4;
	}
	if (y < 3) {
		y = 3;
	}
	else if (y >= dim_y - 3) {
		y = dim_y - 4;
	}
	int *p = neighborhood_48;
	for (int j = -3; j <= 3; j++) {
		for (int i = -3; i <= 3; i++) {
			int pos = (x + i) + (y + j) * dim_x;
			if (pos != self) {
				*p++ = pos;
			}
		}
	}
	*p = -1; // mark list end
	return neighborhood_48;
}

int *get_soft_80_neighborhood(int x, int y)
{
	int self = x + y * dim_x;
	if (x < 4) {
		x = 4;
	}
	else if (x >= dim_x - 4) {
		x = dim_x - 5;
	}
	if (y < 4) {
		y = 4;
	}
	else if (y >= dim_y - 4) {
		y = dim_y - 5;
	}
	int *p = neighborhood_80;
	for (int j = -4; j <= 4; j++) {
		for (int i = -4; i <= 4; i++) {
			int pos = (x + i) + (y + j) * dim_x;
			if (pos != self) {
				*p++ = pos;
			}
		}
	}
	*p = -1; // mark list end
	return neighborhood_80;
}

int *get_soft_120_neighborhood(int x, int y)
{
	int self = x + y * dim_x;
	if (x < 5) {
		x = 5;
	}
	else if (x >= dim_x - 5) {
		x = dim_x - 6;
	}
	if (y < 5) {
		y = 5;
	}
	else if (y >= dim_y - 5) {
		y = dim_y - 6;
	}
	int *p = neighborhood_120;
	for (int j = -5; j <= 5; j++) {
		for (int i = -5; i <= 5; i++) {
			int pos = (x + i) + (y + j) * dim_x;
			if (pos != self) {
				*p++ = pos;
			}
		}
	}
	*p = -1; // mark list end
	return neighborhood_120;
}

float get_distance(Position& a, Position& b)
{
    return sqrtf((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
}

float get_squared_distance(Position& a, Position& b)
{
    return (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y);
}

