/*-------------------------------- INCLUDES --------------------------------*/

#include <climits>
#include <cmath>
#include <iostream>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "common.hpp"

/*-------------------------------- IMPORTED VARIABLES --------------------------------*/

extern int dim_x, dim_y, dim_z;
extern Position *position_matrix;
extern Position3d *position_matrix_3d;
extern int position_counter;

/*-------------------------------- EXPORTED VARIABLES --------------------------------*/

float diagonal_angle = 45;
float diagonal_width = 0.02;

/*-------------------------------- LOCAL VARIABLES --------------------------------*/

gsl_rng *rng = NULL;

/*------------------------ RANDOM NUMBER FUNCTIONS ------------------------*/

void random_init()
{
	rng = gsl_rng_alloc(gsl_rng_mt19937);
}

void random_seed(int seed)
{
	if (rng == NULL) {
		random_init();
	}
	gsl_rng_set(rng, (unsigned long int) seed);
}

void random_done()
{
	gsl_rng_free(rng);
}

float random_range(float min, float max)
{
    return min + (max - min) * ((float) gsl_rng_get(rng) / UINT_MAX);
}

/*------------------------ TEST CASE DISTRIBUTIONS ------------------------*/

void setup_uniform()
{
	for (int row = 0; row < dim_y; row++) {
    	for (int col = 0; col < dim_x; col++) {
    		int pos = col + row * dim_x;
    		position_matrix[pos].id = pos;

    		position_matrix[pos].x = random_range(0, 1);
    		position_matrix[pos].y = random_range(0, 1);
    	}
    }
}

void setup_uniform_3d()
{
	for (int layer = 0; layer < dim_z; layer++) {
		for (int row = 0; row < dim_y; row++) {
			for (int col = 0; col < dim_x; col++) {
				int pos = col + row * dim_x + layer * dim_x * dim_y;
				position_matrix_3d[pos].id = pos;

				position_matrix_3d[pos].x = random_range(0, 1);
				position_matrix_3d[pos].y = random_range(0, 1);
				position_matrix_3d[pos].z = random_range(0, 1);
			}
		}
	}
}

void setup_shuffle()
{
	for (int row = 0; row < dim_y; row++) {
    	for (int col = 0; col < dim_x; col++) {
    		int pos = col + row * dim_x;
    		position_matrix[pos].id = pos;

    		// perfect set, already ordered
    		position_matrix[pos].x = (float) col / (dim_x - 1);
    		position_matrix[pos].y = (float) row / (dim_y - 1);
    		
    		// perfect set, transposed: patological case which is stable until small perturbation
    		//position_matrix[pos].x = (float) row / (dim_y - 1);
    		//position_matrix[pos].y = (float) col / (dim_x - 1);

    		// perfect set, horizontally reflected: sorts to normal order
    		//position_matrix[pos].x = (float) (dim_x - 1 - col) / (dim_x - 1);
    		//position_matrix[pos].y = (float) row / (dim_y - 1);

    		// perfect set, vertically reflected: sorts to normal order
    		//position_matrix[pos].x = (float) col / (dim_x - 1);
    		//position_matrix[pos].y = (float) (dim_y - 1 - row) / (dim_y - 1);

    		// perfect set, rotated 90 degrees clockwise: sorts to transposed
    		//position_matrix[pos].x = (float) (dim_y - 1 - row) / (dim_y - 1);
    		//position_matrix[pos].y = (float) col / (dim_x - 1);

    		// perfect set, rotated 180 degrees clockwise: sorts to normal order
    		//position_matrix[pos].x = (float) (dim_x - 1 - col) / (dim_x - 1);
    		//position_matrix[pos].y = (float) (dim_y - 1 - row) / (dim_y - 1);

    		// perfect set, rotated 270 degrees clockwise: sorts to transposed
    		//position_matrix[pos].x = (float) row / (dim_y - 1);
    		//position_matrix[pos].y = (float) (dim_x - 1 - col) / (dim_x - 1);
    	}
    }
	gsl_ran_shuffle(rng, position_matrix, position_counter, sizeof(Position));
	for (int i = 0; i < position_counter; i++) {
		position_matrix[i].id = i;
	}
}

void setup_gaussian()
{
	for (int row = 0; row < dim_y; row++) {
    	for (int col = 0; col < dim_x; col++) {
    		int pos = col + row * dim_x;
    		position_matrix[pos].id = pos;

    		position_matrix[pos].x = gsl_ran_gaussian(rng, 0.2) + 0.5;
    		position_matrix[pos].y = gsl_ran_gaussian(rng, 0.2) + 0.5;
    	}
    }
}

void setup_inside_circle()
{
	float x, y;
	for (int row = 0; row < dim_y; row++) {
    	for (int col = 0; col < dim_x; col++) {
    		int pos = col + row * dim_x;
    		position_matrix[pos].id = pos;

    		while (true) {
    			x = random_range(0, 1);
    			y = random_range(0, 1);
    			if (((x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5)) <= 0.25) {
    				break;
    			}
    		}
    		position_matrix[pos].x = x;
    		position_matrix[pos].y = y;
    	}
    }
}

void setup_outside_circle()
{
	float x, y;
	for (int row = 0; row < dim_y; row++) {
    	for (int col = 0; col < dim_x; col++) {
    		int pos = col + row * dim_x;
    		position_matrix[pos].id = pos;

    		while (true) {
    			x = random_range(0, 1);
    			y = random_range(0, 1);
    			if (((x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5)) >= 0.25) {
    				break;
    			}
    		}
    		position_matrix[pos].x = x;
    		position_matrix[pos].y = y;
    	}
    }
}

void setup_random_clusters()
{
	int n = 5;
	float cx[n], cy[n], cr[n], x, y;

	int i = 0, j;
	while (i < n) {
		cr[i] = random_range(0.10, 0.20);
		cx[i] = random_range(0 + cr[i], 1 - cr[i]);
		cy[i] = random_range(0 + cr[i], 1 - cr[i]);
		for (j = 0; j < i; j++) {
			if (((cx[i] - cx[j]) * (cx[i] - cx[j]) + (cy[i] - cy[j]) * (cy[i] - cy[j])) <= (cr[i] + cr[j]) * (cr[i] + cr[j])) {
				break; // clusters are not disjoint
			}
		}
		if (j == i) {
			i++; // if all clusters are disjoint, add another one
		}
	}

	for (int row = 0; row < dim_y; row++) {
    	for (int col = 0; col < dim_x; col++) {
    		int pos = col + row * dim_x;
    		position_matrix[pos].id = pos;

    		bool found = false;
    		while (! found) {
    			x = random_range(0, 1);
    			y = random_range(0, 1);
    			for (int i = 0; i < n; i++) {
    				if (((x - cx[i]) * (x - cx[i]) + (y - cy[i]) * (y - cy[i])) <= cr[i] * cr[i]) {
    					found = true;
    					break;
    				}
    			}
    		}
    		position_matrix[pos].x = x;
    		position_matrix[pos].y = y;
    	}
    }
}

void setup_random_holes()
{
	int n = 5;
	float cx[n], cy[n], cr[n], x, y;

	int i = 0, j;
	while (i < n) {
		cr[i] = random_range(0.10, 0.20);
		cx[i] = random_range(0 + cr[i], 1 - cr[i]);
		cy[i] = random_range(0 + cr[i], 1 - cr[i]);
		for (j = 0; j < i; j++) {
			if (((cx[i] - cx[j]) * (cx[i] - cx[j]) + (cy[i] - cy[j]) * (cy[i] - cy[j])) <= (cr[i] + cr[j]) * (cr[i] + cr[j])) {
				break; // clusters are not disjoint
			}
		}
		if (j == i) {
			i++; // if all clusters are disjoint, add another one
		}
	}

	for (int row = 0; row < dim_y; row++) {
    	for (int col = 0; col < dim_x; col++) {
    		int pos = col + row * dim_x;
    		position_matrix[pos].id = pos;

    		bool inside = true;
    		while (inside) {
        		inside = false;
    			x = random_range(0, 1);
    			y = random_range(0, 1);
    			for (int i = 0; i < n; i++) {
    				if (((x - cx[i]) * (x - cx[i]) + (y - cy[i]) * (y - cy[i])) <= cr[i] * cr[i]) {
    					inside = true;
    					break;
    				}
    			}
    		}
    		position_matrix[pos].x = x;
    		position_matrix[pos].y = y;
    	}
    }
}

void setup_streets()
{
	int n = 5;
	float v_pos[n], v_siz[n], h_pos[n], h_siz[n], x, y;

	for (int i = 0; i < n; i++) {
		v_pos[i] = (i + 1.0) / (n + 1.0) + random_range(-0.05, 0.05);
		v_siz[i] = random_range(0.01, 0.05);
		h_pos[i] = (i + 1.0) / (n + 1.0) + random_range(-0.05, 0.05);
		h_siz[i] = random_range(0.01, 0.05);
	}

	for (int row = 0; row < dim_y; row++) {
    	for (int col = 0; col < dim_x; col++) {
    		int pos = col + row * dim_x;
    		position_matrix[pos].id = pos;

    		bool found = false;
    		while (! found) {
    			x = random_range(0, 1);
    			y = random_range(0, 1);
    			for (int i = 0; i < n; i++) {
    				if ((v_pos[i] - v_siz[i] <= x) && (x <= v_pos[i] + v_siz[i])) {
    					found = true;
    					break;
    				}
    				if ((h_pos[i] - h_siz[i] <= y) && (y <= h_pos[i] + h_siz[i])) {
    					found = true;
    					break;
    				}
    			}
    		}
    		position_matrix[pos].x = x;
    		position_matrix[pos].y = y;
    	}
    }
}

void setup_tilted()
{
	float angle = random_range(0, 90) * M_PI / 180;
	float s = sinf(angle);
	float c = cosf(angle);
	int n = 5;
	float v_pos[n], v_siz[n], h_pos[n], h_siz[n], x, y;

	for (int i = 0; i < n; i++) {
		v_pos[i] = (i + 1.0) / (n + 1.0) + random_range(-0.05, 0.05);
		v_siz[i] = random_range(0.01, 0.05);
		h_pos[i] = (i + 1.0) / (n + 1.0) + random_range(-0.05, 0.05);
		h_siz[i] = random_range(0.01, 0.05);
	}

	for (int row = 0; row < dim_y; row++) {
    	for (int col = 0; col < dim_x; col++) {
    		int pos = col + row * dim_x;
    		position_matrix[pos].id = pos;

    		bool found = false;
    		while (! found) {
    			x = random_range(0, 1);
    			y = random_range(0, 1);
    			for (int i = 0; i < n; i++) {
    				if ((v_pos[i] - v_siz[i] <= x) && (x <= v_pos[i] + v_siz[i])) {
    					found = true;
    					break;
    				}
    				if ((h_pos[i] - h_siz[i] <= y) && (y <= h_pos[i] + h_siz[i])) {
    					found = true;
    					break;
    				}
    			}
    		}
    		position_matrix[pos].x = c * (x - 0.5) - s * (y - 0.5) + 0.5;
    		position_matrix[pos].y = s * (x - 0.5) + c * (y - 0.5) + 0.5;
    	}
    }
}

void setup_diagonal()
{
	float angle = diagonal_angle * M_PI / 180;
	float s = sinf(angle);
	float c = cosf(angle);
	float h_pos, h_siz, x, y;

	h_pos = 0.5;
	h_siz = diagonal_width;

	for (int row = 0; row < dim_y; row++) {
    	for (int col = 0; col < dim_x; col++) {
    		int pos = col + row * dim_x;
    		position_matrix[pos].id = pos;

    		while (true) {
        		x = random_range(0, 1);
        		y = random_range(0, 1);
    			if ((h_pos - h_siz <= y) && (y <= h_pos + h_siz)) {
    				break;
    			}
    		}
    		position_matrix[pos].x = c * (x - 0.5) - s * (y - 0.5) + 0.5;
    		position_matrix[pos].y = s * (x - 0.5) + c * (y - 0.5) + 0.5;
    	}
    }
}

void setup_cross()
{
	float angle = 45 * M_PI / 180;
	float s = sinf(angle);
	float c = cosf(angle);
	float v_pos, v_siz, h_pos, h_siz, x, y;

	v_pos = 0.5;
	v_siz = 0.02;
	h_pos = 0.5;
	h_siz = 0.02;

	for (int row = 0; row < dim_y; row++) {
    	for (int col = 0; col < dim_x; col++) {
    		int pos = col + row * dim_x;
    		position_matrix[pos].id = pos;

    		while (true) {
    			x = random_range(0, 1);
    			y = random_range(0, 1);
    			if ((v_pos - v_siz <= x) && (x <= v_pos + v_siz)) {
    				break;
    			}
    			if ((h_pos - h_siz <= y) && (y <= h_pos + h_siz)) {
    				break;
    			}
    		}
    		position_matrix[pos].x = c * (x - 0.5) - s * (y - 0.5) + 0.5;
    		position_matrix[pos].y = s * (x - 0.5) + c * (y - 0.5) + 0.5;
    	}
    }
}

void setup_triangle()
{
	float x, y;
	float angle = 45 * M_PI / 180;
	float s = sinf(angle);
	float c = cosf(angle);

	for (int row = 0; row < dim_y; row++) {
    	for (int col = 0; col < dim_x; col++) {
    		int pos = col + row * dim_x;
    		position_matrix[pos].id = pos;

    		while (true) {
    			x = random_range(0, 1);
    			y = random_range(0, 1);

    			if (s * (x - 0.5) + c * (y - 0.5) + 0.5 <= 0.5) {
    				break;
    			}
    		}
    		position_matrix[pos].x = x;
    		position_matrix[pos].y = y;
    	}
    }
}

void setup_gradient()
{
	float x, y, p;
	for (int row = 0; row < dim_y; row++) {
    	for (int col = 0; col < dim_x; col++) {
    		int pos = col + row * dim_x;
    		position_matrix[pos].id = pos;

    		while (true) {
    			x = random_range(0, 1);
    			y = random_range(0, 1);
    			p = random_range(0, 1);
    			if (p <= (x + y) / 3 - 0.1) {
    				break;
    			}
    		}
    		position_matrix[pos].x = x;
    		position_matrix[pos].y = y;
    	}
    }
}

void setup_compact()
{
	for (int row = 0; row < dim_y; row++) {
    	for (int col = 0; col < dim_x; col++) {
    		int pos = col + row * dim_x;
    		position_matrix[pos].id = pos;

    		position_matrix[pos].x = random_range(0, 1);
    		position_matrix[pos].y = random_range(0, 0.2);
    	}
    }
}

