#include <math.h>
#include "stdlib.h"
#include "cdtw.h"

#include "stdio.h"


double _round(double number) {
    return floor(number + 0.5);
}

/*
 * Function:  min3
 * --------------------
 *  Finds the minimum of 3 doubles

 *  double x:    x
 *  double y:    y
 *  double z:    z

 *  returns:  min3(x,y,z)
 */
double min3(double x, double y, double z) {
    if (x < y && x < z)
        return x;
    else if (y < x && y < z)
        return y;
    else
        return z;
}

/*
 * Function:  min3
 * --------------------
 *  Finds the minimum of n doubles

 *  double *arr:    array
 *  int n: length of the array

 *  returns:  min(arr)
 */
double min_n(double *arr, int n) {
    double min = arr[0];
    int i = 1;
    for (; i < n; i++) {
        if (min > arr[i])
            min = arr[i];
    }
    return min;
}


/*
 * Function:  min2idx
 * --------------------
 *  Finds the minimum of 2 doubles and its position (0 or 1)

 *  double a: a
 *  double b: b

 *  returns:  struct t_item, t_item = { min(a,b), position }
 */
struct t_item min2idx(double a, double b) /*0, 1*/
{
    struct t_item item;
    if (a < b) {
        item.val = a;
        item.idx = 0;
    }
    else {
        item.val = b;
        item.idx = 1;
    }
    return item;
}

/*
 * Function:  min3idx
 * --------------------
 *  Finds the minimum of 3 doubles and its position (0,1 or 2)

 *  double x: x
 *  double y: y
 *  double z: z

 *  returns:  struct t_item, t_item = { min(x,y,z), position }
 */
struct t_item min3idx(double x, double y, double z) /*0, 1, 2*/
{
    struct t_item item;
    if (x < y && x < z) {
        item.val = x;
        item.idx = 0;
    }
    else if (y < x && y < z) {
        item.val = y;
        item.idx = 1;
    }
    else {
        item.val = z;
        item.idx = 2;
    }
    return item;
}

/*
 * Function:  min3idx
 * --------------------
 *  Finds the minimum of n doubles and its position (0..n-1)

 *  double x: array
 *  int n: size of the arr

 *  returns:  struct t_item, t_item = { min(arr), position }
 */
struct t_item min_nidx(double *arr, int n) {
    struct t_item item;
    int i = 0;
    item.val = arr[0];
    item.idx = 0;
    for (i = 1; i < n; i++) {
        if (item.val > arr[i]) {
            item.val = arr[i];
            item.idx = i;
        }
    }
    return item;
}


/*
 * Function:  dpw
 * --------------------

 * NOTE: ALL STEP FUNCTIONS have the same args defined in macro
 *_DP_ARGS(cdtw.h). Step functions like step_pattern_type and
 * step_pattern_typedir are pretty similar, step_pattern_type are used in
 * computing dtw without path(without traceback). step_pattern_typedir are
 * used in computing dtw with path(traceback)

 *  Step patern dpw - weights:
 *  min(
 *      cost_matrix[i][j-1]   +   a*d(r[i],q[j])
        cost_matrix[i-1][j]   +   b*d(r[i],q[j]),
        cost_matrix[i-1][j-1] +   c*d(r[i],q[j])
       )
 * where a,b,c are weights
 *
 * double* ref: reference sequence
 * double* query: query sequence
 * double* cost_matrix: cost matrix
 * int i: row index of the cost matrix
 * int j: column index of the cost matrix
 * int ncols: number of columns (1 means 1 dimensional array)
 * int offset: extra size of the cost matrix
 * int size2:  cost matrix columns count
 * double (*dist)(double a, double b): poiter to distance function

 * returns:  double, value to assign cost_matrix[i][j]
*/
double dpw(_DP_ARGS) {
    double d = dist(ref, i, query, j, ncols);
    return min3(cost_matrix[idx(i, j - 1, size2)] + t_s->weights.a * d,
                cost_matrix[idx(i - 1, j, size2)] + t_s->weights.b * d,
                cost_matrix[idx(i - 1, j - 1, size2)] + t_s->weights.c * d);
}


/*
 * Function:  dp1
 * --------------------
 *  Step patern dp1:
 *  min(
 *      cost_matrix[i][j-1]   +   d(r[i],q[j])
        cost_matrix[i-1][j]   +   d(r[i],q[j]),
        cost_matrix[i-1][j-1] + 2*d(r[i],q[j])
       )
 * see doc for the dpw
 * returns:  double, value to assign cost_matrix[i][j]
*/
double dp1(_DP_ARGS) {
    double d = dist(ref, i, query, j, ncols);
    return min3(cost_matrix[idx(i, j - 1, size2)] + d,
                cost_matrix[idx(i - 1, j, size2)] + d,
                cost_matrix[idx(i - 1, j - 1, size2)] + 2 * d);
}

/*
 * Function:  dp2
 * --------------------
 *  Step patern dp2:
 *  min(
 *      cost_matrix[i][j-1]   +   d(r[i],q[j])
        cost_matrix[i-1][j]   +   d(r[i],q[j]),
        cost_matrix[i-1][j-1] +   d(r[i],q[j])
       )
 * see doc for the dpw
*/
double dp2(_DP_ARGS) {
    double d = dist(ref, i, query, j, ncols);
    return (min3(cost_matrix[idx(i, j - 1, size2)],
                 cost_matrix[idx(i - 1, j, size2)],
                 cost_matrix[idx(i - 1, j - 1, size2)]) + d);
}

/*
 * Function:  dp3
 * --------------------
 *  Step patern dp3:
 *  min(
 *      cost_matrix[i][j-1]   +   d(r[i],q[j])
        cost_matrix[i-1][j]   +   d(r[i],q[j]),
       )
 * see doc for the dpw
*/
double dp3(_DP_ARGS) {
    double d = dist(ref, i, query, j, ncols);
    return (min2(cost_matrix[idx(i, j - 1, size2)],
                 cost_matrix[idx(i - 1, j, size2)]) + d);
}

/*
 * Function:  dp1dir
 * --------------------
 * see doc for the dpw
*/
struct t_item dp1dir(_DP_ARGS) {
    double d = dist(ref, i, query, j, ncols);
    return min3idx(cost_matrix[idx(i, j - 1, size2)] + d,
                   cost_matrix[idx(i - 1, j, size2)] + d,
                   cost_matrix[idx(i - 1, j - 1, size2)] + 2 * d);
}

/*
 * Function:  dp2dir
 * --------------------
 * see doc for the dpw,dp2
*/
struct t_item dp2dir(_DP_ARGS) {
    double d = dist(ref, i, query, j, ncols);
    return min3idx(cost_matrix[idx(i, j - 1, size2)] + d,   //left
                   cost_matrix[idx(i - 1, j, size2)] + d,   //up
                   cost_matrix[idx(i - 1, j - 1, size2)] + d);  //diag
}

/*
 * Function:  dp3dir
 * --------------------
 * see doc for the dpw,dp3
*/
struct t_item dp3dir(_DP_ARGS) {
    double d = dist(ref, i, query, j, ncols);
    return min2idx(cost_matrix[idx(i, j - 1, size2)] + d,
                   cost_matrix[idx(i - 1, j, size2)] + d);
}


/*
 * Function:  p0sym
 * --------------------
 * Sakoe-Chiba classification p = 0, symmetric step pattern
 * This function is alias for the dp1
 * see doc for the dpw
*/
double p0_sym(_DP_ARGS) {
    return dp1(ref, query, cost_matrix, i, j, ncols, t_s, size2, dist);
}

/*
 * Function:  p0asym
 * --------------------
 * Sakoe-Chiba classification p = 0, asymmetric step pattern:
 *  min(
 *      cost_matrix[i][j-1]   +   0
        cost_matrix[i-1][j]   +   d(r[i],q[j]),
        cost_matrix[i-1][j-1] +   d(r[i],q[j]),
       )
 * see doc for the dpw
*/
double p0_asym(_DP_ARGS) {
    double d = dist(ref, i, query, j, ncols);
    return min3(cost_matrix[idx(i, j - 1, size2)],         //0
                cost_matrix[idx(i - 1, j, size2)] + d,   //1
                cost_matrix[idx(i - 1, j - 1, size2)] + d);  //2
}

/*
 * Function:  p1div2_sym
 * --------------------
 * Sakoe-Chiba classification p = 0.5, symmetric step pattern:
 *  min(
 *      cost_matrix[i-1][j-3] + 2d(r[i],q[j-2]) + d(r[i],q[j-1]) + d(r[i],q[j]),
 *      cost_matrix[i-1][j-2] + 2d(r[i],q[j-1]) + d(r[i],q[j]),
 *      cost_matrix[i-1][j-1] + 2d(r[i],q[j]),
 *      cost_matrix[i-2][j-1] + 2d(r[i-1],q[j]) + d(r[i],q[j]),
 *      cost_matrix[i-3][j-1] + 2d(r[i-2],q[j]) + d(r[i-1],q[j]) + d(r[i],q[j])
       )
 * see doc for the dpw
*/
double p1div2_sym(_DP_ARGS) {
    double d00 = dist(ref, i, query, j, ncols);
    double d01 = dist(ref, i, query, j-1, ncols);
    double d02 = dist(ref, i, query, j-2, ncols);
    double d10 = dist(ref, i-1, query, j, ncols);
    double d20 = dist(ref, i-2, query, j, ncols);
    double arr[5] = {
            cost_matrix[idx(i - 1, j - 3, size2)] + 2 * d02 + d01 + d00,
            cost_matrix[idx(i - 1, j - 2, size2)] + 2 * d01 + d00,
            cost_matrix[idx(i - 1, j - 1, size2)] + 2 * d00,
            cost_matrix[idx(i - 2, j - 1, size2)] + 2 * d10 + d00,
            cost_matrix[idx(i - 3, j - 1, size2)] + 2 * d20 + d10 + d00
    };
    return min_n(arr, 5);

}

/*
 * Function:  p1div2_asym
 * --------------------
 * Sakoe-Chiba classification p = 0.5, asymmetric step pattern:
 *  min(
 *   cost_matrix[i-1][j-3] + (d(r[i],q[j-2]) + d(r[i],q[j-1]) + d(r[i],q[j]))/3
 *   cost_matrix[i-1][j-2] + (d(r[i],q[j-1]) + d(r[i],q[j]))/2,
 *   cost_matrix[i-1][j-1] + d(r[i],q[j]),
 *   cost_matrix[i-2][j-1] + d(r[i-1],q[j])  + d(r[i],q[j]),
 *   cost_matrix[i-3][j-1] + d(r[i-2],q[j])  + d(r[i-1],q[j]) + d(r[i],q[j])
       )
 * see doc for the dpw
*/
double p1div2_asym(_DP_ARGS) {
    double d00 = dist(ref, i, query, j, ncols);
    double d01 = dist(ref, i, query, j-1, ncols);
    double d02 = dist(ref, i, query, j-2, ncols);
    double d10 = dist(ref, i-1, query, j, ncols);
    double d20 = dist(ref, i-2, query, j, ncols);

    double arr[5] = {
            cost_matrix[idx(i - 1, j - 3, size2)] + (d02 + d01 + d00) / 3.0,
            cost_matrix[idx(i - 1, j - 2, size2)] + (d01 + d00) / 2,
            cost_matrix[idx(i - 1, j - 1, size2)] + d00,
            cost_matrix[idx(i - 2, j - 1, size2)] + d10 + d00,
            cost_matrix[idx(i - 3, j - 1, size2)] + d20 + d10 + d00
    };
    return min_n(arr, 5);

}

/*
 * Function:  p1_sym
 * --------------------
 * Sakoe-Chiba classification p = 1, symmetric step pattern:
 *  min(
 *      cost_matrix[i-1][j-2] + 2d(r[i],q[j-1]) + d(r[i],q[j]),
 *      cost_matrix[i-1][j-1] + 2d(r[i],q[j]),
 *      cost_matrix[i-2][j-1] + 2d(r[i-1],q[j]) + d(r[i],q[j]),
       )
 * see doc for the dpw
*/
double p1_sym(_DP_ARGS) {
    double d00 = dist(ref, i, query, j, ncols);
    double d01 = dist(ref, i, query, j-1, ncols);
    double d10 = dist(ref, i-1, query, j, ncols);

    return min3(cost_matrix[idx(i - 1, j - 2, size2)] + 2 * d01 + d00,
                cost_matrix[idx(i - 1, j - 1, size2)] + 2 * d00,
                cost_matrix[idx(i - 2, j - 1, size2)] + 2 * d10 + d00
    );
}

/*
 * Function:  p1_asym
 * --------------------
 * Sakoe-Chiba classification p = 1, asymmetric step pattern:
 *  min(
 *      cost_matrix[i-1][j-2] + (d(r[i],q[j-1]) + d(r[i],q[j]))/2,
 *      cost_matrix[i-1][j-1] + d(r[i],q[j]),
 *      cost_matrix[i-2][j-1] + d(r[i-1],q[j]) + d(r[i],q[j]),
       )
 * see doc for the dpw
*/
double p1_asym(_DP_ARGS) {
    double d00 = dist(ref, i, query, j, ncols);
    double d01 = dist(ref, i, query, j-1, ncols);
    double d10 = dist(ref, i-1, query, j, ncols);

    return min3(cost_matrix[idx(i - 1, j - 2, size2)] + (d01 + d00) / 2.0,
                cost_matrix[idx(i - 1, j - 1, size2)] + d00,
                cost_matrix[idx(i - 2, j - 1, size2)] + d10 + d00
    );
}


/*
 * Function:  p2_sym
 * --------------------
 * Sakoe-Chiba classification p = 2, symmetric step pattern:
 *  min(
 *     cost_matrix[i-2][j-3] + 2d(r[i-1],q[j-2]) +
 *                             2d(r[i],q[j-1])   +
 *                             d(r[i],q[j]),
 *
 *     cost_matrix[i-1][j-1] + 2d(r[i],q[j]),
 *
 *     cost_matrix[i-3][j-2] + 2d(r[i-2],q[j-1]) +
 *                             2d(r[i-1],q[j])   +
 *                             d(r[i],q[j])
 *     )
 * see doc for the dpw
*/
double p2_sym(_DP_ARGS) {
    double d00 = dist(ref, i, query, j, ncols);
    double d01 = dist(ref, i, query, j-1, ncols);
    double d10 = dist(ref, i-1, query, j, ncols);
    double d12 = dist(ref, i-1, query, j-2, ncols);
    double d21 = dist(ref, i-2, query, j-1, ncols);

    return min3(cost_matrix[idx(i - 2, j - 3, size2)] + 2 * d12 + 2 * d01 + d00,
                cost_matrix[idx(i - 1, j - 1, size2)] + 2 * d00,
                cost_matrix[idx(i - 3, j - 2, size2)] + 2 * d21 + 2 * d10 + d00
    );

}

/*
 * Function:  p2_asym
 * --------------------
 * Sakoe-Chiba classification p = 2, asymmetric step pattern:
 *  min(
 *     cost_matrix[i-2][j-3] + 2( d(r[i-1],q[j-2]) +
 *                                d(r[i],q[j-1])   +
 *                                d(r[i],q[j]) ),
 *
 *     cost_matrix[i-1][j-1] + d(r[i],q[j]),
 *
 *     cost_matrix[i-3][j-2] + d(r[i-2],q[j-1]) +
 *                             d(r[i-1],q[j])   +
 *                             d(r[i],q[j])
 *     )
 * see doc for the dpw
*/
double p2_asym(_DP_ARGS) {
    double d00 = dist(ref, i, query, j, ncols);
    double d01 = dist(ref, i, query, j-1, ncols);
    double d10 = dist(ref, i-1, query, j, ncols);
    double d12 = dist(ref, i-1, query, j-2, ncols);
    double d21 = dist(ref, i-2, query, j-1, ncols);

    return min3(cost_matrix[idx(i - 2, j - 3, size2)] + 2.0 * (d12 + d01 + d00) / 3.0,
                cost_matrix[idx(i - 1, j - 1, size2)] + d00,
                cost_matrix[idx(i - 3, j - 2, size2)] + d21 + d10 + d00
    );

}


struct t_item dpwdir(_DP_ARGS) {
    double d = dist(ref, i, query, j, ncols);
    return min3idx(cost_matrix[idx(i, j - 1, size2)] + t_s->weights.a * d,
                   cost_matrix[idx(i - 1, j, size2)] + t_s->weights.b * d,
                   cost_matrix[idx(i - 1, j - 1, size2)] + t_s->weights.c * d);
}

/*
 * Function:  p0_symdir
 * --------------------
 * Sakoe-Chiba classification p = 0, symmetric step pattern:
 * see doc for the p0_sym, dp1
*/
struct t_item p0_symdir(_DP_ARGS) {
    return dp1dir(ref, query, cost_matrix, i, j, ncols, t_s, size2, dist);
}

/*
 * Function:  p0_asymdir
 * --------------------
 * Sakoe-Chiba classification p = 0, asymmetric step pattern:
 * see doc for the p0_asym, dp1
*/
struct t_item p0_asymdir(_DP_ARGS) {
    double d = dist(ref, i, query, j, ncols);
    return min3idx(cost_matrix[idx(i, j - 1, size2)],          //0
                   cost_matrix[idx(i - 1, j, size2)] + d,   //1
                   cost_matrix[idx(i - 1, j - 1, size2)] + d);  //2
}

/*
 * Function:  p1div2_symdir
 * --------------------
 * Sakoe-Chiba classification p = 0.5, symmetric step pattern:
 * see doc for the p1div2_sym, dp1
*/
struct t_item p1div2_symdir(_DP_ARGS) {
    int m = i;
    int n = j;
    double d00 = dist(ref, i, query, j, ncols);
    double d01 = dist(ref, i, query, j-1, ncols);
    double d02 = dist(ref, i, query, j-2, ncols);
    double d10 = dist(ref, i-1, query, j, ncols);
    double d20 = dist(ref, i-2, query, j, ncols);
    double arr[5] = {
            cost_matrix[idx(i - 1, j - 3, size2)] + 2.0 * d02 + d01 + d00,
            cost_matrix[idx(i - 1, j - 2, size2)] + 2.0 * d01 + d00,
            cost_matrix[idx(i - 1, j - 1, size2)] + 2.0 * d00,
            cost_matrix[idx(i - 2, j - 1, size2)] + 2.0 * d10 + d00,
            cost_matrix[idx(i - 3, j - 1, size2)] + 2.0 * d20 + d10 + d00
    };
    return min_nidx(arr, 5);

}

/*
 * Function:  p1div2_asymdir
 * --------------------
 * Sakoe-Chiba classification p = 0.5, asymmetric step pattern:
 * see doc for the p1div2_asym, dp1
*/
struct t_item p1div2_asymdir(_DP_ARGS) {
    double d00 = dist(ref, i, query, j, ncols);
    double d01 = dist(ref, i, query, j-1, ncols);
    double d02 = dist(ref, i, query, j-2, ncols);
    double d10 = dist(ref, i-1, query, j, ncols);
    double d20 = dist(ref, i-2, query, j, ncols);
    double arr[5] = {
            cost_matrix[idx(i - 1, j - 3, size2)] + (d02 + d01 + d00) / 3.0,
            cost_matrix[idx(i - 1, j - 2, size2)] + (d01 + d00) / 2,
            cost_matrix[idx(i - 1, j - 1, size2)] + d00,
            cost_matrix[idx(i - 2, j - 1, size2)] + d10 + d00,
            cost_matrix[idx(i - 3, j - 1, size2)] + d20 + d10 + d00
    };
    return min_nidx(arr, 5);

}

/*
 * Function:  p1_symdir
 * --------------------
 * Sakoe-Chiba classification p = 1, symmetric step pattern:
 * see doc for the p1_sym, dp1
*/
struct t_item p1_symdir(_DP_ARGS) {
    double d00 = dist(ref, i, query, j, ncols);
    double d01 = dist(ref, i, query, j-1, ncols);
    double d10 = dist(ref, i-1, query, j, ncols);

    return min3idx(cost_matrix[idx(i - 1, j - 2, size2)] + 2 * d01 + d00,
                   cost_matrix[idx(i - 1, j - 1, size2)] + 2 * d00,
                   cost_matrix[idx(i - 2, j - 1, size2)] + 2 * d10 + d00
    );
}

/*
 * Function:  p1_asymdir
 * --------------------
 * Sakoe-Chiba classification p = 1, asymmetric step pattern:
 * see doc for the p1_asym, dp1
*/
struct t_item p1_asymdir(_DP_ARGS) {
    double d00 = dist(ref, i, query, j, ncols);
    double d01 = dist(ref, i, query, j-1, ncols);
    double d10 = dist(ref, i-1, query, j, ncols);

    return min3idx(cost_matrix[idx(i - 1, j - 2, size2)] + (d01 + d00) / 2.0,
                   cost_matrix[idx(i - 1, j - 1, size2)] + d00,
                   cost_matrix[idx(i - 2, j - 1, size2)] + d10 + d00
    );
}

/*
 * Function:  p2_symdir
 * --------------------
 * Sakoe-Chiba classification p = 2, symmetric step pattern:
 * see doc for the p2_sym, dp1
*/
struct t_item p2_symdir(_DP_ARGS) {
    double d00 = dist(ref, i, query, j, ncols);
    double d01 = dist(ref, i, query, j-1, ncols);
    double d10 = dist(ref, i-1, query, j, ncols);
    double d12 = dist(ref, i-1, query, j-2, ncols);
    double d21 = dist(ref, i-2, query, j-1, ncols);

    return min3idx(cost_matrix[idx(i - 2, j - 3, size2)] + 2 * d12 + 2 * d01 + d00,
                   cost_matrix[idx(i - 1, j - 1, size2)] + 2 * d00,
                   cost_matrix[idx(i - 3, j - 2, size2)] + 2 * d21 + 2 * d10 + d00
    );

}

/*
 * Function:  p2_asymdir
 * --------------------
 * Sakoe-Chiba classification p = 2, asymmetric step pattern:
 * see doc for the p2_asym, dp1
*/
struct t_item p2_asymdir(_DP_ARGS) {
    double d00 = dist(ref, i, query, j, ncols);
    double d01 = dist(ref, i, query, j-1, ncols);
    double d10 = dist(ref, i-1, query, j, ncols);
    double d12 = dist(ref, i-1, query, j-2, ncols);
    double d21 = dist(ref, i-2, query, j-1, ncols);
    return min3idx(cost_matrix[idx(i - 2, j - 3, size2)] + 2.0 * (d12 + d01 + d00) / 3.0,
                   cost_matrix[idx(i - 1, j - 1, size2)] + d00,
                   cost_matrix[idx(i - 3, j - 2, size2)] + d21 + d10 + d00
    );

}

/*
 * Function:  scband
 * --------------------
 *  Sakoe-Chiba band global constraint
 *  NOTE: This function is redundant at the moment(scband is unnecessary, when
 *  there exitsts palival window + there is a faster for-loop implementation)
 *  int i:    row index
 *  int j:    column index
 *  double r: width of the window
 *  double I: length of reference sequence, len_ref, fake arg for this func
 *  double J: length of query sequence, len_query, fake arg for this func
 *  returns:  bool value, if cost[i][j] is legal or not
 */
bool scband(int i, int j, double r, double I, double J) {
    return fabs(i - j) < r;
}

/*
 * Function: palival
 * --------------------
 *  Palival global constraint, it is similar to scband, but but adapts to the
 *  length of sequences
 *  NOTE: This function is redundant at the moment, there is a faster
 *  for-loop implementation.
 *  int i:    row index
 *  int j:    column index
 *  double r: width of the window (abolute value)
 *  double I: length of reference sequence, len_ref
 *  double J: length of query sequence, len_queryc
 *  returns: bool value, if cost[i][j] is legal or not
 */
bool palival(int i, int j, double r, double I, double J) {
    return fabs(i * J / I - j) < r;
}

/*
 * Function: palival_mod
 * --------------------
 *  Palival global constraint, it is similar to scband, but but adapts to the
 *  length of sequences. For the difference palival and palival_mod see the args
 *  NOTE: This function is redundant at the moment, there is a faster
 *  for-loop implementation.
 *  int i:    row index
 *  int j:    column index
 *  double r: width of the window(fraction of the reference sequence length)
 *  double I: length of reference sequence, len_ref
 *  double J: length of query sequence, len_queryc
 *  returns: bool value, if cost[i][j] is legal or not
 */
bool palival_mod(int i, int j, double p, double I, double J) {
    double k = I / J;
    double r = p * I;
    return i > k * j - r && i < k * j + r;
}

/*
 * Function:  itakura
 * --------------------
 *  Itakura global constraints
 *  int i:    row index
 *  int j:    column index
 *  double k: k = 1/s = I/J = len_ref/len_query
 *  double I: length of reference sequence, len_ref
 *  double J: length of query sequence, len_query
 *  returns: bool value, if cost[i][j] is legal or not
 */
bool itakura(int i, int j, double k, double I, double J) {
    return j < 2 * i &&
           i <= 2 * j &&
           i >= I - 1 - 2 * (J - j) &&
           j > J - 1 - 2 * (I - i);

}

/*
 * Function: nowindow
 * --------------------
 * This function is for testing only
 * returns: always true(1)
 */
bool nowindow(int i, int j, double k, double I, double J) {
    return true;
}

/* ------------------------ distance function ------------------------------- */

/*
 * Function:  manhattan
 * --------------------
 *  Euclidean distance
 *  double a:  1 or 2 dimensional point
 *  double b:  1 or 2 dimensional point
 *  int ncols: number of columns of the second dimension (1 means 1 dimensional array)
 *  returns: double, manhattan distance between two dimensional points
 */
double manhattan(double *a, int i, double *b, int j, int ncols) {
    double distance = 0.0;
    double * a_start = a + i*ncols;
    double * b_start = b + j*ncols;
    int k;
    for (k = 0; k < ncols; k++) {
        distance += fabs(a_start[k] - b_start[k]);
    }
    return distance;
}

/*
 * Function:  euclid
 * --------------------
 *  Euclidean distance helper
 *  double a:  1 or 2 dimensional point
 *  double b:  1 or 2 dimensional point
 *  int ncols: number of columns of the second dimension (1 means 1 dimensional array)
 *  returns: double, (a-b)^2
 */
double euclid(double *a, int i, double *b, int j, int ncols) {
    double distance = 0.0;
    double * a_start = a + i*ncols;
    double * b_start = b + j*ncols;
    int k;
    for (k = 0; k < ncols; k++) {
        distance += pow((a_start[k] - b_start[k]), 2);
    }
    return distance;
}



/* ------------------------------- interface -------------------------------- */
/*
 * Function:  choose_dist
 * --------------------
 *  Chooses right distance function(euclid or euclid_squared)
 *  int dist_type: dtw_settings.dist_type [_MANHATTAN, _EUCLID, _EUCLID_SQUARED]
 *  returns: dist_fptr, pointer to a distance function
 */
dist_fptr choose_dist(int dist_type) {
    if (dist_type == _EUCLID || dist_type == _EUCLID_SQUARED) return &euclid;
    else if (dist_type == _MANHATTAN) return &manhattan;
    return &manhattan;
}

/*
 * Function:  choose_dp
 * --------------------
 *  Chooses right step function(without traceback)
 *  int dp_type: dtw_settings.dp_type, step function type
 *  returns: dist_fptr, pointer to a step function without traceback
 */
dp_fptr choose_dp(int dp_type) {

    if (dp_type == _DP1) return &dp1;
    else if (dp_type == _DP2) return &dp2;
    else if (dp_type == _DP3) return &dp3;
    else if (dp_type == _DPW) return &dpw;
    else if (dp_type == _SCP0SYM) return &p0_sym;
    else if (dp_type == _SCP0ASYM) return &p0_asym;
    else if (dp_type == _SCP1DIV2SYM) return &p1div2_sym;
    else if (dp_type == _SCP1DIV2ASYM) return &p1div2_asym;
    else if (dp_type == _SCP1SYM) return &p1_sym;
    else if (dp_type == _SCP1ASYM) return &p1_asym;
    else if (dp_type == _SCP2SYM) return &p2_sym;
    else if (dp_type == _SCP2ASYM) return &p2_asym;
    return &dp2;

}

/*
 * Function:  choose_dpdir
 * --------------------
 *  Chooses right step function(with traceback)
 *  int dp_type: dtw_settings.dp_type, step function type
 *  returns: dpdir_fptr, pointer to a step function with traceback
 */
dpdir_fptr choose_dpdir(int dp_type) {
    if (dp_type == _DP1) return &dp1dir;
    else if (dp_type == _DP2) return &dp2dir;
    else if (dp_type == _DP3) return &dp3dir;
    else if (dp_type == _DPW) return &dpwdir;
    else if (dp_type == _SCP0SYM) return &p0_symdir;
    else if (dp_type == _SCP0ASYM) return &p0_asymdir;
    else if (dp_type == _SCP1DIV2SYM) return &p1div2_symdir;
    else if (dp_type == _SCP1DIV2ASYM) return &p1div2_asymdir;
    else if (dp_type == _SCP1SYM) return &p1_symdir;
    else if (dp_type == _SCP1ASYM) return &p1_asymdir;
    else if (dp_type == _SCP2SYM) return &p2_symdir;
    else if (dp_type == _SCP2ASYM) return &p2_asymdir;
    return &dp2dir;
}

/*
 * Function:  choose_window
 * --------------------
 *  NOTE: function is partly redundant at the moment, it always returns  itakura
 *  Chooses right window function
 *  int dp_type: dtw_settings.win , step function type
 *  returns: dpdir_fptr, pointer to a step function with traceback
 */
window_fptr choose_window(struct t_dtw_settings *dtw_settings) {

    if (dtw_settings->window_type == _SCBAND) { return &scband; }
    else if (dtw_settings->window_type == _PALIVAL) { return &palival; }
    else if (dtw_settings->window_type == _ITAKURA) { return &itakura; }
    else if (dtw_settings->window_type == _PALIVAL_MOD) { return &palival_mod; }
    return &nowindow;
}

/*
 * Function:  choose_window_param
 * --------------------
 *  Computes and return parameter for the window function
 *  struct t_dtw_settings *dtw_settings: structure with dtw settings
 *  int len_ref: length of the reference sequence
 *  int len_query: length of the query sequence
 *  returns: double p, window parameter
 */
double choose_window_param(struct t_dtw_settings *dtw_settings,
                           int len_ref,
                           int len_query) {
    double p = 0;
    if (dtw_settings->window_type == _ITAKURA)
        p = len_ref / (double) len_query;
    else
        p = dtw_settings->window_param;
    return p;

}

/*
 * Function:  choose_path_pattern
 * --------------------
 *  Chooses right path_patterns(2d array) based on the current step function
 *  struct t_dtw_settings *dtw_settings: structure with dtw settings
 *  returns: const int[7][11], path_pattern
 */
const int (*choose_path_pattern(struct t_dtw_settings dtw_settings))[11] {
    if (dtw_settings.dp_type == _DP1 ||
        dtw_settings.dp_type == _DP2 ||
        dtw_settings.dp_type == _SCP0ASYM ||
        dtw_settings.dp_type == _SCP0ASYM)
        return dp2_path_pattern;
    else if (dtw_settings.dp_type == _DP3)
        return dp1_path_pattern;
    else if (dtw_settings.dp_type == _SCP1DIV2SYM ||
             dtw_settings.dp_type == _SCP1DIV2ASYM)
        return p1div2_path_pattern;
    else if (dtw_settings.dp_type == _SCP1SYM ||
             dtw_settings.dp_type == _SCP1ASYM)
        return p1_path_pattern;
    else if (dtw_settings.dp_type == _SCP2SYM ||
             dtw_settings.dp_type == _SCP2ASYM)
        return p2_path_pattern;
    else
        return dp2_path_pattern;
}

/*
 * Function:  extra_size
 * --------------------
 *  Computes extra size for the cost matrix (dtw_settings.offset)
 *  int dp_type: dtw_settings.dp_type, step pattern type
 *  returns: int (1,2 or 3)
 */
int extra_size(int dp_type) {
    if (dp_type == _DP1 ||
        dp_type == _DP2 ||
        dp_type == _DP3 ||
        dp_type == _SCP0SYM ||
        dp_type == _SCP0ASYM)
        return 1;
    else if (dp_type == _SCP1SYM ||
             dp_type == _SCP1ASYM)
        return 2;
    else
        return 3;

}

/*
 * Function:  create_path_from_pattern
 * --------------------
 *  Compute full path from path_points and path_pattern. This function is
 *  necessary for the multistep step functions (like p1div2_sym), and also
 *  to transform "directions" to coordinates(matrix indices).
 *  const int pattern[6][11] : path pattern (see cdtw.h)
 *  int len_ref: length of the reference sequence
 *  int len_query: length of the query sequence
 *  int* path_points: path points as array of directions, where direction is
 *                    is idx returned by dp_dir functions.
 *  int  path_points_count: len of path points = len_ref + len_query
 *  struct t_path_element* path: path array
 *  returns: int, true length of the constructed path.
 */
int create_path_from_pattern(const int pattern[6][11],
                             int len_ref,
                             int len_query,
                             int *path_points,
                             int path_points_count,
                             struct t_path_element *path
) {
    int path_idx = 1;
    int i = 0;
    int j = 0;
    path[0].i = len_ref - 1;
    path[0].j = len_query - 1;

    for (; i < path_points_count; i++) {
        int path_idx_tmp = path_idx;
        for (j = 1; j < 2 * pattern[path_points[i] + 1][0] + 1; j += 2) {
            path[path_idx].i = path[path_idx_tmp - 1].i +
                               pattern[path_points[i] + 1][j];
            path[path_idx].j = path[path_idx_tmp - 1].j +
                               pattern[path_points[i] + 1][j + 1];
            path_idx++;
        }
    }
    return path_idx;
}


/*
 * Function:  direct_matrix_to_path_points
 * --------------------
 * Computes path(direction array) from directions matrix
 * int* dir_matrix: directions matrix
 * int *path_points: direction array(path points)
 * int len_ref: length of the reference sequence
 * int len_query: length of the query sequence
 * const int pattern[6][11]:path pattern
 * returns: path_points length
 */
int direct_matrix_to_path_points(int *dir_matrix,
                                 int *path_points,
                                 int len_ref,
                                 int len_query,
                                 const int pattern[6][11]) {
    long tin_idx = len_ref * len_query - 1;
    int tout_idx = 0;

    int i = 0;
    int j = 0;


    path_points[tout_idx] = dir_matrix[tin_idx];
    for (; tin_idx >= 1;) {
        i = pattern[0][2 * dir_matrix[tin_idx] + 1];
        j = pattern[0][2 * dir_matrix[tin_idx] + 2];
        tin_idx += i * len_query + j;
        tout_idx++;
        if (tin_idx < 1)
            break;
        path_points[tout_idx] = dir_matrix[tin_idx];
    }
    return tout_idx;

}

/*
 * Function:  cdtwnopath
 * --------------------
 * Dynamic Time Warping algorithm(without traceback)
 * double* ref:         reference sequence
 * double* query:       query sequence
 * int len_ref:         length of the reference sequence
 * int len_query:       length of the query sequence
 * int ncols:           number of columns of the second dimension (1 means 1 dimensional array)
 * dist_fptr dist:      pointer to a distance function
 * dp_fptr dp:          pointer to a step function
 * window_fptr window:  pointer to a window function
 * double p:            window parameter
 * double *cost_matrix: cost matrix
 * struct t_dtw_settings dtw_settings: structure with dtw settings
 * returns: doube, distance between reference and query sequences
 */
double cdtwnopath(double *ref,
                  double *query,
                  int len_ref,
                  int len_query,
                  int ncols,
                  dist_fptr dist,
                  dp_fptr dp,
                  window_fptr window,
                  double p, //window param
                  double *cost_matrix,
                  struct t_dtw_settings dtw_settings) {
    int off = dtw_settings.offset;
    /*memory was already allocated*/
    /*extending matrix*/
    int M = len_ref + off;
    int N = len_query + off;
    int i = 0;
    int j = 0;
    double w = 0;
    double s = 0;
    bool fast_glob = (dtw_settings.window_type == _PALIVAL ||
                      dtw_settings.window_type == _PALIVAL_MOD);
    /*no window or fast window case*/
    if (fast_glob || dtw_settings.window_type == 0) {
        if (fast_glob) {
            if (dtw_settings.window_type == _PALIVAL_MOD)
                w = dtw_settings.window_param * (double) len_ref;
            else
                w = dtw_settings.window_param;
            s = len_query / (double) len_ref;
        }
        else {
            w = INFINITY;
            s = 1;
        }

        cost_matrix[idx(off, off, N)] = dist(ref, off, query, off, ncols);
        for (j = max2(off + 1, _round(s * (off - w))); j < min2(N, _round(s * (off + w) + 1)); j++) {
            cost_matrix[idx(off, j, N)] = dp(ref, query, cost_matrix, off, j, ncols, &dtw_settings, N, dist);
        }

        for (i = off + 1; i < M; i++) {
            for (j = max2(off, _round(s * (i - w))); j < min2(N, _round(s * (i + w) + 1)); j++)
                cost_matrix[idx(i, j, N)] = dp(ref, query, cost_matrix, i, j, ncols, &dtw_settings, N, dist);
        }

    }
        /*slow window case*/
    else {
        cost_matrix[idx(off, off, N)] = dist(ref, off, query, off, ncols);
        for (j = off + 1; j < N; j++) {
            if (window(off, j, p, len_ref, len_query))
                cost_matrix[idx(off, j, N)] = dp(ref, query, cost_matrix, off, j, ncols, &dtw_settings, N, dist);
        }
        for (i = off + 1; i < M; i++) {
            for (j = off; j < N; j++) {
                if (window(i, j, p, len_ref, len_query))
                    cost_matrix[idx(i, j, N)] = dp(ref, query, cost_matrix, i, j, ncols, &dtw_settings, N, dist);
            }
        }
    }

    return cost_matrix[idx(M - 1, N - 1, N)];

}

/*
 * Function:  cdtwpath
 * --------------------
 * Dynamic Time Warping algorithm(with traceback)
 * double* ref:         reference sequence
 * double* query:       query sequence
 * int len_ref:         length of the reference sequence
 * int len_query:       length of the query sequence
 * int ncols:           number of columns of the second dimension (1 means 1 dimensional array)
 * dist_fptr dist:      pointer to a distance function
 * dpdir_fptr dp_dir:   pointer to a step function
 * window_fptr window:  pointer to a window function
 * double p:            window parameter
 * double *cost_matrix: cost matrix
 * int *dir_matrix:     direction matrix
 * struct t_dtw_settings dtw_settings: structure with dtw settings
 * returns: doube, distance between reference and query sequences
 */
double cdtwpath(double *ref,
                double *query,
                int len_ref,
                int len_query,
                int ncols,
                dist_fptr dist,
                dpdir_fptr dp_dir,
                window_fptr window,
                double p, //window param
                double *cost_matrix,
                int *dir_matrix,
                struct t_dtw_settings dtw_settings) {
    int off = dtw_settings.offset;

    struct t_item item = {0, 0};
    /*extending matrix*/
    int M = len_ref + off;
    int N = len_query + off;

    int i = 0;
    int j = 0;
    double w = 0;
    double s = 0;
    bool fast_glob = (dtw_settings.window_type == _PALIVAL ||
                      dtw_settings.window_type == _PALIVAL_MOD);
    /*no window or fast window case*/
    if (fast_glob || dtw_settings.window_type == 0) {
        if (fast_glob) {
            if (dtw_settings.window_type == _PALIVAL_MOD)
                w = dtw_settings.window_param * (double) len_ref;
            else
                w = dtw_settings.window_param;
            s = len_query / (double) len_ref;
        }
        else {
            w = INFINITY;
            s = 1;
        }

        cost_matrix[idx(off, off, N)] = dist(ref, off, query, off, ncols);
        for (j = max2(off + 1, _round(s * (off - w))); j < min2(N, _round(s * (off + w) + 1)); j++) {
            item = dp_dir(ref, query, cost_matrix, off, j, ncols, &dtw_settings, N, dist);
            cost_matrix[idx(off, j, N)] = item.val;
            dir_matrix[idx(0, j - off, N - off)] = item.idx;
        }

        for (i = off + 1; i < M; i++) {
            for (j = max2(off, _round(s * (i - w))); j < min2(N, _round(s * (i + w) + 1)); j++) {
                item = dp_dir(ref, query, cost_matrix, i, j, ncols, &dtw_settings, N, dist);
                cost_matrix[idx(i, j, N)] = item.val;
                dir_matrix[idx(i - off, j - off, N - off)] = item.idx;
            }
        }
    }
        /*slow window case*/
    else {
        cost_matrix[idx(off, off, N)] = dist(ref, off, query, off, ncols);
        for (j = off + 1; j < N; j++) {
            if (window(off, j, p, len_ref, len_query)) {
                item = dp_dir(ref, query, cost_matrix, off, j, ncols, &dtw_settings, N, dist);
                cost_matrix[idx(off, j, N)] = item.val;
                dir_matrix[idx(0, j - off, N - off)] = item.idx;
            }
        }
        for (i = off + 1; i < M; i++) {
            for (j = off; j < N; j++) {
                if (window(i, j, p, len_ref, len_query)) {
                    item = dp_dir(ref, query, cost_matrix, i, j, ncols, &dtw_settings, N, dist);
                    cost_matrix[idx(i, j, N)] = item.val;
                    dir_matrix[idx(i - off, j - off, N - off)] = item.idx;
                }
            }
        }
    }
    return cost_matrix[idx(M - 1, N - 1, N)];

}

/*
 * Function:  fill_matrix
 * --------------------
 * Prepares cost matrix. Set extra rows and columns to INFINITY. Set all
                         matrix to INFINITY, if there is a window

 * double *matrix:      pointer to cost matrix
 * double* query:       query sequence
 * int len_ref:         length of the reference sequence
 * int len_query:       length of the query sequence
 * struct t_dtw_settings dtw_settings: structure with dtw settings
 * returns: void
 */
void fill_matrix(double *matrix,
                 int len_ref,
                 int len_query,
                 struct t_dtw_settings dtw_settings) {
    /*
    http://www.cplusplus.com/reference/cstring/memset/
    or for 1 dimension {INFINITY}
    */
    int M = len_ref + dtw_settings.offset;
    int N = len_query + dtw_settings.offset;
    int i = 0;
    int j = 0;
    /*if there is a window or complicated step pattern*/
    if (dtw_settings.window_type != 0 || (dtw_settings.dp_type != _DP1 &&
                                          dtw_settings.dp_type != _DP2 &&
                                          dtw_settings.dp_type != _DP3 &&
                                          dtw_settings.dp_type != _SCP0SYM)
            ) /* + _SCP0ASYM??*/
    {
        for (i = 0; i < M; i++)
            for (j = 0; j < N; j++)
                matrix[idx(i, j, N)] = INFINITY;
    }
    else {
        for (i = 0; i < M; i++)
            for (j = 0; j < dtw_settings.offset; j++)
                matrix[idx(i, j, N)] = INFINITY;

        for (i = 0; i < dtw_settings.offset; i++)
            for (j = 0; j < N; j++)
                matrix[idx(i, j, N)] = INFINITY;
    }

}

/*
 * Function:  cdtw
 * --------------------
 * Dynamic Time Warping
 * This fuction is main entry for the cdtw
 * double* ref:                  reference sequence
 * double* query:               query sequence
 * int len_ref:                 length of the reference sequence
 * int len_query:               length of the query sequence
 * int ncols:                   number of columns of the second dimension (1 means 1 dimensional array)
 * double *cost_matrix:         cost matrix, pointer to an allocated memory for
                                the 2 dimensional matrix. (with extra size)
 * struct t_path_element* path: warping path, pointer to an allocated array,
                                size is maximum possible path length:
                                len_ref + len_query
 * int true_path_len:           length of the computed warping path
 * struct t_dtw_settings dtw_settings: structure with dtw settings
 * returns: doube, distance between reference and query sequences
 */
double cdtw(double *ref,
            double *query,
            int len_ref,
            int len_query,
            int ncols,
            double *cost_matrix,
            struct t_path_element *path,
            int *true_path_len,
            struct t_dtw_settings dtw_settings) {
    /*init distance*/
    double distance = 0;
    /*init window function param*/
    double p = 0;
    dp_fptr dp = NULL;
    /*init pointer to step function*/
    dpdir_fptr dp_dir = NULL;

    /*pointer to window function*/
    window_fptr window = choose_window(&dtw_settings);
    /*pointer to distance function*/
    dist_fptr dist = choose_dist(dtw_settings.dist_type);

    int *dir_matrix = NULL;

    int *path_points = NULL;

    int path_points_count;

    /*assign step function*/
    if (dtw_settings.compute_path)
        dp_dir = choose_dpdir(dtw_settings.dp_type);
    else
        dp = choose_dp(dtw_settings.dp_type);


    /*assign window parameter*/
    p = choose_window_param(&dtw_settings,
                            len_ref,
                            len_query);


    /*lets go !*/

    /*prepare cost matrix*/
    fill_matrix(cost_matrix, len_ref, len_query, dtw_settings);
    /*dtw withou traceback case(only cost_matrix and distance)*/
    if (dtw_settings.compute_path == false) {

        distance = cdtwnopath(ref,
                              query,
                              len_ref,
                              len_query,
                              ncols,
                              dist,
                              dp,
                              window,
                              p,
                              cost_matrix,
                              dtw_settings);
    }
        /*dtw with traceback case*/
    else {
        /*allocate direction matrix*/
        dir_matrix = (int *) malloc(sizeof(int) * len_ref * len_query);
        /*call cdtwpath, computes distance, cost matrix and direction matrix*/
        distance = cdtwpath(ref,
                            query,
                            len_ref,
                            len_query,
                            ncols,
                            dist,
                            dp_dir,
                            window,
                            p,
                            cost_matrix,
                            dir_matrix,
                            dtw_settings);

        /*if distance is INFINITY there is not any path*/
        if (distance == INFINITY) {
            *true_path_len = 0;
            return INFINITY;
        }
        /*allocate path points(direction array)*/
        path_points = (int *) calloc(len_ref + len_query, sizeof(int));


        /*compute path(directions array) from direction matrix*/
        path_points_count = direct_matrix_to_path_points(
                dir_matrix,
                path_points,
                len_ref,
                len_query,
                choose_path_pattern(dtw_settings));


        /*cleaning*/
        free(dir_matrix);
        /*compute warping path(finnally, as array of the path elements*/
        *true_path_len = create_path_from_pattern(
                choose_path_pattern(dtw_settings),
                len_ref,
                len_query,
                path_points,
                path_points_count,
                path);
        /*cleaning*/
        free(path_points);

    }


    /*euclidian metric case*/
    if (dtw_settings.dist_type == _EUCLID)
        distance = sqrt(distance);


    /*normalization case*/
    if (dtw_settings.norm)
        return distance / (double) (len_ref + len_query);  //TOTO: ADD NORM FACTORS
        /*there is no normalization*/
    else
        return distance;

}


void print_matrix(double *matrix, int size1, int size2) {
    int i = 0;
    int j = 0;
    for (i = 0; i < size1; i++) {
        for (j = 0; j < size2; j++)
            printf("%8.2f", (double) matrix[idx(i, j, size2)]);
        printf("\n");
    }
}

void print_intarr(int *arr, int n) {
    int i = 0;
    for (i = 0; i < n; i++)
        printf("%4i", arr[i]);
    printf("\n");
}

void print_floatarr(float *arr, int n) {
    int i = 0;
    for (i = 0; i < n; i++)
        printf("%.2f ", arr[i]);
    printf("\n");
}

void print_doublearr(double *arr, int n) {
    int i = 0;
    for (i = 0; i < n; i++)
        printf("%.2f ", arr[i]);
    printf("\n");
}

//int main() {
//    int i, j;
//
//    int ncols = 2;
//    double r[] = {1, 0, 5, 0, 4, 0, 2, 0};
//    double q[] = {1, 0, 2, 0, 4, 0, 1, 0};
//
////    int ncols = 1;
////    double r[] = {1, 5, 4, 2};
////    double q[] = {1, 2, 4, 1};
//
//    int len_ref = sizeof(r) / sizeof(r[0]) / ncols;
//    int len_query = sizeof(q) / sizeof(q[0]) / ncols;
//    struct t_dtw_settings dtw_settings;
//
//    dtw_settings.compute_path = true;
//    dtw_settings.dist_type = _MANHATTAN;
//    dtw_settings.dp_type = _DP2;
//    dtw_settings.window_type = false;
//    dtw_settings.window_param = 0;
//    dtw_settings.norm = false;
//    dtw_settings.offset = extra_size(dtw_settings.dp_type);
//
//    int offset = dtw_settings.offset * ncols;
//
//    double *expanded_r = malloc(sizeof(double) * (len_ref*ncols + offset));
//    double *expanded_q = malloc(sizeof(double) * (len_query*ncols + offset));
//
//    for (i = 0; i < offset; i++) {
//        expanded_r[i] = expanded_q[i] = 0;
//    }
//
//    for (i = 0; i < len_ref * ncols; i++) {
//        expanded_r[i + offset] = r[i];
//    }
//    for (i = 0; i < len_query * ncols; i++) {
//        expanded_q[i + offset] = q[i];
//    }
//
//    double *cost = (double *) malloc(
//            sizeof(double) * ((len_ref + dtw_settings.offset) * (len_query + dtw_settings.offset)));
//
//
//    struct t_path_element *path = (struct t_path_element *) malloc(
//            sizeof(struct t_path_element) * (len_ref + len_query));
//    int path_len = 0;
//    cdtw(expanded_r, expanded_q, len_ref, len_query, ncols, cost, path, &path_len, dtw_settings);
//    if (cost[14 * 14 - 1] == 4282.0)
//        printf("\nFAIL");
//    print_matrix(cost, len_ref + dtw_settings.offset, len_query + dtw_settings.offset);
//    for (i = 0; i < path_len; i++)
//        printf("i: %d  j: %d\n", path[i].i, path[i].j);
//
//    return 0;
//}
