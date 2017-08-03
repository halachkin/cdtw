#include <math.h>
#include "stdlib.h"
#include "ced.h"

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
 * Function:  max2idx
 * --------------------
 *  Finds the maximum of 2 doubles and its position (0 or 1)

 *  double a: a
 *  double b: b

 *  returns:  struct t_item, t_item = { min(a,b), position }
 */
struct t_item max2idx(double a, double b) /*0, 1*/
{
    struct t_item item;
    if (a >= b) {
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
 *_DP_ARGS(ced.h). Step functions like step_pattern_type and
 * step_pattern_typedir are pretty similar, step_pattern_type are used in
 * computing edit distance without path(without traceback). step_pattern_typedir are
 * used in computing edit distance with path(traceback)

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
 * doubel* args.sigma: maximum difference allowances of all dimension (usable by EDR)
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
    double d = quantise(dist(ref, i, query, j, ncols), args.sigma);
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
    double d = quantise(dist(ref, i, query, j, ncols), args.sigma);
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
    double d = quantise(dist(ref, i, query, j, ncols), args.sigma);
    return (min3(cost_matrix[idx(i, j - 1, size2)],
                 cost_matrix[idx(i - 1, j, size2)],
                 cost_matrix[idx(i - 1, j - 1, size2)]) + d);
}


/*
 * Function:  dp2_edr (step pattern used by EDR)
 * --------------------
 *  Step patern dp2:
 *  min(
 *      cost_matrix[i][j-1]   +   1 # Insert cost = 1
        cost_matrix[i-1][j]   +   1 # Delete cost = 1
        cost_matrix[i-1][j-1] +   d(r[i],q[j]) # substitution cost (0 if similar, 1 otherwise)
       )
 * Source: [L Chen, MT Özsu, V Oria Robust and fast similarity search for moving object trajectories]
*/
double dp2_edr(_DP_ARGS) {
    double d = quantise(dist(ref, i, query, j, ncols), args.sigma);
    return min3(cost_matrix[idx(i, j - 1, size2)] + 1,   //left
                cost_matrix[idx(i - 1, j, size2)] + 1,   //up
                cost_matrix[idx(i - 1, j - 1, size2)] + d);  //diag
}


/*
 * Function:  dp2_erp (step pattern used by ERP)
 * --------------------
 *  Step patern dp2:
 *  min(
 *      cost_matrix[i][j-1]   +   d(r[i], g) # Deletion ( of a args.gap) cost = distance between the element & the args.gap
        cost_matrix[i-1][j]   +   d(g, q[j]) # Insertion cost = distance between the inserted element and the args.gap
        cost_matrix[i-1][j-1] +   d(r[i],q[j]) # substitution cost (0 if similar, 1 otherwise)
       )
 * Source: [L. Chen and R. Ng. On the marriage of edit distance and Lp norms. In Proc. VLDB, 2004.]
*/
double dp2_erp(_DP_ARGS) {
    double del = dist(query, j, args.gap, 0, ncols);
    double ins = dist(ref, i, args.gap, 0, ncols);

    /* The quantisation for ERP should be no-op */
    double d = quantise(dist(ref, i, query, j, ncols), args.sigma);
    return min3(cost_matrix[idx(i, j - 1, size2)] + ins,   //left
                cost_matrix[idx(i - 1, j, size2)] + del,   //up
                cost_matrix[idx(i - 1, j - 1, size2)] + d);  //diag
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
    double d = quantise(dist(ref, i, query, j, ncols), args.sigma);
    return (min2(cost_matrix[idx(i, j - 1, size2)],
                 cost_matrix[idx(i - 1, j, size2)]) + d);
}

/*
 * Function:  dp3_lcss (for LCSS)
 * --------------------
 *  Step patern:
      If the ref[i] and query[j] match, then increase the previous edge (cost_matrix[i-1, j-1]) by 1
      Otherwise, take the max of the up (cost_matrix[i, j-1]) and left (cost_matrix[i-1, j]) elements

 * see  M. Vlachos, M. Hadjieleftheriou, D. Gunopulos, and E. Keogh.
 *      "Indexing multi-dimensional time-series with support for multiple distance measures",
 *      In Proc. SIGKDD, 2003
*/
double dp3_lcss(_DP_ARGS) {
    /* The quantisation function should return 1 for matching pair and 0 for non matching pair */
    double d = quantise(dist(ref, i, query, j, ncols), args.sigma);
    if (d == 1) {
        return cost_matrix[idx(i - 1, j - 1, size2)] + 1;
    }
    return (max2(cost_matrix[idx(i, j - 1, size2)],
                 cost_matrix[idx(i - 1, j, size2)]));
}

/*
 * Function:  dp1dir
 * --------------------
 * see doc for the dpw
*/
struct t_item dp1dir(_DP_ARGS) {
    double d = quantise(dist(ref, i, query, j, ncols), args.sigma);
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
    double d = quantise(dist(ref, i, query, j, ncols), args.sigma);
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
    double d = quantise(dist(ref, i, query, j, ncols), args.sigma);
    return min2idx(cost_matrix[idx(i, j - 1, size2)] + d,
                   cost_matrix[idx(i - 1, j, size2)] + d);
}

/*
 * Function:  dp3lcssdir (for LCSS)
 * --------------------
 *  Step patern:
      If the ref[i] and query[j] match, then increase the previous edge (cost_matrix[i-1, j-1]) by 1
      Otherwise, take the max of the up (cost_matrix[i, j-1]) and left (cost_matrix[i-1, j]) elements

 * see  M. Vlachos, M. Hadjieleftheriou, D. Gunopulos, and E. Keogh.
 *      "Indexing multi-dimensional time-series with support for multiple distance measures",
 *      In Proc. SIGKDD, 2003
*/
struct t_item dp3_lcss_dir(_DP_ARGS) {
    /* The quantisation function should return 1 for matching pair and 0 for non matching pair */
    double d = quantise(dist(ref, i, query, j, ncols), args.sigma);
    if (d == 1) {
        struct t_item item;
        item.val = cost_matrix[idx(i - 1, j - 1, size2)] + 1;
        item.idx = 0;
        return item;
    }
    return max2idx(cost_matrix[idx(i, j - 1, size2)],
                   cost_matrix[idx(i - 1, j, size2)]);
}

/**
 * Step function used in EDR:
 *  min(
 *      cost_matrix[i][j-1]   +   1 # Insert cost = 1
        cost_matrix[i-1][j]   +   1 # Delete cost = 1
        cost_matrix[i-1][j-1] +   d(r[i],q[j]) # substitution cost (0 if similar, 1 otherwise)
       )
 */
struct t_item dp2_edr_dir(_DP_ARGS) {
    double d = quantise(dist(ref, i, query, j, ncols), args.sigma);
    return min3idx(cost_matrix[idx(i, j - 1, size2)] + 1,   //left
                   cost_matrix[idx(i - 1, j, size2)] + 1,   //up
                   cost_matrix[idx(i - 1, j - 1, size2)] + d);  //diag
}

/*
 * Function:  step pattern used by ERP
 * --------------------
 *  Step patern dp2:
 *  min(
 *      cost_matrix[i][j-1]   +   d(r[i], g) # Insert (a args.gap) cost = distance between the element & the args.gap
        cost_matrix[i-1][j]   +   d(g, q[j]) # Delete cost = distance between the deleted element and the args.gap
        cost_matrix[i-1][j-1] +   d(r[i],q[j]) # substitution cost (0 if similar, 1 otherwise)
       )
 * Source: [L. Chen and R. Ng. On the marriage of edit distance and Lp norms. In Proc. VLDB, 2004.]
*/
struct t_item dp2_erp_dir(_DP_ARGS) {
    double del = dist(query, j, args.gap, 0, ncols);
    double ins = dist(ref, i, args.gap, 0, ncols);

    /* The quantisation for ERP should be no-op */
    double d = quantise(dist(ref, i, query, j, ncols), args.sigma);
    return min3idx(cost_matrix[idx(i, j - 1, size2)] + ins,   //left
                   cost_matrix[idx(i - 1, j, size2)] + del,   //up
                   cost_matrix[idx(i - 1, j - 1, size2)] + d);  //diag
}

/*
 * Function:  p0sym
 * --------------------
 * Sakoe-Chiba classification p = 0, symmetric step pattern
 * This function is alias for the dp1
 * see doc for the dpw
*/
double p0_sym(_DP_ARGS) {
    return dp1(ref, query, args, cost_matrix, i, j, ncols, t_s, size2, dist, quantise);
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
    double d = quantise(dist(ref, i, query, j, ncols), args.sigma);
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
    double d00 = quantise(dist(ref, i, query, j, ncols), args.sigma);
    double d01 = quantise(dist(ref, i, query, j - 1, ncols), args.sigma);
    double d02 = quantise(dist(ref, i, query, j - 2, ncols), args.sigma);
    double d10 = quantise(dist(ref, i - 1, query, j, ncols), args.sigma);
    double d20 = quantise(dist(ref, i - 2, query, j, ncols), args.sigma);
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
    double d00 = quantise(dist(ref, i, query, j, ncols), args.sigma);
    double d01 = quantise(dist(ref, i, query, j - 1, ncols), args.sigma);
    double d02 = quantise(dist(ref, i, query, j - 2, ncols), args.sigma);
    double d10 = quantise(dist(ref, i - 1, query, j, ncols), args.sigma);
    double d20 = quantise(dist(ref, i - 2, query, j, ncols), args.sigma);

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
    double d00 = quantise(dist(ref, i, query, j, ncols), args.sigma);
    double d01 = quantise(dist(ref, i, query, j - 1, ncols), args.sigma);
    double d10 = quantise(dist(ref, i - 1, query, j, ncols), args.sigma);

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
    double d00 = quantise(dist(ref, i, query, j, ncols), args.sigma);
    double d01 = quantise(dist(ref, i, query, j - 1, ncols), args.sigma);
    double d10 = quantise(dist(ref, i - 1, query, j, ncols), args.sigma);

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
    double d00 = quantise(dist(ref, i, query, j, ncols), args.sigma);
    double d01 = quantise(dist(ref, i, query, j - 1, ncols), args.sigma);
    double d10 = quantise(dist(ref, i - 1, query, j, ncols), args.sigma);
    double d12 = quantise(dist(ref, i - 1, query, j - 2, ncols), args.sigma);
    double d21 = quantise(dist(ref, i - 2, query, j - 1, ncols), args.sigma);

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
    double d00 = quantise(dist(ref, i, query, j, ncols), args.sigma);
    double d01 = quantise(dist(ref, i, query, j - 1, ncols), args.sigma);
    double d10 = quantise(dist(ref, i - 1, query, j, ncols), args.sigma);
    double d12 = quantise(dist(ref, i - 1, query, j - 2, ncols), args.sigma);
    double d21 = quantise(dist(ref, i - 2, query, j - 1, ncols), args.sigma);

    return min3(cost_matrix[idx(i - 2, j - 3, size2)] + 2.0 * (d12 + d01 + d00) / 3.0,
                cost_matrix[idx(i - 1, j - 1, size2)] + d00,
                cost_matrix[idx(i - 3, j - 2, size2)] + d21 + d10 + d00
    );

}


struct t_item dpwdir(_DP_ARGS) {
    double d = quantise(dist(ref, i, query, j, ncols), args.sigma);
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
    return dp1dir(ref, query, args, cost_matrix, i, j, ncols, t_s, size2, dist, quantise);
}

/*
 * Function:  p0_asymdir
 * --------------------
 * Sakoe-Chiba classification p = 0, asymmetric step pattern:
 * see doc for the p0_asym, dp1
*/
struct t_item p0_asymdir(_DP_ARGS) {
    double d = quantise(dist(ref, i, query, j, ncols), args.sigma);
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
    double d00 = quantise(dist(ref, i, query, j, ncols), args.sigma);
    double d01 = quantise(dist(ref, i, query, j - 1, ncols), args.sigma);
    double d02 = quantise(dist(ref, i, query, j - 2, ncols), args.sigma);
    double d10 = quantise(dist(ref, i - 1, query, j, ncols), args.sigma);
    double d20 = quantise(dist(ref, i - 2, query, j, ncols), args.sigma);
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
    double d00 = quantise(dist(ref, i, query, j, ncols), args.sigma);
    double d01 = quantise(dist(ref, i, query, j - 1, ncols), args.sigma);
    double d02 = quantise(dist(ref, i, query, j - 2, ncols), args.sigma);
    double d10 = quantise(dist(ref, i - 1, query, j, ncols), args.sigma);
    double d20 = quantise(dist(ref, i - 2, query, j, ncols), args.sigma);
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
    double d00 = quantise(dist(ref, i, query, j, ncols), args.sigma);
    double d01 = quantise(dist(ref, i, query, j - 1, ncols), args.sigma);
    double d10 = quantise(dist(ref, i - 1, query, j, ncols), args.sigma);

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
    double d00 = quantise(dist(ref, i, query, j, ncols), args.sigma);
    double d01 = quantise(dist(ref, i, query, j - 1, ncols), args.sigma);
    double d10 = quantise(dist(ref, i - 1, query, j, ncols), args.sigma);

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
    double d00 = quantise(dist(ref, i, query, j, ncols), args.sigma);
    double d01 = quantise(dist(ref, i, query, j - 1, ncols), args.sigma);
    double d10 = quantise(dist(ref, i - 1, query, j, ncols), args.sigma);
    double d12 = quantise(dist(ref, i - 1, query, j - 2, ncols), args.sigma);
    double d21 = quantise(dist(ref, i - 2, query, j - 1, ncols), args.sigma);

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
    double d00 = quantise(dist(ref, i, query, j, ncols), args.sigma);
    double d01 = quantise(dist(ref, i, query, j - 1, ncols), args.sigma);
    double d10 = quantise(dist(ref, i - 1, query, j, ncols), args.sigma);
    double d12 = quantise(dist(ref, i - 1, query, j - 2, ncols), args.sigma);
    double d21 = quantise(dist(ref, i - 2, query, j - 1, ncols), args.sigma);
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

void vector_subtract(double *v1, double *v2, double *result, int n) {
    int i;
    for (i = 0; i < n; i++) {
        result[i] = v1[i] - v2[i];
    }
}

double l1norm(double *a, int n) {
    int i;
    double result = 0;
    for (i = 0; i < n; i++) {
        result += fabs(a[i]);
    }
    return result;
}

double l2norm(double *a, int n) {
    int i;
    double result = 0;
    for (i = 0; i < n; i++) {
        result += a[i] * a[i];
    }
    return result;
}


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
    double *a_start = a + i * ncols;
    double *b_start = b + j * ncols;
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
    double distance = euclid_square(a, i, b, j, ncols);
    return sqrt(distance);
}


/*
 * Function:  euclid squared
 * --------------------
 *  Euclidean distance helper
 *  double a:  1 or 2 dimensional point
 *  double b:  1 or 2 dimensional point
 *  int ncols: number of columns of the second dimension (1 means 1 dimensional array)
 *  returns: double, (a-b)^2
 */
double euclid_square(double *a, int i, double *b, int j, int ncols) {
    double distance = 0.0;
    double *a_start = a + i * ncols;
    double *b_start = b + j * ncols;
    int k;
    for (k = 0; k < ncols; k++) {
        distance += pow((a_start[k] - b_start[k]), 2);
    }
    return distance;
}
//
///**
// * Quantise the value of real distance (d) given the tolerance value (tol)
// * The distance is an Lp-norm_type
// */
//double edr_norm(double* raw_dist, int ncols, double tol, double (*norm_type)(double *, int)){
//    double d = norm_type(raw_dist, ncols);
//    return (d <= tol) ? 0 : 1;
//}
//
///**
// * Quantise the value of the distance between two vectors given the tolerance values (tol)
// * The distance is still a vector, not norm_type. So is the tolerance.
// * The quantised value is 0 (similar) if the distance vector is less than the tolerance vector at all dimension
// */
//double edr_dim(double* raw_dist, int ncols, double * tol){
//    int i;
//    for (i=0; i<ncols; i++) {
//        if (raw_dist[i] > tol[i]) {
//            return 1;
//        }
//    }
//    return 0;
//}

/**
 * Quantise the value of real distance (d) given the tolerance value (s)
 */
double edr(double d, double s) {
    return (d <= s) ? 0 : 1;
}

/**
 * Quantise the value of real distance (d) given the tolerance value (s).
 * For LCSS, the value will be oposite: similar points (d <= s) is quantised as 1 and vice versa.
 */
double lcss(double d, double s) {
    return (d > s) ? 0 : 1;
}

/**
 * Some edit distance algorithms don't quantise the distance (i.e. they use real distance)
 * This dummy function will just pass the distance right back
 */
double no_qtse(double d, double s) {
    return d;
}

/* ------------------------------- interface -------------------------------- */
/*
 * Function:  choose_dist
 * --------------------
 *  Chooses right distance function(euclid or euclid_squared)
 *  int dist_type: settings.dist_type [_MANHATTAN, _EUCLID, _EUCLID_SQUARED]
 *  returns: dist_fptr, pointer to a distance function
 */
dist_fptr choose_dist(int dist_type) {
    switch (dist_type) {
        case _EUCLID:
            return &euclid;
        case _EUCLID_SQUARED:
            return &euclid_square;
        case _MANHATTAN:
            return &manhattan;
        default:
            return NULL;
    }
}


qtse_fptr choose_quantisation(int qtse_type) {
    switch (qtse_type) {
        case _EDR:
            return &edr;
        case _LCSS:
            return &lcss;
        default:
            return &no_qtse;
    }
}


/*
 * Function:  choose_dp
 * --------------------
 *  Chooses right step function(without traceback)
 *  int dp_type: settings.dp_type, step function type
 *  returns: dist_fptr, pointer to a step function without traceback
 */
dp_fptr choose_dp(int dp_type) {

    switch (dp_type) {
        case _DP1:
            return &dp1;
        case _DP2:
            return &dp2;
        case _DP3:
            return &dp3;
        case _DP3_LCSS:
            return &dp3_lcss;
        case _DPW:
            return &dpw;
        case _DP2_EDR:
            return &dp2_edr;
        case _DP2_ERP:
            return &dp2_erp;
        case _SCP0SYM:
            return &p0_sym;
        case _SCP0ASYM:
            return &p0_asym;
        case _SCP1DIV2SYM:
            return &p1div2_sym;
        case _SCP1DIV2ASYM:
            return &p1div2_asym;
        case _SCP1SYM:
            return &p1_sym;
        case _SCP1ASYM:
            return &p1_asym;
        case _SCP2SYM:
            return &p2_sym;
        case _SCP2ASYM:
            return &p2_asym;
        default:
            return &dp2;
    }
}

/*
 * Function:  choose_dpdir
 * --------------------
 *  Chooses right step function(with traceback)
 *  int dp_type: settings.dp_type, step function type
 *  returns: dpdir_fptr, pointer to a step function with traceback
 */
dpdir_fptr choose_dpdir(int dp_type) {
    switch (dp_type) {
        case _DP1:
            return &dp1dir;
        case _DP2:
            return &dp2dir;
        case _DP3:
            return &dp3dir;
        case _DP3_LCSS:
            return &dp3_lcss_dir;
        case _DPW:
            return &dpwdir;
        case _DP2_EDR:
            return &dp2_edr_dir;
        case _DP2_ERP:
            return &dp2_erp_dir;
        case _SCP0SYM:
            return &p0_symdir;
        case _SCP0ASYM:
            return &p0_asymdir;
        case _SCP1DIV2SYM:
            return &p1div2_symdir;
        case _SCP1DIV2ASYM:
            return &p1div2_asymdir;
        case _SCP1SYM:
            return &p1_symdir;
        case _SCP1ASYM:
            return &p1_asymdir;
        case _SCP2SYM:
            return &p2_symdir;
        case _SCP2ASYM:
            return &p2_asymdir;
        default:
            return &dp2dir;
    }
}

/*
 * Function:  choose_window
 * --------------------
 *  NOTE: function is partly redundant at the moment, it always returns  itakura
 *  Chooses right window function
 *  int dp_type: settings.win , step function type
 *  returns: dpdir_fptr, pointer to a step function with traceback
 */
window_fptr choose_window(struct t_settings *settings) {
    switch (settings->window_type) {
        case _SCBAND:
            return &scband;
        case _PALIVAL:
            return &palival;
        case _ITAKURA:
            return &itakura;
        case _PALIVAL_MOD:
            return &palival_mod;
        default:
            return &nowindow;
    }
}

/*
 * Function:  choose_window_param
 * --------------------
 *  Computes and return parameter for the window function
 *  struct t_settings *settings: structure with edit distance settings
 *  int len_ref: length of the reference sequence
 *  int len_query: length of the query sequence
 *  returns: double p, window parameter
 */
double choose_window_param(struct t_settings *settings,
                           int len_ref,
                           int len_query) {
    double p = 0;
    if (settings->window_type == _ITAKURA)
        p = len_ref / (double) len_query;
    else
        p = settings->window_param;
    return p;

}

/*
 * Function:  choose_path_pattern
 * --------------------
 *  Chooses right path_patterns(2d array) based on the current step function
 *  struct t_settings *settings: structure with edit distance settings
 *  returns: const int[7][11], path_pattern
 */
const int (*choose_path_pattern(struct t_settings settings))[11] {
    switch (settings.dp_type) {
        case _DP1:
        case _DP2:
        case _DP2_EDR:
        case _SCP0ASYM:
        case _SCP0SYM:
            return dp2_path_pattern;
        case _DP3:
        case _DP3_LCSS:
            return dp1_path_pattern;
        case _SCP1DIV2SYM:
        case _SCP1DIV2ASYM:
            return p1div2_path_pattern;
        case _SCP1SYM:
        case _SCP1ASYM:
            return p1_path_pattern;
        case _SCP2SYM:
        case _SCP2ASYM:
            return p2_path_pattern;
        default:
            return dp2_path_pattern;
    }
}

/*
 * Function:  extra_size
 * --------------------
 *  Computes extra size for the cost matrix (settings.offset)
 *  int dp_type: settings.dp_type, step pattern type
 *  returns: int (1,2 or 3)
 */
int extra_size(int dp_type) {
    switch (dp_type) {
        case _DP1:
        case _DP2:
        case _DP3:
        case _DP3_LCSS:
        case _DP2_EDR:
        case _DP2_ERP:
        case _SCP0ASYM:
        case _SCP0SYM:
            return 1;
        case _SCP1ASYM:
        case _SCP1SYM:
            return 2;
        default:
            return 3;
    }
}

/*
 * Function:  create_path_from_pattern
 * --------------------
 *  Compute full path from path_points and path_pattern. This function is
 *  necessary for the multistep step functions (like p1div2_sym), and also
 *  to transform "directions" to coordinates(matrix indices).
 *  const int pattern[6][11] : path pattern (see ced.h)
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
 * Function:  cednopath
 * --------------------
 * Edit Distance (e.g. DTW, EDR) algorithm(without traceback)
 * double* ref:         reference sequence
 * double* query:       query sequence
 * doubel* args.sigma:       maximum difference allowances of all dimension (usable by EDR)
 * int len_ref:         length of the reference sequence
 * int len_query:       length of the query sequence
 * int ncols:           number of columns of the second dimension (1 means 1 dimensional array)
 * dist_fptr dist:      pointer to a distance function
 * dp_fptr dp:          pointer to a step function
 * window_fptr window:  pointer to a window function
 * double p:            window parameter
 * double *cost_matrix: cost matrix
 * struct t_settings settings: structure with edit distance settings
 * returns: doube, distance between reference and query sequences
 */
double cednopath(double *ref,
                 double *query,
                 t_extra_args args,
                 int len_ref,
                 int len_query,
                 int ncols,
                 dist_fptr dist,
                 qtse_fptr quantise,
                 dp_fptr dp,
                 window_fptr window,
                 double p, //window param
                 double *cost_matrix,
                 struct t_settings settings) {
    int off = settings.offset;
    /*memory was already allocated*/
    /*extending matrix*/
    int M = len_ref + off;
    int N = len_query + off;
    int i = 0;
    int j = 0;
    double w = 0;
    double s = 0;
    bool fast_glob = (settings.window_type == _PALIVAL ||
                      settings.window_type == _PALIVAL_MOD);
    /*no window or fast window case*/
    if (fast_glob || settings.window_type == 0) {
        if (fast_glob) {
            if (settings.window_type == _PALIVAL_MOD)
                w = settings.window_param * (double) len_ref;
            else
                w = settings.window_param;
            s = len_query / (double) len_ref;
        }
        else {
            w = INFINITY;
            s = 1;
        }

        cost_matrix[idx(off, off, N)] = quantise(dist(ref, off, query, off, ncols), args.sigma);
        for (j = max2(off + 1, _round(s * (off - w))); j < min2(N, _round(s * (off + w) + 1)); j++) {
            cost_matrix[idx(off, j, N)] = dp(ref, query, args, cost_matrix, off, j, ncols, &settings, N, dist,
                                             quantise);
        }

        for (i = off + 1; i < M; i++) {
            for (j = max2(off, _round(s * (i - w))); j < min2(N, _round(s * (i + w) + 1)); j++)
                cost_matrix[idx(i, j, N)] = dp(ref, query, args, cost_matrix, i, j, ncols, &settings, N, dist,
                                               quantise);
        }

    }
        /*slow window case*/
    else {
        cost_matrix[idx(off, off, N)] = quantise(dist(ref, off, query, off, ncols), args.sigma);
        for (j = off + 1; j < N; j++) {
            if (window(off, j, p, len_ref, len_query))
                cost_matrix[idx(off, j, N)] = dp(ref, query, args, cost_matrix, off, j, ncols, &settings, N, dist,
                                                 quantise);
        }
        for (i = off + 1; i < M; i++) {
            for (j = off; j < N; j++) {
                if (window(i, j, p, len_ref, len_query))
                    cost_matrix[idx(i, j, N)] = dp(ref, query, args, cost_matrix, i, j, ncols, &settings, N, dist,
                                                   quantise);
            }
        }
    }

    return cost_matrix[idx(M - 1, N - 1, N)];

}

/*
 * Function:  cedpath
 * --------------------
 * Edit Distance (e.g. DTW, EDR) algorithm(with traceback)
 * double* ref:         reference sequence
 * double* query:       query sequence
 * doubel* args.sigma:       maximum difference allowances of all dimension (usable by EDR)
 * int len_ref:         length of the reference sequence
 * int len_query:       length of the query sequence
 * int ncols:           number of columns of the second dimension (1 means 1 dimensional array)
 * dist_fptr dist:      pointer to a distance function
 * dpdir_fptr dp_dir:   pointer to a step function
 * window_fptr window:  pointer to a window function
 * double p:            window parameter
 * double *cost_matrix: cost matrix
 * int *dir_matrix:     direction matrix
 * struct t_settings settings: structure with edit distance settings
 * returns: doube, distance between reference and query sequences
 */
double cedpath(double *ref,
               double *query,
               t_extra_args args,
               int len_ref,
               int len_query,
               int ncols,
               dist_fptr dist,
               qtse_fptr quantise,
               dpdir_fptr dp_dir,
               window_fptr window,
               double p, //window param
               double *cost_matrix,
               int *dir_matrix,
               struct t_settings settings) {
    int off = settings.offset;

    struct t_item item = {0, 0};
    /*extending matrix*/
    int M = len_ref + off;
    int N = len_query + off;

    int i = 0;
    int j = 0;
    double w = 0;
    double s = 0;
    bool fast_glob = (settings.window_type == _PALIVAL ||
                      settings.window_type == _PALIVAL_MOD);
    /*no window or fast window case*/
    if (fast_glob || settings.window_type == 0) {
        if (fast_glob) {
            if (settings.window_type == _PALIVAL_MOD)
                w = settings.window_param * (double) len_ref;
            else
                w = settings.window_param;
            s = len_query / (double) len_ref;
        }
        else {
            w = INFINITY;
            s = 1;
        }

        cost_matrix[idx(off, off, N)] = quantise(dist(ref, off, query, off, ncols), args.sigma);
        for (j = max2(off + 1, _round(s * (off - w))); j < min2(N, _round(s * (off + w) + 1)); j++) {
            item = dp_dir(ref, query, args, cost_matrix, off, j, ncols, &settings, N, dist, quantise);
            cost_matrix[idx(off, j, N)] = item.val;
            dir_matrix[idx(0, j - off, N - off)] = item.idx;
        }

        for (i = off + 1; i < M; i++) {
            for (j = max2(off, _round(s * (i - w))); j < min2(N, _round(s * (i + w) + 1)); j++) {
                item = dp_dir(ref, query, args, cost_matrix, i, j, ncols, &settings, N, dist, quantise);
                cost_matrix[idx(i, j, N)] = item.val;
                dir_matrix[idx(i - off, j - off, N - off)] = item.idx;
            }
        }
    }
        /*slow window case*/
    else {
        cost_matrix[idx(off, off, N)] = quantise(dist(ref, off, query, off, ncols), args.sigma);
        for (j = off + 1; j < N; j++) {
            if (window(off, j, p, len_ref, len_query)) {
                item = dp_dir(ref, query, args, cost_matrix, off, j, ncols, &settings, N, dist, quantise);
                cost_matrix[idx(off, j, N)] = item.val;
                dir_matrix[idx(0, j - off, N - off)] = item.idx;
            }
        }
        for (i = off + 1; i < M; i++) {
            for (j = off; j < N; j++) {
                if (window(i, j, p, len_ref, len_query)) {
                    item = dp_dir(ref, query, args, cost_matrix, i, j, ncols, &settings, N, dist, quantise);
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
 * doubel* args.sigma:       maximum difference allowances of all dimension (usable by EDR)
 * int len_ref:         length of the reference sequence
 * int len_query:       length of the query sequence
 * struct t_settings settings: structure with edit distance settings
 * returns: void
 */
void fill_matrix(double *matrix,
                 int len_ref,
                 int len_query,
                 struct t_settings settings) {
    /*
    http://www.cplusplus.com/reference/cstring/memset/
    or for 1 dimension {INFINITY}
    */
    int M = len_ref + settings.offset;
    int N = len_query + settings.offset;
    int i = 0;
    int j = 0;

    typedef enum FILL_TYPE {
        FILL_ALL_ZERO,                 ///< If the step pattern is _DP3_LCSS, fill the entire matrix with zero
        FILL_ALL_INFINITY,             ///< If there is a window or complicated step pattern, fill the entire matrix with INF
        FILL_INFINITY_RIM_ZERO_INSIDE  ///< All other case fill the rim with INF and everything inside zero
    } FILL_TYPE;

    FILL_TYPE fill_type;

    switch (settings.dp_type) {
        case _DP1:
        case _DP2:
        case _DP3:
        case _SCP0SYM:
            fill_type = settings.window_type != 0 ? FILL_ALL_INFINITY : FILL_INFINITY_RIM_ZERO_INSIDE;
            break;
        case _DP3_LCSS:
            fill_type = FILL_ALL_ZERO;
            break;
        default:
            fill_type = FILL_ALL_INFINITY;
            break;
    }

    switch (fill_type) {

        case FILL_ALL_ZERO:
            for (i = 0; i < M; i++)
                for (j = 0; j < N; j++)
                    matrix[idx(i, j, N)] = 0;
            break;
        case FILL_ALL_INFINITY:
            for (i = 0; i < M; i++)
                for (j = 0; j < N; j++)
                    matrix[idx(i, j, N)] = INFINITY;
            break;
        case FILL_INFINITY_RIM_ZERO_INSIDE:
            for (i = 0; i < M; i++)
                for (j = 0; j < settings.offset; j++)
                    matrix[idx(i, j, N)] = INFINITY;

            for (i = 0; i < settings.offset; i++)
                for (j = 0; j < N; j++)
                    matrix[idx(i, j, N)] = INFINITY;
            break;
    }
}

/*
 * Function:  ced
 * --------------------
 * Edit Distance (e.g. DTW, EDR)
 * This fuction is main entry for the ced
 * double* ref:                 reference sequence
 * double* query:               query sequence
 * doubel* args.sigma:               maximum difference allowances of all dimension (usable by EDR)
 * doubel* args.sigma:               maximum difference allowances of all dimension (usable by EDR)
 * int len_ref:                 length of the reference sequence
 * int len_query:               length of the query sequence
 * int ncols:                   number of columns of the second dimension (1 means 1 dimensional array)
 * double *cost_matrix:         cost matrix, pointer to an allocated memory for
                                the 2 dimensional matrix. (with extra size)
 * struct t_path_element* path: warping path, pointer to an allocated array,
                                size is maximum possible path length:
                                len_ref + len_query
 * int true_path_len:           length of the computed warping path
 * struct t_settings settings: structure with edit distance settings
 * returns: doube, distance between reference and query sequences
 */
double ced(double *ref,
           double *query,
           t_extra_args args,
           int len_ref,
           int len_query,
           int ncols,
           double *cost_matrix,
           struct t_path_element *path,
           int *true_path_len,
           struct t_settings settings) {
    /*init distance*/
    double distance = 0;
    /*init window function param*/
    double p = 0;
    dp_fptr dp = NULL;
    /*init pointer to step function*/
    dpdir_fptr dp_dir = NULL;

    /*pointer to window function*/
    window_fptr window = choose_window(&settings);
    /*pointer to distance function*/
    dist_fptr dist = choose_dist(settings.dist_type);
    /*pointer to the quantisation function*/
    qtse_fptr quantise = choose_quantisation(settings.qtse_type);

    /* Normalise the tolerance if an array is given */
    if (args.sigmas != NULL) {
        args.sigma = norm(args.sigmas, ncols, dist);
    }

    int *dir_matrix = NULL;

    int *path_points = NULL;

    int path_points_count;

    /*assign step function*/
    if (settings.compute_path)
        dp_dir = choose_dpdir(settings.dp_type);
    else
        dp = choose_dp(settings.dp_type);


    /*assign window parameter*/
    p = choose_window_param(&settings,
                            len_ref,
                            len_query);


    /*lets go !*/

    /*prepare cost matrix*/
    fill_matrix(cost_matrix, len_ref, len_query, settings);

    /*edit distance without traceback case(only cost_matrix and distance)*/
    if (settings.compute_path == false) {

        distance = cednopath(ref,
                             query,
                             args,
                             len_ref,
                             len_query,
                             ncols,
                             dist,
                             quantise,
                             dp,
                             window,
                             p,
                             cost_matrix,
                             settings);
    }
        /*edit distance with traceback case*/
    else {
        /*allocate direction matrix*/
        dir_matrix = (int *) malloc(sizeof(int) * len_ref * len_query);
        /*call cedpath, computes distance, cost matrix and direction matrix*/
        distance = cedpath(ref,
                           query,
                           args,
                           len_ref,
                           len_query,
                           ncols,
                           dist,
                           quantise,
                           dp_dir,
                           window,
                           p,
                           cost_matrix,
                           dir_matrix,
                           settings);

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
                choose_path_pattern(settings));


        /*cleaning*/
        free(dir_matrix);
        /*compute warping path(finnally, as array of the path elements*/
        *true_path_len = create_path_from_pattern(
                choose_path_pattern(settings),
                len_ref,
                len_query,
                path_points,
                path_points_count,
                path);
        /*cleaning*/
        free(path_points);

    }

    /*normalization case*/
    double normed = normalise(distance, len_ref, len_query, settings.norm_type);

    if (settings.dp_type == _DP3_LCSS) {
        return 1.0 - normed;
    }
    return normed;
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

/**
 * Calculate the norm_type of an array given the distance function
 */
double norm(double *arr, int n, double (*dist)(double *, int, double *, int, int)) {
    double *origin = malloc(sizeof(double) * n);
    int i;
    for (i = 0; i < n; i++) {
        origin[i] = 0;
    }
    double result = dist(arr, 0, origin, 0, n);
    free(origin);
    return result;
}

/**
 * Normalise the distance given the normalisation type
 */
double normalise(double d, int len_ref, int len_query, int norm_type) {
    switch (norm_type) {
        case _NORM_BY_MIN_LENGTH:
            return d / (double) min2(len_ref, len_query);
        case _NORM_BY_MAX_LENGTH:
            return d / (double) max2(len_ref, len_query);
        case _NORM_BY_AVG_LENGTH:
            return d / ((len_ref + len_query) / 2.0);
        default:
            return d;
    }
}


//int main() {
//    int i, j;
//
////    int ncols = 2;
////    double r[] = {1, 0, 5, 0, 4, 0, 2, 0};
////    double q[] = {1, 0, 2, 0, 4, 0, 1, 0};
//
//    int ncols = 1;
//    double r[] = {4, 4, 2, 4};
//    double q[] = {5, 4, 5, 6, 4};
//    t_extra_args args;
//
//    args.sigmas = malloc(sizeof(double) * ncols);
//    args.gap = malloc(sizeof(double) * ncols);
//    for (i = 0; i < ncols; i++) {
//        args.gap[i] = 0;
//        args.sigmas[i] = 1;
//    }
//
//    int len_ref = sizeof(r) / sizeof(r[0]) / ncols;
//    int len_query = sizeof(q) / sizeof(q[0]) / ncols;
//    struct t_settings settings;
//
//    settings.compute_path = false;
//    settings.dist_type = _EUCLID;
//    settings.dp_type = _DP3_LCSS;
//    settings.qtse_type = _LCSS;
////    settings.window_type = _PALIVAL;
//    settings.window_param = 1;
//    settings.norm_type = _NORM_BY_AVG_LENGTH;
//    settings.offset = extra_size(settings.dp_type);
//
//    int offset = settings.offset * ncols;
//
//    double *expanded_r = malloc(sizeof(double) * (len_ref * ncols + offset));
//    double *expanded_q = malloc(sizeof(double) * (len_query * ncols + offset));
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
//            sizeof(double) * ((len_ref + settings.offset) * (len_query + settings.offset)));
//
//
//    struct t_path_element *path = (struct t_path_element *) malloc(
//            sizeof(struct t_path_element) * (len_ref + len_query));
//    int path_len = 0;
//    double dist = ced(expanded_r, expanded_q, args, len_ref, len_query, ncols, cost, path, &path_len, settings);
//    if (cost[14 * 14 - 1] == 4282.0)
//        printf("\nFAIL");
//    print_matrix(cost, len_ref + settings.offset, len_query + settings.offset);
//    for (i = 0; i < path_len; i++)
//        printf("i: %d  j: %d\n", path[i].i, path[i].j);
//
//    printf("Dist = %.2f", dist);
//    return 0;
//}

