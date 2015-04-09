#ifndef _CDTW_H_

#define _CDTW_H_

#ifdef _MSC_VER
#define INFINITY ((DBL_MAX+DBL_MAX))
#endif



typedef int bool;

#define true 1
#define false 0

// #define UP 1
// #define LEFT 0
// #define DIAG 2

/*distance types*/
#define _EUCLID 11          /*euclidean distance*/
#define _EUCLID_SQUARED 12  /*squared euclidean distance*/
#define _MANHATTAN 13       /*manhattan distance*/

/*step pattern types, read struct dtw_settings for more info*/
#define _DP1            21  
#define _DP2            22
#define _DP3            23
#define _SCP0SYM        24
#define _SCP0ASYM       25
#define _SCP1DIV2SYM    26
#define _SCP1DIV2ASYM   27
#define _SCP1SYM        28
#define _SCP1ASYM       29
#define _SCP2SYM        210
#define _SCP2ASYM       211

#define _SCBAND         31 
#define _PALIVAL        32 
#define _ITAKURA        33 
#define _PALIVAL_MOD    35 


#define min2(a,b) ((a) < (b) ? a : b)
#define max2(a,b) ((a) > (b) ? a : b)
#define idx(i,j,J) ( (i)*(J) + (j) )


/*step functions args*/
#define _DP_ARGS    double* ref,\
                    double* query,\
                    double* cost_matrix,\
                    int i,\
                    int j,\
                    int offset,\
                    int size2, \
                    double (*dist)(double a, double b)


/*pointer to distance function*/
typedef double (*dist_fptr)(double a, double b);

/*pointer to step function*/
typedef double  (*dp_fptr)(_DP_ARGS);

/*pointer to step function with traceback*/
typedef struct t_item (*dpdir_fptr)(_DP_ARGS);

/*pointer to window function*/
typedef bool (*window_fptr)(int i, int j, double r, double I, double J);

/*warping path*/
struct t_path_element{
    int i;
    int j;
};

/*t_item returns all step functions with traceback*/
struct t_item{
    double val;
    int idx;
};

/*
 * Structure:  t_path_element
 * -------------------- 
 *  Contains dtw settings
 *
 *  int compute_path:   0 - to compute only distance and matrices
 *                      1 - to compute optimal path
 *                      default 0
 *  int dist_type:      distance function:
 *                          11 - euclidean          _EUCLID
 *                          12 - euclidean          _EUCLID_SQUARED
 *                          13 - manhattan          _MANHATTAN
 *                      default euclidean
 *  int dp_type:        step pattern:
 *                          21 - dp1                            _DP1
 *                          22 - dp2                            _DP2
 *                          23 - dp3                            _DP3
 *                          24 - Sakoe-Chiba P0 symmetric       _SCP0SYM        
 *                          25 - Sakoe-Chiba P0 asymmetric      _SCP0ASYM
 *                          26 - Sakoe-Chiba P1/2 symmetric     _SCP1DIV2SYM
 *                          27 - Sakoe-Chiba P1/2 asymmetric    _SCP1DIV2ASYM
 *                          28 - Sakoe-Chiba P1 symmetric       _SCP1SYM
 *                          29 - Sakoe-Chiba P1 asymmetric      _SCP1ASYM
 *                          210 - Sakoe-Chiba P2 symmetric      _SCP2SYM
 *                          211 - Sakoe-Chiba P2 asymmetric     _SCP2ASYM
 *                          default - dp2
 *  int window_type:    global constraint(window)
 *                          31 - scband      _SCBAND
 *                          32 - palival     _PALIVAL 
 *                          33 - itakura     _ITAKURA
 *                          35 - palivan_mod _PALIVAL_MOD
 *                        otherwise - full matrix(without window)
 *
 *  double window_param: param for global constraint
 *  int norm:               normalization 
 *                          0 - without normalization
 *                          1 - with normalization factor
 *  int offset:   number of extra rows and columns for cost matrix
 *
 */
struct t_dtw_settings{

        int compute_path;
        int dist_type;  
        int dp_type;
        int window_type;
        double window_param;
        int norm;
        int offset; 

};

/*path patterns for step functions*/
/*
path pattern is a matrix like:
'Sakoe-Chiba 1/2 sym' 
 [5,-1,-3,-1,-2,-1,-1,-2,-1,-3,-1] <- 5 path_points (-1-,3),(-1,-2),...
 [3, 0,-1, 0,-2,-1,-3, 0, 0, 0, 0] <- path to first target with 
 [2, 0,-1,-1,-2, 0, 0, 0, 0, 0, 0],                3 parts (0,-1),(0,-2),..
 [1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0],
 [2,-1, 0,-2,-1, 0, 0, 0, 0, 0, 0],
 [3,-1, 0,-2, 0,-3,-1, 0, 0, 0, 0],
*/

const int dp1_path_pattern[6][11] =       {{ 2, 0,-1,-1, 0, 0, 0, 0, 0, 0, 0},
                                            {1, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0},
                                            {1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};


const int dp2_path_pattern[6][11] =       {{ 3, 0,-1,-1, 0,-1,-1, 0, 0, 0, 0},
                                            {1, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0},
                                            {1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                            {1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0},
                                            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};

const int p1div2_path_pattern[6][11] =     {{5,-1,-3,-1,-2,-1,-1,-2,-1,-3,-1},
                                            {3, 0,-1, 0,-2,-1,-3, 0, 0, 0, 0},
                                            {2, 0,-1,-1,-2, 0, 0, 0, 0, 0, 0},
                                            {1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0},
                                            {2,-1, 0,-2,-1, 0, 0, 0, 0, 0, 0},
                                            {3,-1, 0,-2, 0,-3,-1, 0, 0, 0, 0}};

const int p1_path_pattern[6][11] =     {{3,-1,-2,-1,-1,-2,-1, 0, 0, 0, 0},
                                        {2, 0,-1,-1,-2, 0, 0, 0, 0, 0 ,0},
                                        {1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0},
                                        {2,-1, 0,-2,-1, 0, 0, 0, 0, 0, 0},
                                        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};


const int p2_path_pattern[6][11] =     {{3,-2,-3,-1,-1,-3,-2, 0, 0, 0, 0},
                                        {3, 0,-1,-1,-2,-2,-3, 0, 0, 0 ,0},
                                        {1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0},
                                        {3,-1, 0,-2,-1,-3,-2, 0, 0, 0, 0},
                                        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};


/*
 * Function:  extra_size
 * --------------------
 *  Computes extra size for the cost matrix (dtw_settings.offset)
 *  int dp_type: dtw_settings.dp_type, step pattern type
 *  returns: int (1,2 or 3)
 */
int extra_size(int dp_type);
/*--------------------finds minimum value----------------------*/
/*
 * Function:  min3
 * -------------------- 
 *  Finds the minimum of 3 doubles

 *  double x:    x
 *  double y:    y
 *  double z:    z

 *  returns:  min3(x,y,z)
 */
double min3(double x, double y, double z);

/*
 * Function:  min3
 * -------------------- 
 *  Finds the minimum of n doubles

 *  double *arr:    array
 *  int n: length of the array

 *  returns:  min(arr)
 */
double min_n(double* arr, int n);

/*--------------------find minimum value + index--------------*/
/*
 * Function:  min2idx
 * -------------------- 
 *  Finds the minimum of 2 doubles and its position (0 or 1)

 *  double a: a 
 *  double b: b

 *  returns:  struct t_item, t_item = { min(a,b), position }
 */
struct t_item min2idx(double a, double b);
/*
 * Function:  min3idx
 * -------------------- 
 *  Finds the minimum of 3 doubles and its position (0,1 or 2)

 *  double x: x 
 *  double y: y
 *  double z: z

 *  returns:  struct t_item, t_item = { min(x,y,z), position }
 */
struct t_item min3idx(double x, double y, double z);

/*
 * Function:  min3idx
 * -------------------- 
 *  Finds the minimum of n doubles and its position (0..n-1)

 *  double x: array
 *  int n: size of the arr

 *  returns:  struct t_item, t_item = { min(arr), position }
 */
struct t_item min_nidx(double* arr, int n);

/*----------------------------------------------------------------------------*/
/* -------------------------------- step patterns ----------------------------*/
/*----------------------------------------------------------------------------*/
/*--------------------------usual------------------------------*/
/*
 * Function:  dp1
 * -------------------- 
 
 * NOTE: ALL STEP FUNCTIONS have the same args, defined in macro 
 *_DP_ARGS(cdtw.h). Step functions like step_pattern_type and 
 * step_pattern_typedir are pretty similar, step_pattern_type are used in 
 * computing dtw without path(without traceback). step_pattern_typedir are 
 * used in computing dtw with path(traceback) 

 *  Step patern dp1:
 *  min(      
 *      cost_matrix[i][j-1]   +   d(r[i],q[j])    
 *       cost_matrix[i-1][j]   +   d(r[i],q[j]),    
 *       cost_matrix[i-1][j-1] + 2*d(r[i],q[j])
 *      )
 *
 * double* ref: reference sequence
 * double* query: query sequence
 * double* cost_matrix: cost matrix
 * int i: row index of the cost matrix
 * int j: column index of the cost matrix
 * int offset: extra size of the cost matrix
 * int size2:  cost matrix columns count 
 * double (*dist)(double a, double b): poiter to distance function

 * returns:  double, value to assign cost_matrix[i][j]
*/
double dp1(_DP_ARGS);

/*
 * Function:  dp2
 * -------------------- 
 *  Step patern dp2:
 *  min(      
 *      cost_matrix[i][j-1]   +   d(r[i],q[j])    
 *      cost_matrix[i-1][j]   +   d(r[i],q[j]),    
 *      cost_matrix[i-1][j-1] +   d(r[i],q[j])
 *      )
 * see doc for the dp1
*/
double dp2(_DP_ARGS);

/*
 * Function:  dp3
 * -------------------- 
 *  Step patern dp3:
 *  min(      
 *      cost_matrix[i][j-1]   +   d(r[i],q[j])    
 *      cost_matrix[i-1][j]   +   d(r[i],q[j]),    
 *      )
 * see doc for the dp1
*/
double dp3(_DP_ARGS);
/* ------------------------usual for traceback-----------------*/
/*
 * Function:  dp1dir
 * -------------------- 
 * see doc for the dp1
*/
struct t_item dp1dir(_DP_ARGS);

/*
 * Function:  dp2dir
 * -------------------- 
 * see doc for the dp1,dp2
*/
struct t_item dp2dir(_DP_ARGS);

/*
 * Function:  dp3dir
 * -------------------- 
 * see doc for the dp1,dp3
*/
struct t_item dp3dir(_DP_ARGS);
/*-------------------Sakoe-Chiba classifikaction---------------*/
/*
 * Function:  p0sym
 * -------------------- 
 * Sakoe-Chiba classification p = 0, symmetric step pattern
 * This function is alias for the dp1
 * see doc for the dp1
*/
double p0_sym(      _DP_ARGS);

/*
 * Function:  p0asym
 * -------------------- 
 * Sakoe-Chiba classification p = 0, asymmetric step pattern:
 *  min(      
 *       cost_matrix[i][j-1]   +   0    
 *       cost_matrix[i-1][j]   +   d(r[i],q[j]), 
 *       cost_matrix[i-1][j-1] +   d(r[i],q[j]),   
 *      )
 * see doc for the dp1
*/
double p0_asym(     _DP_ARGS);

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
 *      )
 * see doc for the dp1
*/
double p1div2_sym(  _DP_ARGS);

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
 *      )
 * see doc for the dp1
*/
double p1div2_asym( _DP_ARGS);

/*
 * Function:  p1_sym
 * -------------------- 
 * Sakoe-Chiba classification p = 1, symmetric step pattern:
 *  min(        
 *      cost_matrix[i-1][j-2] + 2d(r[i],q[j-1]) + d(r[i],q[j]),     
 *      cost_matrix[i-1][j-1] + 2d(r[i],q[j]), 
 *      cost_matrix[i-2][j-1] + 2d(r[i-1],q[j]) + d(r[i],q[j]), 
 *      )
 * see doc for the dp1
*/
double p1_sym(      _DP_ARGS);

/*
 * Function:  p1_asym
 * -------------------- 
 * Sakoe-Chiba classification p = 1, asymmetric step pattern:
 *  min(        
 *      cost_matrix[i-1][j-2] + (d(r[i],q[j-1]) + d(r[i],q[j]))/2,     
 *      cost_matrix[i-1][j-1] + d(r[i],q[j]), 
 *      cost_matrix[i-2][j-1] + d(r[i-1],q[j]) + d(r[i],q[j]), 
 *      )
 * see doc for the dp1
*/
double p1_asym(     _DP_ARGS);

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
 * see doc for the dp1
*/
double p2_sym(      _DP_ARGS);

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
 * see doc for the dp1
*/
double p2_asym(     _DP_ARGS);
/*---------------step patterns for traceback-------------------*/
/*
 * Function:  p0_symdir
 * -------------------- 
 * Sakoe-Chiba classification p = 0, symmetric step pattern:
 * see doc for the p0_sym, dp1
*/
struct t_item p0_symdir(        _DP_ARGS);

/*
 * Function:  p0_asymdir
 * -------------------- 
 * Sakoe-Chiba classification p = 0, asymmetric step pattern:
 * see doc for the p0_asym, dp1
*/
struct t_item p0_asymdir(       _DP_ARGS);

/*
 * Function:  p1div2_symdir
 * -------------------- 
 * Sakoe-Chiba classification p = 0.5, symmetric step pattern:
 * see doc for the p1div2_sym, dp1
*/
struct t_item p1div2_symdir(    _DP_ARGS);

/*
 * Function:  p1div2_asymdir
 * -------------------- 
 * Sakoe-Chiba classification p = 0.5, asymmetric step pattern:
 * see doc for the p1div2_asym, dp1
*/
struct t_item p1div2_asymdir(   _DP_ARGS);

/*
 * Function:  p1_symdir
 * -------------------- 
 * Sakoe-Chiba classification p = 1, symmetric step pattern:
 * see doc for the p1_sym, dp1
*/
struct t_item p1_symdir(        _DP_ARGS);

/*
 * Function:  p1_asymdir
 * -------------------- 
 * Sakoe-Chiba classification p = 1, asymmetric step pattern:
 * see doc for the p1_asym, dp1
*/
struct t_item p1_asymdir(     _DP_ARGS);

/*
 * Function:  p2_symdir
 * -------------------- 
 * Sakoe-Chiba classification p = 2, symmetric step pattern:
 * see doc for the p2_sym, dp1
*/
struct t_item p2_symdir(        _DP_ARGS);

/*
 * Function:  p2_asymdir
 * -------------------- 
 * Sakoe-Chiba classification p = 2, asymmetric step pattern:
 * see doc for the p2_asym, dp1
*/
struct t_item p2_asymdir(       _DP_ARGS);

/*------------------------ window functions ------------------*/

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
bool scband(     int i, int j, double r, double I, double J);

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
bool palival(    int i, int j, double r, double I, double J);

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
bool itakura(    int i, int j, double k, double I, double J);

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
bool palival_mod(int i, int j, double p, double I, double J);

/*
 * Function: nowindow
 * --------------------
 * This function is for testing only
 * returns: always true(1)
 */
bool nowindow(   int i, int j, double k, double I, double J);
 /*------------------------distance function------------------- */

/*
 * Function:  manhattan
 * --------------------
 *  Euclidean distance 
 *  double a:  1 dimensional point  
 *  double b:  1 dimensional point  
 *  returns: double, manhattan distance between two dimensional points 
 */
double manhattan(double a, double b); 

/*
 * Function:  euclid
 * --------------------
 *  Euclidean distance helper
 *  double a:  1 dimensional point  
 *  double b:  1 dimensional point  
 *  returns: double, (a-b)^2 
 */
double euclid(double a, double b);
/*----------------------------------------------------------------------------*/
/* ------------------------------- interface -------------------------------- */
/*----------------------------------------------------------------------------*/
/*
 * Function:  choose_dist
 * --------------------
 *  Chooses right distance function(euclid or euclid_squared) 
 *  int dist_type: dtw_settings.dist_type [_MANHATTAN, _EUCLID, _EUCLID_SQUARED]
 *  returns: dist_fptr, pointer to a distance function
 */
dist_fptr choose_dist(int dist_type);

/*
 * Function:  choose_dp
 * --------------------
 *  Chooses right step function(without traceback)
 *  int dp_type: dtw_settings.dp_type, step function type
 *  returns: dist_fptr, pointer to a step function without traceback
 */
dp_fptr choose_dp(int dp_type);
/*
 * Function:  choose_dpdir
 * --------------------
 *  Chooses right step function(with traceback)
 *  int dp_type: dtw_settings.dp_type, step function type
 *  returns: dpdir_fptr, pointer to a step function with traceback
 */
dpdir_fptr choose_dpdir(int dp_type);

/*
 * Function:  choose_window
 * --------------------
 *  NOTE: function is partly redundant at the moment, it returns always itakura
 *  Chooses right window function
 *  int dp_type: dtw_settings.win , step function type
 *  returns: dpdir_fptr, pointer to a step function with traceback
 */
window_fptr choose_window(struct t_dtw_settings* dtw_settings);

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
                           int len_query);

/*
 * Function:  choose_path_pattern
 * --------------------
 *  Chooses right path_patterns(2d array) based on the current step function 
 *  struct t_dtw_settings *dtw_settings: structure with dtw settings
 *  returns: const int[7][11], path_pattern
 */
const int (*choose_path_pattern(struct t_dtw_settings dtw_settings))[11];

/*----------------------------------------------------------------------------*/
/*--------------------------------- core ------------------------------------ */
/*----------------------------------------------------------------------------*/
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
                                    int* path_points,
                                    int  path_points_count,
                                    struct t_path_element* path
                                    );

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
int direct_matrix_to_path_points(int* dir_matrix, int *path_points, 
                     int len_ref, int len_query,const int pattern[6][11]);

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
void fill_matrix(double *matrix, int len_ref, int len_query, struct t_dtw_settings dtw_settings);

/*---------------------------DTW--------------------------------*/
/*
 * Function:  cdtwnopath
 * --------------------
 * Dynamic Time Warping algorithm(without traceback)
 * double* ref:         reference sequence
 * double* query:       query sequence
 * int len_ref:         length of the reference sequence
 * int len_query:       length of the query sequence
 * dist_fptr dist:      pointer to a distance function
 * dp_fptr dp:          pointer to a step function
 * window_fptr window:  pointer to a window function
 * double p:            window parameter
 * double *cost_matrix: cost matrix 
 * struct t_dtw_settings dtw_settings: structure with dtw settings 
 * returns: doube, distance between reference and query sequences 
 */
double cdtwnopath(double* ref,
                  double* query,
                  int len_ref,
                  int len_query,
                  dist_fptr dist,
                  dp_fptr dp,
                  window_fptr window,
                  double p, //window param
                  double *cost_matrix,
                  struct t_dtw_settings dtw_settings);

/*
 * Function:  cdtwpath
 * --------------------
 * Dynamic Time Warping algorithm(with traceback)
 * double* ref:         reference sequence
 * double* query:       query sequence
 * int len_ref:         length of the reference sequence
 * int len_query:       length of the query sequence
 * dist_fptr dist:      pointer to a distance function
 * dpdir_fptr dp_dir:   pointer to a step function
 * window_fptr window:  pointer to a window function
 * double p:            window parameter
 * double *cost_matrix: cost matrix 
 * int *dir_matrix:     direction matrix
 * struct t_dtw_settings dtw_settings: structure with dtw settings 
 * returns: doube, distance between reference and query sequences 
 */
double cdtwpath(double* ref,
                  double* query,
                  int len_ref,
                  int len_query,
                  dist_fptr dist,
                  dpdir_fptr dp_dir,
                  window_fptr window,
                  double p, //window param
                  double *cost_matrix,
                  int *dir_matrix, 
                  struct t_dtw_settings dtw_settings);


/*
 * Function:  cdtw
 * --------------------
 * Dynamic Time Warping 
 * This fuction is main entry for the cdtw
 * double* ref:                  reference sequence
 * double* query:               query sequence
 * int len_ref:                 length of the reference sequence
 * int len_query:               length of the query sequence
 * double *cost_matrix:         cost matrix, pointer to an allocated memory for 
                                the 2 dimensional matrix. (with extra size)
 * struct t_path_element* path: warping path, pointer to an allocated array,  
                                size is maximum possible path length:
                                len_ref + len_query
 * int true_path_len:           length of the computed warping path
 * struct t_dtw_settings dtw_settings: structure with dtw settings 
 * returns: doube, distance between reference and query sequences 
 */
double cdtw(double* ref,
            double* query,
            int len_ref,
            int len_query,  
            double* cost_matrix,
            struct t_path_element* path, /*path maximum length = len_ref + len_query*/
            int *true_path_len,   /*actual computed path length*/
            struct t_dtw_settings dtw_settings);




/*helpers*/
void print_matrix(double* matrix, int size1, int size2);
void print_intarr(int* arr, int n);



#endif