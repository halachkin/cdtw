stepstr = dict()

stepstr['dp1'] = "\n \
 *  Step patern dp1: \n \
 *  min(\n \
 *       cost_matrix[i][j-1]   +   d(r[i],q[j]) \n \
 *       cost_matrix[i-1][j]   +   d(r[i],q[j]), \n \
 *       cost_matrix[i-1][j-1] + 2*d(r[i],q[j]) \n \
 *      )\n \
 "

stepstr['dp2'] = " \n \
 *  Step patern dp2:\n \
 *  min(\n \
 *      cost_matrix[i][j-1]   +   d(r[i],q[j]) \n \
 *      cost_matrix[i-1][j]   +   d(r[i],q[j]),\n \
 *      cost_matrix[i-1][j-1] +   d(r[i],q[j]) \n \
 *      )\n \
 "

stepstr['dp3'] = "\n \
 *  Step patern dp3:\n \
 *  min(      \n \
 *      cost_matrix[i][j-1]   +   d(r[i],q[j]) \n \
 *      cost_matrix[i-1][j]   +   d(r[i],q[j]),\n \
 *      )\n \
 "

stepstr['p0sym'] = "\n \
 *  Sakoe-Chiba classification p = 0, symmetric step pattern:\n \
 *  min(      \n \
 *      cost_matrix[i][j-1]   +   d(r[i],q[j]) \n \
 *      cost_matrix[i-1][j]   +   d(r[i],q[j]),\n \
 *      )\n \
 "

stepstr['p0asym'] = "\n \
 *  Sakoe-Chiba classification p = 0, asymmetric step pattern: \n \
 *  min(      \n \
 *       cost_matrix[i][j-1]   +   0   \n \
 *       cost_matrix[i-1][j]   +   d(r[i],q[j]), \n \
 *       cost_matrix[i-1][j-1] +   d(r[i],q[j]), \n \
 *      )\n \
 "


stepstr['p05sym'] = "\n \
 * Sakoe-Chiba classification p = 0.5, symmetric step pattern: \n \
 *  min( \n \
 *      cost_matrix[i-1][j-3] + 2d(r[i],q[j-2]) + d(r[i],q[j-1]) + d(r[i],q[j]),\n \
 *      cost_matrix[i-1][j-2] + 2d(r[i],q[j-1]) + d(r[i],q[j]), \n \
 *      cost_matrix[i-1][j-1] + 2d(r[i],q[j]), \n \
 *      cost_matrix[i-2][j-1] + 2d(r[i-1],q[j]) + d(r[i],q[j]), \n \
 *      cost_matrix[i-3][j-1] + 2d(r[i-2],q[j]) + d(r[i-1],q[j]) + d(r[i],q[j]) \n \
 *      )"



stepstr['p05asym'] =  "\n \
 * Sakoe-Chiba classification p = 0.5, asymmetric step pattern:\n \
 *  min( \n \
 *   cost_matrix[i-1][j-3] + (d(r[i],q[j-2]) + d(r[i],q[j-1]) + d(r[i],q[j]))/3 \n \
 *   cost_matrix[i-1][j-2] + (d(r[i],q[j-1]) + d(r[i],q[j]))/2, \n \
 *   cost_matrix[i-1][j-1] + d(r[i],q[j]), \n \
 *   cost_matrix[i-2][j-1] + d(r[i-1],q[j])  + d(r[i],q[j]),\n \
 *   cost_matrix[i-3][j-1] + d(r[i-2],q[j])  + d(r[i-1],q[j]) + d(r[i],q[j]) \n \
 *      )"


stepstr['p1sym'] = "\n \
 * Sakoe-Chiba classification p = 1, symmetric step pattern:\n \
 *  min(        \n \
 *      cost_matrix[i-1][j-2] + 2d(r[i],q[j-1]) + d(r[i],q[j]), \n \
 *      cost_matrix[i-1][j-1] + 2d(r[i],q[j]), \n \
 *      cost_matrix[i-2][j-1] + 2d(r[i-1],q[j]) + d(r[i],q[j]), \n \
 *      )"

stepstr['p1asym'] = "\n \
 * Sakoe-Chiba classification p = 1, asymmetric step pattern:\n \
 *  min( \n \
 *      cost_matrix[i-1][j-2] + (d(r[i],q[j-1]) + d(r[i],q[j]))/2, \n \
 *      cost_matrix[i-1][j-1] + d(r[i],q[j]), \n \
 *      cost_matrix[i-2][j-1] + d(r[i-1],q[j]) + d(r[i],q[j]), \n \
 *      ) \n \
" 

stepstr['p2sym'] = "\n \
 * Sakoe-Chiba classification p = 2, symmetric step pattern: \n \
 *  min(   \n \
 *     cost_matrix[i-2][j-3] + 2d(r[i-1],q[j-2]) + \n \
 *                             2d(r[i],q[j-1])   + \n \
 *                             d(r[i],q[j]), \n \
 * \n \
 *     cost_matrix[i-1][j-1] + 2d(r[i],q[j]), \n \
 * \n \
 *     cost_matrix[i-3][j-2] + 2d(r[i-2],q[j-1]) + \n \
 *                             2d(r[i-1],q[j])   + \n \
 *                             d(r[i],q[j]) \n \
 *     ) \n \
"


stepstr['p2asym'] = " \n \
 * Sakoe-Chiba classification p = 2, asymmetric step pattern: \n \
 *  min(   \n \
 *     cost_matrix[i-2][j-3] + 2( d(r[i-1],q[j-2]) +  \n \
 *                                d(r[i],q[j-1])   +  \n \
 *                                d(r[i],q[j]) ),     \n \
 * \n \
 *     cost_matrix[i-1][j-1] + d(r[i],q[j]), \n \
 * \n \
 *     cost_matrix[i-3][j-2] + d(r[i-2],q[j-1]) +    \n \
 *                             d(r[i-1],q[j])   +    \n \
 *                             d(r[i],q[j])          \n \
 *     )" 
