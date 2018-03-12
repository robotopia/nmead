// From https://github.com/matteotiziano/nelder-mead

#include "nmead.h"

int compare(const void *arg1, const void *arg2) {
    double fx1, fx2;
    fx1 = (((nm_point *)arg1)->fx);
    fx2 = (((nm_point *)arg2)->fx);
    
    if(fx1==fx2) {
        return 0;
    } else {
        return (fx1<fx2)? -1 : 1;
    }
}

void get_centroid(int n, nm_point *points, double *x_bar) {
    int i, j;
    for(j=0; j<n; j++) {
        x_bar[j] = 0;
        for(i=0; i<n; i++) {
            x_bar[j] += points[i].x[j];
        }
        x_bar[j] /= n;
    }
}

double modulus(double x) {
    return (x>0)? x : -x;
}

int continue_minimization(int n, nm_point *points, int eval_count, int iter_count, nm_optimset optimset) {
    int i,j;
    double condx = -1;
    double condf = -1;
    double temp;
    if(eval_count>optimset.max_eval || iter_count>optimset.max_iter) {
        // stop if #evals or #iters are greater than the max allowed
        return 0;
    }
    for(i=1; i<n+1; i++) {
        temp = modulus(points[0].fx-points[i].fx);
        if(condf<temp) {
            condf = temp;
        }
    }
    for(i=1; i<n+1; i++) {
        for(j=0; j<n; j++) {
            temp = modulus(points[0].x[j]-points[i].x[j]);
            if(condx<temp) {
                condx = temp;
            }
        }
    }
    // continue if both tolx or tolf condition is not met
    return condx>optimset.tolx || condf>optimset.tolf;
    
}

void get_point(int n, double *x, double *x_bar, double coeff, double *x_out) {
    int j;
    for(j=0; j<n; j++) {
        x_out[j] = (1+coeff)*x_bar[j] - coeff*x[j];
    }
}

void copy_nm_point(int n, double *x_from, double *x_to, double fx_from, double *fx_to) {
    int j;
    for(j=0; j<n; j++) {
        x_to[j] = x_from[j];
    }
    *fx_to = fx_from;
}

void swap_points(int n, nm_point *p1, nm_point *p2) {
    int j;
    double temp;
    for(j=0; j<n; j++) {
        temp     = p1->x[j];
        p1->x[j] = p2->x[j];
        p2->x[j] = temp;
    }
    temp   = p1->fx;
    p1->fx = p2->fx;
    p2->fx = temp;
}

void print_min(int n, nm_point *points) {
    int j;
    printf("[ ");
    for(j=0; j<n; j++) {
        printf("%.2f ", points[0].x[j]);
    }
    printf("]\t%.2f \n", points[0].fx);
}


void nelder_mead(double *x0, int n, nm_optimset optimset, nm_point *solution, double (*cost_fun)(long int, const double *, void *), void *cost_fun_args) {
    
    int        i, j;
    int        iter_count, eval_count;
    int        shrink;
    double     x_bar[n];
    double     x_r[n], x_e[n], x_c[n];
    double     fx_r, fx_e, fx_c;
    nm_point    points[n+1];
    
    iter_count = 0;
    eval_count = 0;
    
    // initial simplex
    for(i=0; i<n+1; i++) {
        points[i].x = malloc(n*sizeof(double));
        for(j=0; j<n; j++) {
            points[i].x[j] = (i-1==j)? ( x0[j]!=0? 1.05*x0[j] : 0.00025 ) : x0[j];
        }
        points[i].fx = cost_fun(n, points[i].x, cost_fun_args);
        eval_count++;
    }
    qsort((void *)points, n+1, sizeof(nm_point), compare);
    get_centroid(n, points, x_bar);
    iter_count++;
    
    // continue minimization until stop conditions are met
    while(continue_minimization(n, points, eval_count, iter_count, optimset)) {
        shrink = 0;
        
        get_point(n, points[n].x, x_bar, NM_RHO, x_r);
        fx_r = cost_fun(n, x_r, cost_fun_args);
        eval_count++;
        if(fx_r<points[0].fx) {
            get_point(n, points[n].x, x_bar, NM_RHO*NM_CHI, x_e);
            fx_e = cost_fun(n, x_e, cost_fun_args);
            eval_count++;
            if(fx_e<fx_r) {
                // expand
                copy_nm_point(n, x_e, points[n].x, fx_e, &(points[n].fx));
            } else {
                // reflect
                copy_nm_point(n, x_r, points[n].x, fx_r, &(points[n].fx));
            }
        } else {
            if(fx_r<points[n-1].fx) {
                // reflect
                copy_nm_point(n, x_r, points[n].x, fx_r, &(points[n].fx));
            } else {
                if(fx_r<points[n].fx) {
                    get_point(n, points[n].x, x_bar, NM_RHO*NM_GAMMA, x_c);
                    fx_c = cost_fun(n, x_c, cost_fun_args);
                    eval_count++;
                    if(fx_c<=fx_r) {
                        // contract outside
                        copy_nm_point(n, x_c, points[n].x, fx_c, &(points[n].fx));
                    } else {
                        // shrink
                        shrink = 1;
                    }
                } else {
                    get_point(n, points[n].x, x_bar, -NM_GAMMA, x_c);
                    fx_c = cost_fun(n, x_c, cost_fun_args);
                    eval_count++;
                    if(fx_c<=points[n].fx) {
                        // contract inside
                        copy_nm_point(n, x_c, points[n].x, fx_c, &(points[n].fx));
                    } else {
                        // shrink
                        shrink = 1;
                    }
                }
            }
        }
        if(shrink) {
            for(i=1; i<n+1; i++) {
                for(j=0; j<n; j++) {
                    points[i].x[j] = points[0].x[j] + NM_SIGMA*(points[i].x[j]-points[0].x[j]);
                }
                points[i].fx = cost_fun(n, points[i].x, cost_fun_args);
                eval_count++;
            }
            qsort((void *)points, n+1, sizeof(nm_point), compare);
        } else {
            for(i=n-1; i>=0 && points[i+1].fx<points[i].fx; i--) {
                swap_points(n, points+(i+1), points+i);
            }
        }
        get_centroid(n, points, x_bar);
        iter_count++;
    }
    
    // save solution in output argument
    copy_nm_point(n, points[0].x, solution->x, points[0].fx, &(solution->fx));

    // Free memory
    for(i=0; i<n+1; i++)
        free(points[i].x);
    
}


