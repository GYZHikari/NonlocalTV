#include <iostream>
#include <list>
#include <vector>
#include <cmath>
#include <mex.h>
using namespace std;


vector<int> find_index_num(vector<float> &W, float num, int sp_num)
{
    int ind_num = 0;
    vector<int> index(sp_num);
    for(int i = 0; i<sp_num; ++i)
    {
        index[i] = 0;
    }

    for(int i = 0; i<sp_num; ++i)
    {
        if(W[i]!=num)
        {
            index[ind_num] = i;
            ++ind_num;
        }
    }
    vector<int> find_num(ind_num);
    for(int i = 0; i<ind_num; ++i)
    {
        find_num[i] = index[i];
    }
    //find_num[ind_num+1] = ind_num;
    return find_num;
}

int index_num(vector<float> &W, float num, int sp_num)
{
    int ind_num = 0;
    vector<int> index(sp_num);
    for(int i = 0; i<sp_num; ++i)
    {
        index[i] = 0;
    }
    for(int i = 0; i<sp_num; ++i)
    {
        if(W[i]!=num)
        {
            index[ind_num] = i;
            ++ind_num;
        }
    }

    return ind_num;
}

vector<float>  solve_tv_iter(const double *W, double *initial_u, int sp_num, int iter_num,  float lambda)
{
    // initial
    vector<float> new_u(sp_num);
    for(int i = 0; i<sp_num; ++i)
    {
        new_u[i] = initial_u[i];
    }

    vector< vector<float> >  aff_matrix(sp_num, vector<float>(sp_num));

    for(int i = 0; i<sp_num*sp_num; ++i)
    {
        aff_matrix[i/sp_num][i%sp_num] = W[i];
    }

    for (int iter_index = 0; iter_index < iter_num; ++iter_index)
    {
        for(int i = 0; i < sp_num; ++i)
        {
            vector<int> find_num;
            find_num = find_index_num(aff_matrix[i],0,sp_num);
            int size = index_num(aff_matrix[i],0,sp_num);
            float sum_w = 0;
            float sum_j = 0;
            for(int j = 0; j<size; ++j)
            {
                sum_j += aff_matrix[i][find_num[j]]*new_u[find_num[j]];
                sum_w += aff_matrix[i][find_num[j]];
            }
            float sum_tmp = 0;
            sum_tmp = sum_j - sum_w*new_u[i];
            new_u[i] = new_u[i] + sum_tmp - lambda * (new_u[i] - initial_u[i]);
        }

    }

    return new_u;

}
void mexFunction(int numOut, mxArray *pmxOut[],
    int numIn, const mxArray *pmxIn[])
{
    double *W = mxGetPr(pmxIn[0]);
    double *initial_u = mxGetPr(pmxIn[1]);
    double *sp_num = mxGetPr(pmxIn[2]);
    double *iter_num = mxGetPr(pmxIn[3]);
    double *lambda = mxGetPr(pmxIn[4]);

    //int num = sp_num[0];
    vector<float> new_u = solve_tv_iter(W, initial_u, sp_num[0],  iter_num[0], lambda[0]);

    pmxOut[0] = mxCreateDoubleMatrix(1,(int)sp_num[0],mxREAL);//num
    double *solved_u;
    solved_u = mxGetPr(pmxOut[0]);
    for(int i = 0; i<sp_num[0];++i)
    {
        solved_u[i] = new_u[i];
    }

}
