#include <iostream>
#include <list>
#include <vector>
#include <cmath>
#include <mex.h>
using namespace std;
#define max(a, b)  (((a) > (b)) ? (a) : (b))

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

float sum(vector<float> &W, vector<int> find_num, int size)
{
    float sum_num = 0;
    for(int i = 0; i<size; ++i)
    {
        sum_num = sum_num+W[find_num[i]];
    }
    return sum_num;
}

vector<float> sqrt_matrix(vector<float> &u, int sp_num)
{
    vector<float> sqrt_u (sp_num);
    for (int i = 0; i<sp_num; ++i)
    {
        sqrt_u[i] = sqrt(u[i]);
    }
    return sqrt_u;
}


vector<float> solve_tv_c(const double *W, double *initial_u, int sp_num, int iter_num, int iter_u_num, float beta, float lambda, int d_mode)
{
    // initial
    vector<float> new_u(sp_num);
    for(int i = 0; i<sp_num; ++i)
    {
        new_u[i] = initial_u[i];
    }
    vector< vector<float> >  new_b(sp_num, vector<float>(sp_num));
    vector< vector<float> >  new_d(sp_num, vector<float>(sp_num));
    vector< vector<float> >  aff_matrix(sp_num, vector<float>(sp_num));
    vector<float> cost_function(iter_num);
    for(int i = 0; i<sp_num*sp_num; ++i)
    {
        aff_matrix[i/sp_num][i%sp_num] = W[i];
    }
    for(int iter_k = 0; iter_k<iter_num; ++iter_k)
    {
        //u
        for(int iter_u = 0; iter_u< iter_u_num; ++iter_u)
        {
            for(int i = 0; i<sp_num; ++i)
            {
                vector<int> find_num;
                find_num = find_index_num(aff_matrix[i],0,sp_num);
                int size = index_num(aff_matrix[i],0,sp_num);
                float sumw = sum(aff_matrix[i],find_num,size);
                float fDen = 1/(lambda+beta*sumw);
                vector<float> aff_tmp(sp_num);
                for (int j = 0; j<sp_num; ++j)
                {
                    aff_tmp[j] = aff_matrix[i][j]*new_u[j];
                }
                float sumwu = sum(aff_tmp,find_num,size);
                aff_tmp = sqrt_matrix(aff_matrix[i],sp_num);
                float sumdb = 0;
                for(int j = 0; j<size; ++j)
                {
                    sumdb += sqrt(aff_matrix[i][find_num[j]])* (new_d[i][find_num[j]] - new_d[find_num[j]][i] - new_b[i][find_num[j]]+new_b[find_num[j]][i]);
                }
                new_u[i] = fDen*(sumwu*beta + lambda*initial_u[i] - beta/2*sumdb);

            }
        }

        //d
        if (d_mode==1)

            for(int i = 0; i<sp_num; ++i)
            {
                vector<int> find_num;
                find_num = find_index_num(aff_matrix[i],0,sp_num);
                int size = index_num(aff_matrix[i],0,sp_num);

                float max_value_tmp = 0, max_value = 0;
                float sumf1 = 0, sumf2 = 0;
                for(int j = 0; j<size; ++j)
                {
                    //max_value_tmp = max_value_tmp+ aff_matrix[i][find_num[j]]*(new_u[find_num[j]] - new_u[i])*(new_u[find_num[j]] - new_u[i]);
                    //max_value_tmp = max_value_tmp+ aff_matrix[i][find_num[j]]*(new_u[find_num[j]] - new_u[i])*(new_u[find_num[j]] - new_u[i])+ new_b[i][find_num[j]]*new_b[i][find_num[j]];
                    max_value_tmp = max_value_tmp+ (sqrt(aff_matrix[i][find_num[j]])*(new_u[find_num[j]] - new_u[i])+ new_b[i][find_num[j]])*(sqrt(aff_matrix[i][find_num[j]])*(new_u[find_num[j]] - new_u[i])+ new_b[i][find_num[j]]);
                }
                for(int j = 0; j<size; ++j)
                {
                    //sumf2 =  sqrt(max_value_tmp + new_b[i][find_num[j]]*new_b[i][find_num[j]]);
                    sumf2 =  sqrt(max_value_tmp);
                    max_value = sumf2 - beta;
                    max_value = max(max_value,0);
                    sumf1 = sqrt(aff_matrix[i][find_num[j]])*(new_u[find_num[j]] - new_u[i])+ new_b[i][find_num[j]];
		    if (sumf2!=0)
                        new_d[i][find_num[j]] = sumf1/sumf2*max_value;
                }
            }
        else
        {
            // 2����
            for(int i = 0; i<sp_num; ++i)
            {
                vector<int> find_num;
                find_num = find_index_num(aff_matrix[i],0,sp_num);
                int size = index_num(aff_matrix[i],0,sp_num);
                float deta_u = 0;
                for (int j = 0; j<size; ++j)
                {
                    deta_u = sqrt(aff_matrix[i][find_num[j]])* (new_u[find_num[j]] - new_u[i])+new_b[i][find_num[j]];
                    new_d[i][find_num[j]] = beta/(2+beta)*deta_u;
                }
            }
        }

        //b
        for (int i = 0; i<sp_num; ++i)
        {
            vector<int> find_num;

            find_num = find_index_num(aff_matrix[i],0,sp_num);
            int size = index_num(aff_matrix[i],0,sp_num);
            for (int j = 0; j<size; ++j)
            {
                new_b[i][find_num[j]] = new_b[i][find_num[j]] +  sqrt(aff_matrix[i][find_num[j]])*(new_u[find_num[j]] - new_u[i]) - new_d[i][find_num[j]];
            }
        }

        // cost function

        for (int i = 0; i<sp_num ; ++i)
        {
            vector<int> find_num;

            find_num = find_index_num(aff_matrix[i],0,sp_num);
            int size = index_num(aff_matrix[i],0,sp_num);

            float max_value_tmp = 0;
            float sum1 = 0, sum2 = 0;
            for(int j = 0; j<size; ++j)
            {
                max_value_tmp = max_value_tmp+ aff_matrix[i][find_num[j]]*(new_u[find_num[j]] - new_u[i])*(new_u[find_num[j]] - new_u[i]);
            }
            sum1 = sqrt(max_value_tmp);
            sum2 = lambda/2*(new_u[i] - initial_u[i]*initial_u[i]);
            cost_function[iter_k] = cost_function[iter_k]+sum1+sum2;
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
    double *iter_u_num = mxGetPr(pmxIn[4]);
    double *beta = mxGetPr(pmxIn[5]);
    double *lambda = mxGetPr(pmxIn[6]);
    double *d_mode = mxGetPr(pmxIn[7]);
    //int num = sp_num[0];
    vector<float> new_u = solve_tv_c(W, initial_u, sp_num[0],  iter_num[0],  iter_u_num[0], beta[0],  lambda[0],  d_mode[0]);

    pmxOut[0] = mxCreateDoubleMatrix(1,(int)sp_num[0],mxREAL);//num
    double *solved_u;
    solved_u = mxGetPr(pmxOut[0]);
    for(int i = 0; i<sp_num[0]; ++i)
    {
        solved_u[i] = new_u[i];
    }

}
