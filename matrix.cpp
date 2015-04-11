#include "matrix.h"
#include <cstring>
#include <cmath>
#include <iostream>
#include "mem.h"
using namespace std;
const int MAX = 50;

/*求取向量v的最大范数 即绝对值最大的分量*/
double max_value(double * v, int N)
{
    double max = 0;
    for(int i=0;i<N;i++)
    {
        double tmp = fabs(v[i]);
        if(tmp > max)
            max = tmp;
    }
    return max;
}


/*进行矩阵的相减*/
void minus_matrix(double **A,double **B,int size)
{
    for(int i = 0;i<size;i++)
    {
        for(int j = 0;j<size;j++)
        {
            A[i][j]-=B[i][j];
        }

    }
}
/*计算矩阵的转置*/
void trans(double **A,int size)
{
    //double *tmp = new double*[size];
    //double tmp[MAX][MAX];
    double **tmp = get_mem(size);
    for(int i =0;i<size;i++)
        for(int j = 0;j<size;j++)
        {
            tmp[i][j] = A[j][i];
        }
    for(int i =0;i<size;i++)
        for(int j = 0;j<size;j++)
        {
            A[i][j] = tmp[i][j];
        }
    free_mem(tmp,size);

}
void identity_matrix(double **A,int size)
{
    //memset(A,0,sizeof(A));
    for(int i = 0;i<size;i++)
    {
       memset(A[i],0,sizeof(double)*size);
    }
    for(int i = 0;i<size;i++)
    {
        A[i][i] = 1;
    }

}
/*进行矩阵相乘   V = A×u*/
void multiple(double **A,double *u,double *v,int size)
{
    for(int i=0;i<size;i++)
    {
        v[i] = 0;
        for(int j = 0;j<size;j++)
        {
            v[i] += A[i][j] * u[j];
        }

    }
}

/*矩阵和一个常数相乘*/
void multiple_const(double **A,int size,double var)
{
    for(int i = 0;i<size;i++)
        for(int j =0; j<size;j++)
            A[i][j] *= var;

}
/*两个矩阵相乘*/
void multiple_matrix(double **A,double **B,double **C,int size)
{
    for(int i=0;i<size;i++)
        for(int j=0;j<size;j++)
            for(int k = 0;k<size;k++)
                C[i][j] = A[i][k] * B[k][j];

}
/*求矩阵的行列式 默认按照第一行展开*/
double det(double **A,int size)
{

    if(size == 1)
        return A[0][0];

    double **tmp = get_mem(size-1);

    double ret = 0.0;
    for(int i=0;i<size;i++)
    {
        /*构造余子式*/
        for(int j = 0;j<size-1;j++)
        {
            for(int k = 0;k<size-1;k++)
            {
                int x = j+1;
                int y = k>=i?k+1:k;
                tmp[j][k] = A[x][y];
            }
        }
        ret+=A[0][i] * (i%2==1?-1:1)*det(tmp,size-1);

    }
    free_mem(tmp,size-1);
    return ret;
}
/*求解伴随矩阵*/
void cal_adjoint(double **A,double **adj,int size)
{
    double **tmp = get_mem(size-1);
    for(int i = 0;i<size;i++)
    {
        for(int j = 0;j<size;j++)
        {
            for(int k = 0; k<size-1;k++)
                {
                    for(int m = 0;m<size-1;m++)
                    {
                        int x = k>=i?k+1:k;
                        int y = m>=j?m+1:m;
                        tmp[k][m] = A[x][y];

                    }
                }
            adj[j][i] = ((i+j)%2==1?-1:1) * det(tmp,size-1);

        }
    }
    free_mem(tmp,size-1);
}

/*求矩阵的逆*/
void inverse(double **A,double **inver,int size)
{
    double det_value = det(A,size);

    if(fabs(det_value-0.0)<=0.0000001)
    {
        std::cout<<"cannot inverse"<<std::endl;
    }

    //double adj[MAX][MAX];
    double **adj = get_mem(size);
    cal_adjoint(A,adj,size);
    for(int i = 0;i<size;i++)
    {
        for(int j = 0;j<size;j++)
        {
            inver[i][j]= adj[i][j] / det_value;
        }
    }
    free_mem(adj,size);
}

/*输出矩阵*/
void output_matrix(double **A,int size)
{


  for(int i=0;i<size;i++)
    {
        for(int j=0;j<size;j++)
                {
                    std::cout<<A[i][j]<<" ";
                }
    std::cout<<std::endl;
    }
    std::cout<<std::endl;

}

/*根据矩阵A得到对角矩阵，矩阵上的元素对应的是每一个节点的度*/
void get_diag(double **A,double **diag,int size)
{
    for(int i =0;i<size;i++)
        memset(diag[i],0,sizeof(double)*size);
    for(int i=0;i<size;i++)
    {
        int cnt = 0;
        for(int j = 0;j<size;j++)
        {
            if(fabs(A[i][j]-1.0) <0.00001)
                cnt++;
        }
        diag[i][i] = cnt;
    }
}

