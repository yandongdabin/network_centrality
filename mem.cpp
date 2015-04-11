#include "mem.h"

/*开放一个二维数组*/
double ** get_mem(int size)
{
    double ** arr = new double *[size];
    for(int i =0;i<size;i++)
    {
        arr[i] = new double[size];
    }
    return arr;
}

/*释放数组的空间*/
void free_mem(double **A,int size)
{
    for(int i=0;i<size;i++)
    {
        delete(A[i]);
    }
    delete(A);
}
