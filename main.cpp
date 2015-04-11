#include <iostream>
#include <cstring>
#include <math.h>
#include "matrix.h"
#include <fstream>
#include "centrality.h"
#include "mem.h"

const int N = 5;
const double a = 0.3;
const double b = 0.3;
using namespace std;


/*读取矩阵文件到内存*/
double ** read(int size)
{

    double ** A = get_mem(size);
    ifstream f("a.txt");
    if(f.fail())
    {
        cout<<"文件打开失败"<<endl;
        return NULL;
    }
    while(!f.eof())
    {
        for(int i=0;i<N;i++)
            for(int j=0;j<N;j++)
                {
                    f>>A[i][j];
                }

    }
    f.close();
    return A;

}
int main()
{


    double **A = read(N);
    if(A == NULL)
        return -1;
    //output_matrix(A,N);
    //eigenvector_centrality(A,N);

    //cout<<det(A,N);

    //double ** inver = get_mem(N);
    //inverse(A,inver,N);
    //output_matrix(inver,N);
    //free_mem(inver,N);

    //katz_centrality(A,N,a,b);
    //pagerank_centrality(A,N,a,b);
    //betweeness_centrality(A,N);
    //closeness_centrality(A,N);
    free_mem(A,N);
    return 0;
}
