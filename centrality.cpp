#include "centrality.h"
#include "matrix.h"
#include <cmath>
#include <iostream>
#include "mem.h"
#include <algorithm>
#include <numeric>
#include <set>
#include <iterator>
using namespace std;

const int SUP = 1000;


/*幂方法求解方阵的最大特征值及其对应的特征向量*/
void eigenvector_centrality(double **A,int N)
{
    double *v = new double[N];
    for(int i=0;i<N;i++)
        v[i] = 1.0;

    double *u = new double[N];
    int step = 1;
    while(true)
    {
        double max = max_value(v,N);
        for(int i=0;i<N;i++)
        {
            u[i] = v[i] / max;
        }
        //if(fabs(max_value(u,N)- max) > delta) break;
        multiple(A,u,v,N);
        step++;
        if(step > SUP) break;
    }
    std::cout<<max_value(v,N)<<std::endl;
    for(int i =0 ;i<N ;i++)
        std::cout<<u[i]<<" ";
    std::cout<<std::endl;



    delete(u);
    delete(v);

}
/*degree-centrality 求每个节点的度*/
void degree_centrality(double **A,int N)
{
    for(int i=0;i<N;i++)
    {
        int cnt = 0;
        for(int j = 0;j<N;j++)
        {
            if(fabs(A[i][j]-1.0) <0.00001)
                cnt++;
        }
        std::cout<<cnt<<" ";
    }
}
/*katz 方法求centrality*/
void katz_centrality(double **A,int N,double a,double b)
{
    //得到单位矩阵
    double **I = get_mem(N);
    identity_matrix(I,N);

    /*对A求转置并且 与a相乘*/

    trans(A,N);
    multiple_const(A,N,a);

    /*单位矩阵减去上步所求*/
    minus_matrix(I,A,N);
    double **inver = get_mem(N);
    /*求逆乘b*/
    inverse(I,inver,N);
    multiple_const(inver,N,b);
    for(int i=0;i<N;i++)
    {
        double sum = 0.0;
        for(int j =0;j<N;j++)
            sum += inver[i][j];
        std::cout<<sum<<" ";
        //std::cout<<std::accumulate((double *)inver[i],(double *)(inver[i]+N));
    }
    free_mem(inver,N);
    free_mem(I,N);

}
/*pagerank centrality*/
void pagerank_centrality(double **A,int N,double a,double b)
{
    //得到单位矩阵
    double **I = get_mem(N);
    identity_matrix(I,N);
    /*得到D 并且求其逆*/
    double **diag = get_mem(N);
    get_diag(A,diag,N);

    double **diag_inver = get_mem(N);
    inverse(diag,diag_inver,N);



    /*对A求转置并且 与a相乘*/
    trans(A,N);
    multiple_const(A,N,a);
    double **product = get_mem(N);
    multiple_matrix(A,diag_inver,product,N);


    /*单位矩阵减去上步所求*/
    minus_matrix(I,product,N);
    double **inver = get_mem(N);

    /*求逆乘b*/
    inverse(I,inver,N);
    multiple_const(inver,N,b);

    for(int i=0;i<N;i++)
    {
        double sum = 0.0;
        for(int j =0;j<N;j++)
            sum += inver[i][j];
        std::cout<<sum<<" ";
        //std::cout<<std::accumulate((double *)inver[i],(double *)(inver[i]+N));
    }
    free_mem(inver,N);
    free_mem(I,N);
    free_mem(diag,N);
    free_mem(diag_inver,N);
    free_mem(product,N);

}
/*用floyd算法来求betweeness centrality*/
void betweeness_centrality(double **A,int N)
{
    const double INF = 1000000.0;
    double **tmp = get_mem(N);
    double **path = get_mem(N);
    std::multiset<int> path_more[N][N];
    /*记录最短路径的数目*/
    double **cnt = get_mem(N);


    for(int i=0;i<N;i++)
        for(int j = 0;j<N;j++)
        {
            bool is = (fabs(A[i][j]-0) < 0.000001) ;
            tmp[i][j] =is ? INF:A[i][j];
            cnt[i][j] = is ?0:1;
            if(!is)
            {
                path_more[i][j].insert(i+1);
                path_more[i][j].insert(j+1);
            }
            path[i][j] = -1;

        }
    /*floyd算法*/
    for(int k=0;k<N;k++)
    {
        for(int i=0;i<N;i++)
        {
            for(int j=0;j<N;j++)
            {
                if(tmp[i][j] > tmp[i][k] + tmp[k][j])
                {
                    tmp[i][j] = tmp[i][k] + tmp[k][j];
                    path[i][j] = k;
                    cnt[i][j] = 1;
                    path_more[i][j].clear();
                    set_union(path_more[i][k].begin(),path_more[i][k].end(),path_more[k][j].begin(),path_more[k][j].end(),inserter(path_more[i][j],path_more[i][j].begin()));
                }
                else if(tmp[i][j] == tmp[i][k] + tmp[k][j])
                {
                    cnt[i][j] += 1;
                    set_union(path_more[i][k].begin(),path_more[i][k].end(),path_more[k][j].begin(),path_more[k][j].end(),inserter(path_more[i][j],path_more[i][j].end()));
                }
            }
        }
    }
    for(int i=0;i<N;i++)
    {
        double total = 0;
        for(int j=0;j<N;j++)
        {
            for(int k=0;k<N;k++)
            {
                if(k>=j || j==i || k == i) continue;
                double c = 1.0*path_more[j][k].count(i+1);
                if(cnt[j][k] != 0)
                total += c / cnt[j][k];
                //int cnt = 0;
                /*set<int>::iterator beg = path_more[j][k].begin();
                set<int>::iterator end = path_more[j][k].end();
                while(beg!=end)
                {
                   if(beg->find(i)!=end);
                    //if(path_more[j][k].find(i)) ++cnt;
                    ++beg;
                }*/

            }
        }
        total *= 2;
        cout<<total<<" ";

    }
    /*for(int i =0;i<N;i++)
    {
        for(int j =0;j<N;j++)
        {
            if(i>=j) continue;
            cout<<"total count is "<<cnt[i][j]<<endl;
            cout<<"from "<<i+1<<" to "<<j+1<<endl;
            set<int>::iterator it = path_more[i][j].begin();
            while(it!=path_more[i][j].end())
            {
                cout<<*it<<" ";
                ++it;
            }
            cout<<endl;


        }


    }*/

    free_mem(tmp,N);
    free_mem(path,N);
    free_mem(cnt,N);


}

/*用floyd算法来求closeness centrality*/
void closeness_centrality(double **A,int N)
{
    const double INF = 1000000.0;
    double **tmp = get_mem(N);
    double **path = get_mem(N);
    for(int i=0;i<N;i++)
        for(int j = 0;j<N;j++)
        {
            tmp[i][j] =(fabs(A[i][j]-0) < 0.000001) ? INF:A[i][j];
            path[i][j] = -1;

        }
    for(int k=0;k<N;k++)
    {
        for(int i=0;i<N;i++)
        {
            for(int j=0;j<N;j++)
            {
                if(tmp[i][j] > tmp[i][k] + tmp[k][j])
                {
                    tmp[i][j] = tmp[i][k] + tmp[k][j];
                    path[i][j] = k;
                }
            }

        }
    }
    for(int i=0;i<N;i++)
    {
        double sum = 0.0;
        for(int j=0;j<N;j++)
        {
            if(i==j) continue;
            sum+=tmp[i][j];
        }
        sum/=((N-1)*1.0);
        std::cout<<1.0/sum<<" ";


    }
    std::cout<<std::endl;

    free_mem(tmp,N);
    free_mem(path,N);


}
