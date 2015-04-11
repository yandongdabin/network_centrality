#ifndef CENTRALITY_H
#define CENTRALITY_H

/*幂方法求解方阵的最大特征值及其对应的特征向量*/
void eigenvector_centrality(double **A,int N);

/*degree-centrality 求每个节点的度*/
void degree_centrality(double **A,int N);

/*katz 方法求centrality*/
void katz_centrality(double **A,int N,double a,double b);

/*pagerank centrality*/
void pagerank_centrality(double **A,int N,double a,double b);

/*用floyd算法来求betweeness centrality*/
void betweeness_centrality(double **A,int N);

/*用floyd算法来求closeness centrality*/
void closeness_centrality(double **A,int N);
#endif // CENTRALITY_H
