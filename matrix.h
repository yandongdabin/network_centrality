#ifndef MATRIX_H
#define MATRIX_H

/*求取向量v的最大范数 即绝对值最大的分量*/
double max_value(double * v, int N);

/*进行矩阵的相减*/
void minus_matrix(double **A,double **B,int size);

/*计算矩阵的转置*/
void trans(double **A,int size);

/*获得单位矩阵*/
void identity_matrix(double **A,int size);

/*进行矩阵相乘   V = A×u*/
void multiple(double **A,double *u,double *v,int size);

/*矩阵和一个常数相乘*/
void multiple_const(double **A,int size,double var);

/*两个矩阵相乘*/
void multiple_matrix(double **A,double **B,double **C,int size);

/*求矩阵的行列式 默认按照第一行展开*/
double det(double **A,int size);

/*求解伴随矩阵*/
void cal_adjoint(double **A,double **adj,int size);

/*求矩阵的逆*/
void inverse(double **A,double **inver,int size);

/*输出矩阵*/
void output_matrix(double **A,int size);

/*根据矩阵A得到对角矩阵，矩阵上的元素对应的是每一个节点的度*/
void get_diag(double **A,double **diag,int size);


#endif // MATRIX_H
