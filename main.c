#include "Common.h"

int main()
{
    printf("Hello world!\n");
    double x_min = 0;
    double x_max = 1;
    int N = 10;
    /* 然后开辟数组 */
    double *meshes, *us, *KK, *F;
    meshes = (double *)calloc(N+1, sizeof(double));
    us = (double *)calloc(N+1, sizeof(double));
    KK = (double *)calloc((N+1)*(N+1), sizeof(double));
    F = (double *)calloc(N+1, sizeof(double));
    if (meshes == NULL || us == NULL || KK == NULL || F == NULL){
        printf("SOrry! 开辟内存空间失败!\n");
        exit(1);
    }
    /* 首先进行网格划分 */
    Gmesh1D(meshes, N, x_min, x_max);
    printf("打印输出当前的网格节点坐标数据:\n");
    print_d_Array(meshes, N+1);
    /* 生成刚度矩阵 */
    Stiff(KK, N, meshes, p, q);
    printf("打印输出刚度矩阵:\n");
    print_d_Matrix(KK, N+1, N+1);
    right_b(F, N, meshes, f);
    printf("打印输出方程组的右端项:\n");
    print_d_Array(F, N+1);;
    /*添加第一类边界条件 */
    int nodes[2] = {0, N};
    double Bvals[2] = {100, 10};
    int m = 2;
    Boundary_element1(KK, meshes, N, nodes, m);
    printf("打印输出添加了边界条件之后的系数矩阵:\n");
    print_d_Matrix(KK, N+1, N+1);
    /* 进行求解计算 */
    Solver(KK, F, N, nodes, m, Bvals, us);
    printf("打印输出方程组求解结果:\n");
    print_d_Array(us, N+1);
    return 0;
}

double p(double x){
    return 1;
}

double q(double x){
    return 0;
}

double f(double x){
    return 0;
}
