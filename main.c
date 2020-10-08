#include "Common.h"

int main()
{
    printf("Hello world!\n");
    double x_min = 0;
    double x_max = 1;
    int N = 10;
    /* Ȼ�󿪱����� */
    double *meshes, *us, *KK, *F;
    meshes = (double *)calloc(N+1, sizeof(double));
    us = (double *)calloc(N+1, sizeof(double));
    KK = (double *)calloc((N+1)*(N+1), sizeof(double));
    F = (double *)calloc(N+1, sizeof(double));
    if (meshes == NULL || us == NULL || KK == NULL || F == NULL){
        printf("SOrry! �����ڴ�ռ�ʧ��!\n");
        exit(1);
    }
    /* ���Ƚ������񻮷� */
    Gmesh1D(meshes, N, x_min, x_max);
    printf("��ӡ�����ǰ������ڵ���������:\n");
    print_d_Array(meshes, N+1);
    /* ���ɸնȾ��� */
    Stiff(KK, N, meshes, p, q);
    printf("��ӡ����նȾ���:\n");
    print_d_Matrix(KK, N+1, N+1);
    right_b(F, N, meshes, f);
    printf("��ӡ�����������Ҷ���:\n");
    print_d_Array(F, N+1);;
    /*��ӵ�һ��߽����� */
    int nodes[2] = {0, N};
    double Bvals[2] = {100, 10};
    int m = 2;
    Boundary_element1(KK, meshes, N, nodes, m);
    printf("��ӡ�������˱߽�����֮���ϵ������:\n");
    print_d_Matrix(KK, N+1, N+1);
    /* ���������� */
    Solver(KK, F, N, nodes, m, Bvals, us);
    printf("��ӡ��������������:\n");
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
