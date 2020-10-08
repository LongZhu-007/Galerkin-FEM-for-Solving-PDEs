/* 1D微分方程的求解 */
#include "Common.h"
int Gmesh1D(double *meshes, int N, double x_min, double x_max){
    if (x_min >= x_max){
        printf("Sorry! Parameters invalid: x_min >= x_max.\n");
        exit(1);
    }
    double dh = (x_max - x_min) / N;
    int i;
    for (i=0; i<N+1; i++){
        *(meshes+i) = x_min + dh * i;
    }
    return 0;
}

double hat01(double x, double x0, double x1){
    if (x < x0) return 0.0;
    else if (x > x1) return 0.0;
    return (x - x0)/(x1 - x0);
};

double hat02(double x, double x1, double x2){
    if (x < x1) return 0.0;
    else if (x > x2) return x2;
    return (x2 - x) / (x2 - x1);
}

double hat(double x, double x0, double x1, char d){
    if (d == 1) return hat01(x, x0, x1);
    else if (d == 2) return hat02(x, x0, x1);
    printf("Sorry! Parameters invalid: d.\n");
    exit(1);
}

double D_hat01(double x, double x0, double x1){
    if (x < x0) return 0.0;
    else if (x > x1) return 0.0;
    return 1/(x1 - x0);
}

double D_hat02(double x, double x1, double x2){
    if (x < x1) return 0.0;
    else if (x > x2) return 0.0;
    return -1/(x2 - x1);
}

double D_hat(double x, double x0, double x1, char d){
    if (d == 1) return D_hat01(x, x0, x1);
    else if (d == 2) return D_hat02(x, x0, x1);
    printf("Sorry! Parameters invalid: d.\n");
    exit(1);
}

double init_D_hats(double x0, double x1, double (*p)(double), char d1, char d2){
    double c00 = 0.0;
    c00 += (*p)(x0)*D_hat(x0, x0,x1, d1)*D_hat(x0, x0, x1, d2);
    c00 += (*p)(x1)*D_hat(x1, x0, x1, d1)*D_hat(x1, x0,x1, d2);
    c00 *= (x1 - x0) / 2;
    return c00;
}

double init_hats(double x0, double x1, double (*q)(double), char d1, char d2){
    double c00 = 0.0;
    c00 += (*q)(x0)*hat(x0, x0, x1, d1)*hat(x0, x0,x1, d2);
    c00 += (*q)(x1)*hat(x1, x0, x1, d1)*hat(x1, x0, x1, d2);
    c00 *= (x1 - x0) / 2;
    return c00;
}

double init_f_hat(double x0, double x1, double (*f)(double), char d){
    double c00 = 0.0;
    c00 += (*f)(x0)*hat(x0, x0, x1, d);
    c00 += (*f)(x1)*hat(x1, x0, x1, d);
    c00 *= (x1 - x0)/2;
    return c00;
}

int Stiff(double *KK, int N, double *meshes, double (*p)(double), double (*q)(double)){
    /* 网格节点坐标数据 */
    /* 暂时不添加边界条件, 然后刚度矩阵的shape (N+1)*(N+1) */
    int i;
    for (i = 0; i<N; i++){
        *(KK+i*(N+1)+i) += init_D_hats(meshes[i], meshes[i+1], p, 2, 2)+init_hats(meshes[i], meshes[i+1], q, 2, 2);
        *(KK+i*(N+1)+i+1) += init_D_hats(meshes[i], meshes[i+1], p, 2, 1)+init_hats(meshes[i], meshes[i+1], q, 2, 1);
        *(KK+(i+1)*(N+1)+i) += init_D_hats(meshes[i], meshes[i+1], p, 1, 2)+init_hats(meshes[i], meshes[i+1], q, 1, 2);
        *(KK+(i+1)*(N+1)+i+1) += init_D_hats(meshes[i], meshes[i+1], p, 1, 1)+init_hats(meshes[i], meshes[i+1], q, 1, 1);
    }
    return 0;
}

/* 然后我们 来创建一个函数，生成方程组的右端项 */
int right_b(double *F, int N, double *meshes, double (*f)(double)){
    int i;
    for (i = 0; i < N; i++){
        *(F+i) += init_f_hat(meshes[i], meshes[i+1], f, 2);
        *(F+i) += init_f_hat(meshes[i], meshes[i+1], f, 1);
    }
    return 0;
}

/* 然后创建边界元信息 */
int Boundary_element1(double *KK, double *meshes, int N, int *node, int m){
    /* 对于边界元，我们要给出对应的边界元的节点信息 */
    /* 给出对应的节点信息 */
    /* 那么我们接下来给 U中添加固定的边界或者其他边界条件 */
    /* 可是 边界条件是否会作用到 F */
    /*在哪个节点上的设置了1D边界，那么我们便需要来计算出边界元对应的刚度矩阵 */
    double *eK;
    eK = (double *)calloc(4, sizeof(double));
    if (eK == NULL){
        printf("Sorry! 开辟内存失败!\n");
        exit(1);
    }
    /* 这里的边界元的节点， ua = u0; ub = u1 */
    /* 那么第一个边界元是 Element 1; 第二个边界元是 Element 2*/
    *(eK+0) = 0;
    // *(eK+1) = 1/(meshes[1] - meshes[0]);
    *(eK+1) = 0;
    *(eK+2) = 0;
    *(eK+3) = 0;
    *(KK+0) += *(eK+0);
    *(KK+1) += *(eK+1);
    *(KK+(N+1)+0) += *(eK+2);
    *(KK+(N+1)+1) += *(eK+3);
    /* 下面来处理对应的最有一个节点 */
    *(eK+0) = 0;
    *(eK+1) = 0;
    // *(eK+2) = -1/(meshes[N] - meshes[N-1]);
    *(eK+2) = 0;
    *(eK+3) = 0;
    /* 然后添加到刚度矩阵上去 */
    *(KK+(N-1)*(N+1)+(N-1)) += *(eK+0);
    *(KK+(N-1)*(N+1)+N) += *(eK+1);
    *(KK+N*(N+1)+N-1) += *(eK+2);
    *(KK+N*(N+1)+N) += *(eK+3);
    /* 那么当前的边界条件便被作用到刚度矩阵上去了 */
    /* 下面我们划掉已知的位移数据 */
    /* 节点0 对应的行数和列数； 节点N对应的行数和列数 */

    return 0;
}

/* 创建一个函数，进行迭代求解 */
int Solver(double *KK, double *F, int N, int *nodes, int m, double *Bval, double *ui){
    /* 这里的 KK 是系数矩阵Shape=(N+1, N+1)， F 是方程组的右端项(N+1, 1) */
    /* 这里的 N 是要划分的网格数量 */
    /* 这里的 nodes 是指定的参与到边界条件计算的节点编号 */
    /* 这里的 m 是边界节点的个数 */
    /* 这里的ui 是存储了初始数组的数据 Shape = (N+1, 1) */
    /* 首先来删除掉特定行 列上的数据 */
    /* 在进行调整参数数据之前，先来计算出对应的方程组的右端项的新值 */
    int i, j, k;
    double val = 0.0;
    i = 0;
    k = 0;
    double *Fd;
    Fd = (double *)calloc(N+1-m, sizeof(double));
    if (Fd == NULL){
        printf("SOrry! 开辟内存失败!\n");
        exit(1);
    }
    k = 0;
    for (int x=0; x<N+1; x++){
        /* 最外层应该遍历的是有效数据行数 */
        /* 这里取 i 表示 nodes 中的索引 */
        /* 遍历列数， */
        if (i < m && x == *(nodes+i) ){
            /* 转到下一行 */
            i ++;
            continue;
        }
        /* 如果当前是有效行数的话，那么接下来该怎么操作 */
        val = 0.0;
        for (j=0; j < m; j++){
            val += *(KK+x*(N+1)+nodes[j]) * Bval[j];
            /* 作用到F上去 */
        }
        Fd[k++] = val;
    }

    /*
        给出 i 表示数据空位;
             j 表示 有效数据位置
             k 表示 nodes 的遍历索引
    */
    if (m != 0 &&  nodes != NULL){
        k = 0;
        i = *(nodes+k);
        j = i + 1;
        while (j == *(nodes+k+1)){
            j ++;
            k ++;
        }
        k ++;
        /* 然后开始复制 */
        while (j < N+1 && k < m){
            /* memcpy(KK+i, KK+j, (N+1)*sizeof(double)); */
            for (int s=0; s<N+1; s++){
                *(KK+i*(N+1)+s) = *(KK+j*(N+1)+s);
            }
            j ++;
            i ++;
            /* 然后判断对应的j 是否是有效数据 */
            while (k < m && j == *(nodes+k)){
                j ++;
                k ++;
            }
            /* 然后这里的 j 可能超出了N */
        }
        printf("打印输出交换行之后的矩阵:\n");
        print_d_Matrix(KK, N+1, N+1);
        /* 那么到这里对应的行数便交换完毕了， 下面来交换对应的列数 */
        for (int s=0; s<N+1-m; s++){
            /* 这里的 s 是行数, 下面来遍历列数  */
            k = 0;
            i = *(nodes+k);
            j = i + 1;
            while ( j == *(nodes+k+1)){
                j ++;
                k ++;
            }
            k ++;
            while (j < N+1 && k < m){
                *(KK+s*(N+1)+i) = *(KK+s*(N+1)+j);
                i ++;
                j ++;
                while (k < m && j == *(nodes+k)){
                    j ++;
                    k ++;
                }
            }
        }
        /* 下面来划分对应的F */
        k = 0;
        i = *(nodes+k);
        j = i + 1;
        while (j == *(nodes+k+1)){
            j ++;
            k ++;
        }
        k ++;
        while (j < N+1 && k < m){
            *(F+i) = *(F+j);
            i ++;
            j ++;
            while (k < m && j == *(nodes+k)){
                j ++;
                k ++;
            }
        }
        /* 循环结束之后，那么当前的遍历便结束了 */
    }
    printf("打印输出Fd 数组:\n");
    print_d_Array(Fd, N+1-m);
    printf("将Fd作用到F上去:\n");
    for (i = 0; i<N+1-m; i++){
        F[i] -= Fd[i];
    }
    printf("打印输出划掉已知行和列的矩阵：\n");
    print_d_Matrix(KK, N+1, N+1);
    printf("打印输出方程组的右端项:\n");
    print_d_Array(F, N+1);
    int Size = N + 1 - m;
    double *rm;
    rm = (double *)calloc(Size, sizeof(double));
    if (rm == NULL){
        printf("Sorry! 开辟内存失败!\n");
        exit(1);
    }
    double linf, omega = 0.6, sd;
    int itr = 0;
    /* Multipy(rm, KK, ui, Size); */
    for (int y=0; y<Size; y++){
            sd = 0.0;
        for (int z=0; z<Size; z++){
        sd += KK[y*(N+1)+z] * ui[z];
        }
        rm[y] = sd;
    }
    /* 然后计算出误差 */
    for (i=0; i<Size; i++){
        rm[i] = F[i] - rm[i];
    }
    /* 计算出初始化的误差 */
    linf = linfnorm(rm, Size);
    while (itr <= itr_max && linf >= eps){
        /* 然后首先需要把 rm 作用到 ui */
        for (i=0; i<Size; i++){
            ui[i] = ui[i] + omega * rm[i] / KK[i*(N+1)+i];
        }

        /* Multipy(rm, KK, ui, Size); */
        for (int y=0; y<Size; y++){
            sd = 0.0;
            for (int z=0; z<Size; z++){
                sd += KK[y*(N+1)+z] * ui[z];
            }
            rm[y] = sd;
        }
        for (i=0; i<Size; i++){
            rm[i] = F[i] - rm[i];
        }
        linf = linfnorm(rm, Size);
        itr ++;

    }
    printf("我们再回乘系数矩阵:\n");
    for (int y=0; y<Size; y++){
            sd = 0;
        for (int z=0; z<Size; z++){
            sd += KK[y*(N+1)+z] * ui[z];
        }
        rm[y] = sd;
    }
    printf("result rm:\n");
    print_d_Array(rm, Size);
    /* 循环结束之后， 下面我们便得到了对应的方程组的解 */
    return 0;

}

/* 下面来求解出对应的节点上的值 */
int UonNodes(double *res, double *meshes, int N, double *us, int *(nodes), int m, double *(val)){
    /* 这里的 res的长度是 N+1, 保存最后节点上的数据  */
    /* meshes 的长度是 N+1; us 的长度是 N+1-m; nodes的长度为 m; val 的长度是 N+1-m */
    int i, j, k;
    k = 0;
    j = 0;
    for (i = 0; i<N+1; i++){
        if (i == *(nodes+k)){
            *(res+i) = *(val+k);
            k ++;
        }else{
            *(res+i) = *(us+j);
            j ++;
        }
    }
    /* 循环结束之后，得到对应的 res */
    return 0;
}


/* 创建一个函数，进行方程组的求解 */
int Multipy(double *us, double *A, double *b, int Size){
    int i, j;
    double c0;
    for (i = 0; i < Size; i++){
        c0 = 0.0;
        for (j=0; j < Size; j++){
            c0 += A[i*Size+j] * b[j];
        }
        us[i] = c0;
    }
    /* 循环结束之后，便得到了对应的相乘之后的结果 */
    return 0;
}

/* 创建一个函数，计算出对应的 */

double linfnorm(double *v, int Size){
    int i;
    double linf = 0.0;
    for (i = 0; i < Size; i++){
        if (linf < fabs(v[i])){
            linf = fabs(v[i]);
        }
    }
    return linf;
}

/* 创建函数，来打印输出矩阵 */
void print_d_Array(double *Arr, int N){
    int i;
    for (i=0; i<N; i++){
        printf("%.4g\t", Arr[i]);
    }
    printf("\n");
}

void print_d_Matrix(double *Mat, int row, int col){
    int r, c;
    for (r = 0; r < row; r++){
        for (c = 0; c < col; c++){
            printf("%.4g\t", Mat[r*col+c]);
        }
        printf("\n");
    }
}
