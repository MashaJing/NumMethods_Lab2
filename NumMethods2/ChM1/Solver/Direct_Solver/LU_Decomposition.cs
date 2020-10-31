using System;

namespace Com_Methods
{
    /// <summary>
    /// класс, реализующий решатель на базе LU-разложения
    /// </summary>
    class LU_Decomposition
    {
        //матрица-хранилище для LU-разложения
        public Matrix LU { set; get; }

        //конструктор по умолчанию
        public LU_Decomposition() { }

        //реализация LU-разложения
        public LU_Decomposition(Matrix A)
        {
            //хранилище для верхней и нижней треугольных матриц
            LU = new Matrix(A.M, A.N);
            //копирование исходной матрицы
            LU.Copy(A);

            //построение верхней треугольной матрицы
            Gauss_Method.Direct_Way(LU);

            //построение нижней треугольной матрицы
            for (int i = 0; i < A.M; i++)
            {
                for (int j = 0; j < i; j++)
                {
                    double sum_LikUkj = 0;
                    for (int k = 0; k < j; k++)
                    {
                        sum_LikUkj += LU.Elem[i][k] * LU.Elem[k][j];
                    }
                    LU.Elem[i][j] = (A.Elem[i][j] - sum_LikUkj) / LU.Elem[j][j];
                }
            }
        }

        /// <summary>
        /// прямой ход (без деления на диагональ, т.к. L_ii = 1)
        /// </summary>
        void Direct_Way(Matrix A, Vector F, Vector RES)
        {
            RES.Copy(F);
            for (int i = 0; i < A.M; i++)
            {
                for (int j = 0; j < i; j++)
                {
                    RES.Elem[i] -= A.Elem[i][j] * RES.Elem[j];
                }
            }
        }

        /// <summary>
        /// решатель СЛАУ с плотной матрицей на базе LU-разложения
        /// LU разложение предполагается построенным
        /// </summary>
        public Vector Start_Solver(Vector F)
        {
            if (LU == null) throw new Exception("Error: LU-matrix is not built. Constructing of LU-matrix is required...");
            var RES = new Vector(F.N);
            //прямой ход
            Direct_Way(LU, F, RES);
            //обратный ход
            Substitution_Method.Back_Row_Substitution(LU, RES, RES);
            return RES;
        }
    }
}
