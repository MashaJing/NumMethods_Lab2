using System;
using System.Threading;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Com_Methods
{
    class Jacobi_Method : IIteration_Solver 
    {
        //наибольшее число итераций
        public int Max_Iter { set; get; }

        //точность решения
        public double Eps { set; get; }

        //текущая итерация
        public int Iter { set; get; }

        public Jacobi_Method(int MAX_ITER, double EPS)
        {
            Max_Iter = MAX_ITER;
            Eps = EPS;
            Iter = 0;

        }

        //точечный метод Якоби
        public Vector Start_Solver(Matrix A, Vector F)
        {
            //для вычисления норм ||Xnew - Xold||
            double Norm_Xnew_Xold;

            //инициализация вктора результата
            var RES = new Vector(F.N);
            var RES_New = new Vector(F.N);

            for (int k = 0; k < RES.N; k++) RES.Elem[k] = 0.0;

            //итерации метода Якоби
            do
            {
                Norm_Xnew_Xold = 0.0;

                //цикл по неизвестным
                for (int i=0; i<RES.N; i++)
                {
                    //инициализация скобки (F - Ax)
                    double F_Ax = F.Elem[i];

                    //произведение матрицы на старые Х (нижний треугольник)
                    for (int j = 0; j < i; j++)
                    {
                        F_Ax -= A.Elem[i][j] * RES.Elem[j];
                    }

                    //произведение матрицы на старые Х (верхний треугольник)
                    for (int j = i + 1; j < A.N; j++)
                    {
                        F_Ax -= A.Elem[i][j] * RES.Elem[j];
                    }

                    //формируем результат для 1-ой компоненты и квадраты норм
                    RES_New.Elem[i] = F_Ax/A.Elem[i][i];

                    Norm_Xnew_Xold += ((RES.Elem[i] - RES_New.Elem[i])*(RES.Elem[i] - RES_New.Elem[i]));

                }

                //копирование полученного результата
                RES.Copy(RES_New);

                //норма абсолютной погрешности
                Norm_Xnew_Xold = Math.Sqrt(Norm_Xnew_Xold);
                Iter++;


                Console.WriteLine("Iter {0, -10} {1}", Iter, Norm_Xnew_Xold);
            }
            while (Norm_Xnew_Xold > Eps && Iter < Max_Iter);

            return RES;
        }


        public Block_Vector Start_Solver(Block_Matrix A, Block_Vector F)
        {
            //для вычисления норм ||Xnew - Xold||
            double Norm_Xnew_Xold;

            //инициализация вктора результата
            var RES = new Block_Vector(A.N, A.Size_Block);
            var RES_New = new Block_Vector(A.N, A.Size_Block);

            //LU-решатель для обращения диагональных блоков
            var LU_Solver = new LU_Decomposition();

            //вспомогательный вектор для вычисления скобки (F - Ax)
            var F_Ax = new Vector(A.Size_Block);

            //for (int k = 0; k < RES.N; k++) RES.Elem[k] = 0.0;

            //итерации метода Якоби
            do
            {
                Norm_Xnew_Xold = 0.0;

                //цикл по неизвестным
                for (int i = 0; i < RES.N; i++)
                {
                    //инициализация скобки (F - Ax)
                    F_Ax.Copy(F.Block[i]);

                        //произведение блоков матрицы на старые Х (верхний треугольник)
                        for (int j = 0; j < i; j++)
                        {
                            var Current_Matrix_Block = A.Block[i][j];
                            var Current_Vector_Block = RES.Block[j];
                            for (int Row = 0; Row < Current_Matrix_Block.M; Row++)
                                for (int Column = 0; Column < Current_Matrix_Block.N; Column++)
                                    F_Ax.Elem[Row] -= Current_Matrix_Block.Elem[Row][Column] * Current_Vector_Block.Elem[Column];
                        }
                    
                    

                        //произведение блоков матрицы на старые Х (верхний треугольник)
                        for (int j = i + 1; j < A.N; j++)
                        {
                            var Current_Matrix_Block = A.Block[i][j];
                            var Current_Vector_Block = RES.Block[j];
                            for (int Row = 0; Row < Current_Matrix_Block.M; Row++)
                                for (int Column = 0; Column < Current_Matrix_Block.N; Column++)
                                    F_Ax.Elem[Row] -= Current_Matrix_Block.Elem[Row][Column] * Current_Vector_Block.Elem[Column];
                        }                    
                    
                    //Решение СЛАУ с диагональным блоком
                    LU_Solver.LU = A.Block[i][i];
                    F_Ax = LU_Solver.Start_Solver(F_Ax);


                    //формируем результат для i-ой компоненты и квадрат нормы разности решений
                    for (int k = 0; k < RES.Block_Size; k++)
                    {
                        RES_New.Block[i].Elem[k] = F_Ax.Elem[k];
                        Norm_Xnew_Xold += Math.Pow((RES_New.Block[i].Elem[k] - RES.Block[i].Elem[k]), 2);

                    }
                }
                //копирование полученного результата
                for (int i = 0; i < RES.N; i++)
                    RES.Block[i].Copy(RES_New.Block[i]);

                //Норма разности решений
                Norm_Xnew_Xold = Math.Sqrt(Norm_Xnew_Xold);
                Iter++;


                Console.WriteLine("Iter {0, -10} {1}", Iter, Norm_Xnew_Xold);
            }
            while (Norm_Xnew_Xold > Eps && Iter < Max_Iter);

            return RES;
        }
    }
}
