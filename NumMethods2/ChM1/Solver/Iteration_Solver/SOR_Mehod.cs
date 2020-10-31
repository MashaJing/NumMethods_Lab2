using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Com_Methods
{
    class SOR_Method : IIteration_Solver
    {
        //максимальное число итераций
        public int Max_Iter { set; get; }

        //точность решения
        public double Eps { set; get; }

        //текущая итерация
        public int Iter { set; get; }

        public SOR_Method(int MAX_ITER, double EPS)
        {
            Max_Iter = MAX_ITER;
            Eps = EPS;
            Iter = 0;
        }

        //Точечный метод релаксации: реализация
        public Vector Start_Solver(Matrix A, Vector F, double w)
        {
            //норма разности решения на двух итерациях
            double Norm_Xnew_Xold;

            //инициализация вектора результата
            var RES = new Vector(F.N);
            var RES_New = new Vector(F.N);

            //начальное приближение
            for (int k = 0; k < RES.N; k++) RES.Elem[k] = 0.0;

            //итерации метода блочной релаксации
            do
            {
                Norm_Xnew_Xold = 0.0;

                for (int i = 0; i < RES.N; i++)
                {

                    //инициализация скобки (F - Ax)
                    double F_Ax = F.Elem[i];

                    //произведение матрицы на старые Х (нижний треугольник)
                    for (int j = 0; j < i; j++)
                    {
                        F_Ax -= A.Elem[i][j] * RES_New.Elem[j];
                    }

                    //произведение матрицы на старые Х (верхний треугольник)
                    for (int j = i + 1; j < A.N; j++)
                    {
                        F_Ax -= A.Elem[i][j] * RES.Elem[j];
                    }

                    //формируем результат для 1-ой компоненты и квадраты норм

                    F_Ax /= A.Elem[i][i];
                    RES_New.Elem[i] = (1.0 - w)*RES.Elem[i] + F_Ax * w;


                    Norm_Xnew_Xold += Math.Pow((RES.Elem[i] - RES_New.Elem[i]), 2);
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


        //блочный метод релаксации
        public Block_Vector Start_Solver(Block_Matrix A, Block_Vector F, double w)
        {
            //для вычисления норм ||Xnew - Xold||
            double Norm_Xnew_Xold;

            //инициализация вктора результата
            var RES = new Block_Vector(A.N, A.Size_Block);
            var RES_New = new Block_Vector(A.N, A.Size_Block);

            //LU-решатель для обращения диагональных элементов
            var LU_Solver = new LU_Decomposition();

            //Вспомогательный вектор для вычисления скобки (F_Ax)
            var F_Ax = new Vector(A.Size_Block);

            //итерации метода блочной релаксации
            do
            {
                Norm_Xnew_Xold = 0.0;
                //цикл по неизвестным
                for (int i=0; i<RES.N; i++)
                {
                    //инициализация суммы
                    for (int k = 0; k < RES.N; k++)
                        F_Ax.Elem[k] = F.Block[i].Elem[k];

                    //Произведение блоков матрицы на новые Х

                    for (int j = 0; j <  i; j++)
                    {
                        var Current_Matrix_Block = A.Block[i][j];
                        var Current_Vector_Block = RES_New.Block[j]; //RES_NEW????
                        for (int Row = 0; Row < Current_Matrix_Block.M; Row++)
                            for (int Column = 0; Column < Current_Matrix_Block.N; Column++)
                                F_Ax.Elem[Row] -= Current_Matrix_Block.Elem[Row][Column] * Current_Vector_Block.Elem[Column];
                    }

                    //Произведение блоков матрицы на предыдущие Х

                    for (int j = i+1; j < A.N; j++)
                    {
                        var Current_Matrix_Block = A.Block[i][j];
                        var Current_Vector_Block = RES.Block[j];
                        for (int Row = 0; Row < Current_Matrix_Block.M; Row++)
                            for (int Column = 0; Column < Current_Matrix_Block.N; Column++)
                                F_Ax.Elem[Row] -= Current_Matrix_Block.Elem[Row][Column] * Current_Vector_Block.Elem[Column];
                    }

                    //решение СЛАУ с диагональной матрицей
                    LU_Solver.LU = A.Block[i][i];
                    F_Ax = LU_Solver.Start_Solver(F_Ax);

                    //формируем результат для i-ой компоненты
                    for (int k = 0; k < RES.Block_Size; k++)
                    {
                        double X_NEW = (1.0 - w) * RES.Block[i].Elem[k] + F_Ax.Elem[k] * w;
                        Norm_Xnew_Xold += (RES.Block[i].Elem[k] - X_NEW) * (RES.Block[i].Elem[k] - X_NEW);
                        RES.Block[i].Elem[k] = X_NEW;

                    }

                }
            }
            while (Norm_Xnew_Xold > Eps && Iter < Max_Iter);

            return RES;
        }
    }
}
