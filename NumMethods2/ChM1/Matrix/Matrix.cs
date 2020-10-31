using System;
using System.Threading;
using System.IO;

namespace Com_Methods
{
    //интерфейс матрицы
    public interface IMatrix
    {
        //размер матрицы
        int M { set; get; } //строки
        int N { set; get; } //столбцы
    }

    //класс матрицы
    public class Matrix : IMatrix
    {
        //размер матрицы
        public int M { set; get; }
        public int N { set; get; }
        //элементы матрицы
        public double[][] Elem { set; get; }

        //конструктор по умолчанию
        public Matrix() { }

        //конструктор нуль-матрицы m X n
        public Matrix(int m, int n)
        {
            N = n; M = m;
            Elem = new double[m][];
            for (int i = 0; i < m; i++) Elem[i] = new double[n];
        }

        //умножение на скаляр с выделением памяти под новую матрицу
        public static Matrix operator *(Matrix T, double Scal)
        {
            Matrix RES = new Matrix(T.M, T.N);

            for (int i = 0; i < T.M; i++)
            {
                for (int j = 0; j < T.N; j++)
                {
                    RES.Elem[i][j] = T.Elem[i][j] * Scal;
                }
            }
            return RES;
        }

        //умножение на скаляр, результат запишется в исходную матрицу
        public void Dot_Scal(double Scal)
        {
            for (int i = 0; i < M; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    Elem[i][j] *= Scal;
                }
            }
        }

        //умножение матрицы на вектор
        public static Vector operator *(Matrix T, Vector V)
        {
            if (T.N != V.N) throw new Exception("M * V: dim(Matrix) != dim(Vector)...");

            Vector RES = new Vector(T.M);

            for (int i = 0; i < T.M; i++)
            {
                for (int j = 0; j < T.N; j++)
                {
                    RES.Elem[i] += T.Elem[i][j] * V.Elem[j];
                }
            }
            return RES;
        }

        //умножение транспонированной матрицы на вектор
        public Vector Multiplication_Trans_Matrix_Vector(Vector V)
        {
            if (N != V.N) throw new Exception("Mt * V: dim(Matrix) != dim(Vector)...");

            Vector RES = new Vector(M);

            for (int i = 0; i < M; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    RES.Elem[i] += Elem[j][i] * V.Elem[j];
                }
            }
            return RES;
        }

        //умножение матрицы на матрицу
        public static Matrix operator *(Matrix T1, Matrix T2)
        {
            if (T1.N != T2.M) throw new Exception("M * M: dim(Matrix1) != dim(Matrix2)...");

            Matrix RES = new Matrix(T1.M, T2.N);

            for (int i = 0; i < T1.M; i++)
            {
                for (int j = 0; j < T2.N; j++)
                {
                    for (int k = 0; k < T1.N; k++)
                    {
                        RES.Elem[i][j] += T1.Elem[i][k] * T2.Elem[k][j];
                    }
                }
            }
            return RES;
        }

        //сложение матриц с выделением памяти под новую матрицу
        public static Matrix operator +(Matrix T1, Matrix T2)
        {
            if (T1.M != T2.M || T1.N != T2.N) throw new Exception("dim(Matrix1) != dim(Matrix2)...");

            Matrix RES = new Matrix(T1.M, T2.N);

            for (int i = 0; i < T1.M; i++)
            {
                for (int j = 0; j < T2.N; j++)
                {
                    RES.Elem[i][j] = T1.Elem[i][j] + T2.Elem[i][j];
                }
            }
            return RES;
        }

        //сложение матриц без выделения памяти под новую матрицу (добавление в ту же матрицу)
        public void Add(Matrix T2)
        {
            if (M != T2.M || N != T2.N) throw new Exception("dim(Matrix1) != dim(Matrix2)...");

            for (int i = 0; i < M; i++)
            {
                for (int j = 0; j < T2.N; j++)
                {
                    Elem[i][j] += T2.Elem[i][j];
                }
            }
        }

        //транспонирование матрицы
        public Matrix Transpose_Matrix()
        {
            var RES = new Matrix(N, M);

            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < M; j++)
                {
                    RES.Elem[i][j] = Elem[j][i];
                }
            }
            return RES;
        }

        //копирование матрицы A
        public void Copy(Matrix A)
        {
            if (N != A.N || M != A.M) throw new Exception("Copy: dim(matrix1) != dim(matrix2)...");
            for (int i = 0; i < M; i++)
                for (int j = 0; j < N; j++)
                    Elem[i][j] = A.Elem[i][j];
        }

        //вывод матрицы на консоль
        public void Console_Write_Matrix()
        {
            for (int i = 0; i < M; i++)
            {
                for (int j = 0; j < N; j++)
                    Console.Write(String.Format("{0, -22}", Elem[i][j].ToString("E5")));
                Console.WriteLine();
            }
        }

        public void Multiply_By_Noise(double Noise_Param)
        {
            if (N != M) throw new Exception("Multiply_By_Noise: Matrix not squared");
            for (int i = 0; i < N; i++)
                Elem[i][i] *= Noise_Param;
        }

        //вычисление числа обусловленности квадратной матрицы (параллельная версия)
        //делегат на вызов функции решателя СЛАУ
        delegate void Thread_Solver(int Number);
        //число обусловленности с нормой "бесконечность"
        public double Cond_InfinityNorm()
        {
            //проверка на "квадратность" матрицы
            if (M != N) throw new Exception("Cond(A): M != N ...");

            //решатель СЛАУ: A^t = QR и решаем системы A^t * A^(-t) = E
            var QR_Solver = new QR_Decomposition(Transpose_Matrix(), QR_Decomposition.QR_Algorithm.Householder);

            //проверка на невырожденность
            if (Math.Abs(QR_Solver.R.Elem[M - 1][M - 1]) < CONST.EPS)
                throw new Exception("Cond(A): detA = 0 ...");

            //число потоков
            int Number_Threads = Environment.ProcessorCount;

            //семафоры для потоков (по умолчанию false): сигнализируют, что i-ый поток завершился
            var Semaphores = new bool[Number_Threads];

            //максимальные нормы строк (вычисляются на каждом i-ом потоке) 
            var Norma_Row_A = new double[Number_Threads];
            var Norma_Row_A1 = new double[Number_Threads];

            //безымянная функция для решения СЛАУ -> столбцы обратной матрицы
            //Number - номер потока
            var Start_Solver = new Thread_Solver((Number) =>
            {
                //строка обратной матрицы
                var A1 = new Vector(M);
                double S1, S2;
                //первая и последняя обрабатываемые строки для потока
                int Begin = N / Number_Threads * Number;
                int End = Begin + N / Number_Threads;
                //в последний поток добавим остаток
                if (Number + 1 == Number_Threads) End += N % Number_Threads;

                //решаем системы A^t * A^(-t) = E
                for (int i = Begin; i < End; i++)
                {
                    A1.Elem[i] = 1.0;
                    A1 = QR_Solver.Start_Solver(A1);

                    S1 = 0; S2 = 0;
                    for (int j = 0; j < M; j++)
                    {
                        S1 += Math.Abs(Elem[i][j]);
                        S2 += Math.Abs(A1.Elem[j]);
                        A1.Elem[j] = 0.0;
                    }
                    if (Norma_Row_A[Number] < S1) Norma_Row_A[Number] = S1;
                    if (Norma_Row_A1[Number] < S2) Norma_Row_A1[Number] = S2;
                }
                //сигнал о завершении потока
                Semaphores[Number] = true;
            });


            //отцовский поток запускает дочерние
            for (int I = 0; I < Number_Threads - 1; I++)
            {
                int Number = Number_Threads - I - 1;
                ThreadPool.QueueUserWorkItem((Par) => Start_Solver(Number));
            }

            //отцовский поток забирает первую порцию строк
            Start_Solver(0);

            //ожидание отцовским потоком завершения работы дочерних
            while (Array.IndexOf<bool>(Semaphores, false) != -1) ;

            //поиск наибольшей нормы строки
            for (int i = 1; i < Number_Threads; i++)
            {
                if (Norma_Row_A[0] < Norma_Row_A[i]) Norma_Row_A[0] = Norma_Row_A[i];
                if (Norma_Row_A1[0] < Norma_Row_A1[i]) Norma_Row_A1[0] = Norma_Row_A1[i];
            }

            return Norma_Row_A[0] * Norma_Row_A1[0];
        }

        //чтение из бинарного файла

        

        //число обусловленности с нормой 1
        public double Cond_FirstNorm()
        {
            //проверка на "квадратность" матрицы
            if (M != N) throw new Exception("Cond(A): M != N ...");

            //решатель СЛАУ: A = QR и решаем системы A * A^(-1) = E
            var QR_Solver = new QR_Decomposition(this, QR_Decomposition.QR_Algorithm.Householder);

            ///////////////////////////////////////////////////////////////////////////////////////////
            //проверка на невырожденность
            if (Math.Abs(QR_Solver.R.Elem[M - 1][M - 1]) < CONST.EPS)
                throw new Exception("Cond(A): detA = 0 ...");

            //число потоков
            int Number_Threads = Environment.ProcessorCount;

            //семафоры для потоков (по умолчанию false): сигнализируют, что i-ый поток завершился
            var Semaphores = new bool[Number_Threads];

            //максимальные нормы столбцов (вычисляются на каждом i-ом потоке) 
            var Norma_Col_A = new double[Number_Threads];
            var Norma_Col_A1 = new double[Number_Threads];

            //безымянная функция для решения СЛАУ -> столбцы обратной матрицы
            //Number - номер потока
            var Start_Solver = new Thread_Solver((Number) =>
            {
                //строка обратной матрицы
                var A1 = new Vector(N);
                double S1, S2;
                //первая и последняя обрабатываемые строки для потока
                int Begin = N / Number_Threads * Number;
                int End = Begin + N / Number_Threads;
                //в последний поток добавим остаток
                if (Number + 1 == Number_Threads) End += N % Number_Threads;

                ///////////////////////////////////////////////////////////////////////////////////////////
                //решаем системы A * A^(-1) = E
                for (int i = Begin; i < End; i++)
                {
                    A1.Elem[i] = 1.0;
                    A1 = QR_Solver.Start_Solver(A1);

                    S1 = 0; S2 = 0;
                    for (int j = 0; j < N; j++)
                    {
                        S1 += Math.Abs(Elem[i][j]);
                        S2 += Math.Abs(A1.Elem[j]);
                        A1.Elem[j] = 0.0;
                    }

                    if (Norma_Col_A[Number] < S1) Norma_Col_A[Number] = S1;
                    if (Norma_Col_A1[Number] < S2) Norma_Col_A1[Number] = S2;
                }
                //сигнал о завершении потока
                Semaphores[Number] = true;
            });

            //отцовский поток запускает дочерние
            for (int I = 0; I < Number_Threads - 1; I++)
            {
                int Number = Number_Threads - I - 1;
                ThreadPool.QueueUserWorkItem((Par) => Start_Solver(Number));
            }

            //отцовский поток забирает первую порцию строк
            Start_Solver(0);

            //ожидание отцовским потоком завершения работы дочерних
            while (Array.IndexOf<bool>(Semaphores, false) != -1) ;

            double temp1, temp2;
            //поиск наибольшей нормы строки
            for (int i = 0; i < Number_Threads; i++)
            {
                if (Norma_Col_A[0] < Norma_Col_A[i]) Norma_Col_A[0] = Norma_Col_A[i];
                if (Norma_Col_A1[0] < Norma_Col_A1[i]) Norma_Col_A1[0] = Norma_Col_A1[i];
                temp1 = Norma_Col_A[i];
                temp2 = Norma_Col_A1[i];

            }

            return Norma_Col_A[0] * Norma_Col_A1[0];
        }



        public Matrix(String PATH)
        {
            //чтение размера системы
            using (var Reader = new BinaryReader(File.Open(PATH + "Size.bin", FileMode.Open)))
            {
                M = Reader.ReadInt32();
                N = M;
            }

                
            //чтение матрицы
            using (var Reader = new BinaryReader(File.Open(PATH + "Matrix.bin", FileMode.Open)))
            {
                //чтение каждого значения из файла
                try
                {

                    Elem = new double[M][];
                    for (int i = 0; i < M; i++) Elem[i] = new double[N];


                    for (int i = 0; i < M; i++)
                    {
                        for (int j = 0; j < N; j++)
                        {
                            Elem[i][j] = Reader.ReadDouble();
                        }
                    }
                }
                catch { throw new Exception("Block_Vector: data file is not correct..."); }

            }

        }


    }



}