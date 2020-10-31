using System;
using System.Diagnostics;
using System.Threading;

namespace Com_Methods
{
    class CONST
    {
        //точность решения
        public static double EPS = 1e-20;
    }

    class Tools
    {
        //замер времени
        public static string Measurement_Time(Thread thread)
        {
            Stopwatch stopWatch = new Stopwatch();
            stopWatch.Start();
            thread.Start();
            while (thread.IsAlive) ;
            stopWatch.Stop();
            TimeSpan ts = stopWatch.Elapsed;
            string elapsedTime = String.Format("{0:00}:{1:00}:{2:00}.{3:00}", ts.Hours, ts.Minutes, ts.Seconds, ts.Milliseconds / 10);
            return ("RunTime: " + elapsedTime);
        }

        //относительная погрешность
        public static double Relative_Error(Vector X, Vector x)
        {
            double s = 0.0;
            for (int i = 0; i < X.N; i++)
            {
                s += Math.Pow(X.Elem[i] - x.Elem[i], 2);
            }
            return Math.Sqrt(s) / x.Norma();
        }

        //относительная невязка
        public static double Relative_Discrepancy(Matrix A, Vector X, Vector F)
        {
            var x = A * X;
            for (int i = 0; i < X.N; i++) x.Elem[i] -= F.Elem[i];
            return x.Norma() / F.Norma();
        }
    }
    class Program
    {
        static void Main(string[] args)
        {
         
            try
            {
                //double NOISE_PARAM = 0.000001;
                //int SIZE_BLOCK = 1;
                //СЛАУ
                
                var A = new Matrix("Data\\Dense_Format\\System3\\");
                var F = new Vector("Data\\Dense_Format\\System3\\");

                /*Console.WriteLine("Noise_param = " + NOISE_PARAM);

                A.Multiply_By_Noise(NOISE_PARAM);

                Console.WriteLine("||A||1 = " + A.Cond_FirstNorm());
                Console.WriteLine("||A||inf = " + A.Cond_InfinityNorm());
                */

                //решатель СЛАУ
                var Solver = new Jacobi_Method(30000, 1e-12);


                Stopwatch stopWatch = new Stopwatch();
                stopWatch.Start();

                //решение СЛАУ
                var X = Solver.Start_Solver(A, F);

                X.Console_Write_Vector();

                stopWatch.Stop();
                TimeSpan ts = stopWatch.Elapsed;
                string elapsedTime = String.Format("{0:00}:{1:00}:{2:00}.{3:00}", ts.Hours, ts.Minutes, ts.Seconds, ts.Milliseconds / 10);

                X.Console_Write_Vector();
                Console.WriteLine("RunTime: " + elapsedTime);

                /* int k = 500;    ///////////////////////////////////////////////////time
                 Random rnd = new Random();
                 Matrix M = new Matrix(k, k);
                 Vector V =new Vector(k);
                 Vector Fau;
                 Vector X;
                 Vector Y = new Vector(k);

                 for (int i = 0; i < M.M; i++)
                 {
                     for (int j = 0; j < M.N; j++)
                         M.Elem[i][j] = rnd.NextDouble();
                     V.Elem[i] = rnd.NextDouble();
                 }

                 double fault;

                 Stopwatch stopwatch = new Stopwatch();
                 stopwatch.Start();
                 QR_Decomposition qr = new QR_Decomposition(M);
                 X = qr.Solve_Householder(V);
                 stopwatch.Stop();
                 //X.ConsoleWriteVector();

                 TimeSpan ts = stopwatch.Elapsed;
                 string elapsedTime = String.Format("{0:00}:{1:00}:{2:00}.{3:00}",
                 ts.Hours, ts.Minutes, ts.Seconds,
                 ts.Milliseconds / 10);
                 Console.WriteLine("RunTime Householder " + elapsedTime);

                 Fau = X - Y;
                 fault = Fau.Norma() / Y.Norma();
                 Console.WriteLine(fault);
                 Console.WriteLine();
                 Console.WriteLine();


                 Stopwatch stopwatch2 = new Stopwatch();
                 stopwatch2.Start();
                 X = Gauss_Method.Solve(M, V);
                 stopwatch2.Stop();
                 //X.ConsoleWriteVector();

                 ts = stopwatch2.Elapsed;
                 elapsedTime = String.Format("{0:00}:{1:00}:{2:00}.{3:00}",
                 ts.Hours, ts.Minutes, ts.Seconds,
                 ts.Milliseconds / 10);
                 Console.WriteLine("RunTime Gauss " + elapsedTime);

                 Fau = X - Y;
                 fault = Fau.Norma() / Y.Norma();
                 Console.WriteLine(fault);
                 Console.WriteLine();
                 Console.WriteLine();


                 Stopwatch stopwatch1 = new Stopwatch();
                 stopwatch1.Start();
                 QR_Decomposition qr1 = new QR_Decomposition(M);
                 X = qr1.Solve_Givens(V);
                 stopwatch1.Stop();
                 //X.ConsoleWriteVector();

                 ts = stopwatch1.Elapsed;
                 elapsedTime = String.Format("{0:00}:{1:00}:{2:00}.{3:00}",
                 ts.Hours, ts.Minutes, ts.Seconds,
                 ts.Milliseconds / 10);
                 Console.WriteLine("RunTime Givens " + elapsedTime);

                 Fau = X - Y;
                 fault = Fau.Norma() / Y.Norma();
                 Console.WriteLine(fault);
                 Console.WriteLine();
                 Console.WriteLine();


                 Stopwatch stopwatch3 = new Stopwatch();
                 stopwatch3.Start();
                 LU lu = new LU(M, V);
                 X = lu.X;
                 //LU_Decomposition lu = new LU_Decomposition(M);
                 //X = lu.Solve(V);
                 stopwatch3.Stop();
                 //X.ConsoleWriteVector();

                 ts = stopwatch3.Elapsed;
                 elapsedTime = String.Format("{0:00}:{1:00}:{2:00}.{3:00}",
                 ts.Hours, ts.Minutes, ts.Seconds,
                 ts.Milliseconds / 10);
                 Console.WriteLine("RunTime LU  " + elapsedTime);

                 Fau = X - Y;
                 fault = Fau.Norma() / Y.Norma();
                 Console.WriteLine(fault);
                 Console.WriteLine();
                 Console.WriteLine();


                 Stopwatch stopwatch4 = new Stopwatch();
                 stopwatch4.Start();
                 QR_Decomposition qr2 = new QR_Decomposition(M,1);
                 X = qr2.Solve_Gram_Schmidt(V);
                 stopwatch4.Stop();
                 //X.ConsoleWriteVector();

                 ts = stopwatch4.Elapsed;
                 elapsedTime = String.Format("{0:00}:{1:00}:{2:00}.{3:00}",
                 ts.Hours, ts.Minutes, ts.Seconds,
                 ts.Milliseconds / 10);
                 Console.WriteLine("RunTime GH  " + elapsedTime);


                 Fau = X - Y;
                 fault = Fau.Norma() / Y.Norma();
                 Console.WriteLine(fault);
                 Console.WriteLine();
                 Console.WriteLine();

                 /*
                 X = Gauss_Method.Solve(A, F);

                 LU_Decomposition lu = new LU_Decomposition(A);
                 X = lu.Solve(F);////c первой матрицей не работает

                 QR_Decomposition qr = new QR_Decomposition(A);
                 X = qr.Solve(F, 1);
                 Console.WriteLine();
                 Console.WriteLine();
                 QR_Decomposition qr1 = new QR_Decomposition(A);
                 Y = qr1.Solve(F);

                 A.OrthogonalizationMatrix();

                 A.ConsoleWriteMatrix();
                 F.ConsoleWriteVector();
                 X.ConsoleWriteVector();
                 Y.ConsoleWriteVector();
                 */


                Console.ReadKey();
            }
            catch(Exception e)
            {
                Console.Write(e.Message);
                Console.ReadKey();
            }


        }
    }
}
