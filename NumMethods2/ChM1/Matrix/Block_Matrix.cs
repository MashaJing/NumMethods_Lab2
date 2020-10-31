using System;
using System.IO;


namespace Com_Methods
{
    class Block_Matrix : IMatrix
    {
        //размер матрицы
        public int M { set; get; } //строки
        public int N { set; get; } //столбцы
        public int Size_Block { set; get; } //размер блока
        public Matrix[][] Block { set; get; } //элементы матрицы

        //конструктор по умолчанию
        public Block_Matrix() { }
        //конструктор по бинарному файлу
        public Block_Matrix(String PATH, int SIZE_BLOCK)
        {
            //чтение размера системы
            using (var Reader = new BinaryReader(File.Open(PATH + "Size.bin", FileMode.Open)))
            {
                M = Reader.ReadInt32();
                N = M;
            }

            if (M % SIZE_BLOCK != 0) throw new Exception("Block_Matrix: error in the block size");

            //размер матрицы
            N /= SIZE_BLOCK;
            M /= SIZE_BLOCK;
            Size_Block = SIZE_BLOCK;
            Block = new Matrix[M][];

            //чтение матрицы
            using (var Reader = new BinaryReader(File.Open(PATH + "Matrix.bin", FileMode.Open)))
            {
                //чтение каждого значения из файла
                try
                {
                    for (int i = 0; i < M; i++)
                    {
                        //выделили место под блочную матрицу
                        Block[i] = new Matrix[N];

                        //выделили место под каждый блок
                        for (int j = 0; j < N; j++) Block[i][j] = new Matrix(Size_Block, Size_Block);

                        //заполняем строки в блоках i-ой строки блочной матрицы
                        for (int ii = 0; ii < Size_Block; ii++)
                        {
                            for (int j = 0; j < N; j++)
                                for (int k = 0; k < Size_Block; k++)
                                {
                                    Block[i][j].Elem[ii][k] = Reader.ReadDouble();
                                }
                        }

                        //диагональный блок необходимо преобазовать в LU-разложение
                        var LU_Decomp = new LU_Decomposition(Block[i][i]);
                        Block[i][i] = LU_Decomp.LU;


                    }

                }
                catch { throw new Exception("Block_Vector: data file is not correct..."); }

            }

        }
        public void Multiply_By_Noise(double Noise_Param)
        {
            if (N != M) throw new Exception("Multiply_By_Noise: Matrix not squared");
            for (int i = 0; i < N; i++)
                for (int j = 0; j < Size_Block; j++) //в целях экономии времени не будем обращаться к методу multiply_by_noise каждого блока
                    Block[i][i].Elem[j][j]*=Noise_Param;
        }
    }
}
