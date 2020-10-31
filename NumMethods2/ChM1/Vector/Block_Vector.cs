using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Com_Methods
{
    class Block_Vector : IVector
    {
        //размер вектора
        public int N { set; get; }

        //размер блока
        public int Block_Size { set; get; }
        //элементы вектора
        public Vector[] Block { set; get; }

        //конструктор по умолчанию
        public Block_Vector() { }

        //конструктор по размеру
        public Block_Vector(int BLOCK_AMOUNT, int BLOCK_SIZE) 
        {
            N = BLOCK_AMOUNT;
            Block_Size = BLOCK_SIZE;
            Block = new Vector[BLOCK_AMOUNT];
            for (int i = 0; i < N; i++)
                Block[i] = new Vector(Block_Size);
        }

        public Block_Vector(String PATH, int BLOCK_SIZE)
        {
            int Matrix_Size;
            //чтение размера системы
            using (var Reader = new BinaryReader(File.Open(PATH + "Size.bin", FileMode.Open)))
            {
                Matrix_Size = Reader.ReadInt32();
            }

            if (Matrix_Size % BLOCK_SIZE != 0) throw new Exception("Block_Vector: invalid block size");

            //размер вектора
            N = Matrix_Size/ BLOCK_SIZE;

            Block_Size = BLOCK_SIZE;
            Block = new Vector[N];
            
            //чтение вектора правой части
            using (var Reader = new BinaryReader(File.Open(PATH + "F.bin", FileMode.Open)))
            {
                //чтение данных из файла
                try
                {
                    for (int i = 0; i < N; i++)
                    {
                        Block[i] = new Vector(Block_Size);
                        for (int j = 0; j < Block_Size; j++)
                            Block[i].Elem[j] = Reader.ReadDouble();
                    }

                }
                catch { throw new Exception("Block_Vector: data file is not correct..."); }
                
            }

        }
        //Вывод вектора на консоль
        public void Console_Write_Vector()
        {
            for (int i=0; i < N; i++)
                for (int j = 0; j < Block_Size; j++)
                {
                    Console.WriteLine("X[{0}] = {1, -20}", i*Block_Size + j + 1, Block[i].Elem[j]);
                } 
        }

    }
}
