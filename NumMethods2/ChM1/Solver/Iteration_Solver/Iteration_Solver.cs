using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Com_Methods
{
    public interface IIteration_Solver
    {
        //наибольшее число итераций
        int Max_Iter { set; get; }

        //точность решения
        double Eps { set; get; }

        //текущая итерация
        int Iter { set; get; }

    }
}
