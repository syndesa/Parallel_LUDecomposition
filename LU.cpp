#include <iostream>
#include <stdlib.h>
#include <thread>
#include <vector>
#include <math.h>
#include "BS_thread_pool.hpp"


using namespace std;

const auto CPU_COUNT= std::thread::hardware_concurrency();
BS::thread_pool Pool(CPU_COUNT);  // https://github.com/bshoshany/thread-pool#getting-started

class LUDecomposition{
    private:
        vector<vector<float>> Matrix;
        vector<vector<float>> LMatrix; 
    public:
        LUDecomposition(){};
    
    // Define Matrix size with arguments size, and integer range from 1 - max
    void generateMatrix(int size, int max)
    {
        for(int i = 0 ; i < size; i++)
        {
            vector<float> row;
            vector<float> lrow;
            for(int j = 0; j < size; j++)
            {
                row.push_back(rand() % max + 1);
                j == i ? lrow.push_back(1) : lrow.push_back(0);      
            }
            Matrix.push_back(row);
            LMatrix.push_back(lrow);
        }
    };

    // Sequential Implementation of LU Decomposition
    void gaussElimination_s(int col, int row)
    {
        float factor = 0;
        for (int k = 0; k < row; k++){
            for (int i = k + 1; i < row; i++)
            {
                    factor = Matrix[i][k] / Matrix[k][k];
                    LMatrix[i][k] = factor;
                    for(int j = 0; j < col; j++)
                    {
                        Matrix[i][j] =  Matrix[i][j] -  factor*Matrix[k][j] ;
                    }
            }
        }
    };


    void solve_lower_triangle(int pos, int step, int k)
    {
        float factor = 0;
        for (int i = pos; i < pos + step; i++)
            { 
                factor = Matrix[i][k] / Matrix[k][k];
                LMatrix[i][k] = factor;
                for(int j = 0; j < Matrix.size(); j++)
                    {
                        Matrix[i][j] = Matrix[i][j] - factor * Matrix[k][j];
                    }
            }
    }

    // Parallel Implementation of LU Decomposition
    void gaussElimination_p(int col, int row)
    {
        for (int k = 0; k < row; k++){
           parallel_for(k, row);
        }
    };

    void parallel_for(int start, int size)
    {
        std::vector<thread> threads;
        int step = CPU_COUNT;
        for (int i = start + 1; i < size ; i += step)
        {
            if (i + step < size)
            {
                // Pool.submit(&LUDecomposition::solve_lower_triangle, this, i, step, start);
                Pool.push_task(&LUDecomposition::solve_lower_triangle, this, i, step, start);

            }
            else {
                // Pool.submit(&LUDecomposition::solve_lower_triangle, this, i, size - i, start);
                Pool.push_task(&LUDecomposition::solve_lower_triangle, this, i, size - i, start);
            }
        }
    };

    //Print Matrices L and U 
    void printMatrix()
    {       
        cout << "Matrix: " << endl; 

        for (int row=0; row<Matrix.size(); row++)
        {
            for (int col=0; col<Matrix.size(); col++)
            {
                cout << Matrix[row][col] << " ";
            }
            cout << endl << "********************************************" << endl;
        }
        cout << endl << "LMatrix: " << endl; 

        for (int row=0; row<LMatrix.size(); row++)
        {
            for (int col=0; col<LMatrix.size(); col++)
            {
                cout << LMatrix[row][col] << " ";
            }
            cout << endl << "********************************************" << endl;
        }
    }

 

    //Proof of LU Decomposition -- Multiply Matrices L and U to Obtain original Matrix A
    void printMatrixMul() {
        cout << "Matrix A = L*U" << endl;
        for (int row = 0; row < Matrix.size(); row++)
        {

            for (int col = 0; col < Matrix.size(); col++) {
                float total = 0;
                for (int i = 0; i < Matrix.size(); i++) {
                    total += LMatrix[row][i] * Matrix[i][col];
                }
                cout << round(total) << " ";
            }
            cout << endl;
        }
    };


};








int main()
{


    LUDecomposition sLU;	
    LUDecomposition pLU;
    sLU.generateMatrix(1000, 100); // 100 x 100 Matrix with range 1-100
    pLU.generateMatrix(1000, 100); // 100 x 100 Matrix with range 1-100


    auto seqStart = chrono::high_resolution_clock::now();
    sLU.gaussElimination_s(1000, 1000);
    auto seqStop = chrono::high_resolution_clock::now();
    auto seqDuration = chrono::duration_cast<chrono::milliseconds>(seqStop - seqStart);

    auto parStart = chrono::high_resolution_clock::now();
    pLU.gaussElimination_p(1000, 1000);
    auto parStop = chrono::high_resolution_clock::now();
    auto parDuration = chrono::duration_cast<chrono::milliseconds>(parStop - parStart);


    //LU.printMatrix();
    //LU.printMatrixMul();


	cout << "Sequential Implementation: " << seqDuration.count() << " ms" << endl;
	cout << "Parallel Implementation: " << parDuration.count() << " ms";
    
    return 0;
};
