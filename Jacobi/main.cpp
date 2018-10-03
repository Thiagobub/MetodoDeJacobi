#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>

#define EPSILON 0.0001
#define PI 3.14159265359

using namespace std;

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

int main()
{
/** Criando e Preenchendo a Matriz Principal **/
    double** A;
    int tamanho;

    /// Pegando o Tamanho da Matriz
    cout << "Digite o Tamanho da Matriz: ";
    cin >> tamanho;

    /// Alocando a Matriz
    A = new double*[tamanho];
    for(int i = 0; i < tamanho; i++)
        A[i] = new double[tamanho];

    /// Preenchimento da Matriz
    cout << "\nPreencha Aqui a Matriz Simetrica: " << endl;
    for(int i = 0; i < tamanho; i++)
    {
        for(int j = i; j < tamanho; j++)
        {
            cout << i << "." << j << ": ";
            cin >> A[i][j];
            A[j][i] = A[i][j];
        }
        cout << "\n";
    }

    /// Imprimindo a Matriz
    for(int i = 0; i < tamanho; i++)
    {
        for(int j = 0; j < tamanho; j++)
        {
            cout << A[i][j] << "\t";
        }
        cout << "\n";
    }

/** Método de Jacobi **/
    int iteracoes = 0;
    int p = 0, q = 0;
    double max_value = 0;
    double angulo = 0;
    double t = 0;
    double cosseno = 0;
    double seno = 0;

    /// Selecionando o Numero de Iterações
    cout << "Quantas Iteracoes Pretende Fazer: ";
    cin >> iteracoes;

    /// Alocando aMatriz A_1
    double** A_1;

    A_1 = new double*[tamanho];
    for(int i = 0; i < tamanho; i++)
    A_1[i] = new double[tamanho];

    /// Alocando a Matriz A_2
    double** A_2;

    A_2 = new double*[tamanho];
    for(int i = 0; i < tamanho; i++)
    A_2[i] = new double[tamanho];

    while(iteracoes > 0)
    {
        /// Calculando o Maior Elemento em Módulo Fora da Diagonal
        max_value = A[0][1];
        for(int i = 0; i < tamanho; i++)
        {
            for(int j = i; j < tamanho; j++)
            {
                if(i != j && fabs(A[i][j]) >= max_value)
                {
                    max_value = fabs(A[i][j]);
                    p = i;
                    q = j;
                }
            }
        }

        /// Se o Maior Valor da Matriz A For Menor que Epsilon
        /// Então Paramos o Programa
        if (max_value <= EPSILON)
        {
            break;
        }

        /// Calculos Para a Matriz A"
        angulo = (A[q][q] - A[p][p]) / (2 * A[p][q]);
        if(fabs(angulo) <= 0.001)
            t = 1.0;
        else
            t = 1.0 / (angulo + (((angulo > 0) - (angulo < 0)) * sqrt(pow(angulo, 2) + 1)));

        cosseno = 1 / sqrt(pow(t, 2) + 1);
        seno = t / sqrt(pow(t, 2) + 1);

        /// Calculo da Matriz A' = U(transposta) * A
        for(int i = 0; i < tamanho; i++)
        {
            for(int j = 0; j < tamanho; j++)
            {
                if(i == p)
                    A_1[p][j] = A[p][j]*cosseno - A[q][j]*seno;
                else if(i == q)
                    A_1[q][j] = A[p][j]*seno + A[q][j]*cosseno;
                else
                    A_1[i][j] = A[i][j];
            }
        }

        /// Calculo da Matriz A" = A' * U
        for(int i = 0; i < tamanho; i++)
        {
            for(int j = 0; j < tamanho; j++)
            {
                if(j == p)
                    A_2[i][p] = A_1[i][p]*cosseno - A_1[i][q]*seno;
                else if(j == q)
                    A_2[i][q] = A_1[i][p]*seno + A_1[i][q]*cosseno;
                else
                    A_2[i][j] = A_1[i][j];
            }
        }

        A_2[p][p] = A[p][p]*pow(cosseno, 2) - 2*A[p][q]*seno*cosseno + A[q][q]*pow(seno, 2);
        A_2[q][q] = A[p][p]*pow(seno, 2) + 2*A[p][q]*seno*cosseno + A[q][q]*pow(cosseno, 2);
        A_2[p][q] = (A[p][p] - A[q][q])*seno*cosseno + A[p][q]*(pow(cosseno, 2) - pow(seno, 2));
        A_2[q][p] = A_2[p][q];

        /// Transferindo a Matriz A" para a Matriz A
        for(int i = 0; i < tamanho; i++)
            for(int j = 0; j < tamanho; j++)
                A[i][j] = A_2[i][j];

        /* Imprimindo a Matriz
        for(int i = 0; i < tamanho; i++)
        {
            for(int j = 0; j < tamanho; j++)
            {
                cout << A[i][j] << "\t";
            }
            cout << endl;
        }
        cout << endl << endl << endl;*/

        /// Uma Iteração a Menos
        iteracoes--;
    }

    /// Imprimindo os Autovalores
    cout << endl << " - Autovalores - " << endl;
    for(int i = 0; i < tamanho; i++)
    {
        for(int j = 0; j < tamanho; j++)
        {
            if(i == j)
                cout << i << "." << j << ": " << A[i][j] << endl;
        }
    }

/** Desalocando as Matrizes Alocadas Dinamicamente **/
    for(int i = 0; i < tamanho; i++)
        delete [] A[i];
    delete [] A;

    for(int i = 0; i < tamanho; i++)
        delete [] A_1[i];
    delete [] A_1;

    for(int i = 0; i < tamanho; i++)
        delete [] A_2[i];
    delete [] A_2;

    return 0;
}
