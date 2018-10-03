#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>

#define EPSILON 0.0001
#define PI 3.14159265359

using namespace std;


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
    cout << "Matriz Preenchida: " << endl;
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
    cout << endl << "Quantas Iteracoes Pretende Fazer: ";
    cin >> iteracoes;

    cout << endl << "OBS: O Epsilon pode ser alterado no cabecalho do programa." << endl;

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

    /// Alocando a Matriz U
    double** U;

    U = new double*[tamanho];
    for(int i = 0; i < tamanho; i++)
        U[i] = new double[tamanho];

    /// Alocando a Matriz U_1
    double** U_1;

    U_1 = new double*[tamanho];
    for(int i = 0; i < tamanho; i++)
        U_1[i] = new double[tamanho];

    /// Alocando a Matriz U_2
    double** U_2;

    U_2 = new double*[tamanho];
    for(int i = 0; i < tamanho; i++)
        U_2[i] = new double[tamanho];

    for(int i = 0; i < tamanho; i++)
    {
        for(int j = 0; j < tamanho; j++)
        {
            if(i == j)
                U[i][j] = 1;
            else
                U[i][j] = 0;
        }
    }

    /// Opções para o Usuário
    bool matrizA = false;
    bool matrizU = false;

    cout << endl << "Deseja imprimir a Matriz A Durante cada Iteracao? (1/0): ";
    cin >> matrizA;

    cout << endl << "Deseja imprimir a Matriz U Durante cada Iteracao? (1/0): ";
    cin >> matrizU;
    cout << endl << "---------------------------------------------------------------" << endl << endl;

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
            cout << "------------------------------------------------" << endl;
            cout << "Maior Valor da Matriz Ultrapassa o EPSILON" << endl;
            cout << "EPSILON: " << EPSILON << endl;
            cout << "Maior Elemento: " << max_value << endl;
            cout << "------------------------------------------------" << endl;
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

        /// Preenchendo a Matriz U_1
        for(int i = 0; i < tamanho; i++)
        {
            for(int j = 0; j < tamanho; j++)
            {
                if(i == j)
                    U_1[i][j] = 1.0;
                else
                    U_1[i][j] = 0.0;
            }
        }

        U_1[p][p] = cosseno;
        U_1[q][q] = cosseno;
        U_1[p][q] = seno;
        U_1[q][p] = (-1) * seno;

        /// Zerando a Matriz Auxiliar
        for(int i = 0; i < tamanho; i++)
        {
            for(int j = 0; j < tamanho; j++)
            {
                U_2[i][j] = 0;
            }
        }

        /// Efetuando a Multiplicacao da Matriz U
        for(int i = 0; i < tamanho; i++)
        {
            for(int j = 0; j < tamanho; j++)
            {
                for(int k = 0; k < tamanho; k++)
                {
                    U_2[i][j] += U[i][k]*U_1[k][j];
                }
            }
        }

        /// Transferindo a Matriz Auxiliar para a Matriz U
        for(int i = 0; i < tamanho; i++)
        {
            for(int j = 0; j < tamanho; j++)
            {
                U[i][j] = U_2[i][j];
            }
        }


        /// Imprimindo as Matrizes
        if(matrizA)
        {
            cout << "Matriz A, Iteracao: " << iteracoes << endl;
            for(int i = 0; i < tamanho; i++)
            {
                for(int j = 0; j < tamanho; j++)
                {
                    cout << A[i][j] << "\t";
                }
                cout << endl;
            }
            cout << endl << endl;
        }

        if(matrizU)
        {
            cout << "Matriz U, Iteracao: " << iteracoes << endl;
            for(int i = 0; i < tamanho; i++)
            {
                for(int j = 0; j < tamanho; j++)
                {
                    cout << U[i][j] << "\t";
                }
                cout << endl;
            }
            cout << endl << endl;
        }

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

    /// Imprimindo a Matriz de Autovetores
    cout << endl << " - Autovetores - " << endl;
    for(int i = 0; i < tamanho; i++)
    {
        for(int j = 0; j < tamanho; j++)
        {
            cout << U[i][j] << "\t";
        }
        cout << endl;
    }
    cout << endl;

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

    for(int i = 0; i < tamanho; i++)
        delete [] U[i];
    delete [] U;

    for(int i = 0; i < tamanho; i++)
        delete [] U_1[i];
    delete [] U_1;

    for(int i = 0; i < tamanho; i++)
        delete [] U_2[i];
    delete [] U_2;

    return 0;
}
