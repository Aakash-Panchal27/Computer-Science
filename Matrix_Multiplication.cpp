/*
1.Blocked Matrix Multiplication
2.Simple Matrix Multiplciation: 6 types of loop permutation and its performance comparison
*/

#include<bits/stdc++.h>
using namespace std;

int **mult, **a, **b;
double **dmult, **da, **db;
int N;
int B = 8;

void block_mat_mul()
{
  for(int jj = 0; jj < N; jj += B)
	    for(int kk = 0; kk < N; kk += B)
	        for(int i = 0; i < N; i++)
	            for(int j = jj; j < jj + B; j++)
	                for(int k = kk; k < kk + B; k++)
	                    mult[i][j] += a[i][k]*b[k][j];
}

void mat_mul1()
{
  for(int i = 0; i < N; ++i)
        for(int j = 0; j < N; ++j)
            for(int k = 0; k < N; ++k)
            {
                mult[i][j] += a[i][k] * b[k][j];
            }
}
void mat_mul2()
{
  for(int i = 0; i < N; ++i)
        for(int j = 0; j < N; ++j)
            for(int k = 0; k < N; ++k)
            {
                mult[j][i] += a[j][k] * b[k][i];
            }
}
void mat_mul3()
{
  for(int i = 0; i < N; ++i)
        for(int j = 0; j < N; ++j)
            for(int k = 0; k < N; ++k)
            {
                mult[k][j] += a[k][i] * b[i][j];
            }
}
void mat_mul4()
{
  for(int i = 0; i < N; ++i)
        for(int j = 0; j < N; ++j)
            for(int k = 0; k < N; ++k)
            {
                mult[k][i] += a[k][j] * b[j][i];
            }
}
void mat_mul5()
{
  for(int i = 0; i < N; ++i)
        for(int j = 0; j < N; ++j)
            for(int k = 0; k < N; ++k)
            {
                mult[i][k] += a[i][j] * b[j][k];
            }
}
void mat_mul6()
{
  for(int i = 0; i < N; ++i)
        for(int j = 0; j < N; ++j)
            for(int k = 0; k < N; ++k)
            {
                mult[j][k] += a[j][i] * b[i][k];
            }
}

void double_mat_mul()
{
  for(int i = 0; i < N; ++i)
        for(int j = 0; j < N; ++j)
            for(int k = 0; k < N; ++k)
            {
                dmult[i][j] += da[i][k] * db[k][j];
            }
}

int main() 
{

//------------------------------------------Blocked Matrix comparison
  for(int i = 1<<8; i < 1<<12; i<<=1){
    cout << "size is " << i << endl;
    N = i;
    mult = new int*[i];
    a = new int*[i];
    b = new int*[i];
    for(int j = 0; j < i; j++) {
      mult[j] = new int[i];
      a[j] = new int[i];
      b[j] = new int[i];
    }
    
    for(int sz = 0; sz < i; sz++){
        for(int col = 0; col < i; col++){
            mult[sz][col] = 0;
            a[sz][col] = 10;
            b[sz][col] = 10;
        }
    }
    
    cout << "Simple Matrix Multiplication" << endl;
    auto start = chrono::high_resolution_clock::now();
    mat_mul1();
    auto end = chrono::high_resolution_clock::now();
    double total_time =  
      chrono::duration_cast<chrono::nanoseconds>(end - start).count(); 
    total_time *= 1e-9;
    cout << fixed << setprecision(20) << 2*pow(i,3)/total_time << endl;
    
    for(int sz = 0; sz < i; sz++)
        for(int col = 0; col < i; col++)
            mult[sz][col] = 0;
            
    cout << "Blocked Matrix Multiplication" << endl;
    start = chrono::high_resolution_clock::now();
    block_mat_mul();
    end = chrono::high_resolution_clock::now();
    total_time =  
      chrono::duration_cast<chrono::nanoseconds>(end - start).count();
    total_time *= 1e-9;
    cout << fixed << setprecision(20) << 2*pow(i,3)/total_time << endl;
  }

//------------------------------------------1
    cout << endl << endl;
  for(int i = 1<<9; i < 1<<12;i<<=2){
    N = i;
    mult = new int*[i];
    a = new int*[i];
    b = new int*[i];
    for(int j = 0; j < i; j++) {
      mult[j] = new int[i];
      a[j] = new int[i];
      b[j] = new int[i];
    }
    
    for(int sz = 0; sz < i; sz++){
        for(int col = 0; col < i; col++){
            mult[sz][col] = 0;
            a[sz][col] = 10;
            b[sz][col] = 10;
        }
    }
    
    
    auto start = chrono::high_resolution_clock::now();
    mat_mul1();
    auto end = chrono::high_resolution_clock::now();
    double total_time =  
      chrono::duration_cast<chrono::nanoseconds>(end - start).count(); 
    total_time *= 1e-9;
    cout << fixed << setprecision(20) << "i j k: " << total_time << endl;
    
    for(int sz = 0; sz < i; sz++){
        for(int col = 0; col < i; col++){
            mult[sz][col] = 0;
        }
    }
    
    start = chrono::high_resolution_clock::now();
    mat_mul2();
    end = chrono::high_resolution_clock::now();
    total_time =  
      chrono::duration_cast<chrono::nanoseconds>(end - start).count(); 
    total_time *= 1e-9;
    cout << fixed << setprecision(20) << "j i k: " << total_time << endl;
    
    for(int sz = 0; sz < i; sz++){
        for(int col = 0; col < i; col++){
            mult[sz][col] = 0;
        }
    }
    
    start = chrono::high_resolution_clock::now(); ;
    mat_mul3();
    end = chrono::high_resolution_clock::now();
    total_time =  
      chrono::duration_cast<chrono::nanoseconds>(end - start).count(); 
    total_time *= 1e-9;
    cout << fixed << setprecision(20) <<  "k j i: " << total_time << endl;
    
    for(int sz = 0; sz < i; sz++){
        for(int col = 0; col < i; col++){
            mult[sz][col] = 0;
        }
    }
    
    start = chrono::high_resolution_clock::now(); 
    mat_mul4();
    end = chrono::high_resolution_clock::now(); 
    total_time =  
      chrono::duration_cast<chrono::nanoseconds>(end - start).count(); 
    total_time *= 1e-9;
    cout << fixed << setprecision(20) <<  "k i j: " << total_time << endl;
    
    for(int sz = 0; sz < i; sz++){
        for(int col = 0; col < i; col++){
            mult[sz][col] = 0;
        }
    }
    
    start = chrono::high_resolution_clock::now(); 
    mat_mul5();
    end = chrono::high_resolution_clock::now();
    total_time =  
      chrono::duration_cast<chrono::nanoseconds>(end - start).count(); 
    total_time *= 1e-9;
    cout << fixed << setprecision(20) <<  "i k j: " << total_time << endl;

    for(int sz = 0; sz < i; sz++){
        for(int col = 0; col < i; col++){
            mult[sz][col] = 0;
        }
    }

    start = chrono::high_resolution_clock::now(); 
    mat_mul6();
    end = chrono::high_resolution_clock::now();
    total_time =  
      chrono::duration_cast<chrono::nanoseconds>(end - start).count(); 
    total_time *= 1e-9;
    cout << fixed << setprecision(20) <<  "j k i: " << total_time << endl;
    
    cout << endl;
    
  }
  
  //------------------------------------------2
  cout << endl << endl;
  for(int i = 1<<8; i < 1<<12; i<<=1){
    N = i;
    mult = new int*[i];
    a = new int*[i];
    b = new int*[i];
    for(int j = 0; j < i; j++) {
      mult[j] = new int[i];
      a[j] = new int[i];
      b[j] = new int[i];
    }
    
    for(int sz = 0; sz < i; sz++){
        for(int col = 0; col < i; col++){
            mult[sz][col] = 0;
            a[sz][col] = 10;
            b[sz][col] = 10;
        }
    }
    
    auto start = chrono::high_resolution_clock::now();
    mat_mul1();
    auto end = chrono::high_resolution_clock::now();
    double total_time =  
      chrono::duration_cast<chrono::nanoseconds>(end - start).count(); 
    total_time *= 1e-9;
    cout << i << " " << total_time << " " << 2*pow(N,3)/total_time << endl;
  }
  
  //------------------------------------------3
  cout << endl << endl;
  for(int i = 1<<10; i < 1<<12; i<<=1){
    cout << "size is " << i << endl;
    cout << "Single Precision" << endl;
    N = i;
    mult = new int*[i];
    a = new int*[i];
    b = new int*[i];
    for(int j = 0; j < i; j++) {
      mult[j] = new int[i];
      a[j] = new int[i];
      b[j] = new int[i];
    }
    
    for(int sz = 0; sz < i; sz++){
        for(int col = 0; col < i; col++){
            mult[sz][col] = 0;
            a[sz][col] = 10;
            b[sz][col] = 10;
        }
    }
    
    
    auto start = chrono::high_resolution_clock::now();
    mat_mul1();
    auto end = chrono::high_resolution_clock::now();
    double total_time =  
      chrono::duration_cast<chrono::nanoseconds>(end - start).count(); 
    total_time *= 1e-9;
    cout << fixed << setprecision(20) << total_time << " " << 2*pow(N,3)/total_time << endl;
    
    N = i;
    dmult = new double*[i];
    da = new double*[i];
    db = new double*[i];
    for(int j = 0; j < i; j++) {
      dmult[j] = new double[i];
      da[j] = new double[i];
      db[j] = new double[i];
    }
    
    for(int sz = 0; sz < i; sz++){
        for(int col = 0; col < i; col++){
            dmult[sz][col] = (double)0;
            da[sz][col] = (double)10;
            db[sz][col] = (double)10;
        }
    }
    
    start = chrono::high_resolution_clock::now();
    double_mat_mul();
    end = chrono::high_resolution_clock::now();
    total_time =  
      chrono::duration_cast<chrono::nanoseconds>(end - start).count(); 
    total_time *= 1e-9;
    cout << "Double Precision" << endl;
    cout << fixed << setprecision(20) << total_time << " " << 2*pow(N,3)/total_time << endl;
  }
  return 0;
}
