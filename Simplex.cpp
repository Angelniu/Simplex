#include <cstdio>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>
using namespace std;
const int MAXN  = 100;
const double INF = 999999.0;
const double EPS = 1e-9;
double b[MAXN]; //constraint in each equation
double x[MAXN]; //last result
double a[MAXN][MAXN]; //coefficient of each non-basic variable
double c[MAXN]; //coefficient of target function
int n; //denote number of non-basic variables as n
int m; //denote number of basic variables as m
int z; //target function
bool isBasic[MAXN]; //mark the basic

double b1[MAXN]; //constraint in each equation
double a1[MAXN][MAXN]; //coefficient of each non-basic variable
double c1[MAXN]; //coefficient of target function

void InitializeSimplex();
void Simplex();
void GetPivot(int leave, int enter);

int main() {

    //input a linear programming as slack form
    InitializeSimplex();
    Simplex();
    cout << "the optimal solution is " << z << endl;
    return 0;
}

void InitializeSimplex(){
    cin >> n >> m;
    memset(a, 0.0, sizeof(a));
    memset(isBasic, false, sizeof(isBasic));
    for (int i = 1; i <= n; i++) {
        cin >> c[i];
    }
    for (int i = 1; i <= m; i++) {
        cin >> b[i+n];
        for (int j = 1; j <= n; j++) {
            cin >> a[i+n][j];
        }
        isBasic[i+n] = true;
    }
    z = 0;
}

void Simplex() {
    int z = 0;
    while (1) {
        //find s non-basic variable to enter.if all c[] are less than 0, break;
        bool flag = false;
        int enter = 0;
        for (int i = 1; i <= n+m; i++) {
            if (!isBasic[i] && c[i] > 0) {
                flag = true;
                enter = i;
                break;
            }
        }
        if (!flag) break;

        //find a basic variable to leave
        double minimalConstraint = INF;
        int leave = 0;
        for (int i = 1; i <= n+m; i++) {
            if (isBasic[i] && a[i][enter] > 0) {
                if (minimalConstraint > b[i] / a[i][enter]) {
                    minimalConstraint = b[i] / a[i][enter];
                    leave = i;
                }
            }
        }

        //unbounded
        if (fabs(minimalConstraint - INF) < EPS) {
            cout << "UNBOUNDED!" << endl;
            return;
        }
        else {
            // cout <<"***"<<endl;
            // cout <<leave<< " "<<enter<<endl;
            GetPivot(leave, enter);
        }

    }

    int ans = 1;
    for (int i = 1; i <= n+m; i++) {
         if (isBasic[i]) {
            x[ans] = b[i];
            cout <<"x" << ans <<": " << x[ans] << endl;
            ans ++;
         }
    }
}

void GetPivot(int leave, int enter) {
    //compute the coefficients of the equation for new basic variable xe
    b1[enter] = b[leave] / a[leave][enter];
    for (int i = 1; i <= n+m; i++) {
        if (!isBasic[i] && (i != enter)) {
            a1[enter][i] = a[leave][i] / a[leave][enter];
        }
    }
    a1[enter][leave] = 1.0 / a[leave][enter];

    //compute the coefficients of the remaining constraints
    for (int i = 1; i <= n+m; i++) {
        if (isBasic[i] && (i != leave)){
            b1[i] = b[i] - a[i][enter] * b1[enter];
            for (int j = 1; j <= n+m; j++) {
                if (!isBasic[j] && (j != enter)) {
                    a1[i][j] = a[i][j]  - a[i][enter] * a1[enter][j];
                }
            }
            a1[i][leave] = -1 * a[i][enter] * a1[enter][leave];
        }
    }

    //compute the objective function
    z = z + c[enter] * b1[enter];
    for (int i = 1; i <= n+m; i++) {
        if (!isBasic[i] && (i != enter)) {
            c1[i] = c[i] - c[enter] * a1[enter][i];
        }
    }
    c1[leave] = -1 * c[enter] * a1[enter][leave];

    //update new sets of basic and nonbasic variables
    isBasic[leave] = false;
    isBasic[enter] = true;

    for (int i = 1; i <= n+m; i++) {
        c[i] = c1[i];
        b[i] = b1[i];
       // cout <<"c="<<c[i]<<" b="<<b[i]<<endl;
        for (int j = 1; j <= n+m; j++) {
            a[i][j] = a1[i][j];
           // if (a[i][j]!=0)
           // cout <<" a[i][j] = " <<a[i][j]<<endl;
        }
    }
}
/*
3 3
2 -3 3
7 1 1 1
-7 -1 -1 1
4 1 -2 2
// 6 1 0   9

3 3
3 1 2
30 1 1 3
24 2 2 5
36 4 1 2
//8 4 18 28
*/
