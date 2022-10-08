#include <iostream>
#include<bits/stdc++.h>
using namespace std;
using namespace std;
int main() {
   int a[3];
  int x,y;
   a[0]=7;
   a[1]=5;
   int a, b;
    cout << "Enter first number  : ";
    cin >> a;
    cout << "Enter second number : ";
    cin >> b;
     
    int result = a * b;
    cout << "Multiplication : " << result << endl;
   a[2]=a[0] + a[1];
   cout<<"Sum of "<<a[0]<<" and "<<a[1]<<" is "<<a[2];
   return 0;
}
