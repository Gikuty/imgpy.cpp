#include<iostream>
#include <algorithm>
#include<bits/stdc++.h>
using namespace std;
/*int gcd(int a, int b)
{
    int res= min(a,b);
    while(res>0){
    if(a%res==0 && b%res==0){
        break;
    }
    res--;

    
    }
    return res;

}
int main()
{
    int a,b;
    cin>>a>>b;
    cout<<gcd(a,b)<<endl;
    return 0;
}*/
//euclidean method

int gcd(int a, int b)
{
    if(b==0)
    return a;
    else 
    return (b, a%b);
}
int main()
{
    int a,b;
    cin>>a>>b;
    cout<<gcd(a,b)<<endl;
    return 0;
}