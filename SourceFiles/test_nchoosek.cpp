#include <iostream>
using namespace std;
void combinate(int iPos, int iProc, int iTol, int iKey, int data[],int des[])
{
	if(iProc > iTol)
	{
		return;
	}
	if(iPos == iKey)
	{
		for(int i = 0;i < iKey; i++)
		{
			cout<<des[i]<<" ";
		}
		cout<<endl;
		return;
	}
	else
	{
		combinate(iPos,iProc+1,iTol,iKey,data,des);
		des[iPos] = data[iProc];
		combinate(iPos+1,iProc+1,iTol,iKey,data,des);
	}
}
int main()
 {
	int a[6] = {1,2,3,4,5,6}, b[5];
	combinate(0, 0, 5, 3 , a,b);
	cin>>n>>k;
	if(n > k > 0)
	{
		int data[n],temp[n];
		for(i = 0; i < n; i++ )
		{
			data[i] = i+1;
		}
		combinate(0, 0, n, k, data, temp);
	}
	return 0;
}

