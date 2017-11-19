#pragma once

#include <algorithm>
#include <vector>
#include <map>

using namespace std;

class FaulhaberModulo
{
	typedef map<int, long long> BernouilliCache;

public:
	FaulhaberModulo(int maxPower, int modulo);

	long long GetSumOfNSquared(long long n, int p);

private:
	void InitialisePascalsTriangleCache();
	void InitialiseBernoulliCache();

	int ModularInverse(int a) const;
	long long BernoulliNumberMod(int p);
	long long GetBernouilliPolynomial(int p, int power);
	long long SafePower(long long n, int p) const;


	BernouilliCache m_answerCache;
	vector<vector<long long>> m_pascalTriangles;
	const int m_maxPower;
	const int m_modulo;
};