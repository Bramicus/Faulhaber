#include "FaulhaberModulo.h"

#include <assert.h>

FaulhaberModulo::FaulhaberModulo(int maxPower, int modulo) 
: m_maxPower(maxPower) 
, m_modulo(modulo)
{
	assert(m_modulo > 1);

	InitialisePascalsTriangleCache();
	InitialiseBernoulliCache();
}

void FaulhaberModulo::InitialisePascalsTriangleCache()
{
	//We resize to a triangular structure to store all of Pascal's triangle numbers
	const int offset = 2; //Allow +1 (since Faulhabers formula means we go max power+1, AND because we leave 0's for easy arithmatic.
	m_pascalTriangles.resize(m_maxPower + offset);
	
	for (int column = 0; column < m_maxPower + offset; ++column)
	{
		m_pascalTriangles[column].resize(column + offset);
		m_pascalTriangles[column][0] = 1;
	}

	for (int row = 0; row < m_maxPower + offset; ++row)
	{
		for (int column = 1; column <= row; ++column)
		{
			m_pascalTriangles[row][column] = (m_pascalTriangles[row - 1][column - 1] + m_pascalTriangles[row - 1][column]) % m_modulo;
		}
	}
}

void FaulhaberModulo::InitialiseBernoulliCache()
{
	m_answerCache[0] = 1;
	m_answerCache[1] = m_modulo - ModularInverse(2);
}

//Faulhaber's formula - https://en.wikipedia.org/wiki/Faulhaber%27s_formula
//n - the number of positive integers
//p - the pth power we take the numbers to
long long FaulhaberModulo::GetSumOfNSquared(long long n, int p)
{
	assert(p <= m_maxPower);

	long long ans = 0;
	for (int power = 0; power <= p; ++power)
	{
		const long long newAns = (GetBernouilliPolynomial(p + 1, power) * SafePower(n, p + 1 - power)) % m_modulo;
		ans = ((m_modulo + ans) + ((power & 0x1 == 1) ? -newAns : newAns)) % m_modulo;
	}
	ans = (ans * ModularInverse(p + 1)) % m_modulo;
	return ans;
}

//https://en.wikipedia.org/wiki/Modular_multiplicative_inverse
int FaulhaberModulo::ModularInverse(int a) const
{
	int b = m_modulo;
	const int initialB = m_modulo;
	int x0 = 0;
	int x1 = 1;

	while (a > 1)
	{
		const int q = a / b;
		const int oldB = b;
		b = a % b;
		a = oldB;
		const int oldx0 = x0;
		x0 = x1 - (q * x0);
		x1 = oldx0;
	}

	if (x1 < 0)
	{
		x1 += initialB;
	}
	return x1;
}

long long FaulhaberModulo::GetBernouilliPolynomial(int p, int power)
{
	return (m_pascalTriangles[p][power] * BernoulliNumberMod(power)) % m_modulo;
}

long long FaulhaberModulo::BernoulliNumberMod(int p)
{
	if (m_answerCache.find(p) == m_answerCache.end())
	{
		long long sum = (1 + m_pascalTriangles[p + 1][1] * m_answerCache[1]) % m_modulo;
		for (int power = 2; power < p; power += 2)
		{
			sum = (sum + GetBernouilliPolynomial(p + 1, power)) % m_modulo;
		}
		sum = ((sum * ModularInverse(p + 1)) % m_modulo) * -1;
		while (sum < 0)
		{
			sum += m_modulo;
		}
		assert(sum < m_modulo);
		m_answerCache[p] = sum;
	}
	return m_answerCache[p];
}

//Function to safely power a number to another without overflowing (by modulo-ing with a limit)
//n - number to power by
//p - number t power to
long long FaulhaberModulo::SafePower(long long n, int p) const
{
	assert(p > 0);
	long long ans = 1;
	long long ln = n % m_modulo;
	while (p != 0)
	{
		ans = ((p & 0x1) == 1) ? ((ans * ln) % m_modulo) : ans;
		ln = (ln * ln) % m_modulo;
		p >>= 1;
	}
	assert(ans < m_modulo);
	return ans;
}