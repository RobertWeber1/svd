#include "catch.h"
#include <svd/calculator.h>

using vector_t = std::vector<double>;
using matrix_t = std::vector<vector_t>;
using SVD_t    = svd_v1::SVD_t<matrix_t>;

namespace std
{

matrix_t operator*(matrix_t const& matrix1, matrix_t const& matrix2)
{
	matrix_t matrix3;
	matrix3.resize(matrix1.size());
	for (std::uint32_t row = 0; row < matrix1.size(); row++)
	{
		matrix3[row].resize(matrix1[row].size());
		for (std::uint32_t col = 0; col < matrix1[row].size(); col++)
		{
			matrix3[row][col] = 0.00;
			for (std::uint32_t k = 0; k < matrix1[row].size(); k++)
				matrix3[row][col] += matrix1[row][k] * matrix2[k][col];
		}
	}

	return matrix3;
}

std::ostream& operator<<(std::ostream& os, matrix_t const& input)
{
	for (std::size_t row = 0; row < input.size(); row++)
	{
		for (std::size_t col = 0; col < input[row].size(); col++)
		{
			os << std::setprecision(5) << std::fixed << input[row][col] << " ";
		}
		os << "\n";
	}
	return os << "\n";
}

template<class T>
bool are_equal(T f1, T f2, T epsilon = std::numeric_limits<T>::epsilon())
{
	return
		(std::fabs(f1 - f2) <=
		 (epsilon * std::fmax(std::fabs(f1), std::fabs(f2))));
}

bool operator==(matrix_t const& lhs, matrix_t const& rhs)
{
	for (std::size_t row = 0; row < lhs.size(); row++)
	{
		for (std::size_t col = 0; col < lhs[row].size(); col++)
		{
			if(not are_equal(lhs[row][col], rhs[row][col], 0.0000000001))
			{
				return false;
			}
		}
	}

	return true;
}

} //namespace std


matrix_t check(SVD_t const& input)
{
	return input.u * input.s * svd_v1::transpose(input.v);
}

matrix_t generate(std::size_t rows, std::size_t cols)
{
	matrix_t input;
	std::srand(0);
	input.resize(rows);
	for (std::size_t row = 0; row < input.size(); row++)
	{
		input[row].resize(cols);
		for (std::size_t col = 0; col < input[row].size(); col++)
		{
			input[row][col] = std::rand() % 20 - std::rand() % 20;
		}
	}

	return input;
}


TEST_CASE("calculate svd")
{
	auto const input = generate(5,5);
	std::cout << "\ninput = \n" << input;

	auto const ret = svd_v1::compute_svd<vector_t, matrix_t>(input);

	std::cout << "\nS = \n" << ret.s;
	std::cout << "\nU = \n" << ret.u;
	std::cout << "\nV = \n" << ret.v;

	std::cout << "\nA = \n" << check(ret);

	REQUIRE(input == check(ret));
}