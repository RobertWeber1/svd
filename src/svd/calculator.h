#pragma once
#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <limits>


namespace svd_v1
{

template<class Matrix>
Matrix transpose(Matrix const& input)
{
	Matrix result(input.size());

	for (std::uint32_t row = 0; row < input.size(); row++)
	{
		result[row].resize(input[row].size());
		for (std::uint32_t col = 0; col < input[row].size(); col++)
		{
			result[row][col] = input[col][row];
		}
	}

	return result;
}

template<class Matrix>
Matrix reduce_matrix(Matrix input, std::size_t new_size)
{
	Matrix r_matrix(new_size);

	std::size_t index_d = input.size() - new_size;
	std::uint32_t row = index_d, row_n = 0;

	while (row < input.size())
	{
		r_matrix[row_n].resize(new_size);

		std::uint32_t col = index_d;
		std::uint32_t col_n = 0;

		while (col < input.size())
		{
			r_matrix[row_n][col_n++] = input[row][col++];
		}

		row++; row_n++;
	}

	return r_matrix;
}

template<class Matrix>
Matrix inverse_diagonal(Matrix input)
{
	Matrix inv_matrix(input.size());

	for (std::uint32_t index = 0; index < input.size(); index++)
	{
		inv_matrix[index].resize(input[index].size());
		inv_matrix[index][index] = 1.0 / input[index][index];
	}

	return inv_matrix;
}

template<class Vector, class Matrix>
Vector get_eigenvector(Matrix input)
{
	Vector eigenvector;

	const double eps = 0.000001;
	bool eigenv_found = false;

	for (std::uint32_t s = 0; s < input.size() - 1 && !eigenv_found; s++)
	{
		std::uint32_t col = s; double alpha = input[s][s];
		while (col < input[s].size() && alpha != 0 && alpha != 1)
		{
			input[s][col++] /= alpha;
		}

		for (std::uint32_t col = s; col < input[s].size() && !alpha; col++)
		{
			std::swap(input[s][col], input[s + 1][col]);
		}

		for (std::uint32_t row = 0; row < input.size(); row++)
		{
			double gamma = input[row][s];
			for (std::uint32_t col = s; col < input[row].size() && row != s; col++)
			{
				input[row][col] = input[row][col] - input[s][col] * gamma;
			}
		}

		std::uint32_t row = 0;
		while (row < input.size() &&
			(s == input.size() - 1 || std::fabs(input[s + 1][s + 1]) < eps))
		{
			eigenvector.push_back(-input[row++][s + 1]);
		}

		if (eigenvector.size() == input.size())
		{
			eigenv_found = true;
			eigenvector[s + 1] = 1;
			for (std::uint32_t index = s + 1; index < eigenvector.size(); index++)
			{
				eigenvector[index] = (std::fabs(eigenvector[index]) >= eps) ? eigenvector[index] : 0;
			}
		}
	}

	return eigenvector;
}

template<class Matrix, class Vector>
Matrix get_hermitian(Vector const& eigenvector)
{
	Matrix h_matrix(eigenvector.size());

	for (std::uint32_t row = 0; row < eigenvector.size(); row++)
		h_matrix[row].resize(eigenvector.size());

	h_matrix[0][0] = 1 / eigenvector[0];
	for (std::uint32_t row = 1; row < eigenvector.size(); row++)
		h_matrix[row][0] = -eigenvector[row] / eigenvector[0];

	for (std::uint32_t row = 1; row < eigenvector.size(); row++)
		h_matrix[row][row] = 1;

	return h_matrix;
}

template<class Matrix, class Vector>
Matrix get_hermitian_inverse(Vector const& eigenvector)
{
	Matrix ih_matrix(eigenvector.size());

	for (std::uint32_t row = 0; row < eigenvector.size(); row++)
		ih_matrix[row].resize(eigenvector.size());

	ih_matrix[0][0] = eigenvector[0];
	for (std::uint32_t row = 1; row < eigenvector.size(); row++)
		ih_matrix[row][0] = -eigenvector[row];

	for (std::uint32_t row = 1; row < eigenvector.size(); row++)
		ih_matrix[row][row] = 1;

	return ih_matrix;
}

template<class Vector, class Matrix>
struct Eigen_t
{
	Vector values;
	Matrix vectors;
};

template<class Vector, class Matrix>
Eigen_t<Vector, Matrix> compute_evd(
	Matrix const& input,
	Matrix state,
	Eigen_t<Vector, Matrix> eigen,
	std::size_t eig_count)
{
	using Arg = double;
	const std::size_t m_size = input.size();
	Vector vec(m_size);

	std::fill_n(vec.begin(), m_size, 1);

	if (eigen.values.size() == 0 && eigen.vectors.size() == 0)
	{
		eigen.values.resize(m_size);
		eigen.vectors.resize(eigen.values.size());
	}

	Matrix m(m_size);

	for (std::uint32_t row = 0; row < m_size; row++)
	{
		m[row].resize(100);
	}

	Arg lambda_old = 0;

	std::uint32_t index = 0;
	bool is_eval = false;

	while (is_eval == false)
	{
		for (std::uint32_t row = 0; row < m_size && (index % 100) == 0; row++)
		{
			m[row].resize(m[row].size() + 100);
		}

		for (std::uint32_t row = 0; row < m_size; row++)
		{
			m[row][index] = 0;
			for (std::uint32_t col = 0; col < m_size; col++)
			{
				m[row][index] += input[row][col] * vec[col];
			}
		}

		for (std::uint32_t col = 0; col < m_size; col++)
		{
			vec[col] = m[col][index];
		}

		if (index > 0)
		{
			Arg lambda =
				(m[0][index - 1] != 0)
				? (m[0][index] / m[0][index - 1])
				: m[0][index];

			is_eval = (std::fabs(lambda - lambda_old) < 0.0000000001) ? true : false;

			lambda = (std::fabs(lambda) >= 10e-6) ? lambda : 0;
			eigen.values[eig_count] = lambda; lambda_old = lambda;
		}

		index++;
	}

	Matrix matrix_new;

	if (m_size > 1)
	{
		Matrix matrix_target(m_size);

		for (std::uint32_t row = 0; row < m_size; row++)
		{
			matrix_target[row].resize(m_size);
			for (std::uint32_t col = 0; col < m_size; col++)
				matrix_target[row][col] = (row == col) ? \
				(input[row][col] - eigen.values[eig_count]) : input[row][col];
		}

		Vector eigenvector = get_eigenvector<Vector>(matrix_target);

		Matrix hermitian_matrix = get_hermitian<Matrix>(eigenvector);

		Matrix ha_matrix_product = hermitian_matrix * input;

		Matrix inverse_hermitian_matrix = get_hermitian_inverse<Matrix>(eigenvector);

		Matrix iha_matrix_product = ha_matrix_product * inverse_hermitian_matrix;

		matrix_new = reduce_matrix(iha_matrix_product, m_size - 1);
	}

	if (m_size <= 1)
	{
		for (std::uint32_t index = 0; index < eigen.values.size(); index++)
		{
			Arg lambda = eigen.values[index];
			Matrix matrix_target(state.size());

			for (std::uint32_t row = 0; row < state.size(); row++)
			{
				matrix_target[row].resize(state.size());
				for (std::uint32_t col = 0; col < state.size(); col++)
					matrix_target[row][col] = (row == col) ? \
					(state[row][col] - lambda) : state[row][col];
			}

			eigen.vectors.resize(state.size());
			eigen.vectors[index] = get_eigenvector<Vector>(matrix_target);

			Arg eigsum_sq = 0;
			for (std::uint32_t v = 0; v < eigen.vectors[index].size(); v++)
				eigsum_sq += std::pow(eigen.vectors[index][v], 2.0);

			for (std::uint32_t v = 0; v < eigen.vectors[index].size(); v++)
				eigen.vectors[index][v] /= sqrt(eigsum_sq);

			eigen.values[index] = std::sqrt(eigen.values[index]);
		}

		return eigen;
	}

	return
		compute_evd<Vector, Matrix>(
			matrix_new,
			state,
			eigen,
			eig_count + 1);
}

template<class Vector, class Matrix>
Eigen_t<Vector, Matrix> compute_evd(Matrix const& input)
{
	return compute_evd<Vector>(input, input, Eigen_t<Vector, Matrix>{}, 0);
}

template<class Matrix>
struct SVD_t
{
	Matrix s;
	Matrix u;
	Matrix v;
};

template<class Vector, class Matrix>
SVD_t<Matrix> compute_svd(Matrix const& input)
{
	SVD_t<Matrix> result;

	auto const eigen = compute_evd<Vector>(transpose(input) * input);

	result.v = transpose(eigen.vectors);

	result.s.resize(input.size());
	for (std::uint32_t index = 0; index < eigen.values.size(); index++)
	{
		result.s[index].resize(eigen.values.size());
		result.s[index][index] = eigen.values[index];
	}

	result.u = input * result.v * inverse_diagonal(result.s);

	return result;
}

} //namespace svd_v1
