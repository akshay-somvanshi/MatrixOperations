#include "Matrix_h.h"
#include <cassert>
#include <vector>
#include <algorithm>
#include <iterator>
#include <omp.h>
#include <numeric>

using namespace std;

typedef int int32;
typedef unsigned int uint32;
typedef double double64;

namespace adv_prog_cw
{
	// default constructor
	// ---------------------------------------
	template <typename fT>
	Matrix_02508971<fT>::Matrix_02508971()
		: rows(0), cols(0)
	{
	}

	// operator=( DM )
	// ---------------------------------------
	template <typename fT>
	Matrix_02508971<fT> &Matrix_02508971<fT>::operator=(const Matrix_02508971<fT> &mat)
	{
		if (&mat != this)
		{
			rows = mat.rows;
			cols = mat.cols;

			data.resize(rows, vector<fT>(cols));

			for (size_t i = 0U; i < rows; i++)
				for (size_t j = 0U; j < cols; j++)
					data[i][j] = mat.data[i][j];
		}

		return *this;
	}

	template <typename fT>
	Matrix_02508971<fT> &Matrix_02508971<fT>::operator=(fT val)
	{
		for (size_t i = 0U; i < rows; i++)
			for (size_t j = 0U; j < cols; j++)
				data[i][j] = val;

		return *this;
	}

	// copy constructor
	// ---------------------------------------
	template <typename fT>
	Matrix_02508971<fT>::Matrix_02508971(const Matrix_02508971<fT> &mat)
	{
		*this = mat;
	}

	// constructor (i,j, value)
	// ---------------------------------------
	template <typename fT>
	Matrix_02508971<fT>::Matrix_02508971(size_t m, size_t n, fT val)
		: rows(m), cols(n)
	{
		data.resize(rows, vector<fT>(cols));

		for (size_t i = 0U; i < rows; i++)
			for (size_t j = 0U; j < cols; j++)
				data[i][j] = val;
	}

	// destructor
	// ---------------------------------------
	template <typename fT>
	Matrix_02508971<fT>::~Matrix_02508971()
	{
	}

	// operator+=
	// ---------------------------------------
	template <typename fT>
	Matrix_02508971<fT> &Matrix_02508971<fT>::operator+=(const Matrix_02508971<fT> &mat)
	{
		CheckSizes(mat, "Matrix_02508971<fT>::operator+=");
		for (size_t i = 0U; i < rows; i++)
			for (size_t j = 0U; j < cols; j++)
				data[i][j] += mat.data[i][j];

		return *this;
	}

	// operator-=
	// ---------------------------------------
	template <typename fT>
	Matrix_02508971<fT> &Matrix_02508971<fT>::operator-=(const Matrix_02508971<fT> &mat)
	{
		CheckSizes(mat, "Matrix_02508971<fT>::operator-=");
		for (size_t i = 0U; i < rows; i++)
			for (size_t j = 0U; j < cols; j++)
				data[i][j] -= mat.data[i][j];

		return *this;
	}

	// operator +
	// ---------------------------------------
	template <typename fT>
	Matrix_02508971<fT> operator+(const Matrix_02508971<fT> &a, const Matrix_02508971<fT> &b)
	{
		Matrix_02508971<fT> temp(a);
		temp += b;
		return temp;
	}

	// operator -
	// ---------------------------------------
	template <typename fT>
	Matrix_02508971<fT> operator-(const Matrix_02508971<fT> &a, const Matrix_02508971<fT> &b)
	{
		Matrix_02508971<fT> temp(a);
		temp -= b;
		return temp;
	}

	// operator *
	// ---------------------------------------
	template <typename fT>
	Matrix_02508971<fT> operator*(const Matrix_02508971<fT> &a, const Matrix_02508971<fT> &b)
	{
#ifndef NDEBUG
		if (a.Cols() != b.Rows())
		{
			cout << "\nMatrix<" << typeid(fT).name() << "> operator*: Matrices cannot be multiplied ";
			cout << "because of incompatible sizes (A * B, see matrices below): " << endl;
			a.Out(3L);
			b.Out(3L);
			throw length_error("Matrix_02508971<double,mn_max>::operator*");
		}
#endif
		Matrix_02508971<fT> temp(a.Rows(), b.Cols());

		for (size_t i = 0U; i < a.Rows(); i++)
			for (size_t j = 0U; j < b.Cols(); j++)
			{
				temp(i, j) = static_cast<fT>(0.0);
				for (size_t k = 0U; k < b.Rows(); k++)
					temp(i, j) += a(i, k) * b(k, j);
			}

		return temp;
	}

	// -------------------------------------------------------------------------------------------
	// OPERATOR FUNCTIONS
	// -------------------------------------------------------------------------------------------

	// vector^T = Matrix * vector^T
	// ---------------------------------------
	template <typename fT>
	vector<fT> operator*(const Matrix_02508971<fT> &mat, const vector<fT> &vec)
	{
		assert(mat.Cols() == vec.size());

		vector<fT> temp(mat.Rows(), static_cast<fT>(0.0));

		for (size_t i = 0; i < mat.Rows(); i++)
			for (size_t j = 0; j < mat.Cols(); j++)
				temp[i] += mat(i, j) * vec[j];

		return temp;
	} // end operator*

	// Matrix = vector^T * Matrix
	// ---------------------------------------
	template <typename fT>
	Matrix_02508971<fT> operator*(const vector<fT> &vec, const Matrix_02508971<fT> &mat)
	{
		if (vec.size() != mat.Rows())
		{
			cerr << "\noperator*: vector cannot be multiplied with matrix";
			cerr << "because of incompatible sizes (v * M): " << endl;
			throw length_error("Matrix<double,mn_max>::operator*");
		}
		Matrix_02508971<fT> temp(vec.size(), mat.Cols());

		for (size_t i = 0U; i < vec.size(); i++)
			for (size_t j = 0U; j < mat.Cols(); j++)
			{
				temp(i, j) = static_cast<fT>(0.0);
				for (size_t k = 0U; k < mat.Rows(); k++)
					temp(i, j) += vec[k] * mat(k, j);
			}

		return temp;
	} // end operator

	// Resize()
	// ---------------------------------------------
	template <typename fT>
	void Matrix_02508971<fT>::Resize(size_t m, size_t n)
	{
		// resize matrix but keep storage as is if the new matrix
		// is smaller than the old matrix
		if (m <= rows && n <= cols)
		{
			rows = m;
			cols = n;
			return;
		}

		// increase matrix size
		rows = m;
		cols = n;
		data.resize(m);
		for (size_t i = 0; i < m; i++)
			data[i].resize(n);
	}

	// Identity()
	// ---------------------------------------------
	template <typename fT>
	void Matrix_02508971<fT>::Identity()
	{
		for (size_t i = 0U; i < rows; i++)
			for (size_t j = 0U; j < cols; j++)
				if (i == j)
					data[i][j] = static_cast<fT>(1.0);
				else
					data[i][j] = static_cast<fT>(0.0);
	}

	// Zero()
	// ---------------------------------------------
	template <typename fT>
	void Matrix_02508971<fT>::Zero()
	{
		for (size_t i = 0U; i < rows; i++)
			for (size_t j = 0U; j < cols; j++)
				data[i][j] = static_cast<fT>(0.0);
	}

	// Transposed()
	// ---------------------------------------------
	template <typename fT>
	void Matrix_02508971<fT>::Transposed(Matrix_02508971 &M) const
	{
		M.Resize(cols, rows);

		for (size_t i = 0U; i < rows; i++)
			for (size_t j = 0U; j < cols; j++)
				M.data[j][i] = data[i][j];
	}

	// Out( digits )
	// ---------------------------------------
	template <typename fT>
	void Matrix_02508971<fT>::Out(long digits) const
	{
		std::streamsize prec;
		cout << "\nMatrix<" << typeid(fT).name() << ">::Out(): m=" << rows << ", n=" << cols << endl;
		if (digits != 0U)
		{
			cout.setf(ios::scientific);
			prec = cout.precision(digits);
		}
		size_t row_break, split_after(10U);

		for (size_t i = 0; i < rows; i++)
		{
			row_break = 1;
			for (size_t j = 0; j < cols; j++, row_break++)
			{
				if (data[i][j] >= 0.)
					cout << " ";
				cout << data[i][j] << " ";
				if (row_break == split_after)
				{
					cout << endl;
					row_break = 0U;
				}
			}
			cout << endl;
		}

		if (digits != 0U)
		{
			cout.unsetf(ios::scientific);
			cout.precision(prec);
		}

		cout << endl;
	} // end Out()

	// Matrix = Matrix * scalar
	// ---------------------------------------

	/* 
	Method to multiply the matrix by a scalar using the *= operator
	Parameters: 
		fT : scalar to multiply matrix by

	Returns: 
		Matrix_02508971<fT> : Matrix multiplied by the scalar
	*/

	template <typename fT>
	Matrix_02508971<fT> &Matrix_02508971<fT>::operator*=(fT scalar)
	{
		// Calculate this outside the loop
		// We use loop unrolling to optmise the code. 
		size_t unroll_limit = (cols / 2) * 2;

#pragma omp parallel for
		for (size_t i = 0; i < rows; i++)
		{
			for (size_t j = 0; j + 1 < cols; j += 2)
			{
				data[i][j] *= scalar;
				data[i][j + 1] *= scalar;
			}
			for (size_t j = unroll_limit; j < cols; j++)
			{
				data[i][j] *= scalar;
			}
		}

		// Recalculate L,U if matrix changes 
		LUdecom = false; 

		return *this;
	} // end operator*=

	// Matrix = Matrix / scalar
	// ---------------------------------------

	/* 
	Method to divide the matrix by a scalar using the /= operator
	Parameters: 
		fT scalar: scalar to divide matrix by

	Returns: 
		Matrix_02508971<fT> : Matrix divided by the scalar
	*/

	template <typename fT>
	Matrix_02508971<fT> &Matrix_02508971<fT>::operator/=(fT scalar)
	{	
		// Error handling in case of division by 0.
		if (scalar == 0)
		{
			throw logic_error("Cannot divide by 0!");
		}

		// Calculate these outside the loop
		double rec_scalar = 1 / scalar;
		double unroll_limit = (cols / 2) * 2;

#pragma omp parallel for
		for (size_t i = 0; i < rows; i++)
		{
			for (size_t j = 0; j + 1 < cols; j += 2)
			{
				// Perform multiplication as it is a cheaper operation
				data[i][j] *= rec_scalar;
				data[i][j + 1] *= rec_scalar;
			}
			for (double j = unroll_limit; j < cols; j++)
			{
				data[i][j] *= rec_scalar;
			}
		}

		// Recalculate L,U if matrix changes 
		LUdecom = false; 

		return *this;
	} // end operator/=

	// LU Decomposition
	// ----------------------------------------

	/* 
	Performs the LU Decomposition of the matrix for further calculations.
	Changes the Object vector parameters L and U by storing the computed values in those
	vectors. 

	Used for small dense matrices
	*/

	template <typename fT>
	void Matrix_02508971<fT>::LUdecomposition() const
	{
		int n = rows;
	
		// Use a single 1D vector for L and U for cache optimization
		L->assign(n * n, 0.0);
		U->assign(n * n, 0.0);
	
		// Temporary variables to avoid repeated pointer dereferencing
		std::vector<fT>& L_data = *L;
		std::vector<fT>& U_data = *U;
	
		for (int i = 0; i < n; i++)
		{
			// Compute U's i-th row
			#pragma omp parallel for
			for (int k = i; k < n; k++)
			{
				double sum = 0;
				for (int j = 0; j < i; j++)
				{
					sum += L_data[i * n + j] * U_data[j * n + k];
				}
				U_data[i * n + k] = data[i][k] - sum;
			}
	
			// Compute L's i-th column
			#pragma omp parallel for
			for (int k = i; k < n; k++)
			{
				if (i == k)
				{
					L_data[k * n + i] = 1; // Diagonal of L is 1
				}
				else
				{
					double sum = 0;
					for (int j = 0; j < i; j++)
					{
						sum += L_data[k * n + j] * U_data[j * n + i];
					}
					L_data[k * n + i] = (data[k][i] - sum) / U_data[i * n + i];
				}
			}
		}
	
		// Set the flag to indicate that LU decomposition is done
		LUdecom = true;
	}

	// Sparse LU Decomposition
	// ----------------------------------

	/* 
	Performs the Sparse LU Decomposition of the matrix for further calculations.
	Parameters: 
		CSRMatrix A: CSR Matrix form of the matrix being calculated
		CSRMatrix L: Empty CSR Matrix to store the Lower triangular values
		CSRMatrix U: Empty CSR Matrix to store the Upper triangular values

	Returns:

	Used for large sparse matrices
	*/

	void SparseLU(CSRMatrix &A, CSRMatrix &L, CSRMatrix &U)
	{
		int n = A.n;
	
		// Initialize L as the identity matrix
		for (int i = 0; i < n; i++) {
			L.add_value(i, i, 1.0);
		}
	
		// Perform LU decomposition with partial pivoting
		for (int i = 0; i < n; i++) {
			// Partial pivoting: Find the row with the maximum element in the current column
			int max_row = i;
			double max_val = 0.0;
			for (int j_idx = A.row_ptr[i]; j_idx < A.row_ptr[i + 1]; j_idx++) {
				int j = A.col_idx[j_idx];
				if (j == i) {
					double val = std::abs(A.values[j_idx]);
					if (val > max_val) {
						max_val = val;
						max_row = j_idx;
					}
				}
			}
	
			// Compute U(i, j) and L(i, j)
			for (int j_idx = A.row_ptr[i]; j_idx < A.row_ptr[i + 1]; j_idx++) {
				int j = A.col_idx[j_idx];
				double a_ij = A.values[j_idx];
	
				// Compute U(i, j)
				for (int k_idx = L.row_ptr[i]; k_idx < L.row_ptr[i + 1]; k_idx++) {
					int k = L.col_idx[k_idx];
					if (k < j) {
						a_ij -= L.values[k_idx] * U.values[U.row_ptr[k] + j];
					}
				}
	
				// Update U and L
				if (i <= j) {
					U.add_value(i, j, a_ij);
				} else {
					double l_ij = a_ij / U.values[U.row_ptr[j] + j];
					L.add_value(i, j, l_ij);
				}
			}
		}
	
		L.finalize();
		U.finalize();
	} // end SparseLU

	// Determinant
	// ----------------------------------------

	/* 
	Performs the determinant of the matrix.

	Returns:
		fT determinant: The determinant of the matrix
	*/

	template <typename fT>
	fT Matrix_02508971<fT>::Determinant() const
	{	
		// Check for non square matrix
		CheckSquareMatrix(rows, cols);
		int n = rows;

		// In case of large matrices we assume sparsity. 
		if (n > 500000)
		{	
			cout << "Using Sparse LU! \n";
			// Convert the matrix to CSR format
			CSRMatrix A_CSR(n);
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < n; j++)
				{
					if (data[i][j] != 0)
					{
						A_CSR.add_value(i, j, data[i][j]);
					}
				}
			}
			A_CSR.finalize();

			// Perform Sparse LU Decomposition
			CSRMatrix L_CSR(n), U_CSR(n);
			SparseLU(A_CSR, L_CSR, U_CSR);

			// Compute determinant as product of diagonal elements of U
			double det = 1.0;
			for (int i = 0; i < n; i++)
			{
				for (int idx = U_CSR.row_ptr[i]; idx < U_CSR.row_ptr[i + 1]; idx++)
				{
					if (U_CSR.col_idx[idx] == i)
					{
						det *= U_CSR.values[idx];
						break;
					}
				}
			}
			return det;
		}

		// If small matrix, use standard LU
		if (!LUdecom)
		{
			LUdecomposition();
		}
		else
		{
			cout << "Using precomputed LU Decomposition from previous computation for calculating Determinant!\n";
		}

		double det = 1.0;
		for (int i = 0; i < n; i++)
		{
			det *= (*U)[i * n + i];
		}
		return det;
	} // end Determinant

	// Inverse
	// -----------------------------------

	/* 
	Computes the inverse of the matrix.
	Parameters: 
		Matrix_02508971 A: Matrix to store the inverse

	Returns:
		bool: Returns whether the matrix is invertible. 
	*/

	template <typename fT>
	bool Matrix_02508971<fT>::Inverse(Matrix_02508971 &result) const
	{
		result.Resize(rows, cols);

		// Check for singularity
		double det = Determinant();
		CheckSingular(det);

		int n = rows;

		if (!LUdecom)
		{
			LUdecomposition();
		}
		else
		{
			cout << "Using precomputed LU Decomposition from previous computation for calculating Inverse!\n";
		}

		std::vector<double> invL(n * n, 0.0);
		std::vector<double> invU(n * n, 0.0);

// Compute invL using forward substitution
#pragma omp parallel for
		for (int i = 0; i < n; i++)
		{
			invL[i * n + i] = 1.0;
			for (int j = i + 1; j < n; j++)
			{
				double sum = 0.0;
				for (int k = 0; k < j; k++)
				{
					sum += (*L)[j * n + k] * invL[k * n + i];
				}
				invL[j * n + i] = -sum;
			}
		}

// Compute invU using back substitution
#pragma omp parallel for
		for (int i = n - 1; i >= 0; i--)
		{
			invU[i * n + i] = 1.0 / (*U)[i * n + i];
			for (int j = i - 1; j >= 0; j--)
			{
				double sum = 0.0;
				for (int k = j + 1; k < n; k++)
				{
					sum += (*U)[j * n + k] * invU[k * n + i];
				}
				invU[j * n + i] = -sum / (*U)[j * n + j];
			}
		}

// Compute final inverse: result = invU * invL
#pragma omp parallel for collapse(2)
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				double sum = 0.0;
				for (int k = 0; k < n; k++)
				{
					sum += invU[i * n + k] * invL[k * n + j];
				}
				result.data[i][j] = sum;
			}
		}

		return true;
	} // end Inverse

	// ------------------------------------------------------------------------------
	// template instantiations
	// ------------------------------------------------------------------------------

	// double instantiaition
	template class Matrix_02508971<double>;

	template Matrix_02508971<double> operator+(
		const Matrix_02508971<double> &a,
		const Matrix_02508971<double> &b);

	template Matrix_02508971<double> operator-(
		const Matrix_02508971<double> &a,
		const Matrix_02508971<double> &b);

	template Matrix_02508971<double> operator*(
		const Matrix_02508971<double> &a,
		const Matrix_02508971<double> &b);
	template vector<double> operator*(const Matrix_02508971<double> &mat,
									  const vector<double> &vec);

	// int instantiaition
	template class Matrix_02508971<int>;

	template Matrix_02508971<int> operator+(
		const Matrix_02508971<int> &a,
		const Matrix_02508971<int> &b);

	template Matrix_02508971<int> operator-(
		const Matrix_02508971<int> &a,
		const Matrix_02508971<int> &b);

	template Matrix_02508971<int> operator*(
		const Matrix_02508971<int> &a,
		const Matrix_02508971<int> &b);
	template vector<int> operator*(const Matrix_02508971<int> &mat,
								   const vector<int> &vec);

	// extra operators
	template Matrix_02508971<double64>
	operator*(const vector<double64> &vec, const Matrix_02508971<double64> &mat);
}