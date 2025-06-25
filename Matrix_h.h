// Step 1: RENAME THE _ADV_PROG_MATRIX_H_ Label To _ADV_PROG_MATRIX_%Your CID Nr%_H_. So, if your CID is 00112233, then the label should be named: _ADV_PROG_MATRIX_00112233_H_
#ifndef _ADV_PROG_MATRIX_02508971_H_
#define _ADV_PROG_MATRIX_02508971_H_

#include <vector>
#include <iostream>

namespace adv_prog_cw 
{
	// Step 2: RENAME THE MATRIX CLASS To Matrix_%Your CID Nr%. So, if your CID is 00112233, then the class should be named: Matrix_00112233
	template<typename fT>
	class Matrix_02508971 {
	public:
		Matrix_02508971();
		Matrix_02508971(size_t m, size_t n);
		Matrix_02508971(size_t m, size_t n, fT val);
		Matrix_02508971(const Matrix_02508971& M);
		~Matrix_02508971();

		size_t Rows() const;
		size_t Cols() const;
		void Resize(size_t m, size_t n);
		// accessors M(i,j)
		fT& operator()(size_t m, size_t n);
		const fT& operator()(size_t m, size_t n) const;
		// assignment
		Matrix_02508971& operator=(const Matrix_02508971& M);
		Matrix_02508971& operator=(fT val);

		Matrix_02508971& operator+=(const Matrix_02508971& M);
		Matrix_02508971& operator-=(const Matrix_02508971& M);
		void Identity();
		void Zero();
		void Transposed(Matrix_02508971& RES) const;
		void Out(long digits = 5L) const;

		// Step 3: Implement the following methods. These will be tested for speed and accuracy for matrices of increasing size and complexity.
		
		// Step 3.1: Multiplication of the matrix by a scalar
		Matrix_02508971& operator*=(fT scalar);
		// Step 3.2:  Division of the matrix by a scalar
		Matrix_02508971& operator/=(fT scalar);
		// Step 3.3:  A method to compute the determinant of square matrices
		void LUdecomposition() const;
		fT Determinant() const;
		// Step 3.4:  A method to compute the inverse of the matrix
		bool Inverse(Matrix_02508971& result) const;

	private:
		// Store the L, U vectors in heap to avoid overfilling the stack, as they scale with the matrix.
		mutable std::vector<fT>* L = new std::vector<fT>;
		mutable std::vector<fT>* U = new std::vector<fT>;

		std::vector<std::vector<fT> >  data;
		size_t                   rows, cols;

		bool  CheckRange(size_t m, size_t n, const char* originator) const;
		bool  CheckSizes(const Matrix_02508971& mat, const char* originator) const;
		bool  CheckSingular(double& det) const;
		bool  CheckSquareMatrix(int rows, int cols) const;

		// If LU decomposition has been saved, we do not compute it again
		mutable bool LUdecom;
	};

	class SingularMatrixException : public std::runtime_error {
		public:
			SingularMatrixException() : std::runtime_error("Matrix is singular and cannot be inverted!") {}
		};
	
	class SquareMatrixException : public std::runtime_error { 
		public:
			SquareMatrixException() : std::runtime_error("Matrix is not square. Cannot calculate determinant!") {}
		};
	
	// Struct to use with Sparse LU Decomposition
	struct CSRMatrix {
		int n;
		std::vector<double> values;
		std::vector<int> col_idx;
		std::vector<int> row_ptr;
	
		CSRMatrix(int size) : n(size), row_ptr(size + 1, 0) {}
	
		void add_value(int row, int col, double val) {
			values.push_back(val);
			col_idx.push_back(col);
			row_ptr[row + 1]++;
		}
	
		void finalize() {
			for (int i = 1; i <= n; ++i)
				row_ptr[i] += row_ptr[i - 1];
		}
	};

	// associated operators
	template<typename fT>
	Matrix_02508971<fT>  operator+(const Matrix_02508971<fT>& a, const Matrix_02508971<fT>& b);

	template<typename fT>
	Matrix_02508971<fT>  operator-(const Matrix_02508971<fT>& a, const Matrix_02508971<fT>& b);

	template<typename fT>
	Matrix_02508971<fT>  operator*(const Matrix_02508971<fT>& a, const Matrix_02508971<fT>& b);

	template<typename fT>
	std::vector<fT>  operator*(const Matrix_02508971<fT>& mat, const std::vector<fT>& v);

	// v^T = (v^T M)^T
	template<typename fT>
	Matrix_02508971<fT>  operator*(const std::vector<fT>& v, const Matrix_02508971<fT>& M);

	// <!-- ==== class Method documentation format ==================== -->

	// constructor (i,j)
	// ---------------------------------------
	//
	template<typename fT>
	inline Matrix_02508971<fT>::Matrix_02508971(size_t m, size_t n)
		: data(m, std::vector<fT>(n)), rows(m), cols(n)
	{
	}

	// operator (i,j)
	// ---------------------------------------
	//
	template<typename fT>
	inline fT& Matrix_02508971<fT>::operator()(size_t m, size_t n)
	{
#ifndef NDEBUG
		CheckRange(m, n, "Matrix_02508971<fT>::operator()");
#endif
		return data[m][n];
	}

	// operator (i,j) const
	// ---------------------------------------
	//
	template<typename fT>
	inline const fT& Matrix_02508971<fT>::operator()(size_t m, size_t n) const
	{
#ifndef NDEBUG
		CheckRange(m, n, "Matrix<fT>::operator()");
#endif
		return data[m][n];
	}

	// Rows()
	// ---------------------------------------------
	//
	template<typename fT>
	inline size_t Matrix_02508971<fT>::Rows() const { return rows; }

	// Cols()
	// ---------------------------------------------
	//
	template<typename fT>
	inline size_t Matrix_02508971<fT>::Cols() const { return cols; }

	// CheckRange(i,j, message )
	// ---------------------------------------
	//
	template<typename fT>
	inline bool Matrix_02508971<fT>::CheckRange(size_t m, size_t n,
		const char* originator) const
	{
		if (m >= rows) {
			std::cerr << "\n" << originator << " row index violation, index=" << m;
			std::cerr << " versus, row-max=" << rows << std::endl;
			throw std::length_error("Matrix<double,mn_max>::CheckRange");
			return false;
		}
		if (n >= cols) {
			std::cerr << "\n" << originator << " column index violation, index=" << n;
			std::cerr << " versus, column-max=" << cols << std::endl;
			throw std::length_error("Matrix<double,mn_max>::CheckRange");
			return false;
		}
		return true;
	}

	// CheckSizes (i,j,message)
	// ---------------------------------------
	template<typename fT>
	inline bool Matrix_02508971<fT>::CheckSizes(const Matrix_02508971& mat,
		const char* originator) const
	{
		if (rows != mat.rows) {
			std::cerr << "\n" << originator << " matrices have different sizes; rows1=" << rows;
			std::cerr << " versus, rows2=" << mat.rows << std::endl;
			throw std::length_error("Matrix_02508971<double,mn_max>::CheckSizes");
			return false;
		}
		if (cols != mat.cols) {
			std::cerr << "\n" << originator << " matrices have different sizes; columns1=" << cols;
			std::cerr << " versus, columns2=" << mat.cols << std::endl;
			throw std::length_error("Matrix<double,mn_max>::CheckSizes");
			return false;
		}
		return true;
	}

	template<typename fT>
	inline bool Matrix_02508971<fT>::CheckSingular(double& det) const{
		if(det == 0){
			std::cerr << "\n Matrix is has zero determinant (singular), cannot compute inverse.";
			throw SingularMatrixException();
			return false;
		}
		return true;
	}

	template <typename fT>
	inline bool Matrix_02508971<fT>::CheckSquareMatrix(int rows, int cols) const{
		if(rows != cols){
			std::cerr << "\n The matrix is not square! Cannot compute determinant.";
			throw SquareMatrixException();
			return false;
		}
		return true;
	}
} // end scope

#endif
