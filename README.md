# Matrix Class

## Overview

This project extends a provided Matrix class by implementing four additional functions focused on efficiency, accuracy, and scalability for matrix operations.

## Evaluation Criteria

1. **Efficiency** – Optimized code for performance, particularly for large matrices
2. **Accuracy** – Correct implementation of matrix operations  
3. **Scalability** – Ability to handle increasing matrix sizes effectively

## Implementation Requirements

#### 1 Scalar Multiplication
```cpp
Matrix& operator*=(fT scalar);
```
- Multiplies each matrix element by a scalar
- Returns reference to the modified matrix

#### 2 Scalar Division
```cpp
Matrix& operator/=(fT scalar);
```
- Divides each matrix element by a scalar
- Returns reference to the modified matrix
- **Must handle division by zero appropriately**

#### 3 Determinant Calculation
```cpp
fT Determinant() const;
```
- Computes and returns the determinant of a square matrix
- **Must handle non-square matrices appropriately**

#### 4 Matrix Inverse
```cpp
bool Inverse(Matrix& result) const;
```
- Computes matrix inverse and stores result in provided matrix
- Returns `true` if matrix is invertible, `false` otherwise
- **Must handle singular matrices appropriately**

## Features

### Edge Case Handling
- Division by zero
- Singular matrices
- Non-square matrices (where applicable)

### Performance Considerations

Focus on optimization for:
- Large matrix operations
- Memory efficiency
- Computational complexity
- Parallel processing (OpenMP)

---