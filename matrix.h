
#ifndef __MATRIX_H
#define __MATRIX_H

#include <complex>

class matrix
{
	std::complex<double>* element;

public:	
	//--------------------------------------------------------------
	//complex matrix interface
	//--------------------------------------------------------------
	//construtor, specify the size of the matrix
	//m: input column vector length
	//n: input row vector length
	matrix(unsigned int m, unsigned int n);

	//desctructor
	~matrix();

	//move memory start from a to this complex matrix
	//a: pointer to input memory 
	void setMatrix(const double* a);
	//copy this complex matrix to the memory start from a
	//a: pointer to ouput memory 
	void queryMatrix(double* a) const;
	//----------------------------------------------------------------


	const unsigned int columnLen;
	const unsigned int rowLen;
	//fetch one scalar value in the matrix
	//m: input column index
	//n: input row index
	//output: complex<double> value at index (m,n)
	std::complex<double> value(unsigned int m, unsigned int n) const;

	//set one scalar value in the matirx
	//m: input column index
	//n: input row index
	//a: input scalar value
	void setValue(unsigned int m, unsigned int n, std::complex<double> a);

	//set matrix from array, length of the array has to match the size of matrix
	//a: input complex<double> array
	void setValue(const std::complex<double>*);

	//set matrix from another matrix, give a copy of the matrix. Two matrix have to the same size
	//a: input matrix to copy
	void setValue(const matrix* a);
};


namespace matrixOperation
{
	//----------------------------------------------------------
	//complex matrix operation interface
	//----------------------------------------------------------
	//matrix transpose
	//a: input matrix
	//b: output transpose matrix
	int transpose(const matrix* a, matrix* b);

	//matrix conjugate transpose
	//a: input matrix
	//b: output conjugate transpose matrix
	int conjugateTranspose(const matrix* a, matrix* b);

	//set matrix as indentity
	//a: input and output identity matrix
	int identity(matrix* a);

	//matrix multipy
	//a: input left matrix
	//b: input right matrix
	//c: output result matrix
	int multiply(const matrix* a, const matrix* b, matrix* c);
	
	//QR decomposition
	//a: input matrix
	//q: output left unitary matrix Q
	//r: output upper triangular matrix R
	int QRDecomposition(const matrix* a, matrix* q, matrix* r);

	//Eigenvector decomposition
	//a: input matrix
	//lambda: output array for eigenvalues
	//length: output eigenvalue array length
	//q: output eigenvector matrix
	int eigenDecomposition(const matrix* a, double* lambda, unsigned int* length, matrix* q);

	//singluar vector decomposition
	//a: input matrix
	//lambda: output array for singular value array
	//length: output eigenvalue array length
	//u: output left unitary matrix U
	//v: output right unitary matrix V
	int singluarDecomposition(const matrix* a, double* lambda, unsigned int* length, matrix* u, matrix* v);

	//Pseudo-Inverse (Moore–Penrose Inverse)
	//a: input matrix
	//b: output pseudoinverse matrix
	int pseudoInverse(const matrix* a, matrix* b);
	//----------------------------------------------------------
	
	//maximum number of iterations used by QR method
	static const unsigned int g_iterationForQRMethod = 40;

	//threshold check for QR method. 
	static const double g_zeroThresholdForEigenSingluarDecomposition = 1e-12;
	
	//Givens Rotation
	int GivensRotation(const matrix*, unsigned int, unsigned int, unsigned int, unsigned int, matrix*);

	//Eigenvalue decomposition
	//a: input matrix
	//iterationForQRmethod: input maximum number of iterations for QR method
	//lambda: output complex<double> array for eigen values
	int eigenDecomposition(const matrix* a, unsigned int iterationForQRMethod, std::complex<double>* lambda);
	int eigenDecomposition(const matrix* a, std::complex<double>* lambda);

	//Eigenvector decomposition
	//a: input matrix
	//iterationForQRmethod: input maximum number of iterations for QR method
	//lambda: output complex<double> array for eigen values
	//q: output eigenvector matrix
	int eigenDecomposition(const matrix* a, unsigned int iterationForQRMethod, std::complex<double>* lambda, matrix* q);
	int eigenDecomposition(const matrix* a, std::complex<double>* lambda, matrix* q);


	//singluar value decomposition
	//a: input matrix
	//iterationForQRmethod: input maximum number of iterations for QR method
	//lambda: output complex<double> array for singluar
	int singluarDecomposition(const matrix* a, unsigned int iterationForQRMethod, std::complex<double>* lambda);
	int singluarDecomposition(const matrix* a, std::complex<double>* lambda);

	//singluar vector decomposition
	//a: input matrix
	//iterationForQRmethod: input maximum number of iterations for QR method
	//lambda: output complex<double> array for singluar
	//u: output left unitary matrix U
	//v: output right unitary matrix V
	int singluarDecomposition(const matrix* a, unsigned int iterationForQRMethod, std::complex<double>* lambda, matrix* u, matrix* v);
	int singluarDecomposition(const matrix* a, std::complex<double>* lambda, matrix* u, matrix* v);

};
#endif