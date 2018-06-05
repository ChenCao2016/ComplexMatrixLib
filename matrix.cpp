#include "stdafx.h"
#include "matrix.h"
#include <cmath>



static void printMatrix(const matrix* a) {
	std::complex<double> b;
	for (unsigned int m = 0; m < a->columnLen; m++) {
		for (unsigned int n = 0; n < a->rowLen; n++) {
			b = a->value(m, n);
			if (n == 0) {
				printf("\n");
			}
			printf("%f,%fj	", real(b), imag(b));
		}
	}
	printf("\n");
}

static double checkMaxWeight(const matrix* a) {
	double maxWeight = 0;
	for (unsigned int m = 0; m < a->columnLen; m++) {
		for (unsigned int n = 0; n < a->rowLen; n++) {
			if (m != n) {
				double weight = std::norm(a->value(m, n));
				if (weight > maxWeight) {
					maxWeight = weight;
				}
			}
		}
	}
	return maxWeight;
}

matrix::matrix(unsigned int columnLen, unsigned int rowLen):columnLen(columnLen),rowLen(rowLen) {
	element = new std::complex<double>[columnLen*rowLen];
}


matrix::~matrix() {
	delete[] element;
}



void matrix::setMatrix(const double* a) {
	for (unsigned int m = 0; m < columnLen; m++) {
		for (unsigned int n = 0; n < rowLen; n++) {
			setValue(m, n, std::complex<double>(*a,*(a+1)));
			a = a + 2;
		}
	}
}

void matrix::queryMatrix(double* a) const{
	for (unsigned int m = 0; m < columnLen; m++) {
		for (unsigned int n = 0; n < rowLen; n++) {
			*a = std::real(this->value(m, n));
			*(a+1) = std::imag(this->value(m, n));
			a = a + 2;
		}
	}
}

std::complex<double> matrix::value(unsigned int columnInd, unsigned int rowInd) const{
	if (columnInd >= columnLen || rowInd >= rowLen) {
		throw ("ERROR: index excceds matrix size");
	}
	return element[columnInd*rowLen + rowInd];
}

void matrix::setValue(unsigned int columnInd, unsigned int rowInd, std::complex<double> cvalue) {
	if (columnInd >= columnLen || rowInd >= rowLen) {
		throw ("ERROR: index excceds matrix size");
	}
	element[columnInd*rowLen + rowInd] = cvalue;
}

void matrix::setValue(const std::complex<double>* cvalues) {
	for (unsigned int m = 0; m < columnLen; m++) {
		for (unsigned int n = 0; n < rowLen; n++) {
			setValue(m, n, cvalues[m*rowLen + n]);
		}
	}
}


void matrix::setValue(const matrix* a) {

	if (columnLen != a->columnLen || rowLen != a->rowLen) {
		throw("ERROR: matrix size does not match");
	}

	for (unsigned int m = 0; m < columnLen; m++) {
		for (unsigned int n = 0; n < rowLen; n++) {
			setValue(m, n, a->value(m,n));
		}
	}
}


int matrixOperation::transpose(const matrix* a, matrix* b) {
	int err = 0;
	if (a->columnLen != b->rowLen || a->rowLen != b->columnLen) {
		throw ("ERROR: matrix size does not match");
	}
	
	for (unsigned int m = 0; m < a->columnLen; m++) {
		for (unsigned int n = 0; n < a->rowLen; n++){
			b->setValue(n, m, a->value(m, n));
		}
	}
	return err;
}

int matrixOperation::conjugateTranspose(const matrix* a, matrix* b) {
	int err = 0;

	if (a->columnLen != b->rowLen || a->rowLen != b->columnLen) {
		throw ("ERROR: matrix size does not match");
	}

	for (unsigned int m = 0; m < a->columnLen; m++) {
		for (unsigned int n = 0; n < a->rowLen; n++) {
			std::complex<double> temp = a->value(m, n);
			//b->setValue(n, m, std::complex<double>(std::real(temp),-std::imag(temp)));
			b->setValue(n, m, std::conj(temp));
		}
	}

	return err;
}

int matrixOperation::identity(matrix* a) {
	int err = 0;

	for (unsigned int m = 0; m < a->columnLen; m++) {
		for (unsigned int n = 0; n < a->rowLen; n++) {
			if (m == n) {
				a->setValue(m, n, 1);
			}
			else {
				a->setValue(m, n, 0);
			}
		}
	}
	return err;
}

int matrixOperation::multiply(const matrix* a, const matrix* b, matrix* c) {
	int err = 0;

	if (a->rowLen != b->columnLen) {
		throw ("ERROR: matrix size does not match");
	}

	for (unsigned int m = 0; m < a->columnLen; m++) {
		for (unsigned int k = 0; k < b->rowLen; k++) {
			std::complex<double> temp(0, 0);
			for (unsigned int n = 0; n < a->rowLen; n++) {
				temp += (a->value(m, n))*(b->value(n, k));

				//double realPart = std::real(a->value(m, n))*std::real(b->value(n, k))
				//				- std::imag(a->value(m, n))*std::imag(b->value(n, k));
				//double imagPart = std::real(a->value(m, n))*std::imag(b->value(n, k))
				//				+ std::imag(a->value(m, n))*std::real(b->value(n, k));

				//temp = std::complex<double>(std::real(temp) + realPart, std::imag(temp) + imagPart);

			}
			c->setValue(m, k, temp);
		}
	}

	return err;
}

int matrixOperation::GivensRotation(const matrix* a, unsigned int m, unsigned int n, unsigned int k, unsigned int p, matrix* b) {
	int err = 0;

	if (b->columnLen != a->columnLen || b->rowLen != a->columnLen) {
		throw ("ERROR: matrix size does not match");
	}

	if (m>=a->columnLen || k>= a->columnLen || n>= a->rowLen||p>= a->rowLen) {
		throw ("ERROR: index excceds matrix size");
	}

	std::complex<double> a_mn = a->value(m, n);
	std::complex<double> a_kp = a->value(k, p);

	//double magnitude = sqrt(std::real(a_mn)*std::real(a_mn) +
	//						std::imag(a_mn)*std::imag(a_mn) + 
	//						std::real(a_kp)*std::real(a_kp) + 
	//						std::imag(a_kp)*std::imag(a_kp));

	double magnitude = sqrt(std::norm(a_mn) + std::norm(a_kp));

	//std::complex<double> c(std::real(a_mn) / magnitude, -std::imag(a_mn) / magnitude);
	//std::complex<double> s(std::real(a_kp) / magnitude, -std::imag(a_kp) / magnitude);

	std::complex<double> c = std::conj(a_mn) / magnitude;
	std::complex<double> s = std::conj(a_kp) / magnitude;

	identity(b);

	b->setValue(m, m, c);
	b->setValue(k, k, std::complex<double>(std::real(c), -std::imag(c)));

	b->setValue(m, k, s);
	b->setValue(k, m, std::complex<double>(-std::real(s), std::imag(s)));

	return err;
}

int matrixOperation::QRDecomposition(const matrix* a, matrix* q, matrix* r) {
	int err = 0;
	double zero = 0;

	matrix temp(a->columnLen, a->columnLen);
	matrix temp_p0(a->columnLen, a->columnLen);
	matrix temp_p1(a->columnLen, a->columnLen);

	matrixOperation::identity(&temp_p0);
	int flag_temp_p = 0;

	for (unsigned int m = 0; m < a->columnLen - 1 && m < a->rowLen; m++) {
		for (unsigned int n = a->columnLen - 1; n > m; n--) {
			
			if (m == 0 && n == a->columnLen - 1) {

				if (abs(std::real(a->value(n, m)))> zero || 
					abs(std::imag(a->value(n, m)))> zero) 
				{
					matrixOperation::GivensRotation(a, n - 1, m, n, m, &temp);
					matrixOperation::multiply(&temp, &temp_p0, &temp_p1);
					matrixOperation::multiply(&temp_p1, a, r);
					flag_temp_p = 1;
				}

			}
			else {
				if (abs(std::real(r->value(n, m))) > zero ||
					abs(std::imag(r->value(n, m))) > zero) 
				{
					matrixOperation::GivensRotation(r, n - 1, m, n, m, &temp);
					if (flag_temp_p == 0) {
						matrixOperation::multiply(&temp, &temp_p0, &temp_p1);
						matrixOperation::multiply(&temp_p1, a, r);
						flag_temp_p = 1;
					}
					else if (flag_temp_p == 1) {
						matrixOperation::multiply(&temp, &temp_p1, &temp_p0);
						matrixOperation::multiply(&temp_p0, a, r);
						flag_temp_p = 0;
					}
				}
			}
		}
	}

	if (a->columnLen == 1) {
		identity(q);
		r->setValue(a);
	}

	if (flag_temp_p == 0) {
		matrixOperation::conjugateTranspose(&temp_p0, q);
	}

	if (flag_temp_p == 1) {
		matrixOperation::conjugateTranspose(&temp_p1, q);
	}

	return err;
}

int matrixOperation::eigenDecomposition(const matrix* a, unsigned int iteration, std::complex<double>* lambda) {
	int err = 0;

	if (a->columnLen != a->rowLen) {
		throw ("ERROR: matrix is not squared");
	}

	matrix temp_a(a->columnLen, a->columnLen);
	matrix q(a->columnLen, a->columnLen);
	matrix r(a->columnLen, a->columnLen);

	matrixOperation::QRDecomposition(a, &q, &r);
	for (unsigned int i = 1; i < iteration; i++) {
		matrixOperation::multiply(&r, &q, &temp_a);
		matrixOperation::QRDecomposition(&temp_a, &q, &r);

		if (checkMaxWeight(&temp_a) < g_zeroThresholdForEigenSingluarDecomposition) {
			break;
		}
	}
	matrixOperation::multiply(&r, &q, &temp_a);
	

	for (unsigned int m = 0; m < a->columnLen; m++) {
		lambda[m] = temp_a.value(m, m);
	}

	if (checkMaxWeight(&temp_a) > g_zeroThresholdForEigenSingluarDecomposition) {
		err = -1;
	}

	return err;
}

int matrixOperation::eigenDecomposition(const matrix*a, std::complex<double>* lambda) {
	int err = 0;
	
	err = matrixOperation::eigenDecomposition(a, g_iterationForQRMethod, lambda);

	return err;
}


int matrixOperation::eigenDecomposition(const matrix* a, unsigned int iteration, std::complex<double>* lambda, matrix* u) {
	int err = 0;

	if (a->columnLen != a->rowLen) {
		throw ("ERROR: matrix is not squared");
	}

	matrix temp_a(a->columnLen, a->columnLen);
	matrix q(a->columnLen, a->columnLen);
	matrix r(a->columnLen, a->columnLen);
	matrix temp_u0(a->columnLen, a->columnLen);
	matrix temp_u1(a->columnLen, a->columnLen);

	matrixOperation::QRDecomposition(a, &q, &r);
	temp_u0.setValue(&q);
	int flag_temp_u = 0;

	for (unsigned int i = 1; i < iteration; i++) {
		matrixOperation::multiply(&r, &q, &temp_a);
		matrixOperation::QRDecomposition(&temp_a, &q, &r);
		if (flag_temp_u == 0) {
			matrixOperation::multiply(&temp_u0, &q, &temp_u1);
			flag_temp_u = 1;
		}
		else if (flag_temp_u == 1) {
			matrixOperation::multiply(&temp_u1, &q, &temp_u0);
			flag_temp_u = 0;
		}

		if (checkMaxWeight(&temp_a) < g_zeroThresholdForEigenSingluarDecomposition) {
			break;
		}
	}
	matrixOperation::multiply(&r, &q, &temp_a);

	for (unsigned int m = 0; m < a->columnLen; m++) {
		lambda[m] = temp_a.value(m, m);
	}

	if (flag_temp_u == 0) {
		u->setValue(&temp_u0);
	}
	if (flag_temp_u == 1) {
		u->setValue(&temp_u1);
	}


	if (checkMaxWeight(&temp_a) > g_zeroThresholdForEigenSingluarDecomposition) {
		err = -1;
	}

	return err;
}

int matrixOperation::eigenDecomposition(const matrix*a, std::complex<double>* lambda, matrix* u) {
	int err = 0;

	err = matrixOperation::eigenDecomposition(a, g_iterationForQRMethod, lambda,u);

	return err;
}

int matrixOperation::eigenDecomposition(const matrix* a, double* lambda, unsigned int* length, matrix* u) {
	int err = 0;

	if (a->columnLen != a->rowLen) {
		throw ("ERROR: matrix is not squared");
	}

	matrix temp_a(a->columnLen, a->columnLen);
	matrix q(a->columnLen, a->columnLen);
	matrix r(a->columnLen, a->columnLen);
	matrix temp_u0(a->columnLen, a->columnLen);
	matrix temp_u1(a->columnLen, a->columnLen);

	matrixOperation::QRDecomposition(a, &q, &r);
	temp_u0.setValue(&q);
	int flag_temp_u = 0;

	for (unsigned int i = 1; i < g_iterationForQRMethod; i++) {
		matrixOperation::multiply(&r, &q, &temp_a);
		matrixOperation::QRDecomposition(&temp_a, &q, &r);
		if (flag_temp_u == 0) {
			matrixOperation::multiply(&temp_u0, &q, &temp_u1);
			flag_temp_u = 1;
		}
		else if (flag_temp_u == 1) {
			matrixOperation::multiply(&temp_u1, &q, &temp_u0);
			flag_temp_u = 0;
		}

		if (checkMaxWeight(&temp_a) < g_zeroThresholdForEigenSingluarDecomposition) {
			break;
		}
	}
	matrixOperation::multiply(&r, &q, &temp_a);

	*length = 0;
	for (unsigned int m = 0; m < a->columnLen; m++) {
		*lambda = std::real(temp_a.value(m, m));
		*(lambda + 1) = std::imag(temp_a.value(m, m));
		lambda = lambda + 2;
		*length++ ;
	}

	if (flag_temp_u == 0) {
		u->setValue(&temp_u0);
	}
	if (flag_temp_u == 1) {
		u->setValue(&temp_u1);
	}

	if (checkMaxWeight(&temp_a) > g_zeroThresholdForEigenSingluarDecomposition) {
		err = -1;
	}

	return err;
}


int matrixOperation::singluarDecomposition(const matrix*a, unsigned int iteration, std::complex<double>* lambda) {
	int err = 0;

	matrix temp_a_1(a->columnLen, a->rowLen);
	matrix q_1(a->columnLen, a->columnLen);
	matrix r_1(a->columnLen, a->rowLen);

	matrix temp_a_2(a->rowLen, a->columnLen);
	matrix q_2(a->rowLen, a->rowLen);
	matrix r_2(a->rowLen, a->columnLen);

	matrixOperation::QRDecomposition(a, &q_1, &r_1);

	matrixOperation::conjugateTranspose(a, &temp_a_2);
	matrixOperation::QRDecomposition(&temp_a_2, &q_2, &r_2);


	for (unsigned int i = 1; i < iteration; i++) {
		matrixOperation::multiply(&r_1, &q_2, &temp_a_1);
		matrixOperation::multiply(&r_2, &q_1, &temp_a_2);
		matrixOperation::QRDecomposition(&temp_a_1, &q_1, &r_1);
		matrixOperation::QRDecomposition(&temp_a_2, &q_2, &r_2);

		if (checkMaxWeight(&temp_a_1) < g_zeroThresholdForEigenSingluarDecomposition) {
			break;
		}

	}

	matrixOperation::multiply(&r_1, &q_2, &temp_a_1);
	matrixOperation::multiply(&r_2, &q_1, &temp_a_2);


	if (a->rowLen > a->columnLen) {
		for (unsigned int m = 0; m < a->columnLen; m++) {
			lambda[m] = temp_a_1.value(m, m);
		}
	}
	else {
		for (unsigned int m = 0; m < a->rowLen; m++) {
			lambda[m] = temp_a_1.value(m, m);
		}
	}


	if (checkMaxWeight(&temp_a_1) > g_zeroThresholdForEigenSingluarDecomposition) {
		err = -1;
	}

	return err;
}

int matrixOperation::singluarDecomposition(const matrix* a, std::complex<double>* lambda) {
	int err = 0;

	err = matrixOperation::singluarDecomposition(a, g_iterationForQRMethod, lambda);

	return err;
}


int matrixOperation::singluarDecomposition(const matrix*a, unsigned int iteration, std::complex<double>* lambda, matrix* u,matrix* v) {
	int err = 0;

	matrix temp_a_1(a->columnLen, a->rowLen);
	matrix q_1(a->columnLen, a->columnLen);
	matrix r_1(a->columnLen, a->rowLen);

	matrix temp_a_2(a->rowLen, a->columnLen);
	matrix q_2(a->rowLen, a->rowLen);
	matrix r_2(a->rowLen, a->columnLen);

	matrix temp_u0(a->columnLen, a->columnLen);
	matrix temp_u1(a->columnLen, a->columnLen);

	matrix temp_v0(a->rowLen, a->rowLen);
	matrix temp_v1(a->rowLen, a->rowLen);


	matrixOperation::QRDecomposition(a, &q_1, &r_1);
	temp_u0.setValue(&q_1);
	
	
	
	matrixOperation::conjugateTranspose(a, &temp_a_2);
	matrixOperation::QRDecomposition(&temp_a_2, &q_2, &r_2);
	temp_v0.setValue(&q_2);

	//printMatrix(&q_1);
	//printMatrix(&r_2);

	int flag_temp_u = 0;

	for (unsigned int i = 1; i < iteration; i++) {
		matrixOperation::multiply(&r_1, &q_2, &temp_a_1);
		matrixOperation::multiply(&r_2, &q_1, &temp_a_2);

		//printMatrix(&temp_a_1);
		//printMatrix(&temp_a_2);

		matrixOperation::QRDecomposition(&temp_a_1, &q_1, &r_1);
		matrixOperation::QRDecomposition(&temp_a_2, &q_2, &r_2);

		if (flag_temp_u == 0) {
			matrixOperation::multiply(&temp_u0, &q_1, &temp_u1);
			matrixOperation::multiply(&temp_v0, &q_2, &temp_v1);
			flag_temp_u = 1;
		}
		else if (flag_temp_u == 1) {
			matrixOperation::multiply(&temp_u1, &q_1, &temp_u0);
			matrixOperation::multiply(&temp_v1, &q_2, &temp_v0);
			flag_temp_u = 0;
		}

		if (checkMaxWeight(&temp_a_1) < g_zeroThresholdForEigenSingluarDecomposition) {
			break;
		}

	}
	matrixOperation::multiply(&r_1, &q_2, &temp_a_1);
	matrixOperation::multiply(&r_2, &q_1, &temp_a_2);

	//printMatrix(&temp_a_1);
	//printMatrix(&temp_a_2);

	if (a->rowLen > a->columnLen) {
		for (unsigned int m = 0; m < a->columnLen; m++) {
			lambda[m] = temp_a_1.value(m, m);
		}
	}
	else {
		for (unsigned int m = 0; m < a->rowLen; m++) {
			lambda[m] = temp_a_1.value(m, m);
		}
	}

	if (flag_temp_u == 0) {
		u->setValue(&temp_u0);
		v->setValue(&temp_v0);	
	}
	if (flag_temp_u == 1) {
		u->setValue(&temp_u1);
		v->setValue(&temp_v1);
	}

	if (checkMaxWeight(&temp_a_1) > g_zeroThresholdForEigenSingluarDecomposition) {
		err = -1;
	}

	return err;
}

int matrixOperation::singluarDecomposition(const matrix*a, std::complex<double>* lambda, matrix* u, matrix* v) {
	int err = 0;

	err = matrixOperation::singluarDecomposition(a, g_iterationForQRMethod, lambda, u, v);

	return err;
}

int matrixOperation::singluarDecomposition(const matrix* a, double* lambda, unsigned int* length, matrix* u, matrix* v) {
	int err = 0;

	matrix temp_a_1(a->columnLen, a->rowLen);
	matrix q_1(a->columnLen, a->columnLen);
	matrix r_1(a->columnLen, a->rowLen);

	matrix temp_a_2(a->rowLen, a->columnLen);
	matrix q_2(a->rowLen, a->rowLen);
	matrix r_2(a->rowLen, a->columnLen);

	matrix temp_u0(a->columnLen, a->columnLen);
	matrix temp_u1(a->columnLen, a->columnLen);

	matrix temp_v0(a->rowLen, a->rowLen);
	matrix temp_v1(a->rowLen, a->rowLen);


	matrixOperation::QRDecomposition(a, &q_1, &r_1);
	temp_u0.setValue(&q_1);



	matrixOperation::conjugateTranspose(a, &temp_a_2);
	matrixOperation::QRDecomposition(&temp_a_2, &q_2, &r_2);
	temp_v0.setValue(&q_2);

	int flag_temp_u = 0;

	for (unsigned int i = 1; i < g_iterationForQRMethod; i++) {
		matrixOperation::multiply(&r_1, &q_2, &temp_a_1);
		matrixOperation::multiply(&r_2, &q_1, &temp_a_2);

		matrixOperation::QRDecomposition(&temp_a_1, &q_1, &r_1);
		matrixOperation::QRDecomposition(&temp_a_2, &q_2, &r_2);

		if (flag_temp_u == 0) {
			matrixOperation::multiply(&temp_u0, &q_1, &temp_u1);
			matrixOperation::multiply(&temp_v0, &q_2, &temp_v1);
			flag_temp_u = 1;
		}
		else if (flag_temp_u == 1) {
			matrixOperation::multiply(&temp_u1, &q_1, &temp_u0);
			matrixOperation::multiply(&temp_v1, &q_2, &temp_v0);
			flag_temp_u = 0;
		}

		if (checkMaxWeight(&temp_a_1) < g_zeroThresholdForEigenSingluarDecomposition) {
			break;
		}

	}
	matrixOperation::multiply(&r_1, &q_2, &temp_a_1);
	matrixOperation::multiply(&r_2, &q_1, &temp_a_2);

	if (a->rowLen > a->columnLen) {
		*length = 0;
		for (unsigned int m = 0; m < a->columnLen; m++) {
			*lambda = std::real(temp_a_1.value(m, m));
			*(lambda + 1) = std::imag(temp_a_1.value(m, m));
			lambda = lambda + 2;
			*length++;
		}
	}
	else {
		*length = 0;
		for (unsigned int m = 0; m < a->columnLen; m++) {
			*lambda = std::real(temp_a_1.value(m, m));
			*(lambda + 1) = std::imag(temp_a_1.value(m, m));
			lambda = lambda + 2;
			*length++;
		}
	}

	if (flag_temp_u == 0) {
		u->setValue(&temp_u0);
		v->setValue(&temp_v0);
	}
	if (flag_temp_u == 1) {
		u->setValue(&temp_u1);
		v->setValue(&temp_v1);
	}

	if (checkMaxWeight(&temp_a_1) > g_zeroThresholdForEigenSingluarDecomposition) {
		err = -1;
	}

	return err;
}

int matrixOperation::pseudoInverse(const matrix* a, matrix* b) {
	int err = 0;
	
	double zero = 1e-20;

	if (b->columnLen != a->rowLen || b->rowLen != a->columnLen) {
		throw ("ERROR: matrix size does not match");
	}

	matrix u_ct(a->columnLen, a->columnLen);
	matrix v(a->rowLen, a->rowLen);

	unsigned int len = a->rowLen > a->columnLen ? a->columnLen : a->rowLen;

	std::complex<double> *lambda = new std::complex<double>[len];



	matrix u(a->columnLen, a->columnLen);

	err = matrixOperation::singluarDecomposition(a,lambda, &u, &v);

	matrixOperation::conjugateTranspose(&u, &u_ct);


	for (unsigned int m = 0; m < len; m++) {
		
		//double magnitude = std::real(lambda[m])* std::real(lambda[m]) +
		//				   std::imag(lambda[m])* std::imag(lambda[m]);

		double magnitude = std::norm(lambda[m]);

		if (magnitude > zero) {
			//lambda[m] = std::complex<double>(std::real(lambda[m]) / magnitude, -std::imag(lambda[m]) / magnitude);
			lambda[m] = std::conj(lambda[m]) / magnitude;
		}
		else {
			lambda[m] = std::complex<double>(0, 0);
		}

	}

	matrix eigen(a->rowLen, a->columnLen);
	matrixOperation::identity(&eigen);

	for (unsigned int k = 0; k < len; k++) {
		eigen.setValue(k, k, lambda[k]);
	}



	matrix temp(a->rowLen, a->columnLen);

	matrixOperation::multiply(&v, &eigen, &temp);
	matrixOperation::multiply(&temp, &u_ct, b);


	return err;
}