//
// Created by brian on 11/20/18.
//

#include "complex.h"

#include <cmath>
#include <iomanip>

const float PI = 3.14159265358979f;

Complex::Complex() : real(0.0f), imag(0.0f) {}

Complex::Complex(float r) : real(r), imag(0.0f) {}

Complex::Complex(float r, float i) : real(r), imag(i) {}

Complex Complex::operator+(const Complex &b) const {
    
    
        Complex sum(0,0);
      
        sum.real = this->real + b.real;
        sum.imag = this->imag + b.imag;

        return sum;

}

Complex Complex::operator-(const Complex &b) const {

    
        Complex sub(0,0);
      
        sub.real = this->real - b.real;
        sub.imag = this->imag - b.imag;

        return sub;

}

Complex Complex::operator*(const Complex &b) const {
    
      Complex mult(0,0);

      mult.real = this->real * b.real - this->imag*b.imag;
      mult.imag = this->real * b.imag + this->imag*b.real;

      return mult;
}

Complex Complex::mag() const {
    /*
    float Magn;
    Magn = sqrt(pow(this->real,2.0)+pow(this->imag,2.0));
    return Magn;
    */
    return Complex(sqrt(real*real + imag*imag));
}

Complex Complex::angle() const {
    double Ang = atan2(this->imag, this->real); // use the arctangent formula
    return Ang;
}

Complex Complex::conj() const {
    Complex Num(0,0);
    Num.real = this->real;
    Num.imag = -1*imag;
    return Num;
}

std::ostream& operator<< (std::ostream& os, const Complex& rhs) {
    Complex c(rhs);
    if(fabsf(rhs.imag) < 1e-10) c.imag = 0.0f;
    if(fabsf(rhs.real) < 1e-10) c.real = 0.0f;

    os << std::fixed << std::setprecision(4);


    if(c.imag == 0) {
        os << c.real;
    }
    else {
        os << "(" << c.real << "," << c.imag << ")";
    }

    return os;
}

// reading
Complex Complex::operator[] (int index) const
{
    return this->m_elements[index];
}

// writing
Complex &Complex::operator[] (int index)
{
    return this->m_elements[index];
} 

// +=
Complex &Complex::operator+=(const Complex& rhs)
{
        this->real += rhs.real;
        this->imag += rhs.imag;
        return *this;  
}