//
// Created by brian on 11/20/18.
//

#pragma once

#include <iostream>
#include <iomanip>

class Complex {
private:
    //int *array_ptr;
    //int arr_size;
    Complex* m_elements; // random number of indices to have in the array
public:
    Complex();
    Complex(float r, float i);
    Complex(float r);
    Complex operator+(const Complex& b) const;
    Complex operator-(const Complex& b) const;
    Complex operator*(const Complex& b) const;

    Complex mag() const;
    Complex angle() const;
    Complex conj() const;

    Complex operator[](int index) const;
    Complex &operator[](int index);


    Complex &operator+=(const Complex& rhs);

    float real;
    float imag;
};

std::ostream& operator<<(std::ostream& os, const Complex& rhs);


