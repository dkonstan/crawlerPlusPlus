#include <memory>
#include <cassert>

class Complex
{
public:
	// special form, initalizer list
	Complex(double re, double im): real(re), imag(im) {}
	Complex(double re): real(re), imag(0) {}
	Complex(): real(0), imag(0) {}
	Complex(const Complex& c): real(c.real), imag(c.imag) {}  // a copy constructor (made by default to copy all members if not included)
	Complex operator+(const Complex& c2) const
	//                c2 doesn't change   the original Complex object doesn't change
	{
		return Complex(real + c2.real, imag + c2.imag);
	}
	~Complex();
	// friend double& real(Complex& c); // this allows the function "real", defined elsewhere, to access the private members of this class
private:
	double real;
	double imag;
friend double real(Complex c);
friend double imag(Complex c);
};



//structs make everything public by default
//classes make everything private by default

// implement binary operators as free functions
double real(Complex c)
{
	return c.real;
}

double imag(Complex c)
{
	return c.imag;
}
inline Complex operator+(const Complex& c1, const Complex& c2)
{
	return Complex(real(c2) + real(c2), imag(c1) + imag(c2));
}


class Vector
{
public:
	Vector(const Vector& v): my_size(v.my_size), data(new double[my_size]) // this ensures that the copy actually copies the value, not just reference as by default
	{
		for (unsigned i = 0; i < my_size; ++i)
		{
			data[i] = v.data[i];
		}
	}

	// move semantics (how to assign a second variable to the vector object)
	// && means rvalue reference (?)
	Vector(Vector&& v): my_size(v.my_size), data(v.data)  // move constructor steals data and makes rvalue's data 0/null
	{
		v.data = nullptr;  // robbing v.data of its heap data ensures that when the rvalue is cleared nothing happens to the data
		v.my_size = 0;
	}
	// move assignment, v1 = f(), where f returns v2, need to transfer data from v2 to v1
	Vector& operator=(Vector&& src)  // && means reference to an rvalue
	{
		assert(my_size == 0 || my_size == src.my_size);
		std::swap(data, src.data);
		return *this;
	}
	// std::move casts an lvalue to an rvalue
	// Vector x(std::move(w)), in this, x steals the data of w
	// v = std::move(u), this swaps the data of v and u
	~Vector()
	{
		// noexcept by default
		// must not throw exceptions, they will not be caught most likely and the program will crash
		delete[] data;
	}
	const double& operator[](int i) const  // not modifying data (const on left), and 
	{
		assert(i >= 0 && i < my_size);
		return data[i];
	}

private:
	unsigned my_size;
	double *data;
};
