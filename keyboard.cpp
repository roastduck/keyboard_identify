#include <cmath>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include "lodepng.h"

const char filename[] = "input.png";

struct Pixel // Ignore alpha
{
	int r, g, b;

	Pixel() {}
	Pixel(int _r, int _g, int _b) : 
		r(_r), g(_g), b(_b) {}
	Pixel(const Pixel &from) :
		r(from.r), g(from.g), b(from.b) {}

	double norm() const
	{
		return sqrt(r*r+b*b+g*g);
	}
};

inline Pixel operator-(const Pixel &a, const Pixel &b)
{
	return Pixel(a.r-b.r, a.g-b.g, a.b-b.b);
}

template <class T>
class Vec2D : public std::vector<T>
{
	size_t width, height;

public:
	size_t getWidth() const { return width; }
	size_t getHeight() const { return height; }

	T &at(size_t x, size_t y) { return std::vector<T>::at(x*width+y); }
	const T &at(size_t x, size_t y) const { return std::vector<T>::at(x*width+y); }

	Vec2D(size_t _width, size_t _height) :
		std::vector<T>(_width*_height), width(_width), height(_height) {}
};

Vec2D<Pixel> load(const char *filename)
{
	unsigned width, height;
	std::vector<unsigned char> raw;
	unsigned error = lodepng::decode(raw, width, height, filename);
	if (error)
		std::cerr << "Error reading " << filename << " : " << lodepng_error_text(error) << std::endl;
	Vec2D<Pixel> ret(width, height);
	for (size_t i=0; i*4<raw.size(); i++)
		ret[i] = Pixel(raw[i*4], raw[i*4+1], raw[i*4+2]);
	return ret;
}

class ColDiffSum
{
	const Vec2D<Pixel> &image;
	Vec2D<double> diff;

public:
	ColDiffSum(const Vec2D<Pixel> &_image) :
		image(_image), diff(_image.getWidth(), _image.getHeight())
	{
		for (size_t i=0; i<image.getHeight(); i++)
			for (size_t j=0; j<image.getWidth(); j++)
				diff.at(i,j) = (j ? image.at(i,j)-image.at(i,j-1) : image.at(i,j)).norm();
	}
};

#ifdef DEBUG

void printAll(const Vec2D<Pixel> &image)
{
	for (size_t i=0; i<image.getHeight(); i++)
	{
		for (size_t j=0; j<image.getWidth(); j++)
			std::clog << '(' << (int)image.at(i,j).r << ',' << (int)image.at(i,j).g << ',' << (int)image.at(i,j).b << ')';
		std::clog << std::endl;
	}
}

#endif

int main()
{
	Vec2D<Pixel> image = load(filename);
	printAll(image);
	return 0;
}

