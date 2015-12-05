#include <cmath>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include "lodepng.h"

const char filename[] = "input.png";

typedef std::pair<int,double> PID;

template <class T1, class T2>
inline bool cmp_first(const std::pair<T1,T2> &a, const std::pair<T1,T2> &b)
{
	return a.first < b.first;
}

template <class T1, class T2>
inline bool cmp_second(const std::pair<T1,T2> &a, const std::pair<T1,T2> &b)
{
	return a.second < b.second;
}

struct Pixel // Ignore alpha
{
	int r, g, b;

	Pixel() {}
	Pixel(int _r, int _g, int _b) : 
		r(_r), g(_g), b(_b) {}
	Pixel(const Pixel &from) :
		r(from.r), g(from.g), b(from.b) {}
};

inline bool operator==(const Pixel &a, const Pixel &b)
{
	//return a.r==b.r && a.g==b.g && a.b==b.b;
	return abs(a.r-b.r)+abs(a.g-b.g)+abs(a.b-b.b)<15;
}

inline bool operator!=(const Pixel &a, const Pixel &b)
{
	return ! (a==b);
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

Vec2D<Pixel> loadPNG(const char *filename)
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

class DiffSum
{
protected:
	const Vec2D<Pixel> &image;
	Vec2D<int> diff;

	DiffSum(const Vec2D<Pixel> &_image) :
		image(_image), diff(_image.getWidth(), _image.getHeight()) {}
};

class ColDiffSum : public DiffSum
{
public:
	ColDiffSum(const Vec2D<Pixel> &_image) : DiffSum(_image)
	{
		for (size_t i=0; i<image.getHeight(); i++)
			for (size_t j=1; j<image.getWidth(); j++)
			{
				if (i) diff.at(i,j) = diff.at(i-1,j);
				diff.at(i,j) += (image.at(i,j) != image.at(i,j-1));
			}
	}

	int sum(int col, int up, int down)
	{
		return up ? diff.at(down, col)-diff.at(up-1, col) : diff.at(down, col);
	}
};

class RowDiffSum : public DiffSum
{
public:
	RowDiffSum(const Vec2D<Pixel> &_image) : DiffSum(_image)
	{
		for (size_t j=0; j<image.getWidth(); j++)
			for (size_t i=1; i<image.getHeight(); i++)
			{
				if (j) diff.at(i,j) = diff.at(i,j-1);
				diff.at(i,j) += (image.at(i,j) != image.at(i-1,j));
			}
	}

	int sum(int row, int left, int right)
	{
		return left ? diff.at(row, right)-diff.at(row,left-1) : diff.at(row, right);
	}
};

class Seperator
{
private:
	double sharpness(const std::vector<double> &sumFront, const std::vector<double> &sum2Front, int x)
	{
		double sumL = sumFront[x-1], sumR = sumFront.back()-sumL;
		double sum2L = sum2Front[x-1], sum2R = sum2Front.back()-sum2L;
		int nL = x, nR = sumFront.size()-x;
		double avgL = sumL/nL, avgR = sumR/nR;
		double sqrL = sum2L-2*sumL*avgL+nL*avgL*avgL, sqrR = sum2R-2*sumR*avgR+nR*avgR*avgR;
		return sqrL+sqrR;
	}

public:
	std::vector<PID> data;

	void filter()
	{
		sort(data.begin(), data.end(), cmp_second<int,double>);

		std::vector<double> sumFront(data.size());
		std::vector<double> sum2Front(data.size());
		for (size_t i=0; i!=data.size(); i++)
		{
			if (i)
				sumFront[i] = sumFront[i-1], sum2Front[i] = sum2Front[i-1];
			sumFront[i] += data[i].second;
			sum2Front[i] += data[i].second*data[i].second;
		}

		int l(1), r(data.size()-1), ans(0); // seperator on the right
		while (l<=r)
		{
			int midl = l+(r-l)/3, midr = r-(r-l)/3;
			if (sharpness(sumFront, sum2Front, midl) < sharpness(sumFront, sum2Front, midr))
				ans = midl, r = midr-1;
			else
				ans = midr, l = midl+1;
		}

		sort(data.begin()+ans, data.end(), cmp_first<int,double>);
		int full(0);
		for (size_t i=0, j=ans; j<data.size(); i++, j++)
		{
			size_t k = j;
			int sum = data[j].second;
			while (k+1<data.size() && data[k+1].first-data[j].first<5)
				sum += data[++k].second;
			data[i] = PID((data[j].first+data[k].first)/2, sum/(k-j+1));
			j = k;
			full = i;
		}
		data.resize(full+1);
	}
};

struct Square
{
	int l, r, u, d;

	Square(int _l, int _r, int _u, int _d)
		: l(_l), r(_r), u(_u), d(_d) {}
};

typedef std::vector< std::vector<Square> > Squares;

inline double modify_diff_row(double x)
{
	return pow(x,2);
}

inline double modify_diff_col(double x)
{
	return pow(x,3);
}

Squares find_squares(const Vec2D<Pixel> &image)
{
	RowDiffSum mRowDiffSum(image);
	ColDiffSum mColDiffSum(image);

	Seperator lineSep;
	lineSep.data.resize(image.getHeight()-1);
	for (size_t i=0; i<image.getHeight()-1; i++)
		lineSep.data[i] = PID(i+1, modify_diff_row((double)mRowDiffSum.sum(i+1, 0, image.getWidth()-1)/image.getWidth()));
	lineSep.data.push_back(PID(0, 1.0));
	lineSep.data.push_back(PID(image.getHeight(), 1.0));
	lineSep.filter();

	Squares ret(lineSep.data.size()-1);
	for (size_t i=0; i<ret.size(); i++)
	{
		int top(lineSep.data[i].first), bot(lineSep.data[i+1].first-1);

		Seperator colSep;
		colSep.data.resize(image.getWidth()-1);
		for (size_t j=0; j<image.getWidth()-1; j++)
			colSep.data[j] = PID(j+1, modify_diff_col((double)mColDiffSum.sum(j+1, top, bot)/(bot-top+1)));
		colSep.data.push_back(PID(0, 1.0));
		colSep.data.push_back(PID(image.getWidth(), 1.0));
		colSep.filter();

		for (size_t j=0; j<colSep.data.size()-1; j++)
			ret[i].push_back(Square(colSep.data[j].first, colSep.data[j+1].first-1, top, bot));
	}
	return ret;
}

#ifdef DEBUG
void printAll(const Vec2D<Pixel> &image);
void savePNG(const Vec2D<Pixel> &image, const char *filename);
Vec2D<Pixel> drawSquares(const Vec2D<Pixel> &image, const Squares &mSquares);
#endif

int main()
{
	Vec2D<Pixel> image = loadPNG(filename);
	Squares mSquares = find_squares(image);

	savePNG(drawSquares(image,mSquares),"output.png");

	return 0;
}

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

void savePNG(const Vec2D<Pixel> &image, const char *filename)
{
	std::vector<unsigned char> raw(image.getWidth() * image.getHeight() * 4);
	for (size_t i=0; i<image.getWidth()*image.getHeight(); i++)
		raw[i*4]=image[i].r, raw[i*4+1]=image[i].g, raw[i*4+2]=image[i].b, raw[i*4+3]=255;
	lodepng::encode(filename, raw, image.getWidth(), image.getHeight());
}

Vec2D<Pixel> drawSquares(const Vec2D<Pixel> &image, const Squares &mSquares)
{
	Vec2D<Pixel> ret(image);
	for (size_t i=0; i<mSquares.size(); i++)
		for (size_t j=0; j<mSquares[i].size(); j++)
		{
			std::clog << mSquares[i][j].l << ',' << mSquares[i][j].r << ',' << mSquares[i][j].u << ',' << mSquares[i][j].d << std::endl;
			for (int k=mSquares[i][j].l; k<=mSquares[i][j].r; k++)
				ret.at(mSquares[i][j].u,k) = ret.at(mSquares[i][j].d,k) = Pixel(255,0,0);
			for (int k=mSquares[i][j].u; k<=mSquares[i][j].d; k++)
				ret.at(k,mSquares[i][j].l) = ret.at(k,mSquares[i][j].r) = Pixel(255,0,0);
		}
	return ret;
}

#endif
