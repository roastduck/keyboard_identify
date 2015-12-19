#ifndef ANALYZE_H
#define ANALYZE_H

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include "lodepng.h"

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
    return abs(a.r - b.r) + abs(a.g - b.g) + abs(a.b - b.b) < 20;
}

inline bool operator!=(const Pixel &a, const Pixel &b)
{
    return ! (a == b);
}

template <class T>
class Vec2D : public std::vector<T>
{
    size_t width, height;

public:
    size_t getWidth() const { return width; }
    size_t getHeight() const { return height; }

    T &at(size_t x, size_t y) { return std::vector<T>::at(x * width + y); }
    const T &at(size_t x, size_t y) const { return std::vector<T>::at(x * width + y); }

    Vec2D(size_t _width, size_t _height) :
        std::vector<T>(_width * _height), width(_width), height(_height) {}
};

Vec2D<Pixel> loadPNG(std::string filename);

struct Square
{
    int l, r, u, d;

    Square() {}
    Square(int _l, int _r, int _u, int _d)
        : l(_l), r(_r), u(_u), d(_d) {}
};

typedef std::vector< std::vector<Square> > Squares;

Squares find_squares(const Vec2D<Pixel> &image);

struct keyboard
{
    Square a;
    int kind; //1=b 2=w 3=g 4=s
};

typedef std::vector< std::vector<keyboard> > Keyboards;

Keyboards get_keyboard(const Vec2D<Pixel> &image, const Squares &mSquares);

#ifdef DEBUG
void printAll(const Vec2D<Pixel> &image);
void savePNG(const Vec2D<Pixel> &image, const char *filename);
Vec2D<Pixel> drawSquares(const Vec2D<Pixel> &image, const Squares &mSquares);
Vec2D<Pixel> drawMboard(const Vec2D<Pixel> &image, const Keyboards &Mboard);
#endif

#endif
