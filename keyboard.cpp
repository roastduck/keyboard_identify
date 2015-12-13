#include <cmath>
#include <cstdio>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>
#include "lodepng.h"
#include "yaml-cpp/yaml.h"
#include "boost/filesystem.hpp"

typedef std::pair<int, double> PID;

template <class T1, class T2>
inline bool cmp_first(const std::pair<T1, T2> &a, const std::pair<T1, T2> &b)
{
    return a.first < b.first;
}

template <class T1, class T2>
inline bool cmp_second(const std::pair<T1, T2> &a, const std::pair<T1, T2> &b)
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

Vec2D<Pixel> loadPNG(std::string filename)
{
    unsigned width, height;
    std::vector<unsigned char> raw;
    unsigned error = lodepng::decode(raw, width, height, filename);
    if (error)
        std::cerr << "Error reading " << filename << " : " << lodepng_error_text(error) << std::endl;
    Vec2D<Pixel> ret(width, height);
    for (size_t i = 0; i * 4 < raw.size(); i++)
        ret[i] = Pixel(raw[i * 4], raw[i * 4 + 1], raw[i * 4 + 2]);
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
        for (size_t i = 0; i < image.getHeight(); i++)
            for (size_t j = 1; j < image.getWidth(); j++)
            {
                if (i) diff.at(i, j) = diff.at(i - 1, j);
                diff.at(i, j) += (image.at(i, j) != image.at(i, j - 1));
            }
    }

    int sum(int col, int up, int down)
    {
        return up ? diff.at(down, col) - diff.at(up - 1, col) : diff.at(down, col);
    }
};

class RowDiffSum : public DiffSum
{
public:
    RowDiffSum(const Vec2D<Pixel> &_image) : DiffSum(_image)
    {
        for (size_t j = 0; j < image.getWidth(); j++)
            for (size_t i = 1; i < image.getHeight(); i++)
            {
                if (j) diff.at(i, j) = diff.at(i, j - 1);
                diff.at(i, j) += (image.at(i, j) != image.at(i - 1, j));
            }
    }

    int sum(int row, int left, int right)
    {
        return left ? diff.at(row, right) - diff.at(row, left - 1) : diff.at(row, right);
    }
};

class Seperator
{
private:
    double sharpness(const std::vector<double> &sumFront, const std::vector<double> &sum2Front, int x)
    {
        double sumL = sumFront[x - 1], sumR = sumFront.back() - sumL;
        double sum2L = sum2Front[x - 1], sum2R = sum2Front.back() - sum2L;
        int nL = x, nR = sumFront.size() - x;
        double avgL = sumL / nL, avgR = sumR / nR;
        double sqrL = sum2L - 2 * sumL * avgL + nL * avgL * avgL, sqrR = sum2R - 2 * sumR * avgR + nR * avgR * avgR;
        return sqrL + sqrR;
    }

public:
    std::vector<PID> data;

    void filter()
    {
        sort(data.begin(), data.end(), cmp_second<int, double>);

        std::vector<double> sumFront(data.size());
        std::vector<double> sum2Front(data.size());
        for (size_t i = 0; i != data.size(); i++)
        {
            if (i)
                sumFront[i] = sumFront[i - 1], sum2Front[i] = sum2Front[i - 1];
            sumFront[i] += data[i].second;
            sum2Front[i] += data[i].second * data[i].second;
        }

        int l(1), r(data.size() - 1), ans(0); // seperator on the right
        while (l <= r)
        {
            int midl = l + (r - l) / 3, midr = r - (r - l) / 3;
            if (sharpness(sumFront, sum2Front, midl) < sharpness(sumFront, sum2Front, midr))
                ans = midl, r = midr - 1;
            else
                ans = midr, l = midl + 1;
        }

        sort(data.begin() + ans, data.end(), cmp_first<int, double>);
        int full(0);
        for (size_t i = 0, j = ans; j < data.size(); i++, j++)
        {
            size_t k = j;
            int sum = data[j].second;
            while (k + 1 < data.size() && data[k + 1].first - data[j].first < 5)
                sum += data[++k].second;
            data[i] = PID((data[j].first + data[k].first) / 2, sum / (k - j + 1));
            j = k;
            full = i;
        }
        data.resize(full + 1);
    }
};

struct Square
{
    int l, r, u, d;

    Square() {}
    Square(int _l, int _r, int _u, int _d)
        : l(_l), r(_r), u(_u), d(_d) {}
};

typedef std::vector< std::vector<Square> > Squares;

inline double modify_diff_row(double x)
{
    return pow(x, 2);
}

inline double modify_diff_col(double x)
{
    return pow(x, 3);
}

Squares find_squares(const Vec2D<Pixel> &image)
{
    RowDiffSum mRowDiffSum(image);
    ColDiffSum mColDiffSum(image);

    Seperator lineSep;
    lineSep.data.resize(image.getHeight() - 1);
    for (size_t i = 0; i < image.getHeight() - 1; i++)
        lineSep.data[i] = PID(i + 1, modify_diff_row((double)mRowDiffSum.sum(i + 1, 0, image.getWidth() - 1) / image.getWidth()));
    lineSep.data.push_back(PID(0, 1.0));
    lineSep.data.push_back(PID(image.getHeight(), 1.0));
    lineSep.filter();

    Squares ret(lineSep.data.size() - 1);
    for (size_t i = 0; i < ret.size(); i++)
    {
        int top(lineSep.data[i].first), bot(lineSep.data[i + 1].first - 1);

        Seperator colSep;
        colSep.data.resize(image.getWidth() - 1);
        for (size_t j = 0; j < image.getWidth() - 1; j++)
            colSep.data[j] = PID(j + 1, modify_diff_col((double)mColDiffSum.sum(j + 1, top, bot) / (bot - top + 1)));
        colSep.data.push_back(PID(0, 1.0));
        colSep.data.push_back(PID(image.getWidth(), 1.0));
        colSep.filter();

        for (size_t j = 0; j < colSep.data.size() - 1; j++)
            ret[i].push_back(Square(colSep.data[j].first, colSep.data[j + 1].first - 1, top, bot));
    }
    return ret;
}

struct keyboard
{
    Square a;
    int kind; //1=b 2=w 3=g 4=s
};

typedef std::vector< std::vector<keyboard> > Keyboards;

Pixel get_most(const Vec2D<Pixel> &image, const Square &a)
{
    int sum[300], Max = 0, R;
    memset(sum, 0, sizeof(sum));
    for (int i = a.u; i <= a.d; i++)
        for (int j = a.l; j <= a.r; j++)
        {
            int color = image.at(i, j).r;
            sum[color]++;
            if (sum[color] > Max)
            {
                Max = sum[color];
                R = color;
            }
        }

    int G;
    Max = 0;
    memset(sum, 0, sizeof(sum));
    for (int i = a.u; i <= a.d; i++)
        for (int j = a.l; j <= a.r; j++)
        {
            int color = image.at(i, j).g;
            sum[color]++;
            if (sum[color] > Max)
            {
                Max = sum[color];
                G = color;
            }
        }

    int B;
    Max = 0;
    memset(sum, 0, sizeof(sum));
    for (int i = a.u; i <= a.d; i++)
        for (int j = a.l; j <= a.r; j++)
        {
            int color = image.at(i, j).b;
            sum[color]++;
            if (sum[color] > Max)
            {
                Max = sum[color];
                B = color;
            }
        }
    return (Pixel) { R, G, B };
}

int get_kind(const Vec2D<Pixel> &image, const Square &a)
{
    Pixel most = get_most(image, a);
    int tot = 0, all = 0;
    for (int i = a.u + 2; i <= a.d - 2; i++)
        for (int j = a.l + 2; j <= a.r - 2; j++)
        {
            if (most == image.at(i, j)) tot++;
            all++;
        }
    if (tot >= all - 10) return 4;
    if (abs(most.r - most.g) <= 10 && abs(most.r - most.b) <= 10 && abs(most.b - most.g) <= 10 && most.r >= 230) return 2;
    if (abs(most.r - most.g) <= 40 && abs(most.r - most.b) <= 40 && abs(most.b - most.g) <= 40) return 3;
    return 1;
}

bool samecolor(const Vec2D<Pixel> &image, const Square &a, const Square &b)
{
    Pixel mosta = get_most(image, a), mostb = get_most(image, b);
    return mosta == mostb;
}

Keyboards get_keyboard(const Vec2D<Pixel> &image, const Squares &mSquares)
{
    Keyboards Mboard(mSquares.size());
    int now = mSquares.size() - 1, tot = -1;
    for (int k = 1; k <= 3; k++)
    {
        while (mSquares[now].size() == 1)
        {
            Mboard[++tot].push_back((keyboard){mSquares[now][0], 4});
            now--;
        }
        tot++;
        for (int i = 0; i < mSquares[now].size(); i++)
            Mboard[tot].push_back((keyboard){mSquares[now][i], get_kind(image, mSquares[now][i])});
        now--;
    }
    Square rcolor = Mboard[tot][0].a;
    int need = 0,Len=rcolor.d-rcolor.u;
    for (int i = 0; i < Mboard[tot].size(); i++)
        if (Mboard[tot][i].kind != 4) need++;
    need /= 2;
    for (int i = 0; i < Mboard[tot].size(); i++)
        if (Mboard[tot][i].kind != 4)
        {
            need--;
            if (need == 0) rcolor = Mboard[tot][i].a;
        }
    while (now >= 0)
    {
        if (now>=0&&mSquares[now].size() == 1)
        {
            Mboard[++tot].push_back((keyboard){mSquares[now][0], 4});
            now--;
        }
        if (now<0) break;
        Square nowcolor = Mboard[tot][0].a;
        need = 0;
        for (int i = 0; i < mSquares[now].size(); i++)
            if (get_kind(image, mSquares[now][i]) != 4) need++;
		if (need==0) break;
        need /= 2;
        bool have=0;
        for (int i = 0; i < mSquares[now].size(); i++)
            if (get_kind(image, mSquares[now][i]) != 4)
            {
                need--;
                if (need == 0) {nowcolor = mSquares[now][i];have=1;}
            }
        if (!have) break;
        if (have&&!samecolor(image, rcolor, nowcolor)) break;
        tot++;
        for (int i = 0; i < mSquares[now].size(); i++)
            Mboard[tot].push_back((keyboard){mSquares[now][i], get_kind(image, mSquares[now][i])});
        now--;
    }
	while (Mboard.back().size()==0) Mboard.pop_back();
    return Mboard;
}

struct Image
{
	std::vector<std::string> path; // from leaf to curdir
	std::string fullpath, layout, category;
};
typedef std::vector<Image> Files;

inline std::vector<std::string> split(std::string s)
{
	std::vector<std::string> ret(1);
	for (std::string::iterator i=s.begin(); i!=s.end(); i++)
		if ((*i=='_' || *i=='.') && !ret.back().empty())
			ret.push_back("");
		else
			ret.back() += *i;
	if (ret.back().empty()) ret.pop_back();
	return ret;
}

Files search_images(boost::filesystem::path path = boost::filesystem::current_path())
{
	using namespace boost::filesystem;
	if (!exists(path)) return Files();
	Files ret;
	for (directory_iterator i(path); i!=directory_iterator(); i++)
	{
		if (is_directory(i->status()))
		{
			Files deeper = search_images(i->path());
			for (Files::iterator j=deeper.begin(); j!=deeper.end(); j++)
			{
				ret.push_back(*j);
				ret.back().path.push_back(i->path().filename().string());
			} 
		} else
		{
			std::vector<std::string> part = split(i->path().filename().string());
			if (part.back() != "png" || part.size() != 3) continue;
			Image mImage;
			mImage.layout = part[0];
			mImage.category = part[1];
			mImage.fullpath = i->path().string();
			ret.push_back(mImage);
		}
	}
	return ret;
}

struct FrameBody { int x, y, w, h; };
struct Frame { FrameBody body; };
typedef std::vector<Frame> Keys;
struct Row { FrameBody frame; Keys keys; };
typedef std::vector<Row> RowsBody;
struct Rows { RowsBody body; };
typedef std::map< std::string, std::map< std::string, Rows> > Layouts;
struct Result { int width, height; Layouts layouts; };
typedef std::map< std::string, std::map< std::string, Result> > Results;

YAML::Emitter &operator<<(YAML::Emitter &os, const FrameBody &mFrame)
{
	return os << YAML::Flow << YAML::BeginSeq << mFrame.x << mFrame.y << mFrame.w << mFrame.h << YAML::EndSeq;
}

YAML::Emitter &operator<<(YAML::Emitter &os, const Frame &mFrame)
{
	return os << YAML::BeginMap << YAML::Key << "frame" << YAML::Value << mFrame.body << YAML::EndMap;
}

YAML::Emitter &operator<<(YAML::Emitter &os, const Row &mRow)
{
	os << YAML::BeginMap;
	os << YAML::Key << "frame" << YAML::Value << mRow.frame;
	os << YAML::Key << "keys" << YAML::Value << mRow.keys;
	os << YAML::EndMap;
	return os;
}

YAML::Emitter &operator<<(YAML::Emitter &os, const Rows &mRows)
{
	return os << YAML::BeginMap << YAML::Key << "rows" << YAML::Value << mRows.body << YAML::EndMap;
}

YAML::Emitter &operator<<(YAML::Emitter &os, const Result &mResult)
{
	os << YAML::BeginMap;
	os << YAML::Key << "width" << YAML::Value << mResult.width;
	os << YAML::Key << "height" << YAML::Value << mResult.height;
	os << YAML::Key << "layouts" << YAML::Value << mResult.layouts;
	os << YAML::EndMap;
	return os;
}

inline int str_to_int(const std::string &s)
{
	std::istringstream ss(s);
	int ret;
	ss >> ret;
	return ret;
}

#ifdef DEBUG
void printAll(const Vec2D<Pixel> &image);
void savePNG(const Vec2D<Pixel> &image, const char *filename);
Vec2D<Pixel> drawSquares(const Vec2D<Pixel> &image, const Squares &mSquares);
Vec2D<Pixel> drawMboard(const Vec2D<Pixel> &image, const Keyboards &Mboard);
#endif

int main(int argc, char **argv)
{
	std::string path(argc>1 ? argv[1] : ".");
	Files mFiles = search_images(path);
	Results mResults;
	for (Files::iterator i=mFiles.begin(); i!=mFiles.end(); i++)
	{
		Vec2D<Pixel> image = loadPNG(i->fullpath);
		Squares mSquares = find_squares(image);
		Keyboards Mboard = get_keyboard(image, mSquares);
		double point_pixel = (double)str_to_int(i->path.front()) / image.getWidth();

		Result &mResult = mResults[i->path.back()][i->path.front()];
		RowsBody &mRows = mResult.layouts[i->layout][i->category].body;
		
		mResult.width = str_to_int(i->path.front());
		mResult.height = point_pixel*(image.getHeight() - Mboard.back().front().a.u + 1);

		for (Keyboards::reverse_iterator j=Mboard.rbegin(); j!=Mboard.rend(); j++)
			if (!(j->size()==1 && j->front().kind==4))
			{
				Row ret;
				int leftmost = -1, rightmost = -1;
				for (std::vector<keyboard>::iterator k=j->begin(); k!=j->end(); k++)
					if (k->kind != 4)
					{
						if (leftmost == -1) leftmost = k->a.l;
						rightmost = k->a.r;
						
						ret.keys.push_back(Frame());
						ret.keys.back().body = (FrameBody) {
							int(point_pixel*(k->a.l - leftmost)),
							0,
							int(point_pixel*(k->a.r - k->a.l + 1)),
							int(point_pixel*(k->a.d - k->a.u + 1 + 1)) // plus one to match the example
						};
					}
				ret.frame = (FrameBody) {
					int(point_pixel*leftmost),
					int(point_pixel*(j->front().a.u - Mboard.back().front().a.u)),
					int(point_pixel*(rightmost - leftmost + 1)),
					int(point_pixel*(j->front().a.d - j->front().a.u + 1 + 1)) //plus one to match the example
				};
				mRows.push_back(ret);
			}
		std::clog << "Analyzed " << i->fullpath << std::endl;
	}
	for (Results::iterator i=mResults.begin(); i!=mResults.end(); i++)
		for (std::map<std::string,Result>::iterator j=i->second.begin(); j!=i->second.end(); j++)
		{
			YAML::Emitter em;
			em << j->second;
			std::string outputFile = path+"/"+i->first+"/"+j->first+"/layout_"+j->first+".yaml";
			std::ofstream os(outputFile.c_str());
			os << em.c_str();
			os.close();
			std::clog << "Generated " << outputFile << std::endl;
		}
    //savePNG(drawSquares(image,mSquares),"output.png");
    //savePNG(drawMboard(image, Mboard), "output.png");
    return 0;
}

#ifdef DEBUG

void printAll(const Vec2D<Pixel> &image)
{
    for (size_t i = 0; i < image.getHeight(); i++)
    {
        for (size_t j = 0; j < image.getWidth(); j++)
            std::clog << '(' << (int)image.at(i, j).r << ',' << (int)image.at(i, j).g << ',' << (int)image.at(i, j).b << ')';
        std::clog << std::endl;
    }
}

void savePNG(const Vec2D<Pixel> &image, const char *filename)
{
    std::vector<unsigned char> raw(image.getWidth() * image.getHeight() * 4);
    for (size_t i = 0; i < image.getWidth()*image.getHeight(); i++)
        raw[i * 4] = image[i].r, raw[i * 4 + 1] = image[i].g, raw[i * 4 + 2] = image[i].b, raw[i * 4 + 3] = 255;
    lodepng::encode(filename, raw, image.getWidth(), image.getHeight());
}
Vec2D<Pixel> drawSquares(const Vec2D<Pixel> &image, const Squares &mSquares)
{
    Vec2D<Pixel> ret(image);
    for (size_t i = 0; i < mSquares.size(); i++)
        for (size_t j = 0; j < mSquares[i].size(); j++)
        {
            std::clog << mSquares[i][j].l << ',' << mSquares[i][j].r << ',' << mSquares[i][j].u << ',' << mSquares[i][j].d << std::endl;
            for (int k = mSquares[i][j].l; k <= mSquares[i][j].r; k++)
                ret.at(mSquares[i][j].u, k) = ret.at(mSquares[i][j].d, k) = Pixel(255, 0, 0);
            for (int k = mSquares[i][j].u; k <= mSquares[i][j].d; k++)
                ret.at(k, mSquares[i][j].l) = ret.at(k, mSquares[i][j].r) = Pixel(255, 0, 0);
        }
    return ret;
}
Vec2D<Pixel> drawMboard(const Vec2D<Pixel> &image, const Keyboards &Mboard)
{
    Vec2D<Pixel> ret(image);
    for (size_t i = 0; i < Mboard.size(); i++)
        for (size_t j = 0; j < Mboard[i].size(); j++)
        {
            std::clog << Mboard[i][j].a.l << ',' << Mboard[i][j].a.r << ',' << Mboard[i][j].a.u << ',' << Mboard[i][j].a.d << std::endl;
            if (Mboard[i][j].kind == 1)
            {
                for (int k = Mboard[i][j].a.l; k <= Mboard[i][j].a.r; k++)
                    ret.at(Mboard[i][j].a.u + 1, k) = ret.at(Mboard[i][j].a.d - 1, k) = Pixel(124, 252, 0);
                for (int k = Mboard[i][j].a.u; k <= Mboard[i][j].a.d; k++)
                    ret.at(k, Mboard[i][j].a.l + 1) = ret.at(k, Mboard[i][j].a.r - 1) = Pixel(124, 252, 0);
            }
            if (Mboard[i][j].kind == 2)
            {
                for (int k = Mboard[i][j].a.l; k <= Mboard[i][j].a.r; k++)
                    ret.at(Mboard[i][j].a.u + 1, k) = ret.at(Mboard[i][j].a.d - 1, k) = Pixel(184, 134, 11);
                for (int k = Mboard[i][j].a.u; k <= Mboard[i][j].a.d; k++)
                    ret.at(k, Mboard[i][j].a.l + 1) = ret.at(k, Mboard[i][j].a.r - 1) = Pixel(184, 134, 11);
            }
            if (Mboard[i][j].kind == 3)
            {
                for (int k = Mboard[i][j].a.l; k <= Mboard[i][j].a.r; k++)
                    ret.at(Mboard[i][j].a.u + 1, k) = ret.at(Mboard[i][j].a.d - 1, k) = Pixel(255, 0, 0);
                for (int k = Mboard[i][j].a.u; k <= Mboard[i][j].a.d; k++)
                    ret.at(k, Mboard[i][j].a.l + 1) = ret.at(k, Mboard[i][j].a.r - 1) = Pixel(255, 0, 0);
            }
            if (Mboard[i][j].kind == 4)
            {
                for (int k = Mboard[i][j].a.l; k <= Mboard[i][j].a.r; k++)
                    ret.at(Mboard[i][j].a.u + 1, k) = ret.at(Mboard[i][j].a.d - 1, k) = Pixel(25, 25, 112);
                for (int k = Mboard[i][j].a.u; k <= Mboard[i][j].a.d; k++)
                    ret.at(k, Mboard[i][j].a.l + 1) = ret.at(k, Mboard[i][j].a.r - 1) = Pixel(25, 25, 112);
            }
        }
    return ret;
}
#endif
