#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>
#include "yaml-cpp/yaml.h"
#include "boost/filesystem.hpp"
#include "analyze.h"

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
    return 0;
}

