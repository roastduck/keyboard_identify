#include <string>
#include <fstream>
#include <iostream>
#include <iterator>
#include "analyze.h"

const char listfile[] = "list.txt";

static inline char colorname(int x)
{
	switch (x)
	{
		case 1: return 'B';
		case 2: return 'W';
		case 3: return 'G';
		case 4: return 'S';
		default: std::cerr << "colorname error" << std::endl;
	}
}

static inline std::string outputfile(std::string s)
{
	size_t len = s.length();
	s[len-1] = 't';
	s[len-2] = 'x';
	s[len-3] = 't';
	return s;
}

int main()
{
	std::ifstream listfs(listfile);
	std::string PNGfile;
	while (listfs >> PNGfile)
	{
		Vec2D<Pixel> image = loadPNG(PNGfile);
		Squares mSquares = find_squares(image);
		Keyboards Mboard = get_keyboard(image, mSquares);

		std::ofstream outfs(outputfile(PNGfile).c_str());
		outfs << Mboard.size() << std::endl;
		for (Keyboards::reverse_iterator i=Mboard.rbegin(); i!=Mboard.rend(); i++)
		{
			int height = i->front().a.getHeight();
			if (height&1) height += i->size()==1 ? -1 : 1;
			outfs << height;
			if (i+1==Mboard.rend()) outfs << std::endl; else outfs << ' ';
		}
		for (Keyboards::iterator i=Mboard.begin(); i!=Mboard.end(); i++)
			if (i->size() == 1)
				i->clear();
		for (Keyboards::reverse_iterator i=Mboard.rbegin(); i!=Mboard.rend(); i++)
		{
			outfs << i->size();
			if (i+1==Mboard.rend()) outfs << std::endl; else outfs << ' ';
		}
		for (Keyboards::reverse_iterator i=Mboard.rbegin(); i!=Mboard.rend(); i++)
		{
			for (std::vector<keyboard>::iterator j=i->begin(); j!=i->end(); j++)
			{
				outfs << j->a.getWidth() << ' ' << colorname(j->kind);
				if (j+1!=i->end()) outfs << ' ';
			}
			outfs << std::endl;
		}
		outfs.close();
		//std::clog << "Generated " << PNGfile << std::endl;
	}
	return 0;
}

