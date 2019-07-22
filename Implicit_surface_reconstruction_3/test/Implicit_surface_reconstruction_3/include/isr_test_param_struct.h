#include <iostream>

struct Param {
	bool octree;
	bool del_ref;
	bool poisson;
	bool spectral;
	bool march_tets;
	bool make_sm;

	void display(std::ostream &stream) const ;
};

void Param::display(std::ostream &stream) const
{
	bool first = true ;
	stream << "Used features : " ;

	if (octree) 
	{
		if (!first)
			stream << ", " ;
		stream << "octree" ;
		first = false;
	}

	if (del_ref) 
	{
		if (!first)
			stream << ", " ;
		stream << "del_ref" ;
		first = false;
	}

	if (poisson) 
	{
		if (!first)
			stream << ", " ;
		stream << "poisson" ;
		first = false;
	}

	if (spectral) 
	{
		if (!first)
			stream << ", " ;
		stream << "spectral" ;
		first = false;
	}

	if (march_tets) 
	{
		if (!first)
			stream << ", " ;
		stream << "march_tets" ;
		first = false;
	}

	if (make_sm) 
	{
		if (!first)
			stream << ", " ;
		stream << "make_sm" ;
		first = false;
	}
}

std::ostream &operator<<( std::ostream &stream, const Param &p)
{
	p.display(stream);
	return (stream);
}