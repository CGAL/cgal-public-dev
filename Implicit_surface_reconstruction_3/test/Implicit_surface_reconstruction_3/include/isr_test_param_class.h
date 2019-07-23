#ifndef ISR_TEST_PARAM_CLASS_H
#define ISR_TEST_PARAM_CLASS_H

//includes
#include <iostream>
#include <list>

//Param struct
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

//Parameters class ----- /!\ parametres avec poisson + spectral desactives pour l'instant
class Parameters {

	public :
	//Constructor
	Parameters()
	{
		Param p1 = {false,  true, false, true, false, true};
		Param p2 = {false,  true, false, true, true, false};
		Param p3 = {false,  true, true, false, false, true};
		Param p4 = {false,  true, true, false, true, false};
		//Param p5 = {true,  false, false, true, false, true};
		//Param p6 = {true,  false, false, true, true, false};
		Param p7 = {true,  false, true, false, false, true};
		Param p8 = {true,  false, true, false, true, false};

		paramList.push_back(p1);
		paramList.push_back(p2);
		paramList.push_back(p3);
		paramList.push_back(p4);
		//paramList.push_back(p5);
		//paramList.push_back(p6);
		paramList.push_back(p7);
		paramList.push_back(p8);
	}


	//public methods
	void add(const Param p)
	{
		paramList.push_back(p);
	}

	std::list<Param>::const_iterator begin() const
	{
		return (paramList.begin());
	}

	std::list<Param>::const_iterator end() const
	{
		return (paramList.end());
	}

	protected :

	//Attributes
	std::list<Param> paramList;
};

#endif //ISR_TEST_PARAM_CLASS_H