#include <list>
#include "isr_test_param_struct.h"


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