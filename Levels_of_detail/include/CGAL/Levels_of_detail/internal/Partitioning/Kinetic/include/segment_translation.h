#pragma once
#include "defs.h"
#include "defs_cgal.h"


namespace Skippy {

	typedef enum {
		INSTANTANEOUS,
		PROGRESSIVE
	} Translation_Type;


	class Segment_Translation
	{
	public:
		Segment_Translation(const FT & t, const CGAL_Point_2 & _A, const CGAL_Vector_2 & _dA);

		Segment_Translation(const FT & _t, const CGAL_Point_2 & _A, const CGAL_Point_2 & _B);
		
		~Segment_Translation();

		void set_end(const FT & t);

	public:
		Translation_Type type;
		FT t_int_start;
		FT t_int_end;
		CGAL_Point_2 A;
		CGAL_Vector_2 dA;
		CGAL_Point_2 B;
	};

#if 0
	class Segment_Translation
	{
	protected:
		Segment_Translation(const Translation_Type & _type);

	public:
		virtual ~Segment_Translation();

	public:
		virtual const FT & get_t_begin() const = 0;

		virtual const FT & get_t_end() const = 0;

		virtual const CGAL_Point_2 & get_begin() const = 0;

		virtual CGAL_Point_2 get_end() const = 0;

	public:
		const Translation_Type type;
	};

	class Segment_Translation_P : public Segment_Translation
	{
	public:
		Segment_Translation_P(const FT & t, const CGAL_Point_2 & _A, const CGAL_Vector_2 & _dA);

		void set_end(const FT & t);

	public:
		FT t_int_start;
		FT t_int_end;
		CGAL_Point_2 A;
		CGAL_Vector_2 dA;
	};

	class Segment_Translation_I : public Segment_Translation
	{
	public:
		Segment_Translation_I(const FT & _t, const CGAL_Point_2 & _A, const CGAL_Point_2 & _B);

		~Segment_Translation_I();

	public:
		FT t;
		CGAL_Point_2 A;
		CGAL_Point_2 B;
	};
#endif
}