#ifndef CGAL_LEVEL_OF_DETAIL_DATA_STRUCTURE_H
#define CGAL_LEVEL_OF_DETAIL_DATA_STRUCTURE_H

namespace CGAL {

	namespace Level_of_detail {

		template<class InputKernel, class InputData, class InputNormalMap, class InputLabelMap>
		struct Level_of_detail_data_structure {

		public:			
            using Kernel     = InputKernel;
			using Input      = InputData;
			using Normal_map = InputNormalMap;
			using Label_map  = InputLabelMap;

            Level_of_detail_data_structure(const Input &input, const Normal_map &normal_map, const Label_map &label_map) :
            m_input(input),
            m_normal_map(normal_map),
            m_label_map(label_map)
            { }

        private:
            const Input      &m_input;
            const Normal_map &m_normal_map;
            const Label_map  &m_label_map;
        };
    
    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_DATA_STRUCTURE_H