#ifndef CGAL_LEVEL_OF_DETAIL_SELECTION_STRATEGY_H
#define CGAL_LEVEL_OF_DETAIL_SELECTION_STRATEGY_H

// STL includes.
#include <iostream>

// Boost includes.
#include <boost/tuple/tuple.hpp>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class InputContainer>
		class Level_of_detail_selection_strategy {

		public:
			typedef KernelTraits   Traits;
			typedef InputContainer Container;

			virtual void set_input(const Container &input) { m_input = std::make_unique<Container>(input); }
			virtual bool satisfies_condition(const int) const = 0;

			virtual ~Level_of_detail_selection_strategy() { }

		protected:
			using Label     = int;
			using Label_map = typename Container:: template Property_map<Label>;

			std::unique_ptr<Container> m_input;
		};




		template<class KernelTraits, class InputContainer>
		class Level_of_detail_clutter : public Level_of_detail_selection_strategy<KernelTraits, InputContainer> {

		public:
			typedef Level_of_detail_selection_strategy<KernelTraits, InputContainer> Base;
			
			typedef typename Base::Container Container;
			typedef typename Base::Label     Label;
			typedef typename Base::Label_map Label_map;

			void set_input(const Container &input) override {

				Base::set_input(input);
				boost::tie(m_labels, boost::tuples::ignore) = input.template property_map<Label>("label");
			}

			bool satisfies_condition(const int labelIndex) const override {

				const Label vegetation = 3; 
				if (m_labels[labelIndex] == vegetation) return true;
				return false;
			}

		private:
			Label_map m_labels;
		};




		template<class KernelTraits, class InputContainer>
		class Level_of_detail_ground : public Level_of_detail_selection_strategy<KernelTraits, InputContainer> {

		public:
			typedef Level_of_detail_selection_strategy<KernelTraits, InputContainer> Base;
			
			typedef typename Base::Container Container;
			typedef typename Base::Label     Label;
			typedef typename Base::Label_map Label_map;

			void set_input(const Container &input) override {

				Base::set_input(input);
				boost::tie(m_labels, boost::tuples::ignore) = input.template property_map<Label>("label");
			}

			bool satisfies_condition(const int labelIndex) const override {

				const Label ground = 0; 
				if (m_labels[labelIndex] == ground) return true;
				return false;
			}

		private:
			Label_map m_labels;
		};




		template<class KernelTraits, class InputContainer>
		class Level_of_detail_building_boundary : public Level_of_detail_selection_strategy<KernelTraits, InputContainer> {

		public:
			typedef Level_of_detail_selection_strategy<KernelTraits, InputContainer> Base;
			
			typedef typename Base::Container Container;
			typedef typename Base::Label     Label;
			typedef typename Base::Label_map Label_map;

			void set_input(const Container &input) override {

				Base::set_input(input);
				boost::tie(m_labels, boost::tuples::ignore) = input.template property_map<Label>("label");
			}

			bool satisfies_condition(const int labelIndex) const override {

				const Label facade = 1; 
				if (m_labels[labelIndex] == facade) return true;
				return false;
			}

		private:
			Label_map m_labels;
		};




		template<class KernelTraits, class InputContainer>
		class Level_of_detail_building_interior : public Level_of_detail_selection_strategy<KernelTraits, InputContainer> {

		public:
			typedef Level_of_detail_selection_strategy<KernelTraits, InputContainer> Base;
			
			typedef typename Base::Container Container;
			typedef typename Base::Label     Label;
			typedef typename Base::Label_map Label_map;

			void set_input(const Container &input) override {

				Base::set_input(input);
				boost::tie(m_labels, boost::tuples::ignore) = input.template property_map<Label>("label");
			}

			bool satisfies_condition(const int labelIndex) const override {

				const Label roof = 2; 
				if (m_labels[labelIndex] == roof) return true;
				return false;
			}

		private:
			Label_map m_labels;
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_SELECTION_STRATEGY_H