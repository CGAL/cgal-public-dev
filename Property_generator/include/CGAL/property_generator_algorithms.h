/******************************************************************************/
#ifndef CGAL_PROPERTY_GENERATOR_ALGORITHMS_H
#define CGAL_PROPERTY_GENERATOR_ALGORITHMS_H
/******************************************************************************/

#include <cmath>
#include <vector>
#include <iterator>
#include <CGAL/result_of.h>
#include <CGAL/metafunctions.h>

/******************************************************************************/

namespace CGAL
{

/******************************************************************************/
// Doxygen grouping for this file.

/*!

/// \defgroup STLAlgos Statistics Algorithms
/// \ingroup PkgStlExtension
@{

*/

/******************************************************************************/

// Declarations

/******************************************************************************/

/*!

@name Statistics functions
\brief A collection functions to simplify the process of computing properties
of geometric objects in \cgal.

The following functions are provided to operate on iterators and unary
functions. Each function may be applied the iterator to produce a value. When
required, return types of the unary functions are computed using
`CGAL::cpp11::result_of`. 

\cgalAdvancedBegin
These functions have been specifically designed to operate on \cgal containers
and iterators, which often give iterators that can be implicitly converted to a
handle type. For example, `Triangulation_2::Finite_vertices_iterator` has a
value type of `Triangulation_2::Vertex`. The handle for a given iterator is thus
seen as equivalent to the iterator itself. This may cause problems if a user
wishes to use algorithms in the STL, since in general STL algorithms dereference
iterators before applying functions to them. In order to deal with this
technicality, the implementation of the functions in this group will attempt
to match the value type of the iterator to the argument type of the input
function, dereferencing only when necessary. in C++11, this can be done
automatically using type meta-functions. Earlier versions of C++ can be 
used by using function objects which expose the types `argument_type` and 
`result_type` for the unary function arguments.

Users wishing to use other STL algorithms not found below should consider using
the `CGAL::STL_extension::Non_deref_iterator` wrapper class to wrap the
iterators being iterated over.

\cgalAdvancedEnd

@{

*/

/*!
  Compute the mean of values computed by a unary function over the iterator
  range. 

  To make it easier to use this function on \cgal iterator ranges,
  this function will not de-reference the iterators when calling the given
  unary function under the following conditions:
  1) in C++11: if `UnaryFunction` can not be called on
     `InputIterator::value_type`, but can be called on `InputIterator`.
  2) in C++03 if `UnaryFunction::argument_type` exists and is not convertible
     from `InputIterator::value_type`. The result type of the unary function
     must be computable by `CGAL::cpp11::result_of`.

  \param  begin  Start of iterator range.
  \param  end    End of iterator range.
  \param  f      A unary function taking either the value type of the input
                 iterator, or the iterator itself.
  \return The mean value over all the values computed by the given function,
          with double accuracy.
*/

template <typename InputIterator, typename UnaryFunction>
double mean_result(InputIterator begin, InputIterator end, UnaryFunction f);

/******************************************************************************/

/*!
  Compute the max given an iterator range.

  To make it easier to use this function on \cgal iterator ranges,
  this function will not de-reference the iterators when calling the given
  unary function under the following conditions:
  1) in C++11: if `UnaryFunction` can not be called on
     `InputIterator::value_type`, but can be called on `InputIterator`.
  2) in C++03 if `UnaryFunction::argument_type` exists and is not convertible
     from `InputIterator::value_type`. The result type of the unary function
     must be computable by `CGAL::cpp11::result_of`.

  \param begin  Start of iterator range.
  \param end    End of iterator range.
  \param f      A unary function taking either the value type of the input
                iterator, or the iterator itself.
  \return The maximum value computed over the input range.
*/

template <typename InputIterator, typename UnaryFunction>
#ifdef DOXYGEN_RUNNING
unspecified_type
#else
typename internal::get_result_type<UnaryFunction, InputIterator>::type
#endif
    max_result(InputIterator begin, InputIterator end, UnaryFunction f);

/******************************************************************************/

/*!
  Compute the minimum value of a function over an iterator range.

  To make it easier to use this function on \cgal iterator ranges,
  this function will not de-reference the iterators when calling the given
  unary function under the following conditions:
  1) in C++11: if `UnaryFunction` can not be called on
     `InputIterator::value_type`, but can be called on `InputIterator`.
  2) in C++03 if `UnaryFunction::argument_type` exists and is not convertible
     from `InputIterator::value_type`. The result type of the unary function
     must be computable by `CGAL::cpp11::result_of`.

  \param begin  Start of iterator range.
  \param end    End of iterator range.
  \param f      A unary function taking either the value type of the input
                iterator, or the iterator itself.
  \return The minimum value computed over the input range.
*/

template <typename InputIterator, typename UnaryFunction>
#ifdef DOXYGEN_RUNNING
unspecified_type
#else
typename internal::get_result_type<UnaryFunction, InputIterator>::type
#endif
    min_result(InputIterator begin, InputIterator end, UnaryFunction f);

/******************************************************************************/

/*
  Transform iterator range by applying a unary function, then write
  the output to an output iterator.
  This function works exactly as `std::transform`, except that it allows
  users to provide functions that take CGAL handles (in this case the
  iterators are not dereferenced).

  To make it easier to use this function on \cgal iterator ranges,
  this function will not de-reference the iterators when calling the given
  unary function under the following conditions:
  1) in C++11: if `UnaryFunction` can not be called on
     `InputIterator::value_type`, but can be called on `InputIterator`.
  2) in C++03 if `UnaryFunction::argument_type` exists and is not convertible
     from `InputIterator::value_type`. The result type of the unary function
     must be computable by `CGAL::cpp11::result_of`.

  \param begin  Start of iterator range.
  \param end    End of iterator range.
  \param f      A unary function.
  \param output An output iterator to write out the output.

  \return One past the end of the given output iterator.
*/

// template <typename InputIterator,
//           typename OutputIterator,
//           typename UnaryFunction>
// OutputIterator transform(InputIterator begin,
//                          InputIterator end,
//                          OutputIterator output,
//                          UnaryFunction f);

/******************************************************************************/

/*!
  Count all values computed using the input function within a given interval.

  To make it easier to use this function on \cgal iterator ranges,
  this function will not de-reference the iterators when calling the given
  unary function under the following conditions:
  1) in C++11: if `UnaryFunction` can not be called on
     `InputIterator::value_type`, but can be called on `InputIterator`.
  2) in C++03 if `UnaryFunction::argument_type` exists and is not convertible
     from `InputIterator::value_type`. The result type of the unary function
     must be computable by `CGAL::cpp11::result_of`.

  \param begin  Start of iterator range.
  \param end    End of iterator range.
  \param f      Unary function.
  \param min    Smallest property value allowed in the range.
  \param max    Largest property value allowed in the range.
  \tparam Limit A type which is comparable with `InputIterator::value_type`.

  \return Number of values computed by `f` falling in the given range.
 */

template <typename InputIterator, typename UnaryFunction, typename Limit>
typename std::iterator_traits<InputIterator>::difference_type
    count_result_in_interval(InputIterator begin,
                             InputIterator end,
                             UnaryFunction f,
                             Limit min,
                             Limit max);

/******************************************************************************/

/*!
  Compute the Pearson product-moment correlation between the outputs of two
  given unary functions applied to iterator ranges.

  To make it easier to use this function on \cgal iterator ranges,
  this function will not de-reference the iterators when calling the given
  unary function under the following conditions:
  1) in C++11: if `UnaryFunction` can not be called on
     `InputIterator::value_type`, but can be called on `InputIterator`.
  2) in C++03 if `UnaryFunction::argument_type` exists and is not convertible
     from `InputIterator::value_type`. The result type of the unary function
     must be computable by `CGAL::cpp11::result_of`.

  \param begin  Start of iterator range.
  \param end    End of iterator range.
  \param f1     First function, to generate first set of values.
  \param f2     Second function, to generate first set of values.
  \return
  The Pearson product-moment correlation between the values generated by `f1`
  and `f2`.

*/

template <typename InputIterator,
          typename UnaryFunction1,
          typename UnaryFunction2>
double pearson(InputIterator begin,
               InputIterator end,
               UnaryFunction1 f1,
               UnaryFunction2 f2);

/******************************************************************************/
// End of documentation and declarations.

} // namespace STL_extension
} // namespace CGAL

/*!

@}
@}

*/

/******************************************************************************/
// Implementations                                                            //
/******************************************************************************/
// Function to selectively dereference iterators. Limited to this file.

namespace CGAL 
{
namespace internal
{

template <typename Function, typename Iterator>
typename boost::enable_if< 
  boost::mpl::not_< do_dereference<Function,Iterator> >, 
  Iterator
>::type
conditional_dereference(Iterator i)
{
  return i;
}

template <typename Function, typename Iterator>
typename boost::enable_if< 
  do_dereference<Function,Iterator>, 
  typename std::iterator_traits<Iterator>::value_type
>::type
conditional_dereference(Iterator i)
{
  return *i;
}

}  // anonymous namespace.
namespace STL_extension
{

/******************************************************************************/

template <typename Iterator, typename Function>
double mean_result(Iterator begin, Iterator end, Function f)
{
  typename std::iterator_traits<Iterator>::difference_type count = 0;
  double sum = 0.;

  for (; begin != end; ++begin, ++count)
    sum += f( internal::conditional_dereference<Function>(begin) );

  return sum / double(count);
}

/******************************************************************************/

template <typename Iterator, typename Function>
#ifdef DOXYGEN_RUNNING
unspecified_type
#else
typename internal::get_result_type<Function, Iterator>::type
#endif
    max_result(Iterator begin, Iterator end, Function f)
{
  typedef typename internal::get_result_type<Function, Iterator>::type
      result_type;

  // Set max to value of first element.
  result_type max = f( internal::conditional_dereference<Function>(begin) );

  for (++begin; begin != end; ++begin)
  {
    result_type value = f( internal::conditional_dereference<Function>(begin) );
    if (value > max)
      max = value;
  }

  return max;
}

/******************************************************************************/

template <typename Iterator, typename Function>
#ifdef DOXYGEN_RUNNING
unspecified_type
#else
typename internal::get_result_type<Function, Iterator>::type
#endif
    min_result(Iterator begin, Iterator end, Function f)
{
  typedef typename internal::get_result_type<Function, Iterator>::type
      result_type;

  // Set max to value of first element.
  result_type min = f(internal::conditional_dereference<Function>(begin));

  // Can we optimise-out the first function call?
  for (++begin; begin != end; ++begin)
  {
    result_type value = f(internal::conditional_dereference<Function>(begin));
    if (value < min)
      min = value;
  }

  return min;
}

/******************************************************************************/

// template <typename Iterator, typename Output, typename Function>
// Output transform(Iterator begin, Iterator end, Output output, Function f)
// {
//   for (; begin != end; ++begin)
//     *(output++) = f(internal::conditional_dereference<Function>(begin));

//   return output;
// }

/******************************************************************************/

template <typename Iterator, typename Function, typename Limit>
typename std::iterator_traits<Function>::difference_type
    count_result_in_interval(Iterator begin,
                             Iterator end,
                             Function f,
                             Limit min,
                             Limit max)
{
  typedef typename internal::get_result_type<Function, Iterator>::type
      result_type;
  typename std::iterator_traits<Iterator>::difference_type count = 0;

  for (; begin != end; ++begin)
  {
    result_type value = f(internal::conditional_dereference<Function>(begin));
    if (min <= value && value <= max)
      ++count;
  }
  return count;
}

/******************************************************************************/

template <typename Iterator,
          typename Function_1,
          typename Function_2>
double pearson(Iterator begin,
               Iterator end,
               Function_1 f1,
               Function_2 f2)
{
  typedef typename internal::get_result_type<Function_1, Iterator>::type
      result_type_1;
  typedef typename internal::get_result_type<Function_2, Iterator>::type
      result_type_2;

  std::vector<result_type_1> v_1;
  std::vector<result_type_2> v_2;

  // Compute values over all of the statistics we want.
  STL_extension::transform(begin, end, std::back_inserter(v_1), f1);

  // We know how long v_1 should be now.
  v_2.reserve(v_1.size());

  STL_extension::transform(begin, end, std::back_inserter(v_2), f2);

  typename std::vector<result_type_1>::size_type n = v_1.size();
  double mean_1 = 0.;
  double mean_2 = 0.;
  double numerator = 0.;
  double denominator_1 = 0.;
  double denominator_2 = 0.;

  // Compute means
  for (unsigned i = 0; i != n; ++i)
  {
    mean_1 += v_1[i];
    mean_2 += v_2[i];
  }

  mean_1 = mean_1 / n;
  mean_2 = mean_2 / n;

  for (int i = 0; i < n; i++)
  {
    numerator += (v_1[i] - mean_1) * (v_2[i] - mean_2);
    denominator_1 += std::pow((v_1[i] - mean_1), 2);
    denominator_2 += std::pow((v_2[i] - mean_2), 2);
  }

  return numerator / (std::sqrt(denominator_1) * std::sqrt(denominator_2));
}

/******************************************************************************/

}  // namespace CGAL

/******************************************************************************/
#endif
/******************************************************************************/
