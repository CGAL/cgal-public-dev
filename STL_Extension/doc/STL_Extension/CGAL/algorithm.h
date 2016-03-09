/// \defgroup STLAlgos Generic Algorithms
/// \ingroup PkgStlExtension


namespace CGAL {

/*!
\ingroup STLAlgos

\deprecated This function is deprecated, CGAL::cpp11::copy_n should be
used instead.

Copies the first `n` items from `first` to `result`.

\returns the value of `result` after inserting the `n` items.

\note The \stl release June 13, 1997, from SGI contains an equivalent
function, but it is not part of the ISO standard.

\sa `CGAL::Counting_iterator<Iterator, Value>`

copies

*/
template <class InputIterator, class Size, class
OutputIterator> OutputIterator copy_n(InputIterator first, Size n,
OutputIterator result);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup STLAlgos


Computes the minimal and the
maximal element of a range. It is modeled after the STL functions
`std::min_element` and `std::max_element`. The advantage of
`min_max_element` compared to calling both STL functions is that
one only iterates once over the sequence. This is more efficient
especially for large and/or complex sequences.

\cgalHeading{Example}

The following example program computes the minimal and
maximal element of the sequence ` (3,\,6,\,5)`. Hence the output is
`min = 3, max = 6`.

\cgalExample{STL_Extension/min_max_element_example.cpp}

\returns a pair of iterators where
the first component refers to the minimal and the second component
refers to the maximal element in the range [`first`,
`last`). The ordering is defined by `operator<` on the
value type of `ForwardIterator`.
*/
template < class ForwardIterator > std::pair<
ForwardIterator, ForwardIterator > min_max_element(ForwardIterator
first, ForwardIterator last);


/*!
\ingroup STLAlgos

Computes the minimal and the
maximal element of a range. It is modeled after the STL functions
`std::min_element` and `std::max_element`. The advantage of
`min_max_element` compared to calling both STL functions is that
one only iterates once over the sequence. This is more efficient
especially for large and/or complex sequences.


\returns a pair of iterators where the first component refers to the minimal and the
second component refers to the maximal element in the range
[`first`, `last`).

\tparam CompareMin is an adaptable binary
function object: `VT` \f$ \times\f$ `VT` \f$ \rightarrow\f$ `bool` where `VT`
is the value type of `ForwardIterator`.

\tparam CompareMax is an adaptable binary
function object: `VT` \f$ \times\f$ `VT` \f$ \rightarrow\f$ `bool` where `VT`
is the value type of `ForwardIterator`.
*/
template < class ForwardIterator, class CompareMin,
class CompareMax > std::pair< ForwardIterator, ForwardIterator >
min_max_element(ForwardIterator first, ForwardIterator last,
CompareMin comp_min, CompareMax comp_max);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup STLAlgos

\deprecated This function is deprecated. `CGAL::cpp11::prev` should be used
instead.

Returns the previous iterator,
i.e.\ the result of `operator--` on a bidirectional iterator.

\sa `CGAL::successor()`

\returns `--it`.
*/
template <class BidirectionalIterator>
BidirectionalIterator predecessor(BidirectionalIterator it);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup STLAlgos

\deprecated This function is deprecated. `CGAL::cpp11::next` should be used
instead.


Returns the next iterator, i.e.
the result of `operator++` on a forward iterator.


\sa `CGAL::predecessor()`

\returns `++it`.
*/
template <class ForwardIterator>
ForwardIterator successor(ForwardIterator it);

namespace cpp11 {

/*!
\ingroup STLAlgos

The function returns the result of `operator++` on a
`ForwardIterator`. The exact behaviour is described in Paragraph 24.4.4
of the C++ standard draft
<a href="http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2011/n3242.pdf">N3242</a>.

\note There is actually no function in namespace `CGAL::cpp11` with this
name, but a using declaration which imports a function from another
namespace. By order of priority: the one in namespace `std` is used
(provided by C++0x), if not found, then the one in namespace `boost`
is used.



\sa <a href="http://www.boost.org/doc/libs/1_46_1/libs/utility/utility.htm#functions_next_prior">boost::next</a>
\sa `CGAL::cpp11::prev()`

*/
template <typename ForwardIterator>
Iterator next(ForwardIterator it);

/*!
\ingroup STLAlgos

The function returns the result of `operator--` on
a `BidirectionalIterator`. The exact behaviour is described in
Paragraph 24.4.4 of the C++ standard draft
<a href="http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2011/n3242.pdf">N3242</a>.

\note If C++0x is available the function `std::prev` is imported into
the namespace `CGAL::cpp11`, otherwise `CGAL::cpp11::prev` is declared with the
signature as given in Paragraph 24.4.4 of the ISO C++ Standard
and forwarded to `boost::prior`.
*/
template <typename BidirectionalIterator>
Iterator prev(BidirectionalIterator it);


/*!
\ingroup STLAlgos

Copies `n` items from an
input iterator to an output iterator. Its exact behaviour is defined
in Paragraph 25.3.1 of the C++ standard draft
<a href="http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2011/n3242.pdf">N3242</a>.

\note This provides an implementation of the standard function
`copy_n` from the C++0x standard. If `copy_n` is available
in the `std::` namespace a using declaration is used, otherwise
an alternative implementation from \cgal is used.
*/

template< class InputIterator, class Size, class OutputIterator>
OutputIterator copy_n(InputIterator first, Size count, OutputIterator result);

} /* namespace cpp11 */
} /* namespace CGAL */



namespace CGAL {

/*!

\ingroup STLAlgos

@name Statistics Functions
\brief A collection of functions to simplify the process of computing properties
of geometric objects in \cgal. When required, return types of the unary
functions are computed using `CGAL::cpp11::result_of`.
Users wishing to apply the functions below to \cgal handle types should 
consider using the `CGAL::No_deref_iterator` wrapper class to wrap the
iterators being iterated over.

*/

/*!
  \ingroup STLAlgos
  Compute the mean of result of a unary function applied to an iterator
  range using double accuracy. We require that
  `to_double` may be used to convert the return type of the unary 
  function into a double.

  \param  begin  Start of iterator range.
  \param  end    End of iterator range.
  \param  f      A model of `Callable`, taking 
                 the `InputIterator` `value_type` and
                 returning a value which can be converted to double using
                 `to_double`.
  \return The mean value over values computed by the given function,
          with double accuracy.
*/
template <typename InputIterator, typename Callable>
double
mean_result(InputIterator begin, InputIterator end, Callable f);

/*!
  \ingroup STLAlgos
  Compute the maximum of all values computed using the input function over the 
  range `[begin,end)`. We require that the return type of Callable
  be computable using `CGAL::cpp11::result_of`.


  \param begin  Start of iterator range.
  \param end    End of iterator range.
  \param f      A model of `Callable` which operates on 
  the `InputIterator` `value_type` and returning a comparable
  type.
  \return The maximum value computed over the input range.
*/
template <typename InputIterator, typename Callable>
unspecified_type
max_result(InputIterator begin, InputIterator end, Callable f);

/*!
  \ingroup STLAlgos
  Compute the minimum of all values computed using the input function over the 
  range `[begin,end)`. We require that the return type of Callable
  be computable using `CGAL::cpp11::result_of`.


  \param begin  Start of iterator range.
  \param end    End of iterator range.
  \param f      A model of `Callable` which operates on 
  the `InputIterator` `value_type`, returning a comparable
  type.
  \return The minimum value computed over the input range.
*/

template <typename InputIterator, typename Callable>
unspecified_type
min_result(InputIterator begin, InputIterator end, Callable f);

/*!
  \ingroup STLAlgos
  Count all values computed using the input function within a given interval.

  \param begin  Start of iterator range.
  \param end    End of iterator range.
  \param f      A model of `Callable`.
  \param min    Smallest property value allowed in the range.
  \param max    Largest property value allowed in the range.
  \tparam Value A type which is comparable with the `InputIterator` `value_type`.

  \return Number of values computed by `f` falling in the given range.
 */

template <typename InputIterator, typename Callable, typename Value>
unspecified_type
    count_result_in_interval(InputIterator begin,
                             InputIterator end,
                             Callable f,
                             Value min,
                             Value max);

/*!
  \ingroup STLAlgos
  Compute the approximate Pearson product-moment correlation between the outputs
  of two given unary functions applied to iterator ranges. We require that
  `to_double` may be used to convert the return type into a double.

  \param begin  Start of iterator range.
  \param end    End of iterator range.
  \param f1     A model of `Callable`, to generate first set of values.
  \param f2     A model of `Callable`, to generate second set of values.
  \return
  The approximate Pearson product-moment correlation between the values 
  generated by `f1` and `f2`.

*/

template <typename InputIterator,
          typename Callable1,
          typename Callable2>
double pearson(InputIterator begin,
               InputIterator end,
               Callable1 f1,
               Callable2 f2);


} // namespace CGAL


